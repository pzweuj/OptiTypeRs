use std::collections::HashSet;

use anyhow::{Result, ensure};
use serde::{Deserialize, Serialize};

pub(crate) const SEED: usize = 15;
const POSITION_BITS: usize = 14;
const HIT_BITS: usize = 34;
const POSITION_MASK: u64 = (1 << POSITION_BITS) - 1;
const HIT_MASK: u64 = (1 << HIT_BITS) - 1;
const SEED_MASK: u64 = (1 << (2 * SEED)) - 1;

#[derive(Serialize, Deserialize)]
pub(crate) struct Segment {
    pub bases: Vec<u8>,
    pub alleles: Vec<usize>,
}

#[derive(Serialize, Deserialize)]
pub(crate) struct Index {
    pub segments: Vec<Segment>,
    postings: Vec<u64>,
}

fn base_code(base: u8) -> Option<u64> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn seeds(bases: &[u8]) -> impl Iterator<Item = (usize, u64)> + '_ {
    let mut key = 0;
    let mut length = 0;
    bases.iter().enumerate().filter_map(move |(i, &base)| {
        if let Some(code) = base_code(base) {
            key = ((key << 2) | code) & SEED_MASK;
            length += 1;
            if length >= SEED {
                return Some((i + 1 - SEED, key));
            }
        } else {
            length = 0;
        }
        None
    })
}

impl Index {
    pub fn build(segments: Vec<Segment>) -> Result<Self> {
        ensure!(
            segments.len() <= 1 << (HIT_BITS - POSITION_BITS),
            "too many reference segments"
        );
        let mut postings = Vec::new();
        for (id, segment) in segments.iter().enumerate() {
            ensure!(
                segment.bases.len() <= 1 << POSITION_BITS,
                "reference segment is too long"
            );
            let segment_key = (id as u64) << POSITION_BITS;
            postings.extend(
                seeds(&segment.bases)
                    .map(|(pos, key)| (key << HIT_BITS) | segment_key | pos as u64),
            );
        }
        postings.sort_unstable();
        Ok(Self { segments, postings })
    }

    fn postings(&self, key: u64) -> &[u64] {
        let start = self.postings.partition_point(|hit| hit >> HIT_BITS < key);
        let end = self.postings[start..].partition_point(|hit| hit >> HIT_BITS == key);
        &self.postings[start..start + end]
    }

    pub fn best_hits(&self, forward: &[u8]) -> Vec<usize> {
        let reverse: Vec<_> = forward
            .iter()
            .rev()
            .map(|base| match base {
                b'A' => b'T',
                b'T' => b'A',
                b'G' => b'C',
                b'C' => b'G',
                _ => b'N',
            })
            .collect();
        let max_errors = forward.len() * 3 / 100;
        let unknown = forward
            .iter()
            .filter(|&&base| base_code(base).is_none())
            .count();
        let mut hit_segments = HashSet::new();
        let mut candidates = HashSet::new();
        let mut scratch = [Vec::new(), Vec::new()];

        for errors in unknown..=max_errors {
            for read in [forward, reverse.as_slice()] {
                candidates.clear();
                if read.len() < SEED * (errors + 1) {
                    for (id, segment) in self.segments.iter().enumerate() {
                        if hit_segments.contains(&id) {
                            continue;
                        }
                        let min_len = read.len() - errors;
                        if let Some(last_start) = segment.bases.len().checked_sub(min_len) {
                            if (0..=last_start).any(|start| {
                                within_distance(read, &segment.bases[start..], errors, &mut scratch)
                            }) {
                                hit_segments.insert(id);
                            }
                        }
                    }
                    continue;
                }

                // An alignment with d edits leaves one of d+1 disjoint blocks intact.
                for part in 0..=errors {
                    let start = part * read.len() / (errors + 1);
                    let end = (part + 1) * read.len() / (errors + 1);
                    let chosen = seeds(&read[start..end])
                        .step_by(4)
                        .map(|(offset, key)| (start + offset, self.postings(key)))
                        .try_fold(None, |best: Option<(usize, &[u64])>, candidate| {
                            if candidate.1.is_empty() {
                                Err(candidate)
                            } else if best.is_none_or(|best| candidate.1.len() < best.1.len()) {
                                Ok(Some(candidate))
                            } else {
                                Ok(best)
                            }
                        })
                        .unwrap_or_else(Some);
                    let Some((offset, hits)) = chosen else {
                        continue;
                    };
                    for &hit in hits {
                        let id = ((hit & HIT_MASK) >> POSITION_BITS) as usize;
                        if hit_segments.contains(&id) {
                            continue;
                        }
                        let pos = (hit & POSITION_MASK) as usize;
                        let bases = &self.segments[id].bases;
                        let Some(upper) = (pos + errors).checked_sub(offset) else {
                            continue;
                        };
                        let lower = pos.saturating_sub(offset + errors);
                        for alignment_start in lower..=upper.min(bases.len()) {
                            if (errors == 0 || candidates.insert((id, alignment_start)))
                                && within_distance(
                                    read,
                                    &bases[alignment_start..],
                                    errors,
                                    &mut scratch,
                                )
                            {
                                hit_segments.insert(id);
                                break;
                            }
                        }
                    }
                }
            }
            if !hit_segments.is_empty() {
                break;
            }
        }
        let mut alleles: Vec<_> = hit_segments
            .into_iter()
            .flat_map(|id| self.segments[id].alleles.iter().copied())
            .collect();
        alleles.sort_unstable();
        alleles.dedup();
        alleles
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Independent full-matrix oracle with free leading and trailing reference gaps.
    fn distance(read: &[u8], reference: &[u8]) -> usize {
        let mut previous = vec![0; reference.len() + 1];
        for (i, &base) in read.iter().enumerate() {
            let mut row = vec![i + 1];
            for (j, &target) in reference.iter().enumerate() {
                row.push(
                    (previous[j] + usize::from(base != target || base == b'N'))
                        .min(previous[j + 1] + 1)
                        .min(row[j] + 1),
                );
            }
            previous = row;
        }
        *previous.iter().min().unwrap()
    }

    #[test]
    fn matches_unseeded_alignment_on_both_strands() {
        let mut rng = 42_u64;
        let mut dna = |len| {
            (0..len)
                .map(|_| {
                    rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                    b"ACGT"[(rng >> 62) as usize]
                })
                .collect::<Vec<_>>()
        };
        let sequence = dna(180);
        let mut neighbor = sequence.clone();
        neighbor[65] = b'A';
        neighbor.remove(93);
        let index = Index::build(
            [sequence.clone(), neighbor, dna(200)]
                .into_iter()
                .enumerate()
                .map(|(id, bases)| Segment {
                    bases,
                    alleles: vec![id],
                })
                .collect(),
        )
        .unwrap();
        for length in [12, 34, 67, 100, 125] {
            for start in [0, 11, 24, 48] {
                let source = &sequence[start..start + length];
                let mut edited = source.to_vec();
                edited[length / 3] = b'N';
                edited.insert(length / 2, b'T');
                edited.remove(length - 1);
                for read in [source.to_vec(), edited, dna(length)] {
                    let reverse: Vec<_> = read
                        .iter()
                        .rev()
                        .map(|b| match b {
                            b'A' => b'T',
                            b'C' => b'G',
                            b'G' => b'C',
                            b'T' => b'A',
                            _ => b'N',
                        })
                        .collect();
                    let scores: Vec<_> = index
                        .segments
                        .iter()
                        .map(|segment| {
                            distance(&read, &segment.bases).min(distance(&reverse, &segment.bases))
                        })
                        .collect();
                    let best = *scores.iter().min().unwrap();
                    let expected: Vec<_> = scores
                        .iter()
                        .enumerate()
                        .filter(|(_, score)| **score == best && best <= length * 3 / 100)
                        .map(|(id, _)| id)
                        .collect();
                    assert_eq!(index.best_hits(&read), expected, "read {read:?}");
                    assert_eq!(index.best_hits(&reverse), expected, "reverse {reverse:?}");
                }
            }
        }
    }
}

// Align the entire read to any prefix of this reference suffix.
fn within_distance(
    read: &[u8],
    reference: &[u8],
    limit: usize,
    scratch: &mut [Vec<usize>; 2],
) -> bool {
    if limit == 0 {
        return reference.starts_with(read);
    }
    let n = read.len();
    let m = reference.len().min(n + limit);
    if m < n - limit {
        return false;
    }
    let outside = limit + 1;
    let [previous, current] = scratch;
    previous.resize(m + 1, outside);
    current.resize(m + 1, outside);
    for (j, value) in previous[..=limit.min(m)].iter_mut().enumerate() {
        *value = j;
    }
    if limit < m {
        previous[limit + 1] = outside;
    }
    for (i, &base) in read.iter().enumerate() {
        let i = i + 1;
        let start = i.saturating_sub(limit).max(1);
        let end = (i + limit).min(m);
        current[start - 1] = if start == 1 { i.min(outside) } else { outside };
        let mut best = current[start - 1];
        for j in start..=end {
            let mismatch = usize::from(base != reference[j - 1] || base == b'N');
            current[j] = (previous[j - 1] + mismatch)
                .min(previous[j] + 1)
                .min(current[j - 1] + 1);
            best = best.min(current[j]);
        }
        if best > limit {
            return false;
        }
        if end < m {
            current[end + 1] = outside;
        }
        std::mem::swap(previous, current);
    }
    previous[n - limit..=m].iter().any(|&value| value <= limit)
}
