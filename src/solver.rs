use std::collections::{BTreeMap, BTreeSet};

use anyhow::{Result, bail, ensure};
use highs::{HighsModelStatus, RowProblem, Sense};
use serde::Serialize;

use crate::reference::Allele;

pub(crate) type Evidence = BTreeMap<Vec<usize>, f64>;
const RECONSTRUCTION_COST: f64 = 0.01;
pub(crate) const SCORE_TOLERANCE: f64 = 1e-6;

type AlleleNumber<'a> = (usize, &'a str, usize, &'a str, &'a str);

fn allele_number(name: &str) -> AlleleNumber<'_> {
    let (_, fields) = name.split_once('*').unwrap();
    let (first, second) = fields.split_once(':').unwrap();
    let first = first.trim_start_matches('0');
    let second = second
        .trim_end_matches(|c: char| c.is_ascii_uppercase())
        .trim_start_matches('0');
    (first.len(), first, second.len(), second, name)
}

struct Group {
    locus: usize,
    reconstructed: bool,
    members: Vec<usize>,
}

#[derive(Serialize)]
pub(crate) struct Choice {
    pub two_field: Vec<String>,
    pub alleles: Vec<String>,
}

pub(crate) struct Solution {
    pub loci: [Vec<Choice>; 3],
    pub reads: f64,
    pub objective: f64,
}

impl Solution {
    pub(crate) fn genotype(&self) -> [Option<[&str; 2]>; 3] {
        self.loci.each_ref().map(|groups| {
            let first = groups.first()?.two_field[0].as_str();
            let second = groups
                .get(1)
                .map_or(first, |group| group.two_field[0].as_str());
            Some([first, second])
        })
    }

    fn genotype_key(&self) -> [Option<[AlleleNumber<'_>; 2]>; 3] {
        self.genotype()
            .map(|pair| pair.map(|pair| pair.map(allele_number)))
    }
}

pub(crate) struct Results {
    pub solutions: Vec<Solution>,
    pub best_ties_truncated: bool,
}

struct Model {
    groups: Vec<Group>,
    rows: BTreeMap<Vec<usize>, f64>,
    by_locus: BTreeMap<usize, Vec<usize>>,
}

struct TieBreak<'a> {
    ranks: &'a [usize],
    fixed: &'a [(usize, [usize; 2])],
    locus: usize,
}

impl Model {
    fn new(alleles: &[Allele], evidence: &Evidence) -> Self {
        let mut signatures = vec![Vec::new(); alleles.len()];
        for (row, hits) in evidence.keys().enumerate() {
            for &id in hits {
                signatures[id].push(row);
            }
        }
        let mut grouped = BTreeMap::<_, Vec<_>>::new();
        for (id, signature) in signatures.into_iter().enumerate() {
            if !signature.is_empty() {
                let allele = &alleles[id];
                grouped
                    .entry((allele.locus, allele.reconstructed, signature))
                    .or_default()
                    .push(id);
            }
        }
        let mut groups = Vec::new();
        let mut by_locus = BTreeMap::<_, Vec<_>>::new();
        let mut rows = vec![Vec::new(); evidence.len()];
        for ((locus, reconstructed, signature), members) in grouped {
            let group = groups.len();
            by_locus.entry(locus).or_default().push(group);
            for row in signature {
                rows[row].push(group);
            }
            groups.push(Group {
                locus,
                reconstructed,
                members,
            });
        }
        let mut compact = BTreeMap::new();
        for (hits, weight) in rows.into_iter().zip(evidence.values()) {
            *compact.entry(hits).or_insert(0.0) += weight;
        }
        Self {
            groups,
            rows: compact,
            by_locus,
        }
    }

    fn solve(
        &self,
        beta: f64,
        excluded: &[Vec<usize>],
        minimum: Option<f64>,
        tie_break: Option<TieBreak<'_>>,
    ) -> Result<Option<Vec<usize>>> {
        let mut problem = RowProblem::new();
        let mut objective = Vec::new();
        let primary_cost = |cost| if tie_break.is_some() { 0.0 } else { cost };
        let selected: Vec<_> = self
            .groups
            .iter()
            .map(|group| {
                let penalty = if group.reconstructed {
                    RECONSTRUCTION_COST
                } else {
                    0.0
                };
                let column = problem.add_integer_column(primary_cost(-penalty), 0..=1);
                if penalty > 0.0 {
                    objective.push((column, -penalty));
                }
                column
            })
            .collect();
        let loci = self.by_locus.len() as f64;
        for locus in self.by_locus.values() {
            problem.add_row(1..=2, locus.iter().map(|&id| (selected[id], 1.0)));
        }
        let hetero = (beta > 0.0).then(|| {
            let h = problem.add_column(0.0, 0.0..=loci);
            problem.add_row(
                ..=loci,
                selected.iter().map(|&col| (col, 1.0)).chain([(h, -1.0)]),
            );
            h
        });
        for (hits, &weight) in &self.rows {
            let covered = problem.add_integer_column(primary_cost(weight), 0..=1);
            objective.push((covered, weight));
            problem.add_row(
                0..,
                hits.iter()
                    .map(|&id| (selected[id], 1.0))
                    .chain([(covered, -1.0)]),
            );
            if let Some(hetero) = hetero {
                // The negatively weighted variable equals hetero when this row is covered.
                let regularize = problem.add_column(primary_cost(-beta * weight), 0.0..=loci);
                objective.push((regularize, -beta * weight));
                problem.add_row(
                    -loci..,
                    [(regularize, 1.0), (hetero, -1.0), (covered, -loci)],
                );
            }
        }
        if let Some(minimum) = minimum {
            problem.add_row(minimum.., objective);
        }
        for former in excluded {
            problem.add_row(
                (1.0 - former.len() as f64)..,
                selected.iter().enumerate().map(|(id, &column)| {
                    (
                        column,
                        if former.binary_search(&id).is_ok() {
                            -1.0
                        } else {
                            1.0
                        },
                    )
                }),
            );
        }
        if let Some(TieBreak {
            ranks,
            fixed,
            locus,
        }) = tie_break
        {
            for &(locus, [first, second]) in fixed {
                let ids = &self.by_locus[&locus];
                problem.add_row(
                    ..=0,
                    ids.iter()
                        .filter(|&&id| !(first..=second).contains(&ranks[id]))
                        .map(|&id| (selected[id], 1.0)),
                );
                for rank in BTreeSet::from([first, second]) {
                    problem.add_row(
                        1..,
                        ids.iter()
                            .filter(|&&id| ranks[id] == rank)
                            .map(|&id| (selected[id], 1.0)),
                    );
                }
            }
            let ids = &self.by_locus[&locus];
            let max_rank = ids.iter().map(|&id| ranks[id]).max().unwrap() as f64;
            // Pick the first allele from the selected groups. Minimizing the
            // second rank takes the other group's rank, or repeats the first
            // rank for a homozygous call.
            let second = problem.add_column(-1.0, 0.0..=max_rank);
            let first: Vec<_> = ids
                .iter()
                .map(|&id| {
                    let column =
                        problem.add_integer_column(-(max_rank + 1.0) * ranks[id] as f64, 0..=1);
                    problem.add_row(..=0, [(column, 1.0), (selected[id], -1.0)]);
                    (column, ranks[id] as f64)
                })
                .collect();
            problem.add_row(1..=1, first.iter().map(|&(column, _)| (column, 1.0)));
            problem.add_row(..=0.0, first.iter().copied().chain([(second, -1.0)]));
            problem.add_row(
                ..=0.0,
                ids.iter()
                    .map(|&id| (selected[id], ranks[id] as f64))
                    .chain(first.iter().map(|&(column, rank)| (column, -rank)))
                    .chain([(second, -1.0)]),
            );
        }
        let mut model = problem
            .try_optimise(Sense::Maximise)
            .map_err(|err| anyhow::anyhow!("HiGHS model: {err:?}"))?;
        model.set_option("threads", 1);
        model.set_option("random_seed", 0);
        model.set_option("mip_rel_gap", 0.0);
        model.set_option("mip_abs_gap", 0.0);
        if minimum.is_some() {
            model.set_option("mip_feasibility_tolerance", 1e-9);
        }
        let solved = model
            .try_solve()
            .map_err(|err| anyhow::anyhow!("HiGHS solve: {err:?}"))?;
        match solved.status() {
            HighsModelStatus::Optimal => Ok(Some(
                solved.get_solution().columns()[..selected.len()]
                    .iter()
                    .enumerate()
                    .filter_map(|(id, &value)| (value > 0.5).then_some(id))
                    .collect(),
            )),
            HighsModelStatus::Infeasible if !excluded.is_empty() => Ok(None),
            status => bail!("HiGHS did not find an optimal genotype ({status:?})"),
        }
    }

    fn preferred_tie(&self, alleles: &[Allele], beta: f64, minimum: f64) -> Result<Vec<usize>> {
        let representatives: Vec<_> = self
            .groups
            .iter()
            .map(|group| {
                group
                    .members
                    .iter()
                    .map(|&id| allele_number(&alleles[id].two_field))
                    .min()
                    .unwrap()
            })
            .collect();
        let mut ordered = representatives.clone();
        ordered.sort_unstable();
        ordered.dedup();
        let ranks: Vec<_> = representatives
            .iter()
            .map(|name| ordered.binary_search(name).unwrap())
            .collect();
        let mut fixed = Vec::new();
        let mut preferred = Vec::new();
        for &locus in self.by_locus.keys().filter(|&&locus| locus < 3) {
            preferred = self
                .solve(
                    beta,
                    &[],
                    Some(minimum),
                    Some(TieBreak {
                        ranks: &ranks,
                        fixed: &fixed,
                        locus,
                    }),
                )?
                .unwrap();
            let mut pair = preferred
                .iter()
                .filter(|&&id| self.groups[id].locus == locus)
                .map(|&id| ranks[id]);
            let first = pair.next().unwrap();
            let mut pair = [first, pair.next().unwrap_or(first)];
            pair.sort_unstable();
            fixed.push((locus, pair));
        }
        Ok(preferred)
    }

    fn report(&self, selected: &[usize], alleles: &[Allele], beta: f64) -> Solution {
        let mut loci: [Vec<Choice>; 3] = Default::default();
        let reads = self
            .rows
            .iter()
            .filter(|(hits, _)| selected.iter().any(|id| hits.binary_search(id).is_ok()))
            .map(|(_, weight)| weight)
            .sum::<f64>();
        let mut reconstructed = 0;
        for &id in selected {
            let group = &self.groups[id];
            reconstructed += usize::from(group.reconstructed);
            let Some(locus) = loci.get_mut(group.locus) else {
                continue;
            };
            let mut two_field: Vec<_> = group
                .members
                .iter()
                .map(|&id| alleles[id].two_field.clone())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect();
            two_field.sort_by(|a, b| allele_number(a).cmp(&allele_number(b)));
            let names = group
                .members
                .iter()
                .map(|&id| alleles[id].name.clone())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect();
            locus.push(Choice {
                two_field,
                alleles: names,
            });
        }
        let hetero = selected.len() - self.by_locus.len();
        for locus in &mut loci {
            locus.sort_by(|left, right| {
                left.two_field
                    .iter()
                    .map(|name| allele_number(name))
                    .cmp(right.two_field.iter().map(|name| allele_number(name)))
                    .then(left.alleles.cmp(&right.alleles))
            });
        }
        let objective =
            reads * (1.0 - beta * hetero as f64) - RECONSTRUCTION_COST * reconstructed as f64;
        Solution {
            loci,
            reads,
            objective,
        }
    }
}

pub(crate) fn solve(
    alleles: &[Allele],
    evidence: &Evidence,
    beta: f64,
    count: usize,
) -> Result<Results> {
    ensure!(
        !evidence.is_empty(),
        "no HLA read evidence; check the input, reference and paired-end settings"
    );
    let model = Model::new(alleles, evidence);
    let mut excluded = Vec::new();
    let mut solutions = Vec::new();
    const EXTRA_TIES: usize = 9;
    let mut minimum = None;
    let limit = count.saturating_add(EXTRA_TIES);
    for rank in 0..limit {
        if rank == count {
            minimum = Some(
                solutions
                    .first()
                    .map(|best: &Solution| best.objective)
                    .unwrap()
                    - SCORE_TOLERANCE,
            );
            if solutions
                .last()
                .is_some_and(|last| last.objective < minimum.unwrap())
            {
                break;
            }
        }
        let Some(selected) = model.solve(beta, &excluded, minimum, None)? else {
            break;
        };
        solutions.push(model.report(&selected, alleles, beta));
        excluded.push(selected);
    }
    let best_ties_truncated =
        solutions.len() == limit && model.solve(beta, &excluded, minimum, None)?.is_some();
    if best_ties_truncated {
        // The enumeration cap must not affect which genotype is primary.
        let preferred = model.preferred_tie(alleles, beta, minimum.unwrap())?;
        if !excluded.contains(&preferred) {
            solutions.pop();
            solutions.push(model.report(&preferred, alleles, beta));
        }
    }
    solutions.sort_by(|a, b| {
        b.objective
            .total_cmp(&a.objective)
            .then_with(|| a.genotype_key().cmp(&b.genotype_key()))
    });
    if let Some(best) = solutions.first() {
        let cutoff = best.objective - SCORE_TOLERANCE;
        let end = solutions.partition_point(|solution| solution.objective >= cutoff);
        solutions[..end].sort_by(|a, b| a.genotype_key().cmp(&b.genotype_key()));
    }
    Ok(Results {
        solutions,
        best_ties_truncated,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn allele(name: &str, locus: usize, reconstructed: bool) -> Allele {
        Allele {
            accession: name.to_owned(),
            name: name.to_owned(),
            two_field: name.to_owned(),
            locus,
            reconstructed,
        }
    }

    #[test]
    fn jointly_maximizes_weighted_coverage_and_reports_equivalent_calls() {
        let alleles = [
            allele("A*01:01", 0, false),
            allele("A*02:01", 0, false),
            allele("A*03:01", 0, false),
            allele("B*01:01", 1, false),
            allele("B*02:01", 1, true),
        ];
        let evidence = Evidence::from([
            (vec![0, 2, 3, 4], 8.0),
            (vec![0, 2], 2.5),
            (vec![1], 0.5),
            (vec![3, 4], 2.0),
        ]);
        let beta = 0.02;
        let best = solve(&alleles, &evidence, beta, 1)
            .unwrap()
            .solutions
            .remove(0);
        let exhaustive = (1_u32..1 << alleles.len())
            .filter_map(|mask| {
                let selected: Vec<_> = (0..alleles.len())
                    .filter(|&id| mask & (1 << id) != 0)
                    .collect();
                let a = selected
                    .iter()
                    .filter(|&&id| alleles[id].locus == 0)
                    .count();
                let b = selected.len() - a;
                if !(1..=2).contains(&a) || !(1..=2).contains(&b) {
                    return None;
                }
                let covered: f64 = evidence
                    .iter()
                    .filter(|(hits, _)| selected.iter().any(|id| hits.contains(id)))
                    .map(|(_, weight)| weight)
                    .sum();
                let hetero = usize::from(a == 2) + usize::from(b == 2);
                let penalty = selected
                    .iter()
                    .filter(|&&id| alleles[id].reconstructed)
                    .count() as f64
                    * 0.01;
                Some(covered * (1.0 - beta * hetero as f64) - penalty)
            })
            .fold(f64::NEG_INFINITY, f64::max);
        assert!((best.objective - exhaustive).abs() < 1e-9);
        assert_eq!(best.reads, 13.0);
        assert_eq!(best.loci[0][0].two_field, ["A*01:01", "A*03:01"]);
        assert_eq!(best.loci[0][1].two_field, ["A*02:01"]);
        assert_eq!(best.loci[1][0].two_field, ["B*01:01"]);
        assert!(best.loci[2].is_empty());
    }

    #[test]
    fn reports_best_score_ties_separately() {
        let alleles = [allele("A*01:218", 0, false), allele("A*01:01", 0, false)];
        let evidence = Evidence::from([(vec![0], 1.0), (vec![1], 1.0)]);
        let results = solve(&alleles, &evidence, 0.0, 3).unwrap();
        assert_eq!(
            results
                .solutions
                .iter()
                .map(|sol| sol.reads)
                .collect::<Vec<_>>(),
            [2.0, 1.0, 1.0]
        );
        assert!(!results.best_ties_truncated);
        let ties = Evidence::from([(vec![0], 1.0), (vec![1], 1.0), (vec![0, 1], 100.0)]);
        let results = solve(&alleles, &ties, 0.05, 1).unwrap();
        assert_eq!(results.solutions.len(), 2);
        assert!(results.solutions.iter().all(|sol| sol.objective == 101.0));
        assert_eq!(results.solutions[0].genotype()[0], Some(["A*01:01"; 2]));
        assert_eq!(results.solutions[1].genotype()[0], Some(["A*01:218"; 2]));
    }

    #[test]
    fn chooses_lowest_numeric_genotype_when_best_ties_are_truncated() {
        let mut alleles: Vec<_> = (99..=112)
            .rev()
            .map(|number| allele(&format!("A*01:{number}"), 0, false))
            .collect();
        let partner = alleles.len();
        let mut evidence =
            Evidence::from([((0..partner).collect(), 100.0), (vec![partner], 100.0)]);
        evidence.extend((0..partner).map(|id| (vec![id], 1.0)));
        alleles.push(allele("A*11:01", 0, false));

        let results = solve(&alleles, &evidence, 0.009, 1).unwrap();
        assert!(results.best_ties_truncated);
        assert_eq!(results.solutions.len(), 10);
        assert_eq!(
            results.solutions[0].genotype()[0],
            Some(["A*01:99", "A*11:01"])
        );
        assert!(results.solutions.iter().all(|sol| sol.reads == 201.0));
    }
}
