use std::collections::{BTreeMap, BTreeSet};

use anyhow::{Result, bail, ensure};
use highs::{HighsModelStatus, RowProblem, Sense};
use serde::Serialize;

use crate::reference::Allele;

pub(crate) type Evidence = BTreeMap<Vec<usize>, f64>;
const RECONSTRUCTION_COST: f64 = 0.01;
pub(crate) const SCORE_TOLERANCE: f64 = 1e-6;

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

pub(crate) struct Results {
    pub solutions: Vec<Solution>,
    pub best_ties_truncated: bool,
}

struct Model {
    groups: Vec<Group>,
    rows: BTreeMap<Vec<usize>, f64>,
    by_locus: BTreeMap<usize, Vec<usize>>,
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
    ) -> Result<Option<Vec<usize>>> {
        let mut problem = RowProblem::new();
        let mut objective = Vec::new();
        let selected: Vec<_> = self
            .groups
            .iter()
            .map(|group| {
                let penalty = if group.reconstructed {
                    RECONSTRUCTION_COST
                } else {
                    0.0
                };
                let column = problem.add_integer_column(-penalty, 0..=1);
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
            let covered = problem.add_integer_column(weight, 0..=1);
            objective.push((covered, weight));
            problem.add_row(
                0..,
                hits.iter()
                    .map(|&id| (selected[id], 1.0))
                    .chain([(covered, -1.0)]),
            );
            if let Some(hetero) = hetero {
                // The negatively weighted variable equals hetero when this row is covered.
                let regularize = problem.add_column(-beta * weight, 0.0..=loci);
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
        let mut model = problem
            .try_optimise(Sense::Maximise)
            .map_err(|err| anyhow::anyhow!("HiGHS model: {err:?}"))?;
        model.set_option("threads", 1);
        model.set_option("random_seed", 0);
        model.set_option("mip_rel_gap", 0.0);
        model.set_option("mip_abs_gap", 0.0);
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
            let two_field = group
                .members
                .iter()
                .map(|&id| alleles[id].two_field.clone())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect();
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
                    .cmp(&right.two_field)
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
    for rank in 0..count.saturating_add(EXTRA_TIES) {
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
                return Ok(Results {
                    solutions,
                    best_ties_truncated: false,
                });
            }
        }
        let Some(selected) = model.solve(beta, &excluded, minimum)? else {
            return Ok(Results {
                solutions,
                best_ties_truncated: false,
            });
        };
        solutions.push(model.report(&selected, alleles, beta));
        excluded.push(selected);
    }
    let best_ties_truncated = model.solve(beta, &excluded, minimum)?.is_some();
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
        let alleles = [allele("A*01:01", 0, false), allele("A*02:01", 0, false)];
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
    }
}
