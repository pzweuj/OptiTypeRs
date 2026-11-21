mod fastq;
#[cfg(test)]
mod legacy_regression;
mod mapping;
mod reference;
mod solver;

use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    io::Write,
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, ensure};
use clap::{Args, Parser, Subcommand};
use rayon::prelude::*;
use serde_json::json;

#[derive(Parser)]
#[command(version = RELEASE_VERSION, about = "HLA-A, HLA-B and HLA-C typing from short-read FASTQ")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

const RELEASE_VERSION: &str = match option_env!("OPTITYPE_RELEASE_VERSION") {
    Some(version) => version,
    None => env!("CARGO_PKG_VERSION"),
};

#[derive(Subcommand)]
enum Command {
    /// Prepare a versioned IPD-IMGT/HLA reference.
    Db {
        #[command(subcommand)]
        action: DbCommand,
    },
    /// Type one sample against a local reference.
    Run(Run),
}

#[derive(Subcommand)]
enum DbCommand {
    Build {
        #[arg(long, default_value = "latest")]
        release: String,
        /// Directory with hla.dat.zip and release_version.txt for offline builds.
        #[arg(long)]
        source: Option<PathBuf>,
        #[arg(short, long)]
        out: PathBuf,
    },
}

#[derive(Args)]
struct Run {
    #[arg(long, required_unless_present = "rna", conflicts_with = "rna")]
    dna: bool,
    #[arg(long)]
    rna: bool,
    /// One FASTQ, or two mate FASTQs, optionally gzip-compressed.
    #[arg(short, long, num_args = 1..=2, required = true)]
    input: Vec<PathBuf>,
    #[arg(long)]
    db: PathBuf,
    #[arg(short, long)]
    out: PathBuf,
    #[arg(long, default_value_t = 0.009)]
    beta: f64,
    #[arg(long, default_value_t = 0.0)]
    unpaired_weight: f64,
    #[arg(long)]
    use_discordant: bool,
    /// Maximum ranked models, with extra best-score ties reported when found.
    #[arg(long, default_value_t = 1)]
    enumerate: usize,
    /// Include complete allele groups and run metadata in result.json.
    #[arg(long)]
    details: bool,
    /// FASTQ mapping threads, defaulting to the number of available CPUs.
    #[arg(long)]
    threads: Option<usize>,
}

#[derive(Default, serde::Serialize)]
struct Stats {
    fragments: u64,
    paired: u64,
    single_end: u64,
    unpaired: u64,
    discordant: u64,
    unmapped: u64,
}

enum Hit {
    Paired(Vec<usize>),
    Single(Vec<usize>),
    Unpaired(Vec<usize>),
    Discordant(Vec<usize>),
    Unmapped,
}

fn fragment_hits(left: Vec<usize>, right: Option<Vec<usize>>) -> Hit {
    let Some(right) = right else {
        return if left.is_empty() {
            Hit::Unmapped
        } else {
            Hit::Single(left)
        };
    };
    match (left.is_empty(), right.is_empty()) {
        (true, true) => Hit::Unmapped,
        (true, false) => Hit::Unpaired(right),
        (false, true) => Hit::Unpaired(left),
        (false, false) => {
            let overlap: Vec<_> = left
                .iter()
                .copied()
                .filter(|id| right.binary_search(id).is_ok())
                .collect();
            if overlap.is_empty() {
                let mut combined = left;
                combined.extend(right);
                combined.sort_unstable();
                Hit::Discordant(combined)
            } else {
                Hit::Paired(overlap)
            }
        }
    }
}

impl Stats {
    fn observe(&mut self, hit: Hit, args: &Run, evidence: &mut solver::Evidence) {
        self.fragments += 1;
        let (ids, weight) = match hit {
            Hit::Single(ids) => {
                self.single_end += 1;
                (ids, 1.0)
            }
            Hit::Paired(ids) => {
                self.paired += 1;
                (ids, 1.0)
            }
            Hit::Unpaired(ids) => {
                self.unpaired += 1;
                (ids, args.unpaired_weight)
            }
            Hit::Discordant(ids) => {
                self.discordant += 1;
                (
                    ids,
                    if args.use_discordant {
                        args.unpaired_weight
                    } else {
                        0.0
                    },
                )
            }
            Hit::Unmapped => {
                self.unmapped += 1;
                return;
            }
        };
        if weight > 0.0 {
            *evidence.entry(ids).or_default() += weight;
        }
    }
}

fn collect_evidence(args: &Run, index: &mapping::Index) -> Result<(solver::Evidence, Stats)> {
    let mut left = fastq::Reader::open(&args.input[0])?;
    let mut right = args
        .input
        .get(1)
        .map(|path| fastq::Reader::open(path))
        .transpose()?;
    let mut stats = Stats::default();
    let mut evidence = solver::Evidence::new();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads.unwrap_or(0))
        .build()?;
    let mut chunk = Vec::with_capacity(512);
    let mut ordered = true;
    let mut pending: [BTreeMap<Vec<u8>, Vec<usize>>; 2] = Default::default();
    let mut seen: [BTreeSet<Vec<u8>>; 2] = Default::default();
    loop {
        for _ in 0..chunk.capacity() {
            let a = left.next()?;
            let b = right
                .as_mut()
                .map(fastq::Reader::next)
                .transpose()?
                .flatten();
            if a.is_none() && b.is_none() {
                break;
            }
            chunk.push([a, b]);
        }
        if chunk.is_empty() {
            break;
        }
        let hits: Vec<_> = pool.install(|| {
            chunk
                .par_iter()
                .map(|mates| {
                    mates
                        .each_ref()
                        .map(|mate| mate.as_ref().map(|read| index.best_hits(&read.bases)))
                })
                .collect()
        });
        for (mates, [first, second]) in chunk.drain(..).zip(hits) {
            if right.is_none() {
                stats.observe(
                    fragment_hits(first.unwrap_or_default(), None),
                    args,
                    &mut evidence,
                );
                continue;
            }
            if ordered
                && mates[0]
                    .as_ref()
                    .zip(mates[1].as_ref())
                    .is_some_and(|(first, second)| fastq::paired_ids(&first.id, &second.id))
            {
                stats.observe(
                    fragment_hits(first.unwrap_or_default(), second),
                    args,
                    &mut evidence,
                );
                continue;
            }
            ordered = false;
            for (side, (read, hits)) in mates.into_iter().zip([first, second]).enumerate() {
                let Some(read) = read else { continue };
                let id = fastq::pair_key(&read.id).to_owned();
                ensure!(
                    seen[side].insert(id.clone()),
                    "duplicate mate {} FASTQ ID {}",
                    side + 1,
                    String::from_utf8_lossy(&id)
                );
                let hits = hits.unwrap_or_default();
                if let Some(other) = pending[1 - side].remove(&id) {
                    stats.observe(fragment_hits(hits, Some(other)), args, &mut evidence);
                } else {
                    pending[side].insert(id, hits);
                }
            }
        }
    }
    for mates in pending {
        for (_, hits) in mates {
            stats.observe(fragment_hits(hits, Some(Vec::new())), args, &mut evidence);
        }
    }
    Ok((evidence, stats))
}

type Genotype<'a> = BTreeMap<&'static str, Option<[&'a str; 2]>>;

#[derive(serde::Serialize)]
struct ModelSummary<'a> {
    genotype: Genotype<'a>,
    reads: f64,
    objective: f64,
}

impl<'a> From<&'a solver::Solution> for ModelSummary<'a> {
    fn from(solution: &'a solver::Solution) -> Self {
        let genotype = ["A", "B", "C"]
            .into_iter()
            .zip(&solution.loci)
            .map(|(locus, groups)| {
                let call = groups.first().map(|first| {
                    let first = first.two_field[0].as_str();
                    let second = groups
                        .get(1)
                        .map_or(first, |group| group.two_field[0].as_str());
                    let mut pair = [first, second];
                    pair.sort_unstable();
                    pair
                });
                (locus, call)
            })
            .collect();
        Self {
            genotype,
            reads: solution.reads,
            objective: solution.objective,
        }
    }
}

#[derive(serde::Serialize)]
struct Report<'a> {
    #[serde(flatten)]
    primary: ModelSummary<'a>,
    reference: &'a str,
    ambiguous_loci: Vec<&'static str>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    alternatives: Vec<ModelSummary<'a>>,
    best_ties_truncated: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    details: Option<serde_json::Value>,
}

fn write_result(out: &Path, result: &Report<'_>) -> Result<()> {
    let mut file = tempfile::NamedTempFile::new_in(out).context("create temporary result file")?;
    file.write_all(&serde_json::to_vec_pretty(result)?)
        .context("write temporary result file")?;
    file.into_temp_path()
        .persist_noclobber(out.join("result.json"))
        .context("publish result.json")?;
    Ok(())
}

fn run(args: Run) -> Result<()> {
    ensure!(
        (0.0..0.1).contains(&args.beta),
        "--beta must be in [0, 0.1)"
    );
    ensure!(
        (0.0..=1.0).contains(&args.unpaired_weight),
        "--unpaired-weight must be in [0, 1]"
    );
    ensure!(
        args.enumerate > 0 && args.threads != Some(0),
        "--enumerate and --threads must be positive"
    );
    let out = std::path::absolute(&args.out);
    // Linux can still resolve paths against the process's working directory
    // when getcwd fails on a mounted filesystem. Tempfile needs an absolute path.
    #[cfg(target_os = "linux")]
    let out = out.or_else(|error| {
        if error.kind() == std::io::ErrorKind::NotFound && args.out.is_relative() {
            Ok(Path::new("/proc/self/cwd").join(&args.out))
        } else {
            Err(error)
        }
    });
    let out = out.with_context(|| format!("resolve output directory {}", args.out.display()))?;
    fs::create_dir_all(&out)
        .with_context(|| format!("create output directory {}", args.out.display()))?;
    ensure!(
        !out.join("result.json").try_exists()?,
        "{} already exists. Choose a new output directory with -o.",
        args.out.join("result.json").display()
    );
    let (manifest, alleles, index) = reference::load(&args.db, args.dna)?;
    let (mut evidence, stats) = collect_evidence(&args, &index)?;
    drop(index);
    // Fragments that a competing locus explains equally well cannot distinguish
    // a target allele. Apply this after pairing so that an informative mate
    // can rescue a fragment whose other end maps to multiple loci.
    let mut competing_weight = 0.0;
    evidence.retain(|ids, weight| {
        let informative = ids.iter().all(|&id| alleles[id].locus < 3);
        if !informative {
            competing_weight += *weight;
        }
        informative
    });
    eprintln!(
        "Processed {} fragments, {} distinct HLA evidence patterns",
        stats.fragments,
        evidence.len()
    );
    let results = solver::solve(&alleles, &evidence, args.beta, args.enumerate)?;
    let (best, other) = results
        .solutions
        .split_first()
        .context("solver returned no model")?;
    let primary = ModelSummary::from(best);
    let mut seen = BTreeSet::from([primary.genotype.clone()]);
    let alternatives: Vec<_> = other
        .iter()
        .map(ModelSummary::from)
        .filter(|model| seen.insert(model.genotype.clone()))
        .collect();
    let ambiguous_loci = ["A", "B", "C"]
        .into_iter()
        .enumerate()
        .filter_map(|(i, locus)| {
            let ambiguous = best.loci[i].iter().any(|group| group.two_field.len() > 1)
                || alternatives.iter().any(|model| {
                    (model.objective - primary.objective).abs() <= solver::SCORE_TOLERANCE
                        && model.genotype[locus] != primary.genotype[locus]
                })
                || other.iter().any(|model| {
                    (model.objective - primary.objective).abs() <= solver::SCORE_TOLERANCE
                        && model.loci[i].iter().any(|group| group.two_field.len() > 1)
                });
            ambiguous.then_some(locus)
        })
        .collect();
    let details =
        args.details.then(|| {
            let models: Vec<_> = results.solutions.iter().map(|solution| json!({
            "loci": { "A": solution.loci[0], "B": solution.loci[1], "C": solution.loci[2] },
            "reads": solution.reads,
            "objective": solution.objective,
        })).collect();
            json!({
                "commit": manifest.commit,
                "mode": if args.dna { "dna" } else { "rna" },
                "inputs": args.input,
                "beta": args.beta,
                "unpaired_weight": args.unpaired_weight,
                "use_discordant": args.use_discordant,
                "enumerate": args.enumerate,
                "fragment_counts": stats,
            "competing_locus_weight": competing_weight,
                "models": models,
                })
        });
    let result = Report {
        primary,
        reference: &manifest.release,
        ambiguous_loci,
        alternatives,
        best_ties_truncated: results.best_ties_truncated,
        details,
    };
    write_result(&out, &result).context("write typing results")
}

fn main() -> Result<()> {
    match Cli::parse().command {
        Command::Db {
            action:
                DbCommand::Build {
                    release,
                    source,
                    out,
                },
        } => {
            println!(
                "{}",
                reference::build(&release, source.as_deref(), &out)?.display()
            );
            Ok(())
        }
        Command::Run(args) => run(args),
    }
}
