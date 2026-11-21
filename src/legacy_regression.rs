use std::{
    collections::{BTreeMap, BTreeSet},
    fs::{self, File},
    io::{BufRead, BufReader},
    path::Path,
};

use rayon::prelude::*;
use serde::Deserialize;

use crate::{
    Run, collect_evidence, fastq,
    mapping::{Index, Segment},
    reference::Allele,
    solver::{self, Evidence},
};

#[derive(Deserialize)]
struct LegacyAllele {
    name: String,
    two_field: String,
    allowed: bool,
}

fn legacy_index(
    fasta: &Path,
    metadata: &BTreeMap<String, LegacyAllele>,
) -> (Index, Vec<Allele>, Vec<bool>) {
    let mut entries = Vec::new();
    for line in BufReader::new(File::open(fasta).unwrap()).lines() {
        let line = line.unwrap();
        if let Some(header) = line.strip_prefix('>') {
            entries.push((
                header.split_whitespace().next().unwrap().to_owned(),
                Vec::new(),
            ));
        } else {
            let (_, bases): &mut (String, Vec<u8>) = entries.last_mut().unwrap();
            bases.extend(line.bytes().map(|base| match base.to_ascii_uppercase() {
                exact @ (b'A' | b'C' | b'G' | b'T') => exact,
                _ => b'N',
            }));
        }
    }
    let mut alleles = Vec::new();
    let mut allowed = Vec::new();
    let mut segments = BTreeMap::<_, Vec<_>>::new();
    for (accession, bases) in entries {
        let parts = accession.split_once('_');
        let reconstructed = parts.is_some();
        let meta = &metadata[parts.map_or(accession.as_str(), |(first, _)| first)];
        let locus = b"ABCHJG"
            .iter()
            .position(|&locus| meta.two_field.as_bytes()[0] == locus);
        allowed.push(meta.allowed && locus.is_some());
        segments.entry(bases).or_default().push(alleles.len());
        alleles.push(Allele {
            accession,
            name: meta.name.clone(),
            two_field: meta.two_field.clone(),
            locus: locus.unwrap_or(usize::MAX),
            reconstructed,
        });
    }
    let index = Index::build(
        segments
            .into_iter()
            .map(|(bases, alleles)| Segment { bases, alleles })
            .collect(),
    )
    .unwrap();
    (index, alleles, allowed)
}

#[test]
#[ignore = "requires OptiType-master and comparator artifacts for the two upstream samples"]
fn upstream_legacy_reference_comparison() {
    let upstream = std::env::var("OPTITYPE_UPSTREAM").unwrap();
    let baseline = std::env::var("OPTITYPE_BASELINE").unwrap();
    let upstream = Path::new(&upstream);
    let baseline = Path::new(&baseline);
    let metadata = serde_json::from_reader(BufReader::new(
        File::open(baseline.join("legacy.json")).unwrap(),
    ))
    .unwrap();
    for (mode, fixture) in [("rna", "rna"), ("dna", "exome")] {
        let baseline = baseline.join(mode);
        let (index, alleles, allowed) = legacy_index(
            &upstream.join(format!("data/hla_reference_{mode}.fasta")),
            &metadata,
        );
        let mut input: Vec<_> = fs::read_dir(upstream.join("test").join(fixture))
            .unwrap()
            .map(|entry| entry.unwrap().path())
            .filter(|path| path.extension().is_some_and(|ext| ext == "fastq"))
            .collect();
        input.sort();
        let args = Run {
            dna: mode == "dna",
            rna: mode == "rna",
            input,
            db: Default::default(),
            out: Default::default(),
            beta: 0.009,
            unpaired_weight: 0.0,
            use_discordant: false,
            enumerate: 1,
            details: false,
            threads: Some(8),
        };
        let mut upstream_mates = Vec::new();
        for (mate, input) in args.input.iter().enumerate() {
            let mut expected: BTreeMap<_, Vec<_>> = BTreeMap::new();
            for line in
                BufReader::new(File::open(baseline.join(format!("hits_{}.tsv", mate + 1))).unwrap())
                    .lines()
            {
                let line = line.unwrap();
                let (read, hits) = line.split_once('\t').unwrap();
                expected.insert(
                    read.as_bytes().to_vec(),
                    hits.split(',').map(str::to_owned).collect(),
                );
            }
            let mut reads = fastq::Reader::open(input).unwrap();
            let mut total = 0;
            let mut identical = 0;
            loop {
                let chunk: Vec<_> = (0..512).map_while(|_| reads.next().unwrap()).collect();
                if chunk.is_empty() {
                    break;
                }
                total += chunk.len();
                identical += chunk
                    .par_iter()
                    .filter(|read| {
                        let mut ids: Vec<_> = index
                            .best_hits(&read.bases)
                            .into_iter()
                            .map(|id| &alleles[id].accession)
                            .collect();
                        ids.sort_unstable();
                        expected
                            .get(&read.id)
                            .map_or(ids.is_empty(), |expected| ids.into_iter().eq(expected))
                    })
                    .count();
            }
            println!(
                "{mode} mate {}: identical best-hit sets {identical}/{total}",
                mate + 1
            );
            upstream_mates.push(expected);
        }
        let (evidence, stats) = collect_evidence(&args, &index).unwrap();
        let mut filtered = Evidence::new();
        for (ids, weight) in evidence {
            let hits: Vec<_> = ids.into_iter().filter(|&id| allowed[id]).collect();
            if !hits.is_empty() {
                *filtered.entry(hits).or_default() += weight;
            }
        }
        let results = solver::solve(&alleles, &filtered, args.beta, 1).unwrap();
        let expected = fs::read_to_string(baseline.join("upstream_result.tsv")).unwrap();
        let row: Vec<_> = expected.lines().nth(1).unwrap().split('\t').collect();
        let best = &results.solutions[0];
        println!(
            "{mode}: {} Rust proper pairs; Rust reads={} objective={}; upstream reads={} objective={}",
            stats.paired, best.reads, best.objective, row[7], row[8]
        );
        for (locus, expected) in best.loci.iter().zip(row[1..7].chunks_exact(2)) {
            let expected: BTreeSet<_> = expected.iter().copied().collect();
            let reported: BTreeSet<_> = locus
                .iter()
                .flat_map(|choice| &choice.two_field)
                .map(String::as_str)
                .collect();
            assert!(
                expected.is_subset(&reported),
                "upstream {expected:?} absent from Rust candidates {reported:?}"
            );
        }
        let allele_ids: BTreeMap<_, _> = alleles
            .iter()
            .enumerate()
            .filter(|(id, _)| allowed[*id])
            .map(|(id, allele)| (&allele.accession, id))
            .collect();
        let mut upstream_evidence = Evidence::new();
        for (read, hits) in &upstream_mates[0] {
            let Some(other) = upstream_mates[1].get(read) else {
                continue;
            };
            let mut hits: Vec<_> = hits
                .iter()
                .filter(|id| other.binary_search(id).is_ok())
                .filter_map(|id| allele_ids.get(id).copied())
                .collect();
            if !hits.is_empty() {
                hits.sort_unstable();
                *upstream_evidence.entry(hits).or_default() += 1.0;
            }
        }
        let model = solver::solve(&alleles, &upstream_evidence, args.beta, 1).unwrap();
        let best = &model.solutions[0];
        println!(
            "{mode}: Rust solver on upstream RazerS3 evidence reads={} objective={}",
            best.reads, best.objective
        );
        assert!((best.reads - row[7].parse::<f64>().unwrap()).abs() < 1e-6);
        assert!((best.objective - row[8].parse::<f64>().unwrap()).abs() < 1e-6);
    }
}
