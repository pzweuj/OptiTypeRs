# OptiTypeRs

OptiTypeRs calls two-field HLA-A, HLA-B and HLA-C genotypes from single-end or paired-end short-read DNA and RNA FASTQ files. It adapts the joint typing model from the original [OptiType repository](https://github.com/FRED-2/Optitype), with an in-process Rust mapper and a statically linked HiGHS solver.

## Downloads

Download and unpack the Linux x86-64 archive (including for WSL) or the Windows x86-64 ZIP from [GitHub Releases](../../releases). Add the unpacked directory to your PATH or run the executable from that directory as `./optitype-rs` on Linux and `.\optitype-rs.exe` in PowerShell.

```sh
optitype-rs --version
optitype-rs run --help
```

## Install or update the reference database

```sh
optitype-rs db build --release latest --out refs
```

Run the same command when you want to check for a newer IPD-IMGT/HLA release. The command prints the versioned database directory. Use that directory for `--db` in subsequent runs. Existing version directories stay intact, and if the current release is already present the command reports that you already have it. The manifest records the release, source commit and checksums. [IPD-IMGT/HLA releases and files](https://github.com/ANHIG/IMGTHLA)

To pin a release or install it from local files, use the commands below. For the offline command, place `hla.dat.zip` and `release_version.txt` from the same official release in `release-files`.

```sh
optitype-rs db build --release 3.66.0 --out refs
optitype-rs db build --release 3.66.0 --source release-files --out refs
```

All HLA-A, HLA-B and HLA-C alleles with complete, unambiguous exon 2 and 3 annotations are considered. RNA uses concatenated exons 2 and 3. DNA uses those exons with adjacent annotated introns and keeps sequences on separate sides of unsequenced regions apart. The manifest reports exclusions and missing intron 2 annotations by locus. Index format 2 also includes homologous class I loci and pseudogenes as competing references. Their observed sequences can contribute even when exon annotations are missing. After joining mates, a fragment whose best hits still include a competitor contributes no typing weight. A mate specific to A, B or C can resolve an otherwise shared hit. Only A, B and C are reported.

To migrate a format 1 database at the same IPD release, run `optitype-rs db build --release 3.66.0 --out refs-v2` and use `--db refs-v2/3.66.0`. The IPD release number and the local index format are independent.

## Type a sample offline

```sh
optitype-rs run --dna -i sample_R1.fastq.gz sample_R2.fastq.gz --db refs/3.66.0 -o results
optitype-rs run --rna -i sample.fastq.gz --db refs/3.66.0 -o results-rna
```

Mate reads are joined by read ID, including conventional `/1` and `/2` suffixes. Synchronized inputs stream in bounded FASTQ batches. For reordered or independently filtered mate files the join retains pending IDs and hits in memory. Paired fragments support the intersection of the two mates' globally best hits. Use `--unpaired-weight` to give single-mapped mates a fractional weight and add `--use-discordant` to include pairs with disjoint best-hit sets.

- Each sample produces one `result.json`. Start with `genotype`, which shows two representative allele names per locus, or `null` for a locus without evidence. `ambiguous_loci` lists loci with unresolved two-field calls. `reads` is the weighted number of covered fragments.
- `alternatives` lists distinct representative genotypes from other evidence models when available. An alternative with the same objective is tied with the main result. `--enumerate N` reports up to N ranked evidence models and up to nine extra best-score ties. The `best_ties_truncated` flag indicates that more tied models exist.
- Add `--details` to include full allele groups, mapping counts and run parameters in the same JSON file.

The mapper uses full-read edit distance with a 3% threshold. The new database, handling of tied hits and candidate set can change support counts and calls compared with historical OptiType. The verification samples are HLA-enriched, so their run times do not represent whole-genome throughput.

## Validation

With the legacy reference and filter, the Rust run retained the upstream two-field calls for both samples. The Rust mapper matched the RazerS3 best-hit set on 3,080/3,084 RNA reads and 88,627/88,630 DNA reads. On the same RazerS3 evidence, the Rust solver reproduced upstream `Reads` and `Objective` for both samples.

With IPD-IMGT/HLA 3.66.0, the format 2 index retained 30,365 eligible A/B/C records. The table lists the primary two-field calls for paired inputs.

| Sample | HLA-A | HLA-B | HLA-C | Comparison |
| --- | --- | --- | --- | --- |
| NA11995 DNA | 01:01 / 01:01 | 08:01 / 57:01 | 06:02 / 07:01 | Matches the reported cell-line genotype. [GM11995 typing](https://www.cellosaurus.org/CVCL_9600) |
| NA12878 DNA | 01:218 / 11:01 | 08:01 / 56:01 | 01:02 / 07:01 | The published A*01:01 / A*11:01 pair occurs in an equally scoring model. [Reference typing, Table 4](https://pmc.ncbi.nlm.nih.gov/articles/PMC11948951/) |
| CRC_81_N RNA | 31:01 / 68:01 | 40:01 / 51:01 | 03:04 / 15:02 | Matches the upstream OptiType run. [OptiType](https://github.com/FRED-2/Optitype) |
| arcasHLA test BAM, converted to paired FASTQ | 01:01 / 03:01 | 07:02 / 39:39 | 01:02 / 08:01 | Comparison with another caller's test sample. [arcasHLA](https://github.com/RabadanLab/arcasHLA) |

The previous A/B/C-only index called a second NA11995 allele, A*68:293Q. All 15 pairs exactly matching its indexed region also matched HLA-H. Including only exon-annotated competitors left a second spurious candidate, A*29:05. Twelve of its 15 marginal pairs matched both ends of an HLA-U record without the required exon annotations. The format 2 run covers 665 weighted fragments and reports an objective of 653.03. The reads still admit equivalent two-field calls, including a tied C*06:229 alternative to C*06:02. A separate synthetic RNA test recovered A*01:525 from this release. These comparisons exercise a specific mapping failure and provide limited evidence about accuracy across samples.

The independent NA12878 check used 5,951 complete pairs from the 1000 Genomes phase 3 mapped exome BAM. Both primary records were recovered from GRCh37 chromosome 6 intervals 29900000–29920000, 31230000–31245000, or 31315000–31330000 (zero-based, half-open), then restored to forward FASTQ orientation. This restricted input can miss divergent or unmapped mates. The expected genotype and the primary call each cover 282 fragments. More best-score ties exist than were enumerated. [NA12878 exome alignment](https://1000genomes.s3.amazonaws.com/phase3/data/NA12878/exome_alignment/NA12878.mapped.ILLUMINA.bwa.CEU.exome.20121211.bam)

## License

Original OptiTypeRs contributions are available under the [MIT License](LICENSE). Upstream OptiType-derived material remains under [BSD-3-Clause](LICENSE-OptiType), whose copyright notice, conditions and disclaimer are included in the source and binary archives. The reference database is distributed separately by [IPD-IMGT/HLA](https://github.com/ANHIG/IMGTHLA).
