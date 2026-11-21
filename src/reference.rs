use std::{
    collections::{BTreeMap, BTreeSet},
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail, ensure};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use crate::mapping::{Index, Segment};

const FORMAT: u32 = 2;
const LOCI: &[&str] = &[
    "A", "B", "C", "E", "F", "G", "H", "J", "K", "L", "N", "P", "R", "S", "T", "U", "V", "W", "Y",
];
const RAW: &str = "https://raw.githubusercontent.com/ANHIG/IMGTHLA";

#[derive(Clone, Debug, Serialize, Deserialize)]
pub(crate) struct Allele {
    pub accession: String,
    pub name: String,
    pub two_field: String,
    pub locus: usize,
    pub reconstructed: bool,
}

impl Allele {
    fn new(accession: String, name: String) -> Option<Self> {
        let (gene, fields) = name.split_once('*')?;
        let locus = LOCI.iter().position(|&item| item == gene)?;
        let mut fields = fields.split(':');
        let first = fields.next()?;
        let second = fields.next()?;
        let second_digits = second.trim_end_matches(|c: char| c.is_ascii_uppercase());
        if first.is_empty()
            || second_digits.is_empty()
            || !first
                .bytes()
                .chain(second_digits.bytes())
                .all(|base| base.is_ascii_digit())
        {
            return None;
        }
        let two_field = format!("{gene}*{first}:{second}");
        Some(Self {
            accession,
            name,
            two_field,
            locus,
            reconstructed: false,
        })
    }
}

#[derive(Default, Serialize, Deserialize)]
struct Counts {
    included: usize,
    excluded: usize,
    dna_missing_intron_2: usize,
}

#[derive(Serialize, Deserialize)]
pub(crate) struct Manifest {
    format: u32,
    pub release: String,
    pub commit: Option<String>,
    source_sha256: String,
    sha256: BTreeMap<String, String>,
    loci: BTreeMap<String, Counts>,
}

#[derive(Default)]
struct Entry {
    accession: String,
    name: String,
    bases: Vec<u8>,
    feature: Option<(String, String)>,
    intervals: BTreeMap<(String, usize), (usize, usize)>,
    in_sequence: bool,
}

fn interval(location: &str) -> Option<(usize, usize)> {
    let (start, end) = location.split_once("..").unwrap_or((location, location));
    let start = start.parse::<usize>().ok()?.checked_sub(1)?;
    let end = end.parse().ok()?;
    (start < end).then_some((start, end))
}

impl Entry {
    fn read_line(&mut self, line: &str) {
        if let Some(id) = line.strip_prefix("ID   ") {
            self.accession = id.split(';').next().unwrap_or_default().trim().to_owned();
        } else if let Some(name) = line.strip_prefix("DE   HLA-") {
            self.name = name.split([',', ' ']).next().unwrap_or_default().to_owned();
        } else if let Some(feature) = line.strip_prefix("FT   ") {
            let value = feature.trim();
            if !feature.starts_with(' ') {
                self.feature = value
                    .split_once(char::is_whitespace)
                    .map(|(kind, location)| (kind.to_owned(), location.trim().to_owned()));
            } else if let Some(name) = value
                .strip_prefix("/allele=\"HLA-")
                .and_then(|v| v.strip_suffix('"'))
            {
                self.name = name.to_owned();
            } else if let Some(number) = value
                .strip_prefix("/number=\"")
                .and_then(|v| v.strip_suffix('"'))
                && let Ok(number) = number.parse()
                && let Some((kind, location)) = &self.feature
                && let Some(range) = interval(location)
            {
                self.intervals.insert((kind.clone(), number), range);
            }
        } else if line.starts_with("SQ   ") {
            self.in_sequence = true;
        } else if self.in_sequence && line.starts_with(' ') {
            self.bases
                .extend(line.bytes().filter(u8::is_ascii_alphabetic).map(|base| {
                    match base.to_ascii_uppercase() {
                        base @ (b'A' | b'C' | b'G' | b'T') => base,
                        _ => b'N',
                    }
                }));
        }
    }
}

#[derive(Default)]
struct References {
    alleles: Vec<Allele>,
    dna: BTreeMap<Vec<u8>, Vec<usize>>,
    rna: BTreeMap<Vec<u8>, Vec<usize>>,
    loci: BTreeMap<String, Counts>,
}

fn add_segment(index: &mut BTreeMap<Vec<u8>, Vec<usize>>, bases: &[u8], allele: usize) {
    if !bases.is_empty() {
        let ids = index.entry(bases.to_vec()).or_default();
        if ids.last() != Some(&allele) {
            ids.push(allele);
        }
    }
}

impl References {
    fn add_competitor(&mut self, allele: Allele, entry: &Entry) {
        let counts = self.loci.entry(LOCI[allele.locus].to_owned()).or_default();
        let mut exons: Vec<_> = entry
            .intervals
            .iter()
            .filter_map(|((kind, number), &(start, end))| {
                (kind == "exon" && end <= entry.bases.len()).then_some((*number, start, end))
            })
            .collect();
        exons.sort_unstable();
        let mut boundaries = vec![0, entry.bases.len()];
        for pair in exons.windows(2) {
            let (number, _, end) = pair[0];
            let (next, start, _) = pair[1];
            if end <= start && (next != number + 1 || end == start) {
                boundaries.extend([end, start]);
            }
        }
        boundaries.sort_unstable();
        boundaries.dedup();
        let id = self.alleles.len();
        for window in boundaries.windows(2) {
            for segment in entry.bases[window[0]..window[1]].split(|&b| b == b'N') {
                add_segment(&mut self.dna, segment, id);
                add_segment(&mut self.rna, segment, id);
            }
        }
        // Spliced reads may cross an annotated junction. An unannotated
        // pseudogene still contributes its observed contiguous sequence above.
        let mut transcript = Vec::new();
        let mut previous = None;
        for (number, start, end) in exons {
            if previous.is_some_and(|previous| number != previous + 1) {
                for part in transcript.split(|&b| b == b'N') {
                    add_segment(&mut self.rna, part, id);
                }
                transcript.clear();
            }
            transcript.extend_from_slice(&entry.bases[start..end]);
            previous = Some(number);
        }
        for part in transcript.split(|&b| b == b'N') {
            add_segment(&mut self.rna, part, id);
        }
        counts.included += 1;
        self.alleles.push(allele);
    }

    fn add_entry(&mut self, entry: Entry) {
        let Some(allele) = Allele::new(entry.accession.clone(), entry.name.clone()) else {
            return;
        };
        if allele.locus >= 3 && !entry.bases.is_empty() {
            self.add_competitor(allele, &entry);
            return;
        }
        let counts = self.loci.entry(LOCI[allele.locus].to_owned()).or_default();
        let feature = |kind: &str, number| {
            entry
                .intervals
                .get(&(kind.to_owned(), number))
                .copied()
                .filter(|&(_, end)| end <= entry.bases.len())
        };
        let (Some((start2, end2)), Some((start3, end3))) = (feature("exon", 2), feature("exon", 3))
        else {
            counts.excluded += 1;
            return;
        };
        let exon2 = &entry.bases[start2..end2];
        let exon3 = &entry.bases[start3..end3];
        if end2 > start3 || exon2.contains(&b'N') || exon3.contains(&b'N') {
            counts.excluded += 1;
            return;
        }
        let id = self.alleles.len();
        self.alleles.push(allele);
        counts.included += 1;
        add_segment(&mut self.rna, &[exon2, exon3].concat(), id);

        let start = feature("intron", 1)
            .filter(|&(_, end)| end == start2)
            .map_or(start2, |(start, _)| start);
        let end = feature("intron", 3)
            .filter(|&(start, _)| start == end3)
            .map_or(end3, |(_, end)| end);
        let regions = if end2 < start3 && feature("intron", 2) == Some((end2, start3)) {
            vec![&entry.bases[start..end]]
        } else {
            counts.dna_missing_intron_2 += 1;
            vec![&entry.bases[start..end2], &entry.bases[start3..end]]
        };
        for region in regions {
            for segment in region.split(|&base| base == b'N') {
                add_segment(&mut self.dna, segment, id);
            }
        }
    }
}

fn parse(reader: impl BufRead) -> Result<References> {
    let mut refs = References::default();
    let mut entry = Entry::default();
    for line in reader.lines() {
        let line = line?;
        if line == "//" {
            refs.add_entry(std::mem::take(&mut entry));
        } else {
            entry.read_line(&line);
        }
    }
    ensure!(
        entry.accession.is_empty(),
        "truncated hla.dat entry {}",
        entry.accession
    );
    for locus in ["A", "B", "C"] {
        ensure!(
            refs.loci
                .get(locus)
                .is_some_and(|counts| counts.included > 0),
            "no complete exon 2/3 references for HLA-{locus}"
        );
    }
    let mut seen = BTreeSet::new();
    ensure!(
        refs.alleles
            .iter()
            .all(|allele| !allele.accession.is_empty() && seen.insert(&allele.accession)),
        "missing or duplicate HLA accession"
    );
    Ok(refs)
}

fn checksum(path: &Path) -> Result<String> {
    let mut hasher = Sha256::new();
    std::io::copy(&mut File::open(path)?, &mut hasher)?;
    Ok(format!("{:x}", hasher.finalize()))
}

fn version(text: &str) -> Result<String> {
    let release = text
        .lines()
        .find_map(|line| line.trim().strip_prefix("# version: IPD-IMGT/HLA "))
        .context("missing IPD-IMGT/HLA version in release_version.txt")?;
    release_branch(release)?;
    Ok(release.to_owned())
}

fn release_branch(release: &str) -> Result<String> {
    let fields: Vec<_> = release.split('.').collect();
    ensure!(
        fields.len() == 3
            && fields
                .iter()
                .all(|part| !part.is_empty() && part.bytes().all(|base| base.is_ascii_digit())),
        "expected a release such as 3.66.0"
    );
    Ok(fields.concat())
}

pub(crate) fn build(requested: &str, source: Option<&Path>, output: &Path) -> Result<PathBuf> {
    fs::create_dir_all(output)?;
    let staging = tempfile::tempdir_in(output)?;
    let downloaded = staging.path().join("hla.dat.zip");
    let (archive_path, release, commit) = if let Some(path) = source {
        let release = version(
            &fs::read_to_string(path.join("release_version.txt"))
                .context("read release_version.txt from --source")?,
        )?;
        (path.join("hla.dat.zip"), release, None)
    } else {
        let branch = if requested == "latest" {
            "Latest".to_owned()
        } else {
            release_branch(requested)?
        };
        let metadata: serde_json::Value = ureq::get(format!(
            "https://api.github.com/repos/ANHIG/IMGTHLA/commits/{branch}"
        ))
        .header("User-Agent", "optitype-rs")
        .call()
        .context("resolve IPD-IMGT/HLA release commit")?
        .body_mut()
        .read_json()?;
        let sha = metadata["sha"]
            .as_str()
            .context("missing GitHub commit SHA")?;
        ensure!(
            sha.len() == 40 && sha.bytes().all(|b| b.is_ascii_hexdigit()),
            "invalid GitHub commit SHA"
        );
        let release = version(
            &ureq::get(format!("{RAW}/{sha}/release_version.txt"))
                .call()?
                .body_mut()
                .read_to_string()?,
        )?;
        let mut response = ureq::get(format!("{RAW}/{sha}/hla.dat.zip"))
            .call()
            .context("download hla.dat.zip")?;
        std::io::copy(
            &mut response.body_mut().as_reader(),
            &mut File::create(&downloaded)?,
        )?;
        (downloaded.clone(), release, Some(sha.to_owned()))
    };
    ensure!(
        requested == "latest" || requested == release,
        "requested {requested}, but source is {release}"
    );
    let destination = output.join(&release);
    ensure!(
        !destination.exists(),
        "reference already exists at {}",
        destination.display()
    );
    let source_sha256 = checksum(&archive_path).context("hash hla.dat.zip")?;
    let mut archive = zip::ZipArchive::new(File::open(&archive_path)?)?;
    let members: Vec<_> = archive
        .file_names()
        .filter(|name| name.rsplit('/').next() == Some("hla.dat"))
        .map(str::to_owned)
        .collect();
    ensure!(
        members.len() == 1,
        "expected one hla.dat in release archive"
    );
    let refs = parse(BufReader::new(archive.by_name(&members[0])?))?;
    drop(archive);
    if archive_path == downloaded {
        fs::remove_file(downloaded)?;
    }
    let mut sha256 = BTreeMap::new();
    let path = staging.path().join("alleles.json");
    let mut writer = BufWriter::new(File::create(&path)?);
    serde_json::to_writer(&mut writer, &refs.alleles)?;
    writer.flush()?;
    drop(writer);
    sha256.insert("alleles.json".to_owned(), checksum(&path)?);
    for (name, segments) in [("dna.idx", refs.dna), ("rna.idx", refs.rna)] {
        let index = Index::build(
            segments
                .into_iter()
                .map(|(bases, alleles)| Segment { bases, alleles })
                .collect(),
        )?;
        let path = staging.path().join(name);
        let mut writer = BufWriter::new(File::create(&path)?);
        bincode::serialize_into(&mut writer, &index)?;
        writer.flush()?;
        drop(writer);
        sha256.insert(name.to_owned(), checksum(&path)?);
    }
    let manifest = Manifest {
        format: FORMAT,
        release,
        commit,
        source_sha256,
        sha256,
        loci: refs.loci,
    };
    fs::write(
        staging.path().join("manifest.json"),
        serde_json::to_vec_pretty(&manifest)?,
    )?;
    fs::rename(staging.path(), &destination).context("publish reference database")?;
    Ok(destination)
}

pub(crate) fn load(path: &Path, dna: bool) -> Result<(Manifest, Vec<Allele>, Index)> {
    let manifest: Manifest = serde_json::from_reader(BufReader::new(
        File::open(path.join("manifest.json")).context("open reference manifest")?,
    ))?;
    ensure!(
        manifest.format == FORMAT,
        "unsupported database format {}. Rebuild the release with db build --out refs-v2 to include competing HLA loci.",
        manifest.format
    );
    let name = if dna { "dna.idx" } else { "rna.idx" };
    for file in ["alleles.json", name] {
        if manifest.sha256.get(file) != Some(&checksum(&path.join(file))?) {
            bail!("reference checksum mismatch for {file}");
        }
    }
    let alleles = serde_json::from_reader(BufReader::new(File::open(path.join("alleles.json"))?))?;
    let index = bincode::deserialize_from(BufReader::new(File::open(path.join(name))?))?;
    Ok((manifest, alleles, index))
}
