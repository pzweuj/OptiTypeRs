use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use anyhow::{Context, Result, ensure};
use flate2::bufread::MultiGzDecoder;

pub(crate) struct Reader {
    input: Box<dyn BufRead + Send>,
    line: Vec<u8>,
    record: u64,
    file: std::path::PathBuf,
}

pub(crate) struct Read {
    pub id: Vec<u8>,
    pub bases: Vec<u8>,
}

impl Reader {
    pub fn open(path: &Path) -> Result<Self> {
        let mut file =
            BufReader::new(File::open(path).with_context(|| format!("open {}", path.display()))?);
        let input: Box<dyn BufRead + Send> = if file.fill_buf()?.starts_with(&[0x1f, 0x8b]) {
            Box::new(BufReader::new(MultiGzDecoder::new(file)))
        } else {
            Box::new(file)
        };
        Ok(Self {
            input,
            line: Vec::new(),
            record: 0,
            file: path.to_owned(),
        })
    }

    pub fn next(&mut self) -> Result<Option<Read>> {
        self.record += 1;
        self.read_record()
            .with_context(|| format!("{} FASTQ record {}", self.file.display(), self.record))
    }

    fn read_line(&mut self) -> Result<bool> {
        self.line.clear();
        let read = self.input.read_until(b'\n', &mut self.line)? > 0;
        if self.line.pop_if(|b| *b == b'\n').is_some() {
            self.line.pop_if(|b| *b == b'\r');
        }
        Ok(read)
    }

    fn read_record(&mut self) -> Result<Option<Read>> {
        if !self.read_line()? {
            return Ok(None);
        }
        let id = self
            .line
            .strip_prefix(b"@")
            .context("missing @ header")?
            .split(u8::is_ascii_whitespace)
            .next()
            .unwrap_or_default()
            .to_vec();
        ensure!(!id.is_empty(), "empty read ID");
        let mut bases = Vec::new();
        loop {
            ensure!(self.read_line()?, "missing + separator");
            if self.line.starts_with(b"+") {
                break;
            }
            for &base in &self.line {
                let base = match base.to_ascii_uppercase() {
                    base @ (b'A' | b'C' | b'G' | b'T') => base,
                    b'N' | b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B' | b'D' | b'H' | b'V' => {
                        b'N'
                    }
                    _ => anyhow::bail!("invalid nucleotide byte {base}"),
                };
                bases.push(base);
            }
        }
        ensure!(!bases.is_empty(), "empty read sequence");
        let mut quality_len = 0;
        while quality_len < bases.len() {
            ensure!(self.read_line()?, "truncated quality scores");
            ensure!(
                self.line.iter().all(|byte| (33..=126).contains(byte)),
                "invalid quality score"
            );
            quality_len += self.line.len();
        }
        ensure!(
            quality_len == bases.len(),
            "sequence and quality lengths differ"
        );
        Ok(Some(Read { id, bases }))
    }
}

pub(crate) fn paired_ids(left: &[u8], right: &[u8]) -> bool {
    pair_key(left) == pair_key(right)
}

pub(crate) fn pair_key(id: &[u8]) -> &[u8] {
    id.strip_suffix(b"/1")
        .or_else(|| id.strip_suffix(b"/2"))
        .filter(|id| !id.is_empty())
        .unwrap_or(id)
}
