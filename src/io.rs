use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

fn open(path: &Path) -> Result<BufReader<File>> {
    let f = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    Ok(BufReader::new(f))
}

fn parse_f64(s: &str) -> Result<f64> {
    let t = s.trim();
    t.parse::<f64>()
        .map_err(|_| RsomicsError::InvalidInput(format!("non-numeric value '{t}'")))
}

fn split<'a>(line: &'a str, comma: bool) -> std::str::Split<'a, char> {
    line.split(if comma { ',' } else { '\t' })
}

pub struct VarTable {
    pub genes: Vec<String>,
    pub var: Vec<f64>,
    /// present when the variance file carries a third (df) column.
    pub df: Option<Vec<f64>>,
}

/// gene id, variance, and an optional df column. The header is required; the
/// delimiter is `,` when the path ends in .csv, else tab.
pub fn read_var(path: &Path) -> Result<VarTable> {
    let comma = path
        .extension()
        .map(|e| e.eq_ignore_ascii_case("csv"))
        .unwrap_or(false);
    let mut lines = open(path)?.lines();
    let header = lines
        .next()
        .ok_or_else(|| RsomicsError::InvalidInput("empty variance file".into()))?
        .map_err(RsomicsError::Io)?;
    let ncol = split(&header, comma).count();
    let with_df = ncol >= 3;

    let mut genes = Vec::new();
    let mut var = Vec::new();
    let mut df = if with_df { Some(Vec::new()) } else { None };
    for line in lines {
        let line = line.map_err(RsomicsError::Io)?;
        if line.is_empty() {
            continue;
        }
        let mut f = split(&line, comma);
        let gene = f
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput("missing gene id".into()))?;
        let v = f.next().ok_or_else(|| {
            RsomicsError::InvalidInput(format!("gene '{gene}': missing variance"))
        })?;
        genes.push(gene.to_string());
        var.push(parse_f64(v)?);
        if let Some(d) = df.as_mut() {
            let dv = f.next().ok_or_else(|| {
                RsomicsError::InvalidInput(format!("gene '{gene}': missing df column"))
            })?;
            d.push(parse_f64(dv)?);
        }
    }
    if genes.is_empty() {
        return Err(RsomicsError::InvalidInput("no variances in file".into()));
    }
    Ok(VarTable { genes, var, df })
}

/// A per-gene df column file: gene id, df. Returned in file order with no gene
/// matching — callers pair it positionally with the variance vector.
pub fn read_df_column(path: &Path) -> Result<Vec<f64>> {
    let comma = path
        .extension()
        .map(|e| e.eq_ignore_ascii_case("csv"))
        .unwrap_or(false);
    let mut lines = open(path)?.lines();
    lines.next();
    let mut df = Vec::new();
    for line in lines {
        let line = line.map_err(RsomicsError::Io)?;
        if line.is_empty() {
            continue;
        }
        let v = split(&line, comma)
            .nth(1)
            .ok_or_else(|| RsomicsError::InvalidInput("df file: missing df column".into()))?;
        df.push(parse_f64(v)?);
    }
    Ok(df)
}
