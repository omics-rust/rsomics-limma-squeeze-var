//! limma squeezeVar — empirical-Bayes shrinkage of per-gene variances toward a
//! fitted scaled-F prior.
//!
//! Clean-room reimplementation of limma's squeezeVar / fitFDist following Smyth
//! (2004), Stat Appl Genet Mol Biol 3(1):3, doi:10.2202/1544-6115.1027. No limma
//! (GPL) source was consulted; the moment estimator follows the published paper
//! and is validated black-box against the binary.

mod fitfdist;
mod io;
mod special;

use std::io::{BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

pub use fitfdist::{Df, PriorFit, fit_f_dist, squeeze_var};
pub use io::{VarTable, read_df_column, read_var};

pub enum DfSource<'a> {
    Scalar(f64),
    Column(&'a Path),
}

pub struct Options<'a> {
    pub var: &'a Path,
    pub df: DfSource<'a>,
}

pub struct Results {
    pub genes: Vec<String>,
    pub var_post: Vec<f64>,
    pub var_prior: f64,
    pub df_prior: f64,
}

pub fn run(opts: &Options) -> Result<Results> {
    let table = read_var(opts.var)?;
    let n = table.genes.len();

    let df_vec: Option<Vec<f64>> = match opts.df {
        DfSource::Scalar(_) => table.df,
        DfSource::Column(path) => {
            let d = read_df_column(path)?;
            if d.len() != n {
                return Err(RsomicsError::InvalidInput(format!(
                    "df column has {} rows, variance file has {n}",
                    d.len()
                )));
            }
            Some(d)
        }
    };

    let (var_post, fit) = match &df_vec {
        Some(d) => squeeze_var(&table.var, &Df::per_gene(d)),
        None => {
            let DfSource::Scalar(s) = opts.df else {
                return Err(RsomicsError::InvalidInput(
                    "no per-gene df column found; pass --df <number> or --df <file>".into(),
                ));
            };
            squeeze_var(&table.var, &Df::scalar(s))
        }
    };

    Ok(Results {
        genes: table.genes,
        var_post,
        var_prior: fit.var_prior,
        df_prior: fit.df_prior,
    })
}

pub fn write_results(res: &Results, out: &mut dyn Write) -> Result<()> {
    let mut w = BufWriter::with_capacity(1 << 20, out);
    let mut fmt = ryu::Buffer::new();
    let var_prior = fmt.format(res.var_prior).to_string();
    let df_prior = fmt.format(res.df_prior).to_string();
    writeln!(w, "# var.prior={var_prior} df.prior={df_prior}").map_err(RsomicsError::Io)?;
    writeln!(w, "gene\tvar.post").map_err(RsomicsError::Io)?;
    let mut line = String::with_capacity(64);
    for (gene, &v) in res.genes.iter().zip(&res.var_post) {
        line.clear();
        line.push_str(gene);
        line.push('\t');
        line.push_str(fmt.format(v));
        line.push('\n');
        w.write_all(line.as_bytes()).map_err(RsomicsError::Io)?;
    }
    w.flush().map_err(RsomicsError::Io)?;
    Ok(())
}
