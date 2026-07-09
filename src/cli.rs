use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_limma_squeeze_var::{DfSource, Options, run, write_results};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-limma-squeeze-var", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// Variances TSV/CSV: header, col 1 = gene id, col 2 = variance, optional
    /// col 3 = per-gene df.
    #[arg(long)]
    var: PathBuf,
    /// Residual degrees of freedom, shared across genes.
    #[arg(long, default_value_t = 1.0)]
    df: f64,
    /// Per-gene df file (gene id + df); overrides --df and any inline df column.
    #[arg(long)]
    df_column: Option<PathBuf>,
    /// Output TSV; "-" is stdout.
    #[arg(short = 'o', long, default_value = "-")]
    output: String,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let df = match &self.df_column {
            Some(p) => DfSource::Column(p),
            None => DfSource::Scalar(self.df),
        };
        let opts = Options { var: &self.var, df };
        let res = run(&opts)?;

        let mut out: Box<dyn std::io::Write> = if self.output == "-" && self.common.json {
            Box::new(std::io::sink())
        } else if self.output == "-" {
            Box::new(std::io::stdout().lock())
        } else {
            Box::new(std::fs::File::create(&self.output).map_err(RsomicsError::Io)?)
        };
        write_results(&res, &mut out)?;

        if !self.common.quiet {
            eprintln!(
                "{} genes, var.prior={:.6} df.prior={:.6}",
                res.genes.len(),
                res.var_prior,
                res.df_prior
            );
        }
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "Empirical-Bayes shrinkage of gene variances toward a fitted scaled-F prior.",
    origin: Some(Origin {
        upstream: "limma squeezeVar / fitFDist",
        upstream_license: "GPL (>=2)",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.2202/1544-6115.1027"),
    }),
    usage_lines: &["--var vars.tsv --df 4 [--df-column df.tsv] [-o out.tsv]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: None,
                long: "var",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: true,
                default: None,
                description: "Variances TSV/CSV (gene, variance, optional df column).",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "df",
                aliases: &[],
                value: Some("<n>"),
                type_hint: Some("f64"),
                required: false,
                default: Some("1"),
                description: "Residual degrees of freedom shared across genes.",
                why_default: Some("Placeholder; a real residual df should be supplied."),
            },
            FlagSpec {
                short: None,
                long: "df-column",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: false,
                default: None,
                description: "Per-gene df file; overrides --df and any inline df column.",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("String"),
                required: false,
                default: Some("-"),
                description: "Output TSV; \"-\" is stdout.",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Shrink variances with a shared residual df of 4",
            command: "rsomics-limma-squeeze-var --var vars.tsv --df 4 -o post.tsv",
        },
        Example {
            description: "Per-gene degrees of freedom",
            command: "rsomics-limma-squeeze-var --var vars.tsv --df-column df.tsv > post.tsv",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
