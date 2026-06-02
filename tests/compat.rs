//! Differential compat against limma squeezeVar.
//!
//! - `golden_*` always runs: ours vs a committed R-captured output.
//! - `live_r_*` runs only when an Rscript with limma is found; it regenerates
//!   the oracle and diffs against ours (loud-skip otherwise).

use std::path::PathBuf;
use std::process::Command;

const EPS: f64 = 1e-6; // relative, on var.post and the header scalars

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-limma-squeeze-var"))
}

fn golden(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn manifest(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(rel)
}

struct Table {
    var_prior: f64,
    df_prior: f64,
    post: Vec<(String, f64)>,
}

fn parse(text: &str) -> Table {
    let mut var_prior = f64::NAN;
    let mut df_prior = f64::NAN;
    let mut post = Vec::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("# ") {
            for tok in rest.split_whitespace() {
                if let Some(v) = tok.strip_prefix("var.prior=") {
                    var_prior = v.parse().unwrap();
                } else if let Some(v) = tok.strip_prefix("df.prior=") {
                    df_prior = v.parse().unwrap();
                }
            }
            continue;
        }
        if line.starts_with("gene\t") || line.is_empty() {
            continue;
        }
        let mut f = line.split('\t');
        let gene = f.next().unwrap().to_string();
        let v: f64 = f.next().unwrap().trim().parse().unwrap();
        post.push((gene, v));
    }
    Table {
        var_prior,
        df_prior,
        post,
    }
}

fn rel(a: f64, b: f64) -> f64 {
    if a.is_infinite() && b.is_infinite() {
        return 0.0;
    }
    (a - b).abs() / b.abs().max(1e-9)
}

fn assert_close(a: &Table, b: &Table, label: &str) {
    assert!(
        rel(a.var_prior, b.var_prior) < EPS,
        "{label}: var.prior ours={} ref={}",
        a.var_prior,
        b.var_prior
    );
    assert!(
        rel(a.df_prior, b.df_prior) < EPS,
        "{label}: df.prior ours={} ref={}",
        a.df_prior,
        b.df_prior
    );
    assert_eq!(a.post.len(), b.post.len(), "{label}: row count mismatch");
    let mut max_rel = 0.0f64;
    for ((ga, va), (gb, vb)) in a.post.iter().zip(&b.post) {
        assert_eq!(ga, gb, "{label}: gene order mismatch");
        let r = rel(*va, *vb);
        max_rel = max_rel.max(r);
        assert!(r < EPS, "{label}: {ga} ours={va} ref={vb} rel={r:e}");
    }
    eprintln!("{label}: max relative var.post deviation = {max_rel:e}");
}

fn run_ours(vars: &str, df: &str) -> String {
    let mut cmd = Command::new(ours());
    cmd.args(["--var", golden(vars).to_str().unwrap()]);
    if df.parse::<f64>().is_ok() {
        cmd.args(["--df", df]);
    } else {
        cmd.args(["--df-column", golden(df).to_str().unwrap()]);
    }
    let out = cmd.output().unwrap();
    assert!(
        out.status.success(),
        "ours failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

#[test]
fn golden_scalar_df() {
    let ours_out = run_ours("vars.tsv", "6");
    let expected = std::fs::read_to_string(golden("post.scalar.expected.tsv")).unwrap();
    assert_close(&parse(&ours_out), &parse(&expected), "scalar df (golden)");
}

#[test]
fn golden_pergene_df() {
    let ours_out = run_ours("vars_pg.tsv", "df_pg.tsv");
    let expected = std::fs::read_to_string(golden("post.pergene.expected.tsv")).unwrap();
    assert_close(&parse(&ours_out), &parse(&expected), "per-gene df (golden)");
}

fn rscript() -> Option<String> {
    let conda = format!(
        "{}/miniconda3/envs/r-bioc/bin/Rscript",
        std::env::var("HOME").unwrap_or_default()
    );
    for cand in [conda.as_str(), "Rscript"] {
        let ok = Command::new(cand)
            .args([
                "-e",
                "if(!requireNamespace('limma',quietly=TRUE)) quit(status=1)",
            ])
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false);
        if ok {
            return Some(cand.to_string());
        }
    }
    None
}

fn live_diff(rs: &str, vars: &str, df: &str, label: &str) {
    let scratch = std::env::temp_dir();
    let r_out = scratch.join(format!("squeezevar_r_{}_{}.tsv", label, std::process::id()));
    let df_arg = if df.parse::<f64>().is_ok() {
        df.to_string()
    } else {
        golden(df).to_str().unwrap().to_string()
    };
    let oracle = Command::new(rs)
        .arg(manifest("tests/oracle.R"))
        .arg(golden(vars))
        .arg(&df_arg)
        .arg(&r_out)
        .output()
        .unwrap();
    assert!(
        oracle.status.success(),
        "oracle failed: {}",
        String::from_utf8_lossy(&oracle.stderr)
    );
    let ours_out = run_ours(vars, df);
    let r_ref = std::fs::read_to_string(&r_out).unwrap();
    let _ = std::fs::remove_file(&r_out);
    assert_close(&parse(&ours_out), &parse(&r_ref), label);
}

#[test]
fn live_r_scalar_df() {
    let Some(rs) = rscript() else {
        eprintln!("SKIP live_r_scalar_df: no Rscript with limma found");
        return;
    };
    live_diff(&rs, "vars.tsv", "6", "scalar df (live R)");
}

#[test]
fn live_r_pergene_df() {
    let Some(rs) = rscript() else {
        eprintln!("SKIP live_r_pergene_df: no Rscript with limma found");
        return;
    };
    live_diff(&rs, "vars_pg.tsv", "df_pg.tsv", "per-gene df (live R)");
}
