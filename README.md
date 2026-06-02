# rsomics-limma-squeeze-var

Empirical-Bayes shrinkage of a vector of per-gene variances toward a fitted
scaled-F prior — a clean-room Rust reimplementation of limma's `squeezeVar`
(the `fitFDist` moment estimator plus the posterior shrinkage step).

Given per-gene residual variances and their residual degrees of freedom, it
fits a scaled inverse-chi-square prior by method of moments on the trigamma
scale and returns each gene's posterior variance:

```
var.post = (df * var + df.prior * var.prior) / (df + df.prior)
```

with `var.prior` (the prior scale) and `df.prior` (the prior degrees of freedom)
reported in the output header. When the variances are too tight to distinguish
from a fixed prior, `df.prior` is `Inf` and every posterior collapses to
`var.prior`.

## Usage

```
rsomics-limma-squeeze-var --var vars.tsv --df 4 -o post.tsv
rsomics-limma-squeeze-var --var vars.tsv --df-column df.tsv > post.tsv
```

`--var` is a TSV (or `.csv`) with a header, gene id in column 1, variance in
column 2, and an optional per-gene df in column 3. `--df` supplies a shared
residual df; `--df-column` supplies a per-gene df file (overrides `--df` and any
inline column). Output is a `# var.prior=… df.prior=…` header line followed by a
`gene\tvar.post` table in input order.

## Scope

This crate implements the classic Smyth (2004) `squeezeVar` / `fitFDist`
estimator — `var.prior` and `df.prior` are scalars and every gene shares one
fitted prior. It is value-exact against limma for both a shared `--df` and a
per-gene df column.

For a shared df this is limma's default. For *unequal* per-gene df, modern limma
(>= 3.50) defaults to a different estimator, `fitFDistUnequalDF1`, which returns
a per-gene shrunk df; this crate reproduces the legacy estimator
(`squeezeVar(..., legacy = TRUE)`) in that case, not `fitFDistUnequalDF1`. The
`robust = TRUE` path (`fitFDistRobustly`, Phipson 2016) and the covariate/trend
path (spline-fitted prior) are separate, larger estimators and are out of scope
for this crate.

## Origin

This crate is an independent Rust reimplementation of `limma::squeezeVar` /
`limma::fitFDist` based on:

- The published method: Smyth (2004), "Linear models and empirical Bayes methods
  for assessing differential expression in microarray experiments", Statistical
  Applications in Genetics and Molecular Biology 3(1):3,
  doi:[10.2202/1544-6115.1027](https://doi.org/10.2202/1544-6115.1027).
- The documented behaviour of the function (argument semantics, the
  scaled-F moment estimator, the `trigammaInverse` Newton scheme described in
  the paper's appendix).
- Black-box behaviour testing against the limma binary (limma 3.66.0).

No source code from limma (GPL) was used as reference during implementation.
Test fixtures are independently generated. The `digamma`/`trigamma` asymptotic
series and the `trigammaInverse` Newton iteration are standard numerical
primitives reconstructed from their mathematical definitions and validated
bit-for-bit against the oracle.

License: MIT OR Apache-2.0.
Upstream credit: limma (<https://bioconductor.org/packages/limma/>, GPL >= 2).
