//! fitFDist + the squeezeVar posterior (Smyth 2004, doi:10.2202/1544-6115.1027).
//!
//! Method-of-moments on the trigamma scale: the per-gene residual variances are
//! treated as scale * F(df1, df2); df2 (= df.prior) and scale (= var.prior) come
//! from the mean and variance of log(variance) corrected by digamma/trigamma.

use crate::special::{digamma, trigamma, trigamma_inverse};

pub struct Df<'a> {
    scalar: f64,
    per_gene: Option<&'a [f64]>,
}

impl<'a> Df<'a> {
    pub fn scalar(df: f64) -> Self {
        Df {
            scalar: df,
            per_gene: None,
        }
    }
    pub fn per_gene(df: &'a [f64]) -> Self {
        Df {
            scalar: f64::NAN,
            per_gene: Some(df),
        }
    }
    fn get(&self, i: usize) -> f64 {
        match self.per_gene {
            Some(v) => v[i],
            None => self.scalar,
        }
    }
}

pub struct PriorFit {
    pub var_prior: f64,
    pub df_prior: f64,
}

/// The scaled-F prior fitted to `x` with residual degrees of freedom `df1`.
pub fn fit_f_dist(x: &[f64], df1: &Df) -> PriorFit {
    let n = x.len();
    if n == 0 {
        return PriorFit {
            var_prior: f64::NAN,
            df_prior: f64::NAN,
        };
    }
    if n == 1 {
        return PriorFit {
            var_prior: x[0],
            df_prior: 0.0,
        };
    }

    let ok: Vec<usize> = (0..n)
        .filter(|&i| {
            let d = df1.get(i);
            d.is_finite() && d > 1e-15 && x[i].is_finite() && x[i] > -1e-15
        })
        .collect();
    let nok = ok.len();
    if nok == 0 {
        return PriorFit {
            var_prior: f64::NAN,
            df_prior: f64::NAN,
        };
    }
    if nok == 1 {
        return PriorFit {
            var_prior: x[ok[0]],
            df_prior: 0.0,
        };
    }

    let mut xs: Vec<f64> = ok.iter().map(|&i| x[i].max(0.0)).collect();
    let m = median(&xs);
    let m = if m == 0.0 { 1.0 } else { m };
    let floor = 1e-5 * m;
    for v in &mut xs {
        *v = v.max(floor);
    }

    let nf = nok as f64;
    let mut emean = 0.0;
    let mut tri_mean = 0.0;
    let mut e = Vec::with_capacity(nok);
    for (k, &i) in ok.iter().enumerate() {
        let half = df1.get(i) / 2.0;
        let ei = xs[k].ln() - digamma(half) + half.ln();
        emean += ei;
        tri_mean += trigamma(half);
        e.push(ei);
    }
    emean /= nf;
    tri_mean /= nf;

    let evar = e.iter().map(|&v| (v - emean).powi(2)).sum::<f64>() / (nf - 1.0) - tri_mean;

    if evar > 0.0 {
        let df2 = 2.0 * trigamma_inverse(evar);
        let scale = (emean + digamma(df2 / 2.0) - (df2 / 2.0).ln()).exp();
        PriorFit {
            var_prior: scale,
            df_prior: df2,
        }
    } else {
        PriorFit {
            var_prior: xs.iter().sum::<f64>() / nf,
            df_prior: f64::INFINITY,
        }
    }
}

/// var.post: shrink each variance toward the prior by the posterior weights.
pub fn squeeze_var(x: &[f64], df1: &Df) -> (Vec<f64>, PriorFit) {
    // Fewer than three genes carry no information to fit a prior: limma leaves
    // them unshrunk (var.post = var, df.prior = 0) rather than fitting noise.
    if x.len() < 3 {
        return (
            x.to_vec(),
            PriorFit {
                var_prior: x.first().copied().unwrap_or(f64::NAN),
                df_prior: 0.0,
            },
        );
    }
    let fit = fit_f_dist(x, df1);
    let post: Vec<f64> = if fit.df_prior.is_infinite() {
        vec![fit.var_prior; x.len()]
    } else {
        (0..x.len())
            .map(|i| {
                let d = df1.get(i);
                (d * x[i] + fit.df_prior * fit.var_prior) / (d + fit.df_prior)
            })
            .collect()
    };
    (post, fit)
}

fn median(v: &[f64]) -> f64 {
    let mut s: Vec<f64> = v.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = s.len();
    if n % 2 == 1 {
        s[n / 2]
    } else {
        0.5 * (s[n / 2 - 1] + s[n / 2])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn posterior_is_convex_combination() {
        let x = vec![0.5, 1.2, 0.8, 2.0, 0.3, 1.5, 0.9, 1.1, 0.7, 1.8];
        let (post, fit) = squeeze_var(&x, &Df::scalar(4.0));
        if fit.df_prior.is_finite() {
            for (i, &p) in post.iter().enumerate() {
                let lo = x[i].min(fit.var_prior);
                let hi = x[i].max(fit.var_prior);
                assert!(p >= lo - 1e-12 && p <= hi + 1e-12);
            }
        } else {
            assert!(post.iter().all(|&p| (p - fit.var_prior).abs() < 1e-12));
        }
    }

    #[test]
    fn single_gene() {
        let (post, fit) = squeeze_var(&[2.5], &Df::scalar(4.0));
        assert_eq!(fit.df_prior, 0.0);
        assert_eq!(fit.var_prior, 2.5);
        assert_eq!(post[0], 2.5);
    }

    #[test]
    fn two_genes_unshrunk() {
        // limma leaves fewer than three genes untouched: var.post = var, df.prior = 0.
        let x = vec![0.8, 2.4];
        let (post, fit) = squeeze_var(&x, &Df::scalar(4.0));
        assert_eq!(fit.df_prior, 0.0);
        assert_eq!(post, x);
    }

    #[test]
    fn nonpositive_evar_uses_arithmetic_mean() {
        // Equal variances drive evar<=0; limma's Inf-df branch sets var.prior to
        // mean(var) (here 1.5), not exp(mean(log var)).
        let x = vec![1.5, 1.5, 1.5, 1.5, 1.5];
        let (post, fit) = squeeze_var(&x, &Df::scalar(4.0));
        assert!(fit.df_prior.is_infinite());
        assert!((fit.var_prior - 1.5).abs() < 1e-12, "{}", fit.var_prior);
        assert!(post.iter().all(|&p| (p - 1.5).abs() < 1e-12));
    }

    #[test]
    fn nonpositive_evar_mean_differs_from_geometric_mean() {
        // Spread variances with small df still hit evar<=0; the arithmetic mean
        // (1.08) is distinct from the geometric mean that the old code emitted.
        let x = vec![0.5, 1.2, 0.8, 2.0, 0.3, 1.5, 0.9, 1.1, 0.7, 1.8];
        let (_post, fit) = squeeze_var(&x, &Df::scalar(6.0));
        assert!(fit.df_prior.is_infinite());
        let mean = x.iter().sum::<f64>() / x.len() as f64;
        assert!((fit.var_prior - mean).abs() < 1e-12, "{}", fit.var_prior);
    }
}
