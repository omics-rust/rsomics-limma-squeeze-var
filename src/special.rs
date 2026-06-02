//! digamma / trigamma and Smyth's trigamma_inverse (Smyth 2004, appendix).
//!
//! digamma/trigamma use the asymptotic series with a recurrence shift to x >= 20,
//! matching R's mathlib to ~1e-12. trigamma_inverse is the Newton scheme on the
//! 1/x scale from the fitFDist appendix.

pub fn digamma(mut x: f64) -> f64 {
    let mut result = 0.0;
    while x < 20.0 {
        result -= 1.0 / x;
        x += 1.0;
    }
    let inv = 1.0 / x;
    let inv2 = inv * inv;
    result + x.ln()
        - 0.5 * inv
        - inv2
            * (1.0 / 12.0
                - inv2 * (1.0 / 120.0 - inv2 * (1.0 / 252.0 - inv2 * (1.0 / 240.0 - inv2 / 132.0))))
}

pub fn trigamma(mut x: f64) -> f64 {
    let mut result = 0.0;
    while x < 20.0 {
        result += 1.0 / (x * x);
        x += 1.0;
    }
    let inv = 1.0 / x;
    let inv2 = inv * inv;
    result
        + inv
            * (1.0
                + inv
                    * (0.5
                        + inv
                            * (1.0 / 6.0
                                - inv2 * (1.0 / 30.0 - inv2 * (1.0 / 42.0 - inv2 / 30.0)))))
}

/// Solve trigamma(y) = x for y.
pub fn trigamma_inverse(x: f64) -> f64 {
    if x > 1e7 {
        return 1.0 / x.sqrt();
    }
    if x < 1e-6 {
        return 1.0 / x;
    }
    let mut y = 0.5 + 1.0 / x;
    loop {
        let tri = trigamma(y);
        let dif = tri * (1.0 - tri / x) / psigamma2(y);
        y += dif;
        if -dif / y < 1e-8 {
            break;
        }
    }
    y
}

fn psigamma2(mut x: f64) -> f64 {
    let mut result = 0.0;
    while x < 20.0 {
        result -= 2.0 / (x * x * x);
        x += 1.0;
    }
    let inv = 1.0 / x;
    let x2 = x * x;
    -1.0 / x2 - 1.0 / (x2 * x) - 0.5 * inv * inv * inv * inv + result + inv.powi(6) / 6.0
        - inv.powi(8) / 6.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn digamma_known() {
        assert!((digamma(1.0) + 0.577_215_664_901_532_9).abs() < 1e-10);
        assert!((digamma(0.5) + 1.963_510_026_021_423).abs() < 1e-9);
    }

    #[test]
    fn trigamma_known() {
        assert!((trigamma(1.0) - std::f64::consts::PI.powi(2) / 6.0).abs() < 1e-10);
    }

    #[test]
    fn trigamma_inverse_roundtrip() {
        for &y in &[0.7, 1.3, 4.0, 20.0, 100.0] {
            let yi = trigamma_inverse(trigamma(y));
            assert!((yi - y).abs() / y < 1e-6, "y={y} got {yi}");
        }
    }
}
