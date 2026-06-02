use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_limma_squeeze_var::{Df, squeeze_var};

fn bench(c: &mut Criterion) {
    let g = 200_000;
    let mut x = Vec::with_capacity(g);
    let mut s: u64 = 0x9e37_79b9_7f4a_7c15;
    for _ in 0..g {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        x.push(0.2 + (s >> 11) as f64 / (1u64 << 53) as f64 * 4.0);
    }
    c.bench_function("squeeze_var_200k", |b| {
        b.iter(|| squeeze_var(black_box(&x), &Df::scalar(6.0)))
    });
}

criterion_group!(benches, bench);
criterion_main!(benches);
