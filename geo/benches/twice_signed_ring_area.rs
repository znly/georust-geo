#[macro_use]
extern crate criterion;
extern crate geo;
use geo::{Coordinate, LineString};
// use geo::algorithm::winding_order::twice_signed_ring_area;
use geo::algorithm::winding_order::twice_signed_ring_area2 as twice_signed_ring_area;

fn criterion_benchmark(c: &mut criterion::Criterion) {
    c.bench_function("twice_signed_ring_area f32", |bencher| {
        let points = include!("../src/algorithm/test_fixtures/norway_main.rs");

        let ls = LineString::<f32>(
            points
                .iter()
                .map(|e| Coordinate { x: e[0], y: e[1] })
                .collect()
        );

        bencher.iter(|| {
            criterion::black_box(|| {
                let _ = twice_signed_ring_area(&ls);
            });
        });
    });

    c.bench_function("twice_signed_ring_area f64", |bencher| {
        let points = include!("../src/algorithm/test_fixtures/norway_main.rs");

        let ls = LineString::<f64>(
            points
                .iter()
                .map(|e| Coordinate { x: e[0], y: e[1] })
                .collect()
        );

        bencher.iter(|| {
            criterion::black_box(|| {
                let _ = twice_signed_ring_area(&ls);
            });
        });
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
