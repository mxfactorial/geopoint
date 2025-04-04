use criterion::{black_box, criterion_group, criterion_main, Criterion};
use geonum::Multivector;
use geopoint::*;

fn bench_conformal_point_conversion(c: &mut Criterion) {
    let point = Multivector::point(&[1.0, 2.0, 3.0]);

    c.bench_function("point_to_conformal", |b| {
        b.iter(|| black_box(point.to_conformal_point()))
    });

    let conformal_point = point.to_conformal_point();
    c.bench_function("conformal_to_point", |b| {
        b.iter(|| black_box(conformal_point.extract_euclidean()))
    });
}

fn bench_sphere_operations(c: &mut Criterion) {
    let center = Multivector::point(&[0.0, 0.0, 0.0]);

    c.bench_function("create_sphere", |b| {
        b.iter(|| black_box(Multivector::sphere(&center, 2.0)))
    });

    let sphere = Multivector::sphere(&center, 2.0);
    let point = Multivector::point(&[1.0, 1.0, 1.0]);

    c.bench_function("sphere_contains_point", |b| {
        b.iter(|| black_box(sphere.contains_point(&point)))
    });
}

fn bench_high_dimensional_operations(c: &mut Criterion) {
    // create high-dimensional point (10 dimensions)
    let mut coords = Vec::with_capacity(10);
    for i in 0..10 {
        coords.push(i as f64 / 10.0);
    }

    let high_dim_point = Multivector::point(&coords);

    c.bench_function("high_dim_to_conformal", |b| {
        b.iter(|| black_box(high_dim_point.to_conformal_point()))
    });

    c.bench_function("high_dim_sphere", |b| {
        b.iter(|| black_box(Multivector::sphere(&high_dim_point, 1.0)))
    });
}

fn bench_geometric_transformations(c: &mut Criterion) {
    let point = Multivector::point(&[1.0, 2.0, 3.0]);
    let translation = Multivector::point(&[2.0, -1.0, 5.0]);

    c.bench_function("translate_point", |b| {
        b.iter(|| black_box(point.translate(&translation)))
    });

    c.bench_function("dilate_point", |b| b.iter(|| black_box(point.dilate(2.0))));

    // create a rotation axis (z-axis)
    let axis = Multivector::point(&[0.0, 0.0, 1.0]);

    c.bench_function("create_rotor", |b| {
        b.iter(|| black_box(Multivector::rotor(&axis, std::f64::consts::PI / 4.0)))
    });

    let rotor = Multivector::rotor(&axis, std::f64::consts::PI / 4.0);

    c.bench_function("rotate_point", |b| {
        b.iter(|| black_box(point.rotate(&rotor)))
    });
}

criterion_group!(
    benches,
    bench_conformal_point_conversion,
    bench_sphere_operations,
    bench_high_dimensional_operations,
    bench_geometric_transformations
);
criterion_main!(benches);
