use criterion::{black_box, criterion_group, criterion_main, Criterion};
use geonum::{Geonum, Multivector};
use geopoint::*;
use std::f64::consts::{PI, TAU};

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

fn bench_zero_vector_operations(c: &mut Criterion) {
    // Create different zero vectors
    let standard = ZeroVector::standard(3);
    let conformal = ZeroVector::conformal(3);
    let custom = ZeroVector::new(3, PI / 4.0);

    // Benchmark point creation
    c.bench_function("zv_standard_point", |b| {
        b.iter(|| black_box(standard.point(&[1.0, 2.0, 3.0])))
    });

    c.bench_function("zv_conformal_point", |b| {
        b.iter(|| black_box(conformal.point(&[1.0, 2.0, 3.0])))
    });

    c.bench_function("zv_custom_point", |b| {
        b.iter(|| black_box(custom.point(&[1.0, 2.0, 3.0])))
    });

    // Create points and spheres for benchmark
    let center = conformal.point(&[0.0, 0.0, 0.0]);
    let point = conformal.point(&[1.0, 1.0, 1.0]);
    let sphere = conformal.sphere(&center, 2.0);

    // Benchmark sphere creation
    c.bench_function("zv_create_sphere", |b| {
        b.iter(|| black_box(conformal.sphere(&center, 2.0)))
    });

    // Benchmark contains check
    c.bench_function("zv_sphere_contains_point", |b| {
        b.iter(|| black_box(conformal.contains_point(&sphere, &point)))
    });

    // Benchmark distance calculation
    c.bench_function("zv_distance", |b| {
        b.iter(|| black_box(conformal.distance(&center, &point)))
    });

    // Create a translation vector
    let mut translation = Multivector::new();
    translation.0.push(Geonum {
        length: 1.0,
        angle: (conformal.angle + PI / 2.0) % TAU, // x component
    });
    translation.0.push(Geonum {
        length: 1.0,
        angle: (conformal.angle + PI) % TAU, // y component
    });
    translation.0.push(Geonum {
        length: 1.0,
        angle: (conformal.angle + 3.0 * PI / 2.0) % TAU, // z component
    });

    // Benchmark translation
    c.bench_function("zv_translate_point", |b| {
        b.iter(|| black_box(conformal.translate(&point, &translation)))
    });

    // Benchmark custom metric dot product
    let v1 = geonum::Geonum {
        length: 1.0,
        angle: 0.0,
    };
    let v2 = geonum::Geonum {
        length: 1.0,
        angle: PI / 2.0,
    };
    c.bench_function("zv_custom_dot_product", |b| {
        b.iter(|| black_box(custom.dot(&v1, &v2)))
    });
}

fn bench_high_dimensional_zero_vector(c: &mut Criterion) {
    // Create high-dimensional zero vectors (10D)
    let _high_dim_standard = ZeroVector::standard(10);
    let high_dim_conformal = ZeroVector::conformal(10);

    // Create coordinates
    let mut coords = Vec::with_capacity(10);
    for i in 0..10 {
        coords.push(i as f64 / 10.0);
    }

    // Benchmark high-dimensional point creation
    c.bench_function("zv_high_dim_point", |b| {
        b.iter(|| black_box(high_dim_conformal.point(&coords)))
    });

    // Create high-dimensional point and sphere
    let point = high_dim_conformal.point(&coords);
    let origin = high_dim_conformal.point(&vec![0.0; 10]);

    c.bench_function("zv_high_dim_sphere", |b| {
        b.iter(|| black_box(high_dim_conformal.sphere(&origin, 1.0)))
    });

    // Benchmark coordinate extraction
    c.bench_function("zv_high_dim_extract_coords", |b| {
        b.iter(|| black_box(high_dim_conformal.extract_coordinates(&point)))
    });
}

fn bench_approach_comparison(c: &mut Criterion) {
    // Traditional approach
    let trad_point = Multivector::point(&[1.0, 2.0, 3.0]);
    let _trad_conformal = trad_point.to_conformal_point();
    let trad_center = Multivector::point(&[0.0, 0.0, 0.0]);
    let trad_sphere = Multivector::sphere(&trad_center, 2.0);

    // Zero Vector approach
    let zv = ZeroVector::conformal(3);
    let zv_point = zv.point(&[1.0, 2.0, 3.0]);
    let zv_center = zv.point(&[0.0, 0.0, 0.0]);
    let zv_sphere = zv.sphere(&zv_center, 2.0);

    // Compare point creation
    c.bench_function("compare_trad_point_create", |b| {
        b.iter(|| black_box(Multivector::point(&[1.0, 2.0, 3.0])))
    });

    c.bench_function("compare_zv_point_create", |b| {
        b.iter(|| black_box(zv.point(&[1.0, 2.0, 3.0])))
    });

    // Compare conformal conversion
    c.bench_function("compare_trad_to_conformal", |b| {
        b.iter(|| black_box(trad_point.to_conformal_point()))
    });

    // Compare sphere creation
    c.bench_function("compare_trad_sphere_create", |b| {
        b.iter(|| black_box(Multivector::sphere(&trad_center, 2.0)))
    });

    c.bench_function("compare_zv_sphere_create", |b| {
        b.iter(|| black_box(zv.sphere(&zv_center, 2.0)))
    });

    // Compare containment check
    c.bench_function("compare_trad_contains", |b| {
        b.iter(|| black_box(trad_sphere.contains_point(&trad_point)))
    });

    c.bench_function("compare_zv_contains", |b| {
        b.iter(|| black_box(zv.contains_point(&zv_sphere, &zv_point)))
    });
}

criterion_group!(
    benches,
    bench_conformal_point_conversion,
    bench_sphere_operations,
    bench_high_dimensional_operations,
    bench_geometric_transformations,
    bench_zero_vector_operations,
    bench_high_dimensional_zero_vector,
    bench_approach_comparison
);
criterion_main!(benches);
