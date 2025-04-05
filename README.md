[![build](https://github.com/mxfactorial/geopoint/actions/workflows/publish.yaml/badge.svg)](https://github.com/mxfactorial/geopoint/actions)
[![Discord](https://img.shields.io/discord/868565277955203122.svg?label=&logo=discord&logoColor=ffffff&color=7389D8&labelColor=6A7EC2)](https://discord.gg/KQdC65bG)
[![docs](https://docs.rs/geopoint/badge.svg)](https://docs.rs/geopoint)
[![crates.io](https://img.shields.io/crates/v/geopoint.svg)](https://crates.io/crates/geopoint)
[![coverage](https://coveralls.io/repos/github/mxfactorial/geopoint/badge.svg?branch=main)](https://coveralls.io/github/mxfactorial/geopoint?branch=main)
[![contribute](https://img.shields.io/badge/contribute-paypal-brightgreen.svg)](https://www.paypal.com/paypalme/mxfactorial)

# geopoint

conformal geometric algebra built on top of the O(1) [geonum](https://crates.io/crates/geonum) crate

## features

- conformal point representation for euclidean geometry
- geometric primitives: points, lines, planes, spheres, circles
- geometric operations: containment tests, distance computation 
- conformal transformations: translation, rotation, dilation, inversion
- versor implementation for efficient transformations

## use

```
cargo add geopoint
```

### basic

```rust
use geopoint::*;
use geonum::Multivector;

// create a euclidean point
let point = Multivector::point(&[1.0, 2.0, 3.0]);

// convert to conformal representation
let conf_point = point.to_conformal_point();

// extract euclidean point from conformal representation
let euclidean_point = conf_point.extract_euclidean();

// create a sphere
let center = Multivector::point(&[0.0, 0.0, 0.0]);
let sphere = Multivector::sphere(&center, 2.0);

// test if point is inside sphere
let is_inside = sphere.contains_point(&point);

// compute distance to sphere
let distance = sphere.distance_to_point(&point);

// apply translation
let translation = Multivector::point(&[1.0, 0.0, 0.0]);
let translated_point = point.translate(&translation);

// apply uniform scaling
let scaled_point = point.dilate(2.0);
```

### custom standard angle

```rust
use geopoint::ZeroVector;
use std::f64::consts::PI as pi;

// create a zero vector with a custom standard angle (pi/4 angle)
let custom_metric = ZeroVector::new(3, pi/4.0);

// work with this custom standard angle
let point = custom_metric.point(&[1.0, 2.0, 3.0]);

// standard angle (pi/2) matches euclidean geometry
let standard = ZeroVector::standard(3);

// conformal angle (0.0) for conformal transformations
let conformal = ZeroVector::conformal(3);
```

see:
- `tests/lib_test.rs` for examples of basic usage
- `tests/zero_vector_test.rs` for examples with custom metrics
- `tests/conformal_test.rs` for a complete demonstration

### benches

geopoint maintains geonums O(1) computational advantage in high dimensional spaces

#### core operations (nanoseconds)

| operation | time |
|-----------|------|
| point_to_conformal | 125 ns |
| conformal_to_point | 58 ns |
| create_sphere | 254 ns |
| sphere_contains_point | 214 ns |
| translate_point | 136 ns |
| dilate_point | 55 ns |
| create_rotor | 56 ns |
| rotate_point | 996 ns |

#### dimension scaling

| operation | dimensions | time | traditional complexity |
|-----------|------------|------|-----------------------|
| to_conformal (3D) | 3 | 125 ns | O(n) |
| to_conformal (10D) | 10 | 222 ns | O(n) |
| sphere creation (3D) | 3 | 254 ns | O(2^n) |
| sphere creation (10D) | 10 | 409 ns | O(2^n) |

geopoint demonstrates excellent performance with operations completing in hundreds of nanoseconds. high-dimensional operations scale efficiently, showing only a small increase in execution time as dimensions grow - a key advantage of the underlying O(1) representation from geonum

### test

```
cargo fmt --check # format
cargo clippy # lint
cargo test --lib # unit
cargo test --test "*" # feature
cargo bench # bench
cargo llvm-cov # coverage
```

### docs

```
cargo doc --open
```

### todo

- complete coordinate extraction for all dimensions 
- fix rotation with custom angle
- implement meet and join operations
- add blade identification
- optimize high-dimensional operations