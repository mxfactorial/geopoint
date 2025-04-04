[![build](https://github.com/mxfactorial/geopoint/actions/workflows/publish.yaml/badge.svg)](https://github.com/mxfactorial/geopoint/actions)
[![Discord](https://img.shields.io/discord/868565277955203122.svg?label=&logo=discord&logoColor=ffffff&color=7389D8&labelColor=6A7EC2)](https://discord.gg/KQdC65bG)
[![docs](https://docs.rs/geopoint/badge.svg)](https://docs.rs/geopoint)
[![crates.io](https://img.shields.io/crates/v/geopoint.svg)](https://crates.io/crates/geopoint)
[![coverage](https://coveralls.io/repos/github/mxfactorial/geopoint/badge.svg?branch=main)](https://coveralls.io/github/mxfactorial/geopoint?branch=main)
[![contribute](https://img.shields.io/badge/contribute-paypal-brightgreen.svg)](https://www.paypal.com/paypalme/mxfactorial)

# geopoint

conformal geometric algebra built on top of the O(1) [geonum](https://crates.io/crates/geonum) crate

conformal geometric algebra extends euclidean geometric algebra by adding two basis vectors: e₀ (origin) and e∞ (infinity). this allows elegant representation of primitives like:

- points: p = e₀ + x + ½x²e∞ (where x is euclidean position)
- spheres: s = p - ½r²e∞ (where p is center point, r is radius)
- planes: π = n + de∞ (where n is normal vector, d is distance from origin)

all transformations can be expressed through versors (geometric algebra operators) applied via sandwich products

### features

- conformal point representation for euclidean geometry
- geometric primitives: points, lines, planes, spheres, circles
- geometric operations: containment tests, distance computation 
- conformal transformations: translation, rotation, dilation, inversion
- versor implementation for efficient transformations

### use

```
cargo add geopoint
```

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

see `tests/lib_test.rs` for more examples of conformal geometric algebra operations

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

```
cargo bench
```

### test

```
cargo fmt --check # format
cargo clippy # lint
cargo test --lib # unit
cargo test --test lib_test # feature
cargo bench # bench
cargo llvm-cov # coverage
```

### docs

```
cargo doc --open
```

### todo

- **null basis handling**
  - implement null basis elements e+ and e- (e+ = (e₀ + e∞)/2, e- = (e₀ - e∞)/2)
  - create null basis representation satisfying e+² = 0, e-² = 0, e+·e- = -0.5
  - optimize operations using the null basis representation

- **blade identification and extraction**
  - implement `detect_primitive_type()` to identify geometric entity type
  - add detection for spheres, planes, lines, circles, and point pairs
  - extract grade components from multivectors based on angle patterns
  - add `PrimitiveType` enum for classifying geometric objects

- **meet and join operations**
  - implement meet (intersection) operation using the dual approach: X ∧ Y = (X* ∨ Y*)*
  - implement join (union) operation using the wedge product
  - add grade handling and decomposition for outer products
  - implement helpers `highest_grade()` and `decompose_by_grade()`

- **versor operations for transformations**
  - implement sandwich product for versor-based transformations
  - add conjugate method for versors (R⁻¹ = R̃)
  - create geometric product with [length, angle] handling
  - implement optimizations for common transformations

- **optimizations**
  - component compression - merge similar angle components
  - lazy evaluation - delay computations until needed
  - angle pattern caching - build lookup tables for common patterns
  - specialized storage for high-dimensional spaces
  - parallel computation for large multivectors

- **benchmarks and metrics**
  - implement sphere intersection benchmarks
  - add conformal translation performance tests
  - create high-dimensional sphere benchmarks (1000+ dimensions)
  - compare against traditional implementations to show O(1) advantage

- **grade projection and selection** - extract specific grades from multivectors
- **conformal conversions** - optimize conversion between euclidean and conformal
- **motor operations** - unified system for rotation and translation
- **distance and angle calculations** - compute geometric measurements
- **dualization and projection** - implement dual operations and projection operators
- **conformal integration with geonum** - integrate with geonum's O(1) representation