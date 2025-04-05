//! testing conformal geometric algebra with zero vectors
//!
//! this test suite demonstrates the zero vector implementation for
//! conformal geometric algebra operations

use geonum::{Geonum, Multivector};
use geopoint::ZeroVector;
use std::f64::consts::{PI, TAU};

#[test]
fn test_conformal_operations() {
    println!("Conformal Geometric Algebra with Zero Vectors Demo");
    println!("==================================================\n");

    // create a conformal zero vector with standard metric (π/2)
    let standard = ZeroVector::standard(3);

    // create a conformal zero vector with 0.0 metric
    let conformal = ZeroVector::conformal(3);

    // create a custom metric (π/4)
    let custom = ZeroVector::new(3, PI / 4.0);

    println!("Created three different metrics:");
    println!("  - Standard (π/2): Euclidean geometry");
    println!("  - Conformal (0.0): Conformal geometry");
    println!("  - Custom (π/4): Custom metric\n");

    // Example 1: Creating and working with points
    println!("Example 1: Creating and working with points");
    println!("------------------------------------------");

    // create a point at [1.0, 2.0, 3.0]
    let coords = [1.0, 2.0, 3.0];
    let point = conformal.point(&coords);

    // extract the coordinates back
    let extracted = conformal.extract_coordinates(&point);
    println!("Original coordinates: {:?}", coords);
    println!("Extracted coordinates: {:?}", extracted);

    // create a second point
    let point2 = conformal.point(&[4.0, 5.0, 6.0]);

    // compute distance between points
    let distance = conformal.distance(&point, &point2);
    println!("Distance between points: {:.2}", distance);
    println!();

    // Example 2: Creating and working with spheres
    println!("Example 2: Creating and working with spheres");
    println!("------------------------------------------");

    // create a sphere at origin with radius 5.0
    let origin = conformal.point(&[0.0, 0.0, 0.0]);
    let sphere = conformal.sphere(&origin, 5.0);

    // test if original point is in the sphere
    let in_sphere = conformal.contains_point(&sphere, &point);
    println!("Is point [1.0, 2.0, 3.0] in sphere? {}", in_sphere);

    // compute the distance from sphere to point
    let distance_to_sphere = if in_sphere {
        "Inside the sphere"
    } else {
        "Outside the sphere"
    };
    println!("Point is {}", distance_to_sphere);
    println!();

    // Example 3: Translation
    println!("Example 3: Translation");
    println!("---------------------");

    // create a translation vector [1.0, 1.0, 1.0]
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

    // translate the sphere
    let translated_center = conformal.translate(&origin, &translation);
    let translated_sphere = conformal.sphere(&translated_center, 5.0);

    // extract the center coordinates
    let new_center = conformal.extract_coordinates(&translated_center);
    println!("Original sphere center: [0.0, 0.0, 0.0]");
    println!("Translated sphere center: {:?}", new_center);

    // check if a point is in the translated sphere
    let test_point = conformal.point(&[1.0, 1.0, 1.0]);
    let in_translated = conformal.contains_point(&translated_sphere, &test_point);
    println!(
        "Is point [1.0, 1.0, 1.0] in translated sphere? {}",
        in_translated
    );
    println!();

    // Example 4: Custom metrics
    println!("Example 4: Custom metrics");
    println!("------------------------");

    // create vectors
    let v1 = geonum::Geonum {
        length: 1.0,
        angle: 0.0,
    };
    let v2 = geonum::Geonum {
        length: 1.0,
        angle: PI / 2.0,
    };

    // compute dot products in different metrics
    let dot_standard = standard.dot(&v1, &v2);
    let dot_conformal = conformal.dot(&v1, &v2);
    let dot_custom = custom.dot(&v1, &v2);

    println!("Dot product of perpendicular vectors in different metrics:");
    println!("  - Standard (π/2): {:.4}", dot_standard);
    println!("  - Conformal (0.0): {:.4}", dot_conformal);
    println!("  - Custom (π/4): {:.4}", dot_custom);

    println!("\nExpected values:");
    println!("  - Standard: cos(π/2) = 0");
    println!("  - Conformal: cos(0) = 1");
    println!("  - Custom: cos(π/4) ≈ 0.7071");
    println!();

    // Example 5: The power of conformal geometry with zero vectors
    println!("Example 5: The power of conformal geometry with zero vectors");
    println!("-------------------------------------------------------");
    println!("With the ZeroVector approach, we can work with any metric");
    println!("by simply changing the angle of the zero vector.");
    println!();
    println!("This allows us to maintain the simplicity of the approach");
    println!("while still gaining all the benefits of conformal geometry.");
    println!();
    println!("The key advantage is that computations stay O(1) in complexity");
    println!("rather than scaling with the number of dimensions.");
}
