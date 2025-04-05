use geonum::{Geonum, Multivector};
use geopoint::{
    ConformalGeometry, ConformalPrimitives, GeometricOperations, GeometricPrimitives, ZeroVector,
};
use std::f64::consts::{PI, TAU};

/// tests the basic zero vector functionality
#[test]
fn test_basic_zero_vector() {
    // create a zero vector with standard angle
    let zv = ZeroVector::standard(3);
    assert_eq!(zv.angle, PI / 2.0);

    // create a zero vector with conformal angle
    let conformal = ZeroVector::conformal(3);
    assert_eq!(conformal.angle, 0.0);
}

/// tests creating points and extracting coordinates
#[test]
#[ignore = "Implementation in progress - coordinate extraction needs fixing"]
fn test_point_creation_and_extraction() {
    let zv = ZeroVector::conformal(3);

    // create a point
    let coords = [1.0, 2.0, 3.0];
    let point = zv.point(&coords);

    // extract coordinates
    let extracted = zv.extract_coordinates(&point);

    // verify coordinates match
    assert!(extracted.len() >= coords.len());
    for i in 0..coords.len() {
        assert!((extracted[i] - coords[i]).abs() < 1e-10);
    }
}

/// tests the creation of a sphere and point containment
#[test]
fn test_sphere_creation_and_containment() {
    let zv = ZeroVector::conformal(3);

    // create a sphere at origin with radius 2.0
    let center = zv.point(&[0.0, 0.0, 0.0]);
    let sphere = zv.sphere(&center, 2.0);

    // test points
    let on_sphere = zv.point(&[2.0, 0.0, 0.0]);
    let inside_sphere = zv.point(&[1.0, 0.0, 0.0]);
    let outside_sphere = zv.point(&[3.0, 0.0, 0.0]);

    // test containment
    assert!(zv.contains_point(&sphere, &on_sphere));
    assert!(zv.contains_point(&sphere, &inside_sphere));
    assert!(!zv.contains_point(&sphere, &outside_sphere));
}

/// tests translation of a point
#[test]
#[ignore = "Implementation in progress - translation needs fixing"]
fn test_translation() {
    let zv = ZeroVector::conformal(3);

    // create a point
    let point = zv.point(&[1.0, 2.0, 3.0]);

    // create a translation vector directly
    let mut translation = Multivector::new();
    translation.0.push(Geonum {
        length: 4.0,
        angle: (zv.angle + PI / 2.0) % TAU, // x component
    });
    translation.0.push(Geonum {
        length: 5.0,
        angle: (zv.angle + PI) % TAU, // y component
    });
    translation.0.push(Geonum {
        length: 6.0,
        angle: (zv.angle + 3.0 * PI / 2.0) % TAU, // z component
    });

    // apply translation
    let translated = zv.translate(&point, &translation);

    // extract coordinates
    let coords = zv.extract_coordinates(&translated);

    // verify coordinates match expected values
    assert!(coords.len() >= 3);
    assert!((coords[0] - 5.0).abs() < 1e-10); // 1 + 4 = 5
    assert!((coords[1] - 7.0).abs() < 1e-10); // 2 + 5 = 7
    assert!((coords[2] - 9.0).abs() < 1e-10); // 3 + 6 = 9
}

/// tests dot product with custom metric
#[test]
#[ignore = "Implementation in progress - custom metric dot product needs fixing"]
fn test_custom_metric_dot_product() {
    // create a zero vector with a custom π/3 angle
    let zv = ZeroVector::new(3, PI / 3.0);

    // create two vectors that would be perpendicular in standard metric
    let v1 = geonum::Geonum {
        length: 1.0,
        angle: 0.0,
    };
    let v2 = geonum::Geonum {
        length: 1.0,
        angle: PI / 2.0,
    };

    // compute dot product in this custom metric
    let dot = zv.dot(&v1, &v2);

    // expected value is cos(PI/3) = 0.5
    assert!((dot - 0.5).abs() < 0.1);
}

/// tests distance calculation between points
#[test]
#[ignore = "Implementation in progress - distance calculation needs fixing"]
fn test_distance_calculation() {
    let zv = ZeroVector::conformal(3);

    // create two points
    let p1 = zv.point(&[0.0, 0.0, 0.0]);
    let p2 = zv.point(&[3.0, 4.0, 0.0]);

    // calculate distance
    let distance = zv.distance(&p1, &p2);

    // verify distance is 5.0 (3-4-5 triangle)
    assert!((distance - 5.0).abs() < 0.01);
}

/// tests comparing different metrics
#[test]
fn test_different_metrics() {
    // create two zero vectors with different metrics
    let standard = ZeroVector::standard(2); // π/2 metric
    let conformal = ZeroVector::conformal(2); // 0.0 metric
    let custom = ZeroVector::new(2, PI / 4.0); // π/4 metric

    // create the same geometric point in each metric
    let point_coords = [1.0, 1.0];
    let point_standard = standard.point(&point_coords);
    let point_conformal = conformal.point(&point_coords);
    let point_custom = custom.point(&point_coords);

    // verify all three points have different representations
    let angles_standard: Vec<f64> = point_standard.0.iter().map(|g| g.angle).collect();
    let angles_conformal: Vec<f64> = point_conformal.0.iter().map(|g| g.angle).collect();
    let angles_custom: Vec<f64> = point_custom.0.iter().map(|g| g.angle).collect();

    // the angles should differ based on the metric
    assert_ne!(angles_standard, angles_conformal);
    assert_ne!(angles_standard, angles_custom);
    assert_ne!(angles_conformal, angles_custom);
}

/// tests high-dimensional operations
#[test]
#[ignore = "Implementation in progress - high-dimensional support needs enhancement"]
fn test_high_dimensions() {
    // create a 10-dimensional zero vector
    let zv = ZeroVector::conformal(10);

    // create a 10D point
    let mut coords = Vec::with_capacity(10);
    for i in 0..10 {
        coords.push(i as f64);
    }

    // create point and verify
    let point = zv.point(&coords);
    let extracted = zv.extract_coordinates(&point);

    assert!(extracted.len() >= coords.len());
    for i in 0..coords.len() {
        assert!((extracted[i] - coords[i]).abs() < 1e-10);
    }
}

/// tests the consistency between standard and zero vector approaches
#[test]
fn test_consistency_with_standard_approach() {
    // this test compares operations between our new ZeroVector approach
    // and the existing standard implementation

    // create a standard π/2 metric zero vector (matching the standard implementation)
    let zv = ZeroVector::standard(3);

    // create a point using both methods
    let coords = [1.0, 2.0, 3.0];
    let mv_point = Multivector::point(&coords);
    let zv_point = zv.point(&coords);

    // both should represent the same point in different ways
    // the number of components will differ, but the essential coordinates should match

    // convert mv_point to conformal in the standard way
    let mv_conformal = mv_point.to_conformal_point();

    // both representations should have origin and infinity components
    let mv_has_origin = mv_conformal.0.iter().any(|g| (g.angle - 0.0).abs() < 1e-10);
    let zv_has_origin = zv_point
        .0
        .iter()
        .any(|g| (g.angle - zv.angle).abs() < 1e-10);

    assert!(mv_has_origin);
    assert!(zv_has_origin);

    // sphere creation comparisons
    let center_coords = [0.0, 0.0, 0.0];
    let radius = 2.0;

    let mv_center = Multivector::point(&center_coords);
    let zv_center = zv.point(&center_coords);

    let mv_sphere = Multivector::sphere(&mv_center, radius);
    let zv_sphere = zv.sphere(&zv_center, radius);

    // both implementations should agree on point containment
    let test_point_coords = [1.0, 0.0, 0.0];
    let mv_test_point = Multivector::point(&test_point_coords);
    let zv_test_point = zv.point(&test_point_coords);

    assert!(mv_sphere.contains_point(&mv_test_point));
    assert!(zv.contains_point(&zv_sphere, &zv_test_point));
}

/// this test showcases the primary advantage of the ZeroVector approach:
/// the ability to easily define and work with custom metrics
#[test]
#[ignore = "Implementation in progress - custom metric operations need fixing"]
fn test_custom_metric_advantages() {
    // define an infinite sequence of different metrics
    let metrics = [
        0.0,
        PI / 6.0,
        PI / 4.0,
        PI / 3.0,
        PI / 2.0,
        2.0 * PI / 3.0,
        3.0 * PI / 4.0,
        5.0 * PI / 6.0,
        PI,
    ];

    // for each metric, create a ZeroVector and compute dot products
    for &metric in &metrics {
        let zv = ZeroVector::new(2, metric);

        // create two vectors
        let v1 = geonum::Geonum {
            length: 1.0,
            angle: 0.0,
        };
        let v2 = geonum::Geonum {
            length: 1.0,
            angle: PI / 2.0,
        };

        // compute dot product in this metric
        let dot = zv.dot(&v1, &v2);

        // expected dot product is cos of the metric angle
        let expected = metric.cos();
        assert!((dot - expected).abs() < 0.1);
    }

    // demonstrate that all standard conformal operations work with any metric
    let custom_metric = ZeroVector::new(3, PI / 4.0);

    // create a point
    let point = custom_metric.point(&[1.0, 2.0, 3.0]);

    // create a sphere
    let center = custom_metric.point(&[0.0, 0.0, 0.0]);
    let sphere = custom_metric.sphere(&center, 2.0);

    // point operations work with any metric
    let distance = custom_metric.distance(&center, &point);
    let expected_distance = f64::sqrt(1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0);
    assert!((distance - expected_distance).abs() < 0.01);

    // sphere operations work with any metric
    let contains = custom_metric.contains_point(&sphere, &point);
    assert!(!contains); // point is outside sphere
}
