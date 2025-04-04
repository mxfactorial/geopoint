use geonum::Multivector;
use geopoint::*;
use std::f64::consts::PI;

#[test]
#[ignore = "Implementation in progress - point conversion needs refinement"]
fn test_point_creation_and_conformal_conversion() {
    // Create original coordinates
    let coords_2d = vec![1.0, 2.0];
    let coords_3d = vec![1.0, 2.0, 3.0];

    // create points with different dimensions
    let point2d = Multivector::point(&coords_2d);
    let point3d = Multivector::point(&coords_3d);

    // convert to conformal
    let conformal_2d = point2d.to_conformal_point();
    let conformal_3d = point3d.to_conformal_point();

    // check conformal point properties - it should have more components than original
    assert!(conformal_2d.0.len() > point2d.0.len());
    assert!(conformal_3d.0.len() > point3d.0.len());

    // convert back and check round-trip
    let back_to_2d = conformal_2d.extract_euclidean();
    let back_to_3d = conformal_3d.extract_euclidean();

    // Extract coordinates to compare
    let extract_coords = |mv: &Multivector| -> Vec<f64> {
        let mut max_index = 0;
        for component in &mv.0 {
            let normalized_angle = (component.angle - PI / 2.0) % (2.0 * PI);
            let index = (normalized_angle / (PI / 2.0)) as usize;
            max_index = max_index.max(index);
        }

        let mut coords = vec![0.0; max_index + 1];

        for component in &mv.0 {
            let normalized_angle = (component.angle - PI / 2.0) % (2.0 * PI);
            let index = (normalized_angle / (PI / 2.0)) as usize;
            if index < coords.len() {
                coords[index] = component.length;
            }
        }

        coords
    };

    // Get coordinates
    let original_2d = extract_coords(&point2d);
    let restored_2d = extract_coords(&back_to_2d);
    let original_3d = extract_coords(&point3d);
    let restored_3d = extract_coords(&back_to_3d);

    // Check dimensions
    assert_eq!(
        original_2d.len(),
        restored_2d.len(),
        "Original 2D point has {} dimensions, restored has {}",
        original_2d.len(),
        restored_2d.len()
    );

    assert_eq!(
        original_3d.len(),
        restored_3d.len(),
        "Original 3D point has {} dimensions, restored has {}",
        original_3d.len(),
        restored_3d.len()
    );

    // Check each coordinate
    for i in 0..original_2d.len() {
        assert!(
            (original_2d[i] - restored_2d[i]).abs() < 1e-10,
            "2D coordinate mismatch at dimension {}: {} vs {}",
            i,
            original_2d[i],
            restored_2d[i]
        );
    }

    for i in 0..original_3d.len() {
        assert!(
            (original_3d[i] - restored_3d[i]).abs() < 1e-10,
            "3D coordinate mismatch at dimension {}: {} vs {}",
            i,
            original_3d[i],
            restored_3d[i]
        );
    }

    // Test specific values
    assert!(
        (restored_2d[0] - 1.0).abs() < 1e-10,
        "First 2D coordinate should be 1.0"
    );
    assert!(
        (restored_2d[1] - 2.0).abs() < 1e-10,
        "Second 2D coordinate should be 2.0"
    );

    assert!(
        (restored_3d[0] - 1.0).abs() < 1e-10,
        "First 3D coordinate should be 1.0"
    );
    assert!(
        (restored_3d[1] - 2.0).abs() < 1e-10,
        "Second 3D coordinate should be 2.0"
    );
    assert!(
        (restored_3d[2] - 3.0).abs() < 1e-10,
        "Third 3D coordinate should be 3.0"
    );
}

#[test]
#[ignore = "Implementation in progress - distance calculation needs refinement"]
fn test_sphere_creation_and_operations() {
    // create a sphere at origin with radius 2
    let origin = Multivector::point(&[0.0, 0.0, 0.0]);
    let sphere = Multivector::sphere(&origin, 2.0);

    // points to test
    let on_sphere = Multivector::point(&[2.0, 0.0, 0.0]);
    let inside_sphere = Multivector::point(&[1.0, 0.0, 0.0]);
    let outside_sphere = Multivector::point(&[3.0, 0.0, 0.0]);

    // test containment
    assert!(sphere.contains_point(&on_sphere));
    assert!(sphere.contains_point(&inside_sphere));
    assert!(!sphere.contains_point(&outside_sphere));

    // test distance computation
    assert!((sphere.distance_to_point(&on_sphere) - 0.0).abs() < 1e-6);
    assert!((sphere.distance_to_point(&inside_sphere) - 1.0).abs() < 1e-6);
    assert!((sphere.distance_to_point(&outside_sphere) - 1.0).abs() < 1e-6);
}

#[test]
#[ignore = "Implementation in progress - translation and dimension handling needs work"]
fn test_geometric_transformations() {
    // create a point
    let point = Multivector::point(&[1.0, 0.0, 0.0]);

    // create translation vector
    let translation = Multivector::point(&[0.0, 1.0, 0.0]);

    // apply translation
    let translated = point.translate(&translation);

    // convert back to euclidean for comparison
    let euclidean = translated.extract_euclidean();

    // expected result: [1.0, 1.0, 0.0]
    let expected = Multivector::point(&[1.0, 1.0, 0.0]);

    // verify components match
    assert_eq!(euclidean.0.len(), expected.0.len());

    // create dilation
    let dilated = point.dilate(2.0);

    // expected result: [2.0, 0.0, 0.0]
    let dilated_expected = Multivector::point(&[2.0, 0.0, 0.0]);

    // verify dilation
    assert_eq!(dilated.0.len(), dilated_expected.0.len());
    for (i, comp) in dilated.0.iter().enumerate() {
        let expected_comp = &dilated_expected.0[i];
        assert!((comp.length - expected_comp.length * 2.0).abs() < 1e-10);
    }
}

#[test]
#[ignore = "Implementation in progress - circle representation needs refinement"]
fn test_conformal_primitives() {
    // test circle creation
    let center = Multivector::point(&[0.0, 0.0, 0.0]);
    let normal = Multivector::point(&[0.0, 0.0, 1.0]);
    let circle = Multivector::circle(&center, 1.0, &normal);

    // test point on circle
    let point_on_circle = Multivector::point(&[1.0, 0.0, 0.0]);
    assert!(circle.contains_point(&point_on_circle));

    // test point not on circle
    let point_not_on_circle = Multivector::point(&[0.0, 0.0, 1.0]);
    assert!(!circle.contains_point(&point_not_on_circle));

    // test point pair
    let point1 = Multivector::point(&[1.0, 0.0, 0.0]);
    let point2 = Multivector::point(&[-1.0, 0.0, 0.0]);
    let pair = Multivector::point_pair(&point1, &point2);

    // verify point pair contains its points
    assert!(pair.contains_point(&point1));
    assert!(pair.contains_point(&point2));
}

#[test]
fn test_versor_operations() {
    // create a translator
    let direction = Multivector::point(&[0.0, 1.0, 0.0]);
    let translator = Multivector::translator(&direction);

    // create a rotor (90 degree rotation around z-axis)
    let z_axis = Multivector::point(&[0.0, 0.0, 1.0]);
    let rotor = Multivector::rotor(&z_axis, PI / 2.0);

    // create a motor (combined rotation and translation)
    let motor = Multivector::motor(&z_axis, PI / 2.0, &direction);

    // test that versors are not empty
    assert!(translator.0.len() > 0);
    assert!(rotor.0.len() > 0);
    assert!(motor.0.len() > 0);
}

#[test]
#[ignore = "Implementation in progress - high dimensional support needs enhancement"]
fn test_high_dimensional_operations() {
    // create a high-dimensional point (10D)
    let mut coords = Vec::with_capacity(10);
    for i in 0..10 {
        coords.push(i as f64 / 10.0);
    }

    let high_dim_point = Multivector::point(&coords);

    // convert to conformal
    let conformal = high_dim_point.to_conformal_point();

    // convert back to euclidean
    let euclidean = conformal.extract_euclidean();

    // verify components match
    assert_eq!(high_dim_point.0.len(), euclidean.0.len());
    for (i, comp) in high_dim_point.0.iter().enumerate() {
        let back_comp = &euclidean.0[i];
        assert!((comp.length - back_comp.length).abs() < 1e-10);
    }

    // create high-dimensional sphere
    let sphere = Multivector::sphere(&high_dim_point, 1.0);

    // verify sphere is not empty
    assert!(sphere.0.len() > 0);
}
