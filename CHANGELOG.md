# changelog

all notable changes to geopoint will be documented in this file

the format is based on [keep a changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [semantic versioning](https://semver.org/spec/v2.0.0.html)

## [unreleased]

### added
- initial trait architecture for conformal geometry
- implementation of `ConformalGeometry` trait with conversion methods
- basic geometric primitives via `GeometricPrimitives` trait
- conformal primitives via `ConformalPrimitives` trait
- geometric operations for point containment and distance measurement
- basic transformation operations (translation, rotation, dilation)
- versor implementations for transformations
- unit tests for core functionality

## [0.1.0] - 2025-04-04

### added
- initial project setup
- integration with geonums O(1) representation
- core trait definitions and basic implementations
- README with usage examples and mathematical background
- test cases for primitive functionality