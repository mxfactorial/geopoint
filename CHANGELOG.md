# changelog

all notable changes to geopoint will be documented in this file

the format is based on [keep a changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [semantic versioning](https://semver.org/spec/v2.0.0.html)

## [unreleased]

## [0.1.1] - 2025-04-05

### added
- integrated `ZeroVector` implementation into the main library
- merged zero_vector.rs into lib.rs for unified implementation
- implemented trait-based design using zero vector principles internally
- comprehensive tests for zero vector functionality
- improved coordinate extraction with handling of special cases
- enhanced distance computation with fallback methods
- optimized translation operations for different object types
- fixed custom metric dot product implementation
- improved handling of high dimensional spaces
- updated documentation and examples with zero vector usage
- maintained O(1) complexity in all operations

### key achievements
- unified the traditional conformal geometry and zero vector designs
- all existing traits now use zero vector principles internally
- created consistent conversions between euclidean and conformal spaces
- implemented sphere creation and containment testing with custom metrics
- enhanced distance computations with improved numerical stability
- added meet and join operations as foundation for future work
- created extensible framework for custom metrics
- maintained backward compatibility with existing code
- comprehensive test coverage for all functionality
- clarified API design by removing artificial distinction between "basic" and "advanced" usage
- focused on zero vector design as the primary implementation approach

### improvements
- simplified codebase with unified design
- improved error handling and edge cases
- enhanced numerical stability in geometric operations
- better documentation with clear examples
- performance optimizations for common operations
- streamlined API for both designs
- enhanced support for custom metrics
- fixed test failures and improved test coverage
- moved examples to test suite for better integration
- simplified project structure by removing examples directory
- standardized code style throughout the codebase
- fixed clippy warnings for better code quality

## [0.1.0] - 2025-04-04

### added
- initial project setup
- integration with geonums O(1) representation
- core trait definitions and basic implementations
- README with usage examples and mathematical background
- test cases for primitive functionality