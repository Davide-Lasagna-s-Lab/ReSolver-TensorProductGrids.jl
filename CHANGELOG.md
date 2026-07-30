# Changelog

All notable changes to ReSolverRectangularGrids.jl will be documented in this file. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and tagged releases will follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- A concrete `RectangularGrid` for one to three FDGrids directions combined with periodic Fourier coordinates.
- Lazy `RectangularProductWeights` for allocation-free tensor-product quadrature.
- Channel, square two-dimensional cavity, cubic three-dimensional cavity, and square-duct grid constructors; every square or cubic bounded layout reuses one FDGrids discretisation in all symmetric directions.
- Plane Couette, plane Poiseuille, lid-driven-cavity, and pressure-driven square-duct equation factories, including canonical smooth divergence-free moving-lid liftings.
- Reusable Coriolis and nonzero constant-body-force policies, with explicit rejection of inert zero
  body forces.
- Analytical derivative, Laplacian, norm, shift, product-quadrature, and discrete-adjoint tests for every layout.
- Symmetric executable examples, a full Documenter manual, public-repository metadata, and continuous integration.

### Fixed

- Correct Documenter markup for browser-rendered display and inline mathematics.

### Publication status

- The source repository is available for development and review.
- The first tagged release is deferred until compatible NSEBase and FDGrids versions are released and registered.

[Unreleased]: https://github.com/Davide-Lasagna-s-Lab/ReSolver-RectangularGrids.jl/commits/main
