```@raw html
<p style="text-align:center">
  <img src="assets/logo.svg" alt="ReSolverRectangularGrids.jl logo" style="width:min(100%, 760px)">
</p>
```

# ReSolverRectangularGrids.jl

`ReSolverRectangularGrids` provides FDGrids-backed rectangular domains and focused incompressible-flow case constructors for NSEBase. It is designed for wall-bounded simulations, invariant solutions, stability calculations, and adjoint-based analysis in which bounded coordinates use finite differences and periodic coordinates use Fourier series.

The package provides four layouts:

| Layout | Storage | Bounded coordinates | Fourier coordinates |
|:--|:--|:--|:--|
| [Channel](@ref) | `(y,x,z,t)` | `y` | `x,z,t` |
| [2D lid-driven cavity](@ref "Two-dimensional lid-driven-cavity grid") | `(x,y,t)` | `x,y` | `t` |
| [3D lid-driven cavity](@ref "Three-dimensional lid-driven-cavity grid") | `(x,y,z,t)` | `x,y,z` | `t` |
| [Square duct](@ref) | `(x,y,z,t)` | `x,y` | `z,t` |

Each constructor supplies collocation points, compact first- and second-derivative matrices, quadrature-consistent discrete adjoints, lazy product weights, Fourier scales, and resolution growth. The matching case constructors assemble NSEBase primitive-variable equations for plane Couette flow, plane Poiseuille flow, lid-driven cavities, and pressure-driven square ducts.

!!! warning "Boundary conditions stay explicit"
    A grid includes wall points and differential operators but does not impose velocity values. Encode inhomogeneous steady wall data in `base_flow` and require the perturbation basis or residual formulation to satisfy the corresponding homogeneous boundary conditions.

## Start here

- [Installation](@ref) explains the temporary development-branch setup required before compatible upstream releases are registered.
- [Quick start](@ref) builds a channel, verifies a Fourier derivative and an analytical norm, and constructs pressure-driven equations.
- [Coordinates and numerical conventions](@ref) defines storage order, wavenumber scaling, FFT normalization, quadrature, adjoints, dealiasing, and time or phase.
- [Worked examples](@ref) runs the same validation pattern for all four layouts.
- [Extending the package](@ref) constructs a new layout entirely through public FDGrids and NSEBase interfaces.

## Package boundaries

This package deliberately sits between two reusable numerical layers:

- [FDGrids.jl](https://github.com/Davide-Lasagna-s-Lab/FDGrids.jl) provides one-dimensional grids, compact differentiation matrices, quadrature, and weighted adjoints.
- [NSEBase.jl](https://github.com/Davide-Lasagna-s-Lab/NSEBase.jl) provides fields, transforms, physical-coordinate derivatives, norms, shifts, and incompressible equation formulations.

[ReSolver.jl](https://github.com/Davide-Lasagna-s-Lab/ReSolver.jl) and downstream flow packages can consume the concrete grids and equation objects without duplicating geometry-specific infrastructure.

## Public status

The repository is public for development and review. A tagged release is intentionally deferred until compatible NSEBase and FDGrids versions are released in the [organisation registry](https://github.com/Davide-Lasagna-s-Lab/Registry.jl). The current installation instructions pin the exact upstream branches exercised by continuous integration.
