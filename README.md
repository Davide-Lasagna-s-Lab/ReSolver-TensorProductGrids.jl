<p align="center">
  <img src="docs/src/assets/logo.svg" alt="ReSolverRectangularGrids.jl logo" width="720">
</p>

<p align="center">
  <a href="https://github.com/Davide-Lasagna-s-Lab/ReSolver-RectangularGrids.jl/actions/workflows/CI.yml"><img src="https://github.com/Davide-Lasagna-s-Lab/ReSolver-RectangularGrids.jl/actions/workflows/CI.yml/badge.svg?branch=main" alt="CI status"></a>
  <a href="https://Davide-Lasagna-s-Lab.github.io/ReSolver-RectangularGrids.jl/dev/"><img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="Development documentation"></a>
  <a href="https://codecov.io/gh/Davide-Lasagna-s-Lab/ReSolver-RectangularGrids.jl"><img src="https://codecov.io/gh/Davide-Lasagna-s-Lab/ReSolver-RectangularGrids.jl/branch/main/graph/badge.svg" alt="Code coverage"></a>
  <a href="https://julialang.org/"><img src="https://img.shields.io/badge/Julia-1.10%2B-9558B2.svg" alt="Julia 1.10 or later"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-389826.svg" alt="MIT license"></a>
</p>

# ReSolverRectangularGrids.jl

`ReSolverRectangularGrids.jl` provides concrete, FDGrids-backed Cartesian grids and ready-to-use incompressible-flow case constructors for [NSEBase.jl](https://github.com/Davide-Lasagna-s-Lab/NSEBase.jl). It brings the common numerical machinery for channels, lid-driven cavities, and square ducts into one small package: coordinates, compact finite-difference operators, Fourier directions, tensor-product quadrature, type-stable derivatives, discrete adjoints, base-flow liftings, and physical forcing policies.

The grid describes geometry and differentiation. Boundary conditions remain explicit in the basis, projection, or residual formulation, which makes the same grid useful for different wall models and perturbation spaces.

## Included layouts

| Layout | Constructor | Default physical domain | Storage order | Bounded directions | Fourier directions |
|:--|:--|:--|:--|:--|:--|
| Channel | `ChannelGrid(Nx, Ny, Nz)` | `y ∈ [-1,1]`; periodic `x,z` | `(y,x,z,t)` | `y` | `x,z,t` |
| 2D square lid-driven cavity | `LidDrivenCavity2DGrid(N)` | `[0,1]²` | `(x,y,t)` | `x,y` | `t` |
| 3D cubic lid-driven cavity | `LidDrivenCavity3DGrid(N)` | `[0,1]³` | `(x,y,z,t)` | `x,y,z` | `t` |
| Square duct | `SquareDuctGrid(N, Nz)` | `[0,1]²`; periodic `z` | `(x,y,z,t)` | `x,y` | `z,t` |

All Fourier resolutions are positive and odd. `Nt=1` represents a steady field, while larger odd `Nt` values represent a ``2π``-periodic phase coordinate.

## Installation

The package is currently published as a development repository. Compatible NSEBase and FDGrids changes have not yet been released in the organisation registry, so install the tested upstream branches explicitly:

```julia
using Pkg

Pkg.add([
    PackageSpec(url="https://github.com/Davide-Lasagna-s-Lab/NSEBase.jl", rev="dev"),
    PackageSpec(url="https://github.com/Davide-Lasagna-s-Lab/FDGrids.jl", rev="master"),
])
Pkg.develop(url="https://github.com/Davide-Lasagna-s-Lab/ReSolver-RectangularGrids.jl")
```

Julia 1.10 or later is required. A normal registered installation and a tagged `v0.1.0` release will follow compatible NSEBase and FDGrids releases; see the [installation guide](docs/src/installation.md) for the rationale and update path.

## Quick start

This example builds a Fourier–finite-difference channel, differentiates a field in physical `z`, evaluates its weighted norm, and constructs pressure-driven equations at friction Reynolds number `Reτ`:

```julia
using LinearAlgebra
using NSEBase
using ReSolverRectangularGrids

g = ChannelGrid(19, 25, 19; Nt=1, α=0.5, β=1, width=7)
y, x, z, t = points(g) # coordinates are returned in storage order

u = Field(g, (y, x, z, t) -> (1 - y^2) * exp(cos(0.5x)) * exp(cos(z)))
û = FFT(u)
uz = ddz!(FTField(g), û)
energy_norm = norm(û)

Reτ = 180
U = (Reτ / 2) .* plane_poiseuille_base(g)
equations = PlanePoiseuilleFlow(g, Reτ; base_flow=(U, nothing, nothing), f=1, fftw_flags=FFTW.ESTIMATE)
```

The `exp(cos(·))` factors have nonzero coefficients at every Fourier wavenumber, making this a more revealing spectral example than a single trigonometric mode. The complete [getting-started guide](docs/src/quickstart.md) checks the derivative against its analytical value and explains every object above.

## Grids, cases, and boundary conditions

Files under `src/grids/` define numerical layouts. Files under `src/cases/` define base profiles and factories for the primitive incompressible equations:

- `PlaneCouetteFlow` and `PlanePoiseuilleFlow` use `ChannelGrid`;
- `LidDrivenCavity2DFlow` and `LidDrivenCavity3DFlow` provide canonical smooth moving-lid liftings and accept explicit replacements through `base_flow`;
- `SquareDuctFlow` applies a required, nonzero pressure-gradient force in physical streamwise `z`.

The constructors return NSEBase `ProjectedNSE` operators, but they do not silently impose wall values. For an inhomogeneous wall condition, put the prescribed steady contribution in `base_flow` and make the perturbation basis satisfy homogeneous boundary conditions. The manual explains [coordinates and numerical conventions](docs/src/conventions.md), including the precise quadrature-weighted discrete adjoint.

## Worked examples

The four standalone examples have the same structure and are executed by the test suite:

- [channel flow](examples/channel.jl)
- [2D square lid-driven cavity](examples/lid_driven_cavity_2d.jl)
- [3D cubic lid-driven cavity](examples/lid_driven_cavity_3d.jl)
- [square-duct flow](examples/square_duct.jl)

Each example uses an analytical field, checks a derivative and norm, and constructs the relevant equations. See the [worked-examples guide](docs/src/examples.md) for the numerical formulas and expected accuracy.

## Documentation

The [development manual](https://Davide-Lasagna-s-Lab.github.io/ReSolver-RectangularGrids.jl/dev/) covers installation, storage layouts, Fourier conventions, product quadrature, continuous and discrete adjoints, resolution growth, all bundled grids and cases, and extension of the package through public FDGrids and NSEBase interfaces.

Useful local entry points are the [quick start](docs/src/quickstart.md), [rectangular-grid design](docs/src/rectangular.md), [extension guide](docs/src/extending.md), and [API reference](docs/src/api.md).

## ReSolver ecosystem

- [ReSolver.jl](https://github.com/Davide-Lasagna-s-Lab/ReSolver.jl) provides the surrounding reduced-order and flow-analysis ecosystem.
- [NSEBase.jl](https://github.com/Davide-Lasagna-s-Lab/NSEBase.jl) defines fields, transforms, differential operators, and incompressible equation formulations.
- [FDGrids.jl](https://github.com/Davide-Lasagna-s-Lab/FDGrids.jl) supplies compact finite-difference matrices, quadrature, and weighted adjoints.
- [Davide-Lasagna-s-Lab/Registry.jl](https://github.com/Davide-Lasagna-s-Lab/Registry.jl) is the organisation's Julia package registry.

## Development

Run the full local verification with:

```julia
using Pkg
Pkg.test()
```

ReSolverRectangularGrids.jl is available under the [MIT License](LICENSE).
