# Three-dimensional lid-driven-cavity grid

`LidDrivenCavity3DGrid` represents a fully bounded Cartesian cavity with independent FDGrids
discretisations in physical `x`, `y`, and `z`. Its only homogeneous coordinate is periodic time or
phase `t`; setting its resolution to one gives a steady spatial grid.

## Physical layout and storage

| Physical coordinate | Storage dimension | Representation and default domain |
|:--|:--:|:--|
| horizontal `x` | 1 | FDGrids on ``[0, 1]`` |
| vertical `y` | 2 | FDGrids on ``[0, 1]`` |
| spanwise `z` | 3 | FDGrids on ``[0, 1]`` |
| time or phase `t` | 4 | Fourier on ``[0, 1)`` |

The direct layout gives `size(g) == (Nx, Ny, Nz, Nt)` and `points(g) == (x, y, z, t)`, with
broadcast-compatible coordinate arrays. All three bounded directions have their own points,
operators, adjoints, and quadrature weights, allowing cuboidal domains and anisotropic resolution.

## Differentiation in a Fourier direction

For `Nt > 1`, the unit-period `t` coordinate represents a Fourier phase. The profile below uses
``\exp(\cos(2\pi t))`` to exercise every temporal wavenumber while retaining a smooth field that
vanishes on all spatial faces.

```@example cavity_3d_grid_derivative
using NSEBase, ReSolverRectangularGrids

g = LidDrivenCavity3DGrid(13, 13, 13; Nt=19, width=5)
u(x, y, z, t) = x * (1 - x) * y * (1 - y) * z * (1 - z) * exp(cos(2π * t))
ut(x, y, z, t) = -2π * sin(2π * t) * u(x, y, z, t)

û = FFT(Field(g, u))
numerical = ddt!(FTField(g), û)
exact = FFT(Field(g, ut))
derivative_error = maximum(abs, parent(numerical) .- parent(exact)) /
                   maximum(abs, parent(exact))
derivative_error < 5e-7
```

## Constructor choices and defaults

Construct the grid with `LidDrivenCavity3DGrid(Nx, Ny, Nz; ...)`. `Nt=1` selects a steady field.
The shared `lim=(0, 1)`, `dist=FDGrids.UniformGrid()`, and `width=5` keywords provide defaults for
all bounded directions. Override them independently with `xlim`, `ylim`, `zlim`, `xdist`, `ydist`,
`zdist`, `xwidth`, `ywidth`, and `zwidth`. Numerical data use `Float64` unless another real scalar
type is passed as `T`.

The temporal resolution `Nt` must be positive and odd. Each bounded resolution must be large
enough for its selected stencil and quadrature-weighted discrete adjoint; each must exceed twice
the corresponding stencil width.

## Boundary-condition responsibility

The grid includes boundary points when the selected FDGrids distributions include their interval
endpoints. It does not choose the driven face or prescribe lid and no-slip velocities. Supply the
inhomogeneous wall data through a base-flow lifting and enforce homogeneous perturbation
conditions in the basis or residual formulation. Periodicity in `t` is intrinsic to the Fourier
representation.

## Matching physical case

Use this grid with the [Three-dimensional lid-driven-cavity case](@ref). Its
`LidDrivenCavity3DFlow` constructor selects the three-component primitive equations and accepts a
three-component base-flow lifting.

## API reference

```@docs
LID_DRIVEN_CAVITY_3D_AXES
LID_DRIVEN_CAVITY_3D_FFT_ORDER
LID_DRIVEN_CAVITY_3D_INHOMOGENEOUS_DIMS
AbstractLidDrivenCavity3DGrid
LidDrivenCavity3DGrid
```
