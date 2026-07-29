# Two-dimensional lid-driven-cavity grid

`LidDrivenCavity2DGrid` represents a planar cavity with independent bounded discretisations in
physical `x` and `y`. A Fourier time or phase coordinate permits steady, periodic-orbit, and
frequency-domain formulations without adding a third spatial coordinate.

## Physical layout and storage

| Physical coordinate | Storage dimension | Representation and default domain |
|:--|:--:|:--|
| horizontal `x` | 1 | FDGrids on ``[0, 1]`` |
| vertical `y` | 2 | FDGrids on ``[0, 1]`` |
| spanwise `z` | absent | not represented by this two-dimensional grid |
| time or phase `t` | 3 | Fourier on ``[0, 1)`` |

The direct layout gives `size(g) == (Nx, Ny, Nt)` and `points(g) == (x, y, t)`, with each returned
coordinate shaped for broadcasting over a field. The `x` and `y` points, operators, and weights
are independent, so rectangular cavities and direction-specific discretisations are supported.

## Differentiation in a Fourier direction

For `Nt > 1`, `t` is a homogeneous unit-period coordinate. The following example differentiates a
field containing ``\exp(\cos(2\pi t))``, whose Fourier series exercises every temporal wavenumber.

```@example cavity_2d_grid_derivative
using NSEBase, ReSolverRectangularGrids

g = LidDrivenCavity2DGrid(13, 13; Nt=19, width=5)
u(x, y, t) = x * (1 - x) * y * (1 - y) * exp(cos(2π * t))
ut(x, y, t) = -2π * sin(2π * t) * u(x, y, t)

û = FFT(Field(g, u))
numerical = ddt!(FTField(g), û)
exact = FFT(Field(g, ut))
derivative_error = maximum(abs, parent(numerical) .- parent(exact)) /
                   maximum(abs, parent(exact))
derivative_error < 5e-7
```

## Constructor choices and defaults

Construct the grid with `LidDrivenCavity2DGrid(Nx, Ny; ...)`. `Nt=1` selects a steady field. The
shared interval keyword `lim=(0, 1)` supplies the default for both directions, while `xlim` and
`ylim` override it independently. Similarly, `dist=FDGrids.UniformGrid()` and `width=5` are shared
fallbacks; `xdist`, `ydist`, `xwidth`, and `ywidth` select direction-specific distributions and
stencils. Numerical data use `Float64` unless another real scalar type is passed as `T`.

The temporal resolution `Nt` must be positive and odd. Each bounded resolution must also be large
enough for its chosen FDGrids stencil and for construction of the quadrature-weighted discrete
adjoint; in particular, `Nx > 2*xwidth` and `Ny > 2*ywidth`.

## Boundary-condition responsibility

The grid includes boundary points when the selected FDGrids distributions include their interval
endpoints. It does not identify a driven wall or prescribe moving-lid and no-slip values. Supply a
base-flow lifting for inhomogeneous wall data and enforce homogeneous perturbation conditions in
the basis or residual formulation. Periodicity in `t` is intrinsic to the Fourier representation.

## Matching physical case

Use this grid with the [Two-dimensional lid-driven-cavity case](@ref). Its
`LidDrivenCavity2DFlow` constructor selects the two-component primitive equations and accepts the
base-flow lifting that carries the moving-lid data.

## API reference

```@docs
LID_DRIVEN_CAVITY_2D_AXES
LID_DRIVEN_CAVITY_2D_FFT_ORDER
LID_DRIVEN_CAVITY_2D_INHOMOGENEOUS_DIMS
AbstractLidDrivenCavity2DGrid
LidDrivenCavity2DGrid
```
