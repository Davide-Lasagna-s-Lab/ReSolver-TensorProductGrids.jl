# Two-dimensional lid-driven-cavity grid

`LidDrivenCavity2DGrid` represents a square cavity with bounded physical coordinates `x` and `y`. One FDGrids discretisation is shared by both directions, preserving the geometric and numerical symmetry of the square. A Fourier time or phase coordinate permits steady, periodic-orbit, and frequency-domain formulations without adding a third spatial coordinate.

## Physical layout and storage

| Physical coordinate | Storage dimension | Representation and default domain |
|:--|:--:|:--|
| horizontal `x` | 1 | FDGrids on ``[0,1]`` |
| vertical `y` | 2 | the same FDGrids discretisation on ``[0,1]`` |
| spanwise `z` | absent | not represented by this two-dimensional grid |
| time or phase `t` | 3 | Fourier on ``[0,1)`` |

The direct layout gives `size(g) == (N, N, Nt)` and `points(g) == (x, y, t)`, with each coordinate shaped for broadcasting over a field. Sharing is literal in the stored bounded data: `g.xs[1] === g.xs[2]`, and the corresponding derivative matrices, adjoints, and quadrature weights are the same objects.

## Differentiation in a Fourier direction

For `Nt > 1`, `t` is a homogeneous unit-period coordinate. The following example differentiates a field containing ``\exp(\cos(2\pi t))``, whose Fourier series exercises every temporal wavenumber.

```@example cavity_2d_grid_derivative
using NSEBase, ReSolverRectangularGrids

g = LidDrivenCavity2DGrid(13; Nt=19, width=5)
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

Construct the grid with `LidDrivenCavity2DGrid(N; ...)`. The single positional resolution `N` applies to both bounded coordinates. The keyword `lim=(0,1)` selects their common interval, `dist=FDGrids.UniformGrid()` selects their shared point distribution, and `width=5` selects the shared first- and second-derivative stencil width. `T=Float64` selects the numerical scalar type.

`Nt=1` gives a steady field; any larger temporal resolution must be positive and odd. The bounded resolution must be large enough for both the FDGrids stencil and its quadrature-weighted discrete adjoint, requiring `N > 2*width`.

This constructor deliberately represents only the symmetric square-cavity geometry, with the same resolution and numerical method in both directions.

## Boundary-condition responsibility

The grid includes boundary points when the selected FDGrids distribution includes the endpoints of `lim`. It does not prescribe moving-lid or no-slip values. [`LidDrivenCavity2DFlow`](@ref) supplies a canonical smooth moving-lid lifting by default; a basis or residual formulation must still enforce the corresponding homogeneous conditions on perturbations. Periodicity in `t` is intrinsic to the Fourier representation.

## Matching physical case

Use this grid with the [Two-dimensional lid-driven-cavity case](@ref). The equation constructor selects the two-component primitive equations and defaults to [`lid_driven_cavity_2d_base`](@ref); pass another two-component `base_flow` to prescribe different steady wall data.

## API reference

```@docs
LID_DRIVEN_CAVITY_2D_AXES
LID_DRIVEN_CAVITY_2D_FFT_ORDER
LID_DRIVEN_CAVITY_2D_INHOMOGENEOUS_DIMS
AbstractLidDrivenCavity2DGrid
LidDrivenCavity2DGrid
```
