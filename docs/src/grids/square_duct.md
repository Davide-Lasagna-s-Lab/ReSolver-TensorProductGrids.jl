# Square duct

`SquareDuctGrid` represents a bounded square cross-section with periodic streamwise `z` and
periodic time or phase `t`. The two cross-stream directions deliberately share one FDGrids
discretisation, preserving the geometric and numerical symmetry of the square.

## Physical layout and storage

| Physical coordinate | Storage dimension | Representation and domain |
|:--|:--:|:--|
| cross-stream `x` | 1 | FDGrids on ``[0, 1]`` |
| cross-stream `y` | 2 | the same FDGrids discretisation on ``[0, 1]`` |
| streamwise `z` | 3 | Fourier on ``[0, 2\pi/\alpha)`` |
| time or phase `t` | 4 | Fourier on ``[0, 1)`` |

The direct layout gives `size(g) == (N, N, Nz, Nt)` and `points(g) == (x, y, z, t)`. Sharing is
literal in the stored bounded data: `g.xs[1] === g.xs[2]`, and the corresponding derivative
matrices, adjoints, and quadrature weights are also identical objects. Streamwise `z` is the
real-to-complex transform direction.

## Differentiation in a Fourier direction

The streamwise profile ``\exp(\cos(\alpha z))`` excites every Fourier wavenumber. This example
compares `ddz!` with its analytical derivative on a modest steady grid.

```@example square_duct_grid_derivative
using NSEBase, ReSolverRectangularGrids

α = 0.75
g = SquareDuctGrid(13, 19; α, width=5)
u(x, y, z, t) = x * (1 - x) * y * (1 - y) * exp(cos(α * z))
uz(x, y, z, t) = -α * sin(α * z) * u(x, y, z, t)

û = FFT(Field(g, u))
numerical = ddz!(FTField(g), û)
exact = FFT(Field(g, uz))
derivative_error = maximum(abs, parent(numerical) .- parent(exact)) /
                   maximum(abs, parent(exact))
derivative_error < 5e-7
```

## Constructor choices and defaults

Construct the grid with `SquareDuctGrid(N, Nz; ...)`, where one bounded resolution `N` is used for
both sides of the cross-section. `Nt=1` selects a steady field, `α=1` gives streamwise length
`2π`, `dist=FDGrids.GaussLobattoGrid()` selects the shared bounded distribution, `width=5` selects
the shared first- and second-derivative stencil width, and `T=Float64` selects the scalar type.

The Fourier resolutions `Nz` and `Nt` must be positive and odd. The bounded resolution must be
large enough for the FDGrids stencil and construction of its quadrature-weighted discrete
adjoints; in particular, `N > 2*width`. Use a general [`RectangularGrid`](@ref) when the two
cross-stream directions need different points, side lengths, or operators.

## Boundary-condition responsibility

The default Gauss--Lobatto distribution includes the four duct walls, but the grid does not impose
no-slip values or a streamwise pressure gradient. Enforce wall conditions in the basis, projection,
or residual formulation and select pressure forcing through the matching case constructor.
Periodicity in `z` and `t` is intrinsic to the Fourier representation.

## Matching physical case

Use this grid with [Square-duct flow](@ref). Its `SquareDuctFlow` constructor selects the
three-component primitive equations and applies a constant pressure-gradient force to the
streamwise velocity component.

## API reference

```@docs
SQUARE_DUCT_AXES
SQUARE_DUCT_FFT_ORDER
SQUARE_DUCT_INHOMOGENEOUS_DIMS
AbstractSquareDuctGrid
SquareDuctGrid
```
