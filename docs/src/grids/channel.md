# Channel

`ChannelGrid` represents a plane channel with one bounded wall-normal direction and three Fourier
directions. Constructor resolutions are supplied in physical order `(Nx, Ny, Nz)`, but arrays are
stored with wall-normal `y` first so that each finite-difference operation acts on contiguous data.

## Physical layout and storage

| Physical coordinate | Storage dimension | Representation and domain |
|:--|:--:|:--|
| streamwise `x` | 2 | Fourier on ``[0, 2\pi/\alpha)`` |
| wall-normal `y` | 1 | FDGrids on the fixed interval ``[-1, 1]`` |
| spanwise `z` | 3 | Fourier on ``[0, 2\pi/\beta)`` |
| time or phase `t` | 4 | Fourier on ``[0,2\pi)`` |

Thus `size(g) == (Ny, Nx, Nz, Nt)` and `points(g)` returns `(y, x, z, t)` in storage order.
Streamwise `x` is the real-to-complex transform direction; `z` and `t` use complex transforms.
The half-channel height is the length scale because the walls are fixed at ``y=-1`` and ``y=1``.

## Differentiation in a Fourier direction

The function ``\exp(\cos(\alpha x))`` contains every Fourier harmonic, making it a more demanding
check than a single trigonometric mode. Here `ddx!` differentiates with respect to physical
streamwise `x`, even though `x` is the second storage dimension.

```@example channel_grid_derivative
using NSEBase, ReSolverRectangularGrids

α = 0.5
g = ChannelGrid(19, 13, 1; α, width=5)
u(y, x, z, t) = (1 - y^2) * exp(cos(α * x))
ux(y, x, z, t) = -α * sin(α * x) * u(y, x, z, t)

û = FFT(Field(g, u))
numerical = ddx!(FTField(g), û)
exact = FFT(Field(g, ux))
derivative_error = maximum(abs, parent(numerical) .- parent(exact)) /
                   maximum(abs, parent(exact))
derivative_error < 5e-7
```

## Constructor choices and defaults

Two constructor paths are available:

1. `ChannelGrid(Nx, Ny, Nz; ...)` creates the wall-normal points, first and second derivative
   matrices, quadrature-weighted discrete adjoints, and quadrature weights with FDGrids.
2. `ChannelGrid(y, Nx, Nz, Nt, α, β, Dy, Dy2, Dya, Dy2a, wy, T)` accepts a complete precomputed
   wall-normal discretisation. This path preserves correctly typed input objects and uses the
   supplied adjoints without recomputing them.

For the FDGrids-backed constructor, `Nt=1` selects a steady field, `α=β=1`, the default wall-normal
distribution is `FDGrids.GaussLobattoGrid()`, the default stencil width is `5`, and the default
scalar type is `Float64`. Fourier resolutions `Nx`, `Nz`, and `Nt` must be positive and odd. The
precomputed-data constructor requires `y` to span `[-1, 1]`; its derivative matrices, adjoints,
and weights must all be compatible with `length(y)`.

## Boundary-condition responsibility

The default Gauss--Lobatto distribution includes both walls, but a `ChannelGrid` only supplies
coordinates, quadrature, and differentiation operators. It does not impose no-slip, moving-wall,
fixed-flux, or pressure-gradient conditions. Encode velocity boundary conditions in the basis,
projection, or residual formulation; choose forcing and base-flow data through the matching case
constructor. Periodicity in `x`, `z`, and `t` is intrinsic to the Fourier representation.

## Matching physical case

Use this grid with [Channel flows](@ref), which provides the `PlaneCouetteFlow` and
`PlanePoiseuilleFlow` equation constructors together with their canonical base-flow profiles.

## API reference

```@docs
CHANNEL_AXES
CHANNEL_FFT_ORDER
CHANNEL_INHOMOGENEOUS_DIMS
AbstractChannelGrid
ChannelGrid
```
