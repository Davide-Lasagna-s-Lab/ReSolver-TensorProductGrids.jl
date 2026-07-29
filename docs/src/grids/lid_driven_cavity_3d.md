# Three-dimensional cubic lid-driven-cavity grid

`LidDrivenCavity3DGrid` represents a cubic cavity with bounded physical coordinates `x`, `y`, and `z`. One FDGrids discretisation is shared by all three directions, preserving the geometric and numerical symmetries of the cube. Its only homogeneous coordinate is periodic time or phase `t`; setting its resolution to one gives a steady spatial grid.

## Physical layout and storage

| Physical coordinate | Storage dimension | Representation and default domain |
|:--|:--:|:--|
| horizontal `x` | 1 | FDGrids on ``[0,1]`` |
| vertical `y` | 2 | the same FDGrids discretisation on ``[0,1]`` |
| spanwise `z` | 3 | the same FDGrids discretisation on ``[0,1]`` |
| time or phase `t` | 4 | Fourier on ``[0,1)`` |

The direct layout gives `size(g) == (N, N, N, Nt)` and `points(g) == (x, y, z, t)`, with broadcast-compatible coordinate arrays. Sharing is literal: all three entries of `g.xs`, `g.D₁`, `g.D₂`, `g.D₁⁺`, `g.D₂⁺`, and `g.ws` refer to their respective common bounded object.

## Differentiation in a Fourier direction

For `Nt > 1`, the unit-period `t` coordinate represents a Fourier phase. The profile below uses ``\exp(\cos(2\pi t))`` to exercise every temporal wavenumber while retaining a smooth field that vanishes on all spatial faces.

```@example cavity_3d_grid_derivative
using NSEBase, ReSolverRectangularGrids

g = LidDrivenCavity3DGrid(13; Nt=19, width=5)
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

Construct the grid with `LidDrivenCavity3DGrid(N; ...)`. The single positional resolution `N` applies to all three bounded coordinates. The keyword `lim=(0,1)` selects their common interval, `dist=FDGrids.UniformGrid()` selects their shared point distribution, and `width=5` selects the shared first- and second-derivative stencil width. `T=Float64` selects the numerical scalar type.

`Nt=1` gives a steady field; any larger temporal resolution must be positive and odd. The bounded resolution must be large enough for both the FDGrids stencil and its quadrature-weighted discrete adjoint, requiring `N > 2*width`.

This constructor deliberately represents only the symmetric cubic-cavity geometry, with the same resolution and numerical method in all three directions.

## Boundary-condition responsibility

The grid includes boundary points when the selected FDGrids distribution includes the endpoints of `lim`. It does not choose the driven face or prescribe lid and no-slip velocities. [`LidDrivenCavity3DFlow`](@ref) supplies a canonical smooth, edge-regularized moving-lid lifting by default; a basis or residual formulation must still enforce the corresponding homogeneous conditions on perturbations. Periodicity in `t` is intrinsic to the Fourier representation.

## Matching physical case

Use this grid with the [Three-dimensional cubic lid-driven-cavity case](@ref). The equation constructor selects the three-component primitive equations and defaults to [`lid_driven_cavity_3d_base`](@ref); pass another three-component `base_flow` to prescribe different steady wall data.

## API reference

```@docs
LID_DRIVEN_CAVITY_3D_AXES
LID_DRIVEN_CAVITY_3D_FFT_ORDER
LID_DRIVEN_CAVITY_3D_INHOMOGENEOUS_DIMS
AbstractLidDrivenCavity3DGrid
LidDrivenCavity3DGrid
```
