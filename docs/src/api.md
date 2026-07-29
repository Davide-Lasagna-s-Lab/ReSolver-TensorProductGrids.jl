# API reference

The manual renders each public docstring beside the concepts and examples that give it context. This page provides a compact map of that interface.

## Grid core

- [`RectangularGrid`](@ref) — validated Cartesian product grid and low-level extension point.
- [`RectangularProductWeights`](@ref) — lazy one-, two-, or three-dimensional product quadrature.

`RectangularGrid` implements the NSEBase interfaces `points`, `weights`, `wavenumber_scale`, `growto`, and `derivative_matrix`, together with `size`, `eltype`, and `convert`.

## Layout contracts and constructors

| Layout | Abstract contract | Concrete constructor |
|:--|:--|:--|
| Channel | [`AbstractChannelGrid`](@ref) | [`ChannelGrid`](@ref) |
| 2D square cavity | [`AbstractLidDrivenCavity2DGrid`](@ref) | [`LidDrivenCavity2DGrid`](@ref) |
| 3D cubic cavity | [`AbstractLidDrivenCavity3DGrid`](@ref) | [`LidDrivenCavity3DGrid`](@ref) |
| Square duct | [`AbstractSquareDuctGrid`](@ref) | [`SquareDuctGrid`](@ref) |

Each layout page also renders its exported `AXES`, Fourier-order, and inhomogeneous-dimension constants.

## Cases and forcings

- [`PlaneCouetteFlow`](@ref), [`PlanePoiseuilleFlow`](@ref), [`plane_couette_base`](@ref), and [`plane_poiseuille_base`](@ref).
- [`LidDrivenCavity2DFlow`](@ref), [`LidDrivenCavity3DFlow`](@ref), [`lid_driven_cavity_2d_base`](@ref), and [`lid_driven_cavity_3d_base`](@ref).
- [`SquareDuctFlow`](@ref).
- [`CoriolisForce`](@ref) and [`ConstantBodyForce`](@ref).

## Alphabetical index

```@index
Modules = [ReSolverRectangularGrids]
Pages = ["rectangular.md", "grids/channel.md", "grids/lid_driven_cavity_2d.md",
         "grids/lid_driven_cavity_3d.md", "grids/square_duct.md", "forcings.md",
         "cases/channel.md", "cases/lid_driven_cavity_2d.md",
         "cases/lid_driven_cavity_3d.md", "cases/square_duct.md"]
```
