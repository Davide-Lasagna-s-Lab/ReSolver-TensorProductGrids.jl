# =============================================================================================== #
# Three-dimensional lid-driven-cavity grid                                                       #
# =============================================================================================== #
#
# This file defines a fully bounded cubic three-dimensional cavity. Physical `x`, `y`, and `z`
# share one FDGrids discretisation, while `t` is the only Fourier time or phase direction. The grid
# supplies wall points, quadrature, and derivatives without prescribing which wall is driven.
#
# Typical use; the positional resolution `N` is shared by `x`, `y`, and `z`, while the optional
# temporal resolution is the keyword `Nt`:
#
#     g = LidDrivenCavity3DGrid(49; Nt=1, width=7)
#     x, y, z, t = NSEBase.points(g)  # coordinates are returned in storage order
#
# Downstream algorithms can dispatch on `AbstractLidDrivenCavity3DGrid` without depending on this
# concrete FDGrids-backed representation.

# =============================================================================================== #
# Layout constants                                                                               #
# =============================================================================================== #

"""
    LID_DRIVEN_CAVITY_3D_AXES = (1, 2, 3, 4)

Map physical `(x, y, z, t)` directly to 3D-cavity storage order `(x, y, z, t)`. Bounded `x`, `y`,
and `z` occupy the first three dimensions, followed by time or phase `t`.
"""
const LID_DRIVEN_CAVITY_3D_AXES = (1, 2, 3, 4)

"""
    LID_DRIVEN_CAVITY_3D_FFT_ORDER = (4,)

The 3D cavity's only Fourier dimension, corresponding to unit-period time or phase `t`.
"""
const LID_DRIVEN_CAVITY_3D_FFT_ORDER = (4,)

"""
    LID_DRIVEN_CAVITY_3D_INHOMOGENEOUS_DIMS = (1, 2, 3)

The three FDGrids dimensions, corresponding to bounded physical coordinates `(x, y, z)`.
"""
const LID_DRIVEN_CAVITY_3D_INHOMOGENEOUS_DIMS = (1, 2, 3)

# =============================================================================================== #
# Abstract layout contract                                                                       #
# =============================================================================================== #

"""
    AbstractLidDrivenCavity3DGrid{T}

Layout contract for any real scalar type `T` stored as `(x, y, z, t)` with bounded `x`, `y`, and
`z` and a Fourier time or phase direction `t`.

Equal side lengths and shared numerical objects are guarantees of the concrete constructor, not
requirements imposed by this abstract layout contract.

The alias describes a layout rather than one concrete representation. NSEBase's device and
domain-decomposition wrappers preserve the parent's `AbstractGrid` layout parameters, so wrapped
3D-cavity grids still satisfy `AbstractLidDrivenCavity3DGrid`. Downstream algorithms can therefore
use one dispatch for serial, GPU-backed, and decomposed grids with this layout.
"""
const AbstractLidDrivenCavity3DGrid{T} =
    AbstractGrid{T, 4, LID_DRIVEN_CAVITY_3D_AXES, LID_DRIVEN_CAVITY_3D_FFT_ORDER}

# =============================================================================================== #
# Concrete rectangular-grid alias                                                                #
# =============================================================================================== #

"""
    LidDrivenCavity3DGrid

Concrete [`RectangularGrid{3}`](@ref) alias for the 3D cavity layout. Its stored size is
`(N, N, N, Nt)` and its only Fourier scale is the unit-period time or phase scale ``2π``. The
constructor shares its points, derivative matrices, adjoints, and weights among `x`, `y`, and `z`.
"""
const LidDrivenCavity3DGrid{T, S, XS, D1, D2, A1, A2, WS, WP} = RectangularGrid{
    3, T, S, 4, LID_DRIVEN_CAVITY_3D_AXES, LID_DRIVEN_CAVITY_3D_FFT_ORDER,
    XS, D1, D2, A1, A2, WS, WP, 1,
}

# =============================================================================================== #
# FDGrids-backed constructor                                                                     #
# =============================================================================================== #

"""
    LidDrivenCavity3DGrid(N; Nt=1, lim=(0, 1),
                          dist=FDGrids.UniformGrid(), width=5, T=Float64)

Construct a fully bounded cubic 3D cavity stored as `(x, y, z, t)`. All spatial directions use `N`
collocation points on `lim`. One FDGrids discretisation is reused in `x`, `y`, and `z`, so
corresponding entries in `xs`, `D₁`, `D₂`, `D₁⁺`, `D₂⁺`, and `ws` have identical object identity.
`Nt` is a positive odd unit-period time or phase resolution; `Nt=1` represents a steady field.

The grid includes boundary points when its FDGrids distributions do, but it does not prescribe the
driven wall, lid velocity, or no-slip values.

# Arguments

- `N`: number of collocation points along each edge of the cube.

# Keyword arguments

- `Nt`: positive odd temporal Fourier resolution; one represents a steady field.
- `lim`: common two-endpoint interval for all bounded directions.
- `dist`: FDGrids point distribution shared by `x`, `y`, and `z`.
- `width`: finite-difference stencil width shared by `x`, `y`, and `z`; weighted adjoints require
  `N > 2*width`.
- `T`: real scalar type used for points, operators, weights, and the Fourier scale.

# Returns

A [`LidDrivenCavity3DGrid`](@ref) whose bounded directions share the same numerical objects.

# Example

```julia
g = LidDrivenCavity3DGrid(49; Nt=1, width=7)
```
"""
function LidDrivenCavity3DGrid(N::Int; Nt::Int=1, lim::NTuple{2, <:Real}=(0, 1),
                               dist::FDGrids.AbstractGridDistribution=FDGrids.UniformGrid(),
                               width::Int=5, T::Type{<:Real}=Float64)
    x, D₁, D₂, D₁⁺, D₂⁺, w = _fd_direction(N, lim, dist, width, T)
    return RectangularGrid((x, x, x), (D₁, D₁, D₁), (D₂, D₂, D₂),
                           (D₁⁺, D₁⁺, D₁⁺), (D₂⁺, D₂⁺, D₂⁺), (w, w, w), (2π,),
                           (N, N, N, Nt), LID_DRIVEN_CAVITY_3D_AXES,
                           LID_DRIVEN_CAVITY_3D_FFT_ORDER, T)
end
