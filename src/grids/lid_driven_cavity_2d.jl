# =============================================================================================== #
# Two-dimensional lid-driven-cavity grid                                                         #
# =============================================================================================== #
#
# This file defines the numerical layout for a square two-dimensional lid-driven cavity. Physical
# `x` and `y` share one bounded FDGrids discretisation, while `t` is a Fourier time or phase
# direction. The grid provides coordinates, quadrature, and differentiation; it does not encode
# the moving lid or any other velocity boundary condition.
#
# Typical use; the positional resolution `N` is shared by `x` and `y`, while the optional temporal
# resolution is the keyword `Nt`:
#
#     g = LidDrivenCavity2DGrid(65; Nt=1, width=7)
#     x, y, t = NSEBase.points(g)  # coordinates are returned in storage order
#
# Downstream algorithms can dispatch on `AbstractLidDrivenCavity2DGrid` without depending on this
# concrete FDGrids-backed representation.

# =============================================================================================== #
# Layout constants                                                                               #
# =============================================================================================== #

"""
    LID_DRIVEN_CAVITY_2D_AXES = (1, 2, nothing, 3)

Map physical `(x, y, z, t)` to 2D-cavity storage order `(x, y, t)`. Bounded `x` and `y` occupy
the first two dimensions, physical `z` is absent, and time or phase `t` is stored third.
"""
const LID_DRIVEN_CAVITY_2D_AXES = (1, 2, nothing, 3)

"""
    LID_DRIVEN_CAVITY_2D_FFT_ORDER = (3,)

The 2D cavity's only Fourier dimension, corresponding to the ``2π``-periodic phase `t`.
"""
const LID_DRIVEN_CAVITY_2D_FFT_ORDER = (3,)

"""
    LID_DRIVEN_CAVITY_2D_INHOMOGENEOUS_DIMS = (1, 2)

The two FDGrids dimensions, corresponding to bounded physical coordinates `(x, y)`.
"""
const LID_DRIVEN_CAVITY_2D_INHOMOGENEOUS_DIMS = (1, 2)

# =============================================================================================== #
# Abstract layout contract                                                                       #
# =============================================================================================== #

"""
    AbstractLidDrivenCavity2DGrid{T}

Layout contract for any real scalar type `T` stored as `(x, y, t)` with bounded `x` and `y` and a
Fourier time or phase direction `t`.

Equal side lengths and shared numerical objects are guarantees of the concrete constructor, not
requirements imposed by this abstract layout contract.

The alias describes a layout rather than one concrete representation. NSEBase's device and
domain-decomposition wrappers preserve the parent's `AbstractGrid` layout parameters, so wrapped
2D-cavity grids still satisfy `AbstractLidDrivenCavity2DGrid`. Downstream algorithms can therefore
use one dispatch for serial, GPU-backed, and decomposed grids with this layout.
"""
const AbstractLidDrivenCavity2DGrid{T} =
    AbstractGrid{T, 3, LID_DRIVEN_CAVITY_2D_AXES, LID_DRIVEN_CAVITY_2D_FFT_ORDER}

# =============================================================================================== #
# Concrete rectangular-grid alias                                                                #
# =============================================================================================== #

"""
    LidDrivenCavity2DGrid

Concrete [`RectangularGrid{2}`](@ref) alias for the 2D cavity layout. Its stored size is
`(N, N, Nt)` and its only Fourier scale is the phase scale one, with
``t∈[0,2π)``. The
constructor shares its points, derivative matrices, adjoints, and weights between `x` and `y`.
"""
const LidDrivenCavity2DGrid{T, S, XS, D1, D2, A1, A2, WS, WP} = RectangularGrid{
    2, T, S, 3, LID_DRIVEN_CAVITY_2D_AXES, LID_DRIVEN_CAVITY_2D_FFT_ORDER,
    XS, D1, D2, A1, A2, WS, WP, 1,
}

# =============================================================================================== #
# FDGrids-backed constructor                                                                     #
# =============================================================================================== #

"""
    LidDrivenCavity2DGrid(N; Nt=1, lim=(0, 1),
                          dist=FDGrids.UniformGrid(), width=5, T=Float64)

Construct a square bounded 2D cavity stored as `(x, y, t)`. Both spatial directions use `N`
collocation points on `lim`. One FDGrids discretisation is reused in `x` and `y`, so corresponding
entries in `xs`, `D₁`, `D₂`, `D₁⁺`, `D₂⁺`, and `ws` have identical object identity. `Nt` is a
positive odd phase resolution on ``[0,2π)``; `Nt=1` represents a steady field.

The grid includes boundary points when its FDGrids distributions do, but does not prescribe a
moving-lid profile or no-slip values.

# Arguments

- `N`: number of collocation points along each side of the square.

# Keyword arguments

- `Nt`: positive odd temporal Fourier resolution; one represents a steady field.
- `lim`: common two-endpoint interval for both bounded directions.
- `dist`: FDGrids point distribution shared by `x` and `y`.
- `width`: finite-difference stencil width shared by `x` and `y`; weighted adjoints require
  `N > 2*width`.
- `T`: real scalar type used for points, operators, weights, and Fourier scales.

# Returns

A [`LidDrivenCavity2DGrid`](@ref) whose bounded directions share the same numerical objects.

# Example

```julia
g = LidDrivenCavity2DGrid(65; Nt=1, lim=(-1, 1), width=7)
```
"""
function LidDrivenCavity2DGrid(N::Int; Nt::Int=1, lim::NTuple{2, <:Real}=(0, 1),
                               dist::FDGrids.AbstractGridDistribution=FDGrids.UniformGrid(),
                               width::Int=5, T::Type{<:Real}=Float64)
    x, D₁, D₂, D₁⁺, D₂⁺, w = _fd_direction(N, lim, dist, width, T)
    return RectangularGrid((x, x), (D₁, D₁), (D₂, D₂), (D₁⁺, D₁⁺), (D₂⁺, D₂⁺),
                           (w, w), (1,), (N, N, Nt), LID_DRIVEN_CAVITY_2D_AXES,
                           LID_DRIVEN_CAVITY_2D_FFT_ORDER, T)
end
