# =============================================================================================== #
# Square-duct grid                                                                               #
# =============================================================================================== #
#
# This file defines a bounded square cross-section with periodic streamwise and temporal
# directions. Both cross-section axes deliberately share one FDGrids discretisation, preserving
# the symmetry and object-identity assumptions of the original ReSolver-SquareDuct package.
#
# Typical use; positional resolutions are the shared cross-section resolution `N` and streamwise
# resolution `Nz`. The optional temporal resolution and streamwise scale are keywords `Nt` and `α`:
#
#     g = SquareDuctGrid(49, 63; Nt=1, α=0.5, width=7)
#     x, y = g.xs
#
# The pressure-driven equation factory lives separately in `cases/square_duct.jl`. Downstream code
# can dispatch on `AbstractSquareDuctGrid` when defining related duct cases.

# =============================================================================================== #
# Layout constants and contracts                                                                 #
# =============================================================================================== #

"""
    SQUARE_DUCT_AXES = (1, 2, 3, 4)

Map physical `(x, y, z, t)` directly to square-duct storage dimensions. Bounded `x` and `y` form
the cross-section, periodic `z` is streamwise, and `t` is time or phase.
"""
const SQUARE_DUCT_AXES = (1, 2, 3, 4)

"""
    SQUARE_DUCT_FFT_ORDER = (3, 4)

Square-duct Fourier dimensions in transform order. Streamwise `z` is real-to-complex and
unit-period `t` uses a complex transform.
"""
const SQUARE_DUCT_FFT_ORDER = (3, 4)

"""The two FDGrids dimensions spanning the bounded square cross-section."""
const SQUARE_DUCT_INHOMOGENEOUS_DIMS = (1, 2)

"""
    AbstractSquareDuctGrid{T}

Layout contract for any real scalar type `T` stored as `(x, y, z, t)` with Fourier directions
`(z, t)`. Equal side lengths and shared operators are constructor guarantees, not abstract type
requirements.
"""
const AbstractSquareDuctGrid{T} = AbstractGrid{T, 4, SQUARE_DUCT_AXES, SQUARE_DUCT_FFT_ORDER}

"""
    SquareDuctGrid

Concrete [`RectangularGrid{2}`](@ref) alias for a square duct. The constructor deliberately shares
the same points, matrices, adjoints, and weights between bounded `x` and `y` directions.
"""
const SquareDuctGrid{T, S, XS, D1, D2, A1, A2, WS, WP} = RectangularGrid{
    2, T, S, 4, SQUARE_DUCT_AXES, SQUARE_DUCT_FFT_ORDER, XS, D1, D2, A1, A2, WS, WP, 2,
}

# =============================================================================================== #
# FDGrids-backed constructor                                                                     #
# =============================================================================================== #

"""
    SquareDuctGrid(N, Nz; Nt=1, α=1,
                   dist=FDGrids.GaussLobattoGrid(), width=5, T=Float64)

Construct a square-duct grid with bounded cross-section `[0, 1] × [0, 1]`, periodic streamwise
length `Lz = 2π/α`, and unit-period time or phase. The stored size is `(N, N, Nz, Nt)` and the
Fourier scales are `(α, 2π)`.

One FDGrids discretisation is reused in both cross-section directions, so corresponding entries in
`xs`, `D₁`, `D₂`, `D₁⁺`, `D₂⁺`, and `ws` have identical object identity. `Nz` and `Nt` must be
positive odd resolutions.

The grid provides wall points and differentiation operators but imposes no-slip or pressure-driven
flow conditions itself.

# Arguments

- `N`: number of collocation points along each bounded side of the square cross-section.
- `Nz`: positive odd streamwise Fourier resolution.

# Keyword arguments

- `Nt`: positive odd temporal Fourier resolution; one represents a steady field.
- `α`: streamwise wavenumber scale `2π/Lz`.
- `dist`: FDGrids distribution shared by both bounded directions.
- `width`: finite-difference stencil width shared by both bounded directions; weighted adjoints
  require `N > 2*width`.
- `T`: real scalar type used for points, operators, weights, and Fourier scales.

# Returns

A [`SquareDuctGrid`](@ref) whose two bounded directions share the same numerical objects.

# Example

```julia
g = SquareDuctGrid(49, 63; Nt=1, α=0.5, width=7)
```
"""
function SquareDuctGrid(N::Int, Nz::Int; Nt::Int=1, α::Real=1,
                        dist=FDGrids.GaussLobattoGrid(), width=5,
                        T::Type{<:Real}=Float64)
    x, D₁, D₂, D₁⁺, D₂⁺, w = _fd_direction(N, (0, 1), dist, width, T)
    return RectangularGrid((x, x), (D₁, D₁), (D₂, D₂), (D₁⁺, D₁⁺), (D₂⁺, D₂⁺),
                           (w, w), (α, 2π), (N, N, Nz, Nt), SQUARE_DUCT_AXES,
                           SQUARE_DUCT_FFT_ORDER, T)
end
