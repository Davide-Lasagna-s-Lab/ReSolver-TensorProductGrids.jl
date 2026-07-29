# =============================================================================================== #
# Channel grid                                                                                    #
# =============================================================================================== #
#
# This file contains the complete numerical layout for a Fourier–finite-difference channel. The
# physical domain is bounded in wall-normal `y`, periodic in streamwise `x` and spanwise `z`, and
# may include a periodic time or phase coordinate `t`. The matching profiles and equation factories
# live separately in `cases/channel.jl`.
#
# Typical use; the positional resolutions follow physical order `(Nx, Ny, Nz)`, while the optional
# temporal resolution is the keyword `Nt`:
#
#     g = ChannelGrid(63, 65, 63; Nt=1, α=0.5, β=1, width=7)
#     y, x, z, t = NSEBase.points(g)  # coordinates are returned in storage order
#
# To extend the case, construct a compatible `RectangularGrid` with the same `CHANNEL_AXES` and
# `CHANNEL_FFT_ORDER`, or dispatch downstream algorithms on `AbstractChannelGrid`.

# =============================================================================================== #
# Layout constants                                                                               #
# =============================================================================================== #

"""
    CHANNEL_AXES = (2, 1, 3, 4)

Map physical `(x, y, z, t)` to channel storage order `(y, x, z, t)`. Wall-normal `y` is the leading
finite-difference dimension, followed by streamwise `x`, spanwise `z`, and time or phase `t`.
"""
const CHANNEL_AXES = (2, 1, 3, 4)

"""
    CHANNEL_FFT_ORDER = (2, 3, 4)

Channel Fourier dimensions in transform order. Streamwise `x` is the real-to-complex direction;
spanwise `z` and time `t` use complex transforms.
"""
const CHANNEL_FFT_ORDER = (2, 3, 4)

"""
    CHANNEL_INHOMOGENEOUS_DIMS = (1,)

The channel's single FDGrids dimension, corresponding to physical wall-normal `y`.
"""
const CHANNEL_INHOMOGENEOUS_DIMS = (1,)

# =============================================================================================== #
# Abstract layout contract                                                                       #
# =============================================================================================== #

"""
    AbstractChannelGrid{T}

Layout contract for any real scalar type `T` stored as `(y, x, z, t)` with Fourier directions
`(x, z, t)`.

The alias describes a layout rather than one concrete representation. NSEBase's device and
domain-decomposition wrappers preserve the parent's `AbstractGrid` layout parameters, so wrapped
channel grids still satisfy `AbstractChannelGrid`. Downstream algorithms can therefore use one
dispatch for serial, GPU-backed, and decomposed channel grids.
"""
const AbstractChannelGrid{T} = AbstractGrid{T, 4, CHANNEL_AXES, CHANNEL_FFT_ORDER}

# =============================================================================================== #
# Concrete rectangular-grid alias                                                                #
# =============================================================================================== #

"""
    ChannelGrid

Concrete [`RectangularGrid{1}`](@ref) alias for the channel layout. Its stored size is
`(Ny, Nx, Nz, Nt)` and its Fourier scales are `(α, β, 2π)`.
"""
const ChannelGrid{T, S, XS, D1, D2, A1, A2, WS, WP} = RectangularGrid{
    1, T, S, 4, CHANNEL_AXES, CHANNEL_FFT_ORDER, XS, D1, D2, A1, A2, WS, WP, 3,
}

# =============================================================================================== #
# FDGrids-backed constructor                                                                     #
# =============================================================================================== #

"""
    ChannelGrid(Nx, Ny, Nz; Nt=1, α=1, β=1,
                dist=FDGrids.GaussLobattoGrid(), width=5, T=Float64)

Construct a channel grid stored as `(y, x, z, t)`.

`Ny` wall-normal points span the fixed half-height interval `[-1, 1]` and use the selected FDGrids
distribution. `Nx`, `Nz`, and `Nt` are positive odd Fourier resolutions. The streamwise and
spanwise wavenumber scales are `α = 2π/Lx` and `β = 2π/Lz`; the unit-period time or phase scale is
`2π`. `Nt=1` represents a steady or time-independent field.

`width` is the odd FDGrids stencil width used by both first and second derivatives. Weighted
adjoints are constructed from the FDGrids quadrature weights. The default Gauss–Lobatto grid
includes both walls, but no velocity boundary condition is imposed by the grid itself.

# Arguments

- `Nx`: positive odd streamwise Fourier resolution.
- `Ny`: number of bounded wall-normal collocation points.
- `Nz`: positive odd spanwise Fourier resolution.

# Keyword arguments

- `Nt`: positive odd temporal Fourier resolution; one represents a steady field.
- `α`: streamwise wavenumber scale `2π/Lx`.
- `β`: spanwise wavenumber scale `2π/Lz`.
- `dist`: FDGrids distribution on the fixed wall-normal interval `[-1,1]`.
- `width`: finite-difference stencil width for first and second wall-normal derivatives; weighted
  adjoints require `Ny > 2*width`.
- `T`: real scalar type used for points, operators, weights, and Fourier scales.

# Returns

A [`ChannelGrid`](@ref) stored as `(y, x, z, t)`.

# Example

```julia
g = ChannelGrid(63, 65, 63; Nt=1, α=0.5, β=1, width=7)
```
"""
function ChannelGrid(Nx::Int, Ny::Int, Nz::Int; Nt::Int=1, α::Real=1, β::Real=1,
                     dist::FDGrids.AbstractGridDistribution=FDGrids.GaussLobattoGrid(),
                     width::Int=5, T::Type{<:Real}=Float64)
    y, Dy, Dy2, Dya, Dy2a, wy = _fd_direction(Ny, (-1, 1), dist, width, T)
    return RectangularGrid((y,), (Dy,), (Dy2,), (Dya,), (Dy2a,), (wy,), (α, β, 2π),
                           (Ny, Nx, Nz, Nt), CHANNEL_AXES, CHANNEL_FFT_ORDER, T)
end

# =============================================================================================== #
# Constructors from precomputed wall-normal data                                                 #
# =============================================================================================== #

"""
    ChannelGrid(y, Nx, Nz, Nt, α, β, Dy, Dy2, Dya, Dy2a, wy, T=Float64)

Construct a channel from precomputed wall-normal points, differentiation matrices, and quadrature
weights. The points must span the fixed wall interval `[-1, 1]`. `Dya` and `Dy2a` are the
quadrature-weighted discrete adjoints of `Dy` and `Dy2`; they are explicit inputs because a
precomputed discretisation owns the precise adjoint operators used by its tested formulation.

These constructors are useful when a flow package already owns a carefully tuned wall-normal
discretisation. All data are converted to `T` only when needed.

# Arguments

- `y`: wall-normal collocation vector spanning `[-1,1]`.
- `Nx`, `Nz`, `Nt`: positive odd streamwise, spanwise, and temporal Fourier resolutions.
- `α`, `β`: streamwise and spanwise wavenumber scales.
- `Dy`, `Dy2`: first- and second-derivative matrices in `y`.
- `Dya`, `Dy2a`: quadrature-weighted discrete adjoints of `Dy` and `Dy2`.
- `wy`: wall-normal quadrature weights.
- `T`: real scalar type used by the constructed grid.

# Returns

A [`ChannelGrid`](@ref) retaining already correctly typed input objects.
"""
function ChannelGrid(y::AbstractVector, Nx::Int, Nz::Int, Nt::Int, α::Real, β::Real,
                     Dy::AbstractMatrix, Dy2::AbstractMatrix, Dya::AbstractMatrix,
                     Dy2a::AbstractMatrix, wy::AbstractVector, ::Type{T}=Float64) where {T<:Real}
    isempty(y) && throw(ArgumentError("wall-normal points cannot be empty"))
    (minimum(y) ≈ -1 && maximum(y) ≈ 1) ||
        throw(ArgumentError("wall-normal points must span [-1, 1]"))
    Ny = length(y)
    return RectangularGrid((y,), (Dy,), (Dy2,), (Dya,), (Dy2a,), (wy,), (α, β, 2π),
                           (Ny, Nx, Nz, Nt), CHANNEL_AXES, CHANNEL_FFT_ORDER, T)
end
