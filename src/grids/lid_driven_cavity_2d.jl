# =============================================================================================== #
# Two-dimensional lid-driven-cavity grid                                                         #
# =============================================================================================== #
#
# This file defines the numerical layout for a two-dimensional lid-driven cavity. Physical `x`
# and `y` are bounded FDGrids directions, while `t` is a Fourier time or phase direction. The grid
# provides coordinates, quadrature, and differentiation; it does not encode the moving lid or any
# other velocity boundary condition.
#
# Typical use; positional resolutions follow physical order `(Nx, Ny)`, while the optional temporal
# resolution is the keyword `Nt`:
#
#     g = LidDrivenCavity2DGrid(65, 49; Nt=1, width=7)
#     x, y, t = NSEBase.points(g)  # coordinates are returned in storage order
#
# To extend the case, construct a compatible `RectangularGrid` with the same axes and Fourier order,
# or dispatch downstream algorithms on `AbstractLidDrivenCavity2DGrid`.

# =============================================================================================== #
# Two-dimensional layout                                                                         #
# =============================================================================================== #

"""Map physical `(x, y, z, t)` to the 2D cavity storage order `(x, y, t)`."""
const LID_DRIVEN_CAVITY_2D_AXES = (1, 2, nothing, 3)

"""The 2D cavity's unit-period time or phase transform dimension."""
const LID_DRIVEN_CAVITY_2D_FFT_ORDER = (3,)

"""The bounded physical `(x, y)` storage dimensions of a 2D cavity."""
const LID_DRIVEN_CAVITY_2D_INHOMOGENEOUS_DIMS = (1, 2)

"""
    AbstractLidDrivenCavity2DGrid{T}

Layout contract for a 2D cavity stored as `(x, y, t)` with bounded `x` and `y`.
"""
const AbstractLidDrivenCavity2DGrid{T} =
    AbstractGrid{T, 3, LID_DRIVEN_CAVITY_2D_AXES, LID_DRIVEN_CAVITY_2D_FFT_ORDER}

"""
    LidDrivenCavity2DGrid

Concrete [`RectangularGrid{2}`](@ref) alias for the 2D cavity layout. Its stored size is
`(Nx, Ny, Nt)` and its only Fourier scale is the unit-period time or phase scale `2π`.
"""
const LidDrivenCavity2DGrid{T, S, XS, D1, D2, A1, A2, WS, WP} = RectangularGrid{
    2, T, S, 3, LID_DRIVEN_CAVITY_2D_AXES, LID_DRIVEN_CAVITY_2D_FFT_ORDER,
    XS, D1, D2, A1, A2, WS, WP, 1,
}

# =============================================================================================== #
# Two-dimensional constructor                                                                    #
# =============================================================================================== #

"""
    LidDrivenCavity2DGrid(Nx, Ny; Nt=1, lim=(0, 1), xlim=lim, ylim=lim,
                          dist=FDGrids.UniformGrid(), xdist=dist, ydist=dist,
                          width=5, xwidth=width, ywidth=width, T=Float64)

Construct a bounded 2D cavity stored as `(x, y, t)`. `xlim` and `ylim` may differ, as may the
FDGrids distributions and stencil widths in each direction. `Nt` is a positive odd unit-period
time or phase resolution; `Nt=1` represents a steady field.

The grid includes boundary points when its FDGrids distributions do, but does not prescribe a
moving-lid profile or no-slip values.

# Arguments

- `Nx`, `Ny`: numbers of bounded collocation points in physical `x` and `y`.

# Keyword arguments

- `Nt`: positive odd temporal Fourier resolution; one represents a steady field.
- `lim`: common fallback interval for both bounded directions.
- `xlim`, `ylim`: two-endpoint intervals for `x` and `y`.
- `dist`: common fallback FDGrids point distribution.
- `xdist`, `ydist`: FDGrids distributions for `x` and `y`.
- `width`: common fallback finite-difference stencil width.
- `xwidth`, `ywidth`: finite-difference stencil widths for `x` and `y`; weighted adjoints require
  `Nx > 2*xwidth` and `Ny > 2*ywidth`.
- `T`: real scalar type used for points, operators, weights, and Fourier scales.

# Returns

A [`LidDrivenCavity2DGrid`](@ref) with independent FD operators in `x` and `y`.
"""
function LidDrivenCavity2DGrid(Nx::Int, Ny::Int; Nt::Int=1, lim=(0, 1), xlim=lim,
                               ylim=lim, dist=FDGrids.UniformGrid(), xdist=dist, ydist=dist,
                               width=5, xwidth=width, ywidth=width, T::Type{<:Real}=Float64)
    x, Dx, Dx2, Dxa, Dx2a, wx = _fd_direction(Nx, xlim, xdist, xwidth, T)
    y, Dy, Dy2, Dya, Dy2a, wy = _fd_direction(Ny, ylim, ydist, ywidth, T)
    return RectangularGrid((x, y), (Dx, Dy), (Dx2, Dy2), (Dxa, Dya), (Dx2a, Dy2a),
                           (wx, wy), (2π,), (Nx, Ny, Nt), LID_DRIVEN_CAVITY_2D_AXES,
                           LID_DRIVEN_CAVITY_2D_FFT_ORDER, T)
end
