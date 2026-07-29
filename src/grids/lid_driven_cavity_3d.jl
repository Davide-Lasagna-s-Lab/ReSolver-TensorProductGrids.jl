# =============================================================================================== #
# Three-dimensional lid-driven-cavity grid                                                       #
# =============================================================================================== #
#
# This file defines a fully bounded three-dimensional cavity. Physical `x`, `y`, and `z` use
# independent FDGrids discretisations, while `t` is the only Fourier time or phase direction. The
# grid supplies wall points, quadrature, and derivatives without prescribing which wall is driven.
#
# Typical use; positional resolutions follow physical order `(Nx, Ny, Nz)`, while the optional
# temporal resolution is the keyword `Nt`:
#
#     g = LidDrivenCavity3DGrid(49, 49, 49; Nt=1, width=7)
#     x, y, z, t = NSEBase.points(g)  # coordinates are returned in storage order
#
# Compatible wrappers can dispatch on `AbstractLidDrivenCavity3DGrid`; alternative serial layouts
# can be built directly with `RectangularGrid` and the constants below.

# =============================================================================================== #
# Three-dimensional layout                                                                       #
# =============================================================================================== #

"""Map physical `(x, y, z, t)` directly to the 3D cavity storage order `(x, y, z, t)`."""
const LID_DRIVEN_CAVITY_3D_AXES = (1, 2, 3, 4)

"""The 3D cavity's unit-period time or phase transform dimension."""
const LID_DRIVEN_CAVITY_3D_FFT_ORDER = (4,)

"""The bounded physical `(x, y, z)` storage dimensions of a 3D cavity."""
const LID_DRIVEN_CAVITY_3D_INHOMOGENEOUS_DIMS = (1, 2, 3)

"""
    AbstractLidDrivenCavity3DGrid{T}

Layout contract for a 3D cavity stored as `(x, y, z, t)` with bounded `x`, `y`, and `z`.
"""
const AbstractLidDrivenCavity3DGrid{T} =
    AbstractGrid{T, 4, LID_DRIVEN_CAVITY_3D_AXES, LID_DRIVEN_CAVITY_3D_FFT_ORDER}

"""
    LidDrivenCavity3DGrid

Concrete [`RectangularGrid{3}`](@ref) alias for the 3D cavity layout. Its stored size is
`(Nx, Ny, Nz, Nt)` and its only Fourier scale is the unit-period time or phase scale `2π`.
"""
const LidDrivenCavity3DGrid{T, S, XS, D1, D2, A1, A2, WS, WP} = RectangularGrid{
    3, T, S, 4, LID_DRIVEN_CAVITY_3D_AXES, LID_DRIVEN_CAVITY_3D_FFT_ORDER,
    XS, D1, D2, A1, A2, WS, WP, 1,
}

# =============================================================================================== #
# Three-dimensional constructor                                                                  #
# =============================================================================================== #

"""
    LidDrivenCavity3DGrid(Nx, Ny, Nz; Nt=1, lim=(0, 1),
                          xlim=lim, ylim=lim, zlim=lim,
                          dist=FDGrids.UniformGrid(),
                          xdist=dist, ydist=dist, zdist=dist,
                          width=5, xwidth=width, ywidth=width, zwidth=width,
                          T=Float64)

Construct a fully bounded 3D cavity stored as `(x, y, z, t)`. Each spatial direction may use an
independent interval, FDGrids distribution, and stencil width. `Nt` is a positive odd unit-period
time or phase resolution; `Nt=1` represents a steady field.

The grid includes boundary points when its FDGrids distributions do, but it does not prescribe the
driven wall, lid velocity, or no-slip values.

# Arguments

- `Nx`, `Ny`, `Nz`: numbers of bounded collocation points in physical `x`, `y`, and `z`.

# Keyword arguments

- `Nt`: positive odd temporal Fourier resolution; one represents a steady field.
- `lim`: common fallback interval for all bounded directions.
- `xlim`, `ylim`, `zlim`: two-endpoint intervals for `x`, `y`, and `z`.
- `dist`: common fallback FDGrids point distribution.
- `xdist`, `ydist`, `zdist`: FDGrids distributions for `x`, `y`, and `z`.
- `width`: common fallback finite-difference stencil width.
- `xwidth`, `ywidth`, `zwidth`: finite-difference stencil widths for `x`, `y`, and `z`; weighted
  adjoints require each resolution to exceed twice its corresponding width.
- `T`: real scalar type used for points, operators, weights, and the Fourier scale.

# Returns

A [`LidDrivenCavity3DGrid`](@ref) with independent FD operators in all three spatial directions.

# Example

```julia
g = LidDrivenCavity3DGrid(49, 49, 49; Nt=1, width=7)
```
"""
function LidDrivenCavity3DGrid(Nx::Int, Ny::Int, Nz::Int; Nt::Int=1, lim=(0, 1),
                               xlim=lim, ylim=lim, zlim=lim, dist=FDGrids.UniformGrid(),
                               xdist=dist, ydist=dist, zdist=dist, width=5,
                               xwidth=width, ywidth=width, zwidth=width,
                               T::Type{<:Real}=Float64)
    x, Dx, Dx2, Dxa, Dx2a, wx = _fd_direction(Nx, xlim, xdist, xwidth, T)
    y, Dy, Dy2, Dya, Dy2a, wy = _fd_direction(Ny, ylim, ydist, ywidth, T)
    z, Dz, Dz2, Dza, Dz2a, wz = _fd_direction(Nz, zlim, zdist, zwidth, T)
    return RectangularGrid((x, y, z), (Dx, Dy, Dz), (Dx2, Dy2, Dz2),
                           (Dxa, Dya, Dza), (Dx2a, Dy2a, Dz2a), (wx, wy, wz), (2π,),
                           (Nx, Ny, Nz, Nt), LID_DRIVEN_CAVITY_3D_AXES,
                           LID_DRIVEN_CAVITY_3D_FFT_ORDER, T)
end
