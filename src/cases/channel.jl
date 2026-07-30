# =============================================================================================== #
# Channel-flow cases                                                                             #
# =============================================================================================== #
#
# The grid layout lives in `grids/channel.jl`. This file provides canonical wall-normal profiles
# and equation factories for plane Couette and plane Poiseuille flow.
#
#     g = ChannelGrid(63, 65, 63; width=7)
#     couette = PlaneCouetteFlow(g, 400; fftw_flags=FFTW.ESTIMATE)
#     poiseuille = PlanePoiseuilleFlow(g, 400; fftw_flags=FFTW.ESTIMATE)
#
# Pass another `base_flow=(U, V, W)` to change the reference profile without changing the grid.

# =============================================================================================== #
# Canonical wall-normal profiles                                                                 #
# =============================================================================================== #

"""
    plane_couette_base(g::AbstractChannelGrid) -> Vector

Return the canonical streamwise Couette base profile ``U(y)=y`` on the channel's fixed wall-normal
interval ``[-1,1]``.

# Arguments

- `g`: channel grid whose wall-normal points define the returned profile.

# Returns

A newly allocated vector that does not alias the grid's collocation points.

# Example

```julia
g = ChannelGrid(31, 33, 31)
U = plane_couette_base(g)
```
"""
plane_couette_base(g::AbstractChannelGrid) = copy(vec(points(g)[1]))

"""
    plane_poiseuille_base(g::AbstractChannelGrid) -> Vector

Return the canonical streamwise Poiseuille base profile ``U(y)=1-y^2`` on the fixed wall-normal
interval ``[-1,1]``.

# Arguments

- `g`: channel grid whose wall-normal points define the returned profile.

# Returns

A newly allocated vector with centreline value one, wall values zero, and bulk value `2/3`.

# Example

```julia
g = ChannelGrid(31, 33, 31)
U = plane_poiseuille_base(g)
```
"""
plane_poiseuille_base(g::AbstractChannelGrid) = one(eltype(g)) .- plane_couette_base(g) .^ 2

# =============================================================================================== #
# Plane Couette equation constructor                                                             #
# =============================================================================================== #

@doc raw"""
    PlaneCouetteFlow(g::AbstractChannelGrid, Re::Real;
                     base_flow=(plane_couette_base(g), nothing, nothing), Ro::Real=0,
                     mode=AdjointDiscrete(), fftw_flags=FFTW.EXHAUSTIVE,
                     dealias=true) -> ProjectedNSE

Construct the incompressible plane Couette equations

```math
\partial_t \boldsymbol{u} + (\boldsymbol{u}\cdot\nabla)\boldsymbol{u}
= -\nabla p + Re^{-1}\nabla^2\boldsymbol{u} + \boldsymbol{f}_{Ro},
\qquad \nabla\cdot\boldsymbol{u}=0,
```

where ``\boldsymbol{u}=(u,v,w)``. With channel half-height ``h``, reference velocity ``U_{ref}``,
and kinematic viscosity ``\nu``,

```math
Re=\frac{U_{ref}h}{\nu}.
```

For the default symmetric flow, ``U_{ref}`` is the magnitude of either wall velocity and the
nondimensional walls move at ``u=\pm1``. Spanwise rotation
``\boldsymbol{\Omega}=\Omega\boldsymbol{e}_z`` uses

```math
Ro=\frac{2\Omega h}{U_{ref}},\qquad
\boldsymbol{f}_{Ro}=-Ro\,\boldsymbol{e}_z\times\boldsymbol{u}=Ro\,(v,-u,0).
```

`base_flow=(U,V,W)` is added only to the steady zero `(x,z,t)` Fourier mode. With the default
lifting, perturbations that vanish at ``y=\pm1`` recover total wall velocities ``(-1,0,0)`` and
``(1,0,0)``. The grid and equation constructor do not impose those perturbation boundary values.
Every nonzero base-flow component must have the bounded shape `(Ny,)`.

# Arguments

- `g`: channel grid stored as `(y,x,z,t)`, with ``y\in[-1,1]``.
- `Re`: real Reynolds number ``U_{ref}h/\nu``, multiplying viscosity as ``1/Re``.

# Keyword arguments

- `base_flow`: three-component wall-normal tuple added to the steady zero Fourier mode.
- `Ro`: signed spanwise rotation number; zero omits the Coriolis force.
- `mode`: `AdjointDiscrete()` or `AdjointContinuous()` for the linearised adjoint operator.
- `fftw_flags`: FFTW planner flags forwarded to `NSEBase.construct_equations`.
- `dealias`: whether nonlinear products use padded Fourier resolutions.

# Returns

An `NSEBase.ProjectedNSE` for three velocity components.

# Example

```julia
g = ChannelGrid(31, 33, 31; α=0.5, β=1, width=7)
equations = PlaneCouetteFlow(g, 400; Ro=0.1, fftw_flags=FFTW.ESTIMATE)
```
"""
function PlaneCouetteFlow(g::AbstractChannelGrid, Re::Real;
                          base_flow=(plane_couette_base(g), nothing, nothing), Ro::Real=0,
                          mode=AdjointDiscrete(), fftw_flags=FFTW.EXHAUSTIVE,
                          dealias::Bool=true)
    force = iszero(Ro) ? NoForce() : CoriolisForce(eltype(g)(Ro))
    return _plane_channel_flow(g, Re, base_flow, force; mode, fftw_flags, dealias)
end

# =============================================================================================== #
# Plane Poiseuille equation constructor                                                          #
# =============================================================================================== #

@doc raw"""
    PlanePoiseuilleFlow(g::AbstractChannelGrid, Re::Real;
                        base_flow=(plane_poiseuille_base(g), nothing, nothing),
                        f::Real=1, mode=AdjointDiscrete(), fftw_flags=FFTW.EXHAUSTIVE,
                        dealias=true) -> ProjectedNSE

Construct the pressure-driven plane Poiseuille equations

```math
\partial_t \boldsymbol{u} + (\boldsymbol{u}\cdot\nabla)\boldsymbol{u}
= -\nabla p + Re^{-1}\nabla^2\boldsymbol{u} + f\,\boldsymbol{e}_x,
\qquad \nabla\cdot\boldsymbol{u}=0,
```

using the Reynolds-number convention of [`PlaneCouetteFlow`](@ref). For dimensional mean pressure
``\bar p^*``, density ``\rho``, and positive ``x`` downstream, the signed uniform acceleration is

```math
f=-\frac{h}{\rho U_{ref}^2}\frac{d\bar p^*}{dx^*}.
```

`base_flow=(U,V,W)` is added only to the steady zero `(x,z,t)` Fourier mode. The default profile
has unit centreline velocity and vanishes at the walls, but it is a reference profile rather than
an equilibrium for arbitrary `Re` and `f`. Every nonzero component must have shape `(Ny,)`; the
grid and equation constructor leave perturbation boundary conditions to the basis or residual.

# Friction and bulk-flow conventions

For friction-velocity scaling, use

```julia
U = (Reτ / 2) .* plane_poiseuille_base(g)
equations = PlanePoiseuilleFlow(g, Reτ; base_flow=(U, nothing, nothing), f=1)
```

Then the constructor argument is `Reτ`, and the laminar equilibrium is
``U=(Re_\tau/2)(1-y^2)``. More generally,

```math
Re_\tau=Re\sqrt{|f|},
```

or, for positive forcing, choose ``f=(Re_\tau/Re)^2``. The unscaled default profile has bulk
velocity ``2/3``; use ``3(1-y^2)/2`` when the chosen bulk scale is one.

This constructor represents pressure-driven flow and therefore requires `f` to be nonzero. A
prescribed-bulk-velocity formulation instead uses `NSEBase.construct_equations` with `NoForce()`
and constrains the streamwise component of every basis mode in the zero spatial
``(k_x,k_z)=(0,0)`` sector to have zero quadrature-weighted mean:

```math
\sum_j w_j\,\widehat{\phi}_x(0,y_j,0,k_t)=0.
```

Apply the constraint to every temporal mode ``k_t`` when the bulk velocity is fixed at every time.
The uniform pressure gradient is then an unrepresented Lagrange multiplier, and ``Re_\tau`` is
recovered from that gradient or the mean wall shear.

# Arguments

- `g`: channel grid stored as `(y,x,z,t)`, with ``y\in[-1,1]``.
- `Re`: real Reynolds number ``U_{ref}h/\nu``, multiplying viscosity as ``1/Re``.

# Keyword arguments

- `base_flow`: three-component wall-normal tuple added to the steady zero Fourier mode.
- `f`: nonzero signed uniform streamwise forcing; positive values drive flow in ``+x``. Zero throws
  an `ArgumentError`.
- `mode`: `AdjointDiscrete()` or `AdjointContinuous()` for the linearised adjoint operator.
- `fftw_flags`: FFTW planner flags forwarded to `NSEBase.construct_equations`.
- `dealias`: whether nonlinear products use padded Fourier resolutions.

# Returns

An `NSEBase.ProjectedNSE` for three velocity components with constant streamwise forcing.

# Example

```julia
g = ChannelGrid(31, 33, 31; α=0.5, β=1, width=7)
Reτ = 180
U = (Reτ / 2) .* plane_poiseuille_base(g)
equations = PlanePoiseuilleFlow(g, Reτ; base_flow=(U, nothing, nothing), f=1,
                                 fftw_flags=FFTW.ESTIMATE)
```
"""
function PlanePoiseuilleFlow(g::AbstractChannelGrid, Re::Real;
                             base_flow=(plane_poiseuille_base(g), nothing, nothing),
                             f::Real=1, mode=AdjointDiscrete(),
                             fftw_flags=FFTW.EXHAUSTIVE, dealias::Bool=true)
    force = ConstantBodyForce(eltype(g)(f); i=1)
    return _plane_channel_flow(g, Re, base_flow, force; mode, fftw_flags, dealias)
end

# =============================================================================================== #
# Shared channel-equation assembly                                                               #
# =============================================================================================== #

function _plane_channel_flow(g::AbstractChannelGrid, Re::Real, base_flow, force;
                             mode, fftw_flags, dealias::Bool)
    _validate_base_flow(g, base_flow, Val(3), "channel")
    return construct_equations(g, Re, base_flow, CartesianPrimitive3D();
                               force, mode, flags=fftw_flags, dealias)
end
