# =============================================================================================== #
# Channel-flow cases                                                                              #
# =============================================================================================== #
#
# The grid layout lives in `grids/channel.jl`. This file supplies canonical serial wall-normal
# profiles and equation factories for plane Couette and plane Poiseuille flow on that layout.

# =============================================================================================== #
# Canonical wall-normal profiles                                                                  #
# =============================================================================================== #

"""
    plane_couette_base(g::AbstractChannelGrid) -> Vector

Return the canonical streamwise Couette base flow `U(y) = y` on the channel's fixed wall-normal
interval `[-1, 1]`.

# Arguments

- `g`: channel grid whose wall-normal points define the returned profile.

# Returns

A newly allocated vector that does not alias the grid's collocation points.
"""
plane_couette_base(g::AbstractChannelGrid) = copy(vec(points(g)[1]))

"""
    plane_poiseuille_base(g::AbstractChannelGrid) -> Vector

Return the canonical streamwise Poiseuille base flow `U(y) = 1 - y²` on the fixed interval
`[-1, 1]`.

# Arguments

- `g`: channel grid whose wall-normal points define the returned profile.

# Returns

A newly allocated vector with centreline value one, wall values zero, and bulk value `2/3`.
"""
plane_poiseuille_base(g::AbstractChannelGrid) = one(eltype(g)) .- plane_couette_base(g) .^ 2

# =============================================================================================== #
# Plane Couette equations                                                                         #
# =============================================================================================== #

@doc raw"""
    PlaneCouetteFlow(g::AbstractChannelGrid, Re;
                     Ro=0, base_flow=(plane_couette_base(g), nothing, nothing),
                     mode=AdjointDiscrete(), fftw_flags=FFTW.EXHAUSTIVE,
                     dealias=true) -> ProjectedNSE

Construct the incompressible plane Couette equations. With channel half-height `h`, reference
velocity `U_ref`, and kinematic viscosity `ν`, the Reynolds number is

```math
Re = \frac{U_{ref}h}{\nu}.
```

For the default symmetric Couette flow, `U_ref` is the magnitude of either wall velocity and the
nondimensional walls move at `u=±1`. For total velocity `u = (u, v, w)`, the primitive equations
represented by the returned NSEBase operator are

```math
\partial_t \boldsymbol{u} + (\boldsymbol{u}\cdot\nabla)\boldsymbol{u}
= -\nabla p + Re^{-1}\nabla^2\boldsymbol{u} + \boldsymbol{f}_{Ro},
\qquad \nabla\cdot\boldsymbol{u}=0.
```

For spanwise rotation `Ω = Ω e_z`, the signed rotation number and Coriolis term use the convention

```math
Ro = \frac{2\Omega h}{U_{ref}}, \qquad
\boldsymbol{f}_{Ro} = -Ro\,\boldsymbol{e}_z\times\boldsymbol{u}
                    = Ro\,(v,-u,0).
```

# Arguments

- `g`: channel grid stored as `(y, x, z, t)`, with `y∈[-1,1]`.
- `Re`: Reynolds number `U_ref*h/ν` multiplying viscosity as `1/Re`.

# Keyword arguments

- `Ro`: signed rotation number; zero omits the Coriolis force.
- `base_flow`: three-tuple `(U, V, W)` of wall-normal profiles. An identically zero component may be
  `nothing`. The default is `(plane_couette_base(g), nothing, nothing)`.
- `mode`: `AdjointDiscrete()` for the quadrature-consistent discrete adjoint or
  `AdjointContinuous()` for the continuous adjoint of the linearised operator.
- `fftw_flags`: FFTW planner flags forwarded to `NSEBase.construct_equations`.
- `dealias`: whether nonlinear products use padded Fourier directions.

# Boundary conditions

`ChannelGrid` supplies wall points and numerical operators but does not constrain field values.
`ProjectedNSE` adds `base_flow` only to the steady zero `(x,z,t)` Fourier mode. With the default
lifting, perturbation modes that vanish at `y=±1` give total wall velocities `(-1,0,0)` and
`(1,0,0)`. Encode any other prescribed wall motion in `base_flow` and keep the perturbation basis
homogeneous.

# Returns

An `NSEBase.ProjectedNSE` containing the nonlinear and selected linearised-adjoint operators.
"""
function PlaneCouetteFlow(g::AbstractChannelGrid, Re; Ro=0,
                          base_flow=(plane_couette_base(g), nothing, nothing),
                          mode=AdjointDiscrete(), fftw_flags=FFTW.EXHAUSTIVE, dealias=true)
    force = _coriolis_force(eltype(g)(Ro))
    return _plane_channel_flow(g, Re, base_flow, force; mode, fftw_flags, dealias)
end

# =============================================================================================== #
# Plane Poiseuille equations                                                                      #
# =============================================================================================== #

@doc raw"""
    PlanePoiseuilleFlow(g::AbstractChannelGrid, Re;
                        Ro=0, f=1,
                        base_flow=(plane_poiseuille_base(g), nothing, nothing),
                        mode=AdjointDiscrete(), fftw_flags=FFTW.EXHAUSTIVE,
                        dealias=true) -> ProjectedNSE

Construct the pressure-driven plane Poiseuille equations. Using the Reynolds and rotation-number
definitions from [`PlaneCouetteFlow`](@ref), the total velocity satisfies

```math
\partial_t \boldsymbol{u} + (\boldsymbol{u}\cdot\nabla)\boldsymbol{u}
= -\nabla p + Re^{-1}\nabla^2\boldsymbol{u} + f\,\boldsymbol{e}_x
  + \boldsymbol{f}_{Ro},
\qquad \nabla\cdot\boldsymbol{u}=0.
```

Here `f` is the signed dimensionless mean pressure-gradient acceleration. For dimensional mean
pressure `p̄*`, density `ρ`, and positive `x` downstream,

```math
f = -\frac{h}{\rho U_{ref}^2}\frac{d\bar p^*}{dx^*}.
```

# Arguments

- `g`: channel grid stored as `(y, x, z, t)`, with `y∈[-1,1]`.
- `Re`: Reynolds number `U_ref*h/ν` multiplying viscosity as `1/Re`.

# Keyword arguments

- `Ro`: signed rotation number; zero omits the Coriolis force.
- `f`: signed uniform streamwise forcing; positive values drive flow in `+x`.
- `base_flow`: three-tuple `(U, V, W)` of wall-normal profiles. An identically zero component may be
  `nothing`. The default is `(plane_poiseuille_base(g), nothing, nothing)`.
- `mode`, `fftw_flags`, `dealias`: as documented by [`PlaneCouetteFlow`](@ref).

# Friction and bulk-flow conventions

To prescribe the friction Reynolds number directly, use friction-velocity scaling:

```julia
U = (Reτ / 2) .* plane_poiseuille_base(g)
equations = PlanePoiseuilleFlow(g, Reτ; f=1, base_flow=(U, nothing, nothing))
```

Then the constructor argument is exactly `Reτ`, and `U=(Reτ/2)(1-y²)` is the corresponding
laminar equilibrium. More generally,

```math
Re_\tau = Re\sqrt{|f|},
```

or, for positive forcing, choose `f=(Reτ/Re)^2`.

For prescribed bulk velocity, pass `f=0` and constrain the streamwise component of every basis mode
in the zero spatial `(k_x,k_z)=(0,0)` sector to have zero quadrature-weighted wall-normal mean:

```math
\sum_j w_j\,\widehat{\phi}_x(0,y_j,0,k_t)=0.
```

Apply the constraint to every temporal mode `k_t` when the bulk velocity is fixed at every time.
The uniform pressure gradient is then an unrepresented Lagrange multiplier, and `Reτ` is recovered
from that gradient or the mean wall shear rather than prescribed through this constructor. The
default profile `1-y²` has bulk velocity `2/3` and is a convenient reference/lifting profile, not
an equilibrium for arbitrary `Re` and `f`. Use `3/2*(1-y²)` when the chosen bulk scale is one.

# Returns

An `NSEBase.ProjectedNSE` containing the nonlinear and selected linearised-adjoint operators.
"""
function PlanePoiseuilleFlow(g::AbstractChannelGrid, Re; Ro=0, f=1,
                             base_flow=(plane_poiseuille_base(g), nothing, nothing),
                             mode=AdjointDiscrete(), fftw_flags=FFTW.EXHAUSTIVE, dealias=true)
    force = _poiseuille_force(eltype(g)(Ro), eltype(g)(f))
    return _plane_channel_flow(g, Re, base_flow, force; mode, fftw_flags, dealias)
end

# =============================================================================================== #
# Shared channel-equation assembly                                                               #
# =============================================================================================== #

function _plane_channel_flow(g::AbstractChannelGrid, Re, base_flow, force;
                             mode, fftw_flags, dealias)
    # Require one zero-mode lifting profile for each Cartesian velocity component.
    base_flow isa Tuple{Any, Any, Any} ||
        throw(ArgumentError("a channel base flow must contain (U, V, W)"))

    # Build the nonlinear and requested adjoint operators in `(u, v, w)` order.
    return construct_equations(g, Re, base_flow, CartesianPrimitive3D();
                               force, mode, flags=fftw_flags, dealias)
end

# Avoid carrying a no-op Coriolis policy in the nonrotating case.
_coriolis_force(Ro) = iszero(Ro) ? NoForce() : CoriolisForce(Ro)

function _poiseuille_force(Ro, f)
    # Represent the imposed mean pressure gradient by a uniform streamwise DC force.
    pressure = ConstantBodyForce(f; i=1)

    # Compose the pressure and Coriolis contributions only when rotation is enabled.
    return iszero(Ro) ? pressure : CompoundForcing(pressure, CoriolisForce(Ro))
end
