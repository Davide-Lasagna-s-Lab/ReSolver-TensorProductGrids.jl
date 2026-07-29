# Channel flows

The channel cases use the same [`ChannelGrid`](@ref) and differ only in their reference profile
and forcing policy. Physical velocity is ordered as ``\boldsymbol{u}=(u,v,w)``, with ``x``
streamwise, ``y\in[-1,1]`` wall-normal, and ``z`` spanwise. Arrays are stored as `(y,x,z,t)` for
efficient wall-normal finite differences; the physical derivative names remain `ddx!`, `ddy!`,
and `ddz!`.

The constructors return NSEBase primitive-variable operators for

```math
\partial_t\boldsymbol{u}+(\boldsymbol{u}\boldsymbol{\cdot}\nabla)\boldsymbol{u}
=-\nabla p+\frac{1}{Re}\nabla^2\boldsymbol{u}+\boldsymbol{f},
\qquad \nabla\boldsymbol{\cdot}\boldsymbol{u}=0,
```

where ``Re=U_{ref}h/\nu`` and the nondimensional channel half-height is one. The grid includes the
two walls but does not impose velocity boundary conditions.

## Plane Couette flow

[`PlaneCouetteFlow`](@ref) defaults to the canonical streamwise lifting
``\boldsymbol{U}=(y,0,0)``. With homogeneous perturbation boundary conditions, this represents
walls moving at velocities ``(-1,0,0)`` and ``(1,0,0)``. In the default scaling, ``U_{ref}`` is
the magnitude of either wall velocity.

An optional spanwise rotation ``\boldsymbol{\Omega}=\Omega\boldsymbol{e}_z`` adds

```math
Ro=\frac{2\Omega h}{U_{ref}},\qquad
\boldsymbol{f}_{Ro}=-Ro\,\boldsymbol{e}_z\times\boldsymbol{u}
=Ro\,(v,-u,0).
```

The sign of `Ro` therefore determines the sense of rotation. Passing `Ro=0` omits the Coriolis
policy entirely.

## Plane Poiseuille flow

[`PlanePoiseuilleFlow`](@ref) adds the signed, uniform streamwise acceleration
``f\boldsymbol{e}_x``. If ``G=-(1/\rho)d\bar p^*/dx^*`` is the dimensional pressure-gradient
acceleration, then

```math
f=\frac{Gh}{U_{ref}^2},\qquad
\partial_t\boldsymbol{u}+(\boldsymbol{u}\boldsymbol{\cdot}\nabla)\boldsymbol{u}
=-\nabla p+\frac{1}{Re}\nabla^2\boldsymbol{u}
+f\boldsymbol{e}_x+\boldsymbol{f}_{Ro}.
```

For a unidirectional laminar equilibrium ``U=A(1-y^2)``, the streamwise momentum balance gives
``A=Re\,f/2``. Consequently, friction-velocity scaling uses

```math
Re_\tau=Re\sqrt{|f|},\qquad
U(y)=\frac{Re_\tau}{2}(1-y^2)\quad\text{when }Re=Re_\tau,\ f=1.
```

The default [`plane_poiseuille_base`](@ref) is the unit-centreline profile ``1-y^2``. It is a
convenient reference or boundary-data lifting, but it is not the equilibrium for arbitrary `Re`
and `f`. In particular, setting `f=1` and passing `Reτ` requires scaling the profile by `Reτ/2`.
For a different velocity scale, use ``f=(Re_\tau/Re)^2`` for positive downstream forcing.

```@example channel_cases
using NSEBase
using ReSolverRectangularGrids

g = ChannelGrid(7, 13, 7; Nt=1, α=0.5, β=1, width=5)

couette = PlaneCouetteFlow(g, 400; Ro=0.1, fftw_flags=FFTW.ESTIMATE, dealias=false)

Reτ = 180
U = (Reτ / 2) .* plane_poiseuille_base(g)
poiseuille = PlanePoiseuilleFlow(g, Reτ; base_flow=(U, nothing, nothing), f=1,
                                 fftw_flags=FFTW.ESTIMATE, dealias=false)

@assert couette.base[1] == plane_couette_base(g)
@assert poiseuille.base[1] == U
(couette=typeof(couette.nl), poiseuille=typeof(poiseuille.nl))
```

## Base flows and boundary conditions

`base_flow=(U,V,W)` is added only to the steady zero `(x,z,t)` Fourier mode. It is the total-flow
reference about which the equations are evaluated or linearised; an identically zero component
may be written as `nothing`. The tuple does not itself constrain any perturbation value. A basis or
residual formulation must enforce homogeneous wall conditions consistent with the chosen lifting.

For prescribed bulk velocity, use `f=0`, place the prescribed mean profile in `base_flow`, and
require every streamwise component in the spatially homogeneous sector to have zero
quadrature-weighted wall-normal mean:

```math
\sum_j w_j\,\widehat\phi_x(0,y_j,0,k_t)=0.
```

Apply this constraint to every temporal harmonic when bulk velocity is fixed at every time. The
pressure gradient then acts as an eliminated Lagrange multiplier rather than a represented body
force; recover ``Re_\tau`` from that gradient or from the mean wall shear.

`mode=AdjointDiscrete()` selects the quadrature-consistent discrete adjoint of the linearised
operator, while `AdjointContinuous()` selects the discretised continuous-adjoint equations.
`dealias=true` pads all active Fourier directions for nonlinear products; disabling it is useful
for small examples but is generally inappropriate for production nonlinear calculations.

```@docs
plane_couette_base
plane_poiseuille_base
PlaneCouetteFlow
PlanePoiseuilleFlow
```
