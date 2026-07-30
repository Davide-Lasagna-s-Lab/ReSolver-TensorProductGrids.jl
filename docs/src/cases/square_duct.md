# Square-duct flow

[`SquareDuctFlow`](@ref) combines a bounded square cross-section with periodic streamwise ``z`` and
time or phase ``t``. Velocity is ordered as ``\boldsymbol{u}=(u,v,w)``, so the physical streamwise
component is ``w``. The returned NSEBase operator represents

```math
\partial_t\boldsymbol{u}+(\boldsymbol{u}\boldsymbol{\cdot}\nabla)\boldsymbol{u}
=-\nabla p+\frac{1}{Re}\nabla^2\boldsymbol{u}+f\boldsymbol{e}_z,
\qquad \nabla\boldsymbol{\cdot}\boldsymbol{u}=0.
```

Writing ``G=-(1/\rho)d\bar p^*/dz^*>0`` for the dimensional pressure-gradient acceleration and
using the **full square side** ``L`` and a reference velocity ``U_{ref}``, the package convention is

```math
Re=\frac{U_{ref}L}{\nu},\qquad f=\frac{GL}{U_{ref}^2}.
```

The constant force is stored only at zero `(z,t)` wavenumber, which is exactly a spatially and
temporally uniform acceleration in physical space. It acts on component three in forward residual
evaluations. Because it is independent of the state, its linearisation—and hence its continuous
and discrete adjoint action—is zero. `SquareDuctFlow` requires this pressure-gradient amplitude to
be nonzero, so `f=0` throws an `ArgumentError`; use a separately assembled unforced formulation
when no prescribed pressure gradient is intended.

## Friction Reynolds number

A force balance over a duct segment gives the perimeter-averaged wall stress

```math
\bar\tau_w(4L)=\rho G L^2,
\qquad \bar\tau_w=\rho G\frac{L}{4}.
```

With ``u_\tau=\sqrt{\bar\tau_w/\rho}`` and the conventional half-side length ``L/2``, this yields

```math
Re_\tau=\frac{u_\tau(L/2)}{\nu}
=\frac{Re\sqrt{|f|}}{4}.
```

Therefore friction-velocity and half-side scaling is
`SquareDuctFlow(g, 2Reτ; f=4)`: the constructor's side-based Reynolds number is ``2Re_\tau`` and
``f=GL/u_\tau^2=4``. In contrast, `f=1` selects ``U_{ref}=\sqrt{GL}``, the pressure-gradient
velocity; it is not the conventional perimeter-averaged friction velocity.

```@example square_duct_case
using NSEBase
using ReSolverRectangularGrids

g = SquareDuctGrid(9, 5; Nt=1, α=0.5, width=3)
Reτ = 180
equations = SquareDuctFlow(g, 2Reτ; f=4, fftw_flags=FFTW.ESTIMATE, dealias=false)

@assert equations.nl.force isa ConstantBodyForce
@assert equations.nl.force.i == 3
@assert Reτ == (2Reτ) * sqrt(abs(equations.nl.force.value)) / 4
(grid_size=size(g), operator=typeof(equations.nl))
```

## Base flow and wall conditions

`base_flow=(U,V,W)` is added only to the steady zero `(z,t)` Fourier mode. Use
`(nothing,nothing,W)` to evaluate or linearise about a streamwise profile sampled over the bounded
cross-section; any identically zero component may be `nothing`. The default zero tuple describes
the pressure-driven equations without selecting a reference profile.

The grid supplies wall points and operators but does not impose no-slip values. The basis or
residual formulation must enforce the intended homogeneous conditions on perturbations, consistent
with any nonzero lifting in `base_flow`. `dealias=true` pads both Fourier directions for nonlinear
products. `mode=AdjointDiscrete()` selects the quadrature-consistent discrete adjoint of the
linearised operator, while `AdjointContinuous()` selects the discretised continuous-adjoint
equations.

```@docs
SquareDuctFlow
```
