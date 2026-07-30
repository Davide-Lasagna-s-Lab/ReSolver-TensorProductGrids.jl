# Three-dimensional cubic lid-driven-cavity case

[`LidDrivenCavity3DFlow`](@ref) selects NSEBase's three-component Cartesian primitive formulation on a cubic `(x,y,z)` domain. For velocity ``\boldsymbol{u}=(u,v,w)``, it represents

```math
\partial_t\boldsymbol{u}+(\boldsymbol{u}\boldsymbol{\cdot}\nabla)\boldsymbol{u}
=-\nabla p+\frac{1}{Re}\nabla^2\boldsymbol{u},
\qquad \nabla\boldsymbol{\cdot}\boldsymbol{u}=0.
```

For a cavity of side ``L``, lid speed ``U_{lid}``, and kinematic viscosity ``\nu``, the canonical Reynolds number is ``Re=U_{lid}L/\nu``.

## Canonical moving-lid lifting

[`lid_driven_cavity_3d_base`](@ref) returns the default `base_flow=(U,V,nothing)`. Let ``\xi``,
``\eta``, and ``\zeta`` be the common cavity interval normalized to ``[0,1]``, and define

```math
F(\xi)=16\xi^2(1-\xi)^2,
\qquad H(\eta)=\eta^2(\eta-1),
\qquad G(\zeta)=16\zeta^2(1-\zeta)^2.
```

The helper evaluates

```math
U(\xi,\eta,\zeta)=F(\xi)H'(\eta)G(\zeta),
\qquad V(\xi,\eta,\zeta)=-F'(\xi)H(\eta)G(\zeta),
\qquad W=0.
```

The first two terms are divergence free and ``W`` is represented by `nothing`. At the upper wall,
``U=FG``; the polynomial factors taper the moving-wall velocity smoothly to zero at all four lid
edges. The field vanishes on every other wall. This field lifts boundary data; it is not claimed
to solve the steady Navier--Stokes equations.

## Symmetric equation-constructor interface

```julia
LidDrivenCavity3DFlow(g, Re;
    base_flow=lid_driven_cavity_3d_base(g),
    mode=AdjointDiscrete(), fftw_flags=FFTW.EXHAUSTIVE, dealias=true)
```

Every nonzero `base_flow` component must be an array of bounded shape `(N,N,N)`; an identically zero component may be `nothing`. `mode` selects the discrete or continuous adjoint of the linearised equations, `fftw_flags` controls FFT planning, and `dealias` controls temporal padding for nonlinear products. Passing `base_flow` explicitly replaces only the canonical lifting.

```@example cavity_3d_case
using NSEBase
using ReSolverRectangularGrids

g = LidDrivenCavity3DGrid(7; Nt=1, width=3)
base_flow = lid_driven_cavity_3d_base(g)
equations = LidDrivenCavity3DFlow(g, 1000; fftw_flags=FFTW.ESTIMATE, dealias=false)

U, V, W = base_flow
x, z = vec(points(g)[1]), vec(points(g)[3])
Fx = @. 16x^2 * (1 - x)^2
Gz = @. 16z^2 * (1 - z)^2
lid = Fx * transpose(Gz)
@assert maximum(abs, U[:, end, :] .- lid) < 100eps()
@assert isnothing(W) && equations.base == base_flow
(grid_size=size(g), operator=typeof(equations.nl))
```

## Boundary and spectral semantics

The lifting is inserted only in the zero temporal Fourier mode. Neither [`LidDrivenCavity3DGrid`](@ref) nor the equation factory constrains field values. Perturbations must satisfy homogeneous wall conditions through the basis or residual formulation so that lifting plus perturbation has the intended total wall velocity.

Time or phase is the only Fourier direction. `dealias=true` pads that direction during nonlinear products without changing the bounded resolution. `mode=AdjointDiscrete()` gives the quadrature-consistent discrete adjoint of the linearised operator, whereas `AdjointContinuous()` discretises the continuous-adjoint equations.

```@docs
lid_driven_cavity_3d_base
LidDrivenCavity3DFlow
```
