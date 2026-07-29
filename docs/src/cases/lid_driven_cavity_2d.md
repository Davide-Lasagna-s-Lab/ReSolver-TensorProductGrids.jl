# Two-dimensional lid-driven-cavity case

[`LidDrivenCavity2DFlow`](@ref) selects NSEBase's two-component Cartesian primitive formulation on a square `(x,y)` domain. For velocity ``\boldsymbol{u}=(u,v)``, it represents

```math
\partial_t\boldsymbol{u}+(\boldsymbol{u}\boldsymbol{\cdot}\nabla)\boldsymbol{u}
=-\nabla p+\frac{1}{Re}\nabla^2\boldsymbol{u},
\qquad \nabla\boldsymbol{\cdot}\boldsymbol{u}=0.
```

For a cavity of side ``L``, lid speed ``U_{lid}``, and kinematic viscosity ``\nu``, the canonical Reynolds number is ``Re=U_{lid}L/\nu``.

## Canonical moving-lid lifting

[`lid_driven_cavity_2d_base`](@ref) returns the default `base_flow=(U,V)`. Let ``\xi`` and ``\eta`` be the common cavity interval normalized to `[0,1]`, and define

```math
F(\xi)=16\xi^2(1-\xi)^2,
\qquad H(\eta)=\eta^2(\eta-1).
```

The helper evaluates

```math
U(\xi,\eta)=F(\xi)H'(\eta),
\qquad V(\xi,\eta)=-F'(\xi)H(\eta).
```

The two terms in ``\partial_xU+\partial_yV`` cancel because both coordinates use the same interval. At the upper wall, ``U=F`` and ``V=0``; both components vanish on the other walls. The polynomial lid tapers smoothly to zero at the two upper corners, avoiding the discontinuities of an ideal uniform lid. This field lifts boundary data; it is not claimed to solve the steady Navier--Stokes equations.

## Symmetric equation-constructor interface

```julia
LidDrivenCavity2DFlow(g, Re;
    base_flow=lid_driven_cavity_2d_base(g),
    mode=AdjointDiscrete(), fftw_flags=FFTW.EXHAUSTIVE, dealias=true)
```

Every nonzero `base_flow` component must be an array of bounded shape `(N,N)`; an identically zero component may be `nothing`. `mode` selects the discrete or continuous adjoint of the linearised equations, `fftw_flags` controls FFT planning, and `dealias` controls temporal padding for nonlinear products. Passing `base_flow` explicitly replaces only the canonical lifting.

```@example cavity_2d_case
using NSEBase
using ReSolverRectangularGrids

g = LidDrivenCavity2DGrid(9; Nt=1, width=3)
base_flow = lid_driven_cavity_2d_base(g)
equations = LidDrivenCavity2DFlow(g, 1000; fftw_flags=FFTW.ESTIMATE, dealias=false)

U, V = base_flow
x = vec(points(g)[1])
lid = @. 16x^2 * (1 - x)^2
@assert maximum(abs, U[:, end] .- lid) < 100eps()
@assert equations.base == base_flow
(grid_size=size(g), operator=typeof(equations.nl))
```

## Boundary and spectral semantics

The lifting is inserted only in the zero temporal Fourier mode. Neither [`LidDrivenCavity2DGrid`](@ref) nor the equation factory constrains field values. Perturbations must satisfy homogeneous wall conditions through the basis or residual formulation so that lifting plus perturbation has the intended total wall velocity.

Time or phase is the only Fourier direction. `dealias=true` pads that direction during nonlinear products without changing the bounded resolution. `mode=AdjointDiscrete()` gives the quadrature-consistent discrete adjoint of the linearised operator, whereas `AdjointContinuous()` discretises the continuous-adjoint equations.

```@docs
lid_driven_cavity_2d_base
LidDrivenCavity2DFlow
```
