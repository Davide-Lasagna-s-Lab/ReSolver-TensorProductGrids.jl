# Two-dimensional lid-driven-cavity case

[`LidDrivenCavity2DFlow`](@ref) selects NSEBase's two-component Cartesian primitive formulation on
a fully bounded `(x,y)` domain. For velocity ``\boldsymbol{u}=(u,v)``, it represents

```math
\partial_t\boldsymbol{u}+(\boldsymbol{u}\boldsymbol{\cdot}\nabla)\boldsymbol{u}
=-\nabla p+\frac{1}{Re}\nabla^2\boldsymbol{u},
\qquad \nabla\boldsymbol{\cdot}\boldsymbol{u}=0.
```

For a square cavity of side ``L``, lid speed ``U_{lid}``, and kinematic viscosity ``\nu``, the
canonical Reynolds number is ``Re=U_{lid}L/\nu``. More general intervals and reference scales are
valid provided `Re`, the grid coordinates, and the supplied lifting use the same
nondimensionalisation.

## A smooth moving-lid lifting

The required `base_flow=(U,V)` supplies a steady lifting of the inhomogeneous wall data. On the
unit square, the following field is divergence free:

```math
\begin{aligned}
U(x,y)&=16x^2(1-x)^2(3y^2-2y),\\
V(x,y)&=-32x(1-x)(1-2x)y^2(y-1).
\end{aligned}
```

At the upper wall it gives ``U(x,1)=16x^2(1-x)^2`` and ``V(x,1)=0``; it vanishes on the other
walls and regularises the lid velocity to zero at the two upper corners. This avoids encoding the
corner discontinuities of an ideal uniform lid in the lifting. It is boundary data, not a claim
that the field solves the steady Navier--Stokes equations.

```@example cavity_2d_case
using NSEBase
using ReSolverRectangularGrids

g = LidDrivenCavity2DGrid(9, 9; Nt=1, width=3)
X, Y, _ = points(g)
U = @. 16X^2 * (1 - X)^2 * (3Y^2 - 2Y)
V = @. -32X * (1 - X) * (1 - 2X) * Y^2 * (Y - 1)

equations = LidDrivenCavity2DFlow(g, 1000; base_flow=(U, V),
                                  fftw_flags=FFTW.ESTIMATE, dealias=false)

@assert maximum(abs, U[:, end] .- 16 .* X[:, end].^2 .* (1 .- X[:, end]).^2) < 100eps()
@assert equations.base === (U, V)
(grid_size=size(g), operator=typeof(equations.nl))
```

## What the constructor does—and does not do

The lifting is inserted only in the zero temporal Fourier mode. Either entry may be `nothing` if
that component is identically zero. The grid and equation constructor do not impose no-slip or
moving-wall values: perturbation fields must satisfy homogeneous boundary conditions through the
basis or residual formulation. Adding the lifting to such perturbations recovers the prescribed
total-velocity boundary data.

The only homogeneous direction is time or phase. `dealias=true` pads that temporal Fourier
direction during nonlinear products; it does not change either bounded spatial resolution.
`mode=AdjointDiscrete()` gives the quadrature-consistent discrete adjoint of the linearised
operator, whereas `AdjointContinuous()` discretises the continuous-adjoint equations.

```@docs
LidDrivenCavity2DFlow
```
