# Three-dimensional lid-driven-cavity case

[`LidDrivenCavity3DFlow`](@ref) selects NSEBase's three-component Cartesian primitive formulation
on a fully bounded `(x,y,z)` domain. For ``\boldsymbol{u}=(u,v,w)``, it represents

```math
\partial_t\boldsymbol{u}+(\boldsymbol{u}\boldsymbol{\cdot}\nabla)\boldsymbol{u}
=-\nabla p+\frac{1}{Re}\nabla^2\boldsymbol{u},
\qquad \nabla\boldsymbol{\cdot}\boldsymbol{u}=0.
```

For a cubic cavity of side ``L``, lid speed ``U_{lid}``, and kinematic viscosity ``\nu``, the
canonical Reynolds number is ``Re=U_{lid}L/\nu``. Rectangular boxes and other reference scales are
also supported when the grid, `Re`, and base-flow lifting use a consistent nondimensionalisation.

## A divergence-free, edge-regularised lifting

`base_flow=(U,V,W)` carries the steady moving-wall data. On the unit cube, define

```math
\begin{aligned}
G(z)&=16z^2(1-z)^2,\\
U(x,y,z)&=16x^2(1-x)^2(3y^2-2y)G(z),\\
V(x,y,z)&=-32x(1-x)(1-2x)y^2(y-1)G(z),\\
W(x,y,z)&=0.
\end{aligned}
```

Because ``\partial_xU+\partial_yV=0`` and ``W=0``, this lifting is divergence free. On the upper
wall, ``U(x,1,z)=16x^2(1-x)^2G(z)``; the polynomial factors taper the moving-wall velocity smoothly
to zero at all four lid edges. The field vanishes on the remaining walls. As in two dimensions,
the lifting encodes boundary data but is not asserted to be a steady solution.

```@example cavity_3d_case
using NSEBase
using ReSolverRectangularGrids

g = LidDrivenCavity3DGrid(7, 7, 7; Nt=1, width=3)
X, Y, Z, _ = points(g)
G = @. 16Z^2 * (1 - Z)^2
U = @. 16X^2 * (1 - X)^2 * (3Y^2 - 2Y) * G
V = @. -32X * (1 - X) * (1 - 2X) * Y^2 * (Y - 1) * G

equations = LidDrivenCavity3DFlow(g, 1000; base_flow=(U, V, nothing),
                                  fftw_flags=FFTW.ESTIMATE, dealias=false)

lid = @. 16X[:, end, :]^2 * (1 - X[:, end, :])^2 * G[:, end, :]
@assert iszero(U[:, 1, :]) && iszero(U[:, end, :] .- lid)
@assert equations.base === (U, V, nothing)
(grid_size=size(g), operator=typeof(equations.nl))
```

## Boundary and spectral semantics

The lifting is added only to the zero temporal Fourier mode, and an identically zero component may
be `nothing`. Neither [`LidDrivenCavity3DGrid`](@ref) nor the equation factory imposes boundary
values. Perturbation fields must satisfy homogeneous wall conditions through the basis or residual
formulation so that lifting plus perturbation has the desired total wall velocity.

Time or phase is the only Fourier direction. `dealias=true` pads this direction for nonlinear
products while leaving the three bounded spatial resolutions unchanged. `mode=AdjointDiscrete()`
constructs the quadrature-consistent discrete adjoint of the linearised operator;
`AdjointContinuous()` constructs the discretised continuous-adjoint equations.

```@docs
LidDrivenCavity3DFlow
```
