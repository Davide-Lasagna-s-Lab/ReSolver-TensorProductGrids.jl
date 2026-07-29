# Case forcings

Case constructors compose small callable forcing policies and pass them to
`NSEBase.construct_equations`. These policies add physical source terms without coupling them to a
particular grid layout. ReSolverRectangularGrids provides a linear Coriolis coupling and a
state-independent constant body force.

## Coriolis force

[`CoriolisForce`](@ref) couples the first two velocity components. For `Ro` and a three-component
state ``\boldsymbol{u}=(u,v,w)``, its forward action is

```math
\mathcal C\boldsymbol{u}
=Ro\begin{bmatrix}0&1&0\\-1&0&0\\0&0&0\end{bmatrix}\boldsymbol{u}
=Ro\,(v,-u,0).
```

The trailing row and column are absent for a two-component state. Since
``\mathcal C^T=-\mathcal C``, both continuous- and discrete-adjoint evaluations reverse the sign.
The discrete result is the same transpose because all velocity components carry the same spatial
quadrature. The `Ro` convention used by the channel cases is documented under
[`PlaneCouetteFlow`](@ref).

## Constant body force

[`ConstantBodyForce`](@ref) adds a scalar value to component `i` at the zero mode of every Fourier
direction. The selected slice retains every bounded-grid point, so this single spectral mode
represents a uniform physical-space force. For example, channel pressure forcing uses `i=1`, while
square-duct forcing uses `i=3` because physical streamwise velocity is `w` there. Its amplitude
must be nonzero: constructing `ConstantBodyForce(0)` throws an `ArgumentError`. Use `NoForce()`
when the equations have no body force instead of storing an inert constant-force policy.

If ``\mathcal F(\boldsymbol{u})=\boldsymbol{c}`` is independent of the state, then

```math
D\mathcal F[\boldsymbol{u}]\,\delta\boldsymbol{u}=0,
\qquad
\bigl(D\mathcal F[\boldsymbol{u}]\bigr)^*\boldsymbol{v}=0.
```

Accordingly, `ConstantBodyForce` changes only `Forward()` evaluations. Calls in
`AdjointContinuous()` or `AdjointDiscrete()` mode leave the output unchanged; adding the constant
to an adjoint residual would be mathematically incorrect.

```@example forcing_policies
using NSEBase
using ReSolverRectangularGrids

g = ChannelGrid(3, 9, 3; Nt=1, width=3)
u = VectorField(g; N=3)
out = VectorField(g; N=3)
parent(u[1]) .= 1
parent(u[2]) .= 2

rotation = CoriolisForce(0.25)
rotation(out, u, Forward())
@assert all(parent(out[1]) .== 0.5)
@assert all(parent(out[2]) .== -0.25)

out .= 0
pressure = ConstantBodyForce(2.5; i=1)
pressure(out, u, Forward())
@assert all(parent(out[1])[:, 1, 1, 1] .== 2.5)

out .= 0
pressure(out, u, AdjointDiscrete())
@assert iszero(parent(out[1]))
(rotation=typeof(rotation), pressure=typeof(pressure))
```

Forces accumulate into `out`; they do not clear it first. `ConstantBodyForce` rejects a zero
amplitude when it is constructed and validates `i` against the number of components when it is
applied. Use NSEBase's `CompoundForcing` when defining a custom case that requires multiple
independent forcing policies.

```@docs
CoriolisForce
ConstantBodyForce
```
