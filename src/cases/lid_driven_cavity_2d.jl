# =============================================================================================== #
# Two-dimensional lid-driven-cavity case                                                         #
# =============================================================================================== #
#
# The grid layout lives in `grids/lid_driven_cavity_2d.jl`. This file provides a canonical smooth
# divergence-free lifting and the matching two-component equation factory.
#
#     g = LidDrivenCavity2DGrid(65; width=7)
#     equations = LidDrivenCavity2DFlow(g, 1000; fftw_flags=FFTW.ESTIMATE)
#
# Pass another `base_flow=(U, V)` to change the moving-wall data without changing the grid.

# =============================================================================================== #
# Canonical moving-lid lifting                                                                   #
# =============================================================================================== #

@doc raw"""
    lid_driven_cavity_2d_base(g::AbstractLidDrivenCavity2DGrid) -> Tuple

Return the canonical smooth divergence-free moving-lid lifting `(U,V)` for a two-dimensional
cavity.

Let `ξ` and `η` be the physical `x` and `y` coordinates normalized to `[0,1]`. Define

```math
F(\xi)=16\xi^2(1-\xi)^2,\qquad H(\eta)=\eta^2(\eta-1).
```

The returned components are

```math
U=F(\xi)H'(\eta),\qquad
V=-F'(\xi)H(\eta).
```

Because the cavity has equal side lengths, these components satisfy
`\partial_xU+\partial_yV=0`. At the upper wall, `U=F` and `V=0`; both components vanish on the
other walls. The polynomial lid has unit maximum speed and tapers to zero at the upper corners,
avoiding the discontinuities of an ideal uniform lid.

# Arguments

- `g`: square two-dimensional lid-driven-cavity grid; shifted or uniformly scaled intervals are
  supported.

# Returns

A newly allocated bounded-grid tuple `(U,V)` with size `(N,N)` in each nonzero component.

# Example

```julia
g = LidDrivenCavity2DGrid(65; lim=(-1, 1), width=7)
U, V = lid_driven_cavity_2d_base(g)
```
"""
function lid_driven_cavity_2d_base(g::AbstractLidDrivenCavity2DGrid)
    _validate_equal_bounded_lengths(g, Val(2), "2D lid-driven cavity")
    return _lid_driven_cavity_xy_base(g)
end

# =============================================================================================== #
# Equation constructor                                                                           #
# =============================================================================================== #

@doc raw"""
    LidDrivenCavity2DFlow(g::AbstractLidDrivenCavity2DGrid, Re::Real;
                          base_flow=lid_driven_cavity_2d_base(g),
                          mode=AdjointDiscrete(), fftw_flags=FFTW.EXHAUSTIVE,
                          dealias=true) -> ProjectedNSE

Construct the two-dimensional incompressible lid-driven-cavity equations

```math
\partial_t \boldsymbol{u} + (\boldsymbol{u}\cdot\nabla)\boldsymbol{u}
= -\nabla p + Re^{-1}\nabla^2\boldsymbol{u},
\qquad \nabla\cdot\boldsymbol{u}=0,
```

where `\boldsymbol{u}=(u,v)` on a square cavity. With dimensional lid speed `U_{lid}`, side length
`L`, and kinematic viscosity `ν`,

```math
Re=\frac{U_{lid}L}{\nu}.
```

`base_flow=(U,V)` lifts steady moving-wall data into the zero temporal Fourier mode. The default
is [`lid_driven_cavity_2d_base`](@ref); either component may instead be `nothing` when identically
zero. Every nonzero component must have the bounded shape `(N,N)`.

The grid and equation constructor do not impose boundary values. Perturbations must satisfy
homogeneous wall conditions through the basis or residual formulation.

# Arguments

- `g`: square two-dimensional cavity grid stored as `(x,y,t)`.
- `Re`: real Reynolds number `U_{lid}L/ν`, multiplying viscosity as `1/Re`.

# Keyword arguments

- `base_flow`: two-component bounded-grid tuple added to the steady temporal mode.
- `mode`: `AdjointDiscrete()` or `AdjointContinuous()` for the linearised adjoint operator.
- `fftw_flags`: FFTW planner flags forwarded to `NSEBase.construct_equations`.
- `dealias`: whether nonlinear products use a padded temporal Fourier resolution.

# Returns

An `NSEBase.ProjectedNSE` for two velocity components.

# Example

```julia
g = LidDrivenCavity2DGrid(65; width=7)
equations = LidDrivenCavity2DFlow(g, 1000; fftw_flags=FFTW.ESTIMATE)
```
"""
function LidDrivenCavity2DFlow(g::AbstractLidDrivenCavity2DGrid, Re::Real;
                               base_flow=lid_driven_cavity_2d_base(g),
                               mode=AdjointDiscrete(), fftw_flags=FFTW.EXHAUSTIVE,
                               dealias::Bool=true)
    _validate_equal_bounded_lengths(g, Val(2), "2D lid-driven cavity")
    _validate_base_flow(g, base_flow, Val(2), "2D cavity")
    return construct_equations(g, Re, base_flow, CartesianPrimitive2D();
                               force=NoForce(), mode, flags=fftw_flags, dealias)
end
