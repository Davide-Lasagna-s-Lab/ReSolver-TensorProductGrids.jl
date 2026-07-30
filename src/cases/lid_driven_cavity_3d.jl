# =============================================================================================== #
# Three-dimensional lid-driven-cavity case                                                       #
# =============================================================================================== #
#
# The grid layout lives in `grids/lid_driven_cavity_3d.jl`. This file provides a canonical smooth
# divergence-free lifting and the matching three-component equation factory.
#
#     g = LidDrivenCavity3DGrid(49; width=7)
#     equations = LidDrivenCavity3DFlow(g, 1000; fftw_flags=FFTW.ESTIMATE)
#
# Pass another `base_flow=(U, V, W)` to change the moving-wall data without changing the grid.

# =============================================================================================== #
# Canonical moving-lid lifting                                                                   #
# =============================================================================================== #

@doc raw"""
    lid_driven_cavity_3d_base(g::AbstractLidDrivenCavity3DGrid) -> Tuple

Return the canonical smooth divergence-free moving-lid lifting ``(U,V,0)`` for a three-dimensional
cubic cavity. The returned tuple represents its zero third component by `nothing`.

Let ``\xi``, ``\eta``, and ``\zeta`` be physical ``x``, ``y``, and ``z`` normalized to ``[0,1]``,
and define

```math
F(ξ)=16ξ^2(1-ξ)^2,\qquad H(η)=η^2(η-1),\qquad
G(ζ)=16ζ^2(1-ζ)^2.
```

The returned components are

```math
U=F(ξ)H'(η)G(ζ),\qquad
V=-F'(ξ)H(η)G(ζ),\qquad W=0.
```

They satisfy ``\partial_xU+\partial_yV+\partial_zW=0``. At the upper ``y`` wall, the streamwise
velocity is ``U=FG``; all other velocity components and wall values vanish. The lid has unit maximum
speed and tapers smoothly to zero at all four of its edges.

# Arguments

- `g`: cubic three-dimensional lid-driven-cavity grid; shifted or uniformly scaled intervals are
  supported.

# Returns

A newly allocated bounded-grid tuple `(U,V,nothing)`, with `U` and `V` of size `(N,N,N)`.

# Example

```julia
g = LidDrivenCavity3DGrid(49; lim=(-1, 1), width=7)
U, V, W = lid_driven_cavity_3d_base(g)
```
"""
function lid_driven_cavity_3d_base(g::AbstractLidDrivenCavity3DGrid)
    _validate_equal_bounded_lengths(g, Val(3), "3D lid-driven cavity")
    Uxy, Vxy = _lid_driven_cavity_xy_base(g)
    zdim = inhomogeneous_storage_dims(g)[3]
    ζ = _normalized_bounded_coordinate(g, zdim)
    G = @. 16ζ^2 * (1 - ζ)^2
    taper = reshape(G, 1, 1, length(G))
    U = reshape(Uxy, size(Uxy)..., 1) .* taper
    V = reshape(Vxy, size(Vxy)..., 1) .* taper
    return U, V, nothing
end

# =============================================================================================== #
# Equation constructor                                                                           #
# =============================================================================================== #

@doc raw"""
    LidDrivenCavity3DFlow(g::AbstractLidDrivenCavity3DGrid, Re::Real;
                          base_flow=lid_driven_cavity_3d_base(g),
                          mode=AdjointDiscrete(), fftw_flags=FFTW.EXHAUSTIVE,
                          dealias=true) -> ProjectedNSE

Construct the three-dimensional incompressible lid-driven-cavity equations

```math
\partial_t \boldsymbol{u} + (\boldsymbol{u}\cdot\nabla)\boldsymbol{u}
= -\nabla p + Re^{-1}\nabla^2\boldsymbol{u},
\qquad \nabla\cdot\boldsymbol{u}=0,
```

where ``\boldsymbol{u}=(u,v,w)`` on a cubic cavity. With dimensional lid speed ``U_{lid}``, side
length ``L``, and kinematic viscosity ``\nu``,

```math
Re=\frac{U_{lid}L}{\nu}.
```

`base_flow=(U,V,W)` lifts steady moving-wall data into the zero temporal Fourier mode. The default
is [`lid_driven_cavity_3d_base`](@ref); any component may instead be `nothing` when identically
zero. Every nonzero component must have the bounded shape `(N,N,N)`.

The grid and equation constructor do not impose boundary values. Perturbations must satisfy
homogeneous wall conditions through the basis or residual formulation.

# Arguments

- `g`: cubic three-dimensional cavity grid stored as `(x,y,z,t)`.
- `Re`: real Reynolds number ``U_{lid}L/\nu``, multiplying viscosity as ``1/Re``.

# Keyword arguments

- `base_flow`: three-component bounded-grid tuple added to the steady temporal mode.
- `mode`: `AdjointDiscrete()` or `AdjointContinuous()` for the linearised adjoint operator.
- `fftw_flags`: FFTW planner flags forwarded to `NSEBase.construct_equations`.
- `dealias`: whether nonlinear products use a padded temporal Fourier resolution.

# Returns

An `NSEBase.ProjectedNSE` for three velocity components.

# Example

```julia
g = LidDrivenCavity3DGrid(49; width=7)
equations = LidDrivenCavity3DFlow(g, 1000; fftw_flags=FFTW.ESTIMATE)
```
"""
function LidDrivenCavity3DFlow(g::AbstractLidDrivenCavity3DGrid, Re::Real;
                               base_flow=lid_driven_cavity_3d_base(g),
                               mode=AdjointDiscrete(), fftw_flags=FFTW.EXHAUSTIVE,
                               dealias::Bool=true)
    _validate_equal_bounded_lengths(g, Val(3), "3D lid-driven cavity")
    _validate_base_flow(g, base_flow, Val(3), "3D cavity")
    return construct_equations(g, Re, base_flow, CartesianPrimitive3D();
                               force=NoForce(), mode, flags=fftw_flags, dealias)
end
