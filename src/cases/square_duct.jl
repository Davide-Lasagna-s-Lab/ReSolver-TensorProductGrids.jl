# =============================================================================================== #
# Pressure-driven square-duct case                                                               #
# =============================================================================================== #
#
# The grid layout lives in `grids/square_duct.jl`. This file selects the three-component Cartesian
# equations and applies the pressure-gradient force in the physical streamwise component `w`.
#
#     g = SquareDuctGrid(49, 63; width=7)
#     equations = SquareDuctFlow(g, 2000; f=4, fftw_flags=FFTW.ESTIMATE)
#
# Pass another `base_flow=(U, V, W)` to change the reference profile without changing the grid.

# =============================================================================================== #
# Equation constructor                                                                           #
# =============================================================================================== #

@doc raw"""
    SquareDuctFlow(g::AbstractSquareDuctGrid, Re::Real;
                   base_flow=(nothing, nothing, nothing), f::Real=1,
                   mode=AdjointDiscrete(), fftw_flags=FFTW.EXHAUSTIVE,
                   dealias=true) -> ProjectedNSE

Construct the pressure-driven incompressible square-duct equations

```math
\partial_t \boldsymbol{u} + (\boldsymbol{u}\cdot\nabla)\boldsymbol{u}
= -\nabla p + Re^{-1}\nabla^2\boldsymbol{u} + f\,\boldsymbol{e}_z,
\qquad \nabla\cdot\boldsymbol{u}=0,
```

where ``\boldsymbol{u}=(u,v,w)`` and physical ``z`` is streamwise. Writing the positive dimensional
pressure-gradient acceleration as ``G=-(1/\rho)d\bar p^*/dz^*``, and using the full square side
``L`` and velocity ``U_{ref}`` as reference scales,

```math
Re=\frac{U_{ref}L}{\nu},\qquad f=\frac{GL}{U_{ref}^2}.
```

`base_flow=(U,V,W)` is added only to the steady zero `(z,t)` Fourier mode. To evaluate or linearise
about a streamwise profile, pass its cross-section values as `(nothing,nothing,W)`. Every nonzero
component must have bounded shape `(N,N)`. The grid and equation constructor do not impose wall
values; the basis or residual formulation must enforce the required homogeneous conditions.

# Friction scaling

For a full square cross-section of side ``L``, the perimeter-averaged wall stress and conventional
friction Reynolds number are

```math
\bar\tau_w=\rho G\frac{L}{4},\qquad
u_\tau=\sqrt{\bar\tau_w/\rho},\qquad
Re_\tau=\frac{u_\tau(L/2)}{\nu}=\frac{Re\sqrt{|f|}}{4}.
```

Thus friction-velocity and half-side scaling is represented on the unit-side grid by
`SquareDuctFlow(g, 2Reτ; f=4)`. In contrast, `f=1` uses the pressure-gradient velocity ``\sqrt{GL}``;
it is not the conventional perimeter-averaged friction velocity.

# Arguments

- `g`: square-duct grid stored as `(x,y,z,t)`, with physical ``z`` streamwise.
- `Re`: real side-length Reynolds number ``U_{ref}L/\nu``, multiplying viscosity as ``1/Re``.

# Keyword arguments

- `base_flow`: three-component cross-section tuple added to the steady zero Fourier mode.
- `f`: nonzero signed uniform forcing amplitude in the physical streamwise component ``w``. Zero
  throws an `ArgumentError`.
- `mode`: `AdjointDiscrete()` or `AdjointContinuous()` for the linearised adjoint operator.
- `fftw_flags`: FFTW planner flags forwarded to `NSEBase.construct_equations`.
- `dealias`: whether nonlinear products use padded streamwise and temporal Fourier resolutions.

# Returns

An `NSEBase.ProjectedNSE` for three velocity components with constant streamwise forcing.

# Example

```julia
g = SquareDuctGrid(49, 63; Nt=1, α=0.5, width=7)
equations = SquareDuctFlow(g, 2000; f=4, fftw_flags=FFTW.ESTIMATE)
```
"""
function SquareDuctFlow(g::AbstractSquareDuctGrid, Re::Real;
                        base_flow=(nothing, nothing, nothing), f::Real=1,
                        mode=AdjointDiscrete(), fftw_flags=FFTW.EXHAUSTIVE,
                        dealias::Bool=true)
    _validate_equal_bounded_lengths(g, Val(2), "square duct")
    _validate_base_flow(g, base_flow, Val(3), "square duct")
    force = ConstantBodyForce(eltype(g)(f); i=3)
    return construct_equations(g, Re, base_flow, CartesianPrimitive3D();
                               force, mode, flags=fftw_flags, dealias)
end
