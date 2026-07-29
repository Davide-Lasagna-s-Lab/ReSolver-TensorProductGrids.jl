# =============================================================================================== #
# Pressure-driven square-duct equations                                                          #
# =============================================================================================== #
#
# The grid layout lives in `grids/square_duct.jl`. This file chooses the three-component Cartesian
# formulation and applies the pressure-gradient force in the physical streamwise component `w`.

@doc raw"""
    SquareDuctFlow(g::AbstractSquareDuctGrid, Re;
                   base_flow=(nothing, nothing, nothing), f=1,
                   mode=AdjointDiscrete(), fftw_flags=FFTW.EXHAUSTIVE,
                   dealias=true) -> ProjectedNSE

Construct the three-dimensional incompressible square-duct equations at Reynolds number `Re`:

```math
\partial_t \boldsymbol{u} + (\boldsymbol{u}\cdot\nabla)\boldsymbol{u}
= -\nabla p + Re^{-1}\nabla^2\boldsymbol{u} + f\,\boldsymbol{e}_z,
\qquad \nabla\cdot\boldsymbol{u}=0,
```

where the velocity state is `(u, v, w)` and physical `z` is streamwise.

Writing the positive dimensional pressure-gradient acceleration as
`G=-(1/ρ) d p̄*/d z*`, and using the square side length `L` and velocity `U_ref` as reference scales,

```math
Re = \frac{U_{ref}L}{\nu}, \qquad f = \frac{GL}{U_{ref}^2}.
```

`base_flow` must be a Cartesian tuple `(U, V, W)`; each entry may be `nothing`. To linearise about
a streamwise laminar solution, provide its cross-section values as `(nothing, nothing, W)`.
[`ConstantBodyForce`](@ref) applies `f` to component three at zero `(z, t)` wavenumber.

`mode` selects the continuous or quadrature-consistent discrete adjoint. `fftw_flags` controls FFT
planning, and `dealias=true` pads the streamwise and temporal Fourier directions.

# Friction scaling

For the full square cross-section of side `L`, the perimeter-averaged wall stress and conventional
friction Reynolds number are

```math
\bar\tau_w = \rho G\frac{L}{4}, \qquad
u_\tau = \sqrt{\bar\tau_w/\rho}, \qquad
Re_\tau = \frac{u_\tau(L/2)}{\nu} = \frac{Re\sqrt{|f|}}{4}.
```

Thus friction-velocity and half-side scaling is represented on this unit-side grid by
`SquareDuctFlow(g, 2Reτ; f=4)`. In contrast, `f=1` uses the pressure-gradient velocity
`sqrt(G*L)`; it must not be described as the conventional perimeter-averaged friction velocity.

# Arguments

- `g`: square-duct grid stored as `(x,y,z,t)`, with physical `z` streamwise.
- `Re`: side-length Reynolds number `U_ref*L/ν`, which multiplies viscosity as `1/Re`.

# Keyword arguments

- `base_flow`: three-tuple `(U,V,W)` added to the steady zero `(z,t)` Fourier mode; an identically
  zero component may be `nothing`.
- `f`: signed uniform forcing amplitude in the physical streamwise component `w`.
- `mode`: `AdjointDiscrete()` or `AdjointContinuous()` for the linearised adjoint operator.
- `fftw_flags`: FFTW planner flags forwarded to `NSEBase.construct_equations`.
- `dealias`: whether nonlinear products use padded streamwise and temporal Fourier resolutions.

# Returns

An `NSEBase.ProjectedNSE` for three velocity components with streamwise mean forcing.

# Example

```julia
g = SquareDuctGrid(49, 63; Nt=1, α=0.5, width=7)
equations = SquareDuctFlow(g, 2000; f=4, fftw_flags=FFTW.ESTIMATE)
```
"""
function SquareDuctFlow(g::AbstractSquareDuctGrid, Re;
                        base_flow=(nothing, nothing, nothing), f=1,
                        mode=AdjointDiscrete(), fftw_flags=FFTW.EXHAUSTIVE, dealias=true)
    # Require one steady lifting field for each Cartesian velocity component.
    base_flow isa Tuple && length(base_flow) == 3 ||
        throw(ArgumentError("a square-duct base flow must contain (U, V, W)"))

    # Apply the uniform pressure-gradient acceleration in physical streamwise `z`, component three.
    force = ConstantBodyForce(eltype(g)(f); i=3)

    # Build the nonlinear and requested adjoint operators in `(u, v, w)` order.
    return construct_equations(g, Re, base_flow, CartesianPrimitive3D();
                               force, mode, flags=fftw_flags, dealias)
end
