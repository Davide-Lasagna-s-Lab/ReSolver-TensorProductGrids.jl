# =============================================================================================== #
# Three-dimensional lid-driven-cavity equations                                                  #
# =============================================================================================== #
#
# The fully bounded grid layout lives in `grids/lid_driven_cavity_3d.jl`. This file selects
# NSEBase's three-component Cartesian formulation and records the supplied boundary-data lifting.

@doc raw"""
    LidDrivenCavity3DFlow(g::AbstractLidDrivenCavity3DGrid, Re;
                          base_flow, mode=AdjointDiscrete(),
                          fftw_flags=FFTW.EXHAUSTIVE,
                          dealias=true) -> ProjectedNSE

Construct the three-dimensional incompressible lid-driven-cavity equations at Reynolds number
`Re`:

```math
\partial_t \boldsymbol{u} + (\boldsymbol{u}\cdot\nabla)\boldsymbol{u}
= -\nabla p + Re^{-1}\nabla^2\boldsymbol{u},
\qquad \nabla\cdot\boldsymbol{u}=0,
```

where `u = (u, v, w)` on the fully bounded `(x, y, z)` cavity domain.

With dimensional lid speed `U_lid`, reference cavity side length `L`, and kinematic viscosity `ν`,

```math
Re = \frac{U_{lid}L}{\nu}.
```

For a noncanonical scaling, replace `U_lid` and `L` by the velocity and length used to
nondimensionalise the supplied grid and `base_flow`.

The required `base_flow=(U, V, W)` tuple lifts inhomogeneous moving-wall data into the zero temporal
Fourier mode. Any identically zero component may be `nothing`. Perturbations must satisfy
homogeneous wall conditions; neither the grid nor this constructor imposes boundary values.

`mode` selects the continuous or quadrature-consistent discrete adjoint. `fftw_flags` controls FFT
planning, and `dealias=true` pads the temporal Fourier direction used by nonlinear products.

# Arguments

- `g`: fully bounded three-dimensional cavity grid stored as `(x,y,z,t)`.
- `Re`: Reynolds number `U_ref*L_ref/ν`, which multiplies viscosity as `1/Re`.

# Keyword arguments

- `base_flow`: required three-tuple `(U,V,W)` lifting prescribed wall data into the steady temporal
  mode; an identically zero component may be `nothing`.
- `mode`: `AdjointDiscrete()` or `AdjointContinuous()` for the linearised adjoint operator.
- `fftw_flags`: FFTW planner flags forwarded to `NSEBase.construct_equations`.
- `dealias`: whether nonlinear products use a padded temporal Fourier resolution.

# Returns

An `NSEBase.ProjectedNSE` for three velocity components.

# Example

The following divergence-free lifting regularizes the streamwise lid velocity at all four lid
edges:

```julia
g = LidDrivenCavity3DGrid(49, 49, 49; width=7)
X, Y, Z, _ = points(g)
G = @. 16Z^2 * (1 - Z)^2
U = @. 16X^2 * (1 - X)^2 * (3Y^2 - 2Y) * G
V = @. -32X * (1 - X) * (1 - 2X) * Y^2 * (Y - 1) * G
equations = LidDrivenCavity3DFlow(g, 1000; base_flow=(U, V, nothing),
                                  fftw_flags=FFTW.ESTIMATE)
```
"""
function LidDrivenCavity3DFlow(g::AbstractLidDrivenCavity3DGrid, Re; base_flow,
                               mode=AdjointDiscrete(), fftw_flags=FFTW.EXHAUSTIVE, dealias=true)
    # Require one steady lifting field for each Cartesian velocity component.
    base_flow isa Tuple && length(base_flow) == 3 ||
        throw(ArgumentError("a 3D cavity base flow must contain (U, V, W)"))

    # Build the nonlinear and requested adjoint operators in `(u, v, w)` order.
    return construct_equations(g, Re, base_flow, CartesianPrimitive3D();
                               mode, flags=fftw_flags, dealias)
end
