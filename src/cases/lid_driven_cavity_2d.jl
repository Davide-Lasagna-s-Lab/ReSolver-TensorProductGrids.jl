# =============================================================================================== #
# Two-dimensional lid-driven-cavity equations                                                     #
# =============================================================================================== #
#
# The grid layout lives in `grids/lid_driven_cavity_2d.jl`. This file contains only the matching
# two-component equation factory.

@doc raw"""
    LidDrivenCavity2DFlow(g::AbstractLidDrivenCavity2DGrid, Re;
                          base_flow, mode=AdjointDiscrete(),
                          fftw_flags=FFTW.EXHAUSTIVE,
                          dealias=true) -> ProjectedNSE

Construct the two-dimensional incompressible lid-driven-cavity equations at Reynolds number `Re`:

```math
\partial_t \boldsymbol{u} + (\boldsymbol{u}\cdot\nabla)\boldsymbol{u}
= -\nabla p + Re^{-1}\nabla^2\boldsymbol{u},
\qquad \nabla\cdot\boldsymbol{u}=0,
```

where `u = (u, v)` on the bounded `(x, y)` cavity domain.

With dimensional lid speed `U_lid`, reference cavity side length `L`, and kinematic viscosity `ν`,

```math
Re = \frac{U_{lid}L}{\nu}.
```

For a noncanonical scaling, replace `U_lid` and `L` by the velocity and length used to
nondimensionalise the supplied grid and `base_flow`.

The required `base_flow=(U, V)` tuple lifts the inhomogeneous moving-lid data into the zero temporal
Fourier mode. Either component may be `nothing` when it is identically zero. Perturbations must
satisfy homogeneous wall conditions; the grid and equation constructor do not impose them.

`mode` selects the continuous or quadrature-consistent discrete adjoint. `fftw_flags` controls FFT
planning, and `dealias=true` pads the temporal Fourier direction used by nonlinear products.

# Arguments

- `g`: two-dimensional cavity grid stored as `(x,y,t)`.
- `Re`: Reynolds number `U_ref*L_ref/ν`, which multiplies viscosity as `1/Re`.

# Keyword arguments

- `base_flow`: required two-tuple `(U,V)` lifting the prescribed wall data into the steady temporal
  mode; an identically zero component may be `nothing`.
- `mode`: `AdjointDiscrete()` or `AdjointContinuous()` for the linearised adjoint operator.
- `fftw_flags`: FFTW planner flags forwarded to `NSEBase.construct_equations`.
- `dealias`: whether nonlinear products use a padded temporal Fourier resolution.

# Returns

An `NSEBase.ProjectedNSE` for two velocity components.

# Example

```julia
g = LidDrivenCavity2DGrid(65, 65; width=7)
X, Y, _ = points(g)
U = @. 16X^2 * (1 - X)^2 * (3Y^2 - 2Y)
V = @. -32X * (1 - X) * (1 - 2X) * Y^2 * (Y - 1)
equations = LidDrivenCavity2DFlow(g, 1000; base_flow=(U, V), fftw_flags=FFTW.ESTIMATE)
```
"""
function LidDrivenCavity2DFlow(g::AbstractLidDrivenCavity2DGrid, Re; base_flow,
                               mode=AdjointDiscrete(), fftw_flags=FFTW.EXHAUSTIVE, dealias=true)
    # Require one steady lifting field for each Cartesian velocity component.
    base_flow isa Tuple && length(base_flow) == 2 ||
        throw(ArgumentError("a 2D cavity base flow must contain (U, V)"))

    # Build the nonlinear and requested adjoint operators in `(u, v)` order.
    return construct_equations(g, Re, base_flow, CartesianPrimitive2D();
                               mode, flags=fftw_flags, dealias)
end
