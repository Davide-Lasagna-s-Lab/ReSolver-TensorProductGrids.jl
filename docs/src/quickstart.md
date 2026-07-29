# Quick start

This walkthrough constructs a channel, verifies differentiation and quadrature against analytical values, and builds a pressure-driven flow case.

## Grid and field

```@example quickstart
using LinearAlgebra
using NSEBase
using ReSolverRectangularGrids

g = ChannelGrid(19, 25, 19; Nt=1, α=0.5, β=1.25, width=7)
y, x, z, t = points(g)

u(y, x, z, t) = (1 - y^2) * exp(cos(0.5x)) * exp(cos(1.25z))
uz(y, x, z, t) = -1.25sin(1.25z) * u(y, x, z, t)

û = FFT(Field(g, u))
numerical = ddz!(FTField(g), û)
exact = FFT(Field(g, uz))
derivative_error = maximum(abs, parent(numerical) .- parent(exact)) / maximum(abs, parent(exact))

@assert derivative_error < 5e-7
derivative_error
```

The positional resolutions are in physical order `(Nx,Ny,Nz)`, but channel arrays are stored as `(y,x,z,t)` so the bounded coordinate is contiguous. Named derivative functions stay physical: `ddz!` differentiates physical `z` regardless of its storage dimension.

The streamwise and spanwise periods are `2π/α` and `2π/β`. `Nt=1` means that the field is steady. The default Gauss–Lobatto wall-normal grid includes both walls at `y=±1`.

## Analytical norm

NSEBase's spectral norm integrates the bounded coordinate using FDGrids quadrature and averages every periodic coordinate. For the field above,

\[
\|u\|^2 =
\int_{-1}^{1}(1-y^2)^2\,dy
\,
\left(\frac{1}{2\pi}\int_0^{2\pi}e^{2\cos\theta}\,d\theta\right)^2
= \frac{16}{15}I_0(2)^2.
\]

```@example quickstart
I₀₂ = sum(inv(float(factorial(k)))^2 for k in 0:20)
norm_error = abs(norm(û)^2 - (16 / 15) * I₀₂^2)

@assert norm_error < 5e-11
norm_error
```

The `exp(cos(·))` factors have nonzero coefficients at all Fourier wavenumbers, so the derivative checks spectral behavior across the full resolved band.

## Equations

For channel half-height `h`, reference velocity `U_ref`, and viscosity `ν`, the constructor uses `Re=U_ref h/ν`. With friction-velocity scaling, `Re=Reτ` and unit dimensionless forcing, the laminar equilibrium is `U=(Reτ/2)(1-y²)`:

```@example quickstart
Reτ = 180
U = (Reτ / 2) .* plane_poiseuille_base(g)
equations = PlanePoiseuilleFlow(g, Reτ; base_flow=(U, nothing, nothing), f=1,
                                fftw_flags=FFTW.ESTIMATE, dealias=false)

@assert equations isa ProjectedNSE
typeof(equations.nl)
```

`base_flow` is a steady zero-mode lifting. It supplies the total streamwise reference profile; perturbation modes should satisfy homogeneous no-slip conditions at the walls. `dealias=true` is recommended for production nonlinear evaluations. It is disabled here only to keep the documentation example small.

Continue with [Channel flows](@ref) for Couette rotation, Poiseuille pressure-gradient,
prescribed-bulk-flow, and Reynolds-number conventions.
