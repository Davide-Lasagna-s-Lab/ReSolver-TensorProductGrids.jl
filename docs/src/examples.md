# Worked examples

The repository contains one standalone script for each layout. They deliberately use the same sequence:

1. construct a modest grid;
2. define a separable analytical field with `exp(cos(·))` in every active homogeneous direction;
3. compare a Fourier derivative with its exact value;
4. compare the NSEBase norm with an analytical integral;
5. construct the matching incompressible equations.

Because `exp(cos\theta)` has nonzero coefficients at all integer wavenumbers, these examples exercise the full represented Fourier band. The bounded factors are low-degree polynomials whose derivatives and quadrature integrals are exact to the reported tolerance.

## Run all examples

```@example worked-examples
using ReSolverRectangularGrids

root = joinpath(pkgdir(ReSolverRectangularGrids), "examples")
channel = include(joinpath(root, "channel.jl"))
cavity2d = include(joinpath(root, "lid_driven_cavity_2d.jl"))
cavity3d = include(joinpath(root, "lid_driven_cavity_3d.jl"))
duct = include(joinpath(root, "square_duct.jl"))

[(case=:channel, size=channel.grid_size, derivative_error=channel.derivative_error, norm_error=channel.norm_error),
 (case=:cavity2d, size=cavity2d.grid_size, derivative_error=cavity2d.derivative_error, norm_error=cavity2d.norm_error),
 (case=:cavity3d, size=cavity3d.grid_size, derivative_error=cavity3d.derivative_error, norm_error=cavity3d.norm_error),
 (case=:square_duct, size=duct.grid_size, derivative_error=duct.derivative_error, norm_error=duct.norm_error)]
```

The scripts are:

- [`examples/channel.jl`](https://github.com/Davide-Lasagna-s-Lab/ReSolver-RectangularGrids.jl/blob/main/examples/channel.jl)
- [`examples/lid_driven_cavity_2d.jl`](https://github.com/Davide-Lasagna-s-Lab/ReSolver-RectangularGrids.jl/blob/main/examples/lid_driven_cavity_2d.jl)
- [`examples/lid_driven_cavity_3d.jl`](https://github.com/Davide-Lasagna-s-Lab/ReSolver-RectangularGrids.jl/blob/main/examples/lid_driven_cavity_3d.jl)
- [`examples/square_duct.jl`](https://github.com/Davide-Lasagna-s-Lab/ReSolver-RectangularGrids.jl/blob/main/examples/square_duct.jl)

Every script is executed by `Pkg.test()` as well as this documentation page. An example therefore cannot drift away from the current constructors unnoticed.

## Analytical norm

The periodic factor has mean square

\[
\frac{1}{2\pi}\int_0^{2\pi}e^{2\cos\theta}\,d\theta
=I_0(2)
=\sum_{k=0}^{\infty}\frac{1}{(k!)^2}.
\]

The test scripts evaluate the rapidly convergent series without adding a special-function dependency. The bounded polynomial `b(x)=(x-a)(b-x)` satisfies

\[
\int_a^b b(x)^2\,dx=\frac{(b-a)^5}{30}.
\]

Tensor-product separability then gives an exact norm for every layout, matching the lazy product quadrature used in the implementation.
