# Coordinates and numerical conventions

This page defines the conventions shared by every bundled grid. They are part of the public contract: downstream bases, residuals, and analysis code can dispatch on them without inspecting concrete fields.

## Physical coordinates and storage order

Every grid has up to four physical coordinate slots `(x,y,z,t)`. `AXES` maps those slots to storage dimensions, using `nothing` for an absent coordinate. `points(g)` returns broadcast-compatible arrays in storage order, not necessarily physical order.

| Grid | `AXES` | Result of `points(g)` |
|:--|:--|:--|
| Channel | `(2,1,3,4)` | `(y,x,z,t)` |
| 2D cavity | `(1,2,nothing,3)` | `(x,y,t)` |
| 3D cubic cavity | `(1,2,3,4)` | `(x,y,z,t)` |
| Square duct | `(1,2,3,4)` | `(x,y,z,t)` |

Physical derivative wrappers use `AXES`, so `ddx!` always means physical `x`. The unusual channel storage order is therefore an optimization detail rather than a semantic change.

## Fourier order and wavenumbers

`fft_storage_dims(g)` lists transformed storage dimensions. The first is a real-to-complex transform, which stores only nonnegative wavenumbers; later directions are complex transforms with signed modes. Real physical fields recover the omitted modes through Hermitian symmetry.

For an integer mode ``n`` and scale ``\kappa``,

```math
k_{\mathrm{physical}} = n\kappa,\qquad L=\frac{2\pi}{\kappa}.
```

Spatial scales are passed as constructor keywords such as ``\alpha=2\pi/L_x`` and
``\beta=2\pi/L_z``. Every bundled temporal coordinate is a phase ``s\in[0,2\pi)`` with
scale one, so temporal mode ``h`` has derivative ``ih``. A physical or nondimensional period ``P``
is introduced in the governing residual through the fundamental frequency ``\omega=2\pi/P``, as in
ReSolver's residual ``\omega\,\partial_s u - N(u)``.

All Fourier sizes are odd to avoid an unpaired Nyquist mode. `Nt=1` leaves only the temporal mean and represents a steady field.

## Transforms and normalization

`FFT(Field(g,...))` divides the forward transform by the product of the Fourier resolutions; `IFFT` applies the inverse convention. Parseval's identity therefore makes NSEBase's `dot` and `norm` apply quadrature in bounded coordinates and average periodic coordinates:

```math
\langle u,v\rangle =
\sum_{i_1,\ldots,i_N}
\left(\prod_{d=1}^{N}w^{(d)}_{i_d}\right)
\frac{1}{|\Omega_{\mathrm{Fourier}}|}
\int_{\Omega_{\mathrm{Fourier}}}
\operatorname{Re}(\overline{u}v)\,d\boldsymbol{\xi}.
```

In spectral storage the first real-to-complex direction receives Hermitian multiplicity one at the mean mode and two at positive modes. `VectorField` inner products sum this scalar product over velocity components.

## Product quadrature

Each bounded direction owns a one-dimensional FDGrids weight vector. `RectangularProductWeights` exposes their tensor product lazily:

```math
W_{i_1,\ldots,i_N}=\prod_{d=1}^{N}w^{(d)}_{i_d}.
```

Its `AbstractArray` interface lets NSEBase use ordinary Cartesian indexing while avoiding a dense multidimensional allocation. The square and cubic cavity constructors share one bounded discretisation between every spatial direction; the square-duct constructor does the same between cross-stream `x` and `y`. This literal sharing preserves coordinate symmetry and avoids rebuilding identical operators.

## Discrete and continuous adjoints

Let ``D`` be a finite-difference matrix and ``W=\operatorname{diag}(w)`` the positive quadrature
matrix. The discrete inner product is

```math
\langle u,v\rangle_W = u^\ast Wv.
```

The matrix ``D^+`` used by `AdjointDiscrete()` is defined by

```math
D^+ = W^{-1}D^\ast W.
```

It follows directly that

```math
\langle Du,v\rangle_W
=(Du)^\ast Wv
=u^\ast D^\ast Wv
=u^\ast W(W^{-1}D^\ast W)v
=\langle u,D^+v\rangle_W.
```

FDGrids constructs and stores `D₁⁺` and `D₂⁺` from the same quadrature weights as the forward matrices. In multiple bounded directions, both the operator and weights are tensor products. Identity factors and weights in every untouched direction cancel, so the one-dimensional identity above holds for each Cartesian derivative and their Laplacian sum.

For a Fourier first derivative ``D_k=\mathrm{i}n\kappa``, the averaged periodic inner product gives
``D_k^+=-D_k``. NSEBase implements this sign reversal for `AdjointDiscrete()` and applies the
real-FFT Hermitian multiplicities in the inner product. Fourier second derivatives
``-(n\kappa)^2`` are self-adjoint.

At low level, `ddx!`, `ddy!`, `ddz!`, `ddt!`, and `laplacian!` accept only `Forward()` and `AdjointDiscrete()`. `AdjointContinuous()` is an equation-level choice: NSEBase assembles the continuous PDE adjoint analytically from forward derivatives. The continuous and discrete adjoints need not coincide after discretization, especially at boundaries.

!!! important "Boundary terms"
    The weighted matrix identity is exact for the stored arrays, but it does not choose physically admissible wall perturbations. Integration-by-parts boundary terms, pressure compatibility, and homogeneous perturbation conditions remain the responsibility of the basis, projection, or residual formulation.

## Dealiasing and growth

Quadratic nonlinear products can fold unresolved Fourier modes into the represented band. With
`dealias=true`, NSEBase pads every Fourier direction to the nearest odd size at least ``3N/2``,
evaluates products in physical space, and truncates symmetrically on return.

`points(g; dealias=true)` returns coordinates at those padded sizes. `growto(g, sizes)` creates a grid with explicit new Fourier resolutions while sharing every bounded numerical object. This supports spectral interpolation and truncation without rebuilding FDGrids matrices.

## Base-flow lifting

Case constructors pass a tuple `base_flow` to NSEBase. Each non-`nothing` entry contains one Cartesian velocity component over the bounded coordinates and is added only to the zero Fourier mode. This separates prescribed steady wall data or a laminar reference profile from homogeneous perturbations.

The cavity factories default to [`lid_driven_cavity_2d_base`](@ref) and
[`lid_driven_cavity_3d_base`](@ref). Both helpers normalize the common cavity interval `lim` to
``[0,1]``, produce smooth divergence-free liftings, and taper the moving-lid velocity to zero at
its corners or edges. Supplying `base_flow` explicitly replaces that canonical boundary-data
lifting without changing any other constructor keyword.
