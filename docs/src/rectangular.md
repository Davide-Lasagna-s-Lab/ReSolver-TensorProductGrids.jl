# Rectangular-grid design

`RectangularGrid` is the common implementation behind every bundled layout. It combines one to three bounded FDGrids directions with one or more periodic Fourier directions and encodes all layout metadata in its concrete type.

Most users should call a case-specific constructor. The low-level object is documented here to make its behavior inspectable and its extension points precise.

## Type-level layout

An NSEBase grid maps the physical coordinate slots `(x,y,z,t)` to array dimensions through `AXES`. An absent physical coordinate is `nothing`. `FFT_DIMS` lists transformed storage dimensions in transform order; the first entry is stored with a real-to-complex transform.

For a channel,

```text
physical coordinate     x  y  z  t
storage dimension       2  1  3  4
AXES                    (2, 1, 3, 4)
FFT_DIMS                (2, 3, 4)
stored array order      (y, x, z, t)
```

This is why `points(g)` returns `(y,x,z,t)` for a channel even though named wrappers such as `ddx!` and `ddy!` always retain their physical meaning.

## Bounded directions

The tuples `xs`, `D₁`, `D₂`, `D₁⁺`, `D₂⁺`, and `ws` follow increasing inhomogeneous storage dimension. `D₁` and `D₂` are FDGrids first- and second-derivative matrices. Their plus variants are the quadrature-weighted discrete adjoints. Different bounded directions may use different intervals, point distributions, and stencil widths.

`weights(g)` returns a cached `RectangularProductWeights` object. For one-dimensional weights
``w_1,\ldots,w_N``,

```math
W[i_1,\ldots,i_N] = \prod_{j=1}^{N} w_j[i_j].
```

The product is evaluated when indexed; no dense tensor-product array is allocated.

## Fourier directions

`scales` follows `FFT_DIMS` order. An integer mode ``n`` in a direction with scale ``\kappa`` has
physical wavenumber ``n\kappa`` and period ``2\pi/\kappa``. Every scale is positive and finite, and
every Fourier resolution is positive and odd. The first Fourier dimension is real-to-complex;
subsequent dimensions retain signed complex modes.

`points(g; dealias=true)` returns coordinates at NSEBase's odd 3/2-rule padded Fourier sizes. `growto(g, sizes)` changes only the Fourier resolutions and preserves the bounded points, matrices, adjoints, weights, and scales.

## Derivative dispatch

NSEBase maps `ddx!`, `ddy!`, `ddz!`, and `ddt!` through `AXES`. A Fourier derivative multiplies
mode ``n`` by ``\mathrm{i}n\kappa``. A bounded derivative retrieves the correct FDGrids matrix through:

```julia
NSEBase.derivative_matrix(g, storage_dim, Val(order), mode)
```

`Forward()` selects `D₁` or `D₂`, while `AdjointDiscrete()` selects `D₁⁺` or `D₂⁺`. The derivative order is a `Val` because forward and adjoint FDGrids operators can have different concrete types; dispatch therefore remains type stable without a runtime matrix union.

The [adjoint section](@ref "Discrete and continuous adjoints") derives the exact weighted identity and explains why `AdjointContinuous()` belongs at equation level rather than in this matrix-lookup API.

## Accessor example

```@example rectangular-accessors
using NSEBase
using ReSolverRectangularGrids

g = ChannelGrid(7, 17, 5; Nt=3, α=1.25, β=0.75, width=5)
y, x, z, t = points(g)

@assert size(g) == (17, 7, 5, 3)
@assert size.(points(g)) == ((17, 1, 1, 1), (1, 7, 1, 1), (1, 1, 5, 1), (1, 1, 1, 3))
@assert wavenumber_scale(g, 2) == 1.25
@assert weights(g) === g.w

(storage_size=size(g), bounded_size=size(weights(g)), fourier_dimensions=fft_storage_dims(g))
```

## Low-level API

```@docs
RectangularGrid
RectangularProductWeights
```
