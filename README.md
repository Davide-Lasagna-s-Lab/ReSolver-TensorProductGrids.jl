# ReSolverTensorProductGrids.jl

A tensor-product grid for [ReSolverFlowsBase](https://github.com/Davide-Lasagna-s-Lab/ReSolver-FlowsBase.jl):
one to three bounded directions, each given by its points, derivative matrices and quadrature
weights, times one or more periodic Fourier directions.

## The grid

A bounded direction is described by its collocation points `x`, first- and second-derivative
matrices `D₁`, `D₂`, their discrete adjoints `D₁⁺`, `D₂⁺`, and quadrature weights `w`. A Fourier
direction of period `L` is described by its wavenumber scale `2π/L`: the integer mode `n` has
wavenumber `n 2π/L`. The grid is built from these data, one tuple entry per bounded direction:

```julia
using FDGrids, ReSolverTensorProductGrids

# wall-normal direction of a channel: Ny Chebyshev–Lobatto points on [-1, 1]
fd = FDGrids.grid(Ny, -1, 1, GaussLobattoGrid())
y, w = fd.xs, fd.ws
D₁, D₂ = DiffMatrix(y, 7, 1), DiffMatrix(y, 7, 2)

g = TensorProductGrid((y,), (D₁,), (D₂,), (adjoint(D₁, w),), (adjoint(D₂, w),), (w,),
                      (2π/Lx, 2π/Lz, 1),   # wavenumber scales of x, z, t
                      (Ny, Nx, Nz, Nt),    # array size, in storage order
                      (2, 1, 3, 4),        # storage dimension of the coordinates (x1, x2, x3, t)
                      (2, 3, 4))           # Fourier dimensions, the real-to-complex one first
```

The derivative matrices are typically [FDGrids](https://github.com/Davide-Lasagna-s-Lab/FDGrids.jl)
operators, banded and multiplied dimension-wise; any matrix supporting that multiplication works.

## Quadrature weights

Integrals over the bounded directions use the tensor product of the one-dimensional weights,

```
∫∫ f(x, y) dx dy  ≈  Σᵢ Σⱼ wₓ[i] w_y[j] f(xᵢ, yⱼ),
```

returned lazily, without allocation, by `weights(g)` as a `TensorProductWeights` array. The weights
carry any measure of the coordinates; for a radial direction they include `r dr`. The discrete
adjoints are taken with respect to these weights, `D⁺ = W⁻¹DᵀW`, so that
`⟨w, D u⟩ = ⟨D⁺ w, u⟩` exactly. ReSolverFlowsBase combines these weights with an average over the
Fourier directions in its inner products and norms.

## The `AbstractGrid` interface

| ReSolverFlowsBase asks for | `TensorProductGrid` returns |
|---|---|
| `size(g)` | the array size, in storage order |
| `points(g; dealias=false)` | the collocation points in bounded directions, equispaced points on `[0, L)` in Fourier directions, shaped for broadcasting |
| `wavenumber_scale(g, dim)` | `2π/L` for a Fourier direction, one for a bounded one |
| `weights(g)` | the tensor-product quadrature weights |
| `derivative_matrix(g, dim, Val(order), Direct())` | `D₁` or `D₂` of a bounded direction, used for derivatives and Laplacians |
| `derivative_matrix(g, dim, Val(order), DiscreteAdjoint())` | `D₁⁺` or `D₂⁺` |
| `growto(g, sizes)` | the same grid with new Fourier resolutions |
| `convert(T, g)` | the grid with data converted to the scalar type `T`; FDGrids operators keep their banded structure, and other matrix types throw rather than becoming dense |

Concrete grids (channel, cavities, square duct, pipe) and flow cases are in
[ReSolverCases](https://github.com/Davide-Lasagna-s-Lab/ReSolver-Cases.jl).
