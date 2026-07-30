# Extending the package

Bundled constructors hide axes, Fourier order, adjoint matrices, and quadrature details behind a physical API. When a genuinely new rectangular layout is needed, build it through the public FDGrids and `RectangularGrid` interfaces shown below.

## A complete custom layout

The following example defines a two-dimensional physical domain with bounded ``x\in[0,1]``,
periodic ``y``, no ``z`` coordinate, and a unit-period time or phase coordinate. Arrays are stored
as `(x,y,t)`:

```@example custom-layout
using LinearAlgebra
import FDGrids
using NSEBase
using ReSolverRectangularGrids

Nx, Ny, Nt = 17, 9, 5
β = 2
fd = FDGrids.grid(Nx, 0, 1, FDGrids.GaussLobattoGrid())
x, wx = fd.xs, fd.ws

Dx = FDGrids.DiffMatrix(x, 5, 1)
Dx2 = FDGrids.DiffMatrix(x, 5, 2)
Dxa = adjoint(Dx, wx)
Dx2a = adjoint(Dx2, wx)

axes = (1, 2, nothing, 3)
fft_order = (2, 3)
g = RectangularGrid((x,), (Dx,), (Dx2,), (Dxa,), (Dx2a,), (wx,), (β, 2π),
                    (Nx, Ny, Nt), axes, fft_order)

@assert size(g) == (Nx, Ny, Nt)
@assert NSEBase.fft_physical_dims(g) == (:y, :t)
@assert NSEBase.inhomogeneous_physical_dims(g) == (:x,)
g
```

The assignment above intentionally spells out the low-level data:

- each bounded direction contributes points, first and second derivatives, their weighted adjoints, and quadrature weights;
- bounded tuples follow increasing storage dimension;
- Fourier scales follow `fft_order`;
- `axes` always has the four physical slots `(x,y,z,t)`, even when one is absent;
- `grid_size` is in storage order and all Fourier resolutions are positive and odd.

## Add a user-facing constructor

End users should not need the low-level arguments. Put layout constants, an abstract dispatch alias, a concrete alias, and a constructor in one file under `src/grids/`. Follow the bundled pattern:

```julia
const MY_AXES = (1, 2, nothing, 3)
const MY_FFT_ORDER = (2, 3)
const AbstractMyGrid{T} = AbstractGrid{T, 3, MY_AXES, MY_FFT_ORDER}
```

The constructor should accept physical resolutions and meaningful physical keywords, choose documented defaults, create its FDGrids directions, and call `RectangularGrid`. Keep axes and storage details out of the user API. Export an abstract alias when downstream wrappers or equation constructors should dispatch on any compatible layout.

## Add a physical case

Place equation factories separately under `src/cases/`. A case file should:

1. dispatch on the abstract layout contract;
2. state the governing equations and nondimensional parameters in its docstring;
3. validate the number and meaning of `base_flow` components;
4. compose reusable forcing policies;
5. forward `mode`, `fftw_flags`, and `dealias` to `NSEBase.construct_equations`;
6. explain which boundary data remain external.

Keep forcing objects under `src/forcings.jl` when more than one case can reuse them.

## Validation checklist

A new layout should mirror the existing per-case tests:

- constructor defaults, custom intervals, distributions, and widths;
- `AXES`, Fourier order, physical/storage accessors, and `growto`;
- tensor-product weights and analytical scalar/vector norms;
- first derivatives, Laplacian, and homogeneous shifts using a function with all Fourier modes;
- exact quadrature-weighted identities for every discrete-adjoint derivative;
- equation-factory configuration, forcing component, base lifting, and both adjoint modes;
- one standalone, CI-executed example with the same structure as the bundled cases.

See [`test/grids/`](https://github.com/Davide-Lasagna-s-Lab/ReSolver-RectangularGrids.jl/tree/main/test/grids) and [`test/cases/`](https://github.com/Davide-Lasagna-s-Lab/ReSolver-RectangularGrids.jl/tree/main/test/cases) for the maintained templates.
