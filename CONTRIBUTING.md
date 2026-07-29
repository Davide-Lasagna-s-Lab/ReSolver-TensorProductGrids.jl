# Contributing

Contributions to ReSolver-RectangularGrids.jl are welcome. Bug reports, documentation improvements, new tests, and focused extensions to the rectangular-grid or case APIs are all useful.

The package is currently development-only. It depends on functionality from the development branches of NSEBase.jl and FDGrids.jl that has not yet been released, so compatibility may change until those upstream releases are available.

## Before opening a pull request

For a bug or a proposed API change, please open an issue first when discussion would help establish the expected behavior. Keep pull requests focused on one concern and include tests for behavioral changes.

Use Julia 1.10 or later. After cloning the repository, install the required development dependencies and instantiate the project:

```sh
julia --project=. -e 'using Pkg; Pkg.develop(url="https://github.com/Davide-Lasagna-s-Lab/NSEBase.jl", rev="dev"); Pkg.develop(url="https://github.com/Davide-Lasagna-s-Lab/FDGrids.jl", rev="master"); Pkg.instantiate()'
```

Run the package tests from the repository root:

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

Build the documentation when changing public APIs, docstrings, examples, or documentation pages:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.develop(url="https://github.com/Davide-Lasagna-s-Lab/NSEBase.jl", rev="dev"); Pkg.develop(url="https://github.com/Davide-Lasagna-s-Lab/FDGrids.jl", rev="master"); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The product-weight benchmark has its own environment:

```sh
julia --project=benchmark -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=benchmark scratch/benchmark_rectangular_product_weights.jl
```

## Development guidelines

- Preserve the physical coordinate names and storage conventions documented for each grid.
- Prefer small, explicit implementations over generalized machinery when only the supported dimensions are required.
- Extend the generic NSEBase interfaces rather than introducing parallel interfaces.
- Document mathematical conventions, defaults, units, and boundary-condition responsibilities.
- Add analytical derivative and norm checks when introducing or changing a grid.
- Keep examples executable and representative of the public API.
- Avoid unrelated formatting or refactoring in a focused pull request.

## Pull requests

In the pull-request description, explain the motivation and the user-visible result, record the checks you ran, and call out compatibility or migration consequences. A maintainer may ask for revisions before merging; this is part of keeping the numerical and public APIs dependable.
