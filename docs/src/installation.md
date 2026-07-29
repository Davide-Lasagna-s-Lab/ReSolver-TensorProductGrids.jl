# Installation

## Development installation

ReSolverRectangularGrids currently depends on grid and derivative interfaces from unreleased NSEBase and FDGrids branches. The package repository can be used now, but ordinary registry resolution would select older, incompatible upstream versions. Install the tested branches explicitly:

```julia
using Pkg

Pkg.add([
    PackageSpec(url="https://github.com/Davide-Lasagna-s-Lab/NSEBase.jl", rev="dev"),
    PackageSpec(url="https://github.com/Davide-Lasagna-s-Lab/FDGrids.jl", rev="master"),
])
Pkg.develop(url="https://github.com/Davide-Lasagna-s-Lab/ReSolver-RectangularGrids.jl")
```

The package supports Julia 1.10 and later. `Pkg.add` records the selected upstream branches and their exact revisions in the manifest, while `Pkg.develop` tracks the ReSolverRectangularGrids checkout for local development.

!!! warning "No tagged package release yet"
    Do not infer that `Pkg.add("ReSolverRectangularGrids")` is supported. A tagged `v0.1.0` and registry entry will be created only after compatible NSEBase and FDGrids versions are released and the clean registry-resolved test job passes.

## Local development checkout

Clone the repository and develop the same upstream branches in its project:

```bash
git clone https://github.com/Davide-Lasagna-s-Lab/ReSolver-RectangularGrids.jl
cd ReSolver-RectangularGrids.jl
julia --project=. -e 'using Pkg; Pkg.add([PackageSpec(url="https://github.com/Davide-Lasagna-s-Lab/NSEBase.jl", rev="dev"), PackageSpec(url="https://github.com/Davide-Lasagna-s-Lab/FDGrids.jl", rev="master")]); Pkg.instantiate()'
```

Run the package tests:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Build the manual:

```bash
julia --project=docs -e 'using Pkg; Pkg.add([PackageSpec(url="https://github.com/Davide-Lasagna-s-Lab/NSEBase.jl", rev="dev"), PackageSpec(url="https://github.com/Davide-Lasagna-s-Lab/FDGrids.jl", rev="master")]); Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Local `Manifest.toml` files contain checkout paths and are intentionally ignored.

## Updating

Because these are development dependencies, inspect changes before updating a scientific environment. `Pkg.status(; mode=Pkg.PKGMODE_MANIFEST)` records the current revisions; commit the consuming project's manifest when exact reproducibility matters.

Once upstream releases are registered, this page will replace branch-based development with normal versioned `Pkg.add` instructions and document the minimum compatible versions.
