using Documenter
using ReSolverRectangularGrids

DocMeta.setdocmeta!(
    ReSolverRectangularGrids,
    :DocTestSetup,
    :(using LinearAlgebra, NSEBase, ReSolverRectangularGrids);
    recursive=true,
)

makedocs(
    sitename="ReSolverRectangularGrids.jl",
    modules=[ReSolverRectangularGrids],
    authors="Davide Lasagna, Thomas Burton, and contributors",
    doctest=true,
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Davide-Lasagna-s-Lab.github.io/ReSolver-RectangularGrids.jl/dev/",
        edit_link="main",
        repolink="https://github.com/Davide-Lasagna-s-Lab/ReSolver-RectangularGrids.jl",
    ),
    remotes=Dict(normpath(joinpath(@__DIR__, "..")) =>
                 (Documenter.Remotes.GitHub("Davide-Lasagna-s-Lab", "ReSolver-RectangularGrids.jl"), "main")),
    pages=[
        "Home" => "index.md",
        "Getting started" => [
            "Installation" => "installation.md",
            "Quick start" => "quickstart.md",
        ],
        "Concepts" => [
            "Coordinates and conventions" => "conventions.md",
            "Rectangular-grid design" => "rectangular.md",
        ],
        "Grids" => [
            "Channel" => "grids/channel.md",
            "2D lid-driven cavity" => "grids/lid_driven_cavity_2d.md",
            "3D lid-driven cavity" => "grids/lid_driven_cavity_3d.md",
            "Square duct" => "grids/square_duct.md",
        ],
        "Cases" => [
            "Reusable forcings" => "forcings.md",
            "Channel flows" => "cases/channel.md",
            "2D lid-driven cavity" => "cases/lid_driven_cavity_2d.md",
            "3D lid-driven cavity" => "cases/lid_driven_cavity_3d.md",
            "Square-duct flow" => "cases/square_duct.md",
        ],
        "Worked examples" => "examples.md",
        "Extending the package" => "extending.md",
        "API reference" => "api.md",
    ],
    checkdocs=:exports,
    linkcheck=get(ENV, "DOCUMENTER_LINKCHECK", "false") == "true",
    warnonly=false,
)

deploydocs(repo="github.com/Davide-Lasagna-s-Lab/ReSolver-RectangularGrids.jl.git", devbranch="main")
