using Test
using LinearAlgebra

import FDGrids

using NSEBase
using ReSolverRectangularGrids

include("helpers/analytical.jl")

@testset verbose=true "ReSolver rectangular grids                                  " begin
    include("grids/rectangular.jl")
    include("grids/channel.jl")
    include("grids/lid_driven_cavity_2d.jl")
    include("grids/lid_driven_cavity_3d.jl")
    include("grids/square_duct.jl")
    include("forcings.jl")
    include("cases/channel.jl")
    include("cases/lid_driven_cavity_2d.jl")
    include("cases/lid_driven_cavity_3d.jl")
    include("cases/square_duct.jl")
    include("examples.jl")
end
