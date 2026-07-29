using LinearAlgebra

import FDGrids

using NSEBase
using ReSolverRectangularGrids

let
    grid = LidDrivenCavity2DGrid(17; Nt=19, dist=FDGrids.GaussLobattoGrid(), width=5)
    x, y, t = points(grid)

    periodic_norm2 = sum(inv(float(factorial(k)))^2 for k in 0:20)
    u(x, y, t) = x * (1 - x) * y * (1 - y) * exp(cos(2π * t))
    ut(x, y, t) = -2π * sin(2π * t) * u(x, y, t)

    û = FFT(Field(grid, u))
    numerical = ddt!(FTField(grid), û)
    exact = FFT(Field(grid, ut))
    derivative_error = maximum(abs, parent(numerical) .- parent(exact)) / maximum(abs, parent(exact))
    norm_error = abs(norm(û)^2 - periodic_norm2 / 30^2)

    @assert derivative_error < 5e-7
    @assert norm_error < 5e-11

    base_flow = lid_driven_cavity_2d_base(grid)
    equations = LidDrivenCavity2DFlow(grid, 1000; fftw_flags=FFTW.ESTIMATE, dealias=false)
    @assert equations.base == base_flow

    (; grid_size=size(grid), derivative_error, norm_error, equations_type=typeof(equations))
end
