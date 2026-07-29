using LinearAlgebra

import FDGrids

using NSEBase
using ReSolverRectangularGrids

let
    grid = LidDrivenCavity3DGrid(15, 15, 15; Nt=19, dist=FDGrids.GaussLobattoGrid(), width=5)
    x, y, z, t = points(grid)

    periodic_norm2 = sum(inv(float(factorial(k)))^2 for k in 0:20)
    u(x, y, z, t) = x * (1 - x) * y * (1 - y) * z * (1 - z) * exp(cos(2π * t))
    ut(x, y, z, t) = -2π * sin(2π * t) * u(x, y, z, t)

    û = FFT(Field(grid, u))
    numerical = ddt!(FTField(grid), û)
    exact = FFT(Field(grid, ut))
    derivative_error = maximum(abs, parent(numerical) .- parent(exact)) / maximum(abs, parent(exact))
    norm_error = abs(norm(û)^2 - periodic_norm2 / 30^3)

    @assert derivative_error < 5e-7
    @assert norm_error < 5e-11

    X, Y, Z, _ = points(grid)
    G = @. 16Z^2 * (1 - Z)^2
    U = @. 16X^2 * (1 - X)^2 * (3Y^2 - 2Y) * G
    V = @. -32X * (1 - X) * (1 - 2X) * Y^2 * (Y - 1) * G
    equations = LidDrivenCavity3DFlow(grid, 1000; base_flow=(U, V, nothing),
                                      fftw_flags=FFTW.ESTIMATE, dealias=false)

    (; grid_size=size(grid), derivative_error, norm_error, equations_type=typeof(equations))
end
