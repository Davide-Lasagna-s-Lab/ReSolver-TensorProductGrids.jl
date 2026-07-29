using LinearAlgebra

import FDGrids

using NSEBase
using ReSolverRectangularGrids

let
    grid = SquareDuctGrid(17, 19; Nt=1, α=0.75, width=5)
    x, y, z, t = points(grid)

    periodic_norm2 = sum(inv(float(factorial(k)))^2 for k in 0:20)
    u(x, y, z, t) = x * (1 - x) * y * (1 - y) * exp(cos(0.75z))
    uz(x, y, z, t) = -0.75sin(0.75z) * u(x, y, z, t)

    û = FFT(Field(grid, u))
    numerical = ddz!(FTField(grid), û)
    exact = FFT(Field(grid, uz))
    derivative_error = maximum(abs, parent(numerical) .- parent(exact)) / maximum(abs, parent(exact))
    norm_error = abs(norm(û)^2 - periodic_norm2 / 30^2)

    @assert derivative_error < 5e-7
    @assert norm_error < 5e-11

    Reτ = 180
    equations = SquareDuctFlow(grid, 2Reτ; f=4, fftw_flags=FFTW.ESTIMATE, dealias=false)

    (; grid_size=size(grid), derivative_error, norm_error, equations_type=typeof(equations))
end
