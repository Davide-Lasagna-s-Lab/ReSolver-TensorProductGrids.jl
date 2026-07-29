using LinearAlgebra

import FDGrids

using NSEBase
using ReSolverRectangularGrids

let
    grid = ChannelGrid(19, 25, 19; Nt=1, α=0.5, β=1.25, width=7)
    y, x, z, t = points(grid)

    periodic_norm2 = sum(inv(float(factorial(k)))^2 for k in 0:20)
    u(y, x, z, t) = (1 - y^2) * exp(cos(0.5x)) * exp(cos(1.25z))
    uz(y, x, z, t) = -1.25sin(1.25z) * u(y, x, z, t)

    û = FFT(Field(grid, u))
    numerical = ddz!(FTField(grid), û)
    exact = FFT(Field(grid, uz))
    derivative_error = maximum(abs, parent(numerical) .- parent(exact)) / maximum(abs, parent(exact))
    norm_error = abs(norm(û)^2 - (16 / 15) * periodic_norm2^2)

    @assert derivative_error < 5e-7
    @assert norm_error < 5e-11

    Reτ = 180
    U = (Reτ / 2) .* plane_poiseuille_base(grid)
    equations = PlanePoiseuilleFlow(grid, Reτ; f=1, base_flow=(U, nothing, nothing),
                                    fftw_flags=FFTW.ESTIMATE, dealias=false)

    (; grid_size=size(grid), derivative_error, norm_error, equations_type=typeof(equations))
end
