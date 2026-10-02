using Test
using FDGrids
using ReSolverFlowsBase, ReSolverTensorProductGrids
using ReSolverFlowsBase: derivative_matrix

# a channel-like layout: bounded y stored first, Fourier x, z, t
function channel_like(; Ny=17, Nx=9, Nz=7, Nt=3, α=1.5, β=2.5)
    fd = FDGrids.grid(Ny, -1, 1, GaussLobattoGrid())
    y  = Vector(fd.xs)
    w  = Vector(fd.ws)
    D₁ = DiffMatrix(y, 5, 1)
    D₂ = DiffMatrix(y, 5, 2)

    return TensorProductGrid((y,), (D₁,), (D₂,), (adjoint(D₁, w),), (adjoint(D₂, w),), (w,),
                             (α, β, 1), (Ny, Nx, Nz, Nt), (2, 1, 3, 4), (2, 3, 4))
end

@testset "TensorProductGrid" begin
    g = channel_like()

    # ---- layout ----
    @test size(g) == (17, 9, 7, 3)
    @test eltype(g) == Float64
    @test fft_storage_dims(g) == (2, 3, 4)
    @test inhomogeneous_storage_dims(g) == (1,)
    @test wavenumber_scale.(Ref(g), (1, 2, 3, 4)) == (1, 1.5, 2.5, 1)

    # ---- coordinates ----
    y, x, z, t = points(g)
    @test size(y) == (17, 1, 1, 1) && size(x) == (1, 9, 1, 1)
    @test vec(y) == g.xs[1]
    @test x[2] ≈ 2π / 1.5 / 9
    @test map(length, points(g; dealias=true)) == (17, 15, 11, 5)

    # ---- quadrature ----
    @test sum(weights(g)) ≈ 2

    # ---- derivative operators ----
    @test derivative_matrix(g, 1, Val(1), Direct())          === g.D₁[1]
    @test derivative_matrix(g, 1, Val(2), Direct())          === g.D₂[1]
    @test derivative_matrix(g, 1, Val(1), DiscreteAdjoint()) === g.D₁⁺[1]
    @test derivative_matrix(g, 1, Val(2), DiscreteAdjoint()) === g.D₂⁺[1]

    # ---- resolution growth shares all bounded data ----
    h = growto(g, (15, 11, 5))
    @test size(h) == (17, 15, 11, 5)
    @test h.D₁[1] === g.D₁[1] && h.ws[1] === g.ws[1]

    # ---- scalar type conversion ----
    g32 = convert(Float32, g)
    @test convert(Float64, g) === g
    @test eltype(g32) == Float32
    @test eltype(g32.xs[1]) == eltype(g32.ws[1]) == Float32
    @test g32.wavenumber_scales == Float32.((1.5, 2.5, 1))
    @test g32.D₁[1]  isa FDGrids.DiffMatrix{Float32}
    @test g32.D₁⁺[1] isa FDGrids.AdjointDiffMatrix{Float32}

    # ---- no silent dense fallback for operators without a conversion ----
    d = TensorProductGrid((g.xs[1],), (Matrix(g.D₁[1]),), g.D₂, g.D₁⁺, g.D₂⁺, g.ws,
                          (1.5, 2.5, 1), size(g), (2, 1, 3, 4), (2, 3, 4))
    @test_throws ArgumentError convert(Float32, d)
end

@testset "TensorProductWeights" begin
    w = TensorProductWeights(([1.0, 2.0], [3.0, 4.0, 5.0]))
    @test size(w) == (2, 3)
    @test w == [1, 2] * [3, 4, 5]'
end
