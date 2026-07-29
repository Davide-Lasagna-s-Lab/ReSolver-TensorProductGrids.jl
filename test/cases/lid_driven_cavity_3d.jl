@testset verbose=true "Three-dimensional lid-driven-cavity case                    " begin
    @testset verbose=true "Base flow and equation-constructor configuration            " begin
        g = LidDrivenCavity3DGrid(11; lim=(-1, 1), dist=FDGrids.GaussLobattoGrid(),
                                  width=5)
        U, V, W = base_flow = lid_driven_cavity_3d_base(g)
        x = first(g.xs)
        ξ = (x .- first(x)) ./ (last(x) - first(x))
        lid_factor = @. 16ξ^2 * (1 - ξ)^2
        lid = lid_factor * transpose(lid_factor)

        @test size.(base_flow[1:2]) == ((11, 11, 11), (11, 11, 11))
        @test isnothing(W)
        @test U[:, end, :] ≈ lid atol=5e-14
        @test maximum(U[:, end, :]) ≈ 1 atol=5e-14
        @test iszero(U[:, 1, :])
        @test iszero(U[[1, end], :, :])
        @test iszero(U[:, :, [1, end]])
        @test iszero(V[:, [1, end], :])
        @test iszero(V[[1, end], :, :])
        @test iszero(V[:, :, [1, end]])

        Û = ftfield_from_inhomogeneous(g, U)
        V̂ = ftfield_from_inhomogeneous(g, V)
        divergence = parent(ddx!(FTField(g), Û)) .+ parent(ddy!(FTField(g), V̂))
        @test maximum(abs, divergence) < 2e-12

        equations = LidDrivenCavity3DFlow(g, 100; base_flow,
                                           fftw_flags=FFTW.ESTIMATE, dealias=false)
        default_equations = LidDrivenCavity3DFlow(g, 100; fftw_flags=FFTW.ESTIMATE,
                                                  dealias=false)
        @test equations isa ProjectedNSE
        @test equations.base === base_flow
        @test default_equations.base == lid_driven_cavity_3d_base(g)
        @test equations.nl isa CartesianPrimitive3DNSE
        @test equations.ln isa CartesianPrimitive3DLNSE{AdjointDiscrete}
        @test equations.nl.force isa NoForce
        @test equations.ln.force isa NoForce

        continuous = LidDrivenCavity3DFlow(g, 100; base_flow,
                                            mode=AdjointContinuous(),
                                            fftw_flags=FFTW.ESTIMATE, dealias=false)
        @test continuous.ln isa CartesianPrimitive3DLNSE{AdjointContinuous}

        @test_throws ArgumentError LidDrivenCavity3DFlow(g, 100; base_flow=(U, V),
                                                         fftw_flags=FFTW.ESTIMATE)
        @test_throws ArgumentError LidDrivenCavity3DFlow(
            g, 100; base_flow=[U, V, nothing], fftw_flags=FFTW.ESTIMATE)
        @test_throws DimensionMismatch LidDrivenCavity3DFlow(
            g, 100; base_flow=(zeros(10, 11, 11), V, nothing),
            fftw_flags=FFTW.ESTIMATE)
        @test_throws ArgumentError LidDrivenCavity3DFlow(
            g, 100; base_flow=(1, V, nothing), fftw_flags=FFTW.ESTIMATE)
        @test_throws ArgumentError LidDrivenCavity3DFlow(
            g, 100; base_flow, mode=Forward(), fftw_flags=FFTW.ESTIMATE)
    end
end
