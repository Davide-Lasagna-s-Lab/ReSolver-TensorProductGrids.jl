@testset verbose=true "Three-dimensional lid-driven-cavity case                    " begin
    @testset verbose=true "Equation-constructor configuration                          " begin
        g = LidDrivenCavity3DGrid(11, 11, 11; dist=FDGrids.GaussLobattoGrid(), width=5)
        X, Y, Z, _ = points(g)
        G = @. 16Z^2 * (1 - Z)^2
        U = @. 16X^2 * (1 - X)^2 * (3Y^2 - 2Y) * G
        V = @. -32X * (1 - X) * (1 - 2X) * Y^2 * (Y - 1) * G
        base_flow = (U, V, nothing)
        equations = LidDrivenCavity3DFlow(g, 100; base_flow, fftw_flags=FFTW.ESTIMATE)

        @test equations isa ProjectedNSE
        @test equations.base === base_flow
        @test equations.nl isa CartesianPrimitive3DNSE
        @test equations.ln isa CartesianPrimitive3DLNSE{AdjointDiscrete}
        @test equations.nl.force isa NoForce

        continuous = LidDrivenCavity3DFlow(g, 100; base_flow, mode=AdjointContinuous(),
                                           fftw_flags=FFTW.ESTIMATE)
        @test continuous.ln isa CartesianPrimitive3DLNSE{AdjointContinuous}
        @test_throws UndefKeywordError LidDrivenCavity3DFlow(g, 100; fftw_flags=FFTW.ESTIMATE)
        @test_throws ArgumentError LidDrivenCavity3DFlow(g, 100; base_flow=(U, V),
                                                        fftw_flags=FFTW.ESTIMATE)
        @test_throws ArgumentError LidDrivenCavity3DFlow(g, 100; base_flow=[U, V, nothing],
                                                        fftw_flags=FFTW.ESTIMATE)
        @test_throws ArgumentError LidDrivenCavity3DFlow(g, 100; base_flow, mode=Forward(),
                                                        fftw_flags=FFTW.ESTIMATE)
    end
end
