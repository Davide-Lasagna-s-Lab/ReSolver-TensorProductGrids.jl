@testset verbose=true "Two-dimensional lid-driven-cavity case                      " begin
    @testset verbose=true "Equation-constructor configuration                          " begin
        g = LidDrivenCavity2DGrid(13, 11; dist=FDGrids.GaussLobattoGrid(), width=5)
        X, Y, _ = points(g)
        U = @. 16X^2 * (1 - X)^2 * (3Y^2 - 2Y)
        V = @. -32X * (1 - X) * (1 - 2X) * Y^2 * (Y - 1)
        base_flow = (U, V)
        equations = LidDrivenCavity2DFlow(g, 100; base_flow, fftw_flags=FFTW.ESTIMATE)

        @test equations isa ProjectedNSE
        @test equations.base === base_flow
        @test equations.nl isa CartesianPrimitive2DNSE
        @test equations.ln isa CartesianPrimitive2DLNSE{AdjointDiscrete}
        @test equations.nl.force isa NoForce

        continuous = LidDrivenCavity2DFlow(g, 100; base_flow, mode=AdjointContinuous(),
                                           fftw_flags=FFTW.ESTIMATE)
        @test continuous.ln isa CartesianPrimitive2DLNSE{AdjointContinuous}
        @test_throws UndefKeywordError LidDrivenCavity2DFlow(g, 100; fftw_flags=FFTW.ESTIMATE)
        @test_throws ArgumentError LidDrivenCavity2DFlow(g, 100; base_flow=(U,),
                                                        fftw_flags=FFTW.ESTIMATE)
        @test_throws ArgumentError LidDrivenCavity2DFlow(g, 100; base_flow=[U, V],
                                                        fftw_flags=FFTW.ESTIMATE)
        @test_throws ArgumentError LidDrivenCavity2DFlow(g, 100; base_flow, mode=Forward(),
                                                        fftw_flags=FFTW.ESTIMATE)
    end
end
