@testset verbose=true "Square-duct case                                            " begin
    @testset verbose=true "Base flow and equation-constructor configuration            " begin
        g = SquareDuctGrid(13, 5; Nt=3, α=1.25, width=5)
        x, y = g.xs
        W = (x .* (1 .- x)) * transpose(y .* (1 .- y))
        base_flow = (nothing, nothing, W)

        equations = SquareDuctFlow(g, 500; base_flow, f=1.25,
                                   fftw_flags=FFTW.ESTIMATE, dealias=false)
        default_equations = SquareDuctFlow(g, 500; fftw_flags=FFTW.ESTIMATE,
                                           dealias=false)
        @test equations isa ProjectedNSE
        @test equations.base === base_flow
        @test default_equations.base == (nothing, nothing, nothing)
        @test equations.nl isa CartesianPrimitive3DNSE
        @test equations.ln isa CartesianPrimitive3DLNSE{AdjointDiscrete}
        @test equations.nl.force isa ConstantBodyForce
        @test equations.ln.force isa ConstantBodyForce
        @test equations.nl.force.value == 1.25
        @test equations.nl.force.i == 3

        @test_throws ArgumentError SquareDuctFlow(
            g, 500; f=0, fftw_flags=FFTW.ESTIMATE, dealias=false)

        continuous = SquareDuctFlow(g, 500; base_flow, mode=AdjointContinuous(),
                                    fftw_flags=FFTW.ESTIMATE, dealias=false)
        @test continuous.ln isa CartesianPrimitive3DLNSE{AdjointContinuous}

        @test_throws ArgumentError SquareDuctFlow(
            g, 500; base_flow=(nothing, W), fftw_flags=FFTW.ESTIMATE)
        @test_throws ArgumentError SquareDuctFlow(
            g, 500; base_flow=[nothing, nothing, W], fftw_flags=FFTW.ESTIMATE)
        @test_throws DimensionMismatch SquareDuctFlow(
            g, 500; base_flow=(nothing, nothing, zeros(12, 13)),
            fftw_flags=FFTW.ESTIMATE)
        @test_throws ArgumentError SquareDuctFlow(
            g, 500; base_flow=(nothing, nothing, 1), fftw_flags=FFTW.ESTIMATE)
        @test_throws ArgumentError SquareDuctFlow(
            g, 500; base_flow, mode=Forward(), fftw_flags=FFTW.ESTIMATE)
    end
end
