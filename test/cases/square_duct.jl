@testset verbose=true "Square-duct case                                            " begin
    @testset verbose=true "Equation-constructor configuration                          " begin
        g = SquareDuctGrid(13, 5; Nt=3, α=1.25, width=5)
        x, y = g.xs
        W = reshape(x .* (1 .- x), :, 1) .* reshape(y .* (1 .- y), 1, :)
        base_flow = (nothing, nothing, W)
        equations = SquareDuctFlow(g, 500; base_flow, f=1.25,
                                   fftw_flags=FFTW.ESTIMATE, dealias=false)

        @test equations isa ProjectedNSE
        @test equations.base === base_flow
        @test equations.nl isa CartesianPrimitive3DNSE
        @test equations.ln isa CartesianPrimitive3DLNSE{AdjointDiscrete}
        @test equations.nl.force isa ConstantBodyForce
        @test equations.ln.force isa ConstantBodyForce
        @test equations.nl.force.value == 1.25
        @test equations.nl.force.i == 3

        continuous = SquareDuctFlow(g, 500; mode=AdjointContinuous(),
                                    fftw_flags=FFTW.ESTIMATE)
        @test continuous.ln isa CartesianPrimitive3DLNSE{AdjointContinuous}
        @test_throws ArgumentError SquareDuctFlow(g, 500; base_flow=(nothing, nothing),
                                                 fftw_flags=FFTW.ESTIMATE)
        @test_throws ArgumentError SquareDuctFlow(g, 500;
                                                 base_flow=[nothing, nothing, nothing],
                                                 fftw_flags=FFTW.ESTIMATE)
        @test_throws ArgumentError SquareDuctFlow(g, 500; mode=Forward(),
                                                 fftw_flags=FFTW.ESTIMATE)
    end
end
