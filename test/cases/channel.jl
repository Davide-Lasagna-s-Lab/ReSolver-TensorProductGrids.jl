@testset verbose=true "Channel case                                                " begin
    @testset verbose=true "Base flow and equation-constructor configuration            " begin
        g = ChannelGrid(7, 17, 7; Nt=3, α=1.25, β=0.75, width=5)
        y = only(g.xs)
        couette_base = plane_couette_base(g)
        poiseuille_base = plane_poiseuille_base(g)

        @test couette_base == y
        @test couette_base !== y
        @test poiseuille_base ≈ 1 .- y .^ 2
        @test poiseuille_base[[1, end]] ≈ [0, 0] atol=1e-14

        base_flow = (couette_base, nothing, nothing)
        couette = PlaneCouetteFlow(g, 400; base_flow, Ro=0.25,
                                   fftw_flags=FFTW.ESTIMATE, dealias=false)
        default_couette = PlaneCouetteFlow(g, 400; fftw_flags=FFTW.ESTIMATE,
                                           dealias=false)
        @test couette isa ProjectedNSE
        @test couette.base === base_flow
        @test default_couette.base == (plane_couette_base(g), nothing, nothing)
        @test couette.nl isa CartesianPrimitive3DNSE
        @test couette.ln isa CartesianPrimitive3DLNSE{AdjointDiscrete}
        @test couette.nl.force isa CoriolisForce
        @test couette.ln.force isa CoriolisForce
        @test couette.nl.force.Ro == 0.25
        @test default_couette.nl.force isa NoForce

        poiseuille = PlanePoiseuilleFlow(g, 400; f=1.5, fftw_flags=FFTW.ESTIMATE,
                                         dealias=false)
        @test poiseuille.base == (poiseuille_base, nothing, nothing)
        @test poiseuille.nl.force isa ConstantBodyForce
        @test poiseuille.ln.force isa ConstantBodyForce
        @test poiseuille.nl.force.value == 1.5
        @test poiseuille.nl.force.i == 1
        @test poiseuille.ln.force.value == 1.5
        @test poiseuille.ln.force.i == 1
        @test_throws ArgumentError PlanePoiseuilleFlow(
            g, 400; f=0, fftw_flags=FFTW.ESTIMATE, dealias=false)
        @test_throws MethodError PlanePoiseuilleFlow(
            g, 400; Ro=0.25, fftw_flags=FFTW.ESTIMATE, dealias=false)

        continuous = PlaneCouetteFlow(g, 400; base_flow, mode=AdjointContinuous(),
                                      fftw_flags=FFTW.ESTIMATE, dealias=false)
        @test continuous.ln isa CartesianPrimitive3DLNSE{AdjointContinuous}

        @test_throws ArgumentError PlaneCouetteFlow(
            g, 400; base_flow=(couette_base, nothing), fftw_flags=FFTW.ESTIMATE)
        @test_throws ArgumentError PlanePoiseuilleFlow(
            g, 400; base_flow=[poiseuille_base, nothing, nothing],
            fftw_flags=FFTW.ESTIMATE)
        @test_throws DimensionMismatch PlaneCouetteFlow(
            g, 400; base_flow=(zeros(16), nothing, nothing), fftw_flags=FFTW.ESTIMATE)
        @test_throws ArgumentError PlaneCouetteFlow(
            g, 400; base_flow=(1, nothing, nothing), fftw_flags=FFTW.ESTIMATE)
        @test_throws ArgumentError PlaneCouetteFlow(
            g, 400; base_flow, mode=Forward(), fftw_flags=FFTW.ESTIMATE)
    end
end
