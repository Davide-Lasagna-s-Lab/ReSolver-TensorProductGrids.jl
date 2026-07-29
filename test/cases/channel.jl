@testset verbose=true "Channel cases                                               " begin
    @testset verbose=true "Equation-constructor configuration                          " begin
        g = ChannelGrid(7, 17, 7; Nt=3, α=1.25, β=0.75, width=5)
        y = only(g.xs)
        couette_base_flow = plane_couette_base(g)
        poiseuille_base_flow = plane_poiseuille_base(g)

        @test couette_base_flow == y
        @test couette_base_flow !== y
        @test poiseuille_base_flow ≈ 1 .- y .^ 2
        @test poiseuille_base_flow[[1, end]] ≈ [0, 0] atol=1e-14

        base_flow = (couette_base_flow, nothing, nothing)
        couette = PlaneCouetteFlow(g, 400; base_flow, Ro=0.25,
                                   fftw_flags=FFTW.ESTIMATE, dealias=false)
        unrotated = PlaneCouetteFlow(g, 400; fftw_flags=FFTW.ESTIMATE)
        @test couette isa ProjectedNSE
        @test couette.base === base_flow
        @test couette.nl.force isa CoriolisForce
        @test couette.ln.force isa CoriolisForce
        @test couette.nl.force.Ro == 0.25
        @test fieldnames(typeof(couette.nl.force)) == (:Ro,)
        @test unrotated.nl.force isa NoForce
        @test couette.ln isa CartesianPrimitive3DLNSE{AdjointDiscrete}

        poiseuille = PlanePoiseuilleFlow(g, 400; Ro=0.25, f=1.5,
                                         fftw_flags=FFTW.ESTIMATE, dealias=false)
        @test poiseuille.base == (poiseuille_base_flow, nothing, nothing)
        @test poiseuille.nl.force isa CompoundForcing
        pressure, rotation = poiseuille.nl.force.forces
        @test pressure isa ConstantBodyForce
        @test pressure.value == 1.5
        @test pressure.i == 1
        @test rotation isa CoriolisForce
        @test rotation.Ro == 0.25

        continuous = PlaneCouetteFlow(g, 400; mode=AdjointContinuous(),
                                      fftw_flags=FFTW.ESTIMATE)
        @test continuous.ln isa CartesianPrimitive3DLNSE{AdjointContinuous}
        @test_throws ArgumentError PlaneCouetteFlow(g, 400; base_flow=(nothing, nothing),
                                                   fftw_flags=FFTW.ESTIMATE)
        @test_throws ArgumentError PlanePoiseuilleFlow(g, 400;
                                                      base_flow=[nothing, nothing, nothing],
                                                      fftw_flags=FFTW.ESTIMATE)
        @test_throws ArgumentError PlaneCouetteFlow(g, 400; mode=Forward(),
                                                   fftw_flags=FFTW.ESTIMATE)
    end
end
