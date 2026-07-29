@testset verbose=true "Reusable case forcings                                      " begin
    @testset verbose=true "Coriolis forward and adjoint actions                        " begin
        g = ChannelGrid(3, 9, 3; width=3)
        u, out = VectorField(g; N=3), VectorField(g; N=3)
        parent(u[1]) .= 1
        parent(u[2]) .= 2
        parent(u[3]) .= 3

        force = CoriolisForce(0.25)
        force(out, u, Forward())
        @test parent(out[1]) == fill(0.5, size(out[1]))
        @test parent(out[2]) == fill(-0.25, size(out[2]))
        @test iszero(parent(out[3]))

        for mode in (AdjointContinuous(), AdjointDiscrete())
            out .= 0
            force(out, u, mode)
            @test parent(out[1]) == fill(-0.5, size(out[1]))
            @test parent(out[2]) == fill(0.25, size(out[2]))
            @test iszero(parent(out[3]))
        end

        @test force.Ro == 0.25
        @test fieldnames(typeof(force)) == (:Ro,)
        @test sizeof(force) == sizeof(force.Ro)

        cavity = LidDrivenCavity2DGrid(7; width=3)
        u2, out2 = VectorField(cavity; N=2), VectorField(cavity; N=2)
        parent(u2[1]) .= 1
        parent(u2[2]) .= 2
        force(out2, u2, Forward())
        @test parent(out2[1]) == fill(0.5, size(out2[1]))
        @test parent(out2[2]) == fill(-0.25, size(out2[2]))
    end

    @testset verbose=true "Constant forcing of the homogeneous mean mode               " begin
        g = ChannelGrid(3, 9, 3; Nt=3, width=3)
        u, out = VectorField(g; N=3), VectorField(g; N=3)
        force = ConstantBodyForce(2.5; i=1)

        force(out, u, Forward())
        force(out, u, AdjointDiscrete())
        @test parent(out[1])[:, 1, 1, 1] == fill(2.5, size(g, 1))
        parent(out[1])[:, 1, 1, 1] .= 0
        @test iszero(parent(out[1]))
        @test iszero(parent(out[2]))
        @test iszero(parent(out[3]))
        @test_throws ArgumentError ConstantBodyForce(0; i=1)
        @test_throws ArgumentError ConstantBodyForce(1; i=4)(out, u, Forward())
    end
end
