@testset verbose=true "Rectangular grid implementation                             " begin
    @testset verbose=true "Interface accessors and coordinate generation               " begin
        g = ChannelGrid(7, 13, 5; Nt=3, α=1.25, β=0.75, width=3)

        @test g isa RectangularGrid{1, Float64}
        @test eltype(g) === Float64
        @test size(g) == (13, 7, 5, 3)
        @test weights(g) === g.w
        @test weights(g) isa RectangularProductWeights{1}
        @test collect(weights(g)) == only(g.ws)
        @test wavenumber_scale(g, 1) == 1
        @test wavenumber_scale(g, 2) == 1.25
        @test wavenumber_scale(g, 3) == 0.75
        @test wavenumber_scale(g, 4) == 2π
        @test_throws BoundsError wavenumber_scale(g, 0)
        @test_throws BoundsError wavenumber_scale(g, 5)

        y, x, z, t = points(g)
        @test map(size, (y, x, z, t)) == ((13, 1, 1, 1), (1, 7, 1, 1),
                                                  (1, 1, 5, 1), (1, 1, 1, 3))
        @test vec(y) == g.xs[1]
        @test vec(x) ≈ (0:6) .* (2π / (1.25 * 7))
        @test vec(z) ≈ (0:4) .* (2π / (0.75 * 5))
        @test vec(t) ≈ (0:2) ./ 3

        padded = points(g; dealias=true)
        expected_sizes = (13, cld(3 * 7, 2) | 1, cld(3 * 5, 2) | 1, cld(3 * 3, 2) | 1)
        @test map(x -> only(filter(!=(1), size(x))), padded[2:4]) == expected_sizes[2:4]
        @test_throws ArgumentError points(g, (7, 0, 3))
    end

    @testset verbose=true "Growth, conversion, and compact adjoint storage             " begin
        g = ChannelGrid(7, 17, 5; Nt=3, α=1.25, β=0.75, width=5)
        grown = growto(g, (9, 7, 5))

        @test size(grown) == (17, 9, 7, 5)
        @test points(grown) == points(g, (9, 7, 5))
        @test grown.xs === g.xs
        @test grown.D₁ === g.D₁
        @test grown.D₂ === g.D₂
        @test grown.D₁⁺ === g.D₁⁺
        @test grown.D₂⁺ === g.D₂⁺
        @test grown.ws === g.ws
        @test grown.w === g.w
        @test_throws ArgumentError growto(g, (8, 7, 5))
        @test_throws ArgumentError growto(g, (9, 0, 5))

        converted = convert(Float32, g)
        @test converted isa RectangularGrid{1, Float32}
        @test converted.D₁[1] isa FDGrids.DiffMatrix
        @test converted.D₁⁺[1] isa FDGrids.AdjointDiffMatrix
        @test converted.D₂⁺[1] isa FDGrids.AdjointDiffMatrix
        @test eltype(converted.D₁⁺[1]) === Float32

        y = only(converted.xs)
        u = ftfield_from_inhomogeneous(converted, y .^ 3)
        v = ftfield_from_inhomogeneous(converted, (1 .- y .^ 2) .^ 2)
        Du, D⁺v = FTField(converted), FTField(converted)
        ddy!(Du, u)
        ddy!(D⁺v, v, AdjointDiscrete())
        @test dot(Du, v) ≈ dot(u, D⁺v) rtol=2e-5 atol=2e-6
    end

    @testset verbose=true "Operator lookup and multidirectional product weights        " begin
        channel = ChannelGrid(7, 17, 5; width=5)
        forward = @inferred NSEBase.derivative_matrix(channel, 1, Val(1), Forward())
        adjoint_matrix = @inferred NSEBase.derivative_matrix(channel, 1, Val(1), AdjointDiscrete())
        @test forward === channel.D₁[1]
        @test adjoint_matrix === channel.D₁⁺[1]
        @test typeof(forward) !== typeof(adjoint_matrix)
        @test NSEBase.derivative_matrix(channel, 1, Val(2), Forward()) === channel.D₂[1]
        @test_throws ArgumentError NSEBase.derivative_matrix(channel, 1, Val(3), Forward())
        @test_throws ArgumentError NSEBase.derivative_matrix(channel, 2, Val(1), Forward())

        cavity = LidDrivenCavity2DGrid(13, 11; xwidth=3, ywidth=5)
        @test cavity isa RectangularGrid{2}
        @test weights(cavity) isa RectangularProductWeights{2}
        @test weights(cavity) ≈ cavity.ws[1] * transpose(cavity.ws[2])
        @test typeof(cavity.D₁[1]) !== typeof(cavity.D₁[2])

        ws = (collect(1.0:3.0), collect(2.0:4.0), collect(3.0:5.0))
        product3 = RectangularProductWeights(ws)
        expected_weights = reshape(ws[1], :, 1, 1) .* reshape(ws[2], 1, :, 1) .*
                           reshape(ws[3], 1, 1, :)
        @test product3 isa RectangularProductWeights{3}
        @test size(product3) == (3, 3, 3)
        @test product3 ≈ expected_weights
        @test_throws BoundsError product3[0, 1, 1]
        @test_throws BoundsError product3[1, 4, 1]
    end

    @testset verbose=true "Low-level constructor validation                            " begin
        g = ChannelGrid(7, 13, 5; Nt=3, width=3)
        make(grid_size, axes, fft_dims, scales=g.scales) =
            RectangularGrid(g.xs, g.D₁, g.D₂, g.D₁⁺, g.D₂⁺, g.ws, scales,
                            grid_size, axes, fft_dims)

        @test_throws ArgumentError make((13, 7, 5, 3), (2, 2, 3, 4), (2, 3, 4))
        @test_throws ArgumentError make((13, 7, 5, 3), (2, 1, 3), (2, 3, 4))
        @test_throws ArgumentError make((13, 7, 5, 3), CHANNEL_AXES, ())
        @test_throws ArgumentError make((13, 7, 5, 3), CHANNEL_AXES, (2, 2, 4))
        @test_throws ArgumentError make((13, 7, 5, 3), CHANNEL_AXES, (2, 3, 5))
        @test_throws ArgumentError make((13, 7, 5, 3), CHANNEL_AXES, (2, 3), (1,))
        @test_throws ArgumentError make((13, 7, 5, 0), CHANNEL_AXES, (2, 3, 4))
        @test_throws ArgumentError make((13, 8, 5, 3), CHANNEL_AXES, (2, 3, 4))
        @test_throws ArgumentError make((13, 7, 5, 3), CHANNEL_AXES, (2, 3, 4), (1, 0, 2π))
        @test_throws ArgumentError make((13, 7, 5, 3), CHANNEL_AXES, (2, 3, 4), (1, Inf, 2π))
    end
end
