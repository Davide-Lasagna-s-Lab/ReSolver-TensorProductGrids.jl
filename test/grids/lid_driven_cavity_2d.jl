@testset verbose=true "Two-dimensional lid-driven-cavity grid                      " begin
    N, Nt = 17, 19
    lim = (0, 1)
    g = LidDrivenCavity2DGrid(N; Nt, lim, dist=FDGrids.GaussLobattoGrid(), width=7)

    u(x, y, t) = bounded_profile(x, lim) * bounded_profile(y, lim) *
                  periodic_profile(2π * t)
    ux(x, y, t) = bounded_profile_d1(x, lim) * bounded_profile(y, lim) *
                   periodic_profile(2π * t)
    uy(x, y, t) = bounded_profile(x, lim) * bounded_profile_d1(y, lim) *
                   periodic_profile(2π * t)
    ut(x, y, t) = bounded_profile(x, lim) * bounded_profile(y, lim) *
                   2π * periodic_profile_d1(2π * t)
    Δu(x, y, t) = (bounded_profile_d2(x, lim) * bounded_profile(y, lim) +
                    bounded_profile(x, lim) * bounded_profile_d2(y, lim)) *
                   periodic_profile(2π * t)

    @testset verbose=true "Construction, layout, and product quadrature                " begin
        steady = LidDrivenCavity2DGrid(N)
        custom = LidDrivenCavity2DGrid(13; Nt=5, lim=(-1, 1), width=3)

        @test g isa LidDrivenCavity2DGrid
        @test weights(g) isa RectangularProductWeights{2}
        @test g isa AbstractLidDrivenCavity2DGrid
        @test g isa RectangularGrid{2}
        @test size(g) == (N, N, Nt)
        @test LID_DRIVEN_CAVITY_2D_AXES == (1, 2, nothing, 3)
        @test fft_storage_dims(g) == LID_DRIVEN_CAVITY_2D_FFT_ORDER
        @test inhomogeneous_storage_dims(g) == LID_DRIVEN_CAVITY_2D_INHOMOGENEOUS_DIMS
        @test fft_physical_dims(g) == (:t,)
        @test inhomogeneous_physical_dims(g) == (:x, :y)
        @test size(steady) == (N, N, 1)
        @test weights(g) ≈ g.ws[1] * transpose(g.ws[2])
        @test g.scales == (2π,)
        @test size(growto(g, (11,))) == (N, N, 11)
        @test all(object -> object === first(g.xs), g.xs)
        @test all(object -> object === first(g.D₁), g.D₁)
        @test all(object -> object === first(g.D₂), g.D₂)
        @test all(object -> object === first(g.D₁⁺), g.D₁⁺)
        @test all(object -> object === first(g.D₂⁺), g.D₂⁺)
        @test all(object -> object === first(g.ws), g.ws)
        @test extrema(custom.xs[1]) == (-1.0, 1.0)
        @test_throws MethodError LidDrivenCavity2DGrid(9, 9)
    end

    @testset verbose=true "Analytical derivatives and Laplacian                        " begin
        û = FFT(Field(g, u))
        for (derivative!, exact) in ((ddx!, ux), (ddy!, uy))
            @test derivative!(FTField(g), û) ≈ FFT(Field(g, exact)) atol=3e-11 rtol=3e-11
        end
        @test ddt!(FTField(g), û) ≈ FFT(Field(g, ut)) atol=3e-7 rtol=3e-7
        @test laplacian!(FTField(g), û) ≈ FFT(Field(g, Δu)) atol=3e-7 rtol=3e-7
    end

    @testset verbose=true "Analytical norms and homogeneous shifts                     " begin
        û = FFT(Field(g, u))
        exact_norm2 = Float64(bounded_profile_norm2(lim)^2 * PERIODIC_PROFILE_NORM2)
        velocity = VectorField(Field(g, u), Field(g, (x, y, t) -> 2u(x, y, t)),
                               Field(g, (x, y, t) -> 3u(x, y, t)))
        @test norm(û)^2 ≈ exact_norm2 rtol=3e-12
        @test norm(FFT(velocity))^2 ≈ 14exact_norm2 rtol=3e-12

        st = 0.23
        shifted(x, y, t) = bounded_profile(x, lim) * bounded_profile(y, lim) *
                            periodic_profile(2π * (t + st))
        @test shift!(copy(û), (st,)) ≈ FFT(Field(g, shifted)) atol=3e-7 rtol=3e-7
    end

    @testset verbose=true "Quadrature-weighted discrete adjoints                       " begin
        v(x, y, t) = dual_bounded_profile(x, lim) * dual_bounded_profile(y, lim) *
                      periodic_profile(2π * t + 0.4)
        û, v̂ = FFT(Field(g, u)), FFT(Field(g, v))

        for derivative! in (ddx!, ddy!, ddt!)
            Du = derivative!(FTField(g), û)
            D⁺v = derivative!(FTField(g), v̂, AdjointDiscrete())
            @test dot(Du, v̂) ≈ dot(û, D⁺v) atol=5e-12 rtol=5e-12
        end
        Δû = laplacian!(FTField(g), û)
        Δ⁺v = laplacian!(FTField(g), v̂, AdjointDiscrete())
        @test dot(Δû, v̂) ≈ dot(û, Δ⁺v) atol=5e-10 rtol=5e-10
    end
end
