@testset verbose=true "Square-duct grid                                            " begin
    N, Nz, Nt, α = 17, 19, 19, 1.25
    lim = (0, 1)
    g = SquareDuctGrid(N, Nz; Nt, α, dist=FDGrids.GaussLobattoGrid(), width=7)

    u(x, y, z, t) = bounded_profile(x, lim) * bounded_profile(y, lim) *
                     periodic_profile(α * z) * periodic_profile(t)
    ux(x, y, z, t) = bounded_profile_d1(x, lim) * bounded_profile(y, lim) *
                      periodic_profile(α * z) * periodic_profile(t)
    uy(x, y, z, t) = bounded_profile(x, lim) * bounded_profile_d1(y, lim) *
                      periodic_profile(α * z) * periodic_profile(t)
    uz(x, y, z, t) = bounded_profile(x, lim) * bounded_profile(y, lim) *
                      α * periodic_profile_d1(α * z) * periodic_profile(t)
    ut(x, y, z, t) = bounded_profile(x, lim) * bounded_profile(y, lim) *
                      periodic_profile(α * z) * periodic_profile_d1(t)
    Δu(x, y, z, t) = ((bounded_profile_d2(x, lim) * bounded_profile(y, lim) +
                        bounded_profile(x, lim) * bounded_profile_d2(y, lim)) *
                       periodic_profile(α * z) + bounded_profile(x, lim) *
                       bounded_profile(y, lim) * α^2 * periodic_profile_d2(α * z)) *
                      periodic_profile(t)

    @testset verbose=true "Construction, layout, and product quadrature                " begin
        steady = SquareDuctGrid(N, Nz)

        @test g isa SquareDuctGrid
        @test weights(g) isa RectangularProductWeights{2}
        @test g isa AbstractSquareDuctGrid
        @test g isa RectangularGrid{2}
        @test size(g) == (N, N, Nz, Nt)
        @test SQUARE_DUCT_AXES == (1, 2, 3, 4)
        @test fft_storage_dims(g) == SQUARE_DUCT_FFT_ORDER
        @test inhomogeneous_storage_dims(g) == SQUARE_DUCT_INHOMOGENEOUS_DIMS
        @test fft_physical_dims(g) == (:z, :t)
        @test inhomogeneous_physical_dims(g) == (:x, :y)
        @test size(steady) == (N, N, Nz, 1)
        @test weights(g) ≈ g.ws[1] * transpose(g.ws[2])
        @test g.scales == (α, 1)
        @test size(growto(g, (11, 13))) == (N, N, 11, 13)
        @test g.xs[1] === g.xs[2]
        @test g.D₁[1] === g.D₁[2]
        @test g.D₂[1] === g.D₂[2]
        @test g.D₁⁺[1] === g.D₁⁺[2]
        @test g.D₂⁺[1] === g.D₂⁺[2]
        @test g.ws[1] === g.ws[2]
        @test_throws MethodError SquareDuctGrid(9, 9, 9)
    end

    @testset verbose=true "Analytical derivatives and Laplacian                        " begin
        û = FFT(Field(g, u))
        for (derivative!, exact) in ((ddx!, ux), (ddy!, uy))
            @test derivative!(FTField(g), û) ≈ FFT(Field(g, exact)) atol=3e-11 rtol=3e-11
        end
        for (derivative!, exact) in ((ddz!, uz), (ddt!, ut))
            @test derivative!(FTField(g), û) ≈ FFT(Field(g, exact)) atol=3e-7 rtol=3e-7
        end
        @test laplacian!(FTField(g), û) ≈ FFT(Field(g, Δu)) atol=3e-7 rtol=3e-7
    end

    @testset verbose=true "Analytical norms and homogeneous shifts                     " begin
        û = FFT(Field(g, u))
        exact_norm2 = Float64(bounded_profile_norm2(lim)^2 * PERIODIC_PROFILE_NORM2^2)
        velocity = VectorField(Field(g, u), Field(g, (x, y, z, t) -> 2u(x, y, z, t)),
                               Field(g, (x, y, z, t) -> 3u(x, y, z, t)))
        @test norm(û)^2 ≈ exact_norm2 rtol=3e-12
        @test norm(FFT(velocity))^2 ≈ 14exact_norm2 rtol=3e-12

        sz, st = 1.11, 0.23
        shifted(x, y, z, t) = bounded_profile(x, lim) * bounded_profile(y, lim) *
                               periodic_profile(α * (z + sz)) * periodic_profile(t + st)
        @test shift!(copy(û), (sz, st)) ≈ FFT(Field(g, shifted)) atol=3e-7 rtol=3e-7
    end

    @testset verbose=true "Quadrature-weighted discrete adjoints                       " begin
        v(x, y, z, t) = dual_bounded_profile(x, lim) * dual_bounded_profile(y, lim) *
                         periodic_profile(α * z + 0.3) * periodic_profile(t + 0.4)
        û, v̂ = FFT(Field(g, u)), FFT(Field(g, v))

        for derivative! in (ddx!, ddy!, ddz!, ddt!)
            Du = derivative!(FTField(g), û)
            D⁺v = derivative!(FTField(g), v̂, AdjointDiscrete())
            @test dot(Du, v̂) ≈ dot(û, D⁺v) atol=5e-12 rtol=5e-12
        end
        Δû = laplacian!(FTField(g), û)
        Δ⁺v = laplacian!(FTField(g), v̂, AdjointDiscrete())
        @test dot(Δû, v̂) ≈ dot(û, Δ⁺v) atol=5e-10 rtol=5e-10
    end
end
