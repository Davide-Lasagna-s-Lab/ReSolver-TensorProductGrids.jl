@testset verbose=true "Three-dimensional lid-driven-cavity grid                    " begin
    Nx, Ny, Nz, Nt = 17, 17, 17, 19
    xlim, ylim, zlim = (0, 1), (0, 1), (0, 1)
    g = LidDrivenCavity3DGrid(Nx, Ny, Nz; Nt, xlim, ylim, zlim,
                              dist=FDGrids.GaussLobattoGrid(), width=7)

    u(x, y, z, t) = bounded_profile(x, xlim) * bounded_profile(y, ylim) *
                     bounded_profile(z, zlim) * periodic_profile(2π * t)
    ux(x, y, z, t) = bounded_profile_d1(x, xlim) * bounded_profile(y, ylim) *
                      bounded_profile(z, zlim) * periodic_profile(2π * t)
    uy(x, y, z, t) = bounded_profile(x, xlim) * bounded_profile_d1(y, ylim) *
                      bounded_profile(z, zlim) * periodic_profile(2π * t)
    uz(x, y, z, t) = bounded_profile(x, xlim) * bounded_profile(y, ylim) *
                      bounded_profile_d1(z, zlim) * periodic_profile(2π * t)
    ut(x, y, z, t) = bounded_profile(x, xlim) * bounded_profile(y, ylim) *
                      bounded_profile(z, zlim) * 2π * periodic_profile_d1(2π * t)
    Δu(x, y, z, t) = (bounded_profile_d2(x, xlim) * bounded_profile(y, ylim) *
        bounded_profile(z, zlim) + bounded_profile(x, xlim) * bounded_profile_d2(y, ylim) *
        bounded_profile(z, zlim) + bounded_profile(x, xlim) * bounded_profile(y, ylim) *
        bounded_profile_d2(z, zlim)) * periodic_profile(2π * t)

    @testset verbose=true "Construction, layout, and product quadrature                " begin
        independent = LidDrivenCavity3DGrid(13, 11, 15; Nt=5, xlim=(-1, 1), ylim=(0, 2),
                                             zlim=(-2, 2), xwidth=3, ywidth=5, zwidth=7)

        @test g isa LidDrivenCavity3DGrid
        @test g isa AbstractLidDrivenCavity3DGrid
        @test g isa RectangularGrid{3}
        @test size(g) == (Nx, Ny, Nz, Nt)
        @test LID_DRIVEN_CAVITY_3D_AXES == (1, 2, 3, 4)
        @test fft_storage_dims(g) == LID_DRIVEN_CAVITY_3D_FFT_ORDER
        @test inhomogeneous_storage_dims(g) == LID_DRIVEN_CAVITY_3D_INHOMOGENEOUS_DIMS
        @test fft_physical_dims(g) == (:t,)
        @test inhomogeneous_physical_dims(g) == (:x, :y, :z)
        expected_weights = reshape(g.ws[1], :, 1, 1) .* reshape(g.ws[2], 1, :, 1) .*
                           reshape(g.ws[3], 1, 1, :)
        @test weights(g) isa RectangularProductWeights{3}
        @test weights(g) ≈ expected_weights
        @test g.scales == (2π,)
        @test size(growto(g, (11,))) == (Nx, Ny, Nz, 11)
        @test extrema(independent.xs[1]) == (-1.0, 1.0)
        @test extrema(independent.xs[2]) == (0.0, 2.0)
        @test extrema(independent.xs[3]) == (-2.0, 2.0)
        @test allunique(typeof.(independent.D₁))
        @test_throws MethodError LidDrivenCavity3DGrid(9, 9)
    end

    @testset verbose=true "Analytical derivatives and Laplacian                        " begin
        û = FFT(Field(g, u))
        for (derivative!, exact) in ((ddx!, ux), (ddy!, uy), (ddz!, uz))
            @test derivative!(FTField(g), û) ≈ FFT(Field(g, exact)) atol=3e-11 rtol=3e-11
        end
        @test ddt!(FTField(g), û) ≈ FFT(Field(g, ut)) atol=3e-7 rtol=3e-7
        @test laplacian!(FTField(g), û) ≈ FFT(Field(g, Δu)) atol=3e-7 rtol=3e-7
    end

    @testset verbose=true "Analytical norms and homogeneous shifts                     " begin
        û = FFT(Field(g, u))
        exact_norm2 = Float64(bounded_profile_norm2(xlim) * bounded_profile_norm2(ylim) *
                              bounded_profile_norm2(zlim) * PERIODIC_PROFILE_NORM2)
        velocity = VectorField(Field(g, u), Field(g, (x, y, z, t) -> 2u(x, y, z, t)),
                               Field(g, (x, y, z, t) -> 3u(x, y, z, t)))
        @test norm(û)^2 ≈ exact_norm2 rtol=3e-12
        @test norm(FFT(velocity))^2 ≈ 14exact_norm2 rtol=3e-12

        st = 0.23
        shifted(x, y, z, t) = bounded_profile(x, xlim) * bounded_profile(y, ylim) *
                               bounded_profile(z, zlim) * periodic_profile(2π * (t + st))
        @test shift!(copy(û), (st,)) ≈ FFT(Field(g, shifted)) atol=3e-7 rtol=3e-7
    end

    @testset verbose=true "Quadrature-weighted discrete adjoints                       " begin
        v(x, y, z, t) = dual_bounded_profile(x, xlim) * dual_bounded_profile(y, ylim) *
                         dual_bounded_profile(z, zlim) * periodic_profile(2π * t + 0.4)
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
