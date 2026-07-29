@testset verbose=true "Runnable examples                                           " begin
    examples = normpath(joinpath(@__DIR__, "..", "examples"))

    @testset verbose=true "Channel example                                             " begin
        result = include(joinpath(examples, "channel.jl"))
        @test result.grid_size == (25, 19, 19, 1)
        @test result.derivative_error < 5e-7
        @test result.norm_error < 5e-11
    end

    @testset verbose=true "Two-dimensional cavity example                              " begin
        result = include(joinpath(examples, "lid_driven_cavity_2d.jl"))
        @test result.grid_size == (17, 17, 19)
        @test result.derivative_error < 5e-7
        @test result.norm_error < 5e-11
    end

    @testset verbose=true "Three-dimensional cavity example                            " begin
        result = include(joinpath(examples, "lid_driven_cavity_3d.jl"))
        @test result.grid_size == (15, 15, 15, 19)
        @test result.derivative_error < 5e-7
        @test result.norm_error < 5e-11
    end

    @testset verbose=true "Square-duct example                                         " begin
        result = include(joinpath(examples, "square_duct.jl"))
        @test result.grid_size == (17, 17, 19, 1)
        @test result.derivative_error < 5e-7
        @test result.norm_error < 5e-11
    end
end
