module ReSolverRectangularGrids

using LinearAlgebra

import FDGrids
import NSEBase
import NSEBase: growto, points, wavenumber_scale, weights

using NSEBase: AbstractGrid, AdjointContinuous, AdjointDiscrete, CartesianPrimitive2D,
               CartesianPrimitive3D, CompoundForcing,
               FFTW, FTField, Forward, NoForce, VectorField, WaveNumberVector,
               construct_equations, fft_storage_dims, grid, inhomogeneous_storage_dims

export RectangularGrid, RectangularProductWeights

export CHANNEL_AXES, CHANNEL_FFT_ORDER, CHANNEL_INHOMOGENEOUS_DIMS
export AbstractChannelGrid, ChannelGrid

export LID_DRIVEN_CAVITY_2D_AXES, LID_DRIVEN_CAVITY_2D_FFT_ORDER
export LID_DRIVEN_CAVITY_2D_INHOMOGENEOUS_DIMS
export AbstractLidDrivenCavity2DGrid, LidDrivenCavity2DGrid

export LID_DRIVEN_CAVITY_3D_AXES, LID_DRIVEN_CAVITY_3D_FFT_ORDER
export LID_DRIVEN_CAVITY_3D_INHOMOGENEOUS_DIMS
export AbstractLidDrivenCavity3DGrid, LidDrivenCavity3DGrid

export SQUARE_DUCT_AXES, SQUARE_DUCT_FFT_ORDER, SQUARE_DUCT_INHOMOGENEOUS_DIMS
export AbstractSquareDuctGrid, SquareDuctGrid

export ConstantBodyForce, CoriolisForce
export plane_couette_base, plane_poiseuille_base
export PlaneCouetteFlow, PlanePoiseuilleFlow
export LidDrivenCavity2DFlow, LidDrivenCavity3DFlow, SquareDuctFlow

include("helpers.jl")
include("rectangular.jl")
include("forcings.jl")
include("grids/channel.jl")
include("grids/lid_driven_cavity_2d.jl")
include("grids/lid_driven_cavity_3d.jl")
include("grids/square_duct.jl")
include("cases/channel.jl")
include("cases/lid_driven_cavity_2d.jl")
include("cases/lid_driven_cavity_3d.jl")
include("cases/square_duct.jl")

end
