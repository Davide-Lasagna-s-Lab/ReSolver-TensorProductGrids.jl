module ReSolverTensorProductGrids

# Tensor-product grids for ReSolverFlowsBase: bounded directions with given points, derivative
# matrices and weights, times periodic Fourier directions. Concrete grids and flow cases live in
# ReSolverCases.

import FDGrids

import ReSolverFlowsBase: derivative_matrix, growto, points, wavenumber_scale, weights

using ReSolverFlowsBase: AbstractGrid, Direct, DiscreteAdjoint, fft_storage_dims, get_padded_size,
                         inhomogeneous_storage_dims

export TensorProductGrid, TensorProductWeights

include("tensorproductgrid.jl")

end
