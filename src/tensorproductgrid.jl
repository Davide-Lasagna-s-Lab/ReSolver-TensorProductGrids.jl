# Tensor-product grids: one to three bounded finite-difference directions times one or more
# periodic Fourier directions, implementing the ReSolverFlowsBase.AbstractGrid interface.
#
# Bounded data are stored as tuples, one entry per inhomogeneous storage dimension, in increasing
# storage order; wavenumber scales follow ReSolverFlowsBase.spatial_fft_storage_dims, the phase
# direction s ∈ [0, 2π) having scale one by construction. The grid owns
# geometry, quadrature and differentiation; boundary conditions belong to the basis or residual.


# ============================================================================================== #
# Lazy tensor-product quadrature weights                                                         #
# ============================================================================================== #

"""
    TensorProductWeights(ws)

Lazy tensor product of the one-dimensional quadrature vectors in the tuple `ws`: the entry at
`(i, j, …)` is `ws[1][i] * ws[2][j] * …`. Nothing is allocated.
"""
struct TensorProductWeights{N, T, W<:Tuple} <: AbstractArray{T, N}
    ws::W # one-dimensional quadrature vectors

    TensorProductWeights(ws::Tuple) = new{length(ws), eltype(first(ws)), typeof(ws)}(ws)
end

Base.size(A::TensorProductWeights) = map(length, A.ws)

Base.@propagate_inbounds function Base.getindex(A::TensorProductWeights{N},
                                                I::Vararg{Int, N}) where {N}
    @boundscheck checkbounds(A, I...)
    return prod(ntuple(i -> (@inbounds A.ws[i][I[i]]), Val(N)))
end


# ============================================================================================== #
# Grid                                                                                           #
# ============================================================================================== #

"""
    TensorProductGrid(xs, D₁, D₂, D₁⁺, D₂⁺, ws, wavenumber_scales, grid_size, axes, fft_dims)

Tensor-product grid with bounded directions discretised by the tuples of collocation points `xs`,
first- and second-derivative matrices `D₁`, `D₂`, their quadrature-weighted discrete adjoints
`D₁⁺`, `D₂⁺` and quadrature weights `ws`, one entry per bounded direction in increasing storage
order. `wavenumber_scales` holds the scale `2π/L` of each spatial Fourier direction, in `fft_dims`
order; the time phase `s ∈ [0, 2π)` has scale one and no entry. `grid_size` is the array size in
storage order and `axes` maps the coordinates `(x1, x2, x3, s)` to storage dimensions, with
`nothing` for an absent coordinate. The scalar type is that of the
weights; use `convert` to change it.
"""
struct TensorProductGrid{T, S, D, AXES, FFT_DIMS, X, D1, D2, A1, A2, W, C} <:
       AbstractGrid{T, D, AXES, FFT_DIMS}
                   xs::X  # collocation points of each bounded direction
                   D₁::D1 # first-derivative matrices
                   D₂::D2 # second-derivative matrices
                  D₁⁺::A1 # discrete adjoints of D₁
                  D₂⁺::A2 # discrete adjoints of D₂
                   ws::W  # one-dimensional quadrature weights
    wavenumber_scales::C  # scale 2π/L of each spatial Fourier direction, in transform order

    function TensorProductGrid(               xs::Tuple,
                                              D₁::Tuple,
                                              D₂::Tuple,
                                             D₁⁺::Tuple,
                                             D₂⁺::Tuple,
                                              ws::Tuple,
                               wavenumber_scales::Tuple,
                                       grid_size::Tuple,
                                            axes::Tuple,
                                        fft_dims::Tuple)
        T = eltype(first(ws))
        k = map(T, wavenumber_scales)

        return new{T, grid_size, length(grid_size), axes, fft_dims,
                   typeof(xs), typeof(D₁), typeof(D₂), typeof(D₁⁺), typeof(D₂⁺),
                   typeof(ws), typeof(k)}(xs, D₁, D₂, D₁⁺, D₂⁺, ws, k)
    end
end


# ============================================================================================== #
# AbstractGrid interface                                                                         #
# ============================================================================================== #

Base.size(::TensorProductGrid{T, S}) where {T, S} = S

Base.eltype(::Type{<:TensorProductGrid{T}}) where {T} = T
Base.eltype(::TensorProductGrid{T})         where {T} = T

weights(g::TensorProductGrid) = TensorProductWeights(g.ws)

# bounded directions and the phase s: one; spatial Fourier directions: 2π/L
function wavenumber_scale(g::TensorProductGrid, storage_dim::Int)
    i = findfirst(==(storage_dim), spatial_fft_storage_dims(g))
    return isnothing(i) ? one(eltype(g)) : g.wavenumber_scales[i]
end

"""
    points(g::TensorProductGrid; dealias=false)
    points(g::TensorProductGrid, homogeneous_size)

One coordinate array per storage dimension, shaped for broadcasting. Bounded coordinates are the
collocation points; a Fourier coordinate of scale `k` is equispaced on `[0, 2π/k)`, the phase on
`[0, 2π)`. The Fourier
resolutions are those of `g`, their 3/2-rule padded sizes with `dealias=true`, or
`homogeneous_size` in Fourier transform order.
"""
function points(g::TensorProductGrid; dealias::Bool=false)
    storage_size = dealias ? get_padded_size(size(g), fft_storage_dims(g)) : size(g)
    return _points(g, storage_size)
end

points(g::TensorProductGrid, homogeneous_size::Tuple) =
    _points(g, _storage_size(g, homogeneous_size))

"""
    growto(g::TensorProductGrid, homogeneous_size)

The same grid with Fourier resolutions `homogeneous_size`, in transform order; all bounded data
are shared with `g`.
"""
growto(g::TensorProductGrid, homogeneous_size::Tuple) =
    TensorProductGrid(g.xs, g.D₁, g.D₂, g.D₁⁺, g.D₂⁺, g.ws, g.wavenumber_scales,
                      _storage_size(g, homogeneous_size),
                      _axes(g),
                      fft_storage_dims(g))

"""
    convert(T, g::TensorProductGrid)

The grid with points, operators, weights and scales converted to the scalar type `T`.
"""
function Base.convert(::Type{T}, g::TensorProductGrid) where {T<:Real}
    eltype(g) === T && return g

    return TensorProductGrid(map(x -> _convert(T, x), g.xs),
                             map(A -> _convert(T, A), g.D₁),
                             map(A -> _convert(T, A), g.D₂),
                             map(A -> _convert(T, A), g.D₁⁺),
                             map(A -> _convert(T, A), g.D₂⁺),
                             map(w -> _convert(T, w), g.ws),
                             g.wavenumber_scales,
                             size(g),
                             _axes(g),
                             fft_storage_dims(g))
end


# ============================================================================================== #
# Finite-difference operators                                                                    #
# ============================================================================================== #

# Operators of a bounded storage dimension: the derivative matrices, or their discrete adjoints.
derivative_matrix(g::TensorProductGrid, dim::Integer, ::Val{1}, ::Direct)          = g.D₁[_bounded(g, dim)]
derivative_matrix(g::TensorProductGrid, dim::Integer, ::Val{2}, ::Direct)          = g.D₂[_bounded(g, dim)]
derivative_matrix(g::TensorProductGrid, dim::Integer, ::Val{1}, ::DiscreteAdjoint) = g.D₁⁺[_bounded(g, dim)]
derivative_matrix(g::TensorProductGrid, dim::Integer, ::Val{2}, ::DiscreteAdjoint) = g.D₂⁺[_bounded(g, dim)]

# ============================================================================================== #
# Internals                                                                                      #
# ============================================================================================== #

# position of a storage dimension among the bounded ones
_bounded(g::TensorProductGrid, dim::Integer) = findfirst(==(dim), inhomogeneous_storage_dims(g))

# the coordinate map (x1, x2, x3, s) -> storage dimension, from the type
_axes(::AbstractGrid{T, D, AXES}) where {T, D, AXES} = AXES

# storage size with the given Fourier resolutions, in transform order
function _storage_size(g::TensorProductGrid, homogeneous_size::Tuple)
    return ntuple(length(size(g))) do dim
        i = findfirst(==(dim), fft_storage_dims(g))
        isnothing(i) ? size(g, dim) : homogeneous_size[i]
    end
end

# one broadcast-shaped coordinate array per storage dimension
function _points(g::TensorProductGrid, storage_size::Tuple)
    T = eltype(g)
    D = length(size(g))

    return ntuple(D) do dim
        i = findfirst(==(dim), inhomogeneous_storage_dims(g))
        N = storage_size[dim]

        # bounded: collocation points; Fourier: equispaced over one period
        L = T(2π) / wavenumber_scale(g, dim)
        x = isnothing(i) ? range(zero(T), step=L / T(N), length=N) : g.xs[i]

        reshape(x, ntuple(d -> d == dim ? length(x) : 1, D))
    end
end

# Scalar-type conversion of the grid data. Vectors are converted elementwise; FDGrids operators
# keep their banded structure. Any other matrix type has no conversion defined: rather than
# silently falling back to a dense matrix, converting it is an error.
_convert(::Type{T}, x::AbstractVector) where {T} = eltype(x) === T ? x : T.(x)

_convert(::Type{T}, D::FDGrids.DiffMatrix) where {T} = eltype(D) === T ? D : T.(D)

function _convert(::Type{T}, A::FDGrids.AdjointDiffMatrix) where {T}
    eltype(A) === T && return A
    return FDGrids.AdjointDiffMatrix(_convert(T, A.parent), Vector{T}(A.coeffs))
end

function _convert(::Type{T}, A::AbstractMatrix) where {T}
    eltype(A) === T && return A
    throw(ArgumentError("no conversion to $T defined for $(typeof(A)); " *
                        "convert the operators before building the grid"))
end
