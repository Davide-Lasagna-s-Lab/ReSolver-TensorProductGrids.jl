# =============================================================================================== #
# Rectangular-grid weights and construction helpers                                               #
# =============================================================================================== #

# =============================================================================================== #
# Lazy tensor-product quadrature weights                                                          #
# =============================================================================================== #

"""
    RectangularProductWeights{N, T, V} <: AbstractArray{T, N}
    RectangularProductWeights(ws::NTuple{N, V}) where {N, T, V<:AbstractVector{T}}

Lazy tensor product of `N` one-dimensional quadrature vectors. An entry is evaluated as

```text
W[i₁, …, iₙ] = prod(ws[n][iₙ] for n in 1:N)
```

so multidimensional quadrature never allocates a materialized product array. The wrapper implements
the Cartesian `AbstractArray` interface required by NSEBase's weighted inner products and Galerkin
operators. `ws` follows increasing inhomogeneous storage-dimension order.

# Arguments

- `ws`: `N` one-dimensional quadrature vectors with common scalar and vector types.

# Returns

An `N`-dimensional lazy array whose size and axes are inherited from the vectors in `ws`.
"""
struct RectangularProductWeights{N, T, V<:AbstractVector{T}} <: AbstractArray{T, N}
    ws :: NTuple{N, V}
end

Base.IndexStyle(::Type{<:RectangularProductWeights}) = IndexCartesian()
Base.size(A::RectangularProductWeights) = map(length, A.ws)
Base.axes(A::RectangularProductWeights) = map(w -> axes(w, 1), A.ws)

Base.@propagate_inbounds function Base.getindex(
        A::RectangularProductWeights{N}, I::Vararg{Int, N}) where {N}
    @boundscheck checkbounds(A, I...)
    return prod(ntuple(i -> (@inbounds A.ws[i][I[i]]), N))
end

# =============================================================================================== #
# Coordinate and storage helpers                                                                  #
# =============================================================================================== #

_equal_points(N::Int, L::T) where {T} = range(zero(T), step=L / T(N), length=N)

_rectangular_padded_size(g::AbstractGrid{<:Any, D}) where {D} =
    ntuple(dim -> dim in fft_storage_dims(g) ? (cld(3 * size(g, dim), 2) | 1) : size(g, dim), D)

function _rectangular_storage_size(
        g::AbstractGrid{T, D}, homogeneous_size::NTuple{N, Int}) where {T, D, N}
    fft_dims = fft_storage_dims(g)
    N == length(fft_dims) || throw(ArgumentError("expected $(length(fft_dims)) Fourier sizes"))
    return ntuple(dim -> dim in fft_dims ? homogeneous_size[findfirst(==(dim), fft_dims)] :
                  size(g, dim), D)
end

function _rectangular_points(
        g::AbstractGrid{T, D}, homogeneous_size::NTuple{N, Int}, xs::Tuple) where {T, D, N}
    storage_size = _rectangular_storage_size(g, homogeneous_size)
    inhomogeneous_dims = inhomogeneous_storage_dims(g)
    length(xs) == length(inhomogeneous_dims) ||
        throw(ArgumentError("expected one collocation vector per inhomogeneous direction"))

    return ntuple(D) do dim
        index = findfirst(==(dim), inhomogeneous_dims)
        values = isnothing(index) ?
                 _equal_points(storage_size[dim], T(2π) / wavenumber_scale(g, dim)) : xs[index]
        reshape(values, ntuple(d -> d == dim ? length(values) : 1, D))
    end
end

# =============================================================================================== #
# Grid-data conversion and FDGrids construction                                                   #
# =============================================================================================== #

_convert_grid_data(::Type{T}, x) where {T} = eltype(x) === T ? x : T.(x)

function _convert_grid_data(::Type{T}, A::FDGrids.AdjointDiffMatrix) where {T}
    eltype(A) === T && return A
    parent = _convert_grid_data(T, A.parent)
    return FDGrids.AdjointDiffMatrix(parent, Vector{T}(A.coeffs))
end

"""
    _fd_direction(N, lim, dist, width, T)

Build all numerical data for one bounded direction using FDGrids. The `N` collocation points span
the interval `lim` according to `dist`. First- and second-derivative matrices use the common odd
stencil width `width`, and their discrete adjoints are formed with the FDGrids quadrature weights.
Points and weights are materialized as `Vector{T}`; FDGrids constructs derivative matrices with
element type `T` directly.

# Arguments

- `N`: number of bounded collocation points.
- `lim`: two-endpoint physical interval.
- `dist`: FDGrids point distribution on that interval.
- `width`: finite-difference stencil width for both derivative orders.
- `T`: real scalar type for points, operators, and weights.

# Returns

A tuple `(x, D₁, D₂, D₁⁺, D₂⁺, w)` containing the collocation points, first- and
second-derivative matrices, their quadrature-weighted discrete adjoints, and quadrature weights.
"""
function _fd_direction(N::Int, lim::NTuple{2, <:Real}, dist::FDGrids.AbstractGridDistribution,
                       width::Int, ::Type{T}) where {T<:Real}
    fdgrid = FDGrids.grid(N, lim[1], lim[2], dist)
    x, w = Vector{T}(fdgrid.xs), Vector{T}(fdgrid.ws)
    D₁ = FDGrids.DiffMatrix(x, width, 1; eltype=T)
    D₂ = FDGrids.DiffMatrix(x, width, 2; eltype=T)
    return x, D₁, D₂, adjoint(D₁, w), adjoint(D₂, w), w
end
