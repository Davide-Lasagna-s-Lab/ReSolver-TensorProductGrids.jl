# =============================================================================================== #
# Rectangular tensor-product grids                                                               #
# =============================================================================================== #

"""
    RectangularGrid{NI, T, S, D, AXES, FFT_DIMS}

A concrete tensor-product Cartesian grid with `NI` bounded finite-difference directions and one
or more periodic Fourier directions. `RectangularGrid` implements the `NSEBase.AbstractGrid`
interface used by fields, transforms, derivatives, weighted inner products, and resolution growth.

# Layout

`S` is the physical array size in storage order, `AXES` maps physical `(x, y, z, t)` coordinates to
storage dimensions, and `FFT_DIMS` lists homogeneous storage dimensions in transform order. The
first Fourier dimension uses the real-to-complex transform. Direction-dependent tuples follow
`NSEBase.inhomogeneous_storage_dims`, while `scales` follows `NSEBase.fft_storage_dims`.

For example, a channel stored as `(y, x, z, t)` has `AXES == (2, 1, 3, 4)`,
`FFT_DIMS == (2, 3, 4)`, and one inhomogeneous tuple entry for physical `y`.

# Stored data

- `xs`: collocation points for each bounded direction.
- `D₁`, `D₂`: first- and second-derivative operators.
- `D₁⁺`, `D₂⁺`: their quadrature-consistent discrete adjoints.
- `ws`: one-dimensional quadrature weights.
- `w`: the cached tensor product of `ws`, returned by `NSEBase.weights`.
- `scales`: factors mapping integer Fourier modes to physical wavenumbers.

Every operator is an FDGrids matrix and supports dimension-wise multiplication with
`LinearAlgebra.mul!`. The grid owns geometry, quadrature, and differentiation; boundary conditions
remain the responsibility of a basis, projection, or residual formulation.

Most users should construct one of the supplied cases—[`ChannelGrid`](@ref),
[`LidDrivenCavity2DGrid`](@ref), [`LidDrivenCavity3DGrid`](@ref), or [`SquareDuctGrid`](@ref)—instead
of calling the low-level constructor directly.
"""
struct RectangularGrid{NI, T, S, D, AXES, FFT_DIMS, XS, D1, D2, A1, A2, WS, WP, NFFT} <: AbstractGrid{T, D, AXES, FFT_DIMS}
    xs     :: NTuple{NI, XS}
    D₁     :: NTuple{NI, D1}
    D₂     :: NTuple{NI, D2}
    D₁⁺    :: NTuple{NI, A1}
    D₂⁺    :: NTuple{NI, A2}
    ws     :: NTuple{NI, WS}
    w      :: WP
    scales :: NTuple{NFFT, T}
end

"""
    RectangularGrid(xs, D₁, D₂, D₁⁺, D₂⁺, ws, scales,
                    grid_size, axes, fft_dims, T=Float64)

Construct a rectangular grid from existing one-dimensional FDGrids discretisations.

The common length of `xs`, `D₁`, `D₂`, `D₁⁺`, `D₂⁺`, and `ws` becomes `NI`. Those
tuples are ordered by increasing inhomogeneous storage dimension. `scales` and `fft_dims` have the
same length and use transform order. A periodic interval of length ``L`` normally uses scale
``2π/L``.

Every Fourier scale must be a positive finite real number. `grid_size` is the physical array shape
in storage order. `axes` is the four-entry physical map
`(x_dim, y_dim, z_dim, t_dim)`, using `nothing` for an absent coordinate. All present entries must be
a permutation of `1:length(grid_size)`. A present time coordinate must be Fourier transformed.

The constructor validates all sizes, layouts, matrices, and quadrature vectors. Numerical data are
converted to `T` only when required; same-type objects are retained, so deliberately shared
operators remain shared. Supplied adjoints are trusted and are never recomputed.

# Arguments

- `xs`: tuple of collocation vectors, ordered by increasing bounded storage dimension.
- `D₁`, `D₂`: tuples of first- and second-derivative matrices in the same order as `xs`.
- `D₁⁺`, `D₂⁺`: tuples of quadrature-weighted discrete-adjoint derivative matrices.
- `ws`: tuple of one-dimensional quadrature-weight vectors.
- `scales`: positive finite physical wavenumber scale for each entry of `fft_dims`, in transform
  order.
- `grid_size`: physical array shape in storage order.
- `axes`: physical `(x,y,z,t)` to storage-dimension map, with `nothing` for absent coordinates.
- `fft_dims`: homogeneous storage dimensions in transform order.
- `T`: real scalar type used by the constructed grid.

# Returns

A validated [`RectangularGrid`](@ref) with lazy tensor-product quadrature weights.

# Example

```julia
g = RectangularGrid((y,), (Dy,), (Dy2,), (Dya,), (Dy2a,), (wy,),
                    (α, β, 1), (length(y), Nx, Nz, Nt),
                    (2, 1, 3, 4), (2, 3, 4))
```
"""
function RectangularGrid(xs::Tuple{Vararg{Any, NI}}, D₁::Tuple{Vararg{Any, NI}},
                         D₂::Tuple{Vararg{Any, NI}}, D₁⁺::Tuple{Vararg{Any, NI}},
                         D₂⁺::Tuple{Vararg{Any, NI}}, ws::Tuple{Vararg{Any, NI}}, scales::Tuple,
                         grid_size::NTuple{D, Int}, axes::Tuple, fft_dims::Tuple,
                         ::Type{T}=Float64) where {NI, D, T<:Real}
    # Validate the number, storage sizes, and basic container types of all bounded directions.
    1 <= NI <= 3 || throw(ArgumentError("expected one to three inhomogeneous directions"))
    all(>(0), grid_size) || throw(ArgumentError("grid sizes must be positive"))
    all(x -> x isa AbstractVector, xs) || throw(ArgumentError("collocation data must be vectors"))
    all(w -> w isa AbstractVector, ws) || throw(ArgumentError("quadrature data must be vectors"))
    all(A -> A isa AbstractMatrix, (D₁..., D₂..., D₁⁺..., D₂⁺...)) ||
        throw(ArgumentError("derivative operators must be matrices"))

    # Validate the physical-coordinate map and the ordered set of Fourier storage dimensions.
    length(axes) == 4 || throw(ArgumentError("axes must contain the four Cartesian coordinate slots"))
    present_axes = filter(x -> !isnothing(x), axes)
    all(x -> x isa Int, present_axes) && sort(collect(present_axes)) == collect(1:D) ||
        throw(ArgumentError("nonempty axes must be a permutation of the $D storage dimensions"))
    isempty(fft_dims) && throw(ArgumentError("at least one Fourier direction is required"))
    all(dim -> dim isa Int && 1 <= dim <= D, fft_dims) ||
        throw(ArgumentError("Fourier dimensions must lie in 1:$D"))
    allunique(fft_dims) || throw(ArgumentError("Fourier dimensions must be unique"))
    length(scales) == length(fft_dims) ||
        throw(ArgumentError("expected one wavenumber scale per Fourier direction"))
    all(scale -> scale isa Real && isfinite(scale) && scale > 0, scales) ||
        throw(ArgumentError("wavenumber scales must be positive finite real numbers"))
    isnothing(axes[4]) || axes[4] in fft_dims ||
        throw(ArgumentError("a present time coordinate must be Fourier transformed"))
    all(dim -> isodd(grid_size[dim]), fft_dims) ||
        throw(ArgumentError("grid sizes must be odd in Fourier directions"))

    # Match each bounded-data tuple entry to its inhomogeneous storage dimension and matrix size.
    inhomogeneous_dims = Tuple(dim for dim in 1:D if dim ∉ fft_dims)
    length(inhomogeneous_dims) == NI ||
        throw(ArgumentError("layout has $(length(inhomogeneous_dims)) inhomogeneous directions, but data has $NI"))
    for (index, dim) in pairs(inhomogeneous_dims)
        n = grid_size[dim]
        length(xs[index]) == length(ws[index]) == n ||
            throw(ArgumentError("points and weights for storage dimension $dim have incompatible sizes"))
        size(D₁[index]) == size(D₂[index]) == size(D₁⁺[index]) == size(D₂⁺[index]) == (n, n) ||
            throw(ArgumentError("derivative matrices for storage dimension $dim have incompatible sizes"))
    end

    # Convert only data whose scalar type differs, preserving intentionally shared objects.
    converted_xs = map(x -> _convert_grid_data(T, x), xs)
    converted_ws = map(w -> _convert_grid_data(T, w), ws)
    operators = map(values -> map(x -> _convert_grid_data(T, x), values), (D₁, D₂, D₁⁺, D₂⁺))

    # Finite union element types allow direction-specific FDGrids matrices inside explicit NTuples.
    weight_type = Union{typeof.(converted_ws)...}
    operator_types = map(values -> Union{typeof.(values)...}, operators)

    # Cache the lazy tensor product once and encode all layout metadata in the concrete grid type.
    product_weights = RectangularProductWeights{NI, T, weight_type}(converted_ws)

    return RectangularGrid{NI, T, grid_size, D, axes, fft_dims, Union{typeof.(converted_xs)...},
                           operator_types..., Union{typeof.(converted_ws)...},
                           typeof(product_weights), length(scales)}(
        converted_xs, operators..., converted_ws, product_weights, T.(scales))
end

# =============================================================================================== #
# AbstractGrid accessors                                                                         #
# =============================================================================================== #

"""
    size(g::RectangularGrid) -> Tuple

Return the physical array size of `g` in storage order.

# Arguments

- `g`: rectangular grid being queried.

# Returns

The size tuple encoded by the concrete grid type.
"""
Base.size(::RectangularGrid{NI, T, S}) where {NI, T, S} = S

"""
    eltype(::Type{<:RectangularGrid}) -> Type
    eltype(g::RectangularGrid) -> Type

Return the real scalar type used by `g`.

# Arguments

- `g`: rectangular grid instance, or its concrete type.

# Returns

The real scalar type used by its numerical data.
"""
Base.eltype(::Type{<:RectangularGrid{NI, T}}) where {NI, T} = T
Base.eltype(::RectangularGrid{NI, T}) where {NI, T} = T

"""
    weights(g::RectangularGrid) -> RectangularProductWeights

Return the cached lazy tensor-product quadrature weights of `g`.

# Arguments

- `g`: rectangular grid being queried.

# Returns

The grid's cached [`RectangularProductWeights`](@ref).
"""
weights(g::RectangularGrid) = g.w

"""
    wavenumber_scale(g::RectangularGrid, storage_dim) -> Real

Return the physical wavenumber scale for `storage_dim`. An inhomogeneous dimension returns one;
the `i`th Fourier dimension returns `g.scales[i]`. An index outside the grid's storage dimensions
throws `BoundsError`.

# Arguments

- `g`: rectangular grid being queried.
- `storage_dim`: one-based storage-dimension index.

# Returns

The real scale mapping an integer Fourier mode to physical wavenumber in that dimension.
"""
function wavenumber_scale(g::RectangularGrid{NI, T, S, D}, storage_dim::Int) where {NI, T, S, D}
    1 <= storage_dim <= D || throw(BoundsError(size(g), storage_dim))
    index = findfirst(==(storage_dim), fft_storage_dims(g))
    return isnothing(index) ? one(T) : g.scales[index]
end

"""
    points(g::RectangularGrid; dealias=false) -> Tuple

Return one broadcast-compatible coordinate array per storage dimension. Bounded coordinates come
from `g.xs`. For a Fourier direction with scale ``s``, coordinates are equispaced on ``[0,2π/s)``
without a repeated endpoint. With `dealias=true`, Fourier dimensions use NSEBase's odd 3/2-rule
padded sizes.

# Arguments

- `g`: rectangular grid whose coordinates are requested.

# Keyword arguments

- `dealias`: whether to use padded Fourier resolutions.

# Returns

A tuple of coordinate arrays whose singleton dimensions make joint broadcasting allocation-free.
"""
function points(g::RectangularGrid; dealias::Bool=false)
    storage_size = dealias ? _rectangular_padded_size(g) : size(g)
    return points(g, map(dim -> storage_size[dim], fft_storage_dims(g)))
end

"""
    points(g::RectangularGrid, homogeneous_size) -> Tuple

Return broadcast-compatible coordinates at explicit Fourier resolutions. `homogeneous_size` uses
`NSEBase.fft_storage_dims` order; bounded resolutions and coordinates do not change.

# Arguments

- `g`: rectangular grid whose coordinates are requested.
- `homogeneous_size`: positive Fourier resolutions in transform order.

# Returns

A tuple of broadcast-compatible coordinate arrays.
"""
function points(g::RectangularGrid, homogeneous_size::NTuple{N, Int}) where {N}
    all(>(0), homogeneous_size) || throw(ArgumentError("grid sizes must be positive"))
    return _rectangular_points(g, homogeneous_size, g.xs)
end

"""
    growto(g::RectangularGrid, homogeneous_size) -> RectangularGrid

Return the same grid geometry and operators with new positive odd Fourier resolutions. All bounded
data and wavenumber scales retain object identity; only the size type parameter changes.

# Arguments

- `g`: rectangular grid to resize.
- `homogeneous_size`: new positive odd Fourier resolutions in transform order.

# Returns

A grid sharing all bounded data with `g` and carrying the requested Fourier sizes.
"""
function growto(g::RectangularGrid{NI, T, S, D, AXES, FFT_DIMS, XS, D1, D2, A1, A2, WS, WP, NFFT},
                homogeneous_size::NTuple{N, Int}) where {NI, T, S, D, AXES, FFT_DIMS, XS, D1, D2, A1, A2, WS, WP, NFFT, N}
    all(>(0), homogeneous_size) || throw(ArgumentError("grid sizes must be positive"))
    all(isodd, homogeneous_size) || throw(ArgumentError("grid sizes must be odd in Fourier directions"))
    grid_size = _rectangular_storage_size(g, homogeneous_size)
    return RectangularGrid{NI, T, grid_size, D, AXES, FFT_DIMS, XS, D1, D2, A1, A2, WS, WP, NFFT}(
        g.xs, g.D₁, g.D₂, g.D₁⁺, g.D₂⁺, g.ws, g.w, g.scales)
end

"""
    convert(T, g::RectangularGrid) -> RectangularGrid

Convert points, operators, adjoints, weights, and scales to real scalar type `T`. Converting to the
existing type returns `g`; otherwise supplied adjoints are converted directly, not recomputed.

# Arguments

- `T`: target real scalar type.
- `g`: rectangular grid to convert.

# Returns

`g` itself when its scalar type is already `T`, otherwise a converted grid.
"""
function Base.convert(::Type{T}, g::RectangularGrid{NI, T0, S, D, AXES, FFT_DIMS}) where {T<:Real, NI, T0, S, D, AXES, FFT_DIMS}
    T === T0 && return g
    return RectangularGrid(g.xs, g.D₁, g.D₂, g.D₁⁺, g.D₂⁺, g.ws, g.scales, S, AXES, FFT_DIMS, T)
end

# =============================================================================================== #
# Finite-difference operator access                                                              #
# =============================================================================================== #

"""
    NSEBase.derivative_matrix(g::RectangularGrid, storage_dim, Val(order), mode)

Return the forward or discrete-adjoint FDGrids operator for an inhomogeneous storage dimension.
`order` must be one or two. The `Forward()` and `AdjointDiscrete()` tags keep operator lookup type
stable even though FDGrids forward and adjoint matrices have different concrete types.

# Arguments

- `g`: rectangular grid owning the operators.
- `storage_dim`: bounded storage dimension to differentiate.
- `Val(order)`: derivative order, either `Val(1)` or `Val(2)`.
- `mode`: `Forward()` or `AdjointDiscrete()`.

# Returns

The stored FDGrids differentiation matrix selected entirely by dispatch and compile-time order.
"""
function NSEBase.derivative_matrix(g::RectangularGrid, storage_dim::Integer,
                                   ::Val{ORDER}, ::Forward) where {ORDER}
    index = findfirst(==(storage_dim), inhomogeneous_storage_dims(g))
    isnothing(index) && throw(ArgumentError("storage dimension $storage_dim is not inhomogeneous"))
    ORDER == 1 && return g.D₁[index]
    ORDER == 2 && return g.D₂[index]
    throw(ArgumentError("only first- and second-derivative operators are available"))
end

function NSEBase.derivative_matrix(g::RectangularGrid, storage_dim::Integer,
                                   ::Val{ORDER}, ::AdjointDiscrete) where {ORDER}
    index = findfirst(==(storage_dim), inhomogeneous_storage_dims(g))
    isnothing(index) && throw(ArgumentError("storage dimension $storage_dim is not inhomogeneous"))
    ORDER == 1 && return g.D₁⁺[index]
    ORDER == 2 && return g.D₂⁺[index]
    throw(ArgumentError("only first- and second-derivative operators are available"))
end
