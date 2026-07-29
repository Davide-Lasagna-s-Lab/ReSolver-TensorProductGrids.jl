# =============================================================================================== #
# Shared case-constructor validation                                                             #
# =============================================================================================== #

"""Validate a Cartesian base-flow tuple and every non-`nothing` component's bounded shape."""
function _validate_base_flow(g::AbstractGrid, base_flow, ::Val{N}, case_name) where {N}
    base_flow isa Tuple{Vararg{Any, N}} ||
        throw(ArgumentError("a $case_name base flow must be an $N-component tuple"))

    expected_size = map(dim -> size(g, dim), inhomogeneous_storage_dims(g))
    for (i, component) in pairs(base_flow)
        isnothing(component) && continue
        component isa AbstractArray ||
            throw(ArgumentError("base-flow component $i must be an array or nothing"))
        size(component) == expected_size || throw(DimensionMismatch(
            "base-flow component $i has size $(size(component)); expected $expected_size"))
    end
    return base_flow
end

# =============================================================================================== #
# Bounded-domain geometry and canonical cavity liftings                                          #
# =============================================================================================== #

"""Return the underlying global grid used to define a wrapped grid's physical geometry."""
function _global_case_grid(g::AbstractGrid)
    current = g
    while applicable(Base.parent, current)
        next = parent(current)
        next === current && break
        current = next
    end
    return current
end

"""Return the global physical interval of one bounded storage direction."""
function _bounded_interval(g::AbstractGrid, dim::Int)
    x = vec(points(_global_case_grid(g))[dim])
    return extrema(x)
end

"""Require a case's bounded directions to have a common positive physical length."""
function _validate_equal_bounded_lengths(g::AbstractGrid, ::Val{N}, case_name) where {N}
    dims = inhomogeneous_storage_dims(g)
    length(dims) == N || throw(ArgumentError("a $case_name requires $N bounded directions"))

    lengths = map(dims) do dim
        xmin, xmax = _bounded_interval(g, dim)
        xmax - xmin
    end
    all(length -> length > zero(length), lengths) ||
        throw(ArgumentError("a $case_name requires bounded directions of positive length"))
    all(length -> isapprox(length, first(lengths)), lengths) ||
        throw(ArgumentError("a $case_name requires equal side lengths; got $lengths"))
    return nothing
end

"""Return one bounded coordinate normalized from its physical interval to `[0,1]`."""
function _normalized_bounded_coordinate(g::AbstractGrid, dim::Int)
    x = vec(points(g)[dim])
    xmin, xmax = _bounded_interval(g, dim)
    return (x .- xmin) ./ (xmax - xmin)
end

"""Return the divergence-free `(U,V)` polynomial shared by the 2D and 3D cavity liftings."""
function _lid_driven_cavity_xy_base(g::AbstractGrid)
    xdim, ydim = inhomogeneous_storage_dims(g)[1:2]
    ξ = _normalized_bounded_coordinate(g, xdim)
    η = _normalized_bounded_coordinate(g, ydim)

    F = @. 16ξ^2 * (1 - ξ)^2
    dF = @. 32ξ * (1 - ξ) * (1 - 2ξ)
    H = @. η^2 * (η - 1)
    dH = @. 3η^2 - 2η
    return F * transpose(dH), -dF * transpose(H)
end

# =============================================================================================== #
# Reusable forcing selection                                                                     #
# =============================================================================================== #

"""Return no forcing for a zero amplitude and a constant body force otherwise."""
_constant_body_force(value::Real, i::Int) = iszero(value) ? NoForce() : ConstantBodyForce(value; i)
