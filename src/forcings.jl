# =============================================================================================== #
# Reusable case forcings                                                                         #
# =============================================================================================== #
#
# A grid describes geometry; a case constructor chooses physical source terms. The small callable
# policies below can be passed directly to `NSEBase.construct_equations` and remain concretely typed
# inside the nonlinear and linearised operators.

# =============================================================================================== #
# Coriolis coupling                                                                               #
# =============================================================================================== #

@doc raw"""
    CoriolisForce(Ro::Real)

Construct the skew-symmetric Coriolis force coupling streamwise velocity `u` and wall-normal
velocity `v`. For a forward evaluation,

```math
\boldsymbol{f}_{Ro} = Ro\,(v, -u, 0)
```

for a three-component state `(u, v, w)`; the trailing zero is absent for a two-component state.
Continuous- and discrete-adjoint evaluations apply the transpose, `-\boldsymbol{f}_{Ro}`. The
object stores only the dimensionless rotation coefficient `Ro`; the coupled components are always
one and two.

# Arguments

- `Ro`: signed dimensionless rotation coefficient.

# Returns

A lightweight callable force containing only `Ro`.
"""
struct CoriolisForce{T<:Real}
    Ro::T
end

function (force::CoriolisForce)(out::VectorField{N}, u::VectorField{N}, ::Forward) where {N}
    @. out[1] += force.Ro * u[2]
    @. out[2] -= force.Ro * u[1]
    return out
end

function (force::CoriolisForce)(out::VectorField{N}, u::VectorField{N},
                                ::Union{AdjointContinuous, AdjointDiscrete}) where {N}
    @. out[1] -= force.Ro * u[2]
    @. out[2] += force.Ro * u[1]
    return out
end

# =============================================================================================== #
# Constant body forcing                                                                          #
# =============================================================================================== #

"""
    ConstantBodyForce(value=1; i=1)

Construct a spatially constant body force applied to one velocity component. The force is added
only at the zero wavenumber of every Fourier direction; the selected slice still contains every
finite-difference point, so the value is uniform throughout the bounded coordinates.

`i` is the one-based index of the component in the `VectorField` passed to the force. For a channel
state ordered as `(u, v, w)`, streamwise forcing uses `i=1`. A square duct is also stored as
`(u, v, w)`, but its physical streamwise direction is `z`, so [`SquareDuctFlow`](@ref) uses `i=3`.
The state-independent contribution is applied only for `Forward()` evaluation of the nonlinear
equations. Its linearisation is zero, so continuous- and discrete-adjoint calls leave `out`
unchanged.

# Arguments

- `value`: constant force amplitude.

# Keyword arguments

- `i`: one-based index of the forced velocity component.

# Returns

A callable constant forcing policy.
"""
struct ConstantBodyForce{T}
    value::T
    i::Int
end

ConstantBodyForce(value=1; i=1) = ConstantBodyForce(value, i)

function (force::ConstantBodyForce)(out::VectorField{N, <:FTField}, _, ::Forward) where {N}
    # Validate the requested component against the field on which the force is applied.
    1 <= force.i <= N || throw(ArgumentError("force component must lie in 1:$N"))

    # A spatially uniform force occupies the zero mode of every Fourier direction.
    component = out[force.i]
    mean_mode = WaveNumberVector(map(_ -> 0, fft_storage_dims(grid(component))))

    # Leave every nonzero Fourier mode untouched and add the force at all bounded points.
    mean = component[mean_mode]
    mean .+= force.value
    return out
end

# A state-independent force has zero linearisation and therefore zero adjoint action.
(::ConstantBodyForce)(out::VectorField{N, <:FTField}, _,
                      ::Union{AdjointContinuous, AdjointDiscrete}) where {N} = out
