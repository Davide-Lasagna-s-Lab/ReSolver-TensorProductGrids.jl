# `exp(cos(z))` has a nonzero coefficient at every Fourier wavenumber, so the homogeneous
# derivative and norm tests exercise every stored mode rather than a selected band of modes.
periodic_profile(z) = exp(cos(z))
periodic_profile_d1(z) = -sin(z) * exp(cos(z))
periodic_profile_d2(z) = (sin(z)^2 - cos(z)) * exp(cos(z))

# The exact mean square over one period is the modified Bessel function I₀(2). Its defining series
# is evaluated to Float64 precision here to avoid adding SpecialFunctions as a test dependency:
#
#     (1 / 2π) ∫₀²π exp(2cos(z)) dz = I₀(2) = Σₖ₌₀∞ 1 / (k!)².
const PERIODIC_PROFILE_NORM2 = let
    term = 1.0
    value = term
    for k in 1:20
        term /= k^2
        value += term
    end
    value
end

# A common polynomial profile for every bounded direction. Parameterising the interval preserves
# the natural channel domain `[-1, 1]` and cavity/duct domain `[0, 1]` without changing the tests'
# tensor-product structure. Its first two derivatives and squared norm are analytical.
bounded_profile(x, lim) = (x - lim[1]) * (lim[2] - x)
bounded_profile_d1(x, lim) = lim[1] + lim[2] - 2x
bounded_profile_d2(x, lim) = -2one(x)
bounded_profile_norm2(lim) = (lim[2] - lim[1])^5 / 30

# A distinct bounded factor keeps the discrete-adjoint identities nontrivial in every case.
function dual_bounded_profile(x, lim)
    ξ = (x - lim[1]) / (lim[2] - lim[1])
    return 1 + ξ + ξ^2
end

"""Embed values over the bounded directions in the homogeneous zero mode of an `FTField`."""
function ftfield_from_inhomogeneous(g, values)
    u = FTField(g)
    indices = ntuple(dim -> dim in fft_storage_dims(g) ? 1 : Colon(), length(size(g)))
    view(parent(u), indices...) .= values
    return u
end
