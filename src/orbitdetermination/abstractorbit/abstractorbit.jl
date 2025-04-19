@doc raw"""
    AbstractOrbit{T <: Real, U <: Number}

Supertype for the asteroid orbits API.
"""
abstract type AbstractOrbit{T <: Real, U <: Number} end

@doc raw"""
    epoch(::AbstractOrbit)

Return the reference epoch of an orbit in TDB days since J2000.
"""
epoch(orbit::AbstractOrbit) = orbit.bwd.t0

# Number of observations
nobs(orbit::AbstractOrbit) = length(orbit.res)

@doc raw"""
    minmaxdates(::AbstractOrbit)

Return the dates of the earliest and latest observation of an orbit.
"""
function minmaxdates(orbit::AbstractOrbit)
    dates = map(t -> extrema(date, t.radec), orbit.tracklets)
    t = days2dtutc(epoch(orbit))
    return min(t, minimum(first, dates)), max(t, maximum(last, dates))
end

function numberofdays(orbit::AbstractOrbit)
    t0, tf = minmaxdates(orbit)
    return (tf - t0).value / 86_400_000
end

# Evaluation in time method
(orbit::AbstractOrbit)(t = epoch(orbit)) = t <= epoch(orbit) ? orbit.bwd(t) : orbit.fwd(t)

# Target functions
chi2(orbit::AbstractOrbit{T, U}) where {T <: Real, U <: Number} =
    iszero(orbit) ? T(Inf) : chi2(orbit.res)
nms(orbit::AbstractOrbit{T, U}) where {T <: Real, U <: Number} =
    iszero(orbit) ? T(Inf) : nms(orbit.res)
nrms(orbit::AbstractOrbit{T, U}) where {T <: Real, U <: Number} =
    iszero(orbit) ? T(Inf) : nrms(orbit.res)

@doc raw"""
    critical_value(orbit::AbstractOrbit{T, U}) where {T <: Real, U <: Number}

Return the chi-square critical value corresponding to `nms(orbit)`.
"""
function critical_value(orbit::AbstractOrbit{T, U}) where {T <: Real, U <: Number}
    # Unsuccessful orbit
    iszero(orbit) && return one(T)
    # Chi square distribution with N degrees of freedom
    N = 2 * notout(orbit.res)
    χ2 = N * nms(orbit)
    d = Chisq(N)
    # Evaluate cumulative probability at χ2
    return cdf(d, χ2)
end

function init_residuals(::Type{V}, orbit::AbstractOrbit{T, U}) where {T <: Real,
    U <: Number, V <: Number}
    # Initialize vector of residuals
    res = Vector{OpticalResidual{T, V}}(undef, nobs(orbit))
    for i in eachindex(res)
        ξ_α, ξ_δ = zero(V), zero(V)
        @unpack w_α, w_δ, μ_α, μ_δ, outlier = orbit.res[i]
        res[i] = OpticalResidual{T, V}(ξ_α, ξ_δ, w_α, w_δ, μ_α, μ_δ, outlier)
    end

    return res
end