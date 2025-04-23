@doc raw"""
    AbstractOrbit{T <: Real, U <: Number}

Supertype for the orbits API.
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

# Timespan of the observation window in days
function numberofdays(orbit::AbstractOrbit)
    t0, tf = minmaxdates(orbit)
    return (tf - t0).value / 86_400_000
end

# Assemble the vector of optical astrometry
astrometry(orbit::AbstractOrbit) = astrometry(orbit.tracklets)

# Evaluation in time method
(orbit::AbstractOrbit)(t = epoch(orbit)) = t <= epoch(orbit) ? orbit.bwd(t) : orbit.fwd(t)

# Initialize a vector of residuals consistent with orbit
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

@doc raw"""
    sigmas(::AbstractOrbit)

Return the uncertainties in barycentric cartesian coordinates at the reference epoch.
"""
sigmas(orbit::AbstractOrbit) = map(x -> x < 0 ? NaN * x : sqrt(x), variances(orbit))

@doc raw"""
    snr(::AbstractOrbit)

Return the signal-to-noise ratios in barycentric cartesian
coordinates at the reference epoch.
"""
snr(orbit::AbstractOrbit) = abs.(orbit()) ./ sigmas(orbit)

@doc raw"""
    jplcompare(des::String, orbit::AbstractOrbit)

Return the absolute difference between `orbit` and JPL's orbit for object
`des` computed in barycentric cartesian coordinates at `epoch(orbit)` in
units of `sigmas(orbit)`.
"""
function jplcompare(des::String, orbit::AbstractOrbit)
    # Load Solar System ephemerides
    loadjpleph()
    # orbit's initial condition
    q1 = cte.(orbit())
    # Time of first (last) observation
    t0, tf = minmaxdates(orbit)
    # Download JPL ephemerides
    bsp = smb_spk("DES = $(des);", t0 - Minute(10), tf + Minute(10))
    # Object ID
    id = parse(Int, bsp[1:end-4])
    # Load JPL ephemerides
    furnsh(bsp)
    # Delete ephemerides file
    rm(bsp)
    # JPL barycentric state vector
    q2 = kmsec2auday(spkgeo(id, julian2etsecs(epoch(orbit) + PE.J2000),
        "J2000", 0)[1])
    # Absolute difference in sigma units
    return @. abs(q1 - q2) / sigmas(orbit)
end

@doc raw"""
    uncertaintyparameter(od, orbit, params) where {D, T <: Real}

Return the Minor Planet Center uncertainty parameter.

## Arguments

- `od::ODProblem{D, T}`: an orbit determination problem.
- `orbit::AbstractOrbit{T, T}`: a priori orbit.
- `params::Parameters{T}`: see [`Parameters`](@ref).

!!! reference
    https://www.minorplanetcenter.net/iau/info/UValue.html

!!! warning
    This function will change the (global) `TaylorSeries` variables.
"""
function uncertaintyparameter(od::ODProblem{D, T}, orbit::AbstractOrbit{T, T},
    params::Parameters{T}) where {D, T <: Real}
    # Check consistency between od and orbit
    @assert od.tracklets == orbit.tracklets
    # Reference epoch [Julian days TDB]
    t = epoch(orbit)
    jd0 = t + PE.J2000
    # Barycentric initial conditions
    q00 = orbit(t)
    # Scaling factors
    scalings = all(>(0), variances(orbit)) ? sigmas(orbit) : @. abs(q00) / 1e6
    # Jet transport variables
    dq = scaled_variables("dx", scalings; order = 2)
    # Origin
    x0 = zeros(T, 6)
    # Initial conditions
    q0 = q00 + dq
    # Propagation and residuals
    res = init_residuals(TaylorN{T}, orbit)
    propres!(res, od, jd0, q0, params)
    nobs = 2 * notout(res)
    # Covariance matrix
    Q = nms(res)
    C = (nobs/2) * TS.hessian(Q, x0)
    Γ = inv(C)
    # Osculating keplerian elements
    osc = pv2kep(q0 - params.eph_su(t); jd = jd0, frame = :ecliptic)
    # Semimajor axis [au], eccentricity and time of perihelion passage [julian days]
    @unpack a, e, tp = osc
    # Gauss gravitational constant [deg]
    k_0 = 180 * k_gauss / π
    # Orbital period [days]
    P = 2π * sqrt(a^3 / μ_S)
    # Projected covariance matrix
    t_car2kep = TS.jacobian([tp, P], x0)
    Γ_kep = t_car2kep * Γ * t_car2kep'
    # Uncertainties
    dtp, dP = sqrt.(diag(Γ_kep))
    # Convert orbital period to years
    P = P / yr
    # In-orbit longitude runoff [arcsec / decade]
    runoff = (dtp * e(x0) + 10 * dP / P(x0)) * (k_0 / P(x0)) * 3_600 * 3
    # Uncertainty parameter
    C = log(648_000) / 9
    U = floor(Int, log(runoff) / C) + 1

    return clamp(U, 0, 10)
end