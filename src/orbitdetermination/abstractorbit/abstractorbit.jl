@doc raw"""
    AbstractOrbit{T <: Real, U <: Number}

Supertype for the orbits API.
"""
abstract type AbstractOrbit{T <: Real, U <: Number} end

numtypes(::AbstractOrbit{T, U}) where {T, U} = T, U

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
snr(orbit::AbstractOrbit{T, U}) where {T, U} = iszero(orbit) ? Vector{U}(undef, 0) :
    abs.(orbit()) ./ sigmas(orbit)

# Update `orbit` iff `_orbit_` is complete and has a lower nrms
function updateorbit(orbit::AbstractOrbit{T, T}, _orbit_::AbstractOrbit{T, T},
    radec::Vector{RadecMPC{T}}) where {T <: Real}
    L1, L2, L = nobs(orbit), nobs(_orbit_), length(radec)
    # Both orbits are complete
    if L1 == L2 == L
        return nrms(orbit) <= nrms(_orbit_) ? orbit : _orbit_
    # Only orbit is complete
    elseif L1 == L
        return orbit
    # Only _orbit_ is complete
    elseif L2 == L
        return _orbit_
    # Neither orbit is complete
    else
        return orbit
    end
end

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
    keplerian(orbit, params) where {T <: Real}

Return the heliocentric ecliptic keplerian osculating elements of an orbit
at its reference epoch, as well as the projected covariance matrix.

See also [`pv2kep`](@ref).

## Arguments

- `orbit::AbstractOrbit{T, T}`: a priori orbit.
- `params::Parameters{T}`: see [`Parameters`](@ref).

!!! warning
    This function may change the (global) `TaylorSeries` variables.
"""
function keplerian(orbit::AbstractOrbit{T, T}, params::Parameters{T}) where {T <: Real}
    # Set jet transport variables
    set_od_order(T, 2)
    # Reference epoch [Julian days TDB]
    t = epoch(orbit)
    jd0 = t + PE.J2000
    # Jet transport initial condition
    q0 = orbit(t) + diag(orbit.J) .* get_variables(T, 2)
    # Origin
    x0 = zeros(T, 6)
    # Osculating keplerian elements
    osc = pv2kep(q0 - params.eph_su(t); jd = jd0, frame = :ecliptic)
    osc0 = [osc.a, osc.e, osc.i, osc.Ω, osc.ω, mod(osc.M, 360)]
    osc00 = constant_term.(osc0)
    # Projected covariance matrix
    t_car2kep = TS.jacobian(osc0, x0)
    Γ_kep = t_car2kep * covariance(orbit) * t_car2kep'

    return osc00, Γ_kep
end

@doc raw"""
    uncertaintyparameter(orbit, params) where {T <: Real}

Return the Minor Planet Center uncertainty parameter.

## Arguments

- `orbit::AbstractOrbit{T, T}`: a priori orbit.
- `params::Parameters{T}`: see [`Parameters`](@ref).

!!! warning
    This function may change the (global) `TaylorSeries` variables.

!!! reference
    https://www.minorplanetcenter.net/iau/info/UValue.html.
"""
function uncertaintyparameter(orbit::AbstractOrbit{T, T},
    params::Parameters{T}) where {T <: Real}
    # Set jet transport variables
    set_od_order(T, 2)
    # Reference epoch [Julian days TDB]
    t = epoch(orbit)
    jd0 = t + PE.J2000
    # Jet transport initial condition
    q0 = orbit(t) + diag(orbit.J) .* get_variables(T, 2)
    # Origin
    x0 = zeros(T, 6)
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
    Γ_kep = t_car2kep * covariance(orbit) * t_car2kep'
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

function summary(orbit::AbstractOrbit)
    O = nameof(typeof(orbit))
    T, U = numtypes(orbit)
    D = orbit.dynamics
    Nobs, Nout = nobs(orbit), nout(orbit.res)
    Ndays = @sprintf("%.8f", numberofdays(orbit))
    t0 = epoch(orbit) + PE.J2000
    d0 = julian2datetime(t0)
    Q = nrms(orbit)
    q0, σ0 = orbit(), sigmas(orbit)
    sq0 = [rpad(@sprintf("%+.12E", q0[i]), 25) for i in eachindex(q0)]
    sσ0 = [rpad(@sprintf("%+.12E", σ0[i]), 25) for i in eachindex(σ0)]
    s = string(
        "$O with numeric types ($T, $U)\n",
        repeat('-', 68), "\n",
        "Dynamical model: $D\n",
        "Astrometry: $Nobs observations ($Nout outliers) spanning $Ndays days\n",
        "Epoch: $t0 JDTDB ($d0 TDB)\n",
        "NRMS: $Q\n",
        repeat('-', 68), "\n",
        "Variable    Nominal value            Uncertainty              Units\n",
        "x           ", sq0[1], sσ0[1], "au\n",
        "y           ", sq0[2], sσ0[2], "au\n",
        "z           ", sq0[3], sσ0[3], "au\n",
        "vx          ", sq0[4], sσ0[4], "au/day\n",
        "vy          ", sq0[5], sσ0[5], "au/day\n",
        "vz          ", sq0[6], sσ0[6], "au/day\n"
    )
    return s
end