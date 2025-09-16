"""
    AbstractOrbit{D, T <: Real, U <: Number}

Supertype for the orbits interface.

Every orbit has:
- a dynamical function of type `D`,
- a vector of variables of type `Vector{Int}`,
- a vector of optical astrometry of type `O <: AbstractOpticalVector{T}`,
- a vector of optical tracklets of type `TrackletVector{T}`,
- a backward and a forward integration, both of type `DensePropagation2{T, U}`,
- a vector of optical residuals of type `Vector{OpticalResidual{T, U}}`,
- a jacobian representing the transformation from the space of residuals to
    barycentric coordinates, of type `Matrix{T}`.
"""
abstract type AbstractOrbit{D, T <: Real, U <: Number} end

# Numeric types
numtypes(::AbstractOrbit{D, T, U}) where {D, T, U} = T, U

# Degrees of freedom
dof(x::AbstractOrbit) = dof(Val(x.dynamics))

# Variables
variables(x::AbstractOrbit) = x.variables
numvars(x::AbstractOrbit) = length(x.variables)

"""
    epoch(::AbstractOrbit)

Return the reference epoch of an orbit in TDB days since J2000.
"""
epoch(x::AbstractOrbit) = x.bwd.t0

# Number of observations
noptical(x::AbstractOrbit) = length(x.optical)
nradar(x::AbstractOrbit) = hasradar(x) ? length(x.radar) : 0
nobs(x::AbstractOrbit) = noptical(x) + nradar(x)

function notoutobs(x::AbstractOrbit)
    N = notoutobs(x.ores)
    if hasradar(x)
        N += notoutobs(x.rres)
    end
    return N
end

"""
    minmaxdates(::AbstractOrbit)

Return the dates of the earliest and latest observation of an orbit.
"""
function minmaxdates(x::AbstractOrbit)
    t0, tf = minmaxdates(x.optical)
    if hasradar(x)
        _t0_, _tf_ = minmaxdates(x.radar)
        t0, tf = min(t0, _t0_), max(tf, _tf_)
    end
    return t0, tf
end

# Timespan of the observation window in days
function numberofdays(x::AbstractOrbit)
    t0, tf = minmaxdates(x)
    return (tf - t0).value / 86_400_000
end

# Return the optical astrometry in an orbit
optical(x::AbstractOrbit) = x.optical
radar(x::AbstractOrbit) = hasradar(x) ? x.radar : nothing

# Evaluation in time method
(x::AbstractOrbit)(t = epoch(x)) = t <= epoch(x) ? x.bwd(t) : x.fwd(t)

# Initial condition and jet transport perturbation
function initialcondition(q00::AbstractVector{T}, variables::AbstractVector{Int},
                          dof::Int, params::Parameters{T}) where {T <: Real}
    if length(q00) == dof
        return q00
    elseif length(q00) < dof
        q0 = zeros(T, dof)
        q0[variables] .= q00
        for i in 7:9
            q0[i] = iszero(q0[i]) ? params.marsden_coeffs[i-6] : q0[i]
        end
        return q0
    else
        throw(ArgumentError("The number of degrees of freedom `dof` cannot be smaller \
            than the length of the preliminary initial condition `q00`"))
    end
end

initialcondition(orbit::AbstractOrbit, dof::Int, params::Parameters) =
    initialcondition(orbit(), variables(orbit), dof, params)

function jtperturbation(sigmas::AbstractVector{T}, variables::AbstractVector{Int},
                        dof::Int, order::Int, params::Parameters{T}) where {T <: Real}
    scalings = zeros(T, dof)
    jtvariables = get_variables(T, order)
    dq = [zero(jtvariables[1]) for _ in 1:dof]
    dq[variables] .= jtvariables
    for i in 1:6
        scalings[i] = isnan(sigmas[i]) ? 1e-8 : sigmas[i]
    end
    for i in 7:dof
        scalings[i] = params.marsden_scalings[i-6]
    end
    dq .*= scalings
    return dq
end

jtperturbation(orbit::AbstractOrbit, variables::AbstractVector{Int}, dof::Int, order::Int,
    params::Parameters) = jtperturbation(sigmas(orbit), variables, dof, order, params)

# Number of outlier / non-outlier residuals
function nout(x::AbstractOrbit)
    y = nout(x.ores)
    if hasradar(x)
        y += nout(x.rres)
    end
    return y
end

function notout(x::AbstractOrbit)
    y = notout(x.ores)
    if hasradar(x)
        y += notout(x.rres)
    end
    return y
end

# Initialize a vector of optical residuals consistent with orbit
function init_optical_residuals(::Type{V}, orbit::AbstractOrbit{D, T, U}) where {D,
                                T <: Real, U <: Number, V <: Number}
    # Initialize vector of optical residuals
    res = Vector{OpticalResidual{T, V}}(undef, noptical(orbit))
    for i in eachindex(res)
        ra, dec = zero(V), zero(V)
        @unpack wra, wdec, dra, ddec, corr, outlier = orbit.res[i]
        res[i] = OpticalResidual{T, V}(ra, dec, wra, wdec, dra, ddec, corr, outlier)
    end

    return res
end

init_residuals(::Type{U}, od::OpticalODProblem, orbit::AbstractOrbit) where {U <: Number} =
    init_optical_residuals(U, od, orbit)

init_residuals(::Type{U}, od::MixedODProblem, orbit::AbstractOrbit) where {U <: Number} =
    init_optical_residuals(U, od, orbit), init_radar_residuals(U, od, orbit)

# Target functions
function chi2(x::AbstractOrbit{D, T, U}) where {D, T <: Real, U <: Number}
    iszero(x) && return Inf * one(U)
    y = chi2(x.ores)
    if hasradar(x)
        y += chi2(x.rres)
    end
    return y
end

nms(x::AbstractOrbit) = chi2(x) / notoutobs(x)
nrms(x::AbstractOrbit) = sqrt(nms(x))

"""
    critical_value(::AbstractOrbit)

Return the chi-square critical value corresponding to the normlized mean square
residual of an orbit.
"""
function critical_value(x::AbstractOrbit{D, T, U}) where {D, T <: Real, U <: Number}
    # Unsuccessful orbit
    iszero(x) && return one(T)
    # Chi square distribution with N degrees of freedom
    N = 2 * notout(x)
    χ2 = N * nms(x)
    d = Chisq(N)
    # Evaluate cumulative probability at χ2
    return cdf(d, χ2)
end

"""
    sigmas(::AbstractOrbit)

Return the uncertainties in barycentric cartesian coordinates at the reference epoch.
"""
function sigmas(orbit::AbstractOrbit{D, T, U}) where {D, T, U}
    iszero(orbit) && return Vector{T}(undef, 0)
    v = variances(orbit)
    y = NaN * one(T)
    for (i, x) in enumerate(v)
        v[i] = x < 0 ? y : sqrt(x)
    end
    return v
end

"""
    snr(::AbstractOrbit)

Return the signal-to-noise ratios in barycentric cartesian
coordinates at the reference epoch.
"""
function snr(orbit::AbstractOrbit{D, T, U}) where {D, T, U}
    iszero(orbit) && return Vector{U}(undef, 0)
    q0, ss = orbit(), sigmas(orbit)
    for i in eachindex(q0)
        q0[i] = abs(q0[i]) / ss[i]
    end
    return q0
end

# Update `x` iff `y` is complete and has a lower nrms
function updateorbit(x::AbstractOrbit{D, T, T}, y::AbstractOrbit{D, T, T},
                     optical::AbstractOpticalVector{T}) where {D, T <: Real}
    L1, L2, L = noptical(x), noptical(y), length(optical)
    # Both orbits are complete
    if L1 == L2 == L
        return nrms(x) <= nrms(y) ? x : y
    # Only x is complete
    elseif L1 == L
        return x
    # Only y is complete
    elseif L2 == L
        return y
    # Neither orbit is complete
    else
        return x
    end
end

#=
"""
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
=#

"""
    osculating(orbit, params)

Return the heliocentric ecliptic osculating elements of an `orbit`
at its reference epoch.

See also [`cartesian2osculating`](@ref).

!!! warning
    This function may change the (global) `TaylorSeries` variables.
"""
function osculating(orbit::AbstractOrbit{D, T, T},
                    params::Parameters{T}) where {D, T <: Real}
    # Set jet transport variables
    Npar = numvars(orbit)
    set_od_order(T, 2, Npar)
    # Reference epoch [MJD TDB]
    t = epoch(orbit)
    mjd0 = t + MJD2000
    # Jet transport initial condition
    q0 = orbit(t) + diag(orbit.jacobian) .* get_variables(T, 2)
    # Origin
    x0 = zeros(T, Npar)
    # Osculating orbital elements
    osc = cartesian2osculating(q0[1:6] - params.eph_su(t), mjd0; μ = μ_S,
                               frame = :ecliptic, Γ_car = covariance(orbit))

    return evaldeltas(osc, x0)
end

"""
    uncertaintyparameter(orbit, params)

Return the Minor Planet Center uncertainty parameter of an `orbit`.

!!! warning
    This function may change the (global) `TaylorSeries` variables.

!!! reference
    See:
    - https://www.minorplanetcenter.net/iau/info/UValue.html
"""
function uncertaintyparameter(orbit::AbstractOrbit{D, T, T},
                              params::Parameters{T}) where {D, T <: Real}
    # Set jet transport variables
    Npar = numvars(orbit)
    set_od_order(T, 2, Npar)
    # Reference epoch [MJD TDB]
    t = epoch(orbit)
    mjd0 = t + MJD2000
    # Jet transport initial condition
    q0 = orbit(t) + diag(orbit.jacobian) .* get_variables(T, 2)
    # Origin
    x0 = zeros(T, Npar)
    # Osculating keplerian elements
    osc = cartesian2osculating(q0 - params.eph_su(t), mjd0; μ = μ_S, frame = :ecliptic,
                               Γ_car = covariance(orbit))
    # Uncertainty parameter is not defined for hyperbolic orbits
    ishyperbolic(osc) && throw(ArgumentError("Uncertainty parameter is not defined for \
        hyperbolic orbits"))
    # Semimajor axis [au], eccentricity and time of perihelion passage [MJD TDB]
    a = semimajoraxis(osc)
    e = eccentricity(osc)
    tp::TaylorN{T} = timeperipass(osc)
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

# Observed magnitude model
# h = H + 5 * log10(d_BS * d_BO) - 2.5 * log10(qs)
hmodel(t, p) = @. p[1] + t

"""
    absolutemagnitude(::AbstractOrbit, ::Parameters)

Return the absolute magnitude of an orbit, as well as its standard error.
This function uses the linear H and G model for asteroids.

!!! reference
    See:
    - https://minorplanetcenter.net/iau/ECS/MPCArchive/1985/MPC_19851227.pdf
"""
function absolutemagnitude(orbit::AbstractOrbit, params::Parameters)
    # Extract optical astrometry
    optical = orbit.optical
    # Observation times
    ts = @. dtutc2days(date(optical))
    # Observed magnitudes and bands
    hs, bands = @. mag(optical), band.(optical)
    # Convert magnitudes to V band
    hs .+= getindex.(Ref(V_BAND_CORRECTION), bands)
    # Asteroid, Earth and Sun positions
    xa = [orbit(t)[1:3] for t in ts]
    xe = [params.eph_ea(t)[1:3] for t in ts]
    xs = [params.eph_su(t)[1:3] for t in ts]
    # Distances
    d_BS = @. norm(xa - xs)     # Asteroid-Sun
    d_OS = @. norm(xe - xs)     # Earth-Sun
    d_BO = @. norm(xa - xe)     # Asteroid-Earth
    # Phase angles
    αs = @. acos((d_BO^2 + d_BS^2 - d_OS^2) / (2 * d_BO * d_BS))
    # Phase integrals
    ϕ1s = @. exp(-PHASE_INTEGRAL_A1 * tan(αs/2)^PHASE_INTEGRAL_B1)
    ϕ2s = @. exp(-PHASE_INTEGRAL_A2 * tan(αs/2)^PHASE_INTEGRAL_B2)
    qs = (1 - SLOPE_PARAMETER) * ϕ1s + SLOPE_PARAMETER * ϕ2s
    # Absolute magnitude
    tdata = @. 5 * log10(d_BS * d_BO) - 2.5 * log10(qs)
    mask = @. !isnan(tdata) && !isnan(hs)
    H0 = [mean(view(hs, mask))]
    fit = curve_fit(hmodel, view(tdata, mask), view(hs, mask), H0)
    H = fit.param[1]
    dH = stderror(fit)[1]

    return H, dH
end

function summary(orbit::AbstractOrbit)
    O = nameof(typeof(orbit))
    T, U = numtypes(orbit)
    D = orbit.dynamics
    Nobs, Nout = nobs(orbit), nout(orbit)
    Ndays = @sprintf("%.8f", numberofdays(orbit))
    t0 = epoch(orbit) + PE.J2000
    d0 = julian2datetime(t0)
    Q = nrms(orbit)
    q0, σ0 = orbit(), sigmas(orbit)
    sq0 = [rpad(@sprintf("%+.12E", q0[i]), 25) for i in eachindex(q0)]
    sσ0 = [rpad(@sprintf("%+.12E", σ0[i]), 25) for i in eachindex(σ0)]
    names = ["x", "y", "z", "vx", "vy", "vz", "A2", "A1", "A3"]
    units = ["au", "au", "au", "au/day", "au/day", "au/day", "au/day²",
             "au/day²", "au/day²"]
    s1 = string(
        "$O{$T, $U}\n",
        repeat('-', 69), "\n",
        "Dynamical model: $D\n",
        "Astrometry: $Nobs observations ($Nout outliers) spanning $Ndays days\n",
        "Epoch: $t0 JDTDB ($d0 TDB)\n",
        "NRMS: $Q\n",
        repeat('-', 69), "\n",
        "Variable    Nominal value            Uncertainty              Units\n",
    )
    s2 = Vector{String}(undef, numvars(orbit))
    for (i, k) in enumerate(variables(orbit))
        s2[i] = string(rpad(names[k], 12), sq0[i], sσ0[i], units[k], "\n")
    end
    return string(s1, join(s2))
end