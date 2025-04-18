@doc raw"""
    LeastSquaresOrbit{T, U} <: AbstractOrbit{T, U}

An asteroid least squares orbit.

## Fields

- `tracklets::Vector{Tracklet{T}}`: vector of tracklets.
- `bwd/fwd::TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}}`:
    backward (forward) integration.
- `res::Vector{OpticalResidual{T, U}}`: vector of optical residuals.
- `fit::LeastSquaresFit{T}`: least squares fit.
- `jacobian::Matrix{T}`: residuals space to barycentric coordinates jacobian.
"""
@auto_hash_equals struct LeastSquaresOrbit{T, U} <: AbstractOrbit{T, U}
    tracklets::Vector{Tracklet{T}}
    bwd::TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}}
    fwd::TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}}
    res::Vector{OpticalResidual{T, U}}
    fit::LeastSquaresFit{T}
    jacobian::Matrix{T}
    # Inner constructor
    function LeastSquaresOrbit{T, U}(tracklets::Vector{Tracklet{T}},
        bwd::TaylorInterpolant{T, U, 2}, fwd::TaylorInterpolant{T, U, 2},
        res::Vector{OpticalResidual{T, U}}, fit::LeastSquaresFit{T},
        jacobian::Matrix{T}) where {T <: Real, U <: Number}
        @assert bwd.t0 == fwd.t0 "Backward and forward integration initial \
            times must match"
        @assert nobs(tracklets) == length(res) "Number of observations must \
            match number of residuals"
        _bwd_ = TaylorInterpolant(bwd.t0, bwd.t, collect(bwd.x))
        _fwd_ = TaylorInterpolant(fwd.t0, fwd.t, collect(fwd.x))
        return new{T, U}(tracklets, _bwd_, _fwd_, res, fit, jacobian)
    end
end

# Outer constructor
LeastSquaresOrbit(tracklets::Vector{Tracklet{T}}, bwd::TaylorInterpolant{T, U, 2},
    fwd::TaylorInterpolant{T, U, 2}, res::Vector{OpticalResidual{T, U}},
    fit::LeastSquaresFit{T}, jacobian::Matrix{T}) where {T <: Real, U <: Number} =
    LeastSquaresOrbit{T, U}(tracklets, bwd, fwd, res, fit, jacobian)

# Print method for LeastSquaresOrbit
show(io::IO, orbit::LeastSquaresOrbit) = print(io, "Least squares orbit with ",
    length(orbit.res), " residuals")

# Definition of zero LeastSquaresOrbit
function zero(::Type{LeastSquaresOrbit{T, U}}) where {T <: Real, U <: Number}
    tracklets = Vector{Tracklet{T}}(undef, 0)
    bwd = zero(TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}})
    fwd = zero(TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}})
    res = Vector{OpticalResidual{T, U}}(undef, 0)
    fit = zero(LeastSquaresFit{T})
    jacobian = Matrix{T}(undef, 0, 0)
    return LeastSquaresOrbit{T, U}(tracklets, bwd, fwd, res, fit, jacobian)
end

iszero(orbit::LeastSquaresOrbit{T, U}) where {T <: Real, U <: Number} =
    orbit == zero(LeastSquaresOrbit{T, U})

# Evaluate integrations and residuals in fit deltas
function evalfit(orbit::LeastSquaresOrbit{T, TaylorN{T}}) where {T <: Real}
    # Fit δs
    δs = orbit.fit.x
    # Evaluate integrations
    new_bwd_x = map(x -> Taylor1(x.coeffs(δs)), orbit.bwd.x)
    new_bwd = TaylorInterpolant(orbit.bwd.t0, orbit.bwd.t, new_bwd_x)
    new_fwd_x = map(x -> Taylor1(x.coeffs(δs)), orbit.fwd.x)
    new_fwd = TaylorInterpolant(orbit.fwd.t0, orbit.fwd.t, new_fwd_x)
    # Evaluate residuals
    new_res = orbit.res(δs)

    return LeastSquaresOrbit{T, T}(orbit.tracklets, new_bwd, new_fwd, new_res,
        orbit.fit, orbit.jacobian)
end

@doc raw"""
    sigmas(::LeastSquaresOrbit)

Return the uncertainties in barycentric cartesian coordinates
at the reference epoch.
"""
sigmas(orbit::LeastSquaresOrbit) =
    sqrt.(diag(orbit.jacobian * orbit.fit.Γ * orbit.jacobian'))

@doc raw"""
    snr(::LeastSquaresOrbit)

Return the signal-to-noise ratios in barycentric cartesian coordinates
at the reference epoch.
"""
snr(orbit::LeastSquaresOrbit{T, U}) where {T <: Real, U <: Number} =
    iszero(orbit) ? zeros(T, 6) : abs.(orbit()) ./ sigmas(orbit)

@doc raw"""
    jplcompare(des::String, orbit::LeastSquaresOrbit)

Return `abs.(orbit() - R) ./ sigmas(orbit)`, where `R` is JPL's state vector of object
`des` at `orbit`'s  initial epoch.
"""
function jplcompare(des::String, orbit::LeastSquaresOrbit)
    # Load Solar System ephemerides
    loadjpleph()
    # NEOs barycentric state vector
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
    q2 = kmsec2auday(spkgeo(id, julian2etsecs(epoch(orbit) + PE.J2000), "J2000", 0)[1])
    # Absolute difference in sigma units
    return @. abs(q1 - q2) / sigmas(orbit)
end

@doc raw"""
    uncertaintyparameter(od::ODProblem{D, T}, orbit::LeastSquaresOrbit{T, T},
        params::Parameters{T}) where {D, T <: Real}

Return the Minor Planet Center uncertainty parameter.

## Arguments

- `od::ODProblem{D, T}`: an orbit determination problem.
- `orbit::LeastSquaresOrbit{T, T}`: reference orbit.
- `params::Parameters{T}`: see [`Parameters`](@ref).

!!! reference
    https://www.minorplanetcenter.net/iau/info/UValue.html

!!! warning
    This function will change the (global) `TaylorSeries` variables.
"""
function uncertaintyparameter(od::ODProblem{D, T}, orbit::LeastSquaresOrbit{T, T},
    params::Parameters{T}) where {D, T <: Real}
    # Check consistency between od and orbit
    @assert od.tracklets == orbit.tracklets
    # Epoch [Julian days TDB]
    jd0 = epoch(orbit) + PE.J2000
    # Barycentric initial conditions
    q0 = orbit(epoch(orbit))
    # Scaling factors
    scalings = abs.(q0) ./ 10^6
    # Jet transport variables
    dq = scaled_variables("dx", scalings; order = 2)
    # Origin
    x0 = zeros(T, 6)
    # Initial conditions
    q = q0 + dq
    # Propagation and residuals
    _, _, res = propres(od, jd0, q, params)
    res = @. OpticalResidual(ra(res), dec(res), wra(orbit.res), wdec(orbit.res),
        isoutlier(orbit.res))
    nobs = 2 * notout(res)
    # Covariance matrix
    Q = nms(res)
    C = (nobs/2) * TS.hessian(Q, x0)
    Γ = inv(C)
    # Osculating keplerian elements
    osc = pv2kep(q - params.eph_su(epoch(orbit)); jd = jd0, frame = :ecliptic)
    # Eccentricity
    e = osc.e
    # Gauss gravitational constant [deg]
    k_0 = 180 * k_gauss / π
    # Time of perihelion passage [julian days]
    Tp = osc.tp
    # Orbital period [days]
    P = 2π * sqrt(osc.a^3 / μ_S)
    # Projected covariance matrix
    t_car2kep = TS.jacobian([Tp, P], x0)
    Γ = t_car2kep * Γ * t_car2kep'
    # Uncertainties
    dTp, dP = sqrt.(diag(Γ))
    # Convert orbital period to years
    P = P/yr
    # In-orbit longitude runoff [arcsec / decade]
    runoff = (dTp * e(x0) + 10 * dP / P(x0)) * (k_0 / P(x0)) * 3_600 * 3
    # Uncertainty parameter
    C = log(648_000) / 9
    U = floor(Int, log(runoff)/C) + 1

    return U
end