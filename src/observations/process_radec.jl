@doc raw"""
    OpticalResidual{T <: Real, U <: Number}

An astrometric optical observed minus computed residual.

## Fields

- `ξ_α::U`: right ascension residual.
- `ξ_δ::U`: declination residual.
- `w_α::T`: right ascension weight.
- `w_δ::T`: declination weight.
- `outlier::Bool`: whether the residual is an outlier or not.
"""
@auto_hash_equals struct OpticalResidual{T <: Real, U <: Number}
    ξ_α::U
    ξ_δ::U
    w_α::T
    w_δ::T
    outlier::Bool
    # Inner constructor
    function OpticalResidual{T, U}(ξ_α::U, ξ_δ::U, w_α::T, w_δ::T,
        outlier::Bool = false) where {T <: Real, U <:  Number}
        new{T, U}(ξ_α, ξ_δ, w_α, w_δ, outlier)
    end
end

# Outer constructor
OpticalResidual(ξ_α::U, ξ_δ::U, w_α::T, w_δ::T,
    outlier::Bool = false) where {T <: Real, U <: Number} =
    OpticalResidual{T, U}(ξ_α, ξ_δ, w_α, w_δ, outlier)

# Definition of zero OpticalResidual
zero(::Type{OpticalResidual{T, U}}) where {T <: Real, U <: Number} =
    OpticalResidual{T, U}(zero(U), zero(U), zero(T), zero(T), false)

# Evaluate methods
function evaluate(res::OpticalResidual{T, TaylorN{T}}, x::Vector{T}) where {T <: Real}
    OpticalResidual(res.ξ_α(x), res.ξ_δ(x), res.w_α, res.w_δ, res.outlier)
end
(res::OpticalResidual{T, TaylorN{T}})(x::Vector{T}) where {T <: Real} = evaluate(res, x)

function evaluate(res::AbstractVector{OpticalResidual{T, TaylorN{T}}},
    x::Vector{T}) where {T <: Real}
    res_new = Vector{OpticalResidual{T, T}}(undef, length(res))
    for i in eachindex(res)
        res_new[i] = evaluate(res[i], x)
    end
    return res_new
end
(res::AbstractVector{OpticalResidual{T, TaylorN{T}}})(x::Vector{T}) where {T <: Real} =
    evaluate(res, x)

# Print method for OpticalResidual
# Examples:
# α: -138.79801 δ: -89.80025
# α: -134.79450 δ: -91.42509 (outlier)
function show(io::IO, x::OpticalResidual)
    outlier_flag = isoutlier(x) ? " (outlier)" : ""
    print(io, "α: ", @sprintf("%+.5f", cte(x.ξ_α)), " δ: ",
          @sprintf("%+.5f", cte(x.ξ_δ)), outlier_flag)
end

@doc raw"""
    unfold(ξs::AbstractVector{OpticalResidual{T, U}}) where {T <: Real, U <: Number}

Concatenate right ascension and declination residuals for an orbit fit.
"""
function unfold(ξs::AbstractVector{OpticalResidual{T, U}}) where {T <: Real, U <: Number}
    # Number of non outliers
    L = count(x -> !x.outlier, ξs)
    # Vector of residuals
    res = Vector{U}(undef, 2*L)
    # Vector of weights
    w = Vector{T}(undef, 2*L)
    # Global counter
    k = 1
    # Fill residuals and weights
    for i in eachindex(ξs)
        if !ξs[i].outlier
            # Right ascension
            res[k] = ξs[i].ξ_α
            w[k] = ξs[i].w_α
            # Declination
            res[k+L] = ξs[i].ξ_δ
            w[k+L] = ξs[i].w_δ
            # Update global counter
            k += 1
        end
    end

    return res, w
end

# Functions to get specific fields of a OpticalResidual object
ra(res::OpticalResidual) = res.ξ_α
dec(res::OpticalResidual) = res.ξ_δ
wra(res::OpticalResidual) = res.w_α
wdec(res::OpticalResidual) = res.w_δ
isoutlier(res::OpticalResidual) = res.outlier
nout(res::AbstractVector{OpticalResidual{T, U}}) where {T <: Real, U <: Number} =
    count(isoutlier, res)
notout(res::AbstractVector{OpticalResidual{T, U}}) where {T <: Real, U <: Number} =
    count(!isoutlier, res)

euclid3D(x::Vector{T}) where {T <: Real} = sqrt(dot3D(x, x))
function euclid3D(x::Vector{TaylorN{T}}) where {T <: Real}
    z, w = zero(x[1]), zero(x[1])
    @inbounds for i in eachindex(x)
        TS.zero!(w)
        for k in eachindex(x[i])
            TS.mul!(w, x[i], x[i], k)
            TS.add!(z, z, w, k)
        end
    end
    TS.zero!(w)
    for k in eachindex(z)
        TS.sqrt!(w, z, k)
    end
    return w
end

dot3D(x::Vector{T}, y::Vector{T}) where {T <: Real} = x[1]*y[1] + x[2]*y[2] + x[3]*y[3]
function dot3D(x::Vector{TaylorN{T}}, y::Vector{U}) where {T <: Real, U <: Number}
    z, w = zero(x[1]), zero(x[1])
    @inbounds for i in eachindex(x)
        TS.zero!(w)
        for k in eachindex(x[i])
            TS.mul!(w, x[i], y[i], k)
            TS.add!(z, z, w, k)
        end
    end
    return z
end
dot3D(x::Vector{T}, y::Vector{TaylorN{T}}) where {T <: Real} = dot3D(y, x)

@doc raw"""
    compute_radec(observatory::ObservatoryMPC{T}, t_r_utc::DateTime; kwargs...)
    compute_radec(obs::RadecMPC{T}; kwargs...)
    compute_radec(obs::Vector{RadecMPC{T}}; kwargs...) where {T <: Real}

Compute astrometric right ascension and declination [arcsec] for a set of MPC-formatted
observations. Corrections due to Earth orientation, LOD and polar motion are considered
in computations.

## Arguments

- `observatory::ObservatoryMPC{T}`: observation site.
- `t_r_utc::DateTime`: UTC time of astrometric observation.
- `obs::Vector{RadecMPC{T}}`: optical observation(s).

## Keyword arguments

- `niter::Int`: number of light-time solution iterations (default: `5`).
- `xve::EarthEph`: Earth ephemeris (default: `earthposvel`).
- `xvs::SunEph`: Sun ephemeris (default: `sunposvel`).
- `xva::AstEph`: asteroid ephemeris.

All ephemeris must take [et seconds since J2000] and return [barycentric position in km
and velocity in km/sec].
"""
function compute_radec(observatory::ObservatoryMPC{T}, t_r_utc::DateTime; niter::Int = 5,
                       xvs::SunEph = sunposvel, xve::EarthEph = earthposvel,
                       xva::AstEph) where {T <: Real, SunEph, EarthEph, AstEph}
    # Transform receiving time from UTC to TDB seconds since J2000
    et_r_secs = dtutc2et(t_r_utc)
    # Sun barycentric position and velocity at receive time
    rv_s_t_r = xvs(et_r_secs)
    r_s_t_r = rv_s_t_r[1:3]
    # Earth's barycentric position and velocity at receive time
    rv_e_t_r = xve(et_r_secs)
    r_e_t_r = rv_e_t_r[1:3]
    # Asteroid barycentric position and velocity at receive time
    rv_a_t_r = xva(et_r_secs)
    r_a_t_r = rv_a_t_r[1:3]
    # Compute geocentric position/velocity of receiving antenna in inertial frame [km, km/s]
    RV_r = obsposvelECI(observatory, et_r_secs)
    R_r = RV_r[1:3]
    # Receiver barycentric position and velocity at receive time
    r_r_t_r = r_e_t_r + R_r

    # Down-leg iteration
    # τ_D first approximation
    # See equation (1) of https://doi.org/10.1086/116062
    ρ_vec_r = r_a_t_r - r_r_t_r
    ρ_r = euclid3D(ρ_vec_r)
    # -R_b/c, but delay is wrt asteroid Center (Brozovic et al., 2018)
    τ_D = ρ_r/clightkms # (seconds)
    # Bounce time, new estimate
    # See equation (2) of https://doi.org/10.1086/116062
    et_b_secs = et_r_secs - τ_D

    # Allocate memmory for time-delays
    Δτ_D = zero(τ_D)               # Total time delay
    Δτ_rel_D = zero(τ_D)           # Shapiro delay
    # Δτ_corona_D = zero(τ_D)      # Delay due to Solar corona
    # Δτ_tropo_D = zero(τ_D)       # Delay due to Earth's troposphere

    for _ in 1:niter
        # Asteroid barycentric position (in km) and velocity (in km/s) at bounce time (TDB)
        rv_a_t_b = xva(et_b_secs)
        r_a_t_b = rv_a_t_b[1:3]
        v_a_t_b = rv_a_t_b[4:6]
        # Estimated position of the asteroid's center of mass relative to the recieve point
        # See equation (3) of https://doi.org/10.1086/116062
        ρ_vec_r = r_a_t_b - r_r_t_r
        # Magnitude of ρ_vec_r
        # See equation (4) of https://doi.org/10.1086/116062
        ρ_r = euclid3D(ρ_vec_r)

        # Compute down-leg Shapiro delay
        # NOTE: when using PPN, substitute 2 -> 1+γ in expressions for Shapiro delay,
        # Δτ_rel_[D|U]

        # Earth's position at t_r
        e_D_vec  = r_r_t_r - r_s_t_r
        # Heliocentric distance of Earth at t_r
        e_D = euclid3D(e_D_vec)
        # Barycentric position of Sun at estimated bounce time
        rv_s_t_b = xvs(et_b_secs)
        r_s_t_b = rv_s_t_b[1:3]
        # Heliocentric position of asteroid at t_b
        p_D_vec  = r_a_t_b - r_s_t_b
        # Heliocentric distance of asteroid at t_b
        p_D = euclid3D(p_D_vec)
        # Signal path distance (down-leg)
        q_D = ρ_r

        # Shapiro correction to time-delay
        Δτ_rel_D = shapiro_delay(e_D, p_D, q_D)  # seconds
        # Troposphere correction to time-delay
        # Δτ_tropo_D = tropo_delay(R_r, ρ_vec_r) # seconds
        # Solar corona correction to time-delay
        # Δτ_corona_D = corona_delay(constant_term.(r_a_t_b), r_r_t_r, r_s_t_r, F_tx, station_code) # seconds
        # Total time-delay
        Δτ_D = Δτ_rel_D # + Δτ_tropo_D #+ Δτ_corona_D # seconds

        # New estimate
        p_dot_23 = dot3D(ρ_vec_r, v_a_t_b)/ρ_r
        # Time delay correction
        Δt_2 = (τ_D - ρ_r/clightkms - Δτ_rel_D)/(1.0-p_dot_23/clightkms)
        # Time delay new estimate
        τ_D = τ_D - Δt_2
        # Bounce time, new estimate
        # See equation (2) of https://doi.org/10.1086/116062
        et_b_secs = et_r_secs - τ_D
    end

    # Asteroid barycentric position (in km) and velocity (in km/s) at bounce time (TDB)
    rv_a_t_b = xva(et_b_secs)
    r_a_t_b = rv_a_t_b[1:3]
    v_a_t_b = rv_a_t_b[4:6]
    # Estimated position of the asteroid's center of mass relative to the recieve point
    # See equation (3) of https://doi.org/10.1086/116062
    ρ_vec_r = r_a_t_b - r_r_t_r
    # Magnitude of ρ_vec_r
    # See equation (4) of https://doi.org/10.1086/116062
    ρ_r = euclid3D(ρ_vec_r)

    # TODO: add aberration and atmospheric refraction corrections

    # Compute gravitational deflection of light
    # See Explanatory Supplement to the Astronomical Almanac (ESAA) 2014 Section 7.4.1.4
    E_H_vec = r_r_t_r -r_s_t_r    # ESAA 2014, equation (7.104)
    U_vec = ρ_vec_r               # r_a_t_b - r_e_t_r, ESAA 2014, equation (7.112)
    U_norm = ρ_r                  # sqrt(U_vec[1]^2 + U_vec[2]^2 + U_vec[3]^2)
    u_vec = U_vec/U_norm
    # Barycentric position and velocity of Sun at converged bounce time
    rv_s_t_b = xvs(et_b_secs)
    r_s_t_b = rv_s_t_b[1:3]
    Q_vec = r_a_t_b - r_s_t_b     # ESAA 2014, equation (7.113)
    q_vec = Q_vec/ euclid3D(Q_vec)
    E_H = euclid3D(E_H_vec)
    e_vec = E_H_vec/E_H
    # See ESAA 2014, equation (7.115)
    g1 = (2μ_DE430[su]/(c_au_per_day^2))/(E_H/au)
    g2 = 1 + dot3D(q_vec, e_vec)
    # See ESAA 2014, equation (7.116)
    u1_vec = U_norm*(  u_vec + (g1/g2)*( dot3D(u_vec,q_vec)*e_vec - dot3D(e_vec,u_vec)*q_vec )  )
    u1_norm = euclid3D(u1_vec)

    # Compute right ascension, declination angles
    α_rad_ = mod2pi(atan(u1_vec[2], u1_vec[1]))
    α_rad = mod2pi(α_rad_)          # right ascension (rad)
    δ_rad = asin(u1_vec[3]/u1_norm) # declination (rad)

    δ_as = rad2arcsec(δ_rad) # rad -> arcsec + debiasing
    α_as = rad2arcsec(α_rad) # rad -> arcsec + debiasing

    return α_as, δ_as # right ascension, declination both in arcsec
end

compute_radec(obs::RadecMPC{T}; kwargs...) where {T <: Real} =
    compute_radec(obs.observatory, obs.date; kwargs...)

function compute_radec(obs::Vector{RadecMPC{T}}; xva::AstEph,
    kwargs...) where {T <: Real, AstEph}

    # Number of observations
    n_optical_obs = length(obs)
    # UTC time of first astrometric observation
    utc1 = obs[1].date
    # TDB seconds since J2000.0 for first astrometric observation
    et1 = dtutc2et(utc1)
    # Asteroid ephemeris at et1
    a1_et1 = xva(et1)[1]
    # Type of asteroid ephemeris
    S = typeof(a1_et1)

    # Right ascension
    vra = Array{S}(undef, n_optical_obs)
    # Declination
    vdec = Array{S}(undef, n_optical_obs)

    # Iterate over the number of observations
    for i in 1:n_optical_obs
        vra[i], vdec[i] = compute_radec(obs[i]; xva = xva, kwargs...)
    end

    return vra, vdec # arcsec, arcsec
end

# Angle difference taking into account the discontinuity in [0, 2π) -> [0, 2π)
# x and y must be in arcsec
function anglediff(x::T, y::S) where {T, S <: Number}
    # Signed difference
    Δ = x - y
    # Absolute difference
    Δ_abs = abs(Δ)
    # Reflection
    if Δ_abs > 648_000 # Half circle in arcsec
        return -sign(cte(Δ)) * (1_296_000 - Δ_abs)
    else
        return Δ
    end
end

@doc raw"""
    residuals(radec::AbstractVector{RadecMPC{T}}, w8s::AbstractVector{T},
        bias::AbstractVector{Tuple{T, T}}; xva::AstEph,
        kwargs...) where {AstEph, T <: Real}

Compute observed minus computed residuals for optical astrometry `radec`.
Corrections due to Earth orientation, LOD and polar motion are computed
by default.

See also [`OpticalResidual`](@ref) and [`compute_radec`](@ref).

## Arguments

- `radec::AbstractVector{RadecMPC{T}}`: optical astrometry.
- `w8s::AbstractVector{T}`: statistical weights.
- `bias::AbstractVector{Tuple{T, T}}`: debiasing corrections.

## Keyword arguments

- `niter::Int`: number of light-time solution iterations (default: `5`).
- `xve::EarthEph`: Earth ephemeris (default: `earthposvel`).
- `xvs::SunEph`: Sun ephemeris (default: `sunposvel`).
- `xva::AstEph`: asteroid ephemeris.

All ephemeris must take [et seconds since J2000] and return [barycentric position in km
and velocity in km/sec].
"""
function residuals(radec::AbstractVector{RadecMPC{T}}, w8s::AbstractVector{T},
    bias::AbstractVector{Tuple{T, T}}; xva::AstEph, kwargs...) where {AstEph, T <: Real}

    # Check consistency between arrays
    @assert length(radec) == length(w8s) == length(bias)
    # UTC time of first astrometric observation
    utc1 = date(radec[1])
    # TDB seconds since J2000.0 for first astrometric observation
    et1 = dtutc2et(utc1)
    # Asteroid ephemeris at et1
    a1_et1 = xva(et1)[1]
    # Type of asteroid ephemeris
    U = typeof(a1_et1)
    # Vector of residuals
    res = [zero(OpticalResidual{T, U}) for _ in eachindex(radec)]
    residuals!(res, radec, w8s, bias; xva, kwargs...)

    return res
end

function residuals!(res::Vector{OpticalResidual{T, U}},
    radec::AbstractVector{RadecMPC{T}}, w8s::AbstractVector{T},
    bias::AbstractVector{Tuple{T, T}}; xva::AstEph, kwargs...) where {AstEph,
    T <: Real, U <: Number}

    tmap!(res, radec, w8s, bias, isoutlier.(res)) do obs, w8, bias, outlier
        # Observed ra/dec
        α_obs = rad2arcsec(ra(obs))   # arcsec
        δ_obs = rad2arcsec(dec(obs))  # arcsec
        # Computed ra/dec
        α_comp, δ_comp = compute_radec(obs; xva = xva, kwargs...)   # arcsec
        # Debiasing corrections
        α_corr, δ_corr = bias
        # Observed minus computed residual ra/dec
        # Note: ra is multiplied by a metric factor cos(dec) to match the format of
        # debiasing corrections
        return OpticalResidual{T, U}(
            anglediff(α_obs, α_comp) * cos(dec(obs)) - α_corr,
            δ_obs - δ_comp - δ_corr,
            w8,
            w8,
            outlier
        )
    end

    return nothing
end