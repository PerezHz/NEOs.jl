"""
    OpticalResidual{T, U} <: AbstractOpticalResidual{T, U}

An astrometric optical observed minus computed residual.

# Fields

- `ra/dec::U`: right ascension (declination) residual [arcsec].
- `wra/wdec::T`: right ascension (declination) weight [arcsec⁻²].
- `dra/ddec::T`: right ascension (declination) debiasing factor [arcsec].
- `outlier::Bool`: whether the residual is an outlier or not.
"""
@auto_hash_equals struct OpticalResidual{T, U} <: AbstractOpticalResidual{T, U}
    ra::U
    dec::U
    wra::T
    wdec::T
    dra::T
    ddec::T
    outlier::Bool
    # Inner constructor
    function OpticalResidual{T, U}(ra::U, dec::U, wra::T, wdec::T, dra::T, ddec::T,
                                   outlier::Bool = false) where {T <: Real, U <:  Number}
        return new{T, U}(ra, dec, wra, wdec, dra, ddec, outlier)
    end
end

# Definition of zero OpticalResidual
zero(::Type{OpticalResidual{T, U}}) where {T, U} = OpticalResidual{T, U}(
        zero(U), zero(U), zero(T), zero(T), zero(T), zero(T), false)
iszero(x::OpticalResidual{T, U}) where {T, U} = x == zero(OpticalResidual{T, U})

ra(x::OpticalResidual) = x.ra
dec(x::OpticalResidual) = x.dec
wra(x::OpticalResidual) = x.wra
wdec(x::OpticalResidual) = x.wdec
dra(x::OpticalResidual) = x.dra
ddec(x::OpticalResidual) = x.ddec

isoutlier(x::OpticalResidual) = x.outlier
nout(x::AbstractVector{OpticalResidual{T, U}}) where {T, U} = count(isoutlier, x)
notout(x::AbstractVector{OpticalResidual{T, U}}) where {T, U} = count(!isoutlier, x)

# Print method for OpticalResidual
function show(io::IO, x::OpticalResidual)
    outlier_flag = isoutlier(x) ? " (outlier)" : ""
    print(io, "α: ", @sprintf("%+.5f", cte(x.ξ_α)), " δ: ",
        @sprintf("%+.5f", cte(x.ξ_δ)), outlier_flag)
end

# Evaluate methods
evaluate(y::OpticalResidual{T, TaylorN{T}}, x::Vector{T}) where {T <: Real} =
    OpticalResidual{T, T}(y.ra(x), y.dec(x), y.wra, y.wdec, y.dra, y.ddec, y.outlier)

(y::OpticalResidual{T, TaylorN{T}})(x::Vector{T}) where {T <: Real} = evaluate(y, x)

function evaluate(y::AbstractVector{OpticalResidual{T, TaylorN{T}}},
                  x::Vector{T}) where {T <: Real}
    z = Vector{OpticalResidual{T, T}}(undef, length(y))
    for i in eachindex(z)
        z[i] = evaluate(y[i], x)
    end
    return z
end

(y::AbstractVector{OpticalResidual{T, TaylorN{T}}})(x::Vector{T}) where {T <: Real} =
    evaluate(y, x)

"""
    unfold(::AbstractVector{OpticalResidual})

Return three vectors by concatenating the non-outlier right ascension
and declination residuals, weights and debiasing factors.
"""
function unfold(y::AbstractVector{OpticalResidual{T, U}}) where {T <: Real, U <: Number}
    # Number of non outliers
    L = notout(y)
    # Vector of residuals, weights and debiasing factors
    z = Vector{U}(undef, 2L)
    w = Vector{T}(undef, 2L)
    d = Vector{T}(undef, 2L)
    # Global counter
    k = 1
    # Fill residuals, weights and debiasing factors
    for i in eachindex(y)
        isoutlier(y[i]) && continue
        # Right ascension
        z[k], w[k], d[k] = ra(y[i]), wra(y[i]), dra(y[i])
        # Declination
        z[k+L], w[k+L], d[k+L] = dec(y[i]), wdec(y[i]), ddec(y[i])
        # Update global counter
        k += 1
    end

    return z, w, d
end

# Warning: functions euclid3D(x) and dot3D(x) assume length(x) >= 3
euclid3D(x::AbstractVector{T}) where {T <: Real} = sqrt(dot3D(x, x))

function euclid3D(x::AbstractVector{TaylorN{T}}) where {T <: Real}
    z, w = zero(x[1]), zero(x[1])
    @inbounds for i in 1:3
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

dot3D(x::AbstractVector{T}, y::AbstractVector{T}) where {T <: Real} =
    x[1]*y[1] + x[2]*y[2] + x[3]*y[3]

function dot3D(x::Vector{TaylorN{T}}, y::Vector{U}) where {T <: Real, U <: Number}
    z, w = zero(x[1]), zero(x[1])
    @inbounds for i in 1:3
        TS.zero!(w)
        for k in eachindex(x[i])
            TS.mul!(w, x[i], y[i], k)
            TS.add!(z, z, w, k)
        end
    end
    return z
end

dot3D(x::Vector{T}, y::Vector{TaylorN{T}}) where {T <: Real} = dot3D(y, x)

"""
    compute_radec(::AbstractOpticalResidual; xva, kwargs...)

Compute the astrometric right ascension and declination [arcsec]. Corrections
due to Earth orientation, LOD and polar motion are considered.

# Keyword arguments

- `niter::Int`: number of light-time solution iterations (default: `5`).
- `xve::EarthEph`: Earth ephemeris (default: `earthposvel`).
- `xvs::SunEph`: Sun ephemeris (default: `sunposvel`).
- `xva::AstEph`: asteroid ephemeris.

All ephemeris must take [et seconds since J2000] and return [barycentric position
in km and velocity in km/sec].
"""
compute_radec(x::AbstractOpticalResidual; kwargs...) =
    compute_radec(observatory(x), date(x); kwargs...)

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

        # Shapiro correction to time-delay [seconds]
        Δτ_rel_D = shapiro_delay(e_D, p_D, q_D)
        # Troposphere correction to time-delay [seconds]
        # Δτ_tropo_D = tropo_delay(R_r, ρ_vec_r)
        # Solar corona correction to time-delay [seconds]
        # Δτ_corona_D = corona_delay(constant_term.(r_a_t_b),
        #   r_r_t_r, r_s_t_r, F_tx, station_code)
        # Total time-delay [seconds]
        Δτ_D = Δτ_rel_D # + Δτ_tropo_D + Δτ_corona_D

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
    u1_vec = U_norm * ( u_vec + (g1/g2) * ( dot3D(u_vec, q_vec) * e_vec -
        dot3D(e_vec, u_vec) * q_vec ) )
    u1_norm = euclid3D(u1_vec)

    # Compute right ascension, declination angles
    α_rad_ = mod2pi(atan(u1_vec[2], u1_vec[1]))
    α_rad = mod2pi(α_rad_)          # right ascension (rad)
    δ_rad = asin(u1_vec[3]/u1_norm) # declination (rad)

    δ_as = rad2arcsec(δ_rad) # rad -> arcsec + debiasing
    α_as = rad2arcsec(α_rad) # rad -> arcsec + debiasing

    return α_as, δ_as # right ascension, declination both in arcsec
end

function compute_radec(x::AbstractOpticalVector; xva::AstEph, kwargs...) where {AstEph}
    # Number of observations
    N = length(x)
    # UTC time of first astrometric observation
    utc1 = date(x[1])
    # TDB seconds since J2000.0 for first astrometric observation
    et1 = dtutc2et(utc1)
    # Asteroid ephemeris at et1
    a1_et1 = xva(et1)[1]
    # Type of asteroid ephemeris
    S = typeof(a1_et1)
    # Computed right ascension and declination [arcsec]
    ra = Array{S}(undef, N)
    dec = Array{S}(undef, N)
    for i in eachindex(x)
        ra[i], dec[i] = compute_radec(x[i]; xva, kwargs...)
    end

    return ra, dec
end

# Angle difference taking into account the discontinuity in [0, 2π) -> [0, 2π)
# x and y must be in arcsec
function anglediff(x::Number, y::Number)
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

function init_residuals(::Type{U},
                        optical::AbstractOpticalVector{T},
                        w8s::AbstractVector{Tuple{T, T}},
                        bias::AbstractVector{Tuple{T, T}},
                        outliers::AbstractVector{Bool}) where {T <: Real, U <: Number}
    # Check consistency between arrays
    @assert length(optical) == length(w8s) == length(bias) == length(outliers)
    # Initialize vector of residuals
    res = Vector{OpticalResidual{T, U}}(undef, length(optical))
    for i in eachindex(optical)
        ra, dec = zero(U), zero(U)
        wra, wdec = w8s[i]
        dra, ddec = bias[i]
        res[i] = OpticalResidual{T, U}(ra, dec, wra, wdec, dra, ddec, outliers[i])
    end

    return res
end

"""
    residuals(optical, w8s, bias [, outliers]; xva, kwargs...) where {AstEph, T <: Real}

Compute the observed minus computed residuals for a vector of optical astrometry.
Corrections due to Earth orientation, LOD and polar motion are computed by default.

See also [`OpticalResidual`](@ref) and [`compute_radec`](@ref).

# Arguments

- `optical::AbstractOpticalVector{T}`: optical astrometry.
- `w8s::AbstractVector{Tuple{T, T}}`: statistical weights.
- `bias::AbstractVector{Tuple{T, T}}`: debiasing corrections.
- `outliers::AbstractVector{Bool}`: outlier flags (default:
    `fill(false, length(optical))`).

# Keyword arguments

- `niter::Int`: number of light-time solution iterations (default: `5`).
- `xve::EarthEph`: Earth ephemeris (default: `earthposvel`).
- `xvs::SunEph`: Sun ephemeris (default: `sunposvel`).
- `xva::AstEph`: asteroid ephemeris.

All ephemeris must take [et seconds since J2000] and return [barycentric position in km
and velocity in km/sec].
"""
function residuals(optical::AbstractOpticalVector{T}, w8s::AbstractVector{Tuple{T, T}},
                   bias::AbstractVector{Tuple{T, T}}, outliers::AbstractVector{Bool} =
                   fill(false, length(optical)); xva::AstEph, kwargs...) where {AstEph, T <: Real}
    # UTC time of first optical observation
    utc1 = date(optical[1])
    # TDB seconds since J2000.0 for first optical observation
    et1 = dtutc2et(utc1)
    # Asteroid ephemeris at et1
    a1_et1 = xva(et1)[1]
    # Type of asteroid ephemeris
    U = typeof(a1_et1)
    # Vector of residuals
    res = init_residuals(U, optical, w8s, bias, outliers)
    residuals!(res, optical; xva, kwargs...)

    return res
end

function residuals!(res::Vector{OpticalResidual{T, U}}, optical::AbstractOpticalVector{T};
                    xva::AstEph, kwargs...) where {AstEph, T <: Real, U <: Number}

    @allow_boxed_captures tmap!(res, optical, wra.(res), wdec.(res), dra.(res),
        ddec.(res), isoutlier.(res)) do x, w8ra, w8dec, biasra, biasdec, outlier
        # Observed ra/dec [arcsec]
        obsra, obsdec = rad2arcsec.(measure(x))
        # Computed ra/dec [arcsec]
        compra, compdec = compute_radec(x; xva, kwargs...)
        # Observed minus computed residual ra/dec
        # Note: ra is multiplied by a metric factor cos(dec) to match the format of
        # debiasing corrections
        return OpticalResidual{T, U}(
            anglediff(obsra, compra) * cos(dec(x)) - biasra,
            obsdec - compdec - biasdec,
            w8ra,
            w8dec,
            biasra,
            biasdec,
            outlier
        )
    end

    return nothing
end