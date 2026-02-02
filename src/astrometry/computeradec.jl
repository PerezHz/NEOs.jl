"""
    OpticalBuffer{U <: Number} <: AbstractBuffer

Pre-allocated memory for [`compute_radec`](@ref).

# Fields

- `v0::Vector{U}`: array of scalar variables.
- `v1::Vector{Vector{U}}`: array of vector variables.
"""
struct OpticalBuffer{U <: Number} <: AbstractBuffer
    v0::Vector{U}
    v1::Vector{Vector{U}}
end

function OpticalBuffer(x::U) where {U <: Number}
    # Scalars
    ρ_r = zero(x)
    τ_D = zero(x)
    et_b_secs = zero(x)
    Δτ_D = zero(x)
    Δτ_rel_D = zero(x)
    p_D = zero(x)
    q_D = zero(x)
    p_dot_23 = zero(x)
    Δt_2 = zero(x)
    U_norm = zero(x)
    g2 = zero(x)
    u1_norm = zero(x)
    α_rad_ = zero(x)
    α_rad = zero(x)
    δ_rad = zero(x)
    δ_as = zero(x)
    α_as = zero(x)
    # Vectors
    rv_a_t_r = [zero(x) for _ in 1:6]
    r_a_t_r = [zero(x) for _ in 1:3]
    ρ_vec_r = [zero(x) for _ in 1:3]
    rv_a_t_b = [zero(x) for _ in 1:6]
    r_a_t_b = [zero(x) for _ in 1:3]
    v_a_t_b = [zero(x) for _ in 1:3]
    rv_s_t_b = [zero(x) for _ in 1:6]
    r_s_t_b = [zero(x) for _ in 1:3]
    p_D_vec  = [zero(x) for _ in 1:3]
    U_vec = [zero(x) for _ in 1:3]
    u_vec = [zero(x) for _ in 1:3]
    Q_vec = [zero(x) for _ in 1:3]
    q_vec = [zero(x) for _ in 1:3]
    u1_vec = [zero(x) for _ in 1:3]

    return OpticalBuffer{U}(
        [ρ_r, τ_D, et_b_secs, Δτ_D, Δτ_rel_D, p_D, q_D, p_dot_23, Δt_2,
         U_norm, g2, u1_norm, α_rad_, α_rad, δ_rad, δ_as, α_as],
        [rv_a_t_r, r_a_t_r, ρ_vec_r, rv_a_t_b, r_a_t_b, v_a_t_b, rv_s_t_b,
         r_s_t_b, p_D_vec, U_vec, u_vec, Q_vec, q_vec, u1_vec],
    )
end

"""
    compute_radec(::AbstractOpticalAstrometry [, ::OpticalBuffer]; kwargs...)

Compute the astrometric right ascension and declination [arcsec] at the observatory
and date given by an optical observation. An optional buffer can be passed to recycle
memory. Corrections due to Earth orientation, LOD and polar motion are considered.

# Keyword arguments

- `niter::Int`: number of light-time solution iterations (default: `5`).
- `xve`: Earth ephemeris (default: `earthposvel`).
- `xvs`: Sun ephemeris (default: `sunposvel`).
- `xva`: asteroid ephemeris.

All ephemeris must take [et seconds since J2000] and return [barycentric position
in km and velocity in km/sec].
"""
compute_radec(x::AbstractOpticalAstrometry; kwargs...) =
    compute_radec(observatory(x), date(x); kwargs...)

compute_radec(x::AbstractOpticalAstrometry, buffer::OpticalBuffer; kwargs...) =
    compute_radec(observatory(x), date(x), buffer; kwargs...)

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

    # Allocate memory for time-delays
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

function compute_radec(observatory::ObservatoryMPC{T}, t_r_utc::DateTime,
                       buffer::OpticalBuffer{U}; niter::Int = 5,
                       xvs::DensePropagation2{T, T}, xve::DensePropagation2{T, T},
                       xva::NTuple{2, DensePropagation2{T, U}}) where {T <: Real, U <: Number}
    # Unfold
    ρ_r, τ_D, et_b_secs, Δτ_D, Δτ_rel_D, p_D, q_D, p_dot_23, Δt_2, U_norm,
        g2, u1_norm, α_rad_, α_rad, δ_rad, δ_as, α_as = buffer.v0
    rv_a_t_r, r_a_t_r, ρ_vec_r, rv_a_t_b, r_a_t_b, v_a_t_b, rv_s_t_b,
        r_s_t_b, p_D_vec, U_vec, u_vec, Q_vec, q_vec, u1_vec = buffer.v1
    # Transform receiving time from UTC to TDB seconds since J2000
    et_r_secs = dtutc2et(t_r_utc)
    # Sun barycentric position and velocity at receive time
    rv_s_t_r = auday2kmsec(xvs(et_r_secs/daysec))
    r_s_t_r = rv_s_t_r[1:3]
    # Earth's barycentric position and velocity at receive time
    rv_e_t_r = auday2kmsec(xve(et_r_secs/daysec))
    r_e_t_r = rv_e_t_r[1:3]
    # Asteroid barycentric position and velocity at receive time
    evaleph!(rv_a_t_r, et_r_secs, xva[1], xva[2])
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

    # Allocate memory for time-delays
    Δτ_D = zero(τ_D)               # Total time delay
    Δτ_rel_D = zero(τ_D)           # Shapiro delay
    # Δτ_corona_D = zero(τ_D)      # Delay due to Solar corona
    # Δτ_tropo_D = zero(τ_D)       # Delay due to Earth's troposphere

    for _ in 1:niter
        # Asteroid barycentric position (in km) and velocity (in km/s) at bounce time (TDB)
        evaleph!(rv_a_t_b, et_b_secs, xva[1], xva[2])
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
        evaleph!(rv_s_t_b, et_b_secs, xvs)
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
    evaleph!(rv_a_t_b, et_b_secs, xva[1], xva[2])
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
    evaleph!(rv_s_t_b, et_b_secs, xvs)
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