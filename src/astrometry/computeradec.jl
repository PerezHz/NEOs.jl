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
    v0 = [zero(x) for _ in 1:32]
    v1 = [[zero(x) for _ in 1:6] for _ in 1:19]
    return OpticalBuffer{U}(v0, v1)
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

function compute_radec(
        observatory::ObservatoryMPC{T}, t_r_utc::DateTime; niter::Int = 5,
        xvs::SunEph = sunposvel, xve::EarthEph = earthposvel, xva::AstEph
    ) where {T <: Real, SunEph, EarthEph, AstEph}
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

function compute_radec(
        observatory::ObservatoryMPC{T}, t_r_utc::DateTime, buffer::OpticalBuffer{U};
        niter::Int = 5, xvs::DensePropagation2{T, T}, xve::DensePropagation2{T, T},
        xva::NTuple{2, DensePropagation2{T, U}}
    ) where {T <: Real, U <: NumberNotSeries}
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
    E_H_vec = r_r_t_r - r_s_t_r    # ESAA 2014, equation (7.104)
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

# Taylorized version of the function above
function compute_radec(
        observatory::ObservatoryMPC{T}, t_r_utc::DateTime, buffer::OpticalBuffer{U};
        niter::Int = 5, xvs::DensePropagation2{T, T}, xve::DensePropagation2{T, T},
        xva::NTuple{2, DensePropagation2{T, U}}
    ) where {T <: Real, U <: AbstractSeries}
    aux1, aux2, ρ_r, τ_D, et_b_secs, Δτ_D, Δτ_rel_D, p_D, q_D, _p_dot_23_,
    p_dot_23, τ_D_r, τ_D_23, one_minus_τ_D_23, _Δt_2_, Δt_2, U_norm, Q_norm,
    q_dot_e, g2, u_dot_q, e_dot_u, g1_div_g2, u1_norm, u12_div_u11, α_aux,
    u13_div_u1_norm, δ_aux, α_rad,  δ_rad, δ_as, α_as = buffer.v0
    rv_a_t_r, r_a_t_r, ρ_vec_r, rv_a_t_b, r_a_t_b, v_a_t_b, rv_s_t_b,
    r_s_t_b, p_D_vec, U_vec, u_vec, Q_vec, q_vec, u_dot_q_e_vec,
    e_dot_u_q_vec, uqe_minus_euq, g12_uqe_vec, _u1_vec_, u1_vec = buffer.v1
    et_r_secs = dtutc2et(t_r_utc)
    rv_s_t_r = auday2kmsec(xvs(et_r_secs/daysec))
    r_s_t_r = rv_s_t_r[1:3]
    rv_e_t_r = auday2kmsec(xve(et_r_secs/daysec))
    r_e_t_r = rv_e_t_r[1:3]
    evaleph!(rv_a_t_r, et_r_secs, xva[1], xva[2])
    r_a_t_r = rv_a_t_r[1:3]
    order = get_order(r_a_t_r[1])
    RV_r = obsposvelECI(observatory, et_r_secs)
    R_r = RV_r[1:3]
    r_r_t_r = r_e_t_r + R_r
    E_H_vec = e_D_vec  = r_r_t_r - r_s_t_r
    E_H = e_D = euclid3D(e_D_vec)
    e_vec = E_H_vec / E_H
    g1 = (2μ_DE430[su] /(c_au_per_day^2)) / (E_H / au)
    for ord = 0:order
        for i in 1:3
            TS.subst!(ρ_vec_r[i], r_a_t_r[i], r_r_t_r[i], ord)
        end
        euclid3D!(ρ_r, ρ_vec_r, aux1, aux2, ord)
        TS.mul!(τ_D, 1/clightkms, ρ_r, ord)
        TS.subst!(et_b_secs, et_r_secs, τ_D, ord)
    end
    for _ in 1:niter
        evaleph!(rv_a_t_b, et_b_secs, xva[1], xva[2])
        evaleph!(rv_s_t_b, et_b_secs, xvs)
        for ord in 0:order
            for i in 1:3
                TS.identity!(r_a_t_b[i], rv_a_t_b[i], ord)
                TS.identity!(v_a_t_b[i], rv_a_t_b[i+3], ord)
                TS.identity!(r_s_t_b[i], rv_s_t_b[i], ord)
                TS.subst!(ρ_vec_r[i], r_a_t_b[i], r_r_t_r[i], ord)
                TS.subst!(p_D_vec[i], r_a_t_b[i], r_s_t_b[i], ord)
            end
            euclid3D!(ρ_r, ρ_vec_r, aux1, aux2, ord)
            euclid3D!(p_D, p_D_vec, aux1, aux2, ord)
            TS.identity!(q_D, ρ_r, ord)
        end
        Δτ_rel_D = shapiro_delay(e_D, p_D, q_D)
        for ord in 0:order
            TS.identity!(Δτ_D, Δτ_rel_D, ord)
            dot3D!(_p_dot_23_, ρ_vec_r, v_a_t_b, aux1, ord)
            TS.div!(p_dot_23, _p_dot_23_, ρ_r, ord)
            TS.mul!(τ_D_r, 1/clightkms, ρ_r, ord)
            TS.mul!(τ_D_23, 1/clightkms, p_dot_23, ord)
            TS.subst!(one_minus_τ_D_23, 1, τ_D_23, ord)
            TS.subst!(_Δt_2_, τ_D, τ_D_r, ord)
            TS.subst!(_Δt_2_, _Δt_2_, Δτ_rel_D, ord)
            TS.div!(Δt_2, _Δt_2_, one_minus_τ_D_23, ord)
            TS.subst!(τ_D, τ_D, Δt_2, ord)
            TS.subst!(et_b_secs, et_r_secs, τ_D, ord)
        end
    end
    evaleph!(rv_a_t_b, et_b_secs, xva[1], xva[2])
    evaleph!(rv_s_t_b, et_b_secs, xvs)
    for ord in 0:order
        for i in 1:3
            TS.identity!(r_a_t_b[i], rv_a_t_b[i], ord)
            TS.identity!(v_a_t_b[i], rv_a_t_b[i+3], ord)
            TS.identity!(r_s_t_b[i], rv_s_t_b[i], ord)
            TS.subst!(ρ_vec_r[i], r_a_t_b[i], r_r_t_r[i], ord)
            TS.identity!(U_vec[i], ρ_vec_r[i], ord)
            TS.subst!(Q_vec[i], r_a_t_b[i], r_s_t_b[i], ord)
        end
        euclid3D!(ρ_r, ρ_vec_r, aux1, aux2, ord)
        TS.identity!(U_norm, ρ_r, ord)
        euclid3D!(Q_norm, Q_vec, aux1, aux2, ord)
        for i in 1:3
            TS.div!(u_vec[i], U_vec[i], U_norm, ord)
            TS.div!(q_vec[i], Q_vec[i], Q_norm, ord)
        end
        dot3D!(q_dot_e, q_vec, e_vec, aux1, ord)
        TS.add!(g2, 1, q_dot_e, ord)
        TS.div!(g1_div_g2, g1, g2, ord)
        dot3D!(u_dot_q, u_vec, q_vec, aux1, ord)
        dot3D!(e_dot_u, e_vec, u_vec, aux1, ord)
        for i in 1:3
            TS.zero!(e_dot_u_q_vec[i], ord)
            TS.mul!(e_dot_u_q_vec[i], e_dot_u, q_vec[i], ord)
            TS.zero!(u_dot_q_e_vec[i], ord)
            TS.mul!(u_dot_q_e_vec[i], u_dot_q, e_vec[i], ord)
            TS.subst!(uqe_minus_euq[i], u_dot_q_e_vec[i], e_dot_u_q_vec[i], ord)
            TS.zero!(g12_uqe_vec[i], ord)
            TS.mul!(g12_uqe_vec[i], g1_div_g2, uqe_minus_euq[i], ord)
            TS.add!(_u1_vec_[i], u_vec[i], g12_uqe_vec[i], ord)
            TS.zero!(u1_vec[i], ord)
            TS.mul!(u1_vec[i], U_norm, _u1_vec_[i], ord)
        end
        euclid3D!(u1_norm, u1_vec, aux1, aux2, ord)
        TS.div!(u12_div_u11, u1_vec[2], u1_vec[1], ord)
        TS.div!(u13_div_u1_norm, u1_vec[3], u1_norm, ord)
        TS.atan!(α_rad, u12_div_u11, α_aux, ord)
        if ord == 0
            α_rad[0] = mod2pi(mod2pi(atan(cte(u1_vec[2]), cte(u1_vec[1]))))
        end
        TS.mul!(α_as, rad2deg(one(T)), α_rad, ord)
        TS.mul!(α_as, 3_600, α_as, ord)
        TS.asin!(δ_rad, u13_div_u1_norm, δ_aux, ord)
        TS.mul!(δ_as, rad2deg(one(T)), δ_rad, ord)
        TS.mul!(δ_as, 3_600, δ_as, ord)
    end
    return α_as, δ_as
end