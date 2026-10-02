"""
    RadarBuffer{U <: Number} <: AbstractBuffer

Pre-allocated memory for [`compute_delay`](@ref).

# Fields

- `v0::Vector{Taylor1{U}}`: array of scalar variables.
- `v1::Vector{Vector{Taylor1{U}}}`: array of vector variables.

!!! note
    All variables are `Taylor1` expansions in time (around the echo reception
    time) with coefficients of type `U`. The first element of `v0` is the
    independent time variable.
"""
struct RadarBuffer{U <: Number} <: AbstractBuffer
    v0::Vector{Taylor1{U}}
    v1::Vector{Vector{Taylor1{U}}}
end

function RadarBuffer(x::U, tord::Int) where {U <: Number}
    @assert tord ≥ 1 "Order of Taylor expansions must be at least one"
    # Every coefficient must be a different object, since the buffer
    # is modified in place
    zeroT1() = Taylor1([zero(x) for _ in 0:tord], tord)
    # Independent time variable
    tvar = Taylor1([zero(x), one(x), (zero(x) for _ in 2:tord)...], tord)
    v0 = [tvar, (zeroT1() for _ in 1:24)...]
    # Number of components of each vector variable
    nv1 = (3, 3, 3, 3, 3, 3, 6, 3, 3, 3, 3, 6, 3, 3, 3, 3, 3, 3)
    v1 = [[zeroT1() for _ in 1:n] for n in nv1]
    return RadarBuffer{U}(v0, v1)
end

raw"""
    shapiro_delay(e, p, q)

Return the relativistic (Shapiro) time-delay [sec].

See also [`shapiro_doppler`](@ref).

# Arguments

- `e`: heliocentric distance of the Earth.
- `p`: asteroid's heliocentric distance.
- `q`: asteroid's geocentric distance.

!!! reference
    See
    - https://doi.org/10.1103/PhysRevLett.13.789

# Extended help

The Shapiro time-delay is given by:
```math
\Delta\tau[\text{rel}] = \frac{2\mu_\odot}{c^3}\log\left|\frac{d_{E,S} +
    d_{A,S} + d_{A,E}}{d_{E,S}+d_{A,S}-d_{A, E}}\right|,
```
where ``\mu_\odot = GM_\odot`` is the gravitational parameter of the sun,
and ``d_{E,S}``, ``d_{A,S}`` and ``d_{A,E}`` are the heliocentric distance
of the Earth, the asteroid's heliocentric distance, and the asteroid's
geocentric distance, respectively.
"""
function shapiro_delay(e, p, q)
    # Ansatz
    shap = 0.0 # 2μ[1]/(c_au_per_day^2)
    # Shapiro time-delay [days]
    shap_del_days = (2μ_DE430[su] / (c_au_per_day^3)) * log( (e+p+q+shap)/(e+p-q+shap) )
    # Shapiro time-delay [sec]
    return shap_del_days * daysec
end

raw"""
    shapiro_doppler(e, de, p, dp, q, dq, F_tx)

Return the (Shapiro) Doppler shift [units of `F_tx`].

See also [`shapiro_delay`](@ref).

# Arguments

- `e`: heliocentric distance of the Earth.
- `de`: differential of `e`.
- `p`: asteroid's heliocentric distance.
- `dp`: differential of `p`.
- `q`: asteroid's geocentric distance.
- `dq`: differential of `q`.
- `F_tx`: transmitter frequency [MHz].

!!! reference
    See
    - https://doi.org/10.1103/PhysRevLett.17.933

# Extended help

The Doppler shift is given by:
```math
\Delta\nu = -\nu\frac{d\Delta\tau}{dt},
```
where ``\nu`` is the frequency and ``\frac{d\Delta\tau}{dt}`` is the
differential of the Shapiro delay.
"""
function shapiro_doppler(e, de, p, dp, q, dq, F_tx)
    # Differential of Shapiro delay [adim]
    # shap_del_diff = 2μ[1]*( (de+dp+dq)/(e+p+q) - (de+dp-dq)/(e+p-q) )/(c_au_per_day^3)
    shap_del_diff = (4μ_DE430[su] / (c_au_per_day^3)) * (  ( dq*(e+p) - q*(de+dp) ) /
        ( (e+p)^2 - q^2 )  )
    # ν = -F_tx * dτ / dt [units of F_tx]
    # Shapiro Doppler shift [units of F_tx]
    # See footnote 10 of https://doi.org/10.1103/PhysRevLett.17.933
    shap_dop = -F_tx * shap_del_diff
    return shap_dop
end

raw"""
    Ne(p1, p2, r_s_t0, ds, ΔS) where {S<:Number, U<:Number}

Return the density of ionized electrons [electrons/cm^3] in interplanetary medium.

# Arguments

- `p1::Vector{S}`: signal departure point (transmitter/bounce) [au].
- `p2::Vector{S}`: signal arrival point (bounce/receiver) [au].
- `r_s_t0::Vector{S}`: barycentric position [au] of the Sun at initial time of
    propagation of signal path (bounce time for down-leg; transmit time for up-leg).
- `ds::U`: current distance travelled by ray from emission point [au].
- `ΔS::Real`: total distance between p1 and p2 [au].

!!! reference
    See (Explanatory Supplement to the Astronomical Almanac 2014, p. 323, Sec. 8.7.5,
    Eq. 8.22). ESAA 2014 in turn refers to Muhleman and Anderson (1981). Ostro (1993)
    gives a reference to Anderson (1978), where this model is fitted to Mariner 9 ranging
    data. Reading https://gssc.esa.int/navipedia/index.php/Ionospheric_Delay helped a lot
    to clarify things, especially the 40.3, although they talk about Earth's ionosphere.
    Another valuable source is Standish, E.M., Astron. Astrophys. 233, 252-271 (1990).

# Extended help

The density of ionized electrons in interplanetary medium is given by:
```math
N_e = \frac{A}{r^6} + \frac{ab/\sqrt{a^2\sin^2\beta + b^2\cos^2\beta}}{r^2},
```
where ``r`` is the heliocentric distance expressed in units of the solar radius, ``\beta``
is the solar latitude, and ``A``, ``a``, ``b`` are the solar corona parameters.
"""
function Ne(p1::Vector{S}, p2::Vector{S}, r_s_t0::Vector{S},
            ds::U, ΔS::Real) where {S <: Number, U <: Number}
    # s: linear parametrization of ray path, such that:
    # s = 0 -> point on ray path is at p1
    # s = 1 -> point on ray path is at p2
    s = ds/ΔS
    # Rescale Taylor polynomial
    s_p2_p1 = map(x->s*x, Taylor1.(p2-p1, TaylorSeries.order(s)))
    # Heliocentric position [au] of point on ray path at time t_tdb_jul [Julian days]
    r_vec = Taylor1.(p1, TaylorSeries.order(s)) + s_p2_p1 - Taylor1.(r_s_t0, TaylorSeries.order(s))
    # Heliocentric distance [au] of point on ray path at time t_tdb_jul [Julian days]
    r = sqrt( r_vec[1]^2 + r_vec[2]^2 + r_vec[3]^2 )
    # Compute heliocentric position vector of point on ray path wrt Sun's rotation pole
    # and equator (i.e., heliographic)
    α_p_sun_rad = deg2rad(α_p_sun)
    δ_p_sun_rad = deg2rad(δ_p_sun)
    r_vec_heliographic = inv( pole_rotation(α_p_sun_rad, δ_p_sun_rad) ) * r_vec
    # Compute heliographic (ecliptic) solar latitude (Anderson, 1978) of
    # point on ray path [rad]
    β = asin( r_vec_heliographic[3]/r )
    # Heliocentric distance
    r_sr = r/R_sun
    # First term of Ne
    Ne_t1 = (A_sun/r_sr^6)
    # Second term of Ne
    Ne_t2 = ( (a_sun*b_sun)/sqrt((a_sun*sin(β))^2 + (b_sun*cos(β))^2) )/(r_sr^2)
    # Density of ionized electrons
    Ne_val = Ne_t1 + Ne_t2
    return Ne_val
end

# TODO: @taylorize!
"""
    Ne_path_integral(p1, p2, r_s_t) where {S <: Number}

Return the path integral of the density of ionized electrons in interplanetary
medium ``N_e`` [electrons/cm^2], evaluated with `TaylorIntegration`.

# Arguments

- `p1::Vector{S}`: signal departure point (transmitter/bounce) [au].
- `p2::Vector{S}`: signal arrival point (bounce/receiver) [au].
- `r_s_t0::Vector{S}`: barycentric position [au] of the Sun at initial time of
    propagation of signal path (bounce time for down-leg; transmit time for up-leg).
"""
function Ne_path_integral(p1::Vector{S}, p2::Vector{S},
                          r_s_t0::Vector{S}) where {S <: Number}
    # Total distance between p1 and p2, in centimeters
    ΔS = (100_000au) * norm(p2-p1)
    # Kernel of path integral; distance parameter `s` and total distance `ΔS` is in cm
    function int_kernel(x, params, s)
        return Ne(p1, p2, r_s_t0, s, ΔS)
    end
    # Do path integral

    # Initial condition
    i0 = zero(p1[1])
    iT = Taylor1(i0, 24)
    # Independent variable
    tT = Taylor1(24)
    # Integration
    TaylorIntegration.jetcoeffs!(int_kernel, tT, iT, nothing)
    # Evaluate path integral in total distance ΔS
    return iT(ΔS)
end

raw"""
    corona_delay(p1, p2, r_s_t0, F_tx) where {S <: Number, U <: Real}

Return the time-delay [sec] due to thin plasma of solar corona.

# Arguments

- `p1::Vector{S}`: signal departure point (transmitter/bounce for up/down-link,
    resp.) [au].
- `p2::Vector{S}`: signal arrival point (bounce/receiver for up/down-link,
    resp.) [au].
- `r_s_t0::Vector{S}`: barycentric position [au] of the Sun at initial time of
    propagation of signal path (bounce time for down-leg; transmit time for up-leg).
- `F_tx::U`: transmitter frequency [MHz].

!!! reference
    From https://gssc.esa.int/navipedia/index.php/Ionospheric_Delay it seems that ESAA
    2014 text probably should say that in the formula for ``\Delta\tau_\text{cor}``, the
    expression ``40.3 N_e/f^2`` is adimensional, where ``Ne`` is in electrons/cm^3 and
    ``f`` is in Hz therefore, the integral ``(40.3/f^2)\int Ne \ ds`` is in centimeters,
    where ``ds`` is in cm and the expression ``(40.3/(cf^2))\int Ne \ ds``, with ``c`` in
    cm/sec, is in seconds.

# Extended help

The time-delay due to thin plasma of solar corona is given by:
```math
\Delta\tau_\text{cor} = \frac{40.3}{cf^2}\int_{P_1}^{P_2}N_e \ ds,
```math
where ``c`` is the speed of light [cm/sec], ``f`` is the frequency [Hz], ``N_e`` is the
density of ionized electrons in interplanetary medium [electrons/cm^3], and ``s`` is the
linear distance [cm]. ``N_e`` is computed by [`Ne`](@ref) and integrated via
`TaylorIntegration` in [`Ne_path_integral`](@ref).
"""
function corona_delay(p1::Vector{S}, p2::Vector{S}, r_s_t0::Vector{S},
                      F_tx::U) where {S <: Number, U <: Real}
    # For the time being, we're removing the terms associated with higher-order terms in
    # the variationals (ie, Yarkovsky)
    int_path = Ne_path_integral(p1, p2, r_s_t0) # [electrons/cm^2]
    # Time delay due to solar corona [sec]
    Δτ_corona = 40.3e-6int_path/(c_cm_per_sec*(F_tx)^2)
    return Δτ_corona # seconds
end

"""
    zenith_distance(r_antenna, ρ_vec_ae) where {T <: Number, S <: Number}

**VERY** elementary computation of zenith distance.

# Arguments

- `r_antenna::Vector{T}`: position of antenna at receive/transmit time in
    celestial frame wrt geocenter.
- `ρ_vec_ae::Vector{S}`: slant-range vector from antenna to asteroid.
"""
function zenith_distance(r_antenna::Vector{T},
                         ρ_vec_ae::Vector{S}) where {T <: Number, S <: Number}
    # Magnitude of geocentric antenna position
    norm_r_antenna = sqrt(r_antenna[1]^2 + r_antenna[2]^2 + r_antenna[3]^2)
    # Magnitude of slant-range vector from antenna to asteroid
    norm_ρ_vec_ae = sqrt(ρ_vec_ae[1]^2 + ρ_vec_ae[2]^2 + ρ_vec_ae[3]^2)
    # cos( angle between r_antenna and ρ_vec_ae)
    cos_antenna_slant = dot(r_antenna, ρ_vec_ae)/(norm_r_antenna*norm_ρ_vec_ae)
    # zenith distance
    return acos(cos_antenna_slant)
end

raw"""
    tropo_delay(z)

Return the time-delay [sec] due to Earth's troposphere for radio frequencies.

See also [`zenith_distance`](@doc).

# Arguments

- `z`: zenith distance [rad].

# Extended help

The time-delay for radio frequencies due to Earth's troposphere is given by:
```math
\Delta\tau_\text{tropo} = \frac{7 \ \text{nsec}}{\cos z + \frac{0.0014}{0.045 + \cot z}},
```
where ``z`` is the zenith distance at the antenna. This time-delay oscillates between
0.007``\mu``s and 0.225``\mu``s for ``z`` between 0 and ``\pi/2`` rad.
"""
tropo_delay(z) = (7e-9) / ( cos(z) + 0.0014 / (0.045+cot(z)) ) # seconds

"""
    tropo_delay(r_antenna, ρ_vec_ae) where {T <: Number, S <: Number}

Return the time delay [sec] due to Earth's troposphere for radio frequencies.

The function first computes the zenith distance ``z`` via [`zenith_distance`](@ref)
and then substitutes into the first method of [`tropo_delay`](@ref).

# Arguments

- `r_antenna::Vector{T}`: position of antenna at receive/transmit time in celestial
    frame wrt geocenter.
- `ρ_vec_ae::Vector{S}`: slant-range vector from antenna to asteroid.
"""
function tropo_delay(r_antenna::Vector{T},
                     ρ_vec_ae::Vector{S}) where {T <: Number, S <: Number}
    # zenith distance
    zd = zenith_distance(r_antenna, ρ_vec_ae) # rad
    # Time delay due to Earth's troposphere
    return tropo_delay(zd) # seconds
end

"""
    compute_delay(::ObservatoryMPC, ::DateTime [, ::RadarBuffer]; xva, kwargs...)

Compute the Taylor series expansion of the time-delay [us] observable as seen
by an observatory around an UTC echo reception time. An optional buffer can be
passed to recycle memory.

# Keyword arguments

- `tord::Int`: order of Taylor expansions (default: `10`, or the order of the
    buffer).
- `niter::Int`: number of light-time solution iterations (default: `10`).
- `xve::EarthEph`: Earth ephemeris (default: `earthposvel`).
- `xvs::SunEph`: Sun ephemeris (default: `sunposvel`).
- `xva::AstEph`: asteroid ephemeris.

Without a buffer, all ephemeris must take [et seconds since J2000] and return
[barycentric position in km and velocity in km/sec]. With a buffer, `xvs` and
`xve` must be `DensePropagation2` ephemerides [au, au/day] taking TDB days
since J2000, `xva` must be a tuple with the backward and forward propagations
of the asteroid, and the returned time-delay is stored in the buffer, so it is
overwritten by subsequent calls with the same buffer.

!!! reference
    See https://doi.org/10.1086/116062.

# Extended help

This function allows to compute dopplers via automatic differentiation using
```math
\nu = -f\frac{d\tau}{dt},
```
where ``f`` is the transmitter frequency [MHz] and ``\tau`` is the time-delay at
reception time ``t``. Computed values include corrections due to Earth orientation,
LOD and polar motion.

The above works only with dense `TaylorSolution` ephemerides.
"""
function compute_delay(observatory::ObservatoryMPC{T}, t_r_utc::DateTime; tord::Int = 10,
                       niter::Int = 10, xve::EarthEph = earthposvel, xvs::SunEph = sunposvel,
                       xva::AstEph) where {T <: Real, EarthEph, SunEph, AstEph}

    # Transform receiving time from UTC to TDB seconds since j2000
    et_r_secs_0 = dtutc2et(t_r_utc)
    # Auxiliary to evaluate JT ephemeris
    xva1et0 = xva(et_r_secs_0)[1]
    # et_r_secs_0 as a Taylor polynomial
    et_r_secs = Taylor1([et_r_secs_0,one(et_r_secs_0)].*one(xva1et0), tord)
    utc_r_secs = et_r_secs - tdb_utc(et_r_secs)
    utc_r_days = JD_J2000 + utc_r_secs / daysec
    # Compute geocentric position/velocity of receiving antenna in
    # inertial frame [km, km/sec]
    RV_r = obsposvelECI(observatory, utc_r_days)
    R_r = RV_r[1:3]
    # Earth's barycentric position and velocity at receive time
    r_e_t_r = xve(et_r_secs)[1:3]
    # Receiver barycentric position and velocity at receive time
    r_r_t_r = r_e_t_r + R_r
    # Asteroid barycentric position and velocity at receive time
    r_a_t_r = xva(et_r_secs)[1:3]
    # Sun barycentric position and velocity at receive time
    r_s_t_r = xvs(et_r_secs)[1:3]

    # Down-leg iteration
    # τ_D first approximation
    # See equation (1) of https://doi.org/10.1086/116062
    ρ_vec_r = r_a_t_r - r_r_t_r
    ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)
    # -R_b/c, but delay is wrt asteroid Center (Brozovic et al., 2018)
    τ_D = ρ_r/clightkms # [seconds]
    # Bounce time, new estimate
    # See equation (2) of https://doi.org/10.1086/116062
    et_b_secs = et_r_secs - τ_D

    # Allocate memory for time delays
    Δτ_D = zero(τ_D)            # Total time delay
    Δτ_rel_D = zero(τ_D)        # Shapiro delay
    # Δτ_corona_D = zero(τ_D)   # Delay due to Solar corona
    Δτ_tropo_D = zero(τ_D)      # Delay due to Earth's troposphere

    for i in 1:niter
        # Asteroid barycentric position [au] at bounce time (TDB)
        rv_a_t_b = xva(et_b_secs)
        r_a_t_b = rv_a_t_b[1:3]
        v_a_t_b = rv_a_t_b[4:6]
        # Estimated position of the asteroid's center of mass relative to the recieve point
        # See equation (3) of https://doi.org/10.1086/116062.
        ρ_vec_r = r_a_t_b - r_r_t_r
        # Magnitude of ρ_vec_r
        # See equation (4) of https://doi.org/10.1086/116062.
        ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)

        # Compute down-leg Shapiro delay
        # NOTE: when using PPN, substitute 2 -> 1+γ in expressions for Shapiro delay,
        # Δτ_rel_[D|U]

        # Earth's position at t_r
        e_D_vec  = r_r_t_r - r_s_t_r
        # Heliocentric distance of Earth at t_r
        e_D = sqrt(e_D_vec[1]^2 + e_D_vec[2]^2 + e_D_vec[3]^2)
        # Barycentric position of Sun at estimated bounce time
        r_s_t_b = xvs(et_b_secs)[1:3]
        # Heliocentric position of asteroid at t_b
        p_D_vec  = r_a_t_b - r_s_t_b
        # Heliocentric distance of asteroid at t_b
        p_D = sqrt(p_D_vec[1]^2 + p_D_vec[2]^2 + p_D_vec[3]^2)
        # Signal path distance (down-leg)
        q_D = ρ_r

        # Shapiro correction to time-delay [sec]
        Δτ_rel_D = shapiro_delay(e_D, p_D, q_D)
        # Troposphere correction to time-delay [sec]
        Δτ_tropo_D = tropo_delay(R_r, ρ_vec_r)
        # Solar corona correction to time-delay [sec]
        # Δτ_corona_D = corona_delay(constant_term.(r_a_t_b), r_r_t_r,
        #    r_s_t_r, F_tx, station_code)
        # Total time-delay [sec]
        Δτ_D = Δτ_rel_D # + Δτ_tropo_D #+ Δτ_corona_D

        # New estimate
        p_dot_23 = dot(ρ_vec_r, v_a_t_b)/ρ_r
        # Time delay correction
        Δt_2 = (τ_D - ρ_r/clightkms - Δτ_rel_D)/(1.0-p_dot_23/clightkms)
        # Time delay new estimate
        τ_D = τ_D - Δt_2
        # Bounce time, new estimate
        # See equation (2) of https://doi.org/10.1086/116062
        et_b_secs = et_r_secs - τ_D

    end

    # Asteroid's barycentric position and velocity at bounce time t_b
    rv_a_t_b = xva(et_b_secs)
    r_a_t_b = rv_a_t_b[1:3]
    v_a_t_b = rv_a_t_b[4:6]

    # Up-leg iteration
    # τ_U first estimation
    # See equation (5) of https://doi.org/10.1086/116062
    τ_U = τ_D
    # Transmit time, 1st estimate
    # See equation (6) of https://doi.org/10.1086/116062
    et_t_secs = et_b_secs - τ_U
    # Geocentric position and velocity of transmitting antenna in
    # inertial frame [km, km/sec]
    RV_t = RV_r(et_t_secs-et_r_secs_0)
    R_t = RV_t[1:3]
    V_t = RV_t[4:6]
    # Barycentric position and velocity of the Earth at transmit time
    rv_e_t_t = xve(et_t_secs)
    r_e_t_t = rv_e_t_t[1:3]
    v_e_t_t = rv_e_t_t[4:6]
    # Transmitter barycentric position and velocity of at transmit time
    r_t_t_t = r_e_t_t + R_t
    # Up-leg vector at transmit time
    # See equation (7) of https://doi.org/10.1086/116062
    ρ_vec_t = r_a_t_b - r_t_t_t
    # Magnitude of up-leg vector
    ρ_t = sqrt(ρ_vec_t[1]^2 + ρ_vec_t[2]^2 + ρ_vec_t[3]^2)

    # Allocate memory for time delays
    Δτ_U = zero(τ_U)            # Total time delay
    Δτ_rel_U = zero(τ_U)        # Shapiro delay
    # Δτ_corona_U = zero(τ_U)   # Delay due to Solar corona
    Δτ_tropo_U = zero(τ_U)      # Delay due to Earth's troposphere

    for i in 1:niter
        # Geocentric position and velocity of transmitting antenna in
        # inertial frame [km, km/sec]
        RV_t = RV_r(et_t_secs-et_r_secs_0)
        R_t = RV_t[1:3]
        V_t = RV_t[4:6]
        # Earth's barycentric position and velocity at transmit time
        rv_e_t_t = xve(et_t_secs)
        r_e_t_t = rv_e_t_t[1:3]
        v_e_t_t = rv_e_t_t[4:6]
        # Barycentric position and velocity of the transmitter at the transmit time
        r_t_t_t = r_e_t_t + R_t
        v_t_t_t = v_e_t_t + V_t
        # Up-leg vector and its magnitude at transmit time
        # See equation (7) of https://doi.org/10.1086/116062
        ρ_vec_t = r_a_t_b - r_t_t_t
        ρ_t = sqrt(ρ_vec_t[1]^2 + ρ_vec_t[2]^2 + ρ_vec_t[3]^2)


        # Compute up-leg Shapiro delay

        # Sun barycentric position and velocity [km, km/sec] at transmit time (TDB)
        r_s_t_t = xvs(et_t_secs)[1:3]
        # Heliocentric position of Earth at t_t
        e_U_vec = r_t_t_t - r_s_t_t
        # Heliocentric distance of Earth at t_t
        e_U = sqrt(e_U_vec[1]^2 + e_U_vec[2]^2 + e_U_vec[3]^2)
        # Barycentric position/velocity of Sun at bounce time
        r_s_t_b = xvs(et_b_secs)[1:3]
        # Heliocentric position of asteroid at t_b
        p_U_vec = r_a_t_b - r_s_t_b
        # Heliocentric distance of asteroid at t_b
        p_U = sqrt(p_U_vec[1]^2 + p_U_vec[2]^2 + p_U_vec[3]^2)
        # Signal path distance (up-leg)
        q_U = ρ_t

        # Shapiro correction to time-delay [sec]
        Δτ_rel_U = shapiro_delay(e_U, p_U, q_U)
        # Troposphere correction to time-delay [sec]
        Δτ_tropo_U = tropo_delay(R_t, ρ_vec_t)
        # Delay due to Solar corona [sec]
        # Δτ_corona_U = corona_delay(constant_term.(r_t_t_t), constant_term.(r_a_t_b),
        #    constant_term.(r_s_t_b), F_tx, station_code)
        # Total time delay [sec]
        Δτ_U = Δτ_rel_U # + Δτ_tropo_U #+ Δτ_corona_U

        # New estimate
        p_dot_12 = -dot(ρ_vec_t, v_t_t_t)/ρ_t
        # Time-delay correction
        Δt_1 = (τ_U - ρ_t/clightkms - Δτ_rel_U)/(1.0-p_dot_12/clightkms)
        # Time delay new estimate
        τ_U = τ_U - Δt_1
        # Transmit time, new estimate
        # See equation (6) of https://doi.org/10.1086/116062
        et_t_secs = et_b_secs - τ_U
    end

    # Compute TDB-UTC at transmit time
    # Corrections to TT-TDB from Moyer (2003) / Folkner et al. (2014) due to position
    # of measurement station on Earth are of order 0.01μs
    # Δtt_tdb_station_t = - dot(v_e_t_t, r_t_t_t-r_e_t_t)/clightkms^2
    tdb_utc_t = tdb_utc(et_t_secs) # + Δtt_tdb_station_t

    # Compute TDB-UTC at receive time
    # Corrections to TT-TDB from Moyer (2003) / Folkner et al. (2014) due to position
    # of measurement station on Earth are of order 0.01μs
    # Δtt_tdb_station_r = - dot(v_e_t_r, r_r_t_r-r_e_t_r)/clightkms^2
    tdb_utc_r = tdb_utc(et_r_secs) # + Δtt_tdb_station_r

    # Compute total time delay [UTC seconds]; relativistic delay is already included
    # in τ_D, τ_U. See equation (9) of https://doi.org/10.1086/116062
    τ = (τ_D + τ_U) + (Δτ_tropo_D + Δτ_tropo_U) + (tdb_utc_t - tdb_utc_r)

    # Total signal delay [us]
    return 1e6τ
end

# Taylorized version of the function above
function compute_delay(
        observatory::ObservatoryMPC{T}, t_r_utc::DateTime, buffer::RadarBuffer{U};
        tord::Int = TS.order(buffer.v0[1]), niter::Int = 10,
        xvs::DensePropagation2{T, T}, xve::DensePropagation2{T, T},
        xva::NTuple{2, DensePropagation2{T, U}}
    ) where {T <: Real, U <: Number}
    # Unfold
    tvar, et_r_secs, aux1, aux2, auxh, ρ_r, τ_D, et_b_secs, e_D, p_D, _p_dot_, p_dot,
    τ_ρ, τ_p, one_τ_p, _Δt_, Δt, τ_U, et_t_secs, dt_t, ρ_t, e_U, p_U, τ, τ_us = buffer.v0
    rv_e_t_r, rv_s_t_r, rv_a_t_r, r_r_t_r, ρ_vec_r, e_D_vec, rv_a_t_b, rv_s_t_b,
    p_D_vec, R_t, V_t, rv_e_t_t, rv_s_t_t, r_t_t_t, v_t_t_t, ρ_vec_t, e_U_vec,
    p_U_vec = buffer.v1
    order = TS.order(tvar)
    @assert tord == order "Order of Taylor expansions ($tord) does not match \
        the order of the buffer ($order)"
    # Asteroid barycentric velocity at bounce time
    v_a_t_b = view(rv_a_t_b, 4:6)
    # Transform receiving time from UTC to TDB seconds since J2000
    et_r_secs_0 = dtutc2et(t_r_utc)
    for ord in 0:order
        TS.add!(et_r_secs, tvar, et_r_secs_0, ord)
    end
    # TDB-UTC at receive time
    tdb_utc_r = tdb_utc(et_r_secs)
    # Compute geocentric position/velocity of receiving antenna in
    # inertial frame [km, km/sec]
    utc_r_days = JD_J2000 + (et_r_secs - tdb_utc_r) / daysec
    RV_r = obsposvelECI(observatory, utc_r_days)
    R_r = RV_r[1:3]
    # Earth, Sun and asteroid barycentric positions at receive time
    evaleph!(rv_e_t_r, et_r_secs, xve, auxh)
    evaleph!(rv_s_t_r, et_r_secs, xvs, auxh)
    evaleph!(rv_a_t_r, et_r_secs, xva[1], xva[2], auxh)

    # Down-leg iteration
    for ord in 0:order
        for i in 1:3
            # Receiver barycentric position at receive time
            TS.add!(r_r_t_r[i], rv_e_t_r[i], R_r[i], ord)
            # See equation (1) of https://doi.org/10.1086/116062
            TS.subst!(ρ_vec_r[i], rv_a_t_r[i], r_r_t_r[i], ord)
            # Heliocentric position of Earth at receive time
            TS.subst!(e_D_vec[i], r_r_t_r[i], rv_s_t_r[i], ord)
        end
        euclid3D!(ρ_r, ρ_vec_r, aux1, aux2, ord)
        euclid3D!(e_D, e_D_vec, aux1, aux2, ord)
        # τ_D first approximation [seconds]
        TS.mul!(τ_D, c_kms_m1, ρ_r, ord)
        # Bounce time, first estimate
        # See equation (2) of https://doi.org/10.1086/116062
        TS.subst!(et_b_secs, et_r_secs, τ_D, ord)
    end
    # Allocate memory for time delays
    Δτ_tropo_D = zero(τ_D)      # Delay due to Earth's troposphere
    for _ in 1:niter
        # Asteroid barycentric position and velocity, and Sun barycentric
        # position [km, km/sec] at bounce time
        evaleph!(rv_a_t_b, et_b_secs, xva[1], xva[2], auxh)
        evaleph!(rv_s_t_b, et_b_secs, xvs, auxh)
        for ord in 0:order
            for i in 1:3
                # See equation (3) of https://doi.org/10.1086/116062
                TS.subst!(ρ_vec_r[i], rv_a_t_b[i], r_r_t_r[i], ord)
                # Heliocentric position of asteroid at bounce time
                TS.subst!(p_D_vec[i], rv_a_t_b[i], rv_s_t_b[i], ord)
            end
            # See equation (4) of https://doi.org/10.1086/116062
            euclid3D!(ρ_r, ρ_vec_r, aux1, aux2, ord)
            euclid3D!(p_D, p_D_vec, aux1, aux2, ord)
        end
        # Shapiro and troposphere corrections to time delay [seconds]
        Δτ_rel_D = shapiro_delay(e_D, p_D, ρ_r)
        Δτ_tropo_D = tropo_delay(R_r, ρ_vec_r)
        for ord in 0:order
            # New estimate
            dot3D!(_p_dot_, ρ_vec_r, v_a_t_b, aux1, ord)
            TS.div!(p_dot, _p_dot_, ρ_r, ord)
            # Time delay correction
            TS.mul!(τ_ρ, c_kms_m1, ρ_r, ord)
            TS.mul!(τ_p, c_kms_m1, p_dot, ord)
            TS.subst!(one_τ_p, 1, τ_p, ord)
            TS.subst!(_Δt_, τ_D, τ_ρ, ord)
            TS.subst!(_Δt_, _Δt_, Δτ_rel_D, ord)
            TS.div!(Δt, _Δt_, one_τ_p, ord)
            # Time delay new estimate
            TS.subst!(τ_D, τ_D, Δt, ord)
            # Bounce time, new estimate
            # See equation (2) of https://doi.org/10.1086/116062
            TS.subst!(et_b_secs, et_r_secs, τ_D, ord)
        end
    end

    # Asteroid barycentric position and velocity, and Sun barycentric
    # position [km, km/sec] at bounce time
    evaleph!(rv_a_t_b, et_b_secs, xva[1], xva[2], auxh)
    evaleph!(rv_s_t_b, et_b_secs, xvs, auxh)

    # Up-leg iteration
    for ord in 0:order
        # τ_U first estimate
        # See equation (5) of https://doi.org/10.1086/116062
        TS.identity!(τ_U, τ_D, ord)
        # Transmit time, first estimate
        # See equation (6) of https://doi.org/10.1086/116062
        TS.subst!(et_t_secs, et_b_secs, τ_U, ord)
        TS.subst!(dt_t, et_t_secs, et_r_secs_0, ord)
        # Heliocentric position of asteroid at bounce time
        for i in 1:3
            TS.subst!(p_U_vec[i], rv_a_t_b[i], rv_s_t_b[i], ord)
        end
        euclid3D!(p_U, p_U_vec, aux1, aux2, ord)
    end
    # Allocate memory for time delays
    Δτ_tropo_U = zero(τ_U)      # Delay due to Earth's troposphere
    for _ in 1:niter
        # Geocentric position and velocity of transmitting antenna in
        # inertial frame [km, km/sec]
        for i in 1:3
            taylorcompose!(R_t[i], RV_r[i], dt_t, auxh)
            taylorcompose!(V_t[i], RV_r[i+3], dt_t, auxh)
        end
        # Earth's barycentric position and velocity, and Sun barycentric
        # position [km, km/sec] at transmit time
        evaleph!(rv_e_t_t, et_t_secs, xve, auxh)
        evaleph!(rv_s_t_t, et_t_secs, xvs, auxh)
        for ord in 0:order
            for i in 1:3
                # Barycentric position and velocity of the transmitter at transmit time
                TS.add!(r_t_t_t[i], rv_e_t_t[i], R_t[i], ord)
                TS.add!(v_t_t_t[i], rv_e_t_t[i+3], V_t[i], ord)
                # Up-leg vector at transmit time
                # See equation (7) of https://doi.org/10.1086/116062
                TS.subst!(ρ_vec_t[i], rv_a_t_b[i], r_t_t_t[i], ord)
                # Heliocentric position of Earth at transmit time
                TS.subst!(e_U_vec[i], r_t_t_t[i], rv_s_t_t[i], ord)
            end
            euclid3D!(ρ_t, ρ_vec_t, aux1, aux2, ord)
            euclid3D!(e_U, e_U_vec, aux1, aux2, ord)
        end
        # Shapiro and troposphere corrections to time delay [seconds]
        Δτ_rel_U = shapiro_delay(e_U, p_U, ρ_t)
        Δτ_tropo_U = tropo_delay(R_t, ρ_vec_t)
        for ord in 0:order
            # New estimate (p_dot_12 = -p_dot)
            dot3D!(_p_dot_, ρ_vec_t, v_t_t_t, aux1, ord)
            TS.div!(p_dot, _p_dot_, ρ_t, ord)
            # Time delay correction
            TS.mul!(τ_ρ, c_kms_m1, ρ_t, ord)
            TS.mul!(τ_p, c_kms_m1, p_dot, ord)
            TS.add!(one_τ_p, 1, τ_p, ord)
            TS.subst!(_Δt_, τ_U, τ_ρ, ord)
            TS.subst!(_Δt_, _Δt_, Δτ_rel_U, ord)
            TS.div!(Δt, _Δt_, one_τ_p, ord)
            # Time delay new estimate
            TS.subst!(τ_U, τ_U, Δt, ord)
            # Transmit time, new estimate
            # See equation (6) of https://doi.org/10.1086/116062
            TS.subst!(et_t_secs, et_b_secs, τ_U, ord)
            TS.subst!(dt_t, et_t_secs, et_r_secs_0, ord)
        end
    end

    # TDB-UTC at transmit time
    tdb_utc_t = tdb_utc(et_t_secs)
    for ord in 0:order
        # Total time delay [UTC seconds]; relativistic delay is already included
        # in τ_D, τ_U. See equation (9) of https://doi.org/10.1086/116062
        TS.add!(aux1, τ_D, τ_U, ord)
        TS.add!(aux2, Δτ_tropo_D, Δτ_tropo_U, ord)
        TS.add!(τ, aux1, aux2, ord)
        TS.subst!(auxh, tdb_utc_t, tdb_utc_r, ord)
        TS.add!(τ, τ, auxh, ord)
        # Total signal delay [us]
        TS.mul!(τ_us, 1e6, τ, ord)
    end

    return τ_us
end

"""
    radar_astrometry(::AbstractRadarAstrometry; xva, kwargs...)

Return time-delay [us] and Doppler shift [Hz].

# Keyword arguments

- `tc::Real`: time offset wrt echo reception time, to compute Doppler
    shifts by range differences [seconds]. Offsets will be rounded to
    the nearest millisecond (default: `1.0`).
- `autodiff::Bool`: whether to compute Doppler shift via automatic diff
    of [`compute_delay`](@ref) or not (default: `true`).
- `tord::Int`: order of Taylor expansions (default: `10`).
- `niter::Int`: number of light-time solution iterations (default: `10`).
- `xve::EarthEph`: Earth ephemeris (default: `earthposvel`).
- `xvs::SunEph`: Sun ephemeris (default: `sunposvel`).
- `xva::AstEph`: asteroid ephemeris.

All ephemeris must take  [et seconds since J2000] and return [barycentric
position in km and velocity in km/sec].
"""
radar_astrometry(radar::AbstractRadarAstrometry{T}; kwargs...) where {T <: Real} =
    radar_astrometry(observatory(radar), date(radar), frequency(radar); kwargs...)

function radar_astrometry(observatory::ObservatoryMPC, t_r_utc::DateTime, F_tx::Real;
                          tc::Real = 1.0, autodiff::Bool = true, kwargs...)
    # Compute Doppler shift via automatic differentiation of time-delay
    if autodiff
        # Time delay
        τ = compute_delay(observatory, t_r_utc; kwargs...)
        # Time delay [us] and Doppler shift [Hz]
        return τ[0], -F_tx * τ[1]
    # Compute Doppler shift via numerical differentiation of time-delay
    else
        offset = Dates.Millisecond(1000round(tc/2, digits = 3))
        τe = compute_delay(observatory, t_r_utc + offset; kwargs...)
        τn = compute_delay(observatory, t_r_utc       ; kwargs...)
        τs = compute_delay(observatory, t_r_utc - offset; kwargs...)
        # Time delay [us] and Doppler shift [Hz]
        return τn[0], -F_tx * ((τe[0]-τs[0]) / tc)
    end

end
