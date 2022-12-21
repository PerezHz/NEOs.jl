@doc raw"""
    shapiro_delay(e, p, q)

Returns the relativistic (Shapiro) time-delay in seconds
```math
\Delta\tau[\text{rel}] = \frac{2\mu_\odot}{c^3}\log\left|\frac{d_{E,S} + d_{A,S} + d_{A,E}}{d_{E,S}+d_{A,S}-d_{A, E}}\right|,
```
where ``\mu_\odot = GM_\odot`` is the mass parameter of the sun, and ``d_{E,S}``, ``d_{A,S}``
and ``d_{A,E}`` are the heliocentric distance of the Earth, the asteroid's heliocentric
distance, and the asteroid's  geocentric distance, respectively.

See https://doi.org/10.1103/PhysRevLett.13.789.

# Arguments

- `e`: heliocentric distance of the Earth.
- `p`: asteroid's heliocentric distance.
- `q`: asteroid's geocentric distance.
"""
function shapiro_delay(e, p, q)
    shap = 0.0 # 2μ[1]/(c_au_per_day^2)
    shap_del_days = (2PlanetaryEphemeris.μ[su]/(c_au_per_day^3))*log( (e+p+q+shap)/(e+p-q+shap) ) # days
    return shap_del_days*daysec # seconds
end

@doc raw"""
    shapiro_doppler(e, de, p, dp, q, dq, F_tx)

Returns the Doppler shift (in units of `F_tx`)
```math
\Delta\nu = -\nu\frac{d\Delta\tau}{dt},
```
where ``\nu`` is the frequency and ``\frac{d\Delta\tau}{dt}`` is the differential of the
Shapiro delay. See also [`shapiro_delay`](@ref).

See https://doi.org/10.1103/PhysRevLett.17.933.

# Arguments

- `e`: heliocentric distance of the Earth.
- `de`: differential of `e`.
- `p`: asteroid's heliocentric distance.
- `dp`: differential of `p`.
- `q`: asteroid's geocentric distance.
- `dq`: differential of `q`.
- `F_tx`: transmitter frequency (MHz).
"""
function shapiro_doppler(e, de, p, dp, q, dq, F_tx)
    # shap_del_diff = 2μ[1]*( (de+dp+dq)/(e+p+q) - (de+dp-dq)/(e+p-q) )/(c_au_per_day^3) # (adim.)
    shap_del_diff = (4PlanetaryEphemeris.μ[su]/(c_au_per_day^3))*(  ( dq*(e+p) - q*(de+dp) )/( (e+p)^2 - q^2 )  ) # differential of Shapiro delay (adim.)
    # ν = -F_tx*dτ/dt (units of F_tx)
    # See footnote 10 of https://doi.org/10.1103/PhysRevLett.17.933
    shap_dop = -F_tx*shap_del_diff
    return shap_dop # (units of F_tx)
end

@doc raw"""
    Ne(p1::Vector{S}, p2::Vector{S}, r_s_t0::Vector{S}, ds::U, ΔS::Real) where {S<:Number, U<:Number}

Returns the density of ionized electrons (in electrons/cm^3) in interplanetary medium
```math
N_e = \frac{A}{r^6} + \frac{ab/\sqrt{a^2\sin^2\beta + b^2\cos^2\beta}}{r^2},
```
where ``r`` is the heliocentric distance expressed in units of the solar radius, ``\beta`` is
the solar latitude, and ``A``, ``a``, ``b`` are the solar corona parameters.

See (Explanatory Supplement to the Astronomical Almanac 2014, p. 323, Sec. 8.7.5, Eq. 8.22).
ESAA 2014 in turn refers to Muhleman and Anderson (1981). Ostro (1993) gives a reference to
Anderson (1978), where this model is fitted to Mariner 9 ranging data. Reading https://gssc.esa.int/navipedia/index.php/Ionospheric_Delay
helped a lot to clarify things, especially the 40.3, although they talk about Earth's
ionosphere. Another valuable source is Standish, E.M., Astron. Astrophys. 233, 252-271 (1990).

# Arguments
- `p1`: signal departure point (transmitter/bounce) (au).
- `p2`: signal arrival point (bounce/receiver) (au).
- `r_s_t0`: Barycentric position (au) of Sun at initial time of propagation of signal path (bounce time for down-leg; transmit time for up-leg).
- `ds`: current distance travelled by ray from emission point (au).
- `ΔS`: total distance between p1 and p2 (au).
"""
function Ne(p1::Vector{S}, p2::Vector{S}, r_s_t0::Vector{S}, ds::U, ΔS::Real) where {S<:Number, U<:Number}
    # s: linear parametrization of ray path, such that
    # s = 0 -> point on ray path is at p1; s = 1 -> point on ray path is at p2
    s = ds/ΔS
    # Rescale Taylor polynomial
    s_p2_p1 = map(x->s*x, Taylor1.(p2-p1, s.order))
    # Heliocentric position (au) of point on ray path at time t_tdb_jul (Julian days)
    r_vec = Taylor1.(p1, s.order) + s_p2_p1 - Taylor1.(r_s_t0, s.order)
    # Heliocentric distance (au) of point on ray path at time t_tdb_jul (Julian days)
    r = sqrt( r_vec[1]^2 + r_vec[2]^2 + r_vec[3]^2 )
    # Compute heliocentric position vector of point on ray path wrt Sun's rotation pole
    # and equator (i.e., heliographic)
    α_p_sun_rad = deg2rad(α_p_sun)
    δ_p_sun_rad = deg2rad(δ_p_sun)
    r_vec_heliographic = inv( pole_rotation(α_p_sun_rad, δ_p_sun_rad) )*r_vec
    # Compute heliographic solar latitude (Anderson, 1978) of point on ray path
    β = asin( r_vec_heliographic[3]/r ) # Ecliptic solar latitude of point on ray path (rad)
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
@doc raw"""
    Ne_path_integral(p1::Vector{S}, p2::Vector{S}, r_s_t0::Vector{S}) where {S<:Number}

Returns the path integral of the density of ionized electrons in interplanetary medium ``N_e``
in units of electrons/cm^2, evaluated with `TaylorIntegration`.

# Arguments

- `p1`: signal departure point (transmitter/bounce) (au).
- `p2`: signal arrival point (bounce/receiver) (au).
- `r_s_t0`: Barycentric position (au) of Sun at initial time of propagation of signal path (bounce time for down-leg; transmit time for up-leg).
"""
function Ne_path_integral(p1::Vector{S}, p2::Vector{S}, r_s_t0::Vector{S}) where {S<:Number}
    # Total distance between p1 and p2, in centimeters
    ΔS = (100_000au)*norm(p2-p1)
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

@doc raw"""
    corona_delay(p1::Vector{S}, p2::Vector{S}, r_s_t0::Vector{S}, F_tx::U, station_code::Int) where {S<:Number, U<:Real}

Returns the time-delay (in sec) due to thin plasma of solar corona
```math
\Delta\tau_\text{cor} = \frac{40.3}{cf^2}\int_{P_1}^{P_2}N_e \ ds,
```math
where ``c`` is the speed of light (cm/sec), ``f`` is the frequency (Hz), ``N_e`` is the density
of ionized electrons in interplanetary medium (electrons/cm^3), and ``s`` is the linear
distance (cm). ``N_e`` is computed by [`Ne`](@ref) and integrated via `TaylorIntegration` in
[`Ne_path_integral`](@ref).

From [Ionospheric Delay](https://gssc.esa.int/navipedia/index.php/Ionospheric_Delay) it seems
that ESAA 2014 text probably should say that in the formula for ``\Delta\tau_\text{cor}``,
the expression ``40.3 N_e/f^2`` is adimensional, where ``Ne`` is in electrons/cm^3 and ``f`` is
in Hz therefore, the integral ``(40.3/f^2)\int Ne \ ds`` is in centimeters, where ``ds`` is in
cm and the expression ``(40.3/(cf^2))\int Ne \ ds``, with ``c`` in cm/sec, is in seconds.

# Arguments

- `p1`: signal departure point (transmitter/bounce for up/down-link, resp.) (au).
- `p2`: signal arrival point (bounce/receiver for up/down-link, resp.) (au).
- `t_tdb_jd1`, `t_tdb_jd2`: Two-part Julian date (TDB) of signal path (bounce time for downlink; transmit time for uplink).
- `r_s_t0`: Barycentric position (au) of Sun at initial time of propagation of signal path (bounce time for down-leg; transmit time for up-leg).
- `F_tx`: transmitter frequency (MHz).
- `station_code`: observing station identifier (MPC nomenclature).
"""
function corona_delay(p1::Vector{S}, p2::Vector{S}, r_s_t0::Vector{S}, F_tx::U, station_code::Int) where {S<:Number, U<:Real}
    # For the time being, we're removing the terms associated with higher-order terms in the
    # variationals (ie, Yarkovsky)
    int_path = Ne_path_integral(p1, p2, r_s_t0) # (electrons/cm^2)
    # Time delay due to solar corona
    Δτ_corona = 40.3e-6int_path/(c_cm_per_sec*(F_tx)^2) # seconds
    return Δτ_corona # seconds
end

@doc raw"""
    zenith_distance(r_antenna::Vector{T}, ρ_vec_ae::Vector{S}) where {T<:Number, S<:Number}

**VERY** elementary computation of zenith distance.

# Arguments

- `r_antenna`: position of antenna at receive/transmit time in celestial frame wrt geocenter.
- `ρ_vec_ae`: slant-range vector from antenna to asteroid.
"""
function zenith_distance(r_antenna::Vector{T}, ρ_vec_ae::Vector{S}) where {T<:Number, S<:Number}
    # Magnitude of geocentric antenna position
    norm_r_antenna = sqrt(r_antenna[1]^2 + r_antenna[2]^2 + r_antenna[3]^2)
    # Magnitude of slant-range vector from antenna to asteroid
    norm_ρ_vec_ae = sqrt(ρ_vec_ae[1]^2 + ρ_vec_ae[2]^2 + ρ_vec_ae[3]^2)
    # cos( angle between r_antenna and ρ_vec_ae)
    cos_antenna_slant = dot(r_antenna, ρ_vec_ae)/(norm_r_antenna*norm_ρ_vec_ae)
    # zenith distance
    return acos(cos_antenna_slant)
end

@doc raw"""
    tropo_delay(z)

Returns the time delay (in sec) due to Earth's troposphere for radio frequencies
```math
\Delta\tau_\text{tropo} = \frac{7 \ \text{nsec}}{\cos z + \frac{0.0014}{0.045 + \cot z}},
```
where ``z`` is the zenith distance at the antenna. This time delay oscillates between
0.007``\mu``s and 0.225``\mu``s for ``z`` between 0 and ``\pi/2`` rad.

See also [`zenith_distance`](@doc).

# Arguments

- `z`: zenith distance (radians).
"""
tropo_delay(z) = (7e-9)/( cos(z) + 0.0014/(0.045+cot(z)) ) # seconds

@doc raw"""
    tropo_delay(r_antenna::Vector{T}, ρ_vec_ae::Vector{S}) where {T<:Number, S<:Number}

Returns the time delay (in sec) due to Earth's troposphere for radio frequencies. The function
first computes the zenith distance ``z`` via [`zenith_distance`](@ref) and then substitutes
into the fisrt method of [`tropo_delay`](@ref).

# Arguments

- `r_antenna`: position of antenna at receive/transmit time in celestial frame wrt geocenter.
- `ρ_vec_ae`: slant-range vector from antenna to asteroid.
"""
function tropo_delay(r_antenna::Vector{T}, ρ_vec_ae::Vector{S}) where {T<:Number, S<:Number}
    # zenith distance
    zd = zenith_distance(r_antenna, ρ_vec_ae) # rad
    # Time delay due to Earth's troposphere
    return tropo_delay(zd) # seconds
end

@doc raw"""
    tdb_utc(et::T) where {T<:Number}

Auxiliary function to compute (TDB-UTC)
```math
\begin{align*}
TDB-UTC & = (TDB-TAI) + (TAI-UTC) \\
        & = (TDB-TT) + (TT-TAI) + (TAI-UTC) \\
        & = (TDB-TT) + 32.184 s + ΔAT,
\end{align*}
```
where TDB is the Solar System barycentric ephemeris time, TT is the Terrestrial time,
TAI is the International Atomic Time, and UTC is the Coordinated Universal Time.

This function is useful to convert TDB to UTC via UTC + (TDB-UTC) and viceversa. This function
does not include correction due to position of measurement station ``v_E(r_S.r_E)/c^2``
(Folkner et al. 2014; Moyer, 2003).

# Arguments

- `et::T`: TDB seconds since J2000.0.
"""
function tdb_utc(et::T) where {T<:Number}
    # TT-TDB
    tt_tdb_et = ttmtdb(et)
    # TT-TAI
    tt_tai = 32.184

    et_00 = constant_term(constant_term(et))
    # Used only to determine ΔAT; no high-precision needed
    utc_secs = et_00 - deltet(et_00, "ET")
    # ΔAT
    jd_utc = JD_J2000 + utc_secs/daysec
    tai_utc = get_ΔAT(jd_utc)
    # TDB-UTC = (TDB-TT) + (TT-TAI) + (TAI-UTC) = (TDB-TT) + 32.184 s + ΔAT
    return (tt_tai + tai_utc) - tt_tdb_et
end

# TODO: add tdb_utc(utc) method!!!
# strategy: given UTC, do UTC + (TT-TAI) + (TAI-UTC) to get TT
# then, use ttmtdb(et) function iteratively (Newton) to compute TDB

# function dtutc2et(t_utc::DateTime)
#     tt_tai = 32.184
#     jd_utc = datetime2julian(t_utc)
#     fd_utc = (jd_utc+0.5) - floor(jd_utc+0.5)
#     j, tai_utc = iauDat(year(t_utc), month(t_utc), day(t_utc), fd_utc)
#     return et
# end

@doc raw"""
    delay(station_code::Int, t_r_utc::DateTime, t_offset::Real, niter::Int=10; eo::Bool=true,
          xve=earth_pv, xvs=sun_pv, xva=apophis_pv_197)

Compute radar-astrometric round-trip time (``\mu``s) for an asteroid at UTC instant `t_r_utc`
from tracking station with code `station_code`.

See https://doi.org/10.1086/116062.

# Arguments

- `station_code`: observing station identifier (MPC nomenclature).
- `t_r_utc`: UTC time of echo reception (DateTime).
- `t_offset`: time offset wrt echo reception time, to compute Doppler shifts by range differences (seconds).
- `niter`: number of light-time solution iterations.
- `xve`: Earth ephemeris wich takes TDB seconds since J2000 as input and returns Earth barycentric position in km and velocity in km/second.
- `xvs`: Sun ephemeris wich takes TDB seconds since J2000 as input and returns Sun barycentric position in km and velocity in km/second.
- `xva`: asteroid ephemeris wich takes TDB seconds since J2000 as input and returns asteroid barycentric position in km and velocity in km/second.
"""
function delay(station_code::Int, t_r_utc::DateTime, t_offset::Real,
        niter::Int=10; eo::Bool=true, xve=earth_pv, xvs=sun_pv,
        xva=apophis_pv_197)
    # Transform receiving time from UTC to TDB seconds since j2000
    et_r_secs = str2et(string(t_r_utc)) + t_offset
    # Compute geocentric position/velocity of receiving antenna in inertial frame (au, au/day)
    R_r, V_r = obs_pv_ECI(station_code, et_r_secs, eo=eo)
    # Earth's barycentric position and velocity at receive time
    rv_e_t_r = xve(et_r_secs)
    r_e_t_r = rv_e_t_r[1:3]
    v_e_t_r = rv_e_t_r[4:6]
    # Receiver barycentric position and velocity at receive time
    r_r_t_r = r_e_t_r + R_r
    # Asteroid barycentric position and velocity at receive time
    rv_a_t_r = xva(et_r_secs)
    r_a_t_r = rv_a_t_r[1:3]
    # Sun barycentric position and velocity at receive time
    rv_s_t_r = xvs(et_r_secs)
    r_s_t_r = rv_s_t_r[1:3]

    # Down-leg iteration
    # τ_D first approximation
    # See equation (1) of https://doi.org/10.1086/116062
    ρ_vec_r = r_a_t_r - r_r_t_r
    ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)
    # -R_b/c, but delay is wrt asteroid Center (Brozovic et al., 2018)
    τ_D = ρ_r/clightkms # (seconds)
    # Bounce time, new estimate
    # See equation (2) of https://doi.org/10.1086/116062
    et_b_secs = et_r_secs - τ_D

    # Allocate memmory for time delays
    Δτ_D = zero(τ_D)            # Total time delay
    Δτ_rel_D = zero(τ_D)        # Shapiro delay
    # Δτ_corona_D = zero(τ_D)   # Delay due to Solar corona
    Δτ_tropo_D = zero(τ_D)      # Delay due to Earth's troposphere

    for i in 1:niter
        # Asteroid barycentric position (in au) and velocity (in au/day) at bounce time (TDB)
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
        # Barycentric position and velocity of Sun at estimated bounce time
        rv_s_t_b = xvs(et_b_secs)
        r_s_t_b = rv_s_t_b[1:3]
        # Heliocentric position of asteroid at t_b
        p_D_vec  = r_a_t_b - r_s_t_b
        # Heliocentric distance of asteroid at t_b
        p_D = sqrt(p_D_vec[1]^2 + p_D_vec[2]^2 + p_D_vec[3]^2)
        # Signal path distance (down-leg)
        q_D = ρ_r

        # Shapiro correction to time-delay
        Δτ_rel_D = shapiro_delay(e_D, p_D, q_D) # seconds
        # Troposphere correction to time-delay
        Δτ_tropo_D = tropo_delay(R_r, ρ_vec_r)  # seconds
        # Δτ_corona_D = corona_delay(constant_term.(r_a_t_b), r_r_t_r, r_s_t_r, F_tx, station_code) # seconds
        # Total time-delay
        Δτ_D = Δτ_rel_D # + Δτ_tropo_D #+ Δτ_corona_D # seconds

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
    # Geocentric position and velocity of transmitting antenna in inertial frame (au, au/day)
    R_t, V_t = obs_pv_ECI(station_code, et_t_secs, eo=eo)
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

    # Allocate memmory for time delays
    Δτ_U = zero(τ_U)            # Total time delay
    Δτ_rel_U = zero(τ_U)        # Shapiro delay
    # Δτ_corona_U = zero(τ_U)   # Delay due to Solar corona
    Δτ_tropo_U = zero(τ_U)      # Delay due to Earth's troposphere

    for i in 1:niter
        # Geocentric position and velocity of transmitting antenna in inertial frame (au, au/day)
        # TODO: remove `constant_term` to take into account dependency of R_t, V_t wrt initial
        # conditions variations via et_t_secs
        R_t, V_t = obs_pv_ECI(station_code, et_t_secs, eo=eo)
        # Earth's barycentric position and velocity at the transmit time
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

        # Sun barycentric position and velocity (in au, au/day) at transmit time (TDB)
        rv_s_t_t = xvs(et_t_secs)
        r_s_t_t = rv_s_t_t[1:3]
        # Heliocentric position of Earth at t_t
        e_U_vec = r_t_t_t - r_s_t_t
        # Heliocentric distance of Earth at t_t
        e_U = sqrt(e_U_vec[1]^2 + e_U_vec[2]^2 + e_U_vec[3]^2)
        # Barycentric position/velocity of Sun at bounce time
        rv_s_t_b = xvs(et_b_secs)
        r_s_t_b = rv_s_t_b[1:3]
        # Heliocentric position of asteroid at t_b
        p_U_vec = r_a_t_b - r_s_t_b
        # Heliocentric distance of asteroid at t_b
        p_U = sqrt(p_U_vec[1]^2 + p_U_vec[2]^2 + p_U_vec[3]^2)

        q_U_vec = r_a_t_b - r_e_t_t
        # Signal path distance (up-leg)
        q_U = ρ_t
        # Shapiro correction to time-delay
        Δτ_rel_U = shapiro_delay(e_U, p_U, q_U) # seconds
        # Troposphere correction to time-delay
        Δτ_tropo_U = tropo_delay(R_t, ρ_vec_t)  # seconds
        # Delay due to Solar corona
        # Δτ_corona_U = corona_delay(constant_term.(r_t_t_t), constant_term.(r_a_t_b), constant_term.(r_s_t_b), F_tx, station_code) # seconds
        # Total time delay
        Δτ_U = Δτ_rel_U # + Δτ_tropo_U #+ Δτ_corona_U # seconds

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
    # of measurement station on Earth  are of order 0.01μs
    # Δtt_tdb_station_r = - dot(v_e_t_r, r_r_t_r-r_e_t_r)/clightkms^2
    tdb_utc_r = tdb_utc(et_r_secs) # + Δtt_tdb_station_r

    # Compute total time delay (UTC seconds); relativistic delay is already included in τ_D, τ_U
    # See equation (9) of https://doi.org/10.1086/116062
    τ = (τ_D + τ_U) + (Δτ_tropo_D + Δτ_tropo_U) + (tdb_utc_t - tdb_utc_r) # seconds

    # Total signal delay (μs)
    return 1e6τ
end

@doc raw"""
    delay(station_code::Int, t_r_utc::DateTime, niter::Int=10; eo::Bool=true,
          xve::TaylorInterpolant=earth_pv, xvs::TaylorInterpolant=sun_pv,
          xva::TaylorInterpolant=apophis_pv_197, tord::Int=xva.x[1].order)

Compute Taylor series expansion of time-delay observable around echo reception time. This
allows to compute dopplers via automatic differentiation using
```math
\nu = -f\frac{d\tau}{dt},
```
where ``f`` is the transmitter frequency (MHz) and ``\tau`` is the time-delay at reception
time ``t``.

See https://doi.org/10.1086/116062.

**Note:** Works only with `TaylorInterpolant` ephemeris. See
[`PlanetaryEphemeris.TaylorInterpolant`](@ref).

# Arguments

- `station_code`: observing station identifier (MPC nomenclature).
- `t_r_utc`: UTC time of echo reception (DateTime).
- `niter`: number of light-time solution iterations.
- `eo`: wheter to use `EarthOrientation` or not.
- `xve`: Earth ephemeris wich takes TDB seconds since J2000 as input and returns Earth barycentric position in km and velocity in km/second.
- `xvs`: Sun ephemeris wich takes TDB seconds since J2000 as input and returns Sun barycentric position in km and velocity in km/second.
- `xva`: asteroid ephemeris wich takes TDB seconds since J2000 as input and returns asteroid barycentric position in km and velocity in km/second.
- `tord`: order of Taylor expansions.
"""
function delay(station_code::Int, t_r_utc::DateTime,
        niter::Int=10; eo::Bool=true, xve::TaylorInterpolant=earth_pv,
        xvs::TaylorInterpolant=sun_pv, xva::TaylorInterpolant=apophis_pv_197,
        tord::Int=xva.x[1].order)

    # Auxiliary to evaluate JT ephemeris
    q1 = xva.x[1]
    # Transform receiving time from UTC to TDB seconds since j2000
    et_r_secs_0 = str2et(string(t_r_utc))
    # et_r_secs_0 a a Taylor polynomial
    et_r_secs = Taylor1([et_r_secs_0,1.0].*one(q1[0]), tord)
    # Compute geocentric position/velocity of receiving antenna in inertial frame (au, au/day)
    R_r, _ = obs_pv_ECI(station_code, et_r_secs, eo=eo)
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
    τ_D = ρ_r/clightkms # (seconds)
    # Bounce time, new estimate
    # See equation (2) of https://doi.org/10.1086/116062
    et_b_secs = et_r_secs - τ_D

    # Allocate memmory for time delays
    Δτ_D = zero(τ_D)            # Total time delay
    Δτ_rel_D = zero(τ_D)        # Shapiro delay
    # Δτ_corona_D = zero(τ_D)   # Delay due to Solar corona
    Δτ_tropo_D = zero(τ_D)      # Delay due to Earth's troposphere

    for i in 1:niter
        # Asteroid barycentric position (in au) at bounce time (TDB)
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

        # Shapiro correction to time-delay
        Δτ_rel_D = shapiro_delay(e_D, p_D, q_D) # seconds
        # Troposphere correction to time-delay
        Δτ_tropo_D = tropo_delay(R_r, ρ_vec_r)  # seconds
        # Solar corona correction to time-delay
        # Δτ_corona_D = corona_delay(constant_term.(r_a_t_b), r_r_t_r, r_s_t_r, F_tx, station_code) # seconds
        # Total time-delay
        Δτ_D = Δτ_rel_D # + Δτ_tropo_D #+ Δτ_corona_D # seconds

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
    # Geocentric position and velocity of transmitting antenna in inertial frame (au, au/day)
    R_t, V_t = obs_pv_ECI(station_code, et_t_secs, eo=eo)
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

    # Allocate memmory for time delays
    Δτ_U = zero(τ_U)            # Total time delay
    Δτ_rel_U = zero(τ_U)        # Shapiro delay
    # Δτ_corona_U = zero(τ_U)   # Delay due to Solar corona
    Δτ_tropo_U = zero(τ_U)      # Delay due to Earth's troposphere

    for i in 1:niter
        # Geocentric position and velocity of transmitting antenna in inertial frame (au, au/day)
        # TODO: remove `constant_term` to take into account dependency of R_t, V_t wrt initial
        # conditions variations via et_t_secs
        R_t, V_t = obs_pv_ECI(station_code, et_t_secs, eo=eo)
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

        # Sun barycentric position and velocity (in au, au/day) at transmit time (TDB)
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
        q_U_vec = r_a_t_b - r_e_t_t
        # Signal path distance (up-leg)
        q_U = ρ_t

        # Shapiro correction to time-delay
        Δτ_rel_U = shapiro_delay(e_U, p_U, q_U) # seconds
        # Troposphere correction to time-delay
        Δτ_tropo_U = tropo_delay(R_t, ρ_vec_t)  # seconds
        # Delay due to Solar corona
        # Δτ_corona_U = corona_delay(constant_term.(r_t_t_t), constant_term.(r_a_t_b), constant_term.(r_s_t_b), F_tx, station_code) # seconds
        # Total time delay
        Δτ_U = Δτ_rel_U # + Δτ_tropo_U #+ Δτ_corona_U # seconds

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
    # of measurement station on Earth  are of order 0.01μs
    # Δtt_tdb_station_r = - dot(v_e_t_r, r_r_t_r-r_e_t_r)/clightkms^2
    tdb_utc_r = tdb_utc(et_r_secs) # + Δtt_tdb_station_r

    # Compute total time delay (UTC seconds); relativistic delay is already included in τ_D, τ_U
    # See equation (9) of https://doi.org/10.1086/116062
    τ = (τ_D + τ_U) + (Δτ_tropo_D + Δτ_tropo_U) + (tdb_utc_t - tdb_utc_r) # seconds

    # Total signal delay (μs)
    return 1e6τ
end

@doc raw"""
    delay_doppler(station_code::Int, t_r_utc::DateTime, F_tx::Real, niter::Int=10;
                  eo::Bool=true, tc::Real=1.0, xve=earth_pv, xvs=sun_pv, xva=apophis_pv_197,
                  autodiff::Bool=true, tord::Int=10)

Returns time-delay and Doppler shift.

# Arguments

- `station_code`: observing station identifier (MPC nomenclature).
- `t_r_utc`: UTC time of echo reception (DateTime).
- `F_tx`: transmitter frequency (MHz).
- `niter`: number of light-time solution iterations.
- `eo`: wheter to use `EarthOrientation` or not.
- `tc`: time offset wrt echo reception time, to compute Doppler shifts by range differences (seconds).
- `xve`: Earth ephemeris wich takes TDB seconds since J2000 as input and returns Earth barycentric position in km and velocity in km/second.
- `xvs`: Sun ephemeris wich takes TDB seconds since J2000 as input and returns Sun barycentric position in km and velocity in km/second.
- `xva`: asteroid ephemeris wich takes TDB seconds since J2000 as input and returns asteroid barycentric position in km and velocity in km/second.
- `autodiff`: wheter to use the automatic differentiation method of [`delay`](@ref) or not.
- `tord`: order of Taylor expansions.
"""
function delay_doppler(station_code::Int, t_r_utc::DateTime, F_tx::Real,
        niter::Int=10; eo::Bool=true, tc::Real=1.0, xve=earth_pv, xvs=sun_pv,
        xva=apophis_pv_197, autodiff::Bool=true, tord::Int=10)

    # Automatic differentiation method of delay
    if autodiff
        # Time delay
        τ = delay(station_code, t_r_utc, niter, eo=eo, xve=xve, xvs=xvs, xva=xva, tord=tord)
        # Time delay, Doppler shift
        return τ[0], -F_tx*τ[1]
    # No automatic differentiation method of delay
    else
        τe = delay(station_code, t_r_utc,  tc/2, niter, eo=eo, xve=xve, xvs=xvs, xva=xva)
        τn = delay(station_code, t_r_utc,   0.0, niter, eo=eo, xve=xve, xvs=xvs, xva=xva)
        τs = delay(station_code, t_r_utc, -tc/2, niter, eo=eo, xve=xve, xvs=xvs, xva=xva)
        # Time delay, Doppler shift
        return τn, -F_tx*((τe-τs)/tc)
    end

end

@doc raw"""
    delay_doppler(astradarfile::String, niter::Int=10; eo::Bool=true, tc::Real=1.0,
                  xve=earth_pv, xvs=sun_pv, xva=apophis_pv_197, autodiff::Bool=true,
                  tord::Int=10)

Returns a `DataFrame` with the time delays and Doppler shifts associated with the radar
measurements of file `astradarfile`.

See also [`process_radar_data_jpl`](@ref).

# Arguments

- `astradarfile::String`: radar data file.
- `niter::Int`: number of light-time solution iterations.
- `eo`: wheter to use `EarthOrientation` or not.
- `tc`: time offset wrt echo reception time, to compute Doppler shifts by range differences (seconds).
- `xve`: Earth ephemeris wich takes TDB seconds since J2000 as input and returns Earth barycentric position in km and velocity in km/second.
- `xvs`: Sun ephemeris wich takes TDB seconds since J2000 as input and returns Sun barycentric position in km and velocity in km/second.
- `xva`: asteroid ephemeris wich takes TDB seconds since J2000 as input and returns asteroid barycentric position in km and velocity in km/second.
- `autodiff`: wheter to use the automatic differentiation method of [`delay`](@ref) or not.
- `tord`: order of Taylor expansions.
"""
function delay_doppler(astradarfile::String,
        niter::Int=10; eo::Bool=true, tc::Real=1.0, xve=earth_pv, xvs=sun_pv,
        xva=apophis_pv_197, autodiff::Bool=true, tord::Int=10)

    # Read radar measurements
    astradardata = process_radar_data_jpl(astradarfile)
    #
    et1 = str2et(string(astradardata[1].utcepoch))
    #
    a1_et1 = xva(et1)[1]
    #
    S = typeof(a1_et1)

    # Time delays
    vdelay = Array{S}(undef, length(astradardata))
    # Doppler shifts
    vdoppler = Array{S}(undef, length(astradardata))

    # Iterate over the measurements
    for i in eachindex(astradardata)
        # Compute time delay and doppler shift
        vdelay[i], vdoppler[i] = delay_doppler(
            astradardata[i].rcvr,
            astradardata[i].utcepoch,
            astradardata[i].freq,
            niter,
            eo = eo,
            tc = tc,
            xve = xve,
            xvs = xvs,
            xva = xva,
            autodiff = autodiff,
            tord = tord
        )
    end
    # Rows with time delays
    delay_index = map(x->x.delay_units=="us", astradardata)
    # Rows with Doppler shifts
    doppler_index = map(x->x.doppler_units=="Hz", astradardata)
    # Time delays / Doppler shifts table
    radobs_t = DataFrame(
        (
            dt_utc_obs=utcepoch.(astradardata),        # UTD time
            τ_obs=delay.(astradardata),                # Observed time delay
            ν_obs=doppler.(astradardata),              # Observed Doppler shift
            τ_comp=vdelay,                             # Computed time delay
            ν_comp=vdoppler,                           # Computed Doppler shift
            σ_τ=delay_sigma.(astradardata),            # Observed time delay uncertainty
            σ_ν=doppler_sigma.(astradardata),          # Observed Doppler shift uncertainty
            τ_units=delay_units.(astradardata),        # Time delay units
            ν_units=doppler_units.(astradardata),      # Doppler shift units
            freq=freq.(astradardata),                  # Frequency
            rcvr=rcvr.(astradardata),                  # ID of reciever antenna
            xmit=xmit.(astradardata),                  # ID of emission antenna
            bouncepoint=bouncepoint.(astradardata),    # Bounce point
            delay_index=delay_index,                   # Rows with time delays
            doppler_index=doppler_index                # Rows with Doppler shifts
        )
    )
    # Return time delays / Doppler shifts table
    return radobs_t
end

@doc raw"""
    delay_doppler_yeomansetal92(station_code::Int, t_r_utc::DateTime, F_tx::Real,
                                niter::Int=10; eo::Bool=true, xve=earth_pv, xvs=sun_pv,
                                xva=apophis_pv_197)

Compute round-trip time (in ``\mu``s) and Doppler shift (Hz) radar astrometry for an asteroid at
UTC instant `t_r_utc` from tracking station with code `station_code` from Earth, Sun and
asteroid ephemerides. Dopplers are computed following https://doi.org/10.1086/116062.

# Arguments

- `station_code`: observing station identifier (MPC nomenclature).
- `et_r_secs`: time of echo reception (TDB seconds since J2000.0 TDB).
- `F_tx`: transmitter frequency (MHz).
- `niter`: number of light-time solution iterations.
- `xve`: Earth ephemeris wich takes et seconds since J2000 as input and returns Earth barycentric position in au and velocity in au/day.
- `xvs`: Sun ephemeris wich takes et seconds since J2000 as input and returns Sun barycentric position in au and velocity in au/day.
- `xva`: asteroid ephemeris wich takes et seconds since J2000 as input and returns asteroid barycentric position in au and velocity in au/day.
"""
function delay_doppler_yeomansetal92(station_code::Int, t_r_utc::DateTime,
        F_tx::Real, niter::Int=10; eo::Bool=true, xve=earth_pv, xvs=sun_pv,
        xva=apophis_pv_197)
    # Transform receiving time from UTC to TDB seconds since j2000
    et_r_secs = str2et(string(t_r_utc))
    # Compute geocentric position/velocity of receiving antenna in inertial frame (au, au/day)
    R_r, V_r = obs_pv_ECI(station_code, et_r_secs, eo=eo)
    # Earth's barycentric position and velocity at receive time
    rv_e_t_r = xve(et_r_secs)
    r_e_t_r = rv_e_t_r[1:3]
    v_e_t_r = rv_e_t_r[4:6]
    # Receiver barycentric position and velocity at receive time
    r_r_t_r = r_e_t_r + R_r
    v_r_t_r = v_e_t_r + V_r
    # Asteroid barycentric position and velocity at receive time
    rv_a_t_r = xva(et_r_secs)
    r_a_t_r = rv_a_t_r[1:3]
    v_a_t_r = rv_a_t_r[4:6]
    # Sun barycentric position and velocity at receive time
    rv_s_t_r = xvs(et_r_secs)
    r_s_t_r = rv_s_t_r[1:3]
    v_s_t_r = rv_s_t_r[4:6]

    # Down-leg iteration
    # τ_D first approximation
    # See equation (1) of https://doi.org/10.1086/116062
    ρ_vec_r = r_a_t_r - r_r_t_r
    ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)
    # -R_b/c, but delay is wrt asteroid center (Brozovic et al., 2018)
    τ_D = ρ_r/clightkms # (seconds)
    # Bounce time, 1st estimate
    # See equation (2) of https://doi.org/10.1086/116062
    et_b_secs = et_r_secs - τ_D
    # Asteroid barycentric position (in au) at bounce time (TDB)
    rv_a_t_b = xva(et_b_secs)
    r_a_t_b = rv_a_t_b[1:3]
    v_a_t_b = rv_a_t_b[4:6]

    # Allocate memmory for time delays
    Δτ_D = zero(τ_D)            # Total time delay
    Δτ_rel_D = zero(τ_D)        # Shapiro delay
    # Δτ_corona_D = zero(τ_D)   # Delay due to Solar corona
    Δτ_tropo_D = zero(τ_D)      # Delay due to Earth's troposphere

    for i in 1:niter
        # Estimated position of the asteroid's center of mass relative to the recieve point
        # See equation (3) of https://doi.org/10.1086/116062.
        ρ_vec_r = r_a_t_b - r_r_t_r
        ρ_vec_dot_r = v_a_t_b - v_r_t_r
        # Magnitude of ρ_vec_r
        # See equation (4) of https://doi.org/10.1086/116062.
        ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)

        # -R_b/c (COM correction) + Δτ_D (relativistic, tropo, iono...)
        τ_D = ρ_r/clightkms # (seconds)
        # Bounce time, new estimate
        # See equation (2) of https://doi.org/10.1086/116062.
        et_b_secs = et_r_secs - τ_D
        # Asteroid barycentric position (in au) at bounce time (TDB)
        rv_a_t_b = xva(et_b_secs)
        r_a_t_b = rv_a_t_b[1:3]
        v_a_t_b = rv_a_t_b[4:6]

        # Compute down-leg Shapiro delay
        # NOTE: when using PPN, substitute 2 -> 1+γ in expressions for Shapiro delay,
        # Δτ_rel_[D|U]

        # Earth's position at t_r
        e_D_vec  = r_e_t_r - r_s_t_r
        de_D_vec = v_e_t_r - v_s_t_r
        # Heliocentric distance of Earth at t_r
        e_D = sqrt(e_D_vec[1]^2 + e_D_vec[2]^2 + e_D_vec[3]^2)
        # Barycentric position and velocity of Sun at estimated bounce time
        rv_s_t_b = xvs(et_b_secs)
        r_s_t_b = rv_s_t_b[1:3]
        v_s_t_b = rv_s_t_b[4:6]
        # Heliocentric position and velocity of asteroid at t_b
        p_D_vec  = constant_term.(r_a_t_b - r_s_t_b)
        dp_D_vec = constant_term.(v_a_t_b - v_s_t_b)
        # Heliocentric distance of asteroid at t_b
        p_D = sqrt(p_D_vec[1]^2 + p_D_vec[2]^2 + p_D_vec[3]^2)
        # Signal path distance (down-leg)
        q_D_vec  = constant_term.(r_a_t_b - r_e_t_r)
        dq_D_vec = constant_term.(v_a_t_b - v_e_t_r)
        q_D = sqrt(q_D_vec[1]^2 + q_D_vec[2]^2 + q_D_vec[3]^2)

        # Shapiro correction to time-delay
        Δτ_rel_D = shapiro_delay(e_D, p_D, q_D) # seconds
        # Troposphere correction to time-delay
        Δτ_tropo_D = tropo_delay(R_r, ρ_vec_r)  # seconds
        # Solar corona correction to time-delay
        # Δτ_corona_D = corona_delay(constant_term.(r_a_t_b), r_r_t_r, r_s_t_r, F_tx, station_code) # seconds
        # Total time-delay
        Δτ_D = Δτ_rel_D + Δτ_tropo_D #+ Δτ_corona_D # seconds

    end
    # Total time delay
    τ_D = τ_D + Δτ_D
    # Get latest estimates of ρ_vec_r and ρ_r
    # See equation (3) of https://doi.org/10.1086/116062
    ρ_vec_r = r_a_t_b - r_r_t_r
    # See equation (4) of https://doi.org/10.1086/116062
    ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)

    # Up-leg iteration
    # τ_U first estimation
    # See equation (5) of https://doi.org/10.1086/116062
    τ_U = τ_D
    # Transmit time, 1st estimate
    # See equation (6) of https://doi.org/10.1086/116062
    et_t_secs = et_r_secs - (τ_U+τ_D)
    # Geocentric position and velocity of transmitting antenna in inertial frame (au, au/day)
    R_t, V_t = obs_pv_ECI(station_code, constant_term(et_t_secs), eo=eo)
    # Barycentric position and velocity of the Earth at transmit time
    rv_e_t_t = xve(et_t_secs)
    r_e_t_t = rv_e_t_t[1:3]
    v_e_t_t = rv_e_t_t[4:6]
    # Transmitter barycentric position and velocity at transmit time
    r_t_t_t = r_e_t_t + R_t
    v_t_t_t = v_e_t_t + V_t
    # See equation (7) of https://doi.org/10.1086/116062
    ρ_vec_t = r_a_t_b - r_t_t_t
    ρ_vec_dot_t = v_a_t_b - v_t_t_t
    ρ_t = sqrt(ρ_vec_t[1]^2 + ρ_vec_t[2]^2 + ρ_vec_t[3]^2)

    # Allocate memmory for time delays
    Δτ_U = zero(τ_U)            # Total time delay
    Δτ_rel_U = zero(τ_U)        # Shapiro delay
    # Δτ_corona_U = zero(τ_U)   # Delay due to Solar corona
    Δτ_tropo_U = zero(τ_U)      # Delay due to Earth's troposphere

    for i in 1:niter
        # See equation (8) of https://doi.org/10.1086/116062
        # -R_b/c (COM correction) + Δτ_U (relativistic, tropo, iono...)
        τ_U = ρ_t/clightkms # (seconds)
        # Transmit time, new estimate
        # See equation (6) of https://doi.org/10.1086/116062
        et_t_secs = et_r_secs-(τ_U+τ_D)
        # Geocentric position and velocity of transmitting antenna in inertial frame (au, au/day)
        R_t, V_t = obs_pv_ECI(station_code, constant_term(et_t_secs), eo=eo)
        # Earth's barycentric position and velocity at the transmit time
        rv_e_t_t = xve(et_t_secs)
        r_e_t_t = rv_e_t_t[1:3]
        v_e_t_t = rv_e_t_t[4:6]
        # Barycentric position and velocity of the transmitter at the transmit time
        r_t_t_t = r_e_t_t + R_t
        v_t_t_t = v_e_t_t + V_t
        # See equation (7) of https://doi.org/10.1086/116062
        ρ_vec_t = r_a_t_b - r_t_t_t
        ρ_vec_dot_t = v_a_t_b - v_t_t_t
        ρ_t = sqrt(ρ_vec_t[1]^2 + ρ_vec_t[2]^2 + ρ_vec_t[3]^2)

        # Compute up-leg Shapiro delay

        # Sun barycentric position and velocity (in au, au/day) at transmit time (TDB)
        rv_s_t_t = xvs( et_t_secs )
        r_s_t_t = rv_s_t_t[1:3]
        v_s_t_t = rv_s_t_t[4:6]
        # Earth barycentric position and velocity at transmit time
        e_U_vec = constant_term.(r_e_t_t - r_s_t_t)
        de_U_vec = constant_term.(v_e_t_t - v_s_t_t)
        # Heliocentric distance of Earth at t_t
        e_U = sqrt(e_U_vec[1]^2 + e_U_vec[2]^2 + e_U_vec[3]^2)
        # Barycentric position/velocity of Sun at bounce time
        rv_s_t_b = xvs(et_b_secs)
        r_s_t_b = rv_s_t_b[1:3]
        v_s_t_b = rv_s_t_b[4:6]
        # Heliocentric position and velocity of the asteroid at bounce time
        p_U_vec = constant_term.(r_a_t_b - r_s_t_b)
        dp_U_vec = constant_term.(v_a_t_b - v_s_t_b)
        # Heliocentric distance of asteroid at t_b
        p_U = sqrt(p_U_vec[1]^2 + p_U_vec[2]^2 + p_U_vec[3]^2)
        # Signal path (up-leg)
        q_U_vec = constant_term.(r_a_t_b - r_e_t_t)
        dq_U_vec = constant_term.(v_a_t_b - v_e_t_t)
        q_U = sqrt(q_U_vec[1]^2 + q_U_vec[2]^2 + q_U_vec[3]^2)

        # Time-delays
        # Shapiro delay
        Δτ_rel_U = shapiro_delay(e_U, p_U, q_U) # seconds
        # Troposphere correction to time delay
        Δτ_tropo_U = tropo_delay(R_t, ρ_vec_t)  # seconds
        # Solar corona correctino to time delay
        # Δτ_corona_U = corona_delay(constant_term.(r_t_t_t), constant_term.(r_a_t_b), constant_term.(r_s_t_b), F_tx, station_code) # seconds
        # Total time delay
        Δτ_U = Δτ_rel_U + Δτ_tropo_U #+ Δτ_corona_U # seconds
    end
    # Sun barycentric position (in au) at transmit time (TDB)
    r_s_t_t = xvs(et_t_secs)[1:3]
    # Total time delay
    τ_U = τ_U + Δτ_U

    # Compute TDB-UTC at transmit time
    # Corrections to TT-TDB from Moyer (2003) / Folkner et al. (2014) due to position
    # of measurement station on Earth are of order 0.01μs
    # Δtt_tdb_station_t = - dot(v_e_t_t/daysec, r_t_t_t-r_e_t_t)/c_au_per_sec^2
    tdb_utc_t = tdb_utc(constant_term(et_t_secs)) #+ Δtt_tdb_station_t

    # Compute TDB-UTC at receive time
    # Corrections to TT-TDB from Moyer (2003) / Folkner et al. (2014) due to position
    # of measurement station on Earth  are of order 0.01μs
    # Δtt_tdb_station_r = dot(v_e_t_r/daysec, r_r_t_r-r_e_t_r)/c_au_per_sec^2
    tdb_utc_r = tdb_utc(et_r_secs) #+ Δtt_tdb_station_r

    # Compute total time delay (UTC seconds)
    # See equation (9) of https://doi.org/10.1086/116062
    τ = (τ_U + τ_D) + (tdb_utc_t - tdb_utc_r)

    # Compute Doppler shift ν
    # See equation (10) of https://doi.org/10.1086/116062
    ρ_vec_dot_t = v_a_t_b - v_t_t_t
    ρ_vec_dot_r = v_a_t_b - v_r_t_r
    # See equation (11) of https://doi.org/10.1086/116062
    ρ_dot_t = dot(ρ_vec_t, ρ_vec_dot_t)/ρ_t
    ρ_dot_r = dot(ρ_vec_r, ρ_vec_dot_r)/ρ_r

    # See equation (12) of https://doi.org/10.1086/116062
    # Order 1, only term
    doppler_c = (ρ_dot_t+ρ_dot_r)/clightkms
    # Order 1/c^2, first three terms
    p_t = dot(ρ_vec_t, v_t_t_t)/ρ_t
    p_r = dot(ρ_vec_r, v_r_t_r)/ρ_r
    doppler_c2_t1 = ρ_dot_t*p_t - ρ_dot_r*p_r - ρ_dot_t*ρ_dot_r
    # Order 1/c^2, GM_S terms
    r_ts_vec = r_t_t_t - r_s_t_t
    r_ts = sqrt(r_ts_vec[1]^2+r_ts_vec[2]^2+r_ts_vec[3]^2)
    r_rs_vec = r_r_t_r - r_s_t_r
    r_rs = sqrt(r_rs_vec[1]^2+r_rs_vec[2]^2+r_rs_vec[3]^2)
    doppler_c2_t2 = (PlanetaryEphemeris.μ[su]*((au^3)/(daysec^2)))*( (1/r_ts) - (1/r_rs) )
    # Order 1/c^2, last term
    doppler_c2_t3 = (  dot(v_t_t_t, v_t_t_t) - dot(v_r_t_r, v_r_t_r)  )/2

    # Total Doppler shift
    doppler_c2 = (doppler_c2_t1 + doppler_c2_t2 + doppler_c2_t3)/(clightkms^2)
    ν = -F_tx*(doppler_c) # + doppler_c2)

    # Total signal delay (μs) and Doppler shift (Hz)
    return 1e6τ, 1e6ν
end
