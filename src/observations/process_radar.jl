include("catalogue_mpc.jl")
include("observatory_mpc.jl")
include("radec_mpc.jl")
include("radar_jpl.jl")
include("units.jl")
include("jpl_eph.jl")
include("topocentric.jl")
include("process_radec.jl")

@doc raw"""
    shapiro_delay(e, p, q)

Return the relativistic (Shapiro) time-delay in seconds
```math
\Delta\tau[\text{rel}] = \frac{2\mu_\odot}{c^3}\log\left|\frac{d_{E,S} + d_{A,S} + d_{A,E}}{d_{E,S}+d_{A,S}-d_{A, E}}\right|,
```
where ``\mu_\odot = GM_\odot`` is the gravitational parameter of the sun, and ``d_{E,S}``, ``d_{A,S}``
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
    shap_del_days = (2μ_DE430[su]/(c_au_per_day^3))*log( (e+p+q+shap)/(e+p-q+shap) ) # days
    return shap_del_days*daysec # seconds
end

@doc raw"""
    shapiro_doppler(e, de, p, dp, q, dq, F_tx)

Return the Doppler shift (in units of `F_tx`)
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
    shap_del_diff = (4μ_DE430[su]/(c_au_per_day^3))*(  ( dq*(e+p) - q*(de+dp) )/( (e+p)^2 - q^2 )  ) # differential of Shapiro delay (adim.)
    # ν = -F_tx*dτ/dt (units of F_tx)
    # See footnote 10 of https://doi.org/10.1103/PhysRevLett.17.933
    shap_dop = -F_tx*shap_del_diff
    return shap_dop # (units of F_tx)
end

@doc raw"""
    Ne(p1::Vector{S}, p2::Vector{S}, r_s_t0::Vector{S}, ds::U, ΔS::Real) where {S<:Number, U<:Number}

Return the density of ionized electrons (in electrons/cm^3) in interplanetary medium
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
- `p1::Vector{S}`: signal departure point (transmitter/bounce) (au).
- `p2::Vector{S}`: signal arrival point (bounce/receiver) (au).
- `r_s_t0::Vector{S}`: Barycentric position (au) of Sun at initial time of propagation of signal path (bounce time for down-leg; transmit time for up-leg).
- `ds::U`: current distance travelled by ray from emission point (au).
- `ΔS::Real`: total distance between p1 and p2 (au).
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

Return the path integral of the density of ionized electrons in interplanetary medium ``N_e``
in units of electrons/cm^2, evaluated with `TaylorIntegration`.

# Arguments

- `p1::Vector{S}`: signal departure point (transmitter/bounce) (au).
- `p2::Vector{S}`: signal arrival point (bounce/receiver) (au).
- `r_s_t0::Vector{S}`: Barycentric position (au) of Sun at initial time of propagation of signal path (bounce time for down-leg; transmit time for up-leg).
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
    corona_delay(p1::Vector{S}, p2::Vector{S}, r_s_t0::Vector{S}, F_tx::U) where {S<:Number, U<:Real}

Return the time-delay (in sec) due to thin plasma of solar corona
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

- `p1::Vector{S}`: signal departure point (transmitter/bounce for up/down-link, resp.) (au).
- `p2::Vector{S}`: signal arrival point (bounce/receiver for up/down-link, resp.) (au).
- `r_s_t0::Vector{S}`: Barycentric position (au) of Sun at initial time of propagation of signal path (bounce time for down-leg; transmit time for up-leg).
- `F_tx::U`: transmitter frequency (MHz).
"""
function corona_delay(p1::Vector{S}, p2::Vector{S}, r_s_t0::Vector{S}, F_tx::U) where {S<:Number, U<:Real}
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

- `r_antenna::Vector{T}`: position of antenna at receive/transmit time in celestial frame wrt geocenter.
- `ρ_vec_ae::Vector{S}`: slant-range vector from antenna to asteroid.
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

Return the time delay (in sec) due to Earth's troposphere for radio frequencies
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

Return the time delay (in sec) due to Earth's troposphere for radio frequencies. The function
first computes the zenith distance ``z`` via [`zenith_distance`](@ref) and then substitutes
into the first method of [`tropo_delay`](@ref).

# Arguments

- `r_antenna::Vector{T}`: position of antenna at receive/transmit time in celestial frame wrt geocenter.
- `ρ_vec_ae::Vector{S}`: slant-range vector from antenna to asteroid.
"""
function tropo_delay(r_antenna::Vector{T}, ρ_vec_ae::Vector{S}) where {T<:Number, S<:Number}
    # zenith distance
    zd = zenith_distance(r_antenna, ρ_vec_ae) # rad
    # Time delay due to Earth's troposphere
    return tropo_delay(zd) # seconds
end

@doc raw"""
    compute_delay(observatory::ObservatoryMPC{T}, t_r_utc::DateTime, t_offset::Real, niter::Int = 10; eo::Bool = true,
          xve = earthposvel, xvs = sunposvel, xva = apophisposvel197) where {T <: AbstractFloat}
    compute_delay(radar::RadarJPL{T}, t_offset::Real, niter::Int = 10; eo::Bool = true, xve = earthposvel, xvs = sunposvel,
          xva = apophisposvel197) where {T <: AbstractFloat}

Compute radar-astrometric round-trip time (``\mu``s).

See https://doi.org/10.1086/116062.

# Arguments

- `observatory::ObservatoryMPC{T}`: observing station.
- `t_r_utc::DateTime`: UTC time of echo reception.
- `radar::RadarJPL{T}`: radar observation.
- `t_offset::Real`: time offset wrt echo reception time, to compute Doppler shifts by range differences (seconds).
- `niter::Int`: number of light-time solution iterations.
- `eo::Bool`: compute corrections due to Earth orientation, LOD, polar motion.
- `xve`: Earth ephemeris [et seconds since J2000] -> [barycentric position in km and velocity in km/sec].
- `xvs`: Sun ephemeris [et seconds since J2000] -> [barycentric position in km and velocity in km/sec].
- `xva`: asteroid ephemeris [et seconds since J2000] -> [barycentric position in km and velocity in km/sec].
"""
function compute_delay(observatory::ObservatoryMPC{T}, t_r_utc::DateTime, t_offset::Real, niter::Int = 10; eo::Bool = true,
               xve = earthposvel, xvs = sunposvel, xva = apophisposvel197) where {T <: AbstractFloat}
    # Transform receiving time from UTC to TDB seconds since j2000
    et_r_secs = datetime2et(t_r_utc) + t_offset
    # Compute geocentric position/velocity of receiving antenna in inertial frame (au, au/day)
    RV_r = obsposvelECI(observatory, et_r_secs, eo = eo)
    R_r = RV_r[1:3]
    V_r = RV_r[4:6]
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
    RV_t = obsposvelECI(observatory, et_t_secs, eo = eo)
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

    # Allocate memmory for time delays
    Δτ_U = zero(τ_U)            # Total time delay
    Δτ_rel_U = zero(τ_U)        # Shapiro delay
    # Δτ_corona_U = zero(τ_U)   # Delay due to Solar corona
    Δτ_tropo_U = zero(τ_U)      # Delay due to Earth's troposphere

    for i in 1:niter
        # Geocentric position and velocity of transmitting antenna in inertial frame (au, au/day)
        # TODO: remove `constant_term` to take into account dependency of R_t, V_t wrt initial
        # conditions variations via et_t_secs
        RV_t = obsposvelECI(observatory, et_t_secs, eo = eo)
        R_t = RV_t[1:3]
        V_t = RV_t[4:6]
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
function compute_delay(radar::RadarJPL{T}, t_offset::Real, niter::Int = 10; eo::Bool = true, xve = earthposvel, xvs = sunposvel,
               xva = apophisposvel197) where {T <: AbstractFloat}
    return compute_delay(radar.rcvr, radar.date, t_offset, niter; eo = eo, xve = xve, xvs = xvs, xva = xva)
end

@doc raw"""
    compute_delay(observatory::ObservatoryMPC{T}, t_r_utc::DateTime, niter::Int = 10; eo::Bool = true, xve::TaylorInterpolant = earthposvel,
          xvs::TaylorInterpolant = sunposvel, xva::TaylorInterpolant = apophisposvel197, tord::Int = xva.x[1].order) where {T <: AbstractFloat}
    compute_delay(radar::RadarJPL{T}, niter::Int = 10; eo::Bool = true, xve::TaylorInterpolant = earthposvel,
          xvs::TaylorInterpolant = sunposvel, xva::TaylorInterpolant = apophisposvel197, tord::Int = xva.x[1].order) where {T <: AbstractFloat}

Compute Taylor series expansion of time-delay observable around echo reception time. This
allows to compute dopplers via automatic differentiation using
```math
\nu = -f\frac{d\tau}{dt},
```
where ``f`` is the transmitter frequency (MHz) and ``\tau`` is the time-delay at reception
time ``t``.

See https://doi.org/10.1086/116062.

**Note:** Works only with `TaylorInterpolant` ephemeris. See [`PlanetaryEphemeris.TaylorInterpolant`](@ref).

# Arguments

- `observatory::ObservatoryMPC{T}`: observing station.
- `t_r_utc::DateTime`: UTC time of echo reception.
- `radar::RadarJPL{T}`: radar observation.
- `niter::Int`: number of light-time solution iterations.
- `eo::Bool`: compute corrections due to Earth orientation, LOD, polar motion.
- `xve::TaylorInterpolant`: Earth ephemeris [et seconds since J2000] -> [barycentric position in km and velocity in km/sec].
- `xvs::TaylorInterpolant`: Sun ephemeris [et seconds since J2000] -> [barycentric position in km and velocity in km/sec].
- `xva::TaylorInterpolant`: asteroid ephemeris [et seconds since J2000] -> [barycentric position in km and velocity in km/sec].
- `tord::Int`: order of Taylor expansions.
"""
function compute_delay(observatory::ObservatoryMPC{T}, t_r_utc::DateTime, niter::Int = 10; eo::Bool = true, xve::Function = earthposvel,
               xvs::Function = sunposvel, xva::Function = apophisposvel197, tord::Int = xva.x[1].order) where {T <: AbstractFloat}

    # Transform receiving time from UTC to TDB seconds since j2000
    et_r_secs_0 = datetime2et(t_r_utc)
    # Auxiliary to evaluate JT ephemeris
    xva1et0 = xva(et_r_secs_0)[1]
    # et_r_secs_0 as a Taylor polynomial
    et_r_secs = Taylor1([et_r_secs_0,1.0].*one(xva1et0), tord)
    # Compute geocentric position/velocity of receiving antenna in inertial frame (au, au/day)
    RV_r = obsposvelECI(observatory, et_r_secs, eo = eo)
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
    RV_t = obsposvelECI(observatory, et_t_secs, eo = eo)
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

    # Allocate memmory for time delays
    Δτ_U = zero(τ_U)            # Total time delay
    Δτ_rel_U = zero(τ_U)        # Shapiro delay
    # Δτ_corona_U = zero(τ_U)   # Delay due to Solar corona
    Δτ_tropo_U = zero(τ_U)      # Delay due to Earth's troposphere

    for i in 1:niter
        # Geocentric position and velocity of transmitting antenna in inertial frame (au, au/day)
        # TODO: remove `constant_term` to take into account dependency of R_t, V_t wrt initial
        # conditions variations via et_t_secs
        RV_t = obsposvelECI(observatory, et_t_secs, eo = eo)
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

function compute_delay(radar::RadarJPL{T}, niter::Int = 10; eo::Bool = true, xve::TaylorInterpolant = earthposvel,
               xvs::TaylorInterpolant = sunposvel, xva::TaylorInterpolant = apophisposvel197, tord::Int = xva.x[1].order) where {T <: AbstractFloat}
    return compute_delay(radar.rcvr, radar.date, niter; eo = eo, xve = xve, xvs = xvs, xva = xva, tord = tord)
end

@doc raw"""
    radar_astrometry(observatory::ObservatoryMPC{T}, t_r_utc::DateTime, F_tx::Real, niter::Int = 10; eo::Bool = true, tc::Real = 1.0,
                  xve = earthposvel, xvs = sunposvel, xva = apophisposvel197, autodiff::Bool = true, tord::Int = 10) where {T <: AbstractFloat}
    radar_astrometry(radar::RadarJPL{T}, niter::Int = 10; eo::Bool = true, tc::Real = 1.0, xve = earthposvel,
                  xvs = sunposvel, xva = apophisposvel197, autodiff::Bool = true, tord::Int = 10) where {T <: AbstractFloat}
    radar_astrometry(astradarfile::String, niter::Int = 10; eo::Bool = true, tc::Real = 1.0, xve = earthposvel, xvs = sunposvel,
                  xva = apophisposvel197, autodiff::Bool = true, tord::Int=10)
    radar_astrometry(outfilename::String, radarobsfile::String, asteph::TaylorInterpolant, ss16asteph::TaylorInterpolant; tc::Real=1.0, autodiff::Bool=true,
                  tord::Int=10, niter::Int=5)

Return time-delay and Doppler shift.

# Arguments

- `observatory::ObservatoryMPC{T}`: observing station.
- `t_r_utc::DateTime`: UTC time of echo reception.
- `radar::RadarJPL{T}`: radar observation.
- `astradarfile/radarobsfile::String`: file where to retrieve radar observations.
- `outfilename::String`: file where to save radar observations.
- `F_tx::Real`: transmitter frequency (MHz).
- `niter::Int`: number of light-time solution iterations.
- `eo::Bool`: wheter to use `EarthOrientation` or not.
- `tc::Real`: time offset wrt echo reception time, to compute Doppler shifts by range differences (seconds).
- `xve`: Earth ephemeris [et seconds since J2000] -> [barycentric position in km and velocity in km/sec].
- `xvs`: Sun ephemeris [et seconds since J2000] -> [barycentric position in km and velocity in km/sec].
- `xva`: asteroid ephemeris [et seconds since J2000] -> [barycentric position in km and velocity in km/sec].
- `asteph::TaylorInterpolant`: asteroid's ephemeris.
- `ss16asteph::TaylorInterpolant`: solar system ephemeris.
- `autodiff::Bool`: whether to use the automatic differentiation method of [`compute_delay`](@ref) or not.
- `tord::Int`: order of Taylor expansions.
"""
function radar_astrometry(observatory::ObservatoryMPC{T}, t_r_utc::DateTime, F_tx::Real, niter::Int = 10; eo::Bool = true, tc::Real = 1.0,
                       xve = earthposvel, xvs = sunposvel, xva = apophisposvel197, autodiff::Bool = true, tord::Int = 10) where {T <: AbstractFloat}

    # Automatic differentiation method of delay
    if autodiff
        # Time delay
        τ = compute_delay(observatory, t_r_utc, niter, eo=eo, xve=xve, xvs=xvs, xva=xva, tord=tord)
        # Time delay, Doppler shift
        return τ[0], -F_tx*τ[1]
    # No automatic differentiation method of delay
    else
        τe = compute_delay(observatory, t_r_utc,  tc/2, niter, eo=eo, xve=xve, xvs=xvs, xva=xva)
        τn = compute_delay(observatory, t_r_utc,   0.0, niter, eo=eo, xve=xve, xvs=xvs, xva=xva)
        τs = compute_delay(observatory, t_r_utc, -tc/2, niter, eo=eo, xve=xve, xvs=xvs, xva=xva)
        # Time delay, Doppler shift
        return τn, -F_tx*((τe-τs)/tc)
    end

end

function radar_astrometry(radar::RadarJPL{T}, niter::Int = 10; eo::Bool = true, tc::Real = 1.0, xve = earthposvel,
                       xvs = sunposvel, xva = apophisposvel197, autodiff::Bool = true, tord::Int = 10) where {T <: AbstractFloat}
    return radar_astrometry(radar.rcvr, radar.date, radar.freq, niter; eo = eo, tc = tc, xve = xve, xvs = xvs, xva = xva, autodiff = autodiff, tord = tord)
end

function radar_astrometry(astradarfile::String, niter::Int=10; kwargs...)
    # Check that astradarfile is a file
    @assert isfile(astradarfile) "Cannot open file: $astradarfile"
    # Read radar measurements
    astradardata = read_radar_jpl(astradarfile)
    return radar_astrometry(astradardata, niter; kwargs...)
end

function radar_astrometry(astradardata::Vector{RadarJPL{T}}, niter::Int=10; eo::Bool=true, tc::Real=1.0, xve=earthposvel, xvs=sunposvel,
                       xva=apophisposvel197, autodiff::Bool=true, tord::Int=10) where {T <: AbstractFloat}

    # UTC time of first radar observation
    utc1 = astradardata[1].date
    # TDB seconds since J2000.0 for first astrometric observation
    et1 = datetime2et(utc1)
    # Asteroid ephemeris at et1
    a1_et1 = xva(et1)[1]
    # Type of asteroid ephemeris
    S = typeof(a1_et1)

    # Time delays
    vdelay = Array{S}(undef, length(astradardata))
    # Doppler shifts
    vdoppler = Array{S}(undef, length(astradardata))

    # Iterate over the measurements
    for i in eachindex(astradardata)
        # Compute time delay and doppler shift
        vdelay[i], vdoppler[i] = radar_astrometry(astradardata[i], niter; eo = eo, tc = tc, xve = xve, xvs = xvs, xva = xva,
                                               autodiff = autodiff, tord = tord)
    end
    # Rows with time delays
    delay_index = map(x-> x.Δτ_units == "us", astradardata)
    # Rows with Doppler shifts
    doppler_index = map(x -> x.Δν_units == "Hz", astradardata)

    dt_utc_obs = date.(astradardata)           # UTC time
    Δτ_obs = delay.(astradardata)              # Observed time delay
    Δν_obs = doppler.(astradardata)            # Observed Doppler shift
    Δτ_σ = delay_sigma.(astradardata)          # Observed time delay uncertainty
    Δν_σ = doppler_sigma.(astradardata)        # Observed Doppler shift uncertainty
    τ_units = delay_units.(astradardata)       # Time delay units
    ν_units = doppler_units.(astradardata)     # Doppler shift units
    freq_ = freq.(astradardata)                # Frequency
    rcvr_ = rcvr.(astradardata)                # Receiver antenna
    xmit_ = xmit.(astradardata)                # Emission antenna
    bouncepoint_ = bouncepoint.(astradardata)   # Bounce point

    # Return time delays / Doppler shifts table
    return dt_utc_obs, Δτ_obs, Δν_obs, vdelay, vdoppler, Δτ_σ, Δν_σ, τ_units, ν_units, freq_, rcvr_, xmit_, bouncepoint_, delay_index,
           doppler_index
end

@doc raw"""
    radar_astrometry_yeomansetal92(observatory::ObservatoryMPC{T}, t_r_utc::DateTime, F_tx::Real, niter::Int = 10; eo::Bool = true,
                                xve = earthposvel, xvs = sunposvel, xva = apophisposvel197) where {T <: AbstractFloat}
    radar_astrometry_yeomansetal92(radar::RadarJPL{T}, niter::Int = 10; eo::Bool = true, xve = earthposvel, xvs = sunposvel,
                                xva = apophisposvel197) where {T <: AbstractFloat}

Compute round-trip time (in ``\mu``s) and Doppler shift (Hz) radar astrometry. Dopplers are computed following https://doi.org/10.1086/116062.

# Arguments

- `observatory::ObservatoryMPC{T}`: observing station.
- `t_r_utc::DateTime`: UTC time of echo reception.
- `radar::RadarJPL{T}`: radar observation.
- `F_tx::Real`: transmitter frequency (MHz).
- `niter::Int`: number of light-time solution iterations.
- `eo::Bool`: wheter to use `EarthOrientation` or not.
- `xve`: Earth ephemeris [et seconds since J2000] -> [barycentric position in km and velocity in km/sec].
- `xvs`: Sun ephemeris [et seconds since J2000] -> [barycentric position in km and velocity in km/sec].
- `xva`: asteroid ephemeris [et seconds since J2000] -> [barycentric position in km and velocity in km/sec].
"""
function radar_astrometry_yeomansetal92(observatory::ObservatoryMPC{T}, t_r_utc::DateTime, F_tx::Real, niter::Int = 10; eo::Bool = true,
                                     xve = earthposvel, xvs = sunposvel, xva = apophisposvel197) where {T <: AbstractFloat}
    # Transform receiving time from UTC to TDB seconds since j2000
    et_r_secs = datetime2et(t_r_utc)
    # Compute geocentric position/velocity of receiving antenna in inertial frame (au, au/day)
    RV_r = obsposvelECI(observatory, et_r_secs, eo = eo)
    R_r = RV_r[1:3]
    V_r = RV_r[4:6]
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
    RV_t = obsposvelECI(observatory, constant_term(et_t_secs), eo = eo)
    R_t = RV_t[1:3]
    V_t = RV_t[4:6]
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
        RV_t = obsposvelECI(observatory, constant_term(et_t_secs), eo = eo)
        R_t = RV_t[1:3]
        V_t = RV_t[4:6]
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
    doppler_c2_t2 = (μ_DE430[su]*((au^3)/(daysec^2)))*( (1/r_ts) - (1/r_rs) )
    # Order 1/c^2, last term
    doppler_c2_t3 = (  dot(v_t_t_t, v_t_t_t) - dot(v_r_t_r, v_r_t_r)  )/2

    # Total Doppler shift
    doppler_c2 = (doppler_c2_t1 + doppler_c2_t2 + doppler_c2_t3)/(clightkms^2)
    ν = -F_tx*(doppler_c) # + doppler_c2)

    # Total signal delay (μs) and Doppler shift (Hz)
    return 1e6τ, 1e6ν
end

function radar_astrometry_yeomansetal92(radar::RadarJPL{T}, niter::Int = 10; eo::Bool = true, xve = earthposvel, xvs = sunposvel,
                                     xva = apophisposvel197) where {T <: AbstractFloat}
    return radar_astrometry_yeomansetal92(radar.observatory, radar.date, radar.freq, niter, eo = eo, xve = xve, xvs = xvs, xva = xva)
end

function radar_astrometry(outfilename::String, radarobsfile::String, asteph::TaylorInterpolant, ss16asteph::TaylorInterpolant;
                           tc::Real = 1.0, autodiff::Bool = true, tord::Int = 10, niter::Int = 5)
    # Check that radarobsfile is a file
    @assert isfile(radarobsfile) "Cannot open file: $radarobsfile"
    # Read optical observations
    radar = read_radar_jpl(radarobsfile)

    # Check that first and last observation times are within interpolation interval
    Δt_0 = datetime2julian(radar[1].date) - JD_J2000 - asteph.t0
    Δt_f = datetime2julian(radar[end].date) - JD_J2000 - asteph.t0
    t_min, t_max = minmax(asteph.t[1], asteph.t[end])
    @assert t_min ≤ Δt_0 ≤ Δt_f ≤ t_max "First and/or last observation times are outside interpolation interval"

    # Number of massive bodies
    Nm1 = (size(ss16asteph.x)[2]-13) ÷ 6
    # Number of bodies, including NEA
    N = Nm1 + 1

    # NEO
    # Change t, x, v units, resp., from days, au, au/day to sec, km, km/sec
    asteph_et(et) = auday2kmsec(asteph(et/daysec)[1:6])
    # Earth
    # Change x, v units, resp., from au, au/day to km, km/sec
    eph_ea = selecteph(ss16asteph, ea)
    xve(et) = auday2kmsec(eph_ea(et))
    # Sun
    # Change x, v units, resp., from au, au/day to km, km/sec
    eph_su = selecteph(ss16asteph, su)
    xvs(et) = auday2kmsec(eph_su(et))

    # Compute radar astrometry
    dt_utc_obs, Δτ_obs, Δν_obs, vdelay, vdoppler, Δτ_σ, Δν_σ, τ_units, ν_units, freq, rcvr, xmit, bouncepoint, delay_index,
           doppler_index = radar_astrometry(radarobsfile, niter, xve = earth_et, xvs = sun_et, xva = asteph_et, tc = tc, tord = tord, autodiff = autodiff)
    # Save data to file
    println("Saving data to file: $outfilename")
    JLD2.jldopen(outfilename, "w") do file
        # Write variables to jld2 file
        JLD2.write(file,
        "dt_utc_obs", dt_utc_obs,
        "Δτ_obs", Δτ_obs,
        "Δν_obs", Δν_obs,
        "vdelay", vdelay,
        "vdoppler", vdoppler,
        "Δτ_σ", Δτ_σ,
        "Δν_σ", Δν_σ,
        "τ_units", τ_units,
        "ν_units", ν_units,
        "freq", freq,
        "rcvr", rcvr,
        "xmit", xmit,
        "bouncepoint", bouncepoint,
        "delay_index", delay_index,
        "doppler_index", doppler_index,
        )
    end

    return nothing
end

@doc raw"""
    residuals(obs::Vector{RadarJPL{T}}, niter::Int = 10; eo::Bool = true, debias_table::String = "2018",
              xvs::Function = sunposvel, xve::Function = earthposvel, xva::Function = apophisposvel197) where {T <: AbstractFloat}

Compute O-C residuals for radar astrometry.

See also [`compute_delay`](@ref) and [`radar_astrometry`](@ref).

# Arguments

- `obs::Vector{RadarJPL{T}}`: vector of observations.
- `niter::Int`: number of light-time solution iterations.
- `eo::Bool`: compute corrections due to Earth orientation, LOD, polar motion.
- `tc::Real`: time offset wrt echo reception time, to compute Doppler shifts by range differences (seconds).
- `autodiff::Bool`: whether to use the automatic differentiation method of [`compute_delay`](@ref) or not.
- `tord::Int`: order of Taylor expansions.
- `xvs::Function`: Sun ephemeris [et seconds since J2000] -> [barycentric position in km and velocity in km/sec].
- `xve::Function`: Earth ephemeris [et seconds since J2000] -> [barycentric position in km and velocity in km/sec].
- `xva::Function`: asteroid ephemeris [et seconds since J2000] -> [barycentric position in km and velocity in km/sec].
"""
function residuals(obs::Vector{RadarJPL{T}}, niter::Int = 10; eo::Bool = true, tc::Real = 1.0, autodiff::Bool = true, tord::Int = 10,
                   xvs::Function = sunposvel, xve::Function = earthposvel, xva::Function = apophisposvel197) where {T <: AbstractFloat}
    # Radar delay/Doppler astrometric data
    x_jt = radar_astrometry(obs, niter; eo = eo, tc=tc, autodiff = autodiff, tord = tord, xvs = xvs, xve = xve, xva = xva)
    # Time-delay residuals
    res_τ = x_jt[2][x_jt[14]] .- x_jt[4][x_jt[14]]
    # Doppler-shift residuals
    res_ν = x_jt[3][x_jt[15]] .- x_jt[5][x_jt[15]]
    # Total residuals
    res = vcat(res_τ, res_ν)
    # Weights
    w_τ = repeat(1 ./ x_jt[6][x_jt[14]].^2, 2)
    w_ν = repeat(1 ./ x_jt[7][x_jt[15]].^2, 2)
    w = vcat(w_τ, w_ν)

    return res, w
end
