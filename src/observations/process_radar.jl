include("catalogue_mpc.jl")
include("observatory_mpc.jl")
include("radec_mpc.jl")
include("radar_jpl.jl")
include("units.jl")
include("jpl_eph.jl")
include("topocentric.jl")
include("observation_night.jl")
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

# Arguments

- `e`: heliocentric distance of the Earth.
- `p`: asteroid's heliocentric distance.
- `q`: asteroid's geocentric distance.

!!! reference
    See https://doi.org/10.1103/PhysRevLett.13.789.
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

# Arguments

- `e`: heliocentric distance of the Earth.
- `de`: differential of `e`.
- `p`: asteroid's heliocentric distance.
- `dp`: differential of `p`.
- `q`: asteroid's geocentric distance.
- `dq`: differential of `q`.
- `F_tx`: transmitter frequency (MHz).

!!! reference
    See https://doi.org/10.1103/PhysRevLett.17.933.
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

# Arguments
- `p1::Vector{S}`: signal departure point (transmitter/bounce) (au).
- `p2::Vector{S}`: signal arrival point (bounce/receiver) (au).
- `r_s_t0::Vector{S}`: Barycentric position (au) of Sun at initial time of propagation of signal path (bounce time for down-leg; transmit time for up-leg).
- `ds::U`: current distance travelled by ray from emission point (au).
- `ΔS::Real`: total distance between p1 and p2 (au).

!!! reference
    See (Explanatory Supplement to the Astronomical Almanac 2014, p. 323, Sec. 8.7.5, Eq. 8.22).
    ESAA 2014 in turn refers to Muhleman and Anderson (1981). Ostro (1993) gives a reference to
    Anderson (1978), where this model is fitted to Mariner 9 ranging data. Reading 
    https://gssc.esa.int/navipedia/index.php/Ionospheric_Delay helped a lot to clarify things,
    especially the 40.3, although they talk about Earth's ionosphere. Another valuable source is
    Standish, E.M., Astron. Astrophys. 233, 252-271 (1990).
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

# Arguments

- `p1::Vector{S}`: signal departure point (transmitter/bounce for up/down-link, resp.) (au).
- `p2::Vector{S}`: signal arrival point (bounce/receiver for up/down-link, resp.) (au).
- `r_s_t0::Vector{S}`: Barycentric position (au) of Sun at initial time of propagation of signal path (bounce time for down-leg; transmit time for up-leg).
- `F_tx::U`: transmitter frequency (MHz).

!!! reference
    From [Ionospheric Delay](https://gssc.esa.int/navipedia/index.php/Ionospheric_Delay) it seems
    that ESAA 2014 text probably should say that in the formula for ``\Delta\tau_\text{cor}``,
    the expression ``40.3 N_e/f^2`` is adimensional, where ``Ne`` is in electrons/cm^3 and ``f`` is
    in Hz therefore, the integral ``(40.3/f^2)\int Ne \ ds`` is in centimeters, where ``ds`` is in
    cm and the expression ``(40.3/(cf^2))\int Ne \ ds``, with ``c`` in cm/sec, is in seconds.
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
    compute_delay(observatory::ObservatoryMPC{T}, t_r_utc::DateTime; kwargs...) where {T <: AbstractFloat}

Compute Taylor series expansion of time-delay observable around echo reception time. This
allows to compute dopplers via automatic differentiation using
```math
\nu = -f\frac{d\tau}{dt},
```
where ``f`` is the transmitter frequency (MHz) and ``\tau`` is the time-delay at reception
time ``t``. Computed values include corrections due to Earth orientation, LOD and polar
motion.

**Note:** Works only with `TaylorInterpolant` ephemeris. See [`PlanetaryEphemeris.TaylorInterpolant`](@ref).

# Arguments

- `observatory::ObservatoryMPC{T}`: observing station.
- `t_r_utc::DateTime`: UTC time of echo reception.

# Keyword arguments 

- `tord::Int = 5`: order of Taylor expansions.
- `niter::Int = 10`: number of light-time solution iterations.
- `xve::EarthEph = earthposvel`: Earth ephemeris.
- `xvs::SunEph = sunposvel`: Sun ephemeris.
- `xva::AstEph`: asteroid ephemeris.

All ephemeris must take  [et seconds since J2000] and return [barycentric position in km
and velocity in km/sec].

!!! reference
    See https://doi.org/10.1086/116062.
"""
function compute_delay(observatory::ObservatoryMPC{T}, t_r_utc::DateTime; tord::Int = 5,
                       niter::Int = 10, xve::EarthEph = earthposvel, xvs::SunEph = sunposvel,
                       xva::AstEph) where {T <: AbstractFloat, EarthEph, SunEph, AstEph}

    # Transform receiving time from UTC to TDB seconds since j2000
    et_r_secs_0 = datetime2et(t_r_utc)
    # Auxiliary to evaluate JT ephemeris
    xva1et0 = xva(et_r_secs_0)[1]
    # et_r_secs_0 as a Taylor polynomial
    et_r_secs = Taylor1([et_r_secs_0,one(et_r_secs_0)].*one(xva1et0), tord)
    # Compute geocentric position/velocity of receiving antenna in inertial frame (km, km/sec)
    RV_r = obsposvelECI(observatory, et_r_secs)
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
    # Geocentric position and velocity of transmitting antenna in inertial frame (km, km/sec)
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

    # Allocate memmory for time delays
    Δτ_U = zero(τ_U)            # Total time delay
    Δτ_rel_U = zero(τ_U)        # Shapiro delay
    # Δτ_corona_U = zero(τ_U)   # Delay due to Solar corona
    Δτ_tropo_U = zero(τ_U)      # Delay due to Earth's troposphere

    for i in 1:niter
        # Geocentric position and velocity of transmitting antenna in inertial frame (km, km/sec)
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

        # Sun barycentric position and velocity (in km, km/sec) at transmit time (TDB)
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
    # of measurement station on Earth are of order 0.01μs
    # Δtt_tdb_station_r = - dot(v_e_t_r, r_r_t_r-r_e_t_r)/clightkms^2
    tdb_utc_r = tdb_utc(et_r_secs) # + Δtt_tdb_station_r

    # Compute total time delay (UTC seconds); relativistic delay is already included in τ_D, τ_U
    # See equation (9) of https://doi.org/10.1086/116062
    τ = (τ_D + τ_U) + (Δτ_tropo_D + Δτ_tropo_U) + (tdb_utc_t - tdb_utc_r) # seconds

    # Total signal delay (μs)
    return 1e6τ
end

@doc raw"""
    radar_astrometry(observatory::ObservatoryMPC{T}, t_r_utc::DateTime, F_tx::Real; kwargs...) where {T <: AbstractFloat}
    radar_astrometry(radar::RadarJPL{T}; kwargs...) where {T <: AbstractFloat}
    radar_astrometry(astradarfile::String; kwargs...)
    radar_astrometry(astradardata::Vector{RadarJPL}; kwargs...) where {T <: AbstractFloat} 

Return time-delay and Doppler shift.

# Arguments

- `observatory::ObservatoryMPC{T}`: observing station.
- `t_r_utc::DateTime`: UTC time of echo reception.
- `F_tx::Real`: transmitter frequency (MHz).
- `radar::RadarJPL{T}` / `astradardata::Vector{RadarJPL{T}}`: radar observation(s).
- `astradarfile::String`: file where to retrieve radar observations.

# Keyword arguments

- `tc::Real = 1.0`: time offset wrt echo reception time, to compute Doppler shifts by range differences (seconds). Offsets will be rounded to the nearest millisecond.
- `autodiff::Bool = true`: whether to compute Doppler shift via automatic diff of [`compute_delay`](@ref) or not.
- `tord::Int = 5`: order of Taylor expansions.
- `niter::Int = 10`: number of light-time solution iterations.
- `xve::EarthEph = earthposvel`: Earth ephemeris.
- `xvs::SunEph = sunposvel`: Sun ephemeris.
- `xva::AstEph`: asteroid ephemeris.

All ephemeris must take  [et seconds since J2000] and return [barycentric position in km
and velocity in km/sec].
"""
function radar_astrometry(observatory::ObservatoryMPC{T}, t_r_utc::DateTime, F_tx::Real;
                          tc::Real = 1.0, autodiff::Bool = true, kwargs...) where {T <: AbstractFloat}
    # Compute Doppler shift via automatic differentiation of time-delay
    if autodiff
        # Time delay
        τ = compute_delay(observatory, t_r_utc; kwargs...)
        # Time delay, Doppler shift
        return τ[0], -F_tx*τ[1]
    # Compute Doppler shift via numerical differentiation of time-delay
    else
        offset = Dates.Millisecond(1000round(tc/2, digits=3))
        τe = compute_delay(observatory, t_r_utc+offset; kwargs...)
        τn = compute_delay(observatory, t_r_utc       ; kwargs...)
        τs = compute_delay(observatory, t_r_utc-offset; kwargs...)
        # Time delay, Doppler shift
        return τn[0], -F_tx*((τe[0]-τs[0])/tc)
    end

end

function radar_astrometry(radar::RadarJPL{T}; kwargs...) where {T <: AbstractFloat}
    return radar_astrometry(radar.rcvr, radar.date, radar.freq; kwargs...)
end

function radar_astrometry(astradarfile::String; kwargs...)
    # Check that astradarfile is a file
    @assert isfile(astradarfile) "Cannot open file: $astradarfile"
    # Read radar measurements
    astradardata = read_radar_jpl(astradarfile)
    return radar_astrometry(astradardata; kwargs...)
end

function radar_astrometry(astradardata::Vector{RadarJPL{T}}; xva::AstEph, kwargs...) where {T <: AbstractFloat, AstEph}

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
    Threads.@threads for i in eachindex(astradardata)
        # Compute time delay and doppler shift
        vdelay[i], vdoppler[i] = radar_astrometry(astradardata[i]; xva, kwargs...)
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
    residuals(obs::Vector{RadarJPL{T}}; kwargs...)  where {T <: AbstractFloat}

Compute O-C residuals for radar astrometry, including corrections due to Earth orientation,
LOD and polar motion.

See also [`compute_delay`](@ref) and [`radar_astrometry`](@ref).

# Arguments

- `obs::Vector{RadarJPL{T}}`: vector of observations.

# Keyword arguments

- `tc::Real = 1.0`: time offset wrt echo reception time, to compute Doppler shifts by range differences (seconds). Offsets will be rounded to the nearest millisecond.
- `autodiff::Bool = true`: whether to compute Doppler shift via automatic diff of [`compute_delay`](@ref) or not.
- `tord::Int = 5`: order of Taylor expansions.
- `niter::Int = 10`: number of light-time solution iterations.
- `xve::EarthEph = earthposvel`: Earth ephemeris.
- `xvs::SunEph = sunposvel`: Sun ephemeris.
- `xva::AstEph`: asteroid ephemeris.

All ephemeris must take [et seconds since J2000] and return [barycentric position in km
and velocity in km/sec].
"""
function residuals(obs::Vector{RadarJPL{T}}; kwargs...) where {T <: AbstractFloat}
    # Radar delay/Doppler astrometric data
    x_jt = radar_astrometry(obs; kwargs...)
    # Time-delay residuals
    res_τ = x_jt[2][x_jt[14]] .- x_jt[4][x_jt[14]]
    # Doppler-shift residuals
    res_ν = x_jt[3][x_jt[15]] .- x_jt[5][x_jt[15]]
    # Weights
    w_τ = 1 ./ (x_jt[6][x_jt[14]].^2)
    w_ν = 1 ./ (x_jt[7][x_jt[15]].^2)

    return res_τ, w_τ, res_ν, w_ν
end
