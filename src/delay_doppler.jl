# construct path of JPL ephemerides
const jplephpath = joinpath(dirname(pathof(Apophis)), "../jpleph")

# read JPL ephemerides (Apophis, Solar System, TT-TDB)
function loadjpleph()
    furnsh.( joinpath.(jplephpath, ["a99942.bsp", "de430_1850-2150.bsp", "TTmTDB.de430.19feb2015.bsp"]) )
end

# this is an auxiliary function which converts a [x,y,z,vx,vy,vz] "state" vector from km,km/sec units to au,au/day
function kmsec2auday(pv)
    pv /= au # (km, km/sec) -> (au, au/sec)
    pv[4:6] *= daysec # (au, au/sec) -> (au, au/day)
    return pv
end

# get [x,y,z,vx,vy,vz] geometric "state" vector at TDB instant `et` from
# SPK-formatted ephemeris file wrt J2000 frame
function getpv(target::Int, observer::Int, et)
    return spkgeo(target, et, "J2000", observer)[1] # units: km,km/sec
end

# NAIF IDs:
#0: Solar System Barycenter
#10: Sun (heliocenter)
#2099942: Apophis
#399: Earth (geocenter)
#301: Moon
#1000000001 from body 1000000000: TT-TDB
# Here, we follow the convention from the CSPICE, library, that the ephemeris
# time is referred to the J2000 frame epoch:
# https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html#Terminology
# argument `et` represents "ephemeris seconds" (TDB seconds) since J2000.0 TDB epoch
# position and velocity are assumed to be returned in km, km/sec, resp., by spkgeo
# https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkgeo_c.html (see: Detailed ouput section)
apophis_pv(et) = kmsec2auday( getpv(2099942, 0, et) ) # units: au, au/day
sun_pv(et) = kmsec2auday( getpv(10, 0, et) ) # units: au, au/day
earth_pv(et) = kmsec2auday( getpv(399, 0, et) ) # units: au, au/day
moon_pv(et) = kmsec2auday( getpv(301, 0, et) ) # units: au, au/day
tt_m_tdb(et) = getpv(1000000001, 1000000000, et) # units: seconds

# Convert julian days to ephemeris seconds since J2000
function julian2etsecs(jd)
    return (jd-J2000)*daysec
end

# Convert ephemeris seconds since J2000 to julian days
function etsecs2julian(et)
    return J2000 + et/daysec
end

# Standard formula for relativistic (Shapiro) delay
function shapiro_delay(e, p, q)
    shap_del_days = 2μ[1]*log( abs((e+p+q)/(e+p-q)) )/(c_au_per_day^3) # days
    return shap_del_days*daysec # seconds
end

function shapiro_doppler(e, de, p, dp, q, dq, F_tx)
    # shap_del_diff = 2μ[1]*( (de+dp+dq)/(e+p+q) - (de+dp-dq)/(e+p-q) )/(c_au_per_day^3) # (adim.)
    shap_del_diff = 4μ[1]*(  ( dq*(e+p) - q*(de+dp) )/( (e+p)^2 - q^2 )  )/(c_au_per_day^3) # differential of Shapiro delay (adim.)
    shap_dop = -F_tx*shap_del_diff # ν = -F_tx*dτ/dt (units of F_tx) <-- Shapiro, Ash, Tausner (1966), footnote 10
    return shap_dop # (units of F_tx)
end

# Density of ionized electrons in interplanetary medium (ESAA 2014, p. 323, Sec. 8.7.5, Eq. 8.22)
# ESAA 2014 text probably has an error: Maybe it should say that transmitter frequency is in MHz instead of Hz
# But it doesn't specify output units of correction to time delay either
# Currently there is no errata reported about this (Apr 3, 2019)
# Errata URL: (https://aa.usno.navy.mil/publications/docs/exp_supp_errata.pdf)
# ESAA 2014 in turn refers to Muhleman and Anderson (1981)
# Ostro (1993) gives a reference to Anderson (1978), where this model is fitted to Mariner 9 ranging data
# Reading https://gssc.esa.int/navipedia/index.php/Ionospheric_Delay
# Helped a lot to clarify things, especially the 40.3, although they talk about Earth's ionosphere
# Another valuable source is Standish, E.M., Astron. Astrophys. 233, 252-271 (1990)
# Ne: ionized electrons density: electrons/cm^3
# p1: signal departure point (transmitter/bounce) (au)
# p2: signal arrival point (bounce/receiver) (au)
# r_s_t0, v_s_t0: Barycentric position (au) and velocity (au/day) of Sun at initial time of propagation of signal path (bounce time for downlink; transmit time for uplink)
# ds: current distance travelled by ray from emission point (au)
# ΔS: total distance between p1 and p2 (au)
function Ne(p1::Vector{S}, p2::Vector{S}, r_s_t0::Vector{S}, ds::U, ΔS::Real) where {S<:Number, U<:Number}
    # s: linear parametrization of ray path, such that
    # s=0 -> point on ray path is at p1; s=1 -> point on ray path is at p2
    # @show p1, p2, r_s_t0
    s = ds/ΔS
    # @show s
    s_p2_p1 = map(x->s*x, Taylor1.(p2-p1, s.order))
    r_vec = Taylor1.(p1, s.order) + s_p2_p1 - Taylor1.(r_s_t0, s.order) # heliocentric position (au) of point on ray path at time t_tdb_jul (Julian days)
    # @show r_vec
    r = sqrt( r_vec[1]^2 + r_vec[2]^2 + r_vec[3]^2 ) # heliocentric distance (au) of point on ray path at time t_tdb_jul (Julian days)
    # compute heliocentric position vector of point on ray path wrt Sun's rotation pole and equator (i.e., heliographic)
    α_p_sun_rad = deg2rad(α_p_sun)
    δ_p_sun_rad = deg2rad(δ_p_sun)
    r_vec_heliographic = inv( pole_rotation(α_p_sun_rad, δ_p_sun_rad) )*r_vec
    # compute heliographic solar latitude (Anderson, 1978) of point on ray path
    β = asin( r_vec_heliographic[3]/r ) # ecliptic solar latitude of point on ray path (rad)
    # β.order == 1 && @show β
    r_sr = r/R_sun
    # r_sr.order == 2 && @show r_sr
    Ne_t1 = (A_sun/r_sr^6)
    Ne_t2 = ( (a_sun*b_sun)/sqrt((a_sun*sin(β))^2 + (b_sun*cos(β))^2) )/(r_sr^2)
    Ne_val = Ne_t1 + Ne_t2
    # @show Ne_val
    return Ne_val
end

# Integral of Ne, evaluated with TaylorIntegration.jl
# TODO: @taylorize!
# p1: signal departure point (transmitter/bounce) (au)
# p2: signal arrival point (bounce/receiver) (au)
# t_tdb_jd1, t_tdb_jd2: Two-part Julian date (TDB) of signal path (bounce time for downlink; transmit time for uplink)
# current_delay: total time-delay for current signal path (TDB days)
# output has units of electrons/cm^2
function Ne_path_integral(p1::Vector{S}, p2::Vector{S}, r_s_t0::Vector{S}) where {S<:Number}
    ΔS = (100_000au)*norm(p2-p1) # total distance between p1 and p2, in centimeters
    # kernel of path integral; distance parameter `s` and total distance `ΔS` is in cm
    function int_kernel(s, x)
        return Ne(p1, p2, r_s_t0, s, ΔS)
    end
    # do path integral
    i0 = zero(p1[1])
    tT = Taylor1(24)
    x0T = Taylor1(i0, 24)
    # @show i0, tT, x0T
    TaylorIntegration.jetcoeffs!(int_kernel, tT, x0T)
    # @show x0T, x0T(ΔS)
    return x0T(ΔS)
end

# Empirical numerical factor correction to solar corona delay
# F_tx: transmitter frequency (MHz)
function corona_k(F_tx::T) where {T<:Real}
    if issband(F_tx)
        return 5
    elseif isxband(F_tx)
        return 400
    else
        @warn "Frequency is neither X or S band. Returning k=1"
        return 1
    end
end

# Time-delay generated by thin plasma of solar corona
# Ne: ionized electrons density: electrons/cm^3
# p1: signal departure point (transmitter/bounce for up/down-link, resp.) (au)
# p2: signal arrival point (bounce/receiver for up/down-link, resp.) (au)
# t_tdb_jd1, t_tdb_jd2: Two-part Julian date (TDB) of signal path (bounce time for downlink; transmit time for uplink)
# F_tx: transmitter frequency (MHz)
# From https://gssc.esa.int/navipedia/index.php/Ionospheric_Delay it seems that
# the expression 40.3*Ne/f^2 is adimensional, where Ne [electrons/centimeters^3] and f is in Hz
# therefore, the integral (40.3/f^2)*∫Ne*ds is in centimeters, where ds is in centimeters
# and the expression (40.3/(c*f^2))*∫Ne*ds, with c in centimeters/second, is in seconds
function corona_delay(p1::Vector{S}, p2::Vector{S}, r_s_t0::Vector{S}, F_tx::U, station_code::Int) where {S<:Number, U<:Real}
    # for the time being, we're removing the terms associated with higher-order terms in the variationals (ie, Yarkovsky)
    # int_path = Ne_path_integral(map(x->constant_term.(x), (p1, p2, r_s_t0))...) # (electrons/cm^2)
    int_path = Ne_path_integral(p1, p2, r_s_t0) # (electrons/cm^2)
    Δτ_corona = corona_k(F_tx)*40.3e-6int_path/(c_cm_per_sec*(F_tx)^2) # seconds
    return Δτ_corona # seconds
end

# *VERY* elementary computation of zenith distance
# r_antenna: position of antenna at receive/transmit time in celestial frame wrt geocenter
# ρ_vec_ae: slant-range vector from antenna to asteroid
function zenith_distance(r_antenna::Vector{T}, ρ_vec_ae::Vector{S}) where {T<:Number, S<:Number}
    norm_r_antenna = sqrt(r_antenna[1]^2 + r_antenna[2]^2 + r_antenna[3]^2)
    norm_ρ_vec_ae = sqrt(ρ_vec_ae[1]^2 + ρ_vec_ae[2]^2 + ρ_vec_ae[3]^2)
    cos_antenna_slant = dot(r_antenna, ρ_vec_ae)/(norm_r_antenna*norm_ρ_vec_ae)
    return acos(cos_antenna_slant)
end

# Time delay from the Earth's troposphere for radio frequencies
# oscillates between 0.007μs and 0.225μs for z between 0 and π/2 rad
# z: zenith distance (radians)
tropo_delay(z) = (7e-9)/( cos(z) + 0.0014/(0.045+cot(z)) ) # seconds

function tropo_delay(r_antenna::Vector{T}, ρ_vec_ae::Vector{S}) where {T<:Number, S<:Number}
    zd = zenith_distance(r_antenna, ρ_vec_ae)
    # @show rad2deg(zd)
    return tropo_delay(zd) # seconds
end

# Alternate version of delay_doppler, using JPL ephemerides
# station_code: observing station identifier (MPC nomenclature)
# t_r_utc_julian: time of echo reception (UTC)
# F_tx: transmitter frequency (MHz)
# xve: Earth ephemeris wich takes et seconds since J2000 as input and returns Earth barycentric position in au and velocity in au/day
# xvs: Sun ephemeris wich takes et seconds since J2000 as input and returns Sun barycentric position in au and velocity in au/day
# xva: Apophis ephemeris wich takes et seconds since J2000 as input and returns Apophis barycentric position in au and velocity in au/day
function delay_doppler(station_code::Int, t_r_utc::DateTime, F_tx::Real,
        niter::Int=10; xve::Function=earth_pv, xvs::Function=sun_pv,
        xva::Function=apophis_pv)
    # Transform receiving time from UTC to TDB seconds since j2000
    et_r_secs = seconds( AstroTime.j2000(TDB, UTCEpoch(t_r_utc)) ).Δt
    # Compute geocentric position/velocity of receiving antenna in inertial frame (au, au/day)
    R_r, V_r = observer_position(station_code, t_r_utc)
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
    vr_s_t_r = xvs(et_r_secs)
    r_s_t_r = vr_s_t_r[1:3]
    v_s_t_r = vr_s_t_r[4:6]
    # down-leg iteration
    # τ_D first approximation: Eq. (1) Yeomans et al. (1992)
    ρ_vec_r = r_a_t_r - r_r_t_r
    ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)
    τ_D = ρ_r/c_au_per_sec # (seconds) -R_b/c, but delay is wrt asteroid Center (Brozovic et al., 2018)
    # bounce time, 1st estimate, Eq. (2) Yeomans et al. (1992)
    et_b_secs = et_r_secs - τ_D
    # asteroid barycentric position (in au) at bounce time (TDB)
    rv_a_t_b = xva(et_b_secs)
    r_a_t_b = rv_a_t_b[1:3]
    v_a_t_b = rv_a_t_b[4:6]

    Δτ_D = zero(τ_D)
    Δτ_rel_D = zero(τ_D)
    Δν_rel_D = zero(τ_D)
    Δτ_corona_D = zero(τ_D)
    Δτ_tropo_D = zero(τ_D)
    # Δν_corona_D = zero(τ_D)
    Δν_tropo_D = zero(τ_D)

    for i in 1:niter
        # Eq. (3) Yeomans et al. (1992)
        ρ_vec_r = r_a_t_b - r_r_t_r
        ρ_vec_dot_r = v_a_t_b - v_r_t_r
        # Eq. (4) Yeomans et al. (1992)
        ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)
        τ_D = ρ_r/c_au_per_sec # (seconds) -R_b/c (COM correction) + Δτ_D (relativistic, tropo, iono...)

        # compute down-leg Shapiro delay
        # NOTE: when using PPN, substitute 2 -> 1+γ in expressions for Shapiro delay, Δτ_rel_[D|U]
        r_rs_t_r = r_r_t_r-r_s_t_r
        v_rs_t_r = v_r_t_r-v_s_t_r
        r_rs_t_rT = r_rs_t_r + map(x->x*Taylor1(1), v_rs_t_r)
        e_DT = sqrt(r_rs_t_rT[1]^2 + r_rs_t_rT[2]^2 + r_rs_t_rT[3]^2)
        e_D = e_DT[0] # heliocentric distance of station at t_r
        de_D = e_DT[1]
        rv_s_t_b = xvs(et_b_secs) # barycentric position and velocity of Sun at estimated bounce time
        r_s_t_b = rv_s_t_b[1:3]
        v_s_t_b = rv_s_t_b[4:6]
        r_as_t_b = constant_term.(r_a_t_b-r_s_t_b)
        v_as_t_b = constant_term.(v_a_t_b-v_s_t_b)
        r_as_t_bT = r_as_t_b + map(x->x*Taylor1(1), v_as_t_b)
        p_DT = sqrt(r_as_t_bT[1]^2 + r_as_t_bT[2]^2 + r_as_t_bT[3]^2)
        p_D = p_DT[0] # heliocentric distance of asteroid at t_b
        dp_D = p_DT[1]
        r_ae_t_b = constant_term.(ρ_vec_r)
        v_ae_t_b = constant_term.(ρ_vec_dot_r)
        r_ae_t_bT = r_ae_t_b + map(x->x*Taylor1(1), v_ae_t_b)
        q_DT = sqrt(r_ae_t_bT[1]^2 + r_ae_t_bT[2]^2 + r_ae_t_bT[3]^2)
        q_D = q_DT[0] #signal path distance (down-leg)
        dq_D = q_DT[1]
        Δτ_rel_D_T = shapiro_delay(e_DT, p_DT, q_DT)
        # Δτ_rel_D = shapiro_delay(e_D, p_D, q_D) # seconds
        # @show shapiro_delay(e_DT, p_DT, q_DT)
        Δτ_rel_D = Δτ_rel_D_T[0] # seconds
        # @show Δτ_rel_D
        # Δν_rel_D = shapiro_doppler(e_D, de_D, p_D, dp_D, q_D, dq_D, F_tx) # MHz
        Δν_rel_D = -F_tx*Δτ_rel_D_T[1]/daysec # MHz
        # @show Δν_rel_D
        # r_a_t_bT = constant_term.(r_a_t_b) + map(x->x*Taylor1(1), constant_term.(v_a_t_b))
        # r_r_t_rT = constant_term.(r_r_t_r) + map(x->x*Taylor1(1), constant_term.(v_r_t_r))
        # r_s_t_rT = constant_term.(r_s_t_r) + map(x->x*Taylor1(1), constant_term.(v_s_t_r))
        # Δτ_corona_DT = corona_delay(r_a_t_bT, r_r_t_rT, r_s_t_rT, F_tx, station_code) # seconds
        # Δτ_corona_D = Δτ_corona_DT[0] # seconds
        Δτ_corona_D = corona_delay(constant_term.(r_a_t_b), r_r_t_r, r_s_t_r, F_tx, station_code) # seconds
        # @show Δτ_corona_D
        # @show Δν_corona_D
        # # R_rT = R_r + map(x->x*Taylor1(1), V_r)
        # # ρ_vec_rT = ρ_vec_r + map(x->x*Taylor1(1), ρ_vec_dot_r)
        # # @show R_rT, ρ_vec_rT
        # # Δν_tropo_D = F_tx*daysec*tropo_delay(R_rT, ρ_vec_rT)[1]
        # # @show Δν_tropo_D
        # Δτ_tropo_D = tropo_delay(constant_term.(R_r), constant_term.(ρ_vec_r)) # -F_tx*Δτ_corona_DT[1]/daysec
        R_rT = constant_term.(R_r) + map(x->x*Taylor1(1), constant_term.(V_r))
        ρ_vec_rT = constant_term.(ρ_vec_r) + map(x->x*Taylor1(1), constant_term.(ρ_vec_dot_r))
        Δτ_tropo_DT = tropo_delay(R_rT, ρ_vec_rT)
        Δτ_tropo_D = Δτ_tropo_DT[0] # seconds
        Δν_tropo_D = -F_tx*Δτ_tropo_DT[1]/daysec # MHz
        # @show Δτ_tropo_D
        # @show Δτ_tropo_DT
        Δτ_D = Δτ_rel_D + Δτ_corona_D + Δτ_tropo_D # seconds TODO: convert formulas to seconds
        # # @show τ_D
        # # @show Δτ_rel_D
        # # @show Δτ_corona_D
        # # @show Δτ_tropo_D

        # bounce time, new estimate Eq. (2) Yeomans et al. (1992)
        et_b_secs = et_r_secs - τ_D
        # asteroid barycentric position (in au) at bounce time (TDB)
        rv_a_t_b = xva(et_b_secs)
        r_a_t_b = rv_a_t_b[1:3]
        v_a_t_b = rv_a_t_b[4:6]
    end
    τ_D = τ_D + Δτ_D
    # @show Δν_corona_D
    # get latest estimates of ρ_vec_r and ρ_r
    # Eq. (3) Yeomans et al. (1992)
    ρ_vec_r = r_a_t_b - r_r_t_r
    # Eq. (4) Yeomans et al. (1992)
    ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)

    # up-leg iteration
    # τ_U first estimation: Eq. (5) Yeomans et al. (1992)
    τ_U = τ_D
    # @show τ_U
    # transmit time, 1st estimate Eq. (6) Yeomans et al. (1992)
    et_t_secs = et_r_secs - (τ_U+τ_D)
    # convert transmit time from TDB to UTC time
    t_t_utc = DateTime(UTCEpoch(TDBEpoch(constant_term(et_t_secs)/daysec, origin=:j2000)))
    # Geocentric position and velocity of transmitting antenna in inertial frame (au, au/day)
    R_t, V_t = observer_position(station_code, t_t_utc)
    rv_e_t_t = xve(et_t_secs)
    r_e_t_t = rv_e_t_t[1:3]
    v_e_t_t = rv_e_t_t[4:6]
    # Transmitter barycentric position and velocity of at transmit time
    r_t_t_t = r_e_t_t + R_t
    v_t_t_t = v_e_t_t + V_t
    # Eq. (7) Yeomans et al. (1992)
    ρ_vec_t = r_a_t_b - r_t_t_t
    ρ_vec_dot_t = v_a_t_b - v_t_t_t
    ρ_t = sqrt(ρ_vec_t[1]^2 + ρ_vec_t[2]^2 + ρ_vec_t[3]^2)

    Δτ_U = zero(τ_U)
    Δτ_rel_U = zero(τ_U)
    Δν_rel_U = zero(τ_U)
    Δτ_corona_U = zero(τ_U)
    Δτ_tropo_U = zero(τ_U)
    # Δν_corona_U = zero(τ_U)
    Δν_tropo_U = zero(τ_U)

    for i in 1:niter
        # Eq. (8) Yeomans et al. (1992)
        τ_U = ρ_t/c_au_per_sec # (seconds) -R_b/c (COM correction) + Δτ_U (relativistic, tropo, iono...)

        # compute up-leg Shapiro delay
        # Sun barycentric position and velocity (in au, au/day) at transmit time (TDB)
        rv_s_t_t = xvs( et_r_secs-(τ_U+τ_D) )
        r_s_t_t = rv_s_t_t[1:3]
        v_s_t_t = rv_s_t_t[4:6]
        r_ts_t_t = constant_term.(r_t_t_t - r_s_t_t)
        v_ts_t_t = constant_term.(v_t_t_t - v_s_t_t)
        r_ts_t_tT = r_ts_t_t + map(x->x*Taylor1(1), v_ts_t_t)
        e_UT = sqrt(r_ts_t_tT[1]^2 + r_ts_t_tT[2]^2 + r_ts_t_tT[3]^2)
        e_U = e_UT[0] # heliocentric distance of station at t_t
        de_U = e_UT[1]
        rv_s_t_b = xvs(et_b_secs) # barycentric position/velocity of Sun at bounce time
        r_s_t_b = rv_s_t_b[1:3]
        v_s_t_b = rv_s_t_b[4:6]
        r_as_t_b = constant_term.(r_a_t_b-r_s_t_b)
        v_as_t_b = constant_term.(v_a_t_b-v_s_t_b)
        r_as_t_bT = r_as_t_b + map(x->x*Taylor1(1), v_as_t_b)
        p_UT = sqrt(r_as_t_bT[1]^2 + r_as_t_bT[2]^2 + r_as_t_bT[3]^2)
        p_U = p_UT[0] # heliocentric distance of asteroid at t_b
        dp_U = p_UT[1]
        # rv_e_t_b = xve(et_b_secs) # barycentric position/velocity of Earth at bounce time
        # r_e_t_b = rv_e_t_b[1:3]
        # v_e_t_b = rv_e_t_b[4:6]
        r_ae_t_b = constant_term.(r_a_t_b-r_t_t_t)
        v_ae_t_b = constant_term.(v_a_t_b-v_t_t_t)
        r_ae_t_bT = r_ae_t_b + map(x->x*Taylor1(1), v_ae_t_b)
        q_UT = sqrt(r_ae_t_bT[1]^2 + r_ae_t_bT[2]^2 + r_ae_t_bT[3]^2)
        q_U = q_UT[0] #signal path (up-leg)
        dq_U = q_UT[1]
        Δτ_rel_U_T = shapiro_delay(e_UT, p_UT, q_UT)
        # Δτ_rel_U = shapiro_delay(e_U, p_U, q_U) # seconds
        # @show shapiro_delay(e_UT, p_UT, q_UT)
        Δτ_rel_U = Δτ_rel_U_T[0] # seconds
        # @show Δτ_rel_U
        # Δν_rel_U = shapiro_doppler(e_U, de_U, p_U, dp_U, q_U, dq_U, F_tx)
        Δν_rel_U = -F_tx*Δτ_rel_U_T[1]/daysec # MHz
        # @show Δν_rel_U
        # r_t_t_tT = constant_term.(r_t_t_t) + map(x->x*Taylor1(1), constant_term.(v_t_t_t))
        # r_a_t_bT = constant_term.(r_a_t_b) + map(x->x*Taylor1(1), constant_term.(v_a_t_b))
        # r_s_t_bT = constant_term.(r_s_t_b) + map(x->x*Taylor1(1), constant_term.(v_s_t_b))
        # Δτ_corona_UT = corona_delay(r_t_t_tT, r_a_t_bT, r_s_t_bT, F_tx, station_code) # seconds
        Δτ_corona_U = corona_delay(constant_term.(r_t_t_t), constant_term.(r_a_t_b), constant_term.(r_s_t_b), F_tx, station_code) # seconds
        # Δτ_corona_U = Δτ_corona_UT[0] # seconds
        #@show Δτ_corona_UT
        # @show Δτ_corona_U
        # @show Δν_corona_U
        # # R_tT = R_t + map(x->x*Taylor1(1), V_t)
        # # ρ_vec_tT = ρ_vec_t + map(x->x*Taylor1(1), ρ_vec_dot_t)
        # # @show R_tT, ρ_vec_tT
        # # Δν_tropo_U = F_tx*daysec*tropo_delay(R_tT, ρ_vec_tT)[1]
        # # @show Δν_tropo_U
        # Δτ_tropo_U = tropo_delay(constant_term.(R_t), constant_term.(ρ_vec_t)) # seconds # -F_tx*Δτ_corona_UT[1]/daysec
        R_tT = constant_term.(R_t) + map(x->x*Taylor1(1), constant_term.(V_t))
        ρ_vec_tT = constant_term.(ρ_vec_t) + map(x->x*Taylor1(1), constant_term.(ρ_vec_dot_t))
        Δτ_tropo_UT = tropo_delay(R_tT, ρ_vec_tT)
        Δτ_tropo_U = Δτ_tropo_UT[0] #seconds
        Δν_tropo_U = -F_tx*Δτ_tropo_UT[1]/daysec # MHz
        # @show Δτ_tropo_U
        # @show Δτ_tropo_UT
        Δτ_U = Δτ_rel_U + Δτ_corona_U + Δτ_tropo_U
        # # @show τ_U
        # # @show Δτ_rel_U
        # # @show Δτ_corona_U
        # # @show Δτ_tropo_U
        # # transmit time, 1st estimate Eq. (6) Yeomans et al. (1992)

        et_t_secs = et_r_secs-(τ_U+τ_D)
        # convert transmit time to UTC time
        # TDB -> UTC
        t_t_utc = DateTime(UTCEpoch(TDBEpoch(constant_term(et_t_secs)/daysec, origin=:j2000)))
        # Geocentric position/velocity of receiving antenna in inertial frame (au, au/day)
        R_t, V_t = observer_position(station_code, t_t_utc)
        # Earth's barycentric position and velocity at the transmit time
        rv_e_t_t = xve(et_t_secs)
        r_e_t_t = rv_e_t_t[1:3]
        v_e_t_t = rv_e_t_t[4:6]
        # Barycentric position and velocity of the transmitter at the transmit time
        r_t_t_t = r_e_t_t + R_t
        v_t_t_t = v_e_t_t + V_t
        # Eq. (7) Yeomans et al. (1992)
        ρ_vec_t = r_a_t_b - r_t_t_t
        ρ_vec_dot_t = v_a_t_b - v_t_t_t
        ρ_t = sqrt(ρ_vec_t[1]^2 + ρ_vec_t[2]^2 + ρ_vec_t[3]^2)
    end
    # Sun barycentric position (in au) at transmit time (TDB)
    r_s_t_t = xvs(et_t_secs)[1:3]
    τ_U = τ_U + Δτ_U

    # compute TDB-UTC at transmit time
    # transmit time: TDB -> UTC
    t_t_utc = UTCEpoch(TDBEpoch(constant_term(et_t_secs)/daysec, origin=:j2000))
    # ut1_utc_t = EarthOrientation.getΔUT1( julian(t_t_utc).Δt ) # seconds
    utc_jd1_t, utc_jd2_t = julian_twopart(t_t_utc) # days, days
    # @show utc_jd2_t.Δt
    j, tai_utc_t = iauDat(year(t_t_utc), month(t_t_utc), day(t_t_utc), utc_jd2_t.Δt)
    j != 0 && @warn "iauDat: j = $j. See SOFA.jl docs for more detail."
    tt_utc_t = 32.184 + tai_utc_t # seconds
    ##TT-TDB (transmit time)
    tt_tdb_t = tt_m_tdb(constant_term(et_t_secs))[1] # seconds
    tdb_utc_t = tt_utc_t - tt_tdb_t #seconds
    # @show tdb_utc_t

    # compute TDB-UTC at receive time
    # UTC two-part "quasi" Julian date (see http://www.iausofa.org/sofa_ts_c.pdf)
    utc_jd1_r, utc_jd2_r = julian_twopart(UTCEpoch(t_r_utc)) # days, days
    # ΔAT = TAI - UTC
    j, tai_utc_r = iauDat(year(t_r_utc), month(t_r_utc), day(t_r_utc), utc_jd2_r.Δt)
    j != 0 && @warn "iauDat: j = $j. See SOFA.jl docs for more detail."
    tt_utc_r = 32.184 + tai_utc_r
    # compute `tt_tdb_r` recursively
    tt_tdb_r = 0.0
    for i in 1:niter+1 #TODO: check that convergence always holds
        utc_j2000_r = seconds( AstroTime.j2000(UTCEpoch(t_r_utc)) ).Δt # days since J2000 (UTC)
        tt_tdb_r = tt_m_tdb(utc_j2000_r + (tt_utc_r + tt_tdb_r))[1] # seconds
        # @show tt_tdb_r
    end
    tdb_utc_r = tt_utc_r - tt_tdb_r #seconds
    # @show tdb_utc_r

    # compute total time delay (UTC seconds)
    # Eq. (9) Yeomans et al. (1992)
    τ = (τ_U + τ_D) + (tdb_utc_t - tdb_utc_r)

    # compute Doppler shift ν
    # Eq. (10) Yeomans et al. (1992)
    ρ_vec_dot_t = v_a_t_b - v_t_t_t
    ρ_vec_dot_r = v_a_t_b - v_r_t_r
    # Eq. (11) Yeomans et al. (1992)
    ρ_dot_t = dot(ρ_vec_t, ρ_vec_dot_t)/ρ_t
    ρ_dot_r = dot(ρ_vec_r, ρ_vec_dot_r)/ρ_r
    # Eq. (12) Yeomans et al. (1992)
    doppler_c = -F_tx*(ρ_dot_t+ρ_dot_r)/c_au_per_day
    p_t = dot(ρ_vec_t, v_t_t_t)/ρ_t
    p_r = dot(ρ_vec_r, v_r_t_r)/ρ_r
    doppler_c2_t1 = ρ_dot_t*p_t - ρ_dot_r*p_r - ρ_dot_t*ρ_dot_r # order c^-2, 1st term
    r_ts = r_t_t_t - r_s_t_t
    r_rs = r_r_t_r - r_s_t_r
    ϕ1 = μ[1]/sqrt(r_ts[1]^2+r_ts[2]^2+r_ts[3]^2)
    ϕ3 = μ[1]/sqrt(r_rs[1]^2+r_rs[2]^2+r_rs[3]^2)
    doppler_c2_t2 = ϕ1 - ϕ3 # order c^-2, 2nd term
    doppler_c2_t3 = (  dot(v_t_t_t, v_t_t_t) - dot(v_r_t_r, v_r_t_r)  )/2  # order c^-2, 3rd term
    doppler_c2 = -F_tx*(doppler_c2_t1 + doppler_c2_t2 + doppler_c2_t3)/(c_au_per_day^2)
    # Corrections of order c^-3 to Doppler shift (Moyer, 1971, p. 56, Eq. 343)
    doppler_c3_t1 = ρ_dot_t*p_t^2 - ρ_dot_r*p_r^2 - ρ_dot_t*ρ_dot_r*(p_r-p_t)
    doppler_c3_t2 = -(ρ_dot_t+ρ_dot_r)*(doppler_c2_t2+doppler_c2_t3)
    # Note that c^3 corrections to Doppler include a "Shapiro-like" term, which has already been multiplied by the transmitter frequency
    doppler_c3 = -F_tx*(doppler_c3_t1+doppler_c3_t2)/(c_au_per_day^3) + (Δν_rel_D+Δν_rel_U)
    ν = (doppler_c + doppler_c2) + doppler_c3  + (Δν_tropo_D+Δν_tropo_U)

    return 1e6τ, 1e6ν, (Δτ_corona_D+Δτ_corona_U) # total signal delay (μs) and Doppler shift (Hz)
end
