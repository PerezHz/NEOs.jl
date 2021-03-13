# load ttmtdb as a TaylorInterpolant saved in .jld file
const jldephpath = joinpath(pkgdir(Apophis), "jldeph")
const ttmtdb = JLD.load(joinpath(jldephpath, "ttmtdb_DE430_2003_2030.jld"), "ttmtdb")

# read JPL ephemerides (Apophis, Solar System, TT-TDB)
function loadjpleph()
    furnsh(
        joinpath(artifact"naif0012", "naif0012.tls"),
        joinpath(artifact"TTmTDBde430", "TTmTDB.de430.19feb2015.bsp"),
        joinpath(artifact"de430", "de430_1850-2150.bsp"),
        joinpath(artifact"a99942", "a99942_s197.bsp"),
        joinpath(artifact"a99942", "a99942_s199.bsp"),
    )
end

# this is an auxiliary function which converts a [x,y,z,vx,vy,vz] "state" vector from km,km/sec units to au,au/day
function kmsec2auday(pv)
    pv /= au # (km, km/sec) -> (au, au/sec)
    pv[4:6] *= daysec # (au, au/sec) -> (au, au/day)
    return pv
end

# this is an auxiliary function which converts a [x,y,z,vx,vy,vz] "state" vector from au,au/day units to km,km/sec
function auday2kmsec(pv)
    pv *= au # (au, au/day) -> (km, km/day)
    pv[4:6] /= daysec # (km, km/day) -> (km, km/sec)
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
apophis_pv_197(et) = getpv(9904406, 0, constant_term(et)) # units: km, km/second
apophis_pv_199(et) = getpv(2099942, 0, constant_term(et)) # units: km, km/second
sun_pv(et) = getpv(10, 0, constant_term(et)) # units: km, km/second
earth_pv(et) = getpv(399, 0, constant_term(et)) # units: km, km/second
moon_pv(et) = getpv(301, 0, constant_term(et)) # units: km, km/second
tt_tdb(et) = getpv(1000000001, 1000000000, constant_term(et))[1] # units: seconds
dtt_tdb(et) = getpv(1000000001, 1000000000, constant_term(et))[4] # units: seconds/seconds

# Convert julian days to ephemeris seconds since J2000
function julian2etsecs(jd)
    return (jd-JD_J2000)*daysec
end

# Convert ephemeris seconds since J2000 to julian days
function etsecs2julian(et)
    return JD_J2000 + et/daysec
end

# Standard formula for relativistic (Shapiro) delay
function shapiro_delay(e, p, q)
    shap = 0.0 #2μ[1]/(c_au_per_day^2)
    shap_del_days = (2PlanetaryEphemeris.μ[su]/(c_au_per_day^3))*log( (e+p+q+shap)/(e+p-q+shap) ) # days
    return shap_del_days*daysec # seconds
end

function shapiro_doppler(e, de, p, dp, q, dq, F_tx)
    # shap_del_diff = 2μ[1]*( (de+dp+dq)/(e+p+q) - (de+dp-dq)/(e+p-q) )/(c_au_per_day^3) # (adim.)
    shap_del_diff = (4PlanetaryEphemeris.μ[su]/(c_au_per_day^3))*(  ( dq*(e+p) - q*(de+dp) )/( (e+p)^2 - q^2 )  ) # differential of Shapiro delay (adim.)
    shap_dop = -F_tx*shap_del_diff # ν = -F_tx*dτ/dt (units of F_tx) <-- Shapiro, Ash, Tausner (1966), footnote 10
    return shap_dop # (units of F_tx)
end

# Density of ionized electrons in interplanetary medium (ESAA 2014, p. 323, Sec. 8.7.5, Eq. 8.22)
# ESAA 2014 in turn refers to Muhleman and Anderson (1981)
# Ostro (1993) gives a reference to Anderson (1978), where this model is fitted to Mariner 9 ranging data
# Reading https://gssc.esa.int/navipedia/index.php/Ionospheric_Delay
# Helped a lot to clarify things, especially the 40.3, although they talk about Earth's ionosphere
# Another valuable source is Standish, E.M., Astron. Astrophys. 233, 252-271 (1990)
# Ne: ionized electrons density: electrons/cm^3
# p1: signal departure point (transmitter/bounce) (au)
# p2: signal arrival point (bounce/receiver) (au)
# r_s_t0: Barycentric position (au) of Sun at initial time of propagation of signal path (bounce time for down-leg; transmit time for up-leg)
# ds: current distance travelled by ray from emission point (au)
# ΔS: total distance between p1 and p2 (au)
function Ne(p1::Vector{S}, p2::Vector{S}, r_s_t0::Vector{S}, ds::U, ΔS::Real) where {S<:Number, U<:Number}
    # s: linear parametrization of ray path, such that
    # s=0 -> point on ray path is at p1; s=1 -> point on ray path is at p2
    s = ds/ΔS
    s_p2_p1 = map(x->s*x, Taylor1.(p2-p1, s.order))
    r_vec = Taylor1.(p1, s.order) + s_p2_p1 - Taylor1.(r_s_t0, s.order) # heliocentric position (au) of point on ray path at time t_tdb_jul (Julian days)
    r = sqrt( r_vec[1]^2 + r_vec[2]^2 + r_vec[3]^2 ) # heliocentric distance (au) of point on ray path at time t_tdb_jul (Julian days)
    # compute heliocentric position vector of point on ray path wrt Sun's rotation pole and equator (i.e., heliographic)
    α_p_sun_rad = deg2rad(α_p_sun)
    δ_p_sun_rad = deg2rad(δ_p_sun)
    r_vec_heliographic = inv( pole_rotation(α_p_sun_rad, δ_p_sun_rad) )*r_vec
    # compute heliographic solar latitude (Anderson, 1978) of point on ray path
    β = asin( r_vec_heliographic[3]/r ) # ecliptic solar latitude of point on ray path (rad)
    r_sr = r/R_sun
    Ne_t1 = (A_sun/r_sr^6)
    Ne_t2 = ( (a_sun*b_sun)/sqrt((a_sun*sin(β))^2 + (b_sun*cos(β))^2) )/(r_sr^2)
    Ne_val = Ne_t1 + Ne_t2
    return Ne_val
end

# Integral of Ne, evaluated with TaylorIntegration.jl
# TODO: @taylorize!
# p1: signal departure point (transmitter/bounce) (au)
# p2: signal arrival point (bounce/receiver) (au)
# r_s_t0: Barycentric position (au) of Sun at initial time of propagation of signal path (bounce time for down-leg; transmit time for up-leg)
# output is in units of electrons/cm^2
function Ne_path_integral(p1::Vector{S}, p2::Vector{S}, r_s_t0::Vector{S}) where {S<:Number}
    ΔS = (100_000au)*norm(p2-p1) # total distance between p1 and p2, in centimeters
    # kernel of path integral; distance parameter `s` and total distance `ΔS` is in cm
    function int_kernel(x, params, s)
        return Ne(p1, p2, r_s_t0, s, ΔS)
    end
    # do path integral
    i0 = zero(p1[1])
    tT = Taylor1(24)
    iT = Taylor1(i0, 24)
    TaylorIntegration.jetcoeffs!(int_kernel, tT, iT, nothing)
    return iT(ΔS)
end

# Time-delay due to thin plasma of solar corona
# Ne: ionized electrons density: electrons/cm^3
# p1: signal departure point (transmitter/bounce for up/down-link, resp.) (au)
# p2: signal arrival point (bounce/receiver for up/down-link, resp.) (au)
# t_tdb_jd1, t_tdb_jd2: Two-part Julian date (TDB) of signal path (bounce time for downlink; transmit time for uplink)
# F_tx: transmitter frequency (MHz)
# From https://gssc.esa.int/navipedia/index.php/Ionospheric_Delay it seems that
# ESAA 2014 text probably should say that in the formula for Δτ_corona,
# the expression 40.3*Ne/f^2 is adimensional, where Ne [electrons/centimeters^3] and f is in Hz
# therefore, the integral (40.3/f^2)*∫Ne*ds is in centimeters, where ds is in centimeters
# and the expression (40.3/(c*f^2))*∫Ne*ds, with c in centimeters/second, is in seconds
function corona_delay(p1::Vector{S}, p2::Vector{S}, r_s_t0::Vector{S}, F_tx::U, station_code::Int) where {S<:Number, U<:Real}
    # for the time being, we're removing the terms associated with higher-order terms in the variationals (ie, Yarkovsky)
    int_path = Ne_path_integral(p1, p2, r_s_t0) # (electrons/cm^2)
    Δτ_corona = 40.3e-6int_path/(c_cm_per_sec*(F_tx)^2) # seconds
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
    return tropo_delay(zd) # seconds
end

# auxiliary function to compute (TDB-UTC); useful to convert TDB to UTC via UTC + (TDB-UTC) and viceversa
# TDB: Solar System barycentric ephemeris time
# et: TDB seconds since J2000.0 (TDB)
# TT: Terrestrial time
# TAI: International Atomic Time
# UTC: Coordinated Universal Time
# TDB-UTC = (TDB-TAI) + (TAI-UTC)
#         = (TDB-TT) + (TT-TAI) + (TAI-UTC)
#         = (TDB-TT) + (32.184 s) + ΔAT
# does not include correction due to position of measurement station v_E*(r_S.r_E)/c^2 (Folkner et al. 2014; Moyer, 2003)
function tdb_utc(et::T) where {T<:Number}
    tt_tdb_et = ttmtdb(et)
    tt_tai = 32.184
    et_00 = constant_term(constant_term(et))
    utc_secs = et_00 - deltet(et_00, "ET") # used only to determine ΔAT; no high-precision needed
    jd_utc = JD_J2000 + utc_secs/daysec
    dt_utc = julian2datetime(jd_utc)
    fd_utc = (jd_utc+0.5) - floor(jd_utc+0.5)
    tai_utc = get_ΔAT(jd_utc)
    return (tt_tai + tai_utc) - tt_tdb_et # TDB-UTC = (TDB-TT) + (TT-TAI) + (TAI-UTC) = (TDB-TT) + 32.184 s + ΔAT
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

# Compute radar-astrometric round-trip time for an asteroid at
# UTC instant `t_r_utc` from tracking station with code `station_code`.
# station_code: observing station identifier (MPC nomenclature)
# t_r_utc: UTC time of echo reception (DateTime)
# t_offset: time offset wrt echo reception time, to compute Doppler shifts by range differences (seconds)
# niter: number of light-time solution iterations
# xve: Earth ephemeris wich takes TDB seconds since J2000 as input and returns Earth barycentric position in km and velocity in km/second
# xvs: Sun ephemeris wich takes TDB seconds since J2000 as input and returns Sun barycentric position in km and velocity in km/second
# xva: asteroid ephemeris wich takes TDB seconds since J2000 as input and returns asteroid barycentric position in km and velocity in km/second
function delay(station_code::Int, t_r_utc::DateTime, t_offset::Real,
        niter::Int=10; eo::Bool=true, xve=earth_pv, xvs=sun_pv,
        xva=apophis_pv_197)
    et_r_secs = str2et(string(t_r_utc)) + t_offset
    # Compute geocentric position/velocity of receiving antenna in inertial frame (au, au/day)
    R_r, V_r = observer_position(station_code, et_r_secs, eo=eo)
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
    # down-leg iteration
    # τ_D first approximation: Eq. (1) Yeomans et al. (1992)
    ρ_vec_r = r_a_t_r - r_r_t_r
    ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)
    τ_D = ρ_r/clightkms # (seconds) -R_b/c, but delay is wrt asteroid Center (Brozovic et al., 2018)
    # bounce time, new estimate Eq. (2) Yeomans et al. (1992)
    et_b_secs = et_r_secs - τ_D

    Δτ_D = zero(τ_D)
    Δτ_rel_D = zero(τ_D)
    # Δτ_corona_D = zero(τ_D)
    Δτ_tropo_D = zero(τ_D)

    for i in 1:niter
        # asteroid barycentric position (in au) at bounce time (TDB)
        rv_a_t_b = xva(et_b_secs)
        r_a_t_b = rv_a_t_b[1:3]
        v_a_t_b = rv_a_t_b[4:6]
        # Eq. (3) Yeomans et al. (1992)
        ρ_vec_r = r_a_t_b - r_r_t_r
        # Eq. (4) Yeomans et al. (1992)
        ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)
        # compute down-leg Shapiro delay
        # NOTE: when using PPN, substitute 2 -> 1+γ in expressions for Shapiro delay, Δτ_rel_[D|U]
        e_D_vec  = r_r_t_r - r_s_t_r
        e_D = sqrt(e_D_vec[1]^2 + e_D_vec[2]^2 + e_D_vec[3]^2) # heliocentric distance of Earth at t_r
        rv_s_t_b = xvs(et_b_secs) # barycentric position and velocity of Sun at estimated bounce time
        r_s_t_b = rv_s_t_b[1:3]
        p_D_vec  = r_a_t_b - r_s_t_b
        p_D = sqrt(p_D_vec[1]^2 + p_D_vec[2]^2 + p_D_vec[3]^2) # heliocentric distance of asteroid at t_b
        q_D = ρ_r #signal path distance (down-leg)
        # Shapiro correction to time-delay
        Δτ_rel_D = shapiro_delay(e_D, p_D, q_D)
        # troposphere correction to time-delay
        Δτ_tropo_D = tropo_delay(R_r, ρ_vec_r) # seconds
        # Δτ_corona_D = corona_delay(constant_term.(r_a_t_b), r_r_t_r, r_s_t_r, F_tx, station_code) # seconds
        Δτ_D = Δτ_rel_D # + Δτ_tropo_D #+ Δτ_corona_D # seconds
        p_dot_23 = dot(ρ_vec_r, v_a_t_b)/ρ_r
        Δt_2 = (τ_D - ρ_r/clightkms - Δτ_rel_D)/(1.0-p_dot_23/clightkms)
        τ_D = τ_D - Δt_2
        et_b_secs = et_r_secs - τ_D
    end
    rv_a_t_b = xva(et_b_secs)
    r_a_t_b = rv_a_t_b[1:3]
    v_a_t_b = rv_a_t_b[4:6]

    # up-leg iteration
    # τ_U first estimation: Eq. (5) Yeomans et al. (1992)
    τ_U = τ_D
    # transmit time, 1st estimate Eq. (6) Yeomans et al. (1992)
    et_t_secs = et_b_secs - τ_U
    # Geocentric position and velocity of transmitting antenna in inertial frame (au, au/day)
    R_t, V_t = observer_position(station_code, et_t_secs, eo=eo)
    rv_e_t_t = xve(et_t_secs)
    r_e_t_t = rv_e_t_t[1:3]
    v_e_t_t = rv_e_t_t[4:6]
    # Transmitter barycentric position and velocity of at transmit time
    r_t_t_t = r_e_t_t + R_t
    # Eq. (7) Yeomans et al. (1992)
    ρ_vec_t = r_a_t_b - r_t_t_t
    ρ_t = sqrt(ρ_vec_t[1]^2 + ρ_vec_t[2]^2 + ρ_vec_t[3]^2)

    Δτ_U = zero(τ_U)
    Δτ_rel_U = zero(τ_U)
    # Δτ_corona_U = zero(τ_U)
    Δτ_tropo_U = zero(τ_U)

    # println("   *** DOWNLEG LOOP ***")
    for i in 1:niter
        # Geocentric position and velocity of transmitting antenna in inertial frame (au, au/day)
        # TODO: remove `constant_term` to take into account dependency of R_t, V_t wrt initial conditions variations via et_t_secs
        R_t, V_t = observer_position(station_code, et_t_secs, eo=eo)
        # Earth's barycentric position and velocity at the transmit time
        rv_e_t_t = xve(et_t_secs)
        r_e_t_t = rv_e_t_t[1:3]
        v_e_t_t = rv_e_t_t[4:6]
        # Barycentric position and velocity of the transmitter at the transmit time
        r_t_t_t = r_e_t_t + R_t
        v_t_t_t = v_e_t_t + V_t
        # Eq. (7) Yeomans et al. (1992)
        ρ_vec_t = r_a_t_b - r_t_t_t
        ρ_t = sqrt(ρ_vec_t[1]^2 + ρ_vec_t[2]^2 + ρ_vec_t[3]^2)
        # compute up-leg Shapiro delay
        # Sun barycentric position and velocity (in au, au/day) at transmit time (TDB)
        rv_s_t_t = xvs(et_t_secs)
        r_s_t_t = rv_s_t_t[1:3]
        e_U_vec = r_t_t_t - r_s_t_t
        e_U = sqrt(e_U_vec[1]^2 + e_U_vec[2]^2 + e_U_vec[3]^2) # heliocentric distance of Earth at t_t
        rv_s_t_b = xvs(et_b_secs) # barycentric position/velocity of Sun at bounce time
        r_s_t_b = rv_s_t_b[1:3]
        p_U_vec = r_a_t_b - r_s_t_b
        p_U = sqrt(p_U_vec[1]^2 + p_U_vec[2]^2 + p_U_vec[3]^2) # heliocentric distance of asteroid at t_b
        q_U_vec = r_a_t_b - r_e_t_t
        q_U = ρ_t # signal path distance (up-leg)
        Δτ_rel_U = shapiro_delay(e_U, p_U, q_U) # seconds
        Δτ_tropo_U = tropo_delay(R_t, ρ_vec_t) # seconds
        # Δτ_corona_U = corona_delay(constant_term.(r_t_t_t), constant_term.(r_a_t_b), constant_term.(r_s_t_b), F_tx, station_code) # seconds
        Δτ_U = Δτ_rel_U # + Δτ_tropo_U #+ Δτ_corona_U # seconds
        p_dot_12 = -dot(ρ_vec_t, v_t_t_t)/ρ_t
        Δt_1 = (τ_U - ρ_t/clightkms - Δτ_rel_U)/(1.0-p_dot_12/clightkms)
        τ_U = τ_U - Δt_1
        # transmit time, new estimate
        et_t_secs = et_b_secs - τ_U
    end

    # compute TDB-UTC at transmit time
    # corrections to TT-TDB from Moyer (2003) / Folkner et al. (2014) due to position of measurement station on Earth are of order 0.01μs
    # Δtt_tdb_station_t = - dot(v_e_t_t, r_t_t_t-r_e_t_t)/clightkms^2
    tdb_utc_t = tdb_utc(et_t_secs) # + Δtt_tdb_station_t
    # compute TDB-UTC at receive time
    # corrections to TT-TDB from Moyer (2003) / Folkner et al. (2014) due to position of measurement station on Earth  are of order 0.01μs
    # Δtt_tdb_station_r = - dot(v_e_t_r, r_r_t_r-r_e_t_r)/clightkms^2
    tdb_utc_r = tdb_utc(et_r_secs) # + Δtt_tdb_station_r

    # compute total time delay (UTC seconds); relativistic delay is already included in τ_D, τ_U
    # Eq. (9) Yeomans et al. (1992)
    τ = (τ_D + τ_U) + (Δτ_tropo_D + Δτ_tropo_U) + (tdb_utc_t - tdb_utc_r) # seconds

    return 1e6τ # total signal delay (μs)
end

# Compute Taylor series expansion of time-delay observable around echo reception
# time. This allows to compute dopplers via autodiff using ν = -F_tx*dτ/dt
# NOTE: Works only with TaylorInterpolant ephemeris
function delay(station_code::Int, t_r_utc::DateTime,
        niter::Int=10; eo::Bool=true, xve::TaylorInterpolant=earth_pv,
        xvs::TaylorInterpolant=sun_pv, xva::TaylorInterpolant=apophis_pv_197,
        tord::Int=xva.x[1].order)
    q1 = xva.x[1] # auxiliary to evaluate JT ephemeris
    et_r_secs_0 = str2et(string(t_r_utc))
    et_r_secs = Taylor1([et_r_secs_0,1.0].*one(q1[0]), tord)
    # Compute geocentric position/velocity of receiving antenna in inertial frame (au, au/day)
    R_r, _ = observer_position(station_code, et_r_secs, eo=eo)
    # Earth's barycentric position and velocity at receive time
    r_e_t_r = xve(et_r_secs)[1:3]
    # Receiver barycentric position and velocity at receive time
    r_r_t_r = r_e_t_r + R_r
    # Asteroid barycentric position and velocity at receive time
    r_a_t_r = xva(et_r_secs)[1:3]
    # Sun barycentric position and velocity at receive time
    r_s_t_r = xvs(et_r_secs)[1:3]
    # down-leg iteration
    # τ_D first approximation: Eq. (1) Yeomans et al. (1992)
    ρ_vec_r = r_a_t_r - r_r_t_r
    ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)
    # @show ρ_r
    τ_D = ρ_r/clightkms # (seconds) -R_b/c, but delay is wrt asteroid Center (Brozovic et al., 2018)
    # @show τ_D
    # bounce time, new estimate Eq. (2) Yeomans et al. (1992)
    # @show et_r_secs
    et_b_secs = et_r_secs - τ_D
    # @show et_b_secs et_r_secs

    Δτ_D = zero(τ_D)
    Δτ_rel_D = zero(τ_D)
    # Δτ_corona_D = zero(τ_D)
    Δτ_tropo_D = zero(τ_D)
    # @show Δτ_D

    for i in 1:niter
        # asteroid barycentric position (in au) at bounce time (TDB)
        rv_a_t_b = xva(et_b_secs)
        r_a_t_b = rv_a_t_b[1:3]
        v_a_t_b = rv_a_t_b[4:6]
        # Eq. (3) Yeomans et al. (1992)
        ρ_vec_r = r_a_t_b - r_r_t_r
        # Eq. (4) Yeomans et al. (1992)
        ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)
        # @show ρ_r
        # compute down-leg Shapiro delay
        # NOTE: when using PPN, substitute 2 -> 1+γ in expressions for Shapiro delay, Δτ_rel_[D|U]
        e_D_vec  = r_r_t_r - r_s_t_r
        # @show e_D_vec
        e_D = sqrt(e_D_vec[1]^2 + e_D_vec[2]^2 + e_D_vec[3]^2) # heliocentric distance of Earth at t_r
        # barycentric position and velocity of Sun at estimated bounce time
        r_s_t_b = xvs(et_b_secs)[1:3]
        # @show r_s_t_b
        p_D_vec  = r_a_t_b - r_s_t_b
        p_D = sqrt(p_D_vec[1]^2 + p_D_vec[2]^2 + p_D_vec[3]^2) # heliocentric distance of asteroid at t_b
        q_D = ρ_r #signal path distance (down-leg)
        # Shapiro correction to time-delay
        Δτ_rel_D = shapiro_delay(e_D, p_D, q_D)
        # @show Δτ_rel_D
        # troposphere correction to time-delay
        Δτ_tropo_D = tropo_delay(R_r, ρ_vec_r) # seconds
        # @show Δτ_tropo_D
        # Δτ_corona_D = corona_delay(constant_term.(r_a_t_b), r_r_t_r, r_s_t_r, F_tx, station_code) # seconds
        Δτ_D = Δτ_rel_D # + Δτ_tropo_D #+ Δτ_corona_D # seconds
        # @show Δτ_D
        p_dot_23 = dot(ρ_vec_r, v_a_t_b)/ρ_r
        # @show p_dot_23
        Δt_2 = (τ_D - ρ_r/clightkms - Δτ_rel_D)/(1.0-p_dot_23/clightkms)
        # @show Δt_2
        τ_D = τ_D - Δt_2
        et_b_secs = et_r_secs - τ_D
    end
    # @show τ_D
    rv_a_t_b = xva(et_b_secs)
    r_a_t_b = rv_a_t_b[1:3]
    v_a_t_b = rv_a_t_b[4:6]

    # up-leg iteration
    # τ_U first estimation: Eq. (5) Yeomans et al. (1992)
    τ_U = τ_D
    # transmit time, 1st estimate Eq. (6) Yeomans et al. (1992)
    et_t_secs = et_b_secs - τ_U
    # @show τ_U et_t_secs et_b_secs et_r_secs
    # Geocentric position and velocity of transmitting antenna in inertial frame (au, au/day)
    R_t, V_t = observer_position(station_code, et_t_secs, eo=eo)
    # @show R_t V_t
    # Earth's barycentric position and velocity at transmit time
    rv_e_t_t = xve(et_t_secs)
    r_e_t_t = rv_e_t_t[1:3]
    v_e_t_t = rv_e_t_t[4:6]
    # @show r_e_t_t norm(r_e_t_t)/au norm(v_e_t_t)
    # Transmitter barycentric position and velocity of at transmit time
    r_t_t_t = r_e_t_t + R_t
    # # Eq. (7) Yeomans et al. (1992)
    ρ_vec_t = r_a_t_b - r_t_t_t
    ρ_t = sqrt(ρ_vec_t[1]^2 + ρ_vec_t[2]^2 + ρ_vec_t[3]^2)
    # @show ρ_t

    Δτ_U = zero(τ_U)
    Δτ_rel_U = zero(τ_U)
    # Δτ_corona_U = zero(τ_U)
    Δτ_tropo_U = zero(τ_U)
    # @show Δτ_U

    # println("   *** DOWNLEG LOOP ***")
    for i in 1:niter
        # Geocentric position and velocity of transmitting antenna in inertial frame (au, au/day)
        R_t, V_t = observer_position(station_code, et_t_secs, eo=eo)
        # Earth's barycentric position and velocity at transmit time
        rv_e_t_t = xve(et_t_secs)
        r_e_t_t = rv_e_t_t[1:3]
        v_e_t_t = rv_e_t_t[4:6]
        # Barycentric position and velocity of the transmitter at the transmit time
        r_t_t_t = r_e_t_t + R_t
        v_t_t_t = v_e_t_t + V_t
        # Eq. (7) Yeomans et al. (1992)
        ρ_vec_t = r_a_t_b - r_t_t_t
        ρ_t = sqrt(ρ_vec_t[1]^2 + ρ_vec_t[2]^2 + ρ_vec_t[3]^2)
        # @show ρ_t
        # compute up-leg Shapiro delay
        # Sun barycentric position and velocity (in au, au/day) at transmit time (TDB)
        r_s_t_t = xvs(et_t_secs)[1:3]
        e_U_vec = r_t_t_t - r_s_t_t
        # @show e_U_vec
        e_U = sqrt(e_U_vec[1]^2 + e_U_vec[2]^2 + e_U_vec[3]^2) # heliocentric distance of Earth at t_t
        # @show e_U
        r_s_t_b = xvs(et_b_secs)[1:3] # barycentric position/velocity of Sun at bounce time
        p_U_vec = r_a_t_b - r_s_t_b
        p_U = sqrt(p_U_vec[1]^2 + p_U_vec[2]^2 + p_U_vec[3]^2) # heliocentric distance of asteroid at t_b
        # @show p_U
        q_U_vec = r_a_t_b - r_e_t_t
        q_U = ρ_t # signal path distance (up-leg)
        # @show q_U
        Δτ_rel_U = shapiro_delay(e_U, p_U, q_U) # seconds
        Δτ_tropo_U = tropo_delay(R_t, ρ_vec_t) # seconds
        # Δτ_corona_U = corona_delay(constant_term.(r_t_t_t), constant_term.(r_a_t_b), constant_term.(r_s_t_b), F_tx, station_code) # seconds
        Δτ_U = Δτ_rel_U # + Δτ_tropo_U #+ Δτ_corona_U # seconds
        # @show Δτ_U
        p_dot_12 = -dot(ρ_vec_t, v_t_t_t)/ρ_t
        Δt_1 = (τ_U - ρ_t/clightkms - Δτ_rel_U)/(1.0-p_dot_12/clightkms)
        # @show Δt_1
        τ_U = τ_U - Δt_1
        # transmit time, new estimate
        et_t_secs = et_b_secs - τ_U
    end

    # compute TDB-UTC at transmit time
    # corrections to TT-TDB from Moyer (2003) / Folkner et al. (2014) due to position of measurement station on Earth are of order 0.01μs
    # Δtt_tdb_station_t = - dot(v_e_t_t, r_t_t_t-r_e_t_t)/clightkms^2
    tdb_utc_t = tdb_utc(et_t_secs) # + Δtt_tdb_station_t
    # compute TDB-UTC at receive time
    # corrections to TT-TDB from Moyer (2003) / Folkner et al. (2014) due to position of measurement station on Earth  are of order 0.01μs
    # Δtt_tdb_station_r = - dot(v_e_t_r, r_r_t_r-r_e_t_r)/clightkms^2
    tdb_utc_r = tdb_utc(et_r_secs) # + Δtt_tdb_station_r
    # @show tdb_utc_t tdb_utc_r

    # # compute total time delay (UTC seconds); relativistic delay is already included in τ_D, τ_U
    # # Eq. (9) Yeomans et al. (1992)
    τ = (τ_D + τ_U) + (Δτ_tropo_D + Δτ_tropo_U) + (tdb_utc_t - tdb_utc_r) # seconds

    return 1e6τ # total signal delay (μs)
end

function delay_doppler(station_code::Int, t_r_utc::DateTime, F_tx::Real,
        niter::Int=10; eo::Bool=true, tc::Real=1.0, xve=earth_pv, xvs=sun_pv,
        xva=apophis_pv_197, autodiff::Bool=true, tord::Int=10)
    if autodiff
        τ = delay(station_code, t_r_utc, niter, eo=eo, xve=xve, xvs=xvs, xva=xva, tord=tord)
        return τ[0], -F_tx*τ[1]
    else
        τe = delay(station_code, t_r_utc,  tc/2, niter, eo=eo, xve=xve, xvs=xvs, xva=xva)
        τn = delay(station_code, t_r_utc,   0.0, niter, eo=eo, xve=xve, xvs=xvs, xva=xva)
        τs = delay(station_code, t_r_utc, -tc/2, niter, eo=eo, xve=xve, xvs=xvs, xva=xva)
        return τn, -F_tx*((τe-τs)/tc)
    end
end

function delay_doppler(astradarfile::String,
        niter::Int=10; eo::Bool=true, tc::Real=1.0, xve=earth_pv, xvs=sun_pv,
        xva=apophis_pv_197, autodiff::Bool=true, tord::Int=10)

    astradardata = process_radar_data_jpl(astradarfile)

    et1 = str2et(string(astradardata[1].utcepoch))
    a1_et1 = xva(et1)[1]
    S = typeof(a1_et1)

    vdelay = Array{S}(undef, length(astradardata))
    vdoppler = Array{S}(undef, length(astradardata))

    for i in eachindex(astradardata)
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

    delay_index = map(x->x.delay_units=="us", astradardata)
    doppler_index = map(x->x.doppler_units=="Hz", astradardata)

    radobs_t = table(
        (
            dt_utc_obs=utcepoch.(astradardata),
            τ_obs=delay.(astradardata),
            ν_obs=doppler.(astradardata),
            τ_comp=vdelay,
            ν_comp=vdoppler,
            σ_τ=delay_sigma.(astradardata),
            σ_ν=doppler_sigma.(astradardata),
            τ_units=delay_units.(astradardata),
            ν_units=doppler_units.(astradardata),
            freq=freq.(astradardata),
            rcvr=rcvr.(astradardata),
            xmit=xmit.(astradardata),
            bouncepoint=bouncepoint.(astradardata),
            delay_index=delay_index,
            doppler_index=doppler_index
        )
    )

    return radobs_t
end

# Compute round-trip time and Doppler shift radar astrometry for an asteroid at
# UTC instant `t_r_utc` from tracking station with code `station_code` from Earth,
# Sun and asteroid ephemerides. Dopplers are computed following Yeomans et al. (1992)
# station_code: observing station identifier (MPC nomenclature)
# et_r_secs: time of echo reception (TDB seconds since J2000.0 TDB)
# F_tx: transmitter frequency (MHz)
# niter: number of light-time solution iterations
# xve: Earth ephemeris wich takes et seconds since J2000 as input and returns Earth barycentric position in au and velocity in au/day
# xvs: Sun ephemeris wich takes et seconds since J2000 as input and returns Sun barycentric position in au and velocity in au/day
# xva: asteroid ephemeris wich takes et seconds since J2000 as input and returns asteroid barycentric position in au and velocity in au/day
function delay_doppler_yeomansetal92(station_code::Int, t_r_utc::DateTime,
        F_tx::Real, niter::Int=10; eo::Bool=true, xve=earth_pv, xvs=sun_pv,
        xva=apophis_pv_197)
    # Transform receiving time from UTC to TDB seconds since j2000
    et_r_secs = str2et(string(t_r_utc))
    # Compute geocentric position/velocity of receiving antenna in inertial frame (au, au/day)
    R_r, V_r = observer_position(station_code, et_r_secs, eo=eo)
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
    # down-leg iteration
    # τ_D first approximation: Eq. (1) Yeomans et al. (1992)
    ρ_vec_r = r_a_t_r - r_r_t_r
    ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)
    τ_D = ρ_r/clightkms # (seconds) -R_b/c, but delay is wrt asteroid Center (Brozovic et al., 2018)
    # bounce time, 1st estimate, Eq. (2) Yeomans et al. (1992)
    et_b_secs = et_r_secs - τ_D
    # asteroid barycentric position (in au) at bounce time (TDB)
    rv_a_t_b = xva(et_b_secs)
    r_a_t_b = rv_a_t_b[1:3]
    v_a_t_b = rv_a_t_b[4:6]

    Δτ_D = zero(τ_D)
    Δτ_rel_D = zero(τ_D)
    # Δτ_corona_D = zero(τ_D)
    Δτ_tropo_D = zero(τ_D)

    for i in 1:niter
        # Eq. (3) Yeomans et al. (1992)
        ρ_vec_r = r_a_t_b - r_r_t_r
        ρ_vec_dot_r = v_a_t_b - v_r_t_r
        # Eq. (4) Yeomans et al. (1992)
        ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)
        τ_D = ρ_r/clightkms # (seconds) -R_b/c (COM correction) + Δτ_D (relativistic, tropo, iono...)
        # bounce time, new estimate Eq. (2) Yeomans et al. (1992)
        et_b_secs = et_r_secs - τ_D
        # asteroid barycentric position (in au) at bounce time (TDB)
        rv_a_t_b = xva(et_b_secs)
        r_a_t_b = rv_a_t_b[1:3]
        v_a_t_b = rv_a_t_b[4:6]
        # compute down-leg Shapiro delay
        # NOTE: when using PPN, substitute 2 -> 1+γ in expressions for Shapiro delay, Δτ_rel_[D|U]
        e_D_vec  = r_e_t_r - r_s_t_r
        de_D_vec = v_e_t_r - v_s_t_r
        e_D = sqrt(e_D_vec[1]^2 + e_D_vec[2]^2 + e_D_vec[3]^2) # heliocentric distance of Earth at t_r
        rv_s_t_b = xvs(et_b_secs) # barycentric position and velocity of Sun at estimated bounce time
        r_s_t_b = rv_s_t_b[1:3]
        v_s_t_b = rv_s_t_b[4:6]
        p_D_vec  = constant_term.(r_a_t_b - r_s_t_b)
        dp_D_vec = constant_term.(v_a_t_b - v_s_t_b)
        p_D = sqrt(p_D_vec[1]^2 + p_D_vec[2]^2 + p_D_vec[3]^2) # heliocentric distance of asteroid at t_b
        q_D_vec  = constant_term.(r_a_t_b - r_e_t_r)
        dq_D_vec = constant_term.(v_a_t_b - v_e_t_r)
        q_D = sqrt(q_D_vec[1]^2 + q_D_vec[2]^2 + q_D_vec[3]^2) #signal path distance (down-leg)
        # Shapiro correction to time-delay
        Δτ_rel_D = shapiro_delay(e_D, p_D, q_D)
        # troposphere correction to time-delay
        Δτ_tropo_D = tropo_delay(R_r, ρ_vec_r) # seconds
        # Δτ_corona_D = corona_delay(constant_term.(r_a_t_b), r_r_t_r, r_s_t_r, F_tx, station_code) # seconds
        Δτ_D = Δτ_rel_D + Δτ_tropo_D #+ Δτ_corona_D # seconds
    end
    τ_D = τ_D + Δτ_D
    # get latest estimates of ρ_vec_r and ρ_r
    # Eq. (3) Yeomans et al. (1992)
    ρ_vec_r = r_a_t_b - r_r_t_r
    # Eq. (4) Yeomans et al. (1992)
    ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)

    # up-leg iteration
    # τ_U first estimation: Eq. (5) Yeomans et al. (1992)
    τ_U = τ_D
    # transmit time, 1st estimate Eq. (6) Yeomans et al. (1992)
    et_t_secs = et_r_secs - (τ_U+τ_D)
    # Geocentric position and velocity of transmitting antenna in inertial frame (au, au/day)
    R_t, V_t = observer_position(station_code, constant_term(et_t_secs), eo=eo)
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
    # Δτ_corona_U = zero(τ_U)
    Δτ_tropo_U = zero(τ_U)

    for i in 1:niter
        # Eq. (8) Yeomans et al. (1992)
        τ_U = ρ_t/clightkms # (seconds) -R_b/c (COM correction) + Δτ_U (relativistic, tropo, iono...)
        # transmit time, new estimate
        et_t_secs = et_r_secs-(τ_U+τ_D)
        # Geocentric position and velocity of transmitting antenna in inertial frame (au, au/day)
        R_t, V_t = observer_position(station_code, constant_term(et_t_secs), eo=eo)
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
        # compute up-leg Shapiro delay
        # Sun barycentric position and velocity (in au, au/day) at transmit time (TDB)
        rv_s_t_t = xvs( et_t_secs )
        r_s_t_t = rv_s_t_t[1:3]
        v_s_t_t = rv_s_t_t[4:6]
        e_U_vec = constant_term.(r_e_t_t - r_s_t_t)
        de_U_vec = constant_term.(v_e_t_t - v_s_t_t)
        e_U = sqrt(e_U_vec[1]^2 + e_U_vec[2]^2 + e_U_vec[3]^2) # heliocentric distance of Earth at t_t
        rv_s_t_b = xvs(et_b_secs) # barycentric position/velocity of Sun at bounce time
        r_s_t_b = rv_s_t_b[1:3]
        v_s_t_b = rv_s_t_b[4:6]
        p_U_vec = constant_term.(r_a_t_b - r_s_t_b)
        dp_U_vec = constant_term.(v_a_t_b - v_s_t_b)
        p_U = sqrt(p_U_vec[1]^2 + p_U_vec[2]^2 + p_U_vec[3]^2) # heliocentric distance of asteroid at t_b
        q_U_vec = constant_term.(r_a_t_b - r_e_t_t)
        dq_U_vec = constant_term.(v_a_t_b - v_e_t_t)
        q_U = sqrt(q_U_vec[1]^2 + q_U_vec[2]^2 + q_U_vec[3]^2) # signal path (up-leg)
        Δτ_rel_U = shapiro_delay(e_U, p_U, q_U) # seconds
        Δτ_tropo_U = tropo_delay(R_t, ρ_vec_t) # seconds
        # Δτ_corona_U = corona_delay(constant_term.(r_t_t_t), constant_term.(r_a_t_b), constant_term.(r_s_t_b), F_tx, station_code) # seconds
        Δτ_U = Δτ_rel_U + Δτ_tropo_U #+ Δτ_corona_U # seconds
    end
    # Sun barycentric position (in au) at transmit time (TDB)
    r_s_t_t = xvs(et_t_secs)[1:3]
    τ_U = τ_U + Δτ_U

    # compute TDB-UTC at transmit time
    # corrections to TT-TDB from Moyer (2003) / Folkner et al. (2014) due to position of measurement station on Earth are of order 0.01μs
    # Δtt_tdb_station_t = - dot(v_e_t_t/daysec, r_t_t_t-r_e_t_t)/c_au_per_sec^2
    tdb_utc_t = tdb_utc(constant_term(et_t_secs)) #+ Δtt_tdb_station_t
    # compute TDB-UTC at receive time
    # corrections to TT-TDB from Moyer (2003) / Folkner et al. (2014) due to position of measurement station on Earth  are of order 0.01μs
    # Δtt_tdb_station_r = dot(v_e_t_r/daysec, r_r_t_r-r_e_t_r)/c_au_per_sec^2
    tdb_utc_r = tdb_utc(et_r_secs) #+ Δtt_tdb_station_r

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
    doppler_c = (ρ_dot_t+ρ_dot_r)/clightkms
    p_t = dot(ρ_vec_t, v_t_t_t)/ρ_t
    p_r = dot(ρ_vec_r, v_r_t_r)/ρ_r
    doppler_c2_t1 = ρ_dot_t*p_t - ρ_dot_r*p_r - ρ_dot_t*ρ_dot_r # order c^-2, 1st term
    r_ts_vec = r_t_t_t - r_s_t_t
    r_ts = sqrt(r_ts_vec[1]^2+r_ts_vec[2]^2+r_ts_vec[3]^2)
    r_rs_vec = r_r_t_r - r_s_t_r
    r_rs = sqrt(r_rs_vec[1]^2+r_rs_vec[2]^2+r_rs_vec[3]^2)
    doppler_c2_t2 = (PlanetaryEphemeris.μ[su]*((au^3)/(daysec^2)))*( (1/r_ts) - (1/r_rs) ) # order c^-2, 2nd term
    doppler_c2_t3 = (  dot(v_t_t_t, v_t_t_t) - dot(v_r_t_r, v_r_t_r)  )/2  # order c^-2, 3rd term
    doppler_c2 = (doppler_c2_t1 + doppler_c2_t2 + doppler_c2_t3)/(clightkms^2)
    ν = -F_tx*(doppler_c) # + doppler_c2)

    return 1e6τ, 1e6ν # total signal delay (μs) and Doppler shift (Hz)
end
