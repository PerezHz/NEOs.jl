# x: Solar System + Apophis pos/vel at receiving time (au, au/day)
# station_code: observing station identifier (MPC nomenclature)
# t: time of echo reception (TDB)
# f_T: transmitter frequency (Hz)
function delay_doppler(x, station_code, t_r_utc_julian, f_T, niter::Int=10)
    # UTC -> TT
    t_r_tt = TTEpoch( UTCEpoch(t_r_utc_julian, origin=:julian) )
    # TT -> TDB
    t_r_tt_jul = julian(t_r_tt).Δt
    t = t_r_tt_jul - tt_m_tdb(t_r_tt_jul)[1]/86400 # TDB ≈ TT - f(TT), where f(t) = (TT-TDB)(t)

    # Compute Taylor expansion at receiving time (TDB)
    t_r = Taylor1([t, one(t)], Apophis.order)
    x_r = Taylor1.(x, Apophis.order)
    dx_r = similar(x_r)
    @time TaylorIntegration.jetcoeffs!(Val(Apophis.RNp1BP_pN_A_J234E_J2S_ng!), t_r, x_r, dx_r)

    # Earth-fixed geocentric position/velocity of receiving antenna (au, au/day)
    R_r, V_r = observer_position(station_code, t_r_utc_julian)

    # Earth's barycentric position and velocity at the receive time
    #rv_e_t_r = earth_pv(t_r)
    r_e_t_r = x_r[(3ea-2):3ea]() #rv_e_t_r[1:3]
    v_e_t_r = x_r[(3(N+ea)-2):3(N+ea)]() #rv_e_t_r[4:6]

    # Barycentric position and velocity of the receiver at the receive time
    r_r_t_r = r_e_t_r + R_r
    v_r_t_r = v_e_t_r + V_r

    # asteroid barycentric position (in a.u.) at receiving time (TDB)
    #rv_a_t_r = apophis_pv(t_r)
    r_a_t_r = x_r[apophisdofs[1:3]]() #rv_a_t_r[1:3]
    v_a_t_r = x_r[apophisdofs[4:6]]() #rv_a_t_r[4:6]

    # down-leg iteration
    # τ_D first approximation: Eq. (1) Yeomans et al. (1992)
    ρ_vec_r = r_a_t_r - r_r_t_r
    ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)
    τ_D = ρ_r/c_au_per_day # (days) -R_b/c, but delay is wrt asteroid Center (Brozovic et al., 2018)
    # @show τ_D
    # bounce time, 1st estimate, Eq. (2) Yeomans et al. (1992)
    t_b = t - τ_D
    # asteroid barycentric position (in a.u.) at bounce time (TDB)
    #rv_a_t_b = apophis_pv(t_b)
    r_a_t_b = x_r[apophisdofs[1:3]](-τ_D) #rv_a_t_b[1:3]
    v_a_t_b = x_r[apophisdofs[4:6]](-τ_D) #rv_a_t_b[4:6]
    for i in 1:niter
        # Eq. (3) Yeomans et al. (1992)
        ρ_vec_r = r_a_t_b - r_r_t_r
        # Eq. (4) Yeomans et al. (1992)
        ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)
        τ_D = ρ_r/c_au_per_day # (days) -R_b/c (COM correction) + Δτ_D (relativistic, tropo, iono...)
        # @show τ_D
        # bounce time, new estimate Eq. (2) Yeomans et al. (1992)
        t_b = t - τ_D
        # asteroid barycentric position (in a.u.) at bounce time (TDB)
        #rv_a_t_b = apophis_pv(t_b)
        r_a_t_b = x_r[apophisdofs[1:3]](-τ_D) #rv_a_t_b[1:3]
        v_a_t_b = x_r[apophisdofs[4:6]](-τ_D) #rv_a_t_b[4:6]
    end
    # @show τ_D

    # up-leg iteration
    # τ_U first estimation: Eq. (5) Yeomans et al. (1992)
    τ_U = τ_D
    # @show τ_U
    # transmit time, 1st estimate Eq. (6) Yeomans et al. (1992)
    t_t = t_b - τ_U
    # convert transmit time to UTC time
    t_t_utc_julian = julian(AstroTime.UTC, TDBEpoch(constant_term(t_t), origin=:julian)).Δt
    # Earth-fixed geocentric position/velocity of receiving antenna (au, au/day)
    R_t, V_t = observer_position(station_code, t_t_utc_julian)
    # Earth's barycentric position and velocity at the transmit time
    #rv_e_t_t = earth_pv(t_t)
    r_e_t_t = x_r[(3ea-2):3ea](t_t-t) #rv_e_t_t[1:3]
    v_e_t_t = x_r[(3(N+ea)-2):3(N+ea)](t_t-t) #rv_e_t_t[4:6]
    # Barycentric position and velocity of the transmitter at the transmit time
    r_t_t_t = r_e_t_t + R_t
    v_t_t_t = v_e_t_t + V_t
    # Eq. (7) Yeomans et al. (1992)
    ρ_vec_t = r_a_t_b - r_t_t_t
    ρ_t = sqrt(ρ_vec_t[1]^2 + ρ_vec_t[2]^2 + ρ_vec_t[3]^2)
    for i in 1:niter
        # Eq. (8) Yeomans et al. (1992)
        τ_U = ρ_t/c_au_per_day # (days) -R_b/c (COM correction) + Δτ_U (relativistic, tropo, iono...)
        # @show τ_U
        # transmit time, 1st estimate Eq. (6) Yeomans et al. (1992)
        t_t = t_b - τ_U
        # convert transmit time to UTC time
        t_t_utc_julian = julian(AstroTime.UTC, TDBEpoch(constant_term(t_t), origin=:julian)).Δt
        # Earth-fixed geocentric position/velocity of receiving antenna (au, au/day)
        R_t, V_t = observer_position(station_code, t_t_utc_julian)
        # Earth's barycentric position and velocity at the transmit time
        #rv_e_t_t = earth_pv(t_t)
        r_e_t_t = x_r[(3ea-2):3ea](t_t-t) #rv_e_t_t[1:3]
        v_e_t_t = x_r[(3(N+ea)-2):3(N+ea)](t_t-t) #rv_e_t_t[4:6]
        # Barycentric position and velocity of the transmitter at the transmit time
        r_t_t_t = r_e_t_t + R_t
        v_t_t_t = v_e_t_t + V_t
        # Eq. (7) Yeomans et al. (1992)
        ρ_vec_t = r_a_t_b - r_t_t_t
        ρ_t = sqrt(ρ_vec_t[1]^2 + ρ_vec_t[2]^2 + ρ_vec_t[3]^2)
    end
    # @show τ_U

    # Eq. (9) Yeomans et al. (1992)
    #t_r_UTC = UTCEpoch(TDBEpoch(constant_term(t_r), origin=:julian))
    #t_t_UTC = UTCEpoch(TDBEpoch(constant_term(t_t), origin=:julian))
    #total_time_delay = ( t_r_UTC - t_t_UTC ).Δt
    #τ = τ_D + τ_U
    t_r_TDB = TDBEpoch(t, origin=:julian)
    t_r_UTC = UTCEpoch(t_r_TDB)
    TDB_minus_UTC_r = t - julian(t_r_UTC).Δt
    t_t_TDB = TDBEpoch(constant_term(t_t), origin=:julian)
    t_t_UTC = UTCEpoch(t_t_TDB)
    TDB_minus_UTC_t = constant_term(t_t) - julian(t_t_UTC).Δt
    total_time_delay = (τ_U + τ_D) + (TDB_minus_UTC_t - TDB_minus_UTC_r)
    #@show TDB_minus_UTC_t - TDB_minus_UTC_r
    #@show total_time_delay

    # Eq. (10) Yeomans et al. (1992)
    ρ_vec_dot_t = v_a_t_b - v_t_t_t
    ρ_vec_dot_r = v_a_t_b - v_r_t_r

    # Eq. (11) Yeomans et al. (1992)
    ρ_dot_t = dot(ρ_vec_t, ρ_vec_dot_t)/ρ_t
    ρ_dot_r = dot(ρ_vec_r, ρ_vec_dot_r)/ρ_r

    # @show ρ_dot_t
    # @show ρ_dot_r

    # Eq. (12) Yeomans et al. (1992)
    doppler1 = -f_T*(ρ_dot_t+ρ_dot_r)/c_au_per_day
    doppler2 = (ρ_dot_t/ρ_t)*dot(ρ_vec_t, v_t_t_t) - (ρ_dot_r/ρ_r)*dot(ρ_vec_r, v_r_t_r) - ρ_dot_t*ρ_dot_r
    r_s_t_r = x_r[sundofs[1:3]]() #sun_pv(t_r)[1:3]
    r_s_t_t = x_r[sundofs[1:3]](t_t-t) #sun_pv(t_t)[1:3]
    r_ts = r_t_t_t - r_s_t_t
    r_rs = r_r_t_r - r_s_t_r
    factor = 1/sqrt(r_ts[1]^2+r_ts[2]^2+r_ts[3]^2) - 1/sqrt(r_rs[1]^2+r_rs[2]^2+r_rs[3]^2)
    doppler3 = μ[1]*factor
    doppler4 = (  dot(v_t_t_t, v_t_t_t) - dot(v_r_t_r, v_r_t_r)  )/2
    doppler_234 = -f_T*(doppler2 + doppler3 + doppler4)/(c_au_per_day^2)
    f_D = doppler1+doppler_234

    # compute Shapiro delay
    # down-leg
    r_es_t_r = r_e_t_r-r_s_t_r
    e_D = sqrt( r_es_t_r[1]^2 + r_es_t_r[2]^2 + r_es_t_r[3]^2 ) # heliocentric distance of Earth at t_r
    r_s_t_b = x_r[sundofs[1:3]](t_b-t)
    r_as_t_b = r_a_t_b-r_s_t_b
    p_D = sqrt( r_as_t_b[1]^2 + r_as_t_b[2]^2 + r_as_t_b[3]^2 ) # heliocentric distance of asteroid at t_b
    signal_path_d = r_a_t_b-r_e_t_r
    q_D = sqrt( signal_path_d[1]^2 + signal_path_d[2]^2 + signal_path_d[3]^2 ) #signal path (down-leg)
    Δτ_rel_D = (1+1)*μ[1]*log( (e_D+p_D+q_D)/(e_D+p_D-q_D) )/(c_au_per_day^3)
    # up-leg
    r_es_t_t = r_e_t_r-r_s_t_t
    e_U = sqrt( r_es_t_t[1]^2 + r_es_t_t[2]^2 + r_es_t_t[3]^2 ) # heliocentric distance of Earth at t_t
    p_U = p_D # heliocentric distance of asteroid at t_b
    signal_path_u = r_a_t_b-r_e_t_t
    q_U = sqrt( signal_path_u[1]^2 + signal_path_u[2]^2 + signal_path_u[3]^2 ) #signal path (up-leg)
    Δτ_rel_U = (1+1)*μ[1]*log( (e_U+p_U+q_U)/(e_U+p_U-q_U) )/(c_au_per_day^3)
    # total
    Δτ_rel = Δτ_rel_D + Δτ_rel_U # seconds
    # @show Δτ_rel

    return (1e6*86400)*(total_time_delay+Δτ_rel), f_D # total signal delay (μs) and Doppler shift (Hz)
end

const eph = Ephem(["jpleph/a99942.bsp", "jpleph/de430_1850-2150.bsp", "jpleph/TTmTDB.de430.19feb2015.bsp"])
prefetch(eph)

# CALCEPH index conventions:
#12: Solar System Barycenter
#11: Sun (heliocenter)
#2099942: Apophis
#3: Earth (geocenter)
#10: Moon
#1000000001 from body 1000000000: TT-TDB
apophis_pv(t) = compute(eph, t, 0.0, 2099942   , 12        , unitKM+unitDay, 1)/au # units: au, au/day
sun_pv(t)     = compute(eph, t, 0.0, 11        , 12        , unitKM+unitDay, 1)/au # units: au, au/day
earth_pv(t)   = compute(eph, t, 0.0, 3         , 12        , unitKM+unitDay, 1)/au # units: au, au/day
moon_pv(t)    = compute(eph, t, 0.0, 10        , 12        , unitKM+unitDay, 1)/au # units: au, au/day
tt_m_tdb(t)   = compute(eph, t, 0.0, 1000000001, 1000000000, unitSec       , 1) # units: seconds

# Shapiro delay
function shapiro_delay(e, p, q)
    return 2μ[1]*log( (e+p+q)/(e+p-q) )/(c_au_per_day^3) # days
end

# Density of ionized electrons in interplanetary medium (ESAA 2014, p. 323, Sec. 8.7.5, Eq. 8.22)
# At time `t0_tdb_jul` (days) the signal is at p1 (au),
# and propagates from p1 to p2 at the speed of light
# ESAA 2014 in turn refers to Muhleman and Anderson (1981)
# Not sure if the ESAA formula has some errors (?), currently there is no
# errata reported about this
# But it seems as if Eq. (8.23) is for the Earth ionosphere delay, instead of for the solar corona
# Errata URL: (https://aa.usno.navy.mil/publications/docs/exp_supp_errata.pdf)
# Ostro (1993) gives a reference to Anderson (1978), where this model is fitted to Mariner 9 ranging data
# Reading https://gssc.esa.int/navipedia/index.php/Ionospheric_Delay
# Helped a lot to clarify things, especially the 40.3, although they talk about Earth's ionosphere
# Ne: ionized electrons density: electrons/cm^3
# p1: signal departure point (transmitter/bounce) (au)
# p2: signal arrival point (bounce/receiver) (au)
# t0_tdb_jul: time at which the delay is being evaluated (TDB Julian days)
# ds: distance travelled by ray (au)
function Ne(p1::Vector{S}, p2::Vector{S}, t0_tdb_jul::T, ds::U) where {T<:Real, S<:Number, U<:Number}
    ΔS = norm(p2-p1)/R_sun
    # s: linear parametrization of ray path, such that
    # s=0 -> point on ray path is at p1; s=1 -> point on ray path is at p2
    s = ds/ΔS
    s_p2_p1 = map(x->s*x, (p2-p1))
    dt = constant_term(ds/c_au_per_day) # (TDB days)
    # Linear interpolation of sun position: r_sun(t0+dt) ≈ r_sun(t0)+dt*v_sun(t0)
    vr_s_t0 = sun_pv(t0_tdb_jul)
    r_s_t0 = vr_s_t0[1:3]
    v_s_t0 = vr_s_t0[4:6]
    dt_v_s_t0 = map(x->dt*x, v_s_t0)
    r_s_t0_dt = r_s_t0 + dt_v_s_t0
    r_vec = p1+s_p2_p1 - r_s_t0_dt # heliocentric position (au) of point on ray path at time t_tdb_jul (Julian days)
    r = sqrt( r_vec[1]^2 + r_vec[2]^2 + r_vec[3]^2 ) # heliocentric distance (au) of point on ray path at time t_tdb_jul (Julian days)
    # compute ecliptic solar latitude of point on ray path
    # NOTE: actually, (Anderson, 1978) says it should be heliographic
    # (i.e., wrt solar equator), but diff is ~7 deg and things don't seem to change a lot
    # compute heliographic latitude (ecliptic plane is not equal to heliographic but almost by ~7deg)
    r_vec_ecliptic = Rx(deg2rad(23.43929))*r_vec
    β = asin( r_vec_ecliptic[3]/r ) # ecliptic solar latitude of point on ray path (rad)
    r_sr = r/R_sun
    return (A_sun/r_sr^6) + ( (a_sun*b_sun)/sqrt((a_sun*sin(β))^2 + (b_sun*cos(β))^2) )/(r_sr^2)
end

# Integral of Ne, evaluated with TaylorIntegration.jl
# TODO: @taylorize!
# p1: signal departure point (transmitter/bounce) (au)
# p2: signal arrival point (bounce/receiver) (au)
# t0_tdb_jul: time at which the delay is being evaluated (TDB Julian days)
# current_delay: total time-delay for current signal path (TDB days)
# output has units of electrons/cm^2
function Ne_path_integral(p1::Vector{S}, p2::Vector{S}, t0_tdb_jul::T, current_delay::T) where {T<:Real, S<:Number}
    # kernel of path integral; distance is in au
    function int_kernel(s, x)
        return Ne(p1, p2, t0_tdb_jul, s)
    end
    # do path integral; distance vector sv is in au; integral vector iv is in (electrons/cm^3)*au
    sv, iv = taylorinteg(int_kernel, 0.0, 0.0, c_au_per_day*current_delay/R_sun, 28, 1e-20)
    return iv[end] # (electrons/cm^3)*solar-radii
end

# Time-delay generated by thin plasma of solar corona (but ESAA's formula seems to be for ionosphere delay)
# Ne: ionized electrons density: electrons/cm^3
# p1: signal departure point (transmitter/bounce for up/down-link, resp.) (au)
# p2: signal arrival point (bounce/receiver for up/down-link, resp.) (au)
# t_1_tdb: initial time (TDB) of signal path (bounce time for downlink; transmit time for uplink)
# current_delay: total time-delay of of current signal path (TDB days)
# f_T: transmitter frequency (Hz)
# From https://gssc.esa.int/navipedia/index.php/Ionospheric_Delay it seems that
# the expression 40.3*Ne/f^2 is adimensional, where Ne [electrons/meters^3] and f is in Hz
# therefore, the integral (40.3/f^2)*∫Ne*ds is in meters, where ds is in meters
# and the expression (40.3/(c*f^2))*∫Ne*ds, with c in meters/second, is in seconds
function corona_delay(p1::Vector{S}, p2::Vector{S}, t_1_tdb::T, current_delay::T, f_T::T) where {T<:Real, S<:Number}
    f_T_MHz = (1e-6f_T)
    int_path = Ne_path_integral(p1, p2, t_1_tdb, current_delay) # (electrons/cm^3)*solar-radii
    Δτ_corona_D_μs = 187.0int_path/(f_T_MHz^2) # μseconds
    return 1e-6Δτ_corona_D_μs/86400 # μseconds -> days
end

# *VERY* elementary computation of zenith distance
# r_antenna: position of antenna at receive/transmit time in celestial frame wrt geocenter
# ρ_vec_ae: slant-range vector from antenna to asteroid
function zenith_distance(r_antenna::Vector{S}, ρ_vec_ae::Vector{S}) where {S<:Number}
    cos_antenna_slant = dot(r_antenna, ρ_vec_ae)/(norm(r_antenna)*norm(ρ_vec_ae))
    return acos(cos_antenna_slant)
end

# Time delay from the Earth's troposphere for radio frequencies
# oscillates between 0.007μs and 0.225μs for z between 0 and π/2 rad
tropo_delay(z) = (7e-9/86400)/( cos(z) + 0.0014/(0.045+cot(z)) ) # days

function tropo_delay(r_antenna::Vector{S}, ρ_vec_ae::Vector{S}) where {S<:Number}
    zd = zenith_distance(r_antenna, ρ_vec_ae)
    # @show rad2deg(zd)
    return tropo_delay(zd)
end

# Alternate version of delay_doppler, using JPL ephemerides
# eph: Solar System + Apophis ephemerides pos/vel (au, au/day)
# station_code: observing station identifier (MPC nomenclature)
# t_r_utc_julian: time of echo reception (UTC)
# f_T: transmitter frequency (Hz)
function delay_doppler_jpleph(station_code, t_r_utc_julian, f_T, niter::Int=10)
    # UTC -> TDB
    t_r = julian(TDB, UTCEpoch(t_r_utc_julian, origin=:julian)).Δt

    # Geocentric position/velocity of receiving antenna in inertial frame (au, au/day)
    R_r, V_r = observer_position(station_code, t_r_utc_julian)

    # Earth's barycentric position and velocity at the receive time
    rv_e_t_r = earth_pv(t_r)
    r_e_t_r = rv_e_t_r[1:3]
    v_e_t_r = rv_e_t_r[4:6]

    # Barycentric position and velocity of the receiver at the receive time
    r_r_t_r = r_e_t_r + R_r
    v_r_t_r = v_e_t_r + V_r

    # asteroid barycentric position (in a.u.) at receiving time (TDB)
    rv_a_t_r = apophis_pv(t_r)
    r_a_t_r = rv_a_t_r[1:3]
    v_a_t_r = rv_a_t_r[4:6]

    # Sun barycentric position (in a.u.) at receiving time (TDB)
    r_s_t_r = sun_pv(t_r)[1:3]

    # down-leg iteration
    # τ_D first approximation: Eq. (1) Yeomans et al. (1992)
    ρ_vec_r = r_a_t_r - r_r_t_r
    ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)
    τ_D = ρ_r/c_au_per_day # (days) -R_b/c, but delay is wrt asteroid Center (Brozovic et al., 2018)
    # bounce time, 1st estimate, Eq. (2) Yeomans et al. (1992)
    t_b = t_r - τ_D
    # asteroid barycentric position (in a.u.) at bounce time (TDB)
    rv_a_t_b = apophis_pv(t_b)
    r_a_t_b = rv_a_t_b[1:3]
    v_a_t_b = rv_a_t_b[4:6]
    for i in 1:niter
        # Eq. (3) Yeomans et al. (1992)
        ρ_vec_r = r_a_t_b - r_r_t_r
        # Eq. (4) Yeomans et al. (1992)
        ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)
        τ_D = ρ_r/c_au_per_day # (days) -R_b/c (COM correction) + Δτ_D (relativistic, tropo, iono...)
        # compute down-leg Shapiro delay
        # NOTE: when using PPN, substitute 2 -> 1+γ in expressions for Shapiro delay, Δτ_rel_[D|U]
        e_D = norm(r_e_t_r-r_s_t_r) # heliocentric distance of Earth at t_r
        r_s_t_b = sun_pv(t_r - τ_D)[1:3] # barycentric position of Sun at estimated bounce time
        p_D = norm(r_a_t_b-r_s_t_b) # heliocentric distance of asteroid at t_b
        q_D = norm(ρ_vec_r) #signal path (down-leg)
        Δτ_rel_D = shapiro_delay(e_D, p_D, q_D) # days
        Δτ_corona_D = corona_delay(r_a_t_b, r_e_t_r, t_b, τ_D, f_T)
        Δτ_tropo_D = tropo_delay(R_r, ρ_vec_r)
        τ_D = τ_D + Δτ_rel_D # + Δτ_corona_D + Δτ_tropo_D
        # @show τ_D
        @show Δτ_rel_D
        @show Δτ_corona_D
        @show Δτ_tropo_D
        # bounce time, new estimate Eq. (2) Yeomans et al. (1992)
        t_b = t_r - τ_D
        # asteroid barycentric position (in a.u.) at bounce time (TDB)
        rv_a_t_b = apophis_pv(t_b)
        r_a_t_b = rv_a_t_b[1:3]
        v_a_t_b = rv_a_t_b[4:6]
    end
    # get latest estimates of ρ_vec_r and ρ_r
    # Eq. (3) Yeomans et al. (1992)
    ρ_vec_r = r_a_t_b - r_r_t_r
    # Eq. (4) Yeomans et al. (1992)
    ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)

    # up-leg iteration
    # τ_U first estimation: Eq. (5) Yeomans et al. (1992)
    τ_U = τ_D
    # transmit time, 1st estimate Eq. (6) Yeomans et al. (1992)
    t_t = t_b - τ_U
    # convert transmit time to UTC time
    # TDB -> UTC
    t_t_utc_julian = julian(AstroTime.UTC, TDBEpoch(t_t, origin=:julian)).Δt
    # Geocentric position/velocity of receiving antenna in inertial frame (au, au/day)
    R_t, V_t = observer_position(station_code, t_t_utc_julian)
    rv_e_t_t = earth_pv(t_t)
    r_e_t_t = rv_e_t_t[1:3]
    v_e_t_t = rv_e_t_t[4:6]
    # Barycentric position and velocity of the transmitter at the transmit time
    r_t_t_t = r_e_t_t + R_t
    v_t_t_t = v_e_t_t + V_t
    # Eq. (7) Yeomans et al. (1992)
    ρ_vec_t = r_a_t_b - r_t_t_t
    ρ_t = sqrt(ρ_vec_t[1]^2 + ρ_vec_t[2]^2 + ρ_vec_t[3]^2)
    # Sun barycentric position (in a.u.) at transmit time (TDB)
    r_s_t_t = sun_pv(t_b - τ_U)[1:3]
    for i in 1:niter
        # Eq. (8) Yeomans et al. (1992)
        τ_U = ρ_t/c_au_per_day # (days) -R_b/c (COM correction) + Δτ_U (relativistic, tropo, iono...)
        # Sun barycentric position (in a.u.) at transmit time (TDB)
        r_s_t_t = sun_pv(t_b - τ_U)[1:3]
        # compute up-leg Shapiro delay
        e_U = norm(r_e_t_t-r_s_t_t) # heliocentric distance of Earth at t_t
        r_s_t_b = sun_pv(t_b)[1:3] # barycentric position of Sun at bounce time
        p_U = norm(r_a_t_b-r_s_t_b) # heliocentric distance of asteroid at t_b
        q_U = norm(ρ_vec_t) #signal path (up-leg)
        Δτ_rel_U = shapiro_delay(e_U, p_U, q_U) # days
        Δτ_corona_U = corona_delay(r_e_t_t, r_a_t_b, t_t, τ_U, f_T) # days
        Δτ_tropo_U = tropo_delay(R_t, ρ_vec_t) # days
        τ_U = τ_U + Δτ_rel_U # + Δτ_corona_U + Δτ_tropo_U
        # @show τ_U
        # @show Δτ_rel_U
        # @show Δτ_corona_U
        # @show Δτ_tropo_U
        # transmit time, 1st estimate Eq. (6) Yeomans et al. (1992)
        t_t = t_b - τ_U
        # convert transmit time to UTC time
        # TDB -> UTC
        t_t_utc_julian = julian(AstroTime.UTC, TDBEpoch(t_t, origin=:julian)).Δt
        # Earth-fixed geocentric position/velocity of receiving antenna (au, au/day)
        R_t, V_t = observer_position(station_code, t_t_utc_julian)
        # Earth's barycentric position and velocity at the transmit time
        rv_e_t_t = earth_pv(t_t)
        r_e_t_t = rv_e_t_t[1:3]
        v_e_t_t = rv_e_t_t[4:6]
        # Barycentric position and velocity of the transmitter at the transmit time
        r_t_t_t = r_e_t_t + R_t
        v_t_t_t = v_e_t_t + V_t
        # Eq. (7) Yeomans et al. (1992)
        ρ_vec_t = r_a_t_b - r_t_t_t
        ρ_t = sqrt(ρ_vec_t[1]^2 + ρ_vec_t[2]^2 + ρ_vec_t[3]^2)
    end

    # compute total time delay (UTC seconds)
    # Eq. (9) Yeomans et al. (1992)
    #τ = τ_D + τ_U + (TDB-UTC)_t - (TDB-UTC)_r
    # TDB -> UTC
    t_r_UTC = UTCEpoch( TDBEpoch(t_r, origin=:julian) )
    t_t_UTC = UTCEpoch( TDBEpoch(t_t, origin=:julian) )
    total_time_delay = ( t_r_UTC - t_t_UTC ).Δt

    # compute Doppler shift ν
    # Eq. (10) Yeomans et al. (1992)
    ρ_vec_dot_t = v_a_t_b - v_t_t_t
    ρ_vec_dot_r = v_a_t_b - v_r_t_r
    # Eq. (11) Yeomans et al. (1992)
    ρ_dot_t = dot(ρ_vec_t, ρ_vec_dot_t)/ρ_t
    ρ_dot_r = dot(ρ_vec_r, ρ_vec_dot_r)/ρ_r
    # Eq. (12) Yeomans et al. (1992)
    doppler1 = -f_T*(ρ_dot_t+ρ_dot_r)/c_au_per_day
    doppler2 = (ρ_dot_t/ρ_t)*dot(ρ_vec_t, v_t_t_t) - (ρ_dot_r/ρ_r)*dot(ρ_vec_r, v_r_t_r) - ρ_dot_t*ρ_dot_r
    r_ts = r_t_t_t - r_s_t_t
    r_rs = r_r_t_r - r_s_t_r
    factor = 1/sqrt(r_ts[1]^2+r_ts[2]^2+r_ts[3]^2) - 1/sqrt(r_rs[1]^2+r_rs[2]^2+r_rs[3]^2)
    doppler3 = μ[1]*factor
    doppler4 = (  dot(v_t_t_t, v_t_t_t) - dot(v_r_t_r, v_r_t_r)  )/2
    doppler_234 = -f_T*(doppler2 + doppler3 + doppler4)/(c_au_per_day^2)
    ν = doppler1+doppler_234

    return 1e6total_time_delay, ν # total signal delay (μs) and Doppler shift (Hz)
end
