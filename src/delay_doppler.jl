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

# Density of ionized electrons in interplanetary medium (ESAA 2014, p. 323, Sec. 8.7.5, Eq. 8.22)
# At time `t0_tdb_jul` (days) the signal is at p1 (au),
# and is somewhere between p1 and p2 (au) after `dt_tdb` (days)
# No sure if the ESAA formula has some errors (?), currently there is no
# errata reported about this
# Errata URL: (https://aa.usno.navy.mil/publications/docs/exp_supp_errata.pdf)
# Ostro (1993) gives a reference to Anderson (1978)
#
# Ne: ionized electrons density: electrons/cm^3
# p1: signal departure point (transmitter/bounce) (au)
# p2: signal arrival point (bounce/receiver) (au)
# t0_tdb_jul: time at which the delay is being evaluated (TDB Julian days)
# dt_tdb: elapsed time of signal travel (TDB days)
function Ne(p1::Vector{S}, p2::Vector{S}, t0_tdb_jul::T, dt_tdb::T) where {T<:Real, S<:Number}
    ΔS = norm(p2-p1)
    # s, the parameter of the signal path, goes from 0 to 1 (adimensional)
    s = c_au_per_day*dt_tdb/ΔS
    # @show ΔS
    # @show c_au_per_day*dt_tdb
    r_vec = p1+s*(p2-p1) - sun_pv(t0_tdb_jul+dt_tdb)[1:3] # heliocentric position (au) of point on ray path at time t_tdb_jul (Julian days)
    r = sqrt( r_vec[1]^2 + r_vec[2]^2 + r_vec[3]^2 ) # heliocentric distance (au) of point on ray path at time t_tdb_jul (Julian days)
    β = atan( r_vec[3], sqrt(r_vec[1]^2 + r_vec[2]^2) ) # solar latitude of point on ray path (rad)
    r_sr = r/R_sun
    return (A_sun/r_sr^6) + ( (a_sun*b_sun)/sqrt((a_sun*sin(β))^2 + (b_sun*cos(β))^2) )/(r_sr^2)
end

# Time-delay generated by thin plasma of solar corona
# Ne: ionized electrons density: electrons/cm^3
# p1: signal departure point (transmitter/bounce for up/down-link, resp.) (au)
# p2: signal arrival point (bounce/receiver for up/down-link, resp.) (au)
# t_1_tdb: initital time (TDB) of signal path (bounce time for downlink; transmit time for uplink)
# current_delay: total time-delay of of current signal path (TDB days)
# f_T: transmitter frequency (Hz)
# ***WARNING*** path integral is evaluated by the most naive approach
function delay_corona(p1::Vector{S}, p2::Vector{S}, t_1_tdb::T, current_delay::T, f_T::T) where {T<:Real, S<:Number}
    # vector of s parameter values along path
    sv = collect(0:0.01:1)
    nev = map(x->Ne(p1, p2, t_1_tdb, x*current_delay), sv)
    s_D_cm = c_au_per_day*current_delay*au*100000 # cm
    #ESAA
    f_T_MHz = (1e-6f_T)
    Δτ_corona_D = 40.3mean(nev)*s_D_cm/(c_cm_per_sec*f_T_MHz) # WHAT UNITS??? μsec????
    return 3e-6Δτ_corona_D/86400 # Xseconds -> days, TODO: CHECK UNITS!!!
end

# Time delay from the Earth's troposphere for radio frequencies
# oscillates between 0.007μs and 0.225μs for z between 0 and π/2 rad
delay_tropo(z) = (7e-9)/( cos(z) + 0.0014/(0.045+cot(z)) ) # seconds

# Alternate version of delay_doppler, using JPL ephemerides
# eph: Solar System + Apophis ephemerides pos/vel (au, au/day)
# station_code: observing station identifier (MPC nomenclature)
# t_r_utc_julian: time of echo reception (UTC)
# f_T: transmitter frequency (Hz)
function delay_doppler_jpleph(station_code, t_r_utc_julian, f_T, niter::Int=10)
    # UTC -> TDB
    t_r = julian(TDB, UTCEpoch(t_r_utc_julian, origin=:julian)).Δt

    # Earth-fixed geocentric position/velocity of receiving antenna (au, au/day)
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
        q_D = norm(r_a_t_b-r_e_t_r) #signal path (down-leg)
        Δτ_rel_D = 2μ[1]*log( (e_D+p_D+q_D)/(e_D+p_D-q_D) )/(c_au_per_day^3) # days
        # @show τ_D
        # bounce time, new estimate Eq. (2) Yeomans et al. (1992)
        @show τ_D
        @show Δτ_rel_D
        Δτ_corona_D = delay_corona(r_a_t_b, r_e_t_r, t_r, τ_D, f_T)
        @show Δτ_corona_D
        τ_D = τ_D + Δτ_rel_D + Δτ_corona_D
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
    # Earth-fixed geocentric position/velocity of receiving antenna (au, au/day)
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
        q_U = norm(r_a_t_b-r_e_t_t) #signal path (up-leg)
        Δτ_rel_U = 2μ[1]*log( (e_U+p_U+q_U)/(e_U+p_U-q_U) )/(c_au_per_day^3) # days
        @show τ_U
        @show Δτ_rel_U
        τ_U = τ_U + Δτ_rel_U
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
