# x: Solar System + Apophis pos/vel at receiving time (au, au/day)
# station_code: observing station identifier (MPC nomenclature)
# t: time of echo reception (TDB)
# f_T: transmitter frequency (Hz)
function delay_doppler(x, station_code::Int, t_r_utc::DateTime, f_T, niter::Int=10)
    # UTC -> TDB (receiving time)
    t_r = julian(TDB, UTCEpoch(t_r_utc)).Δt
    @show t_r

    # Compute Taylor expansion at receiving time (TDB)
    t_rT1 = Taylor1([t_r, one(t_r)], Apophis.order)
    x_r = Taylor1.(x, Apophis.order)
    dx_r = similar(x_r)
    @time TaylorIntegration.jetcoeffs!(Val(Apophis.RNp1BP_pN_A_J234E_J2S_ng!), t_rT1, x_r, dx_r)

    # Geocentric position/velocity of receiving antenna in inertial frame (au, au/day)
    R_r, V_r = observer_position(station_code, t_r_utc)

    # Earth's barycentric position and velocity at the receive time
    #rv_e_t_r = earth_pv(t_r)
    r_e_t_r = x_r[earthdofs[1:3]]() # rv_e_t_r[1:3]
    v_e_t_r = x_r[earthdofs[4:6]]() # rv_e_t_r[4:6]

    # Barycentric position and velocity of the receiver at the receive time
    r_r_t_r = r_e_t_r + R_r
    v_r_t_r = v_e_t_r + V_r

    # asteroid barycentric position (in au) at receiving time (TDB)
   # rv_a_t_r = apophis_pv(t_r)
    r_a_t_r = x_r[apophisdofs[1:3]]() # rv_a_t_r[1:3]
    v_a_t_r = x_r[apophisdofs[4:6]]() # rv_a_t_r[4:6]

    # Sun barycentric position (in au) at receiving time (TDB)
    r_s_t_r = x_r[sundofs[1:3]]() # sun_pv(t_r)[1:3]

    # down-leg iteration
    # τ_D first approximation: Eq. (1) Yeomans et al. (1992)
    ρ_vec_r = r_a_t_r - r_r_t_r
    ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)
    τ_D = ρ_r/c_au_per_day # (days) -R_b/c, but delay is wrt asteroid Center (Brozovic et al., 2018)
    # bounce time, 1st estimate, Eq. (2) Yeomans et al. (1992)
    t_b = t_r - τ_D
    # asteroid barycentric position (in au) at bounce time (TDB)
    #rv_a_t_b = apophis_pv(t_r, - τ_D)
    r_a_t_b = x_r[apophisdofs[1:3]](-τ_D) #rv_a_t_b[1:3]
    v_a_t_b = x_r[apophisdofs[4:6]](-τ_D) # rv_a_t_b[4:6]
    for i in 1:niter
        # Eq. (3) Yeomans et al. (1992)
        ρ_vec_r = r_a_t_b - r_r_t_r
        # Eq. (4) Yeomans et al. (1992)
        ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)
        τ_D = ρ_r/c_au_per_day # (days) -R_b/c (COM correction) + Δτ_D (relativistic, tropo, iono...)
        # compute down-leg Shapiro delay
        # NOTE: when using PPN, substitute 2 -> 1+γ in expressions for Shapiro delay, Δτ_rel_[D|U]
        e_D = norm(r_e_t_r-r_s_t_r) # heliocentric distance of Earth at t_r
        r_s_t_b = x_r[sundofs[1:3]](-τ_D) # sun_pv(t_r, -τ_D)[1:3] # barycentric position of Sun at estimated bounce time
        v_s_t_b = x_r[sundofs[4:6]](-τ_D) # sun_pv(t_r, -τ_D)[1:3] # barycentric velocity of Sun at estimated bounce time
        p_D = norm(r_a_t_b-r_s_t_b) # heliocentric distance of asteroid at t_b
        q_D = norm(ρ_vec_r) #signal path (down-leg)
        Δτ_rel_D = shapiro_delay(e_D, p_D, q_D) # days
        Δτ_corona_D = corona_delay(r_a_t_b, r_e_t_r, r_s_t_b, v_s_t_b, f_T, station_code)
        Δτ_tropo_D = tropo_delay(R_r, ρ_vec_r)
        Δτ_D = Δτ_corona_D + (Δτ_rel_D + Δτ_tropo_D)
        τ_D = τ_D + Δτ_D
        @show τ_D
        @show Δτ_rel_D*(1e6*daysec)
        @show Δτ_corona_D*(1e6*daysec)
        @show Δτ_tropo_D*(1e6*daysec)
        # bounce time, new estimate Eq. (2) Yeomans et al. (1992)
        t_b = t_r - τ_D
        # asteroid barycentric position (in au) at bounce time (TDB)
        # rv_a_t_b = apophis_pv(t_r, -τ_D)
        r_a_t_b = x_r[apophisdofs[1:3]](-τ_D) # rv_a_t_b[1:3]
        v_a_t_b = x_r[apophisdofs[4:6]](-τ_D) # rv_a_t_b[4:6]
    end
    # get latest estimates of ρ_vec_r and ρ_r
    # Eq. (3) Yeomans et al. (1992)
    ρ_vec_r = r_a_t_b - r_r_t_r
    # Eq. (4) Yeomans et al. (1992)
    ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)

    # up-leg iteration
    # τ_U first estimation: Eq. (5) Yeomans et al. (1992)
    τ_U = τ_D
    @show τ_U
    # transmit time, 1st estimate Eq. (6) Yeomans et al. (1992)
    t_t = t_b - τ_U
    # convert transmit time to UTC time
    # TDB -> UTC
    t_t_utc = DateTime(UTCEpoch(TDBEpoch(constant_term(t_t), origin=:julian)))
    # Geocentric position/velocity of receiving antenna in inertial frame (au, au/day)
    R_t, V_t = observer_position(station_code, t_t_utc)
    # rv_e_t_t = earth_pv(t_r, -τ_U-τ_D)
    r_e_t_t = x_r[earthdofs[1:3]](-τ_U-τ_D) # rv_e_t_t[1:3]
    v_e_t_t = x_r[earthdofs[4:6]](-τ_U-τ_D) # rv_e_t_t[4:6]
    # Barycentric position and velocity of the transmitter at the transmit time
    r_t_t_t = r_e_t_t + R_t
    v_t_t_t = v_e_t_t + V_t
    # Eq. (7) Yeomans et al. (1992)
    ρ_vec_t = r_a_t_b - r_t_t_t
    ρ_t = sqrt(ρ_vec_t[1]^2 + ρ_vec_t[2]^2 + ρ_vec_t[3]^2)
    # Sun barycentric position (in au) at transmit time (TDB)
    r_s_t_t = x_r[sundofs[1:3]](-τ_U-τ_D) # sun_pv(t_r, -τ_U-τ_D)[1:3]
    for i in 1:niter
        # Eq. (8) Yeomans et al. (1992)
        τ_U = ρ_t/c_au_per_day # (days) -R_b/c (COM correction) + Δτ_U (relativistic, tropo, iono...)
        # Sun barycentric position and velocity (in au, au/day) at transmit time (TDB)
        r_s_t_t = x_r[sundofs[1:3]](-τ_U-τ_D) #sun_pv(t_r, -τ_U-τ_D)[1:3]
        v_s_t_t = x_r[sundofs[4:6]](-τ_U-τ_D) #sun_pv(t_r, -τ_U-τ_D)[4:6] # barycentric position of Sun at estimated bounce time
        # compute up-leg Shapiro delay
        e_U = norm(r_e_t_t-r_s_t_t) # heliocentric distance of Earth at t_t
        r_s_t_b = x_r[sundofs[1:3]](-τ_D) # sun_pv(t_r, -τ_D)[1:3] # barycentric position of Sun at estimated bounce time
        p_U = norm(r_a_t_b-r_s_t_b) # heliocentric distance of asteroid at t_b
        q_U = norm(ρ_vec_t) #signal path (up-leg)
        Δτ_rel_U = shapiro_delay(e_U, p_U, q_U) # days
        # corona_delay(r_e_t_t, r_a_t_b, t_t_jd1, t_t_jd2, f_T, station_code)
        Δτ_corona_U = corona_delay(r_e_t_t, r_a_t_b, r_s_t_t, v_s_t_t, f_T, station_code) # days
        Δτ_tropo_U = tropo_delay(R_t, ρ_vec_t) # days
        Δτ_U = Δτ_corona_U + (Δτ_rel_U + Δτ_tropo_U)
        τ_U = τ_U + Δτ_U
        @show τ_U
        @show Δτ_rel_U*(1e6*daysec)
        @show Δτ_corona_U*(1e6*daysec)
        @show Δτ_tropo_U*(1e6*daysec)
        # transmit time, 1st estimate Eq. (6) Yeomans et al. (1992)
        t_t = t_b - τ_U
        # convert transmit time to UTC time
        # TDB -> UTC
        t_t_utc = DateTime(UTCEpoch(TDBEpoch(constant_term(t_t), origin=:julian)))
        # Geocentric position/velocity of receiving antenna in inertial frame (au, au/day)
        R_t, V_t = observer_position(station_code, t_t_utc)
        # Earth's barycentric position and velocity at the transmit time
        # rv_e_t_t = earth_pv(t_r, -τ_U-τ_D)
        r_e_t_t = x_r[earthdofs[1:3]](-τ_U-τ_D) #rv_e_t_t[1:3]
        v_e_t_t = x_r[earthdofs[4:6]](-τ_U-τ_D) #rv_e_t_t[4:6]
        # Barycentric position and velocity of the transmitter at the transmit time
        r_t_t_t = r_e_t_t + R_t
        v_t_t_t = v_e_t_t + V_t
        # Eq. (7) Yeomans et al. (1992)
        ρ_vec_t = r_a_t_b - r_t_t_t
        ρ_t = sqrt(ρ_vec_t[1]^2 + ρ_vec_t[2]^2 + ρ_vec_t[3]^2)
    end
    @show τ_U

    # TDB -> UTC
    t_t_utc = UTCEpoch(TDBEpoch(constant_term(t_r), constant_term(-τ_D-τ_U), origin=:julian))
    # ut1_utc_t = EarthOrientation.getΔUT1( julian(t_t_utc).Δt ) # seconds
    utc_jd1_t, utc_jd2_t = julian_twopart(t_t_utc) # days, days
    j, tai_utc_t = iauDat(year(t_t_utc), month(t_t_utc), day(t_t_utc), utc_jd2_t.Δt)
    j != 0 && @warn "iauDat: j = $j. See SOFA.jl docs for more detail."
    tt_utc_t = 32.184 + tai_utc_t # seconds
    ##TT-TDB (transmit time)
    tt_tdb_t = tt_m_tdb(constant_term(t_t))[1] # seconds
    tdb_utc_t = tt_utc_t - tt_tdb_t #seconds
    @show tdb_utc_t

    # UTC two-part "quasi" Julian date (see http://www.iausofa.org/sofa_ts_c.pdf)
    utc_jd1_r, utc_jd2_r = julian_twopart(UTCEpoch(t_r_utc)) # days, days
    # ΔAT = TAI - UTC
    j, tai_utc_r = iauDat(year(t_r_utc), month(t_r_utc), day(t_r_utc), utc_jd2_r.Δt)
    j != 0 && @warn "iauDat: j = $j. See SOFA.jl docs for more detail."
    tt_utc_r = 32.184 + tai_utc_r
    tt_tdb_r = 0.0
    for i in 1:niter+1
        tt_tdb_r = tt_m_tdb(utc_jd1_r.Δt + utc_jd2_r.Δt + tt_utc_r/daysec + tt_tdb_r)[1] # seconds
        @show tt_tdb_r
    end
    tdb_utc_r = tt_utc_r - tt_tdb_r #seconds
    @show tdb_utc_r

    # compute total time delay (UTC seconds)
    # Eq. (9) Yeomans et al. (1992)
    #τ = τ_D + τ_U + (TDB-UTC)_t - (TDB-UTC)_r
    τ = daysec*(τ_U + τ_D) + (tdb_utc_t - tdb_utc_r)
    @show tdb_utc_t - tdb_utc_r
    @show τ

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

    @show ν
    return 1e6τ, 1e6ν # total signal delay (μs) and Doppler shift (Hz)
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
apophis_pv(jd1, jd2) = compute(eph, jd1, jd2, 2099942   , 12        , unitKM+unitDay, 1)/au # units: au, au/day
sun_pv(jd1, jd2)     = compute(eph, jd1, jd2, 11        , 12        , unitKM+unitDay, 1)/au # units: au, au/day
earth_pv(jd1, jd2)   = compute(eph, jd1, jd2, 3         , 12        , unitKM+unitDay, 1)/au # units: au, au/day
moon_pv(jd1, jd2)    = compute(eph, jd1, jd2, 10        , 12        , unitKM+unitDay, 1)/au # units: au, au/day
tt_m_tdb(jd1, jd2)   = compute(eph, jd1, jd2, 1000000001, 1000000000, unitSec       , 1) # units: seconds

tt_m_tdb(t) = tt_m_tdb(t, 0.0) # units: seconds

# Standard formula for relativistic (Shapiro) delay
function shapiro_delay(e, p, q)
    return 2μ[1]*log( abs((e+p+q)/(e+p-q)) )/(c_au_per_day^3) # days
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
function Ne(p1::Vector{S}, p2::Vector{S}, r_s_t0::Vector{T}, v_s_t0::Vector{T}, ds::U, ΔS::Real) where {T<:Number, S<:Number, U<:Number}
    # s: linear parametrization of ray path, such that
    # s=0 -> point on ray path is at p1; s=1 -> point on ray path is at p2
    s = ds/ΔS
    s_p2_p1 = map(x->s*x, (p2-p1))
    dt = constant_term(ds/c_cm_per_sec/daysec) # (TDB days)
    # Linear interpolation of sun position: r_sun(t0+dt) ≈ r_sun(t0)+dt*v_sun(t0)
    dt_v_s_t0 = map(x->dt*x, v_s_t0)
    r_s_t0_dt = r_s_t0 + dt_v_s_t0
    r_vec = p1+s_p2_p1 - r_s_t0_dt # heliocentric position (au) of point on ray path at time t_tdb_jul (Julian days)
    r = sqrt( r_vec[1]^2 + r_vec[2]^2 + r_vec[3]^2 ) # heliocentric distance (au) of point on ray path at time t_tdb_jul (Julian days)
    # compute ecliptic solar latitude of point on ray path
    # NOTE: actually, (Anderson, 1978) says it should be heliographic
    # (i.e., wrt solar equator), but diff is ~7 deg and things don't seem to change a lot
    # compute heliographic latitude (ecliptic plane is not equal to heliographic but almost by ~7deg)
    r_vec_ecliptic = r_vec # Rx(deg2rad(23.43929))*r_vec
    β = asin( r_vec_ecliptic[3]/r ) # ecliptic solar latitude of point on ray path (rad)
    r_sr = r/R_sun
    return (A_sun/r_sr^6) + ( (a_sun*b_sun)/sqrt((a_sun*sin(β))^2 + (b_sun*cos(β))^2) )/(r_sr^2)
end

# Integral of Ne, evaluated with TaylorIntegration.jl
# TODO: @taylorize!
# p1: signal departure point (transmitter/bounce) (au)
# p2: signal arrival point (bounce/receiver) (au)
# t_tdb_jd1, t_tdb_jd2: Two-part Julian date (TDB) of signal path (bounce time for downlink; transmit time for uplink)
# current_delay: total time-delay for current signal path (TDB days)
# output has units of electrons/cm^2
function Ne_path_integral(p1::Vector{S}, p2::Vector{S}, r_s_t0::Vector{T}, v_s_t0::Vector{T}) where {T<:Number, S<:Number}
    ΔS = (100000au)*norm(p2-p1) # total distance between p1 and p2, in centimeters
    # kernel of path integral; distance is in cm
    function int_kernel(s, x)
        return Ne(p1, p2, r_s_t0, v_s_t0, s, ΔS)
    end
    # do path integral; distance vector sv is in au; integral vector iv is in (electrons/cm^3)*au
    sv, iv = taylorinteg(int_kernel, 0.0, 0.0, ΔS, 28, 1e-20, maxsteps=5)
    @show size(iv)
    return iv[end] # (electrons/cm^2)
end

# Time-delay generated by thin plasma of solar corona
# Ne: ionized electrons density: electrons/cm^3
# p1: signal departure point (transmitter/bounce for up/down-link, resp.) (au)
# p2: signal arrival point (bounce/receiver for up/down-link, resp.) (au)
# t_tdb_jd1, t_tdb_jd2: Two-part Julian date (TDB) of signal path (bounce time for downlink; transmit time for uplink)
# f_T: transmitter frequency (MHz)
# From https://gssc.esa.int/navipedia/index.php/Ionospheric_Delay it seems that
# the expression 40.3*Ne/f^2 is adimensional, where Ne [electrons/meters^3] and f is in Hz
# therefore, the integral (40.3/f^2)*∫Ne*ds is in meters, where ds is in meters
# and the expression (40.3/(c*f^2))*∫Ne*ds, with c in meters/second, is in seconds
function corona_delay(p1::Vector{S}, p2::Vector{S}, r_s_t0::Vector{T}, v_s_t0::Vector{T}, f_T::U, station_code::Int) where {T<:Number, S<:Number, U<:Real}
    # for the time being, we're removing the terms associated with higher-order terms in the variationals (ie, Yarkovsky)
    int_path = Ne_path_integral(map(x->constant_term.(x), (p1, p2, r_s_t0, v_s_t0))...) # (electrons/cm^2)
    Δτ_corona = 40.3int_path/(c_cm_per_sec*f_T^2) # seconds
    # numerical factor k dependent on transmitter frequency (S-band, X-band) (Muhleman & Anderson, 1981 ApJ 247-1093)
    if issband(f_T)
        k = 5e-7
    elseif isxband(f_T)
        k = 5e-4
    end
    return k*Δτ_corona/daysec # seconds -> days
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
    return tropo_delay(zd)/daysec # days
end

# Alternate version of delay_doppler, using JPL ephemerides
# eph: Solar System + Apophis ephemerides pos/vel (au, au/day)
# station_code: observing station identifier (MPC nomenclature)
# t_r_utc_julian: time of echo reception (UTC)
# f_T: transmitter frequency (MHz)
function delay_doppler_jpleph(station_code::Int, t_r_utc::DateTime, f_T, niter::Int=10)
    # UTC -> TDB (receiving time)
    __t_r_jd1, __t_r_jd2 = julian_twopart(TDB, UTCEpoch(t_r_utc))
    t_r_jd1, t_r_jd2 = __t_r_jd1.Δt, __t_r_jd2.Δt
    @show t_r_jd1+t_r_jd2

    # Geocentric position/velocity of receiving antenna in inertial frame (au, au/day)
    R_r, V_r = observer_position(station_code, t_r_utc)

    # Earth's barycentric position and velocity at the receive time
    rv_e_t_r = earth_pv(t_r_jd1, t_r_jd2)
    r_e_t_r = rv_e_t_r[1:3]
    v_e_t_r = rv_e_t_r[4:6]

    # Barycentric position and velocity of the receiver at the receive time
    r_r_t_r = r_e_t_r + R_r
    v_r_t_r = v_e_t_r + V_r

    # asteroid barycentric position (in au) at receiving time (TDB)
    rv_a_t_r = apophis_pv(t_r_jd1, t_r_jd2)
    r_a_t_r = rv_a_t_r[1:3]
    v_a_t_r = rv_a_t_r[4:6]

    # Sun barycentric position (in au) at receiving time (TDB)
    r_s_t_r = sun_pv(t_r_jd1, t_r_jd2)[1:3]

    # down-leg iteration
    # τ_D first approximation: Eq. (1) Yeomans et al. (1992)
    ρ_vec_r = r_a_t_r - r_r_t_r
    ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)
    τ_D = ρ_r/c_au_per_day # (days) -R_b/c, but delay is wrt asteroid Center (Brozovic et al., 2018)
    # bounce time, 1st estimate, Eq. (2) Yeomans et al. (1992)
    t_b_jd1, t_b_jd2 = t_r_jd1, t_r_jd2-τ_D
    # asteroid barycentric position (in au) at bounce time (TDB)
    rv_a_t_b = apophis_pv(t_b_jd1, t_b_jd2)
    r_a_t_b = rv_a_t_b[1:3]
    v_a_t_b = rv_a_t_b[4:6]
    Δτ_rel_D = zero(τ_D)
    for i in 1:niter
        # Eq. (3) Yeomans et al. (1992)
        ρ_vec_r = r_a_t_b - r_r_t_r
        # Eq. (4) Yeomans et al. (1992)
        ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)
        τ_D = ρ_r/c_au_per_day # (days) -R_b/c (COM correction) + Δτ_D (relativistic, tropo, iono...)
        # compute down-leg Shapiro delay
        # NOTE: when using PPN, substitute 2 -> 1+γ in expressions for Shapiro delay, Δτ_rel_[D|U]
        e_D = norm(r_r_t_r-r_s_t_r) # norm(r_e_t_r-r_s_t_r) # heliocentric distance of station at t_r
        rv_s_t_b = sun_pv(t_b_jd1, t_b_jd2) # barycentric position and velocity of Sun at estimated bounce time
        r_s_t_b = rv_s_t_b[1:3]
        v_s_t_b = rv_s_t_b[4:6]
        p_D = norm(r_a_t_b-r_s_t_b) # norm(r_a_t_b-r_s_t_b) # heliocentric distance of asteroid at t_b
        q_D = norm(ρ_vec_r) #signal path (down-leg)
        Δτ_rel_D = shapiro_delay(e_D, p_D, q_D) # days
        Δτ_corona_D = corona_delay(r_a_t_b, r_e_t_r, r_s_t_b, v_s_t_b, f_T, station_code)
        Δτ_tropo_D = tropo_delay(R_r, ρ_vec_r)
        Δτ_D = Δτ_corona_D + (Δτ_rel_D + Δτ_tropo_D)
        τ_D = τ_D + Δτ_D
        @show τ_D
        @show Δτ_rel_D*(1e6*daysec)
        @show Δτ_corona_D*(1e6*daysec)
        @show Δτ_tropo_D*(1e6*daysec)
        # bounce time, new estimate Eq. (2) Yeomans et al. (1992)
        t_b_jd2 = t_r_jd2-τ_D
        # asteroid barycentric position (in au) at bounce time (TDB)
        rv_a_t_b = apophis_pv(t_b_jd1, t_b_jd2)
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
    @show τ_U
    # transmit time, 1st estimate Eq. (6) Yeomans et al. (1992)
    t_t_jd1, t_t_jd2 = t_r_jd1, t_r_jd2-(τ_U+τ_D)
    # convert transmit time to UTC time
    # TDB -> UTC
    t_t_utc = DateTime(UTCEpoch(TDBEpoch(t_t_jd1, t_t_jd2, origin=:julian)))
    # Geocentric position/velocity of receiving antenna in inertial frame (au, au/day)
    R_t, V_t = observer_position(station_code, t_t_utc)
    rv_e_t_t = earth_pv(t_r_jd1, t_r_jd2-(τ_U+τ_D))
    r_e_t_t = rv_e_t_t[1:3]
    v_e_t_t = rv_e_t_t[4:6]
    # Barycentric position and velocity of the transmitter at the transmit time
    r_t_t_t = r_e_t_t + R_t
    v_t_t_t = v_e_t_t + V_t
    # Eq. (7) Yeomans et al. (1992)
    ρ_vec_t = r_a_t_b - r_t_t_t
    ρ_t = sqrt(ρ_vec_t[1]^2 + ρ_vec_t[2]^2 + ρ_vec_t[3]^2)
    # Sun barycentric position (in au) at transmit time (TDB)
    r_s_t_t = sun_pv(t_r_jd1, t_r_jd2-(τ_U+τ_D))[1:3]
    Δτ_rel_U = zero(τ_U)
    for i in 1:niter
        # Eq. (8) Yeomans et al. (1992)
        τ_U = ρ_t/c_au_per_day # (days) -R_b/c (COM correction) + Δτ_U (relativistic, tropo, iono...)
        # Sun barycentric position and velocity (in au, au/day) at transmit time (TDB)
        rv_s_t_t = sun_pv(t_r_jd1, t_r_jd2-(τ_U-τ_D))
        r_s_t_t = rv_s_t_t[1:3]
        v_s_t_t = rv_s_t_t[4:6]
        # compute up-leg Shapiro delay
        e_U = norm(r_t_t_t-r_s_t_t) # heliocentric distance of station at t_t
        r_s_t_b = sun_pv(t_b_jd1, t_b_jd2)[1:3] # barycentric position of Sun at bounce time
        p_U = norm(r_a_t_b-r_s_t_b) # heliocentric distance of asteroid at t_b
        q_U = norm(ρ_vec_t) #signal path (up-leg)
        Δτ_rel_U = shapiro_delay(e_U, p_U, q_U) # days
        Δτ_corona_U = corona_delay(r_e_t_t, r_a_t_b, r_s_t_t, v_s_t_t, f_T, station_code) # days
        Δτ_tropo_U = tropo_delay(R_t, ρ_vec_t) # days
        Δτ_U = Δτ_corona_U + (Δτ_rel_U + Δτ_tropo_U)
        τ_U = τ_U + Δτ_U
        @show τ_U
        @show Δτ_rel_U*(1e6*daysec)
        @show Δτ_corona_U*(1e6*daysec)
        @show Δτ_tropo_U*(1e6*daysec)
        # transmit time, 1st estimate Eq. (6) Yeomans et al. (1992)
        t_t_jd1, t_t_jd2 = t_r_jd1, t_r_jd2-(τ_U+τ_D)
        # convert transmit time to UTC time
        # TDB -> UTC
        t_t_utc = DateTime(UTCEpoch(TDBEpoch(t_t_jd1, t_t_jd2, origin=:julian)))
        # Geocentric position/velocity of receiving antenna in inertial frame (au, au/day)
        R_t, V_t = observer_position(station_code, t_t_utc)
        # Earth's barycentric position and velocity at the transmit time
        rv_e_t_t = earth_pv(t_r_jd1, t_r_jd2-(τ_U+τ_D))
        r_e_t_t = rv_e_t_t[1:3]
        v_e_t_t = rv_e_t_t[4:6]
        # Barycentric position and velocity of the transmitter at the transmit time
        r_t_t_t = r_e_t_t + R_t
        v_t_t_t = v_e_t_t + V_t
        # Eq. (7) Yeomans et al. (1992)
        ρ_vec_t = r_a_t_b - r_t_t_t
        ρ_t = sqrt(ρ_vec_t[1]^2 + ρ_vec_t[2]^2 + ρ_vec_t[3]^2)
    end

    # compute TDB-UTC at transmit time
    # TDB -> UTC
    t_t_utc = UTCEpoch(TDBEpoch(t_t_jd1, t_t_jd2, origin=:julian))
    # ut1_utc_t = EarthOrientation.getΔUT1( julian(t_t_utc).Δt ) # seconds
    utc_jd1_t, utc_jd2_t = julian_twopart(t_t_utc) # days, days
    @show utc_jd2_t.Δt
    j, tai_utc_t = iauDat(year(t_t_utc), month(t_t_utc), day(t_t_utc), utc_jd2_t.Δt)
    j != 0 && @warn "iauDat: j = $j. See SOFA.jl docs for more detail."
    tt_utc_t = 32.184 + tai_utc_t # seconds
    ##TT-TDB (transmit time)
    tt_tdb_t = tt_m_tdb(t_t_jd1, t_t_jd2)[1] # seconds
    tdb_utc_t = tt_utc_t - tt_tdb_t #seconds
    @show tdb_utc_t

    # compute TDB-UTC at receive time
    # UTC two-part "quasi" Julian date (see http://www.iausofa.org/sofa_ts_c.pdf)
    utc_jd1_r, utc_jd2_r = julian_twopart(UTCEpoch(t_r_utc)) # days, days
    # ΔAT = TAI - UTC
    j, tai_utc_r = iauDat(year(t_r_utc), month(t_r_utc), day(t_r_utc), utc_jd2_r.Δt)
    j != 0 && @warn "iauDat: j = $j. See SOFA.jl docs for more detail."
    tt_utc_r = 32.184 + tai_utc_r
    tt_tdb_r = 0.0
    for i in 1:niter+1
        tt_tdb_r = tt_m_tdb(utc_jd1_r.Δt, utc_jd2_r.Δt + tt_utc_r/daysec + tt_tdb_r)[1] # seconds
        @show tt_tdb_r
    end
    tdb_utc_r = tt_utc_r - tt_tdb_r #seconds
    @show tdb_utc_r

    # compute total time delay (UTC seconds)
    # Eq. (9) Yeomans et al. (1992)
    #τ = τ_D + τ_U + (TDB-UTC)_t - (TDB-UTC)_r
    @show daysec*(τ_U + τ_D) + (tdb_utc_t - tdb_utc_r)
    τ = (τ_U + τ_D) + (tdb_utc_t - tdb_utc_r)/daysec
    τ *= daysec
    @show tdb_utc_t - tdb_utc_r
    @show τ

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

    return 1e6τ, 1e6ν # total signal delay (μs) and Doppler shift (Hz)
end
