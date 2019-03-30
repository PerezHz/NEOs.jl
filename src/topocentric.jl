function observer_position(station_code, t_utc::T) where {T<:Real}
    # east_long: East longitude (deg)
    # r_cos_phi: distance from spin axis (km), taken from Yeomans et al. (1992)
    # r_sin_phi: height above equatorial plane (km), taken from Yeomans et al. (1992)
    # East long data is consistent with MPC database (2019-Mar-7)
    # Goldstone site numbers 13, 14 seem to be interchanged in Yeomans et al. (1992) paper
    # TODO: add more radar stations
    if station_code == 251 # Arecibo
        east_long = 293.24692 #deg
        r_cos_phi = 6056.525 #km
        r_sin_phi = 1994.665 #km
    elseif station_code == 252 # Goldstone DSS 13 (Venus site), Fort Irwin
        east_long = 243.20512 #deg
        r_cos_phi = 5215.484 #km
        r_sin_phi = 3660.957 #km
    elseif station_code == 253 # Goldstone DSS 14 (Mars site), Fort Irwin
        east_long = 243.11047 #deg
        r_cos_phi = 5203.997 #km
        r_sin_phi = 3677.052 #km
    elseif station_code == 254 # Haystack, Westford, MA
        east_long = 288.51128 #deg
        r_cos_phi = 4700.514 #km
        r_sin_phi = 4296.900 #km
    else
        @error "Unknown station."
    end
    # cartesian components of Earth-fixed position of observer
    east_long_rad = deg2rad(east_long)
    x_gc = r_cos_phi*cos(east_long_rad) #km
    y_gc = r_cos_phi*sin(east_long_rad) #km
    z_gc = r_sin_phi*one(east_long_rad) #km

    pos_geo = [x_gc, y_gc, z_gc]/au #au

    # Apply rotation from geocentric, Earth-fixed frame to inertial (celestial) frame
    return c2t_rotation_iau_00_06(t_utc, pos_geo)
end

# conversion of micro-arcseconds to radians
mas2rad(x) = deg2rad(x/3.6e6) # mas/1000 -> arcsec; arcsec/3600 -> deg; deg2rad(deg) -> rad

# Celestial-to-terrestrial rotation matrix (including polar motion)
# Reproduction of Section 5.3 of SOFA Tools for Earth Attitude
# "IAU 2000A, CIO based, using classical angles"
# found at SOFA website, Mar 27, 2019
# Some modifications were applied, using TaylorSeries.jl, in order to compute
# the geocentric velocity of the observer, following the guidelines from
# ESAA 2014, Sec 7.4.3.3 (page 295)
function c2t_rotation_iau_00_06(t_utc::T, pos_geo::Vector{S}) where {T<:Real, S<:Real}
    # UTC
    t0_utc = UTCEpoch(t_utc, origin=:julian)
    t0_utc_jul = julian(t0_utc)

    # UT1
    # dut1 = EarthOrientation.getΔUT1(t0_utc_jul.Δt) # UT1-UTC (seconds)
    t0_ut1 = UT1Epoch(t0_utc)
    t0_ut1_jd1, t0_ut1_jd2 = julian_twopart(t0_ut1)
    # Earth rotation angle ERA = ERA0 + ω*(ERA - ERA0)
    era = iauEra00( t0_ut1_jd1.Δt, t0_ut1_jd2.Δt ) #rad
    # this trick allows us to compute the whole Celestial->Terrestrial matrix and its first derivative
    # For more details, see ESAA 2014, p. 295, Sec. 7.4.3.3, Eqs. 7.137-7.140
    eraT1 = era + ω*Taylor1(1) #rad/day
    # Rz(-ERA)
    Rz_minus_era_T1 = [cos(eraT1) sin(-eraT1) zero(eraT1);
        sin(eraT1) cos(eraT1) zero(eraT1);
        zero(eraT1) zero(eraT1) one(eraT1)
    ]
    Rz_minus_ERA = Rz_minus_era_T1()
    # dRz(-ERA)/dt
    dRz_minus_era_T1 = differentiate.(Rz_minus_era_T1)
    dRz_minus_ERA = dRz_minus_era_T1()

    # TT
    t0_tt = TTEpoch(t0_utc)
    t0_tt_jd1, t0_tt_jd2 = julian_twopart(t0_tt)
    # Polar motion (arcsec->radians)
    xp_arcsec, yp_arcsec = EarthOrientation.polarmotion(t0_utc_jul.Δt)
    xp = deg2rad(xp_arcsec/3600)
    yp = deg2rad(yp_arcsec/3600)
    # Polar motion matrix (TIRS->ITRS, IERS 2003)
    sp = iauSp00( t0_tt_jd1.Δt, t0_tt_jd2.Δt )
    W = iauPom00( xp, yp, sp)
    W_inv = inv(W)

    # CIP and CIO, IAU 2000A
    x, y, s = iauXys00a( t0_tt_jd1.Δt, t0_tt_jd2.Δt )
    # CIP offsets wrt IAU 2000A (mas->radians)
    dx00_mas, dy00_mas = EarthOrientation.precession_nutation00(t0_utc_jul.Δt)
    dx00 = mas2rad(dx00_mas)
    dy00 = mas2rad(dy00_mas)
    # Add CIP corrections
    x += dx00
    y += dy00
    # GCRS to CIRS matrix
    C = iauC2ixys( x, y, s)
    C_inv = inv(C) # CIRS -> GCRS

    # g(t), \dot g(t) ESAA vectors
     g_vec_ESAA =  Rz_minus_ERA*(W_inv*pos_geo)
    dg_vec_ESAA = dRz_minus_ERA*(W_inv*pos_geo)

    # G(t), \dot G(t) ESAA vectors
     G_vec_ESAA = C_inv* g_vec_ESAA
    dG_vec_ESAA = C_inv*dg_vec_ESAA

     return G_vec_ESAA, dG_vec_ESAA
end
