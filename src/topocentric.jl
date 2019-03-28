function observer_position(station_code, t_utc::T) where {T<:Number}
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
    x_gc = r_cos_phi*cos(east_long_rad)
    y_gc = r_cos_phi*sin(east_long_rad)
    z_gc = r_sin_phi*one(east_long_rad)

    # Apply rotation from geocentric, Earth-fixed frame to inertial (celestial) frame
    W = inv(sofa_c2t_rotation(t_utc))

    vec_pos_pod = [x_gc, y_gc, z_gc]
    vec_pos_frame = W*vec_pos_pod

    return vec_pos_frame
end

# conversion of micro-arcseconds to radians
mas2rad(x) = deg2rad(x/3.6e6) # mas/1000 -> arcsec; arcsec/3600 -> deg; deg2rad(deg) -> rad

# Celestial-to-terrestrial rotation matrix (including polar motion)
# Reproduction of Section 5.3 of SOFA Tools for Earth Attitude
# "IAU 2000A, CIO based, using classical angles"
# found at SOFA website, Mar 27, 2019
# Some modifications were applied, using TaylorSeries.jl, in order to compute
# the derivative of the observer position, following the guidelines from
# ESAA 2014, Sec 7.4.3.3 (page 295)
function sofa_c2t_rotation(t_utc::Taylor1{T}) where {T<:Real}
    # UTC
    t0_utc = UTCEpoch(constant_term(t_utc), origin=:julian)
    t0_utc_jul = julian(t0_utc)

    # Polar motion (arcsec->radians)
    xp_arcsec, yp_arcsec = EarthOrientation.polarmotion(t0_utc_jul.Δt)
    xp = deg2rad(xp_arcsec/3600)
    yp = deg2rad(yp_arcsec/3600)

    # UT1-UTC (s)
    dut1 = EarthOrientation.getΔUT1(t0_utc_jul.Δt) # seconds

    # CIP offsets wrt IAU 2000A (mas->radians)
    dx00_mas, dy00_mas = EarthOrientation.precession_nutation00(t0_utc_jul.Δt)
    dx00 = mas2rad(dx00_mas)
    dy00 = mas2rad(dy00_mas)

    # TT
    t0_tt = TTEpoch(t0_utc)
    t0_tt_jd1, t0_tt_jd2 = julian_twopart(t0_tt)

    # UT1
    t0_ut1 = UT1Epoch(t0_utc)
    t0_ut1_jd1, t0_ut1_jd2 = julian_twopart(t0_ut1)

    # CIP and CIO, IAU 2000A
    x, y, s = iauXys00a( t0_tt_jd1.Δt, t0_tt_jd2.Δt )

    # Add CIP corrections
    x += dx00
    y += dy00

    # GCRS to CIRS matrix
    rc2i = iauC2ixys( x, y, s)

    # Earth rotation angle
    era = iauEra00( t0_ut1_jd1.Δt, t0_ut1_jd2.Δt ) #rad
    eraT1 = era + ω*Taylor1(t_utc.order) #rad/day

    # Form celestial-terrestrial matrix (no polar motion yet)
    rc2ti = iauCr(rc2i)
    Rz_eraT1 = [cos(eraT1) sin(eraT1) zero(eraT1);
        -sin(eraT1) cos(eraT1) zero(eraT1);
        zero(eraT1) zero(eraT1) one(eraT1)
    ]
    rc2tiT1 = Rz_eraT1*rc2ti

    # Polar motion matrix (TIRS->ITRS, IERS 2003)
    sp = iauSp00( t0_tt_jd1.Δt, t0_tt_jd2.Δt )
    rpom = iauPom00( xp, yp, sp)

    # Form celestial-terrestrial matrix (including polar motion)
    rc2itT1 = rpom*rc2tiT1
    return rc2itT1
end
