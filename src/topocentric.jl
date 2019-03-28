#eastlong must be in degrees; result is in degrees
function localsiderealtime(jd::Real, east_long)
    J0 = floor(jd)-0.5 # Julian date at 0 h UT
    T0 = (J0-2.451545e6)/36525 # Julian centuries since J2000.0 = 2.451545e6
    # Greenwich sidereal time at T0 (degrees)
    θG0 = mod(100.4606184 + 36000.77004T0 + 0.000387933T0^2 - 2.583E-8T0^3, 360.0)
    UT = (jd-J0) # UT (in units of days)
    θG = θG0 + 360.98564724UT # Greenwich sidereal time at jd (degrees)
    return mod(θG + east_long, 360.0)
end

function localsiderealtime(jd::Taylor1{T}, east_long) where {T<:Real}
    J0 = Taylor1([floor(constant_term(jd))-0.5, one(eltype(jd))], jd.order) # Julian date at 0 h UT
    T0 = (J0-2.451545e6)/36525 # Julian centuries since J2000.0 = 2.451545e6
    # Greenwich sidereal time at T0 (degrees)
    θG0 = mod(100.4606184 + 36000.77004T0 + 0.000387933T0^2 - 2.583E-8T0^3, 360.0)
    UT = (jd-J0) # UT (in units of days)
    θG = θG0 + 360.98564724UT # Greenwich sidereal time at jd (degrees)
    return mod(θG + east_long, 360.0)
end

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
    # local sidereal time of observation
    #lmst = deg2rad( localsiderealtime(t_utc, east_long) )
    # cartesian components of topocentric position of observer
    east_long_rad = deg2rad(east_long)
    x_gc = r_cos_phi*cos(east_long_rad) #*cos(lmst)
    y_gc = r_cos_phi*sin(east_long_rad) #*sin(lmst)
    z_gc = r_sin_phi*one(east_long_rad) #*one(lmst)

    # rotate observer topocentric position from Earth frame to inertial frame
    # t_tdb = julian(TDBEpoch(UTCEpoch(t_utc, origin=:julian))).Δt
    # W = inv(earth_pole_rotation(t_tdb-J2000))
    W = inv(sofa_c2t_rotation(t_utc))

    vec_pos_pod = [x_gc, y_gc, z_gc]
    vec_pos_frame = W*vec_pos_pod
    # @show vec_pos_pod
    # @show vec_pos_frame

    return vec_pos_frame
end

function observer_position(station_code, t_utc::Taylor1{T}) where {T<:Number}
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
    # local sidereal time of observation
    #lmst = deg2rad( localsiderealtime(t_utc, east_long) )
    # cartesian components of topocentric position of observer
    east_long_rad = deg2rad(east_long)
    x_gc = r_cos_phi*cos(east_long_rad) #*cos(lmst)
    y_gc = r_cos_phi*sin(east_long_rad) #*sin(lmst)
    z_gc = r_sin_phi*one(east_long_rad) #*one(lmst)

    # rotate observer topocentric position from Earth frame to inertial frame
    #t_tdb = julian(TDBEpoch(UTCEpoch(constant_term(t_utc), origin=:julian))).Δt
    #W = inv( earth_pole_rotation(Taylor1([t_tdb, one(t_tdb)], t_utc.order)-J2000) )
    #W2 = inv(sofa_c2t_rotation(constant_term(t_utc)))
    W = inv(sofa_c2t_rotation(t_utc))
    #@show constant_term.(W)
    #@show W2

    vec_pos_pod = [x_gc, y_gc, z_gc]
    vec_pos_frame = W*vec_pos_pod

    return vec_pos_frame
end

# conversion of micro-arcseconds to radians
mas2rad(x) = (5/18)*deg2rad(x) # 5/18 = 1000/3600

# Celestial-to-terrestrial rotation matrix (including polar motion)
# Reproduction of Section 5.3 of SOFA Tools for Earth Attitude
# "IAU 2000A, CIO based, using classical angles"
# found at SOFA website, Mar 27, 2019
function sofa_c2t_rotation(t_utc::T) where {T<:Real}
    # UTC
    # iy = 2007
    # im = 4
    # id = 5
    # ih = 12
    # min = 0
    # sec = 0.0
    t0_utc = UTCEpoch(t_utc, origin=:julian)
    t0_utc_jul = julian(t0_utc)

    # Polar motion (arcsec->radians)
    xp_arcsec, yp_arcsec = polarmotion(t0_utc_jul.Δt)
    xp = deg2rad(xp_arcsec/3600)
    yp = deg2rad(yp_arcsec/3600)

    # # UT1-UTC (s)
    # dut1 = getΔUT1(t0_utc_jul.Δt)

    # # Nutation corrections wrt IAU 1976/1980 (mas->radians)
    # ddp80_mas, dde80_mas = precession_nutation80(t0_utc_jul.Δt)
    # ddp80 = mas2rad(ddp80_mas) # -55.0655 * DMAS2R;
    # dde80 = mas2rad(dde80_mas) # -6.3580 * DMAS2R;

    # CIP offsets wrt IAU 2000A (mas->radians)
    dx00_mas, dy00_mas = precession_nutation00(t0_utc_jul.Δt)
    dx00 = mas2rad(dx00_mas) # 0.1725 * DMAS2R;
    dy00 = mas2rad(dy00_mas) # -0.2650 * DMAS2R;

    # CIP offsets wrt IAU 2006/2000A (mas->radians)
    ### ***NOT CURRENTLY AVAILABLE FROM EarthOrientation.jl*** ###

    # TT (MJD)
    # t0_tt = TTEpoch(t0_utc)
    djmjd0, tt = julian_twopart(TT, t0_utc)

    # UT1
    t0_ut1 = UT1Epoch(t0_utc)

    # CIP and CIO, IAU 2000A
    #x, y, s_rad = iauXys00a( djmjd0.Δt, tt.Δt )
    #s = rad2deg(s_rad)*3600 #arcsec
    x, y, s = iauXys00a( djmjd0.Δt, tt.Δt )

    # Add CIP corrections
    #x += dx00
    #y += dy00

    # GCRS to CIRS matrix
    rc2i = iauC2ixys( x, y, s)

    # Earth rotation angle
    t0_ut1_jd1, t0_ut1_jd2 = julian_twopart(t0_ut1)
    #era_rad = iauEra00( t0_ut1_jd1.Δt, t0_ut1_jd2.Δt )
    #era = rad2deg(era_rad)
    era = iauEra00( t0_ut1_jd1.Δt, t0_ut1_jd2.Δt ) #rad
    # eraT1 = era + ω*Taylor1(t_utc.order)

    # Form celestial-terrestrial matrix (no polar motion yet)
    rc2ti = iauCr(rc2i)
    # Rz_eraT1 = [cos(eraT1) sin(eraT1) zero(eraT1); -sin(eraT1) cos(eraT1) zero(eraT1); zero(eraT1) zero(eraT1) one(eraT1)]
    # rc2tiT1 = Rz_eraT1*rc2ti
    rc2ti = iauRz( era, rc2ti )

    # Polar motion matrix (TIRS->ITRS, IERS 2003)
    sp = iauSp00( djmjd0.Δt, tt.Δt )
    rpom = iauPom00( xp, yp, sp)

    # Form celestial-terrestrial matrix (including polar motion)
    # rc2itT1 = rpom*rc2tiT1
    rc2it = iauRxr( rpom, rc2ti )
    return rc2it
end

function sofa_c2t_rotation(t_utc::Taylor1{T}) where {T<:Real}
    # UTC
    # iy = 2007
    # im = 4
    # id = 5
    # ih = 12
    # min = 0
    # sec = 0.0
    t0_utc = UTCEpoch(constant_term(t_utc), origin=:julian)
    t0_utc_jul = julian(t0_utc)

    # Polar motion (arcsec->radians)
    xp_arcsec, yp_arcsec = polarmotion(t0_utc_jul.Δt)
    xp = deg2rad(xp_arcsec/3600)
    yp = deg2rad(yp_arcsec/3600)

    # # UT1-UTC (s)
    # dut1 = getΔUT1(t0_utc_jul.Δt)

    # # Nutation corrections wrt IAU 1976/1980 (mas->radians)
    # ddp80_mas, dde80_mas = precession_nutation80(t0_utc_jul.Δt)
    # ddp80 = mas2rad(ddp80_mas) # -55.0655 * DMAS2R;
    # dde80 = mas2rad(dde80_mas) # -6.3580 * DMAS2R;

    # CIP offsets wrt IAU 2000A (mas->radians)
    dx00_mas, dy00_mas = precession_nutation00(t0_utc_jul.Δt)
    dx00 = mas2rad(dx00_mas) # 0.1725 * DMAS2R;
    dy00 = mas2rad(dy00_mas) # -0.2650 * DMAS2R;

    # CIP offsets wrt IAU 2006/2000A (mas->radians)
    ### ***NOT CURRENTLY AVAILABLE FROM EarthOrientation.jl*** ###

    # TT (MJD)
    # t0_tt = TTEpoch(t0_utc)
    djmjd0, tt = julian_twopart(TT, t0_utc)

    # UT1
    t0_ut1 = UT1Epoch(t0_utc)

    # CIP and CIO, IAU 2000A
    #x, y, s_rad = iauXys00a( djmjd0.Δt, tt.Δt )
    #s = rad2deg(s_rad)*3600 #arcsec
    x, y, s = iauXys00a( djmjd0.Δt, tt.Δt )

    # Add CIP corrections
    #x += dx00
    #y += dy00

    # GCRS to CIRS matrix
    rc2i = iauC2ixys( x, y, s)

    # Earth rotation angle
    t0_ut1_jd1, t0_ut1_jd2 = julian_twopart(t0_ut1)
    #era_rad = iauEra00( t0_ut1_jd1.Δt, t0_ut1_jd2.Δt )
    #era = rad2deg(era_rad)
    era = iauEra00( t0_ut1_jd1.Δt, t0_ut1_jd2.Δt ) #rad
    eraT1 = era + ω*Taylor1(t_utc.order)

    # Form celestial-terrestrial matrix (no polar motion yet)
    rc2ti = iauCr(rc2i)
    Rz_eraT1 = [cos(eraT1) sin(eraT1) zero(eraT1); -sin(eraT1) cos(eraT1) zero(eraT1); zero(eraT1) zero(eraT1) one(eraT1)]
    # rc2ti = iauRz( era, rc2ti )
    rc2tiT1 = Rz_eraT1*rc2ti

    # Polar motion matrix (TIRS->ITRS, IERS 2003)
    sp = iauSp00( djmjd0.Δt, tt.Δt )
    rpom = iauPom00( xp, yp, sp)

    # Form celestial-terrestrial matrix (including polar motion)
    # rc2it = iauRxr( rpom, rc2ti )
    rc2itT1 = rpom*rc2tiT1
    return rc2itT1
end