function observer_position(station_code::Int, t_utc::DateTime, pm::Bool=true)
    # λ_deg: East longitude (deg)
    # u: distance from spin axis (km), taken from Yeomans et al. (1992) u = r*cos(ϕ)
    # v: height above equatorial plane (km), taken from Yeomans et al. (1992) v = r*sin(ϕ)
    # East long data is consistent with MPC database (2019-Mar-7)
    # TODO: add more radar stations
    if station_code == 251 # Arecibo
        λ_deg = 293.24692 #deg
        u = 6056.525 #km #6056.497 #km
        v = 1994.665 #km #1994.649 #km
    elseif station_code == 252 # Goldstone DSS 13 (Venus site), Fort Irwin
        λ_deg = 243.20512 #deg
        u = 5215.484 #km
        v = 3660.957 #km
    elseif station_code == 253 # Goldstone DSS 14 (Mars site), Fort Irwin
        λ_deg = 243.11047 #deg
        u = 5203.997 #km
        v = 3677.052 #km
    elseif station_code == 254 # Haystack, Westford, MA
        λ_deg = 288.51128 #deg
        u = 4700.514 #km
        v = 4296.900 #km
    else
        @error "Unknown station."
    end
    # cartesian components of Earth-fixed position of observer
    λ_rad = deg2rad(λ_deg) # rad
    x_gc = u*cos(λ_rad) #km
    y_gc = u*sin(λ_rad) #km
    z_gc = v #km

    pos_geo = [x_gc, y_gc, z_gc]/au #au

    G_vec_ESAA, dG_vec_ESAA, gast = t2c_rotation_iau_76_80(t_utc, pos_geo, pm)
    # G_vec_ESAA, dG_vec_ESAA, era = t2c_rotation_iau_00_06(t_utc, pos_geo, pm)

    # Apply rotation from geocentric, Earth-fixed frame to inertial (celestial) frame
    return G_vec_ESAA, dG_vec_ESAA
end

# conversion of milli-arcseconds to radians
mas2rad(x) = deg2rad(x/3.6e6) # mas/1000 -> arcsec; arcsec/3600 -> deg; deg2rad(deg) -> rad

# Terrestrial-to-celestial rotation matrix (including polar motion)
# Reproduction of Section 5.3 of SOFA Tools for Earth Attitude
# "IAU 2000A, CIO based, using classical angles"
# found at SOFA website, Mar 27, 2019
# Some modifications were applied, using TaylorSeries.jl, in order to compute
# the geocentric velocity of the observer, following the guidelines from
# ESAA 2014, Sec 7.4.3.3 (page 295)
function t2c_rotation_iau_00_06(t_utc::DateTime, pos_geo::Vector, pm::Bool=true)
    # UTC
    t0_utc = UTCEpoch(t_utc)
    t0_utc_jul = datetime2julian(t_utc)

    # UT1
    # dut1 = EarthOrientation.getΔUT1(t_utc) # UT1-UTC (seconds)
    t0_ut1 = UT1Epoch(t0_utc)
    t0_ut1_jd1, t0_ut1_jd2 = julian_twopart(t0_ut1)
    # Earth rotation angle
    era = iauEra00( t0_ut1_jd1.Δt, t0_ut1_jd2.Δt ) #rad
    # this trick allows us to compute the whole Celestial->Terrestrial matrix and its first derivative
    # For more details, see ESAA 2014, p. 295, Sec. 7.4.3.3, Eqs. 7.137-7.140
    eraT1 = era + ω*Taylor1(1) #rad/day
    # eraT1 = era + omega(getlod(t_utc))*Taylor1(1) #rad/day
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
    W = Array{Float64}(I, 3, 3)
    if pm
        xp_arcsec, yp_arcsec = EarthOrientation.polarmotion(t_utc)
        xp = deg2rad(xp_arcsec/3600)
        yp = deg2rad(yp_arcsec/3600)
        # Polar motion matrix (TIRS->ITRS, IERS 2003)
        sp = iauSp00( t0_tt_jd1.Δt, t0_tt_jd2.Δt )
        W = iauPom00( xp, yp, sp)
    end
    W_inv = inv(W)

    # CIP and CIO, IAU 2000A
    x, y, s = iauXys00a( t0_tt_jd1.Δt, t0_tt_jd2.Δt )
    # CIP offsets wrt IAU 2000A (mas->radians)
    dx00_mas, dy00_mas = EarthOrientation.precession_nutation00(t_utc)
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
     G_vec_ESAA = convert(Vector{Float64}, C_inv* g_vec_ESAA)
    dG_vec_ESAA = convert(Vector{Float64}, C_inv*dg_vec_ESAA)

    return G_vec_ESAA, dG_vec_ESAA, era
end

# Terrestrial-to-celestial rotation matrix (including polar motion)
# Using 1976/1980 Earth orientation/rotation model
# Reproduction of Section 5.2 of SOFA Tools for Earth Attitude
# "IAU 1976/1980/1982/1994, equinox based"
# found at SOFA website, Mar 27, 2019
# Some modifications were applied, using TaylorSeries.jl, in order to compute
# the geocentric velocity of the observer, following the guidelines from
# ESAA 2014, Sec 7.4.3.3 (page 295)
function t2c_rotation_iau_76_80(t_utc::DateTime, pos_geo::Vector, pm::Bool=true)
    # UTC
    t0_utc = UTCEpoch(t_utc)
    # TT
    t0_tt = TTEpoch(t0_utc)
    djmjd0 = julian(J2000_EPOCH).Δt # J2000.0 (TT) epoch (days)
    tt = AstroTime.j2000(t0_tt).Δt
    # IAU 1976 precession matrix, J2000.0 to date
    # https://github.com/sisl/SOFA.jl/blob/be9ddfd412c5ab77b291b17decfd369041ef365b/src/pmat76.jl#L14
    rp = iauPmat76(djmjd0, tt)
    # IAU 1980 nutation angles Δψ (nutation in longitude), Δϵ (nutation in obliquity)
    # Output of `SOFA.iauNut80` is in radians:
    # https://github.com/sisl/SOFA.jl/blob/dc911b990dba79435399e1af0206e4acfc94c630/src/nut80.jl#L14
    dp80, de80  = iauNut80(djmjd0, tt)
    #Nutation corrections wrt IAU 1976/1980
    # Output of `EarthOrientation.precession_nutation80` is in mas:
    # https://github.com/JuliaAstro/EarthOrientation.jl/blob/529f12425a6331b133f989443aeb3fbbafd8f324/src/EarthOrientation.jl#L413
    ddp80_mas, dde80_mas = EarthOrientation.precession_nutation80(t_utc)
    # Convert mas -> radians
    ddp80 = mas2rad(ddp80_mas)
    dde80 = mas2rad(dde80_mas)
    # Add adjustments: frame bias, precession-rates, geophysical
    dpsi = dp80 + ddp80 #rad
    deps = de80 + dde80 #rad
    # Mean obliquity (output in radians)
    epsa = iauObl80(djmjd0, tt) # rad
    # Form the matrix of nutation
    rn = iauNumat(epsa, dpsi, deps)
    # Combine the matrices:  PN = N x P (with frame-bias B included)
    C = rn*rp

    # Equation of the equinoxes `ee = GAST - GMST`, including nutation correction
    ee = iauEqeq94(djmjd0, tt) + ddp80*cos(epsa)
    # UT1
    # dut1 = EarthOrientation.getΔUT1(t0_utc_jul.Δt) # UT1-UTC (seconds)
    t0_ut1 = UT1Epoch(t0_utc)
    djmjd0_plus_date, tut = julian_twopart(t0_ut1)
    # Greenwich apparent sidereal time (IAU 1982/1994)
    gast = iauAnp( iauGmst82( djmjd0_plus_date.Δt, tut.Δt ) + ee ) #rad
    # this trick allows us to compute the whole Celestial->Terrestrial matrix and its first derivative
    # For more details, see ESAA 2014, p. 295, Sec. 7.4.3.3, Eqs. 7.137-7.140
    gastT1 = gast + ω*Taylor1(1) #rad/day
    # gastT1 = gast + omega(EarthOrientation.getlod(t_utc))*Taylor1(1) #rad/day
    # Rz(-GAST)
    Rz_minus_gast_T1 = [cos(gastT1) -sin(gastT1) zero(gastT1);
        sin(gastT1) cos(gastT1) zero(gastT1);
        zero(gastT1) zero(gastT1) one(gastT1)]
    # Rz_minus_gast_T1 = Rz(-gastT1)
    Rz_minus_GAST = Rz_minus_gast_T1()
    # dRz(-GAST)/dt
    dRz_minus_gast_T1 = differentiate.(Rz_minus_gast_T1)
    dRz_minus_GAST = dRz_minus_gast_T1()
    # Rz_minus_GAST = [cos(gast) -sin(gast) 0.0;
    #     sin(gast) cos(gast) 0.0;
    #     0.0 0.0 1.0]
    # dRz_minus_GAST = ω*[-sin(gast) -cos(gast) 0.0;
    #     cos(gast) -sin(gast) 0.0;
    #     0.0 0.0 0.0]

    # # Form celestial-terrestrial matrix (no polar motion yet)
    # rc2ti = deepcopy(C)
    # rc2ti = iauRz(gast, rc2ti)

    # Polar motion matrix (TIRS->ITRS, IERS 1996)
    W = Array{Float64}(I, 3, 3)
    # Polar motion (arcsec->radians)
    if pm
        xp_arcsec, yp_arcsec = EarthOrientation.polarmotion(t_utc)
        xp = deg2rad(xp_arcsec/3600)
        yp = deg2rad(yp_arcsec/3600)
        W = iauRx(-yp, W)
        W = iauRy(-xp, W)
    end

    # # Form celestial-terrestrial matrix (including polar motion)
    # rc2it = W*rc2ti

    W_inv = transpose(W)
    C_inv = transpose(C)

    # g(t), \dot g(t) ESAA vectors
    g_vec_ESAA =  Rz_minus_GAST*(W_inv*pos_geo)
    dg_vec_ESAA = dRz_minus_GAST*(W_inv*pos_geo)
    # G(t), \dot G(t) ESAA vectors
     G_vec_ESAA = convert(Vector{Float64}, C_inv* g_vec_ESAA)
    dG_vec_ESAA = convert(Vector{Float64}, C_inv*dg_vec_ESAA)
    # return g(t), \dot g(t), GAST
    return G_vec_ESAA, dG_vec_ESAA, gast
end
