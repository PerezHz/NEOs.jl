### FWF reader, due to @aplavin
### https://gist.github.com/aplavin/224a31ea457b6e0ef0f4c1a20bd28850
convert_val(::Type{String}, val::String) = val
convert_val(::Type{Symbol}, val::String) = Symbol(val)
convert_val(typ::Type{<:Integer}, val::String) = parse(typ, val)
convert_val(typ::Type{<:AbstractFloat}, val::String) = parse(typ, replace(val, "D" => "E"))  # tables output from Fortran often have floats as 1D+5 instead of 1E+5
# # usage:
# read_fwf(
#   "file.tab",
#   (
#     colname1=(1, 20, String),
#     colname2=(100, 120, Int),
#     # ...
#   ),
#   skiprows=[1, 2, 3, 7]
# )
function readfwf(io, colspecs; skiprows=[], missingstrings=[])
    cols = Dict(
        k => Vector{Union{typ, Missing}}(undef, 0)
        for (k, (from, to, typ)) in pairs(colspecs)
    )
    for (irow, line) in eachline(io) |> enumerate
        if irow ∈ skiprows continue end
        for (k, (from, to, typ)) in pairs(colspecs)
            s_val = from <= length(line) ? line[from:min(length(line), to)] : ""
            # @show irow line k typ s_val
            f_val = s_val in missingstrings ? missing : convert_val(typ, s_val)
            push!(cols[k], f_val)
        end
    end
    DataFrame([k => identity.(cols[k]) for k in keys(colspecs)])
end

mpc_format_obscode = (Code=(1,3,String),
Long=(5,13,Float64),
cos=(14,21,Float64),
sin=(22,30,Float64),
Name=(31,80,String)
)

# MPC minor planet observatory code reader
readmpcobs(mpcfile::String=joinpath(dirname(pathof(Apophis)), "ObsCodes.txt")) = table(readfwf(mpcfile, mpc_format_obscode, skiprows=union([1,247,249,251,252,260],1222:1230)))

const mpcobscodes = readmpcobs()

# et: ephemeris time (TDB seconds since J2000.0 epoch)
function observer_position(station_code::Union{Int, String}, et::T;
        pm::Bool=true, lod::Bool=true, eocorr::Bool=true) where {T<:Number}
    # λ_deg: East longitude (deg)
    # u: distance from spin axis (km), taken from Yeomans et al. (1992) u = r*cos(ϕ)
    # v: height above equatorial plane (km), taken from Yeomans et al. (1992) v = r*sin(ϕ)
    # East long data is consistent with MPC database (2019-Mar-7)
    st_code_str = lpad(string(station_code), 3, "0")
    # @show st_code_str filter(i->i.Code==st_code_str, mpcobscodes)
    obscode = filter(i->i.Code==st_code_str, mpcobscodes)
    if length(obscode) == 0
        @error "Unknown station code: $station_code"
    elseif length(obscode) > 1
        @warn "More than one observatory assigned to code: $station_code"
    end
    λ_deg = select(obscode, :Long)[1]
    u = select(obscode, :cos)[1]*RE
    v = select(obscode, :sin)[1]*RE
    # cartesian components of Earth-fixed position of observer
    λ_rad = deg2rad(λ_deg) # rad
    x_gc = u*cos(λ_rad) #km
    y_gc = u*sin(λ_rad) #km
    z_gc = v #km

    pos_geo = [x_gc, y_gc, z_gc] #km

    # G_vec_ESAA, dG_vec_ESAA, era = t2c_rotation_iau_00_06(et, pos_geo, pm=pm)
    mt2c, dmt2c, gast = t2c_rotation_iau_76_80(et, pm=pm, lod=lod, eocorr=eocorr)

    # Apply rotation from geocentric, Earth-fixed frame to inertial (celestial) frame
    return mt2c*pos_geo, dmt2c*pos_geo
end

# conversion of radians to arcseconds
rad2arcsec(x) = 3600rad2deg(x) # rad2deg(rad) -> deg; 3600deg -> arcsec

# conversion of arcseconds to radians
arcsec2rad(x) = deg2rad(x/3600) # arcsec/3600 -> deg; deg2rad(deg) -> rad

# conversion of milli-arcseconds to radians
mas2rad(x) = arcsec2rad(x/1000) # mas/1000 -> arcsec; arcsec2rad(arcsec) -> rad

# Terrestrial-to-celestial rotation matrix (including polar motion)
# Reproduction of Section 5.3 of SOFA Tools for Earth Attitude
# "IAU 2000A, CIO based, using classical angles"
# found at SOFA website, Mar 27, 2019
# Some modifications were applied, using TaylorSeries.jl, in order to compute
# the geocentric velocity of the observer, following the guidelines from
# ESAA 2014, Sec 7.4.3.3 (page 295)
# et: ephemeris time (TDB seconds since J2000.0 epoch)
function t2c_rotation_iau_00_06(et::Float64, pos_geo::Vector; pm::Bool=true)
    # UTC
    t_utc = DateTime(et2utc(constant_term(et), "ISOC", 3))
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

# IAU 1976/1980 nutation-precession matrix
# tt: days since J2000.0 (TT)
# IAU 1980 nutation angles Δψ (nutation in longitude), Δϵ (nutation in obliquity), both in radians
# dde80, ddp80: nutation corrections wrt IAU 1976/1980 (in radians)
# eocorr: EarthOrientation.jl corrections
function nupr7680mat(tt, Δϵ_1980, ΔΨ_1980, dde80, ddp80; eocorr::Bool=true)
    # IAU 1976 precession matrix, J2000.0 to date
    rp = Rz(-PlanetaryEphemeris.zeta(tt))*Ry(PlanetaryEphemeris.Theta(tt))*Rz(-PlanetaryEphemeris.Zeta(tt))
    # Add nutation corrections
    dpsi = ΔΨ_1980 + ddp80 #rad
    deps = Δϵ_1980 + dde80 #rad
    # IAU 1980 nutation matrix (ESAA 1992)
    # Eq. (5-152) from Moyer, 2003
    # ϵ0 (rad): mean obliquity obliquity of the ecliptic
    ϵ0 = PlanetaryEphemeris.ϵ̄(tt)
    # Δϵ (rad): nutation in obliquity
    Δϵ = deps
    # Δψ (rad): nutation in longitude
    Δψ = dpsi
    # ϵ (rad): true obliquity of date
    ϵ = ϵ0 + Δϵ
    rn = Rx(-ϵ)*Rz(-Δψ)*Rx(ϵ0)
    # Combine the matrices:  PN = N x P (with frame-bias B included)
    return rn*rp
end

# Polar motion matrix (TIRS->ITRS, IERS 1996)
# tt: days since J2000.0
function polarmotionmat(tt; pm::Bool=true)
    # Polar motion (arcsec->radians)
    if pm
        xp_arcsec, yp_arcsec = EarthOrientation.polarmotion(J2000+tt)
        xp = deg2rad(xp_arcsec/3600)
        yp = deg2rad(yp_arcsec/3600)
    else
        xp = 0.0
        yp = 0.0
    end
    return Ry(-xp)*Rx(-yp)
end

# Terrestrial-to-celestial rotation matrix (including polar motion)
# Using 1976/1980 Earth orientation/rotation model
# Reproduction of Section 5.2 of SOFA Tools for Earth Attitude
# "IAU 1976/1980/1982/1994, equinox based"
# found at SOFA website, Mar 27, 2019
# et: ephemeris time (TDB seconds since J2000.0 epoch)
function t2c_rotation_iau_76_80(et::T; pm::Bool=true, lod::Bool=true,
        eocorr::Bool=true) where {T<:Number}
    # UTC (JD)
    utc_secs = et - tdb_utc(et)
    t_utc = J2000 + utc_secs/daysec
    # TT
    jd0 = datetime2julian(DateTime(2008,9,24))
    t0_tt = et + ttmtdb(et)
    tt = t0_tt/daysec

    #Nutation corrections wrt IAU 1976/1980
    # Output of `EarthOrientation.precession_nutation80` is in mas:
    # https://github.com/JuliaAstro/EarthOrientation.jl/blob/529f12425a6331b133f989443aeb3fbbafd8f324/src/EarthOrientation.jl#L413
    if eocorr
        ddp80_mas, dde80_mas = EarthOrientation.precession_nutation80(constant_term(t_utc))
    else
        ddp80_mas, dde80_mas = 0.0, 0.0
    end
    # Convert mas -> radians
    ddp80 = mas2rad(ddp80_mas)
    dde80 = mas2rad(dde80_mas)

    et_days = et/daysec # TDB days since J2000.0
    # Mean obliquity (output in radians)
    epsa = PlanetaryEphemeris.ϵ̄(et_days) # rad
    # Lunar longitude of ascending node (measured from ecliptic)
    Ω_M = PlanetaryEphemeris.Ω(et_days) # rad
    # IAU 1980 nutation angles
    _, Δϵ_1980, ΔΨ_1980 = nutation_fk5(J2000 + et_days)
    # Equation of the equinoxes `ee = GAST - GMST`, including nutation correction
    # Expression from Eq. (5-184) of Moyer (2003): Δθ = ∆ψ cosϵ̄ + 0''.00264 sinΩ + 0''.000063 sin2Ω
    # where the nutation in longitude includes corrections from IERS eop file
    ee = (ΔΨ_1980+ddp80)*cos(epsa) + arcsec2rad( 0.00264sin(Ω_M) + 0.000063sin(2Ω_M) )
    # ΔUT1 = UT1-UTC (seconds)
    if eocorr
        dut1 = EarthOrientation.getΔUT1(constant_term(t_utc))
    else
        dut1 = 0.0
    end
    # UT1
    ut1 = utc_secs + dut1 # elapsed UT1 secs since J2000.0
    ut1_days = ut1/daysec # elapsed UT1 days since J2000.0
    # Greenwich apparent sidereal time (IAU 1982/1994)
    # Eq. (5-173) from Moyer (2003)
    gmst82 = J2000toGMST(ut1_days)
    gast = mod2pi( gmst82 + ee )

    # For more details, see ESAA 2014, p. 295, Sec. 7.4.3.3, Eqs. 7.137-7.140
    Rz_minus_GAST = [
            cos(gast) -sin(gast) 0.0;
            sin(gast) cos(gast) 0.0;
            0.0 0.0 1.0
        ]
    if lod
        β_dot = omega(EarthOrientation.getlod(constant_term(t_utc)))
    else
        β_dot = ω
    end
    dRz_minus_GAST = (β_dot)*[
            -sin(gast) -cos(gast) 0.0;
            cos(gast) -sin(gast) 0.0;
            0.0 0.0 0.0
        ]

    # Polar motion matrix (TIRS->ITRS, IERS 1996)
    W = polarmotionmat(constant_term(tt), pm=pm)

    # IAU 76/80 nutation-precession matrix
    C = nupr7680mat(et_days, Δϵ_1980, ΔΨ_1980, dde80, ddp80, eocorr=eocorr)

    W_inv = convert(Matrix{Float64}, transpose(W))
    C_inv = transpose(C)

    # eop_IAU1980 = get_iers_eop();
    # _mt2c = rECEFtoECI(DCM, ITRF(), GCRF(), t_utc, eop_IAU1980)
    # mt2c = convert(Matrix{eltype(_mt2c)}, _mt2c)
    mt2c = C_inv*Rz_minus_GAST*W_inv
    dmt2c = C_inv*dRz_minus_GAST*W_inv

    return mt2c, dmt2c, gast
end
