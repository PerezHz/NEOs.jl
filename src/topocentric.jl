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
    d = [k => identity.(cols[k]) for k in keys(colspecs)] # form Dict from data
    table((;d...)) # Dict -> NamedTuple -> IndexedTable
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
# et: ephemeris time (TDB seconds since J2000.0 epoch)
function t2c_rotation_iau_76_80(et::T; pm::Bool=true, lod::Bool=true,
        eocorr::Bool=true) where {T<:Number}
    # UTC (JD)
    utc_secs = et - tdb_utc(et)
    t_utc = J2000 + utc_secs/daysec
    t_utc_00 = constant_term(constant_term(t_utc))
    # TT
    t0_tt = et + ttmtdb(et)
    tt = t0_tt/daysec

    #Nutation corrections wrt IAU 1976/1980
    # Output of `EarthOrientation.precession_nutation80` is in mas:
    # https://github.com/JuliaAstro/EarthOrientation.jl/blob/529f12425a6331b133f989443aeb3fbbafd8f324/src/EarthOrientation.jl#L413
    if eocorr
        ddp80_mas, dde80_mas = EarthOrientation.precession_nutation80(t_utc_00)
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
        dut1 = EarthOrientation.getΔUT1(t_utc_00)
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
    β_dot = lod ? omega(EarthOrientation.getlod(t_utc_00)) : ω
    # if isa(gast, Taylor1)
    #     gast[1] = β_dot*one(gast[1])
    # end

    # For more details, see ESAA 2014, p. 295, Sec. 7.4.3.3, Eqs. 7.137-7.140
    Rz_minus_GAST = [
            cos(gast) -sin(gast) zero(gast);
            sin(gast) cos(gast) zero(gast);
            zero(gast) zero(gast) one(gast)
        ]
    dRz_minus_GAST = (β_dot)*[
            -sin(gast) -cos(gast) zero(gast);
            cos(gast) -sin(gast) zero(gast);
            zero(gast) zero(gast) zero(gast)
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
    # Velocity transformation may be retrieved also from: SatelliteToolbox.svECEFtoECI
    mt2c = C_inv*Rz_minus_GAST*W_inv
    dmt2c = C_inv*dRz_minus_GAST*W_inv

    return mt2c, dmt2c, gast
end
