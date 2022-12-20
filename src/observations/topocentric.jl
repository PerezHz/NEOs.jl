# Path to leap seconds kernel 
const naif_path = joinpath(observations_path, "naif0012.tls")

@doc raw"""
    load_leap_seconds()

Loads leap seconds kernel. 
"""
function load_leap_seconds()
    furnsh(naif_path)
end

@doc raw"""
    datetime2et(x::DateTime)
    datetime2et(x::RadecMPC{T}) where {T <: AbstractFloat}
    
Retuns the  TDB seconds past the J2000 epoch.
"""
function datetime2et(x::DateTime)
    return str2et(string(x))
end     

datetime2et(x::RadecMPC{T}) where {T <: AbstractFloat} = datetime2et(x.date)

@doc raw"""
    rad2arcsec(x)

Converts radians to arcseconds. 

See also [`arcsec2rad`](@ref) and [`mas2rad`](@ref).
"""
rad2arcsec(x) = 3600 * rad2deg(x) # rad2deg(rad) -> deg; 3600 * deg -> arcsec

@doc raw"""
    arcsec2rad(x)

Converts arcseconds to radians. 

See also [`rad2arcsec`](@ref) and [`mas2rad`](@ref).
"""
arcsec2rad(x) = deg2rad(x / 3600) # arcsec/3600 -> deg; deg2rad(deg) -> rad

@doc raw"""
    mas2rad(x)

Converts milli-arcseconds to radians. 

See also [`rad2arcsec`](@ref) and [`arcsec2rad`](@ref).
"""
mas2rad(x) = arcsec2rad(x / 1000) # mas/1000 -> arcsec; arcsec2rad(arcsec) -> rad

# Functions get_eop_iau1980, get_eop_iau2000a were adapted from SatelliteToolbox.jl; MIT-licensed
# these functions avoid the use of @eval
# https://github.com/JuliaSpace/SatelliteToolbox.jl/blob/b95e7f54f85c26744c64270841c874631f5addf1/src/transformations/eop.jl#L87

@doc raw"""
    get_eop_iau1980(; force=false)

Adaptation of [`SatelliteToolbox.get_iers_eop_iau_1980`](@ref) that avoids the use of `@eval`.
See https://github.com/JuliaSpace/SatelliteToolbox.jl/blob/master/src/transformations/eop.jl.
"""
function get_eop_iau1980(; force = false)
    # Remote file 
    eop_iau1980_rf = RemoteFiles.@RemoteFile(
        EOP_IAU1980_RF, 
        "https://datacenter.iers.org/data/csv/finals.all.csv", 
        file = "EOP_IAU1980.TXT",
        dir = observations_path, 
        updates = :weekly
    )
    # Download remote file 
    RemoteFiles.download(EOP_IAU1980_RF, force = force)
    # Read delimited file 
    eop, ~ = readdlm(RemoteFiles.path(eop_iau1980_rf), ';'; header = true)
    # Parse eop 
    SatelliteToolbox._parse_iers_eop_iau_1980(eop)
end

@doc raw"""
    get_eop_iau2000a(; force=false)

Adaptation of [`SatelliteToolbox.get_iers_eop_iau_2000A`](@ref) that avoids the use of `@eval`.
See https://github.com/JuliaSpace/SatelliteToolbox.jl/blob/master/src/transformations/eop.jl.
"""
function get_eop_iau2000a(; force = false)
    # Remote file 
    eop_iau2000a_rf = RemoteFiles.@RemoteFile(
        EOP_IAU2000A_RF, 
        "https://datacenter.iers.org/data/csv/finals2000A.all.csv", 
        file = "EOP_IAU2000A.TXT", 
        dir = observations_path, 
        updates = :weekly
    )
    # Download remote file 
    RemoteFiles.download(EOP_IAU2000A_RF, force = force)
    # Read delimited file 
    eop, ~ = readdlm(RemoteFiles.path(eop_iau2000a_rf), ';'; header = true)
    # Parse eop 
    SatelliteToolbox._parse_iers_eop_iau_2000A(eop)
end

# Earth orientation parameters (eop) 1980
const eop_IAU1980 = get_eop_iau1980()
# Earth orientation parameters (eop) 2000
const eop_IAU2000A = get_eop_iau2000a()

@doc raw"""
    observer_position(obs::RadecMPC{T}; eo::Bool=true, eop::Union{EOPData_IAU1980, EOPData_IAU2000A} = eop_IAU1980) where {T <: AbstractFloat}

Returns the observer's position and velocity in Earth-Centered Inertial (ECI) reference frame.

See also [`SatelliteToolbox.satsv`](@ref) and [`SatelliteToolbox.svECEFtoECI`](@ref).

# Arguments 

- `obs::RadecMPC{T}`: observation instant.
- `eo::Bool=true`: wheter to use Earth Orientation Parameters (eop) or not. 
- `eop::Union{EOPData_IAU1980, EOPData_IAU2000A}`: Earth Orientation Parameters (eop).
"""
function observer_position(obs::RadecMPC{T}; eo::Bool=true, eop::Union{EOPData_IAU1980, EOPData_IAU2000A} = eop_IAU1980) where {T <: AbstractFloat}

    observatory = obs.observatory

    @assert !isunknown(observatory) "Cannot compute position for unknown observatory."
    @assert hascoord(observatory) "Cannot compute position for observatory without coordinates"

    # λ_deg: longitude [degrees east of Greenwich]
    # u: distance from spin axis [km], u = ρ*cos(ϕ')
    # v: height above equatorial plane [km], v = ρ*sin(ϕ'),
    #  where ϕ' is the geocentric latitude and ρ is the geocentric distance in km

    # Cilindrical components of Earth-fixed position of observer
    λ_deg = observatory.long
    u = observatory.cos * RE
    v = observatory.sin * RE
    # Cartesian components of Earth-fixed position of observer
    λ_rad = deg2rad(λ_deg)       # rad
    x_gc = u * cos(λ_rad)        # km
    y_gc = u * sin(λ_rad)        # km
    z_gc = v                     # km
    # Earth-fixed position of observer 
    pos_geo = [x_gc, y_gc, z_gc] # km

    # UTC seconds 
    et = datetime2et(obs)
    utc_secs = et - tdb_utc(et)
    # Julian days UTC 
    jd_utc = JD_J2000 + utc_secs/daysec
    # State vector 
    sv_geo = SatelliteToolbox.satsv(jd_utc, pos_geo, zeros(3), zeros(3))
    # Transform position/velocity from Earth-Centered, Earth-fixed frame to Earth-Centered Inertial frame
    # ITRF: International Terrestrial Reference Frame
    # GCRF: Geocentric Celestial Reference Frame
    # Use earth orientation parameters
    if eo
        sv_c = SatelliteToolbox.svECEFtoECI(sv_geo, Val(:ITRF), Val(:GCRF), jd_utc, eop)
    # Not use earth orientation parameters        
    else
        sv_c = SatelliteToolbox.svECEFtoECI(sv_geo, Val(:ITRF), Val(:GCRF), jd_utc)
    end
    # Inertial position 
    r_c = convert(Vector{eltype(sv_c.r)}, sv_c.r)
    # Inertial velocity
    v_c = convert(Vector{eltype(sv_c.v)}, sv_c.v)

    return r_c, v_c
end

@doc raw"""
    nupr7680mat(tt, Δϵ_1980, ΔΨ_1980, dde80, ddp80) 

Returns IAU 1976/1980 nutation-precession matrix. 

# Arguments 

- `tt`: days since J2000.0 (TT). 
- `Δϵ_1980, ΔΨ_1980`: IAU 1980 nutation angles Δψ (nutation in longitude), Δϵ (nutation in obliquity), both in radians.
- `dde80, ddp80`: nutation corrections wrt IAU 1976/1980 (in radians).
"""
function nupr7680mat(tt, Δϵ_1980, ΔΨ_1980, dde80, ddp80)
    # IAU 1976 precession matrix, J2000.0 to date
    rp = Rz(-PlanetaryEphemeris.zeta(tt))*Ry(PlanetaryEphemeris.Theta(tt))*Rz(-PlanetaryEphemeris.Zeta(tt))
    # Add nutation corrections
    dpsi = ΔΨ_1980 + ddp80   # rad
    deps = Δϵ_1980 + dde80   # rad
    # IAU 1980 nutation matrix 
    # See Explanatory Supplement to the Astronomical Almanac 1992
    # See equation (5-152) in page (5-60) of https://doi.org/10.1002/0471728470
    # ϵ0 (rad): mean obliquity of the ecliptic
    ϵ0 = PlanetaryEphemeris.ϵ̄(tt)
    # Δϵ (rad): nutation in obliquity
    Δϵ = deps
    # Δψ (rad): nutation in longitude
    Δψ = dpsi
    # ϵ (rad): true obliquity of date
    ϵ = ϵ0 + Δϵ
    # Nutation matrix 
    rn = Rx(-ϵ)*Rz(-Δψ)*Rx(ϵ0)
    # Combine the matrices:  PN = N x P (with frame-bias B included)
    return rn*rp
end

@doc raw"""
    polarmotionmat(tt; eo::Bool=true)

Returns the polar motion matrix in radians (TIRS->ITRS, IERS 1996).

See also [`EarthOrientation.polarmotion`](@ref). 

# Arguments 

- `tt`: days since J2000.0.
- `eo::Bool`: wheter to use Earth orientation corrections from IERS or not.
"""
function polarmotionmat(tt; eo::Bool=true)
    # Polar motion (arcsec->radians)
    if eo
        xp_arcsec, yp_arcsec = EarthOrientation.polarmotion(JD_J2000+tt)
        xp = deg2rad(xp_arcsec/3600)
        yp = deg2rad(yp_arcsec/3600)
    else
        xp = 0.0
        yp = 0.0
    end
    return Ry(-xp)*Rx(-yp)
end

@doc raw"""
    t2c_rotation_iau_76_80(et::T; eo::Bool=true) where {T<:Number}

Returns the terrestrial-to-celestial rotation matrix (including polar motion) 
using 1976/1980 Earth orientation/rotation model.

See also [`EarthOrientation.precession_nutation80`](@ref), [`SatelliteToolbox.nutation_fk5`](@ref),
[`EarthOrientation.getΔUT1`](@ref), [`SatelliteToolbox.J2000toGMST`](@ref),
[`EarthOrientation.getlod`](@ref), [`polarmotionmat`](@ref) and [`nupr7680mat`](@ref). 

# Arguments

- `et`: ephemeris time (TDB seconds since J2000.0 epoch).
- `eo::Bool`: wheter to use Earth orientation corrections from IERS or not.
"""
function t2c_rotation_iau_76_80(et::T; eo::Bool=true) where {T<:Number}
    # UTC (JD)
    utc_secs = et - tdb_utc(et)                        # seconds
    t_utc = JD_J2000 + utc_secs/daysec                 # days 
    t_utc_00 = constant_term(constant_term(t_utc))     # days 
    # TT
    t0_tt = et + ttmtdb(et)                            # seconds
    tt = t0_tt/daysec                                  # days 

    # Nutation corrections wrt IAU 1976/1980
    # Output of `EarthOrientation.precession_nutation80` is in milli-arcseconds:
    # https://github.com/JuliaAstro/EarthOrientation.jl/blob/529f12425a6331b133f989443aeb3fbbafd8f324/src/EarthOrientation.jl#L413
    
    # Use Earth orientation corrections from IERS
    if eo
        ddp80_mas, dde80_mas = EarthOrientation.precession_nutation80(t_utc_00)
    # Not use Earth orientation corrections
    else
        ddp80_mas, dde80_mas = 0.0, 0.0
    end
    # Convert milli-arcseconds -> radians
    ddp80 = mas2rad(ddp80_mas)
    dde80 = mas2rad(dde80_mas)
    # TDB days since J2000.0
    et_days = et/daysec 
    # Mean obliquity (output in radians)
    epsa = PlanetaryEphemeris.ϵ̄(et_days) # rad
    # Lunar longitude of ascending node (measured from ecliptic)
    Ω_M = PlanetaryEphemeris.Ω(et_days)  # rad
    # IAU 1980 nutation angles
    _, Δϵ_1980, ΔΨ_1980 = nutation_fk5(JD_J2000 + et_days)
    # Equation of the equinoxes `ee = GAST - GMST`, including nutation correction
    # See equation (5-184) in page (5-70) of https://doi.org/10.1002/0471728470:
    # Δθ = ∆ψ cosϵ̄ + 0''.00264 sinΩ + 0''.000063 sin2Ω
    # where the nutation in longitude includes corrections from IERS eop file
    ee = (ΔΨ_1980+ddp80)*cos(epsa) + arcsec2rad( 0.00264sin(Ω_M) + 0.000063sin(2Ω_M) )
    # ΔUT1 = UT1-UTC (seconds)

    # Use Earth orientation corrections from IERS
    if eo
        dut1 = EarthOrientation.getΔUT1(t_utc_00)
        # dut1 = eop_IAU1980.UT1_UTC(t_utc)
    # Not use Earth orientation corrections
    else
        dut1 = 0.0
    end

    # UT1
    ut1 = utc_secs + dut1 # elapsed UT1 secs since J2000.0
    ut1_days = ut1/daysec # elapsed UT1 days since J2000.0
    # Greenwich apparent sidereal time (IAU 1982/1994)
    # See equation (5-173) in page (5-67) of https://doi.org/10.1002/0471728470:
    gmst82 = J2000toGMST(ut1_days)
    gast = mod2pi( gmst82 + ee )
    β_dot = eo ? omega(EarthOrientation.getlod(t_utc_00)) : ω

    # For more details, see Explanatory Supplement to the Astronomical Almanac 2014, 
    # p. 295, Sec. 7.4.3.3, Eqs. 7.137-7.140
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
    W = polarmotionmat(constant_term(constant_term(tt)), eo=eo)
    W_inv = convert(Matrix{Float64}, transpose(W))
    # _W_inv = SatelliteToolbox.rECEFtoECEF(DCM, ITRF(), PEF(), t_utc, eop_IAU1980)
    # W_inv = convert(Matrix{eltype(_W_inv)}, _W_inv)

    # IAU 76/80 nutation-precession matrix
    C = nupr7680mat(et_days, Δϵ_1980, ΔΨ_1980, dde80, ddp80)
    C_inv = transpose(C)
    # _C_inv = SatelliteToolbox.rECItoECI(DCM, TOD(), GCRF(), t_utc, eop_IAU1980)
    # C_inv = convert(Matrix{eltype(_C_inv)}, _C_inv)

    # Velocity transformation may be retrieved also from: SatelliteToolbox.svECEFtoECI
    mt2c = C_inv*Rz_minus_GAST*W_inv
    dmt2c = C_inv*dRz_minus_GAST*W_inv

    return mt2c, dmt2c, gast
end
