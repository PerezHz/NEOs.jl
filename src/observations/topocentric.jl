# Earth orientation parameters (eop) 1980
# const eop_IAU1980 = ST.get_iers_eop_iau_1980()
# Earth orientation parameters (eop) 2000
const eop_IAU2000A::EopIau2000A = fetch_iers_eop(Val(:IAU2000A))

@doc raw"""
    obs_pos_ECEF(observatory::ObservatoryMPC{T}) where {T <: AbstractFloat}
    obs_pos_ECEF(x::RadecMPC{T}) where {T <: AbstractFloat}
    obs_pos_ECEF(x::RadarJPL{T}) where {T <: AbstractFloat}

Return the observer's geocentric `[x, y, z]` position vector in Earth-Centered Earth-Fixed (ECEF) reference frame.
"""
function obs_pos_ECEF(observatory::ObservatoryMPC{T}) where {T <: AbstractFloat}

    # Make sure observatory has coordinates
    @assert !isunknown(observatory) "Cannot compute position for unknown observatory."
    @assert hascoord(observatory) "Cannot compute position for observatory [$(observatory.code)] without coordinates"

    # λ_deg: longitude [degrees east of Greenwich]
    # u: distance from spin axis [km], u = ρ*cos(ϕ')
    # v: height above equatorial plane [km], v = ρ*sin(ϕ'),
    # where ϕ' is the geocentric latitude and ρ is the geocentric distance in km

    # Cilindrical components of Earth-Centered Earth-Fixed position of observer
    λ_deg = observatory.long     # deg
    u = observatory.cos * RE     # km
    v = observatory.sin * RE     # km

    # Cartesian components of Earth-Centered Earth-Fixed position of observer
    λ_rad = deg2rad(λ_deg)       # rad
    x_gc = u * cos(λ_rad)        # km
    y_gc = u * sin(λ_rad)        # km
    z_gc = v                     # km

    # Earth-Centered Earth-Fixed position position of observer
    pos_ECEF = [x_gc, y_gc, z_gc] # km

    return pos_ECEF
end

obs_pos_ECEF(x::RadecMPC{T}) where {T <: AbstractFloat} = obs_pos_ECEF(x.observatory)
obs_pos_ECEF(x::RadarJPL{T}) where {T <: AbstractFloat} = obs_pos_ECEF(x.rcvr)

function sv_ecef_to_eci(
    sv::OrbitStateVector,
    T_ECEF::Union{Val{:PEF}, Val{:TIRS}},
    T_ECI::Union{T_ECIs, T_ECIs_IAU_2006},
    jd_utc::Taylor1{TaylorN{Float64}},
    eop_data::Union{Nothing, EopIau1980, EopIau2000A} = nothing
)
    # Get the matrix that converts the ECEF to the ECI.
    if eop_data === nothing
        D = r_ecef_to_eci(DCM, T_ECEF, T_ECI, jd_utc)
    else
        D = r_ecef_to_eci(DCM, T_ECEF, T_ECI, jd_utc, eop_data)
    end

    # Since the ECI and ECEF frames have a relative velocity between them, then
    # we must account from it when converting the velocity and acceleration. The
    # angular velocity between those frames is computed using `we` and corrected
    # by the length of day (LOD) parameter of the EOP data, if available.
    ω  = EARTH_ANGULAR_SPEED * (1 - (eop_data !== nothing ? eop_data.lod(jd_utc) / 86400000 : 0))
    vω = [0, 0, ω]

    # Compute the position in the ECI frame.
    r_eci::Vector{Taylor1{TaylorN{Float64}}} = D * sv.r

    # Compute the velocity in the ECI frame.
    vω_x_r = vω × sv.r
    v_eci::Vector{Taylor1{TaylorN{Float64}}} = D * (sv.v + vω_x_r )

    # Compute the acceleration in the ECI frame.
    a_eci::Vector{Taylor1{TaylorN{Float64}}} = D * (sv.a + vω × vω_x_r + 2vω × sv.v)

    return orbsv(sv.t, r_eci, v_eci, a_eci)
end

@doc raw"""
    obsposvelECI(observatory::ObservatoryMPC{T}, et::T;
               eop::Union{EopIau1980, EopIau2000A} = eop_IAU2000A) where {T <: AbstractFloat}
    obsposvelECI(x::RadecMPC{T}; eop::Union{EopIau1980, EopIau2000A} = eop_IAU2000A) where {T <: AbstractFloat}
    obsposvelECI(x::RadarJPL{T}; eop::Union{EopIau1980, EopIau2000A} = eop_IAU2000A) where {T <: AbstractFloat}

Return the observer's geocentric `[x, y, z, v_x, v_y, v_z]` "state" vector in Earth-Centered
Inertial (ECI) reference frame. By default, the IAU200A Earth orientation model is used to
transform from Earth-centered, Earth-fixed (ECEF) frame to ECI frame. Other Earth orientation
models, such as the IAU1976/80 model can be used by importing the `SatelliteToolbox.EopIau1980`
type and passing it to the `eop` keyword argument in the function call.

See also [`SatelliteToolbox.orbsv`](@ref) and [`SatelliteToolbox.sv_ecef_to_eci`](@ref).

# Arguments

- `observatory::ObservatoryMPC{T}`: observation site.
- `et::T`: ephemeris time (TDB seconds since J2000.0 epoch).
- `eop::Union{EopIau1980, EopIau2000A}`: Earth Orientation Parameters (eop).
"""
function obsposvelECI(observatory::ObservatoryMPC{T}, et::ET;
        eop::Union{EopIau1980, EopIau2000A} = eop_IAU2000A) where {T <: AbstractFloat, ET<:Union{T,Taylor1{T},Taylor1{TaylorN{T}}}}
    # Earth-Centered Earth-Fixed position position of observer
    pos_ECEF = obs_pos_ECEF(observatory)

    # UTC seconds
    utc_secs = et - tdb_utc(et)
    # Julian days UTC
    jd_utc = JD_J2000 + utc_secs/daysec
    # State vector
    pv_ECEF = orbsv(jd_utc, pos_ECEF, zeros(3), zeros(3))

    # Transform position/velocity from Earth-Centered Earth-fixed (ECEF) frame to Earth-Centered Inertial (ECI) frame
    # ITRF: International Terrestrial Reference Frame
    # GCRF: Geocentric Celestial Reference Frame
    pv_ECI = sv_ecef_to_eci(pv_ECEF, Val(:ITRF), Val(:GCRF), jd_utc, eop)

    # Inertial position
    p_ECI = convert(Vector{eltype(pv_ECI.r)}, pv_ECI.r)
    # Inertial velocity
    v_ECI = convert(Vector{eltype(pv_ECI.v)}, pv_ECI.v)

    # Concat position and velocity
    return vcat(p_ECI, v_ECI)
end

function obsposvelECI(x::RadecMPC{T}; eop::Union{EopIau1980, EopIau2000A} = eop_IAU2000A) where {T <: AbstractFloat}
    return obsposvelECI(x.observatory, datetime2et(x.date); eop)
end
function obsposvelECI(x::RadarJPL{T}; eop::Union{EopIau1980, EopIau2000A} = eop_IAU2000A) where {T <: AbstractFloat}
    return obsposvelECI(x.rcvr, datetime2et(x.date); eop)
end

@doc raw"""
    nupr7680mat(tt, Δϵ_1980, ΔΨ_1980, dde80, ddp80)

Return IAU 1976/1980 nutation-precession matrix.

# Arguments

- `tt`: days since J2000.0 (TT).
- `Δϵ_1980, ΔΨ_1980`: IAU 1980 nutation angles Δψ (nutation in longitude), Δϵ (nutation in obliquity), both in radians.
- `dde80, ddp80`: nutation corrections wrt IAU 1976/1980 (in radians).
"""
function nupr7680mat(tt, Δϵ_1980, ΔΨ_1980, dde80, ddp80)
    # IAU 1976 precession matrix, J2000.0 to date
    rp = Rz(-PE.zeta(tt))*Ry(PE.Theta(tt))*Rz(-PE.Zeta(tt))
    # Add nutation corrections
    dpsi = ΔΨ_1980 + ddp80   # rad
    deps = Δϵ_1980 + dde80   # rad
    # IAU 1980 nutation matrix
    # See Explanatory Supplement to the Astronomical Almanac 1992
    # See equation (5-152) in page (5-60) of https://doi.org/10.1002/0471728470
    # ϵ0 (rad): mean obliquity of the ecliptic
    ϵ0 = PE.ϵ̄(tt)
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
    polarmotionmat(tt; eo::Bool = true)

Return the polar motion matrix in radians (TIRS->ITRS, IERS 1996).

See also [`EarthOrientation.polarmotion`](@ref).

# Arguments

- `tt`: days since J2000.0.
- `eo::Bool`: wheter to use Earth orientation corrections from IERS or not.
"""
function polarmotionmat(tt; eo::Bool = true)
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
    omega(lod)

Return the angular velocity of the earth in units of rad/sec
```math
\omega = (72921151.467064 - 0.843994809\text{LOD})\times 10^{-12},
```
where LOD is the length of the day in milliseconds.

See https://www.iers.org/IERS/EN/Science/EarthRotation/UT1LOD.html.
"""
omega(lod) = (1e-12)*(72921151.467064 - 0.843994809lod)

@doc raw"""
    t2c_rotation_iau_76_80(et::T; eo::Bool = true) where {T <: Number}

Return the terrestrial-to-celestial rotation matrix (including polar motion)
using 1976/1980 Earth orientation/rotation model.

See also [`EarthOrientation.precession_nutation80`](@ref), [`SatelliteToolbox.nutation_fk5`](@ref),
[`EarthOrientation.getΔUT1`](@ref), [`SatelliteToolbox.J2000toGMST`](@ref),
[`EarthOrientation.getlod`](@ref), [`polarmotionmat`](@ref) and [`nupr7680mat`](@ref).

# Arguments

- `et`: ephemeris time (TDB seconds since J2000.0 epoch).
- `eo::Bool`: wheter to use Earth orientation corrections from IERS or not.
"""
function t2c_rotation_iau_76_80(et::T; eo::Bool = true) where {T <: Number}
    # UTC (JD)
    utc_secs = et - tdb_utc(et)                        # seconds
    t_utc = JD_J2000 + utc_secs/daysec                 # days
    t_utc_00 = constant_term(constant_term(t_utc))     # days
    # TT
    t0_tt = et + ttmtdb(et/daysec)                            # seconds
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
    epsa = PE.ϵ̄(et_days) # rad
    # Lunar longitude of ascending node (measured from ecliptic)
    Ω_M = PE.Ω(et_days)  # rad
    # IAU 1980 nutation angles
    _, Δϵ_1980, ΔΨ_1980 = ST.nutation_fk5(JD_J2000 + et_days)
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
    gmst82 = ST.J2000toGMST(ut1_days)
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
