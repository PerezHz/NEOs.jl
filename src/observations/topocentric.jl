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

function sv_ecef_to_ecef(
    sv::OrbitStateVector,
    T_ECEF1::Val{:ITRF},
    T_ECEF2::Val{:TIRS},
    jd_utc::Taylor1{TaylorN{Float64}},
    eop_data::Union{Nothing, EopIau1980, EopIau2000A} = nothing
)
    D = r_ecef_to_ecef(DCM, T_ECEF1, T_ECEF2, jd_utc, eop_data)

    # Since both frames does not have a significant angular velocity between
    # them, then we just need to convert the representations.
    r_ecef::Vector{Taylor1{TaylorN{Float64}}} = D * sv.r
    v_ecef::Vector{Taylor1{TaylorN{Float64}}} = D * sv.v
    a_ecef::Vector{Taylor1{TaylorN{Float64}}} = D * sv.a
    return OrbitStateVector(sv.t, r_ecef, v_ecef, a_ecef)
end

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

    return OrbitStateVector(sv.t, r_eci, v_eci, a_eci)
end

@doc raw"""
    obsposvelECI(observatory::ObservatoryMPC{T}, et::T;
               eop::Union{EopIau1980, EopIau2000A} = eop_IAU2000A) where {T <: AbstractFloat}
    obsposvelECI(x::RadecMPC{T}; eop::Union{EopIau1980, EopIau2000A} = eop_IAU2000A) where {T <: AbstractFloat}
    obsposvelECI(x::RadarJPL{T}; eop::Union{EopIau1980, EopIau2000A} = eop_IAU2000A) where {T <: AbstractFloat}

Return the observer's geocentric `[x, y, z, v_x, v_y, v_z]` "state" vector in Earth-Centered
Inertial (ECI) reference frame. By default, the IAU200A Earth orientation model is used to
transform from Earth-centered, Earth-fixed (ECEF) frame to ECI frame. Other Earth orientation
models, such as the IAU1976/80 model, can be used by importing the
`SatelliteToolboxTransformations.EopIau1980` type and passing it to the `eop` keyword
argument in the function call.

See also [`SatelliteToolboxBase.OrbitStateVector`](@ref) and [`SatelliteToolboxTransformations.sv_ecef_to_eci`](@ref).

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
    pv_ECEF = OrbitStateVector(jd_utc, pos_ECEF, zeros(3), zeros(3))

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
