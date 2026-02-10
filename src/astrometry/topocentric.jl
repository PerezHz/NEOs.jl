"""
    TimeOfDay <: AbstractAstrometryTime

Day/night at a particular timezone.

# Fields

- `light::Symbol`:
    - for geocentric observatories: `:geocentric`,
    - for space-based observatories: `:satellite`,
    - for ground-based observatories: `:day` or `:night`.
- `start/stop::DateTime`: beginning/end of the `:day` or `:night`.
- `utc::Int`: time difference from UTC [hours].
"""
@auto_hash_equals struct TimeOfDay <: AbstractAstrometryTime
    light::Symbol
    start::DateTime
    stop::DateTime
    utc::Int
end

isday(x::TimeOfDay) = x.light == :day
isnight(x::TimeOfDay) = x.light == :night
issatellite(x::TimeOfDay) = x.light == :satellite
isgeocentric(x::TimeOfDay) = x.light == :geocentric

# Print method for TimeOfDay
function show(io::IO, m::TimeOfDay)
    print(io, uppercasefirst(string(m.light)), " from ", Date(m.start), " to ",
        Date(m.stop), " at UTC", @sprintf("%+d", m.utc))
end

# Constructors
TimeOfDay(x::AbstractOpticalAstrometry; eop::EOPIAU = EOP_IAU2000A) =
    TimeOfDay(observatory(x), date(x); eop)

function TimeOfDay(observatory::ObservatoryMPC, date::DateTime; eop::EOPIAU = EOP_IAU2000A)
    if isgeocentric(observatory)
        return TimeOfDay(:geocentric, Date(date), Date(date), 0)
    elseif issatellite(observatory) || isoccultation(observatory)
        return TimeOfDay(:satellite, date, date, 0)
    end
    # Hours from UTC
    utc = hours_from_UTC(observatory, dtutc2et(date); eop)
    # Today's sunrise / sunset
    today = sunriseset(observatory, date; eop)
    # Yesterday's sunrise / sunset
    yesterday = sunriseset(observatory, date - Day(1); eop)
    # Tomorrow's sunrise / sunset
    tomorrow = sunriseset(observatory, date + Day(1); eop)
    # Selection
    if yesterday[2] <= date <= today[1]
        return TimeOfDay(:night, Date(yesterday[2]), Date(today[1]), utc)
    elseif today[1] <= date <= today[2]
        return TimeOfDay(:day, Date(today[1]), Date(today[2]), utc)
    elseif today[2] <= date <= tomorrow[1]
        return TimeOfDay(:night, Date(today[2]), Date(tomorrow[1]), utc)
    else
        throw(ArgumentError("Cannot find the sunrise and sunset for date"))
    end
end

# Return the naive hour difference between longitude `lon` [rad] and UTC.
hours_from_UTC(lon::Real) = ceil(Int, 12*lon/π - 0.5)

function hours_from_UTC(observatory::ObservatoryMPC, et::Number;
                        eop::EOPIAU = EOP_IAU2000A)
    lon, _ = lonlat(observatory, et; eop)
    return hours_from_UTC(lon)
end

# Return the geocentric longitude and latitude [rad] of an observatory.
function lonlat(observatory::ObservatoryMPC, et::Number; eop::EOPIAU = EOP_IAU2000A)
    # Observer's position in ITRF (ECEF) frame [m]
    posECEF = obsposECEF(observatory, et; eop) * 1_000
    # Transform from ITRF (ECEF) [m] to Geocentric [m]
    lat, lon, _ = ecef_to_geocentric(posECEF)

    return lon, lat
end

# Return the `DateTime` of sunrise and sunset at a particular date and location.
# See "General Solar Position Calculations" by NOAA at:
# - https://gml.noaa.gov/grad/solcalc/solareqns.PDF.
sunriseset(x::AbstractOpticalAstrometry; eop::EOPIAU = EOP_IAU2000A) =
    sunriseset(observatory(x), date(x); eop)

function sunriseset(observatory::ObservatoryMPC, date::DateTime;
                    eop::EOPIAU = EOP_IAU2000A)
    # Fractional year [rad]
    γ = 2π * (dayofyear(date) - 1 + (hour(date)-12)/24) / daysinyear(date)
    # Equation of time [min]
    eqtime = 229.18 * (0.000075 + 0.001868*cos(γ) - 0.032077*sin(γ) - 0.014615*cos(2γ)
             - 0.040849*sin(2γ) )
    # Solar declination angle [rad]
    δ = 0.006918 - 0.399912*cos(γ) + 0.070257*sin(γ) - 0.006758*cos(2γ) + 0.000907*sin(2γ)
        - 0.002697*cos(3γ) + 0.00148*sin(3γ)
    # Longitude and latitude [rad]
    lon, lat = lonlat(observatory, dtutc2et(date); eop)
    # Solar hour angle [deg]
    ha_sunrise = acosd(cosd(90.833)/cos(lat)/cos(δ) - tan(lat)*tan(δ))
    ha_sunset = -ha_sunrise
    # Day [DateTime]
    day = DateTime(Date(date))
    # UTC time of sunrise [min]
    sunrise_min = 720 - 4*(rad2deg(lon) + ha_sunrise) - eqtime
    # UTC time of sunrise [DateTime]
    sunrise = day + Microsecond(round(Int, sunrise_min*6e7))
    # UTC time of sunset [min]
    sunset_min = 720 - 4*(rad2deg(lon) + ha_sunset) - eqtime
    # UTC time of sunset [DateTime]
    sunset = day + Microsecond(round(Int, sunset_min*6e7))

    return sunrise, sunset
end

"""
    obsposECEF(::AbstractAstrometryObservation; kwargs...)

Return the observer's geocentric cartesian position vector [km] in
Earth-Centered Earth-Fixed (ECEF) reference frame.

# Keyword argument

- `eop::Union{EopIau1980, EopIau2000A}`: Earth Orientation Parameters
    (default: `EOP_IAU2000A`).
"""
obsposECEF(x::AbstractAstrometryObservation; eop::EOPIAU = EOP_IAU2000A) =
    obsposECEF(observatory(x), dtutc2et(x); eop)

obsposECEF(obs::ObservatoryMPC, et::Number; eop::EOPIAU = EOP_IAU2000A) =
    obsposECEF(Val(Symbol(obs.frame)), obs.coords, et; eop)

function obsposECEF(::Val{:MPC}, coords::SVector{3, T}, et::U;
                    eop::EOPIAU = EOP_IAU2000A) where {T <: Real, U <: Number}
    # One with correct type (and order)
    oneU = one(et)
    # Cilindrical coordinates of the observer's position in ITRF (ECEF) frame
    # λ: longitude East of the prime meridian [rad]
    # u: distance from spin axis ( ρ * cos(ϕ')) [km]
    # v: height above equatorial plane ( ρ * sin(ϕ') ) [km]
    # where ϕ' is the geocentric latitude and ρ is the geocentric distance [km]
    λ = deg2rad(coords[1])
    u = coords[2] * RE
    v = coords[3] * RE
    # Cartesian coordinates of the observer's position in ITRF (ECEF) frame [km]
    x = u * cos(λ)
    y = u * sin(λ)
    z = v
    posECEF = SVector{3, U}(x * oneU, y * oneU, z * oneU)

    return posECEF
end

function obsposECEF(::Val{:ICRF_KM}, coords::SVector{3, T}, et::U;
                    eop::EOPIAU = EOP_IAU2000A) where {T <: Real, U <: Number}
    # Zero and one with correct type (and order)
    zeroU, oneU = zero(et), one(et)
    # UTC seconds
    utc_secs = et - tdb_utc(et)
    # Julian days UTC
    jd_utc = JD_J2000 + utc_secs/daysec
    # Cartesian coordinates of the observer's state vector in GCRF (ECI) frame [km]
    posECI = coords * oneU
    velECI = SVector{3, U}(zeroU, zeroU, zeroU)
    accECI = SVector{3, U}(zeroU, zeroU, zeroU)
    posvelECI = OrbitStateVector(jd_utc, posECI, velECI, accECI)
    # Transform state vector from:
    # GCRF: Geocentric Celestial Reference Frame (ECI)
    # to:
    # ITRF: International Terrestrial Reference Frame (ECEF)
    posvelECEF = sv_eci_to_ecef(posvelECI, Val(:GCRF), Val(:ITRF), jd_utc, eop)

    return posvelECEF.r
end

function obsposECEF(::Val{:ICRF_AU}, coords::SVector{3, T}, et::U;
                    eop::EOPIAU = EOP_IAU2000A) where {T <: Real, U <: Number}
    # Zero and one with correct type (and order)
    zeroU, oneU = zero(et), one(et)
    # UTC seconds
    utc_secs = et - tdb_utc(et)
    # Julian days UTC
    jd_utc = JD_J2000 + utc_secs/daysec
    # Cartesian coordinates of the observer's state vector in GCRF (ECI) frame [km]
    posECI = coords * oneU / au
    velECI = SVector{3, U}(zeroU, zeroU, zeroU)
    accECI = SVector{3, U}(zeroU, zeroU, zeroU)
    posvelECI = OrbitStateVector(jd_utc, posECI, velECI, accECI)
    # Transform state vector from:
    # GCRF: Geocentric Celestial Reference Frame (ECI)
    # to:
    # ITRF: International Terrestrial Reference Frame (ECEF)
    posvelECEF = sv_eci_to_ecef(posvelECI, Val(:GCRF), Val(:ITRF), jd_utc, eop)

    return posvelECEF.r
end

function obsposECEF(::Val{:WGS84}, coords::SVector{3, T}, et::U;
                    eop::EOPIAU = EOP_IAU2000A) where {T <: Real, U <: Number}
    # One with correct type (and order)
    oneU = one(et)
    # Geodetic coordinates of the observer's position in WGS84 (ECEF) frame
    # lon: East longitude [rad]
    # lat: latitude [rad]
    # alt: altitude [m]
    lon = deg2rad(coords[1]) * oneU
    lat = deg2rad(coords[2]) * oneU
    alt = coords[3] * oneU
    # Transform position vector from WGS84 (ECEF) to ITRF (ECEF) [km]
    posECEF = geodetic_to_ecef(lat, lon, alt) / 1_000

    return posECEF
end

# TODO: avoid sv_ecef_to_ecef and sv_ecef_to_eci overloads by defining proper product
# between DCMs and Taylor1/TaylorN. The method below has been adapted from
# SatelliteToolboxTransformations.jl, MIT-licensed
# https://github.com/JuliaSpace/SatelliteToolboxTransformations.jl
for EOP in (:Nothing, :EopIau1980, :EopIau2000A)
    @eval begin

        function sv_ecef_to_ecef(
            sv::OrbitStateVector,
            T_ECEF1::Val{:ITRF},
            T_ECEF2::Val{:TIRS},
            jd_utc::Taylor1{TaylorN{Float64}},
            eop_data::$EOP
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
            eop_data::$EOP
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

    end
end

"""
    obsposvelECI(::AbstractOpticalAstrometry; kwargs...)

Return the observer's geocentric cartesian state vector [km, km/sec] in
Earth-Centered Inertial (ECI) reference frame.

# Keyword argument

- `eop::Union{EopIau1980, EopIau2000A}`: Earth Orientation Parameters
    (default: `EOP_IAU2000A`).

# Extended help

By default, the IAU200A Earth orientation model is used to transform
from  Earth-centered, Earth-fixed (ECEF) frame to ECI frame. Other Earth
orientation  models, such as the IAU1976/80 model, can be used by importing
the `SatelliteToolboxTransformations.EopIau1980` type and passing it to the
`eop` keyword argument in the function call.
"""
obsposvelECI(x::AbstractAstrometryObservation; eop::EOPIAU = EOP_IAU2000A) =
    obsposvelECI(observatory(x), dtutc2et(x); eop)

obsposvelECI(obs::ObservatoryMPC, et::Number; eop::EOPIAU = EOP_IAU2000A) =
    obsposvelECI(Val(Symbol(obs.frame)), obs.coords, et; eop)

function obsposvelECI(::Val{:MPC}, coords::SVector{3, T}, et::U;
                      eop::EOPIAU = EOP_IAU2000A) where {T <: Real, U <: Number}
    # Zero of correct type (and order)
    zeroU = zero(et)
    # UTC seconds
    utc_secs = et - tdb_utc(et)
    # Julian days UTC
    jd_utc = JD_J2000 + utc_secs/daysec
    # Cartesian coordinates of the observer's state vector in ITRF (ECEF) frame [km]
    posECEF = obsposECEF(Val(:MPC), coords, et; eop)
    velECEF = SVector{3, U}(zeroU, zeroU, zeroU)
    accECEF = SVector{3, U}(zeroU, zeroU, zeroU)
    posvelECEF = OrbitStateVector(jd_utc, posECEF, velECEF, accECEF)
    # Transform state vector from:
    # ITRF: International Terrestrial Reference Frame (ECEF)
    # to:
    # GCRF: Geocentric Celestial Reference Frame (ECI)
    posvelECI = sv_ecef_to_eci(posvelECEF, Val(:ITRF), Val(:GCRF), jd_utc, eop)

    return vcat(posvelECI.r, posvelECI.v)
end

function obsposvelECI(::Val{:ICRF_KM}, coords::SVector{3, T}, et::U;
                      eop::EOPIAU = EOP_IAU2000A) where {T <: Real, U <: Number}
    # One with correct type (and order)
    zeroU, oneU = zero(et), one(et)
    # Cartesian coordinates of the observer's state vector in GCRF (ECI) frame [km]
    posECI = coords * oneU
    velECI = SVector{3, U}(zeroU, zeroU, zeroU)

    return vcat(posECI, velECI)
end

function obsposvelECI(::Val{:ICRF_AU}, coords::SVector{3, T}, et::U;
                      eop::EOPIAU = EOP_IAU2000A) where {T <: Real, U <: Number}
    # One with correct type (and order)
    zeroU, oneU = zero(et), one(et)
    # Cartesian coordinates of the observer's state vector in GCRF (ECI) frame [km]
    posECI = coords * oneU / au
    velECI = SVector{3, U}(zeroU, zeroU, zeroU)

    return vcat(posECI, velECI)
end

function obsposvelECI(::Val{:WGS84}, coords::SVector{3, T}, et::U;
                      eop::EOPIAU = EOP_IAU2000A) where {T <: Real, U <: Number}
    # Zero of correct type (and order)
    zeroU = zero(et)
    # UTC seconds
    utc_secs = et - tdb_utc(et)
    # Julian days UTC
    jd_utc = JD_J2000 + utc_secs/daysec
    # Cartesian coordinates of the observer's state vector in ITRF (ECEF) frame [km]
    posECEF = obsposECEF(Val(:WGS84), coords, et; eop)
    velECEF = SVector{3, U}(zeroU, zeroU, zeroU)
    accECEF = SVector{3, U}(zeroU, zeroU, zeroU)
    posvelECEF = OrbitStateVector(jd_utc, posECEF, velECEF, accECEF)
    # Transform state vector from:
    # ITRF: International Terrestrial Reference Frame (ECEF)
    # to:
    # GCRF: Geocentric Celestial Reference Frame (ECI)
    posvelECI = sv_ecef_to_eci(posvelECEF, Val(:ITRF), Val(:GCRF), jd_utc, eop)

    return vcat(posvelECI.r, posvelECI.v)
end