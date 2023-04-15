@doc raw"""
    julian2etsecs(jd)

Convert `jd` julian days to ephemeris seconds since J2000.

See also [`etsecs2julian`](@ref).
"""
function julian2etsecs(jd)
    return (jd-JD_J2000)*daysec
end

@doc raw"""
    etsecs2julian(et)

Convert `et` ephemeris seconds since J2000 to julian days.

See also [`julian2etsecs`](@ref).
"""
function etsecs2julian(et)
    return JD_J2000 + et/daysec
end

@doc raw"""
    datetime2et(x::DateTime)
    datetime2et(x::T) where {T <: AbstractObservation}
    
Retun the TDB seconds past the J2000 epoch.

See also [`SPICE.str2et`](@ref).
"""
function datetime2et(x::DateTime)
    return str2et(string(x))
end

datetime2et(x::T) where {T <: AbstractObservation} = datetime2et(x.date)

@doc raw"""
    et_to_200X(et::T) where {T <: Number}

Convert `et` ephemeris seconds since J2000 to years `200X`. 
"""
et_to_200X(et::T) where {T <: Number} = 2000 + et/daysec/yr

@doc raw"""
    days_to_200X(jd::T) where {T <: Number}

Convert `jd` julian days since J2000 to years `200X`. 
"""
days_to_200X(jd::T) where {T <: Number} = 2000 + jd/yr

@doc raw"""
    datetime_to_200X(x::DateTime)

Convert `x` to years `200X`. 
"""
datetime_to_200X(x::DateTime) = et_to_200X(datetime2et(x))

@doc raw"""
    datetime2days(x::DateTime)

Convert `x` to julian days since 2000. 
"""
datetime2days(x::DateTime) = datetime2julian(x) - JD_J2000

@doc raw"""
    days2datetime(x::DateTime)

Convert `jd` julian days since 2000 to `DateTime`. 
"""
days2datetime(jd::T) where {T <: Number} = julian2datetime(jd + JD_J2000)

@doc raw"""
    tdb_utc(et::T) where {T<:Number}

Auxiliary function to compute (TDB-UTC)
```math
\begin{align*}
TDB-UTC & = (TDB-TAI) + (TAI-UTC) \\
        & = (TDB-TT) + (TT-TAI) + (TAI-UTC) \\
        & = (TDB-TT) + 32.184 s + ΔAT,
\end{align*}
```
where TDB is the Solar System barycentric ephemeris time, TT is the Terrestrial time,
TAI is the International Atomic Time, and UTC is the Coordinated Universal Time.

This function is useful to convert TDB to UTC via UTC + (TDB-UTC) and viceversa. This function
does not include correction due to position of measurement station ``v_E(r_S.r_E)/c^2``
(Folkner et al. 2014; Moyer, 2003).

# Arguments

- `et::T`: TDB seconds since J2000.0.
"""
function tdb_utc(et::T) where {T<:Number}
    # TT-TDB
    tt_tdb_et = ttmtdb(et)
    # TT-TAI
    tt_tai = 32.184

    et_00 = constant_term(constant_term(et))
    # Used only to determine ΔAT; no high-precision needed
    utc_secs = et_00 - deltet(et_00, "ET")
    # ΔAT
    jd_utc = JD_J2000 + utc_secs/daysec
    tai_utc = get_ΔAT(jd_utc)
    # TDB-UTC = (TDB-TT) + (TT-TAI) + (TAI-UTC) = (TDB-TT) + 32.184 s + ΔAT
    return (tt_tai + tai_utc) - tt_tdb_et
end

# TODO: add tdb_utc(utc) method!!!
# strategy: given UTC, do UTC + (TT-TAI) + (TAI-UTC) to get TT
# then, use ttmtdb(et) function iteratively (Newton) to compute TDB

# function dtutc2et(t_utc::DateTime)
#     tt_tai = 32.184
#     jd_utc = datetime2julian(t_utc)
#     fd_utc = (jd_utc+0.5) - floor(jd_utc+0.5)
#     j, tai_utc = iauDat(year(t_utc), month(t_utc), day(t_utc), fd_utc)
#     return et
# end

@doc raw"""
    rad2arcsec(x)

Convert radians to arcseconds. 

See also [`arcsec2rad`](@ref) and [`mas2rad`](@ref).
"""
rad2arcsec(x) = 3600 * rad2deg(x) # rad2deg(rad) -> deg; 3600 * deg -> arcsec

@doc raw"""
    arcsec2rad(x)

Convert arcseconds to radians. 

See also [`rad2arcsec`](@ref) and [`mas2rad`](@ref).
"""
arcsec2rad(x) = deg2rad(x / 3600) # arcsec/3600 -> deg; deg2rad(deg) -> rad

@doc raw"""
    mas2rad(x)

Convert milli-arcseconds to radians. 

See also [`rad2arcsec`](@ref) and [`arcsec2rad`](@ref).
"""
mas2rad(x) = arcsec2rad(x / 1000) # mas/1000 -> arcsec; arcsec2rad(arcsec) -> rad