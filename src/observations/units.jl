@doc raw"""
    kmsec2auday(pv)

Converts a `[x, y, z, v_x, v_y, v_z]` state vector from km, km/sec to au, au/day.

See also [`auday2kmsec`](@ref).
"""
function kmsec2auday(pv)
    pv /= au          # (km, km/sec) -> (au, au/sec)
    pv[4:6] *= daysec # (au, au/sec) -> (au, au/day)
    return pv
end

@doc raw"""
    auday2kmsec(pv)

Converts a `[x, y, z, v_x, v_y, v_z]` state vector from au, au/day to km, km/sec.

See also [`kmsec2auday`](@ref).
"""
function auday2kmsec(pv)
    pv *= au          # (au, au/day) -> (km, km/day)
    pv[4:6] /= daysec # (km, km/day) -> (km, km/sec)
    return pv
end

@doc raw"""
    julian2etsecs(jd)

Converts `jd` julian days to ephemeris seconds since J2000.

See also [`etsecs2julian`](@ref).
"""
function julian2etsecs(jd)
    return (jd-JD_J2000)*daysec
end

@doc raw"""
    etsecs2julian(et)

Converts `et` ephemeris seconds since J2000 to julian days.

See also [`julian2etsecs`](@ref).
"""
function etsecs2julian(et)
    return JD_J2000 + et/daysec
end

@doc raw"""
    datetime2et(x::DateTime)
    datetime2et(x::RadecMPC{T}) where {T <: AbstractFloat}
    
Retuns the TDB seconds past the J2000 epoch.

See also [`SPICE.str2et`](@ref).
"""
function datetime2et(x::DateTime)
    return str2et(string(x))
end

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