@doc raw"""
    julian2etsecs(jd)

Convert the Julian date `jd` to ephemeris seconds since J2000.

See also [`etsecs2julian`](@ref).
"""
function julian2etsecs(jd)
    return (jd-JD_J2000)*daysec
end

@doc raw"""
    etsecs2julian(et)

Convert `et` ephemeris seconds since J2000 to Julian date.

See also [`julian2etsecs`](@ref).
"""
function etsecs2julian(et)
    return JD_J2000 + et/daysec
end

@doc raw"""
    datetime2et(x::DateTime)
    datetime2et(x::T) where {T <: AbstractAstrometry}

Convert a UTC `DateTime` to TDB seconds past the J2000 epoch.
"""
function datetime2et(x::DateTime)
    # UTC seconds since J2000.0 epoch
    utc_seconds = Millisecond(x - DateTime(2000,1,1,12)).value/1000
    # TAI - UTC
    tai_utc = get_Δat(datetime2julian(x))
    # TT - TAI
    tt_tai = 32.184
    # TT - UTC = (TT-TAI) + (TAI-UTC)
    tt_utc = tt_tai +tai_utc
    # TT seconds  = UTC seconds + (TT-UTC)
    tt_seconds = utc_seconds + tt_utc
    # TDB seconds = TT seconds + (TDB-TT)
    return tt_seconds - ttmtdb_tt(tt_seconds)
end

datetime2et(x::AbstractAstrometry) = datetime2et(x.date)

@doc raw"""
    et_to_200X(et::T) where {T <: Number}

Convert `et` ephemeris seconds since J2000 to years `200X`.
"""
et_to_200X(et::Number) = 2000 + et/daysec/yr

@doc raw"""
    days_to_200X(d::T) where {T <: Number}

Convert `d` days since J2000 to years `200X`.
"""
days_to_200X(d::Number) = 2000 + d/yr

@doc raw"""
    datetime_to_200X(x::DateTime)

Convert `DateTime` `x` to years `200X`.
"""
datetime_to_200X(x::DateTime) = et_to_200X(datetime2et(x))

@doc raw"""
    datetime2days(x::DateTime)

Convert `DateTime` `x` to days since J2000.
"""
datetime2days(x::DateTime) = datetime2julian(x) - JD_J2000

@doc raw"""
    days2datetime(d::T) where {T <: Number}

Convert `d` days since J2000 to `DateTime`.
"""
days2datetime(d::Number) = julian2datetime(d + JD_J2000)

@doc raw"""
    tdb_utc(et::Number)

Given `et`, a number of TDB seconds past J2000.0 epoch, compute
the difference (TDB-UTC)
```math
\begin{align*}
TDB-UTC & = (TDB-TAI) + (TAI-UTC) \\
        & = (TDB-TT) + (TT-TAI) + (TAI-UTC) \\
        & = (TDB-TT) + 32.184 s + ΔAT,
\end{align*}
```
where TDB is the Barycentric Dynamical Time (Temps Dynamique Barycentrique), TT is the Terrestrial Time,
TAI is the International Atomic Time, and UTC is the Coordinated Universal Time.

This function is useful to convert TDB to UTC via UTC + (TDB-UTC) and viceversa. It does
not include the correction due to the position of the measurement station ``v_E.(r_S-r_E)/c^2``
(Folkner et al. 2014; Moyer, 2003).

# Arguments

- `et::Number`: TDB seconds since J2000.0.
"""
function tdb_utc(et::Number)
    # TT-TDB
    tt_tdb_et = ttmtdb(et/daysec)
    # TT-TAI
    tt_tai = 32.184
    # TDB - TAI = (TT-TAI) + (TDB-TT) = (TDB-TT) + 32.184 s
    tdb_tai = tt_tai - tt_tdb_et

    # TAI seconds since J2000.0; used to determine ΔAT = TAI - UTC
    tai_secs = cte(cte(et)) - cte(cte(tdb_tai))
    # Julian date corresponding to tai_secs
    jd_tai = JD_J2000 + tai_secs/daysec
    # ΔAT = TAI - UTC
    tai_utc = get_Δat(jd_tai)
    # TDB-UTC = (TDB-TAI) + (TAI-UTC) = (TDB-TT) + 32.184 s + ΔAT
    return tdb_tai + tai_utc
end

@doc raw"""
    tt_tdb_tt(tt::Real)

Given `tt`, a number of TT seconds past J2000.0 epoch, compute
the difference (TT-TDB), where TDB is the Barycentric Dynamical Time (Temps Dynamique
Barycentrique) and TT is the Terrestrial Time.

This function is useful during the reduction of observations, when the TDB instant is
computed from a known UTC instant. The computed value does not include the
correction due to the position of the measurement station ``v_E.(r_S-r_E)/c^2``
(Folkner et al. 2014; Moyer, 2003).

# Arguments

- `tt::Real`: tt seconds since J2000.0.
"""
function ttmtdb_tt(tt::Real; niter=5)
    # Ansatz: TDB - TT = 0
    ttmtdb_order = ttmtdb.x[1].order
    tdb = Taylor1([tt,one(tt)], ttmtdb_order)
    for _ in 1:niter
        ttmtdb_tdb = ttmtdb(tdb/daysec)
        # we look for TDB* such that TT-TDB* = (TT-TDB)(TDB*)
        y = tt - tdb - ttmtdb_tdb
        dy = - 1 - TaylorSeries.differentiate(ttmtdb_tdb)
        # perform Newton iteration
        tdb[0] -= cte(y/dy)
    end
    return ttmtdb(tdb[0]/daysec)
end

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