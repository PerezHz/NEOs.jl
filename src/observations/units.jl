@doc raw"""
    julian2etsecs(jdtdb::T) where {T <: Number}

Convert `jdtdb`, a Julian date in TDB time scale, to TDB seconds since J2000.

See also [`etsecs2julian`](@ref).
"""
function julian2etsecs(jdtdb::T) where {T <: Number}
    return (jdtdb - JD_J2000) * daysec
end

@doc raw"""
    etsecs2julian(et::T) where {T <: Number}

Convert `et` TDB seconds since J2000 to Julian date in TDB time scale.

See also [`julian2etsecs`](@ref).
"""
function etsecs2julian(et::T) where {T <: Number}
    return JD_J2000 + et / daysec
end

@doc raw"""
    dtutc2et(dtutc::DateTime)
    dtutc2et(dtutc::T) where {T <: AbstractAstrometry}

Convert a UTC `DateTime` to TDB seconds past the J2000 epoch.
"""
function dtutc2et(dtutc::DateTime)
    # UTC seconds since J2000.0 epoch
    utc_seconds = Millisecond(dtutc - DateTime(2000,1,1,12)).value/1000
    # TAI - UTC
    tai_utc = get_Δat(datetime2julian(dtutc))
    # TT - TAI
    tt_tai = 32.184
    # TT - UTC = (TT-TAI) + (TAI-UTC)
    tt_utc = tt_tai +tai_utc
    # TT seconds  = UTC seconds + (TT-UTC)
    tt_seconds = utc_seconds + tt_utc
    # TDB seconds = TT seconds + (TDB-TT)
    return tt_seconds - ttmtdb_tt(tt_seconds)
end

dtutc2et(astrometry::T) where {T <: AbstractAstrometry} = dtutc2et(astrometry.date)

@doc raw"""
    dtutc2jdtdb(dtutc::DateTime)

Convert a UTC `DateTime` to Julian date in TDB time scale.
"""
function dtutc2jdtdb(dtutc::DateTime)
    et = dtutc2et(dtutc) # TDB seconds since J2000.0
    return etsecs2julian(et) # JDTDB
end

@doc raw"""
    et2dtutc(et::T) where {T <: Number}

Convert `et` TDB seconds past the J2000 epoch to a UTC `DateTime`.
"""
function et2dtutc(et::T) where {T <: Number}
    # UTC seconds past J2000
    utc = et - tdb_utc(et)
    # UTC milliseconds past J2000
    mill = Millisecond(round(Int, utc * 1_000)).value
    # UTC DateTime
    return epochms2datetime(EPOCHMSJ2000 + mill)
end

@doc raw"""
    jdtdb2dtutc(jdtdb::T) where {T <: Number}

Convert `jdtdb`, a Julian date in TDB time scale, to a UTC `DateTime`.
"""
function jdtdb2dtutc(jdtdb::T) where {T <: Number}
    et = julian2etsecs(jdtdb) # TDB seconds since J2000.0
    return et2dtutc(et) # UTC DateTime
end

@doc raw"""
    et_to_200X(et::T) where {T <: Number}

Convert `et` TDB seconds since J2000 to TDB years `200X`.
"""
et_to_200X(et::T) where {T <: Number} = 2000 + et/daysec/yr

@doc raw"""
    days_to_200X(d::T) where {T <: Number}

Convert `d` TDB days since J2000 to TDB years `200X`.
"""
days_to_200X(d::T) where {T <: Number} = 2000 + d/yr

@doc raw"""
    datetime_to_200X(x::DateTime)

Convert `x`, a UTC `DateTime`, to TDB years `200X`.
"""
datetime_to_200X(x::DateTime) = et_to_200X(dtutc2et(x))

@doc raw"""
    datetime2days(x::DateTime)

Convert `x`, a UTC `DateTime`, to TDB days since J2000.
"""
datetime2days(x::DateTime) = dtutc2jdtdb(x) - JD_J2000

@doc raw"""
    days2datetime(d::T) where {T <: Number}

Convert `d` TDB days since J2000 to a UTC `DateTime`.
"""
days2datetime(d::T) where {T <: Number} = jdtdb2dtutc(d + JD_J2000)

@doc raw"""
    tdb_utc(et::T) where {T <: Number}

Given `et`, a number of TDB seconds past J2000.0 epoch, compute
the difference (TDB-UTC)
```math
\begin{align*}
TDB-UTC & = (TDB-TAI) + (TAI-UTC) \\
        & = (TDB-TT) + (TT-TAI) + (TAI-UTC) \\
        & = (TDB-TT) + 32.184 s + ΔAT,
\end{align*}
```
where TDB is the Barycentric Dynamical Time (Temps Dynamique Barycentrique),
TT is the Terrestrial Time, TAI is the International Atomic Time, and UTC is
the Coordinated Universal Time.

# Arguments

- `et::Number`: TDB seconds since J2000.0.

!!! reference
    This function is useful to convert TDB to UTC via UTC + (TDB-UTC) and viceversa.
    It does not include the correction due to the position of the measurement station
    ``v_E.(r_S-r_E)/c^2`` (Folkner et al. 2014; Moyer, 2003).
"""
function tdb_utc(et::T) where {T <: Number}
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
    ttmtdb_tt(tt::T; niter::Int = 5) where {T <: Real}

Given `tt`, a number of TT seconds past J2000.0 epoch, compute
the difference (TT-TDB), where TDB is the Barycentric Dynamical Time (Temps Dynamique
Barycentrique) and TT is the Terrestrial Time.

# Arguments

- `tt::T`: tt seconds since J2000.0.

!!! reference
    This function is useful during the reduction of observations, when the TDB instant
    is computed from a known UTC instant. The computed value does not include the
    correction due to the position of the measurement station ``v_E.(r_S-r_E)/c^2``
    (Folkner et al. 2014; Moyer, 2003).
"""
function ttmtdb_tt(tt::T; niter::Int = 5) where {T <: Real}
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
    rad2arcsec(x::T) where {T <: Number}

Convert radians to arcseconds.

See also [`arcsec2rad`](@ref) and [`mas2rad`](@ref).
"""
rad2arcsec(x::T) where {T <: Number} = 3600 * rad2deg(x)

@doc raw"""
    arcsec2rad(x::T) where {T <: Number}

Convert arcseconds to radians.

See also [`rad2arcsec`](@ref) and [`mas2rad`](@ref).
"""
arcsec2rad(x::T) where {T <: Number} = deg2rad(x / 3600)

@doc raw"""
    mas2rad(x::T) where {T <: Number}

Convert milli-arcseconds to radians.

See also [`rad2arcsec`](@ref) and [`arcsec2rad`](@ref).
"""
mas2rad(x::T) where {T <: Number} = arcsec2rad(x / 1000)