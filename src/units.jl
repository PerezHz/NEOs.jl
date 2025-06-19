"""
    julian2etsecs(::Number)

Convert a Julian date in TDB time scale to TDB seconds since J2000.

See also [`etsecs2julian`](@ref).
"""
julian2etsecs(jdtdb::Number) = (jdtdb - JD_J2000) * daysec

"""
    etsecs2julian(::Number)

Convert TDB seconds since J2000 to Julian date in TDB time scale.

See also [`julian2etsecs`](@ref).
"""
etsecs2julian(et::Number) = JD_J2000 + et / daysec

"""
    dtutc2et(::DateTime)

Convert a UTC date to TDB seconds past the J2000 epoch.

See also [`et2dtutc`](@ref).
"""
function dtutc2et(dtutc::DateTime)
    # UTC seconds since J2000.0 epoch
    utc_seconds = Millisecond(dtutc - DateTime(2000,1,1,12)).value / 1000
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

"""
    dtutc2jdtdb(::DateTime)

Convert a UTC date to Julian date in TDB time scale.

See also [`jdtdb2dtutc`](@ref).
"""
function dtutc2jdtdb(dtutc::DateTime)
    et = dtutc2et(dtutc) # TDB seconds since J2000.0
    return etsecs2julian(et) # JDTDB
end

"""
    et2dtutc(::Number)

Convert TDB seconds past the J2000 epoch to a UTC date.

See also [`dtutc2et`](@ref).
"""
function et2dtutc(et::Number)
    # UTC seconds past J2000
    utc = et - tdb_utc(et)
    # UTC milliseconds past J2000
    mill = Millisecond(round(Int, utc * 1_000)).value
    # UTC DateTime
    return epochms2datetime(EPOCHMSJ2000 + mill)
end

"""
    jdtdb2dtutc(::Number)

Convert a Julian date in TDB time scale to a UTC date.

See also [`dtutc2jdtdb`](@ref).
"""
function jdtdb2dtutc(jdtdb::Number)
    et = julian2etsecs(jdtdb) # TDB seconds since J2000.0
    return et2dtutc(et) # UTC DateTime
end

"""
    et_to_200X(::Number)

Convert TDB seconds since J2000 to TDB years `200X`.
"""
et_to_200X(et::Number) = 2000 + et/daysec/yr

"""
    days_to_200X(::Number)

Convert TDB days since J2000 to TDB years `200X`.
"""
days_to_200X(d::Number) = 2000 + d/yr

"""
    dtutc_to_200X(::DateTime)

Convert a UTC date to TDB years `200X`.
"""
dtutc_to_200X(x::DateTime) = et_to_200X(dtutc2et(x))

"""
    dtutc2days(::DateTime)

Convert a UTC date to TDB days since J2000.

See also [`days2dtutc`](@ref).
"""
dtutc2days(x::DateTime) = dtutc2jdtdb(x) - JD_J2000

"""
    days2dtutc(::Number)

Convert TDB days since J2000 to a UTC date.

See also [`dtutc2days`](@ref).
"""
days2dtutc(d::Number) = jdtdb2dtutc(d + JD_J2000)

"""
    tdb_utc(::Number)

Return the difference (TDB-UTC) given a number of TDB seconds past
the J2000.0 epoch

!!! reference
    This function is useful to convert TDB to UTC via UTC + (TDB-UTC)
    and viceversa. It does not include the correction due to the position
    of the measurement station ``v_E.(r_S-r_E)/c^2`` (Folkner et al. 2014;
    Moyer, 2003).

# Extended help

The (TDB-UTC) difference is given by:
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

"""
    ttmtdb_tt(::Real; kwargs...)

Return the difference (TT-TDB), where TDB is the Barycentric Dynamical Time
(Temps Dynamique Barycentrique) and TT is the Terrestrial Time, given a number
of TT seconds past the J2000.0 epoch.

# Keyword arguments

- `niter::Int`: number of iterations (default: `5`).

!!! reference
    This function is useful during the reduction of observations, when the TDB
    instant is computed from a known UTC instant. The computed value does not
    include the correction due to the position of the measurement station
    ``v_E.(r_S-r_E)/c^2`` (Folkner et al. 2014; Moyer, 2003).
"""
function ttmtdb_tt(tt::Real; niter::Int = 5)
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

"""
    rad2arcsec(::Number)

Convert radians to arcseconds.

See also [`arcsec2rad`](@ref) and [`mas2rad`](@ref).
"""
rad2arcsec(x::Number) = 3600 * rad2deg(x)

"""
    arcsec2rad(::Number)

Convert arcseconds to radians.

See also [`rad2arcsec`](@ref) and [`mas2rad`](@ref).
"""
arcsec2rad(x::Number) = deg2rad(x / 3600)

"""
    mas2rad(::Number)

Convert milli-arcseconds to radians.

See also [`rad2arcsec`](@ref) and [`arcsec2rad`](@ref).
"""
mas2rad(x::Number) = arcsec2rad(x / 1000)

"""
    range2delay(::Number)

Convert radar range [km] to time-delay [us].
"""
range2delay(ρ::Number) = 2 * ρ / c_km_per_us

raw"""
    rangerate2doppler(::Number, ::Number)

Convert radar range rate [km/day] and transmitter frequency [MHz]
to Doppler shift [Hz].

# Extended help

The Doppler shift is given by:
```math
f_d = -\frac{2\dot{\rho}f_t}{c - \dot{\rho}},
```
where ``\dot{\rho}`` is the range rate and ``f_t`` the transmitter
frequency. However, this function uses only the first term in the
expansion of ``(c - \dot{\rho})^{-1}``.
"""
rangerate2doppler(v::Number, f_t::Number) = -2 * v * f_t * 1e6 / c_km_per_day