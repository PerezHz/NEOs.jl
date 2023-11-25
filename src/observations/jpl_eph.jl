# Load TT-TDB (ttmtdb) from jld2 file
@doc raw"""
    loadjpleph()

Load JPL ephemerides (NAIF IDs, DE430 TT-TDB and ephemerides, #197 and #199 solutions for
Apophis).

See also [`SPICE.furnsh`](@ref).
"""
function loadjpleph()
    furnsh(
        # JPL DE430 TT-TDB
        joinpath(artifact"TTmTDBde430", "TTmTDB.de430.19feb2015.bsp"),
        # JPL DE430 ephemerides
        joinpath(artifact"de430", "de430_1850-2150.bsp"),
        # JPL #197 solution for Apophis
        joinpath(artifact"a99942", "a99942_s197.bsp"),
        # JPL #199 solution for Apophis
        joinpath(artifact"a99942", "a99942_s199.bsp"),
    )
end

@doc raw"""
    getposvel(target::Int, observer::Int, et)

Return the `[x, y, z, v_x, v_y, v_z]` state vector (in units of km, km/sec) at
TDB instant `et` from SPK-formatted ephemeris file with respect to J2000 frame.

See also [`SPICE.spkgeo`](@ref).
"""
function getposvel(target::Int, observer::Int, et)
    return spkgeo(target, et, "J2000", observer)[1] # units: km,km/sec
end

# NAIF IDs:
# 0: Solar System Barycenter
# 10: Sun (heliocenter)
# 2099942: Apophis
# 399: Earth (geocenter)
# 301: Moon
# 1000000001 from body 1000000000: TT-TDB
# Here, we follow the convention from the CSPICE, library, that the ephemeris
# time is referred to the J2000 frame epoch:
# https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html#Terminology
# argument `et` represents "ephemeris seconds" (TDB seconds) since J2000.0 TDB epoch
# position and velocity are assumed to be returned in km, km/sec, resp., by spkgeo
# https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkgeo_c.html (see: Detailed ouput section)

@doc raw"""
    sunposvel(et)

Return the `[x, y, z, v_x, v_y, v_z]` state vector (in units of km, km/sec)
of the Sun at TDB instant `et` with respect to J2000 frame.

See also [`getposvel`](@ref).
"""
sunposvel(et) = getposvel(10, 0, cte(et)) # units: km, km/second

@doc raw"""
    earthposvel(et)

Return the `[x, y, z, v_x, v_y, v_z]` state vector (in units of km, km/sec)
of the Earth at TDB instant `et` with respect to J2000 frame.

See also [`getposvel`](@ref).
"""
earthposvel(et) = getposvel(399, 0, cte(et)) # units: km, km/second

@doc raw"""
    moonposvel(et)

Return the `[x, y, z, v_x, v_y, v_z]` state vector (in units of km, km/sec)
of the Moon at TDB instant `et` with respect to J2000 frame.

See also [`getposvel`](@ref).
"""
moonposvel(et) = getposvel(301, 0, cte(et)) # units: km, km/second

@doc raw"""
    apophisposvel197(et)

Return the `[x, y, z, v_x, v_y, v_z]` state vector (in units of km, km/sec)
of Apophis at TDB instant `et` from JPL #197 solution with respect to J2000 frame.

See also [`getposvel`](@ref).
"""
apophisposvel197(et) = getposvel(9904406, 0, cte(et)) # units: km, km/second

@doc raw"""
    apophisposvel199(et)

Return the `[x, y, z, v_x, v_y, v_z]` state vector (in units of km, km/sec)
of Apophis at TDB instant `et` from JPL #199 solution with respect to J2000 frame.

See also [`getposvel`](@ref).
"""
apophisposvel199(et) = getposvel(2099942, 0, cte(et)) # units: km, km/second

@doc raw"""
    tt_tdb(et)

Return the difference TT-TDB (in units of sec) at TDB instant `et` with respect to J2000
frame.

See also [`getposvel`](@ref).
"""
tt_tdb(et) = getposvel(1000000001, 1000000000, cte(et))[1] # units: seconds

@doc raw"""
    dtt_tdb(et)

Return the rate of change of TT-TDB (in units of sec/sec) at TDB instant `et` with respect
to J2000 frame.

See also [`getposvel`](@ref).
"""
dtt_tdb(et) = getposvel(1000000001, 1000000000, cte(et))[4] # units: seconds/seconds

@doc raw"""
    loadpeeph(eph::TaylorInterpolant{Float64, Float64, 2} = sseph, t_0::T = sseph.t0,
              t_f::S = sseph.t0 + sseph.t[end]) where {T, S <: Real}

Load ephemeris produced by `PlanetaryEphemeris.jl` in timerange `[t_0, t_f] ⊆ [0.0, 36525.0]`
where `t` must have units of TDB days since J2000. The available options for `eph` are:

- `NEOs.sseph`: Solar system ephemeris.
- `NEOs.acceph`: accelerations ephemeris.
- `NEOs.poteph`: newtonian potentials ephemeris.

!!! warning
    Running this function for the first time will download the `sseph_p100` artifact (885 MB)
    which can take several minutes.
"""
function loadpeeph(eph::TaylorInterpolant{Float64, Float64, 2} = sseph, t_0::T = sseph.t0,
                   t_f::S = sseph.t0 + sseph.t[end]) where {T, S <: Real}
    @assert 0.0 ≤ t_0 ≤ t_f ≤ 36525.0
    i_0 = searchsortedlast(eph.t, t_0)
    i_f = searchsortedfirst(eph.t, t_f)
    return TaylorInterpolant(eph.t0, eph.t[i_0:i_f], eph.x[i_0:i_f-1, :])
end

@doc raw"""
    bwdfwdeph(et::Union{T, Taylor1{T}, TaylorN{T}, Taylor1{TaylorN{T}}},
              bwd::TaylorInterpolant{T, U, 2},
              fwd::TaylorInterpolant{T, U, 2}
              ) where {T <: AbstractFloat, U <: Union{T, TaylorN{T}}}

Paste a backward and a forward integration, evaluate at `et` and convert from 
[au, au/day] -> [km, km/sec].
"""
function bwdfwdeph(et::Union{T, Taylor1{T}, TaylorN{T}, Taylor1{TaylorN{T}}},
                   bwd::TaylorInterpolant{T,U,2}, 
                   fwd::TaylorInterpolant{T,U,2}
                   ) where {T <: AbstractFloat, U <: Union{T, TaylorN{T}}}
    @assert bwd.t0 == fwd.t0 "Backward and forward TaylorInterpolant initial times must match"
    t = et/daysec
    t0 = bwd.t0
    if t <= t0
        return auday2kmsec(bwd(t))
    else
        return auday2kmsec(fwd(t))
    end
end
