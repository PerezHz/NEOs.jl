# Load TT-TDB (ttmtdb) as a TaylorInterpolant saved in .jld file
const ttmtdb_artifact_path = joinpath(artifact"ttmtdb_DE430_1995_2030", "ttmtdb_DE430_1995_2030_20221103.jld")
const ttmtdb_t0 = JLD.load(ttmtdb_artifact_path, "t0")
const ttmtdb_t = JLD.load(ttmtdb_artifact_path, "t")
const ttmtdb_x_coeffs = JLD.load(ttmtdb_artifact_path, "x_coeffs")
const ttmtdb = TaylorInterpolant(ttmtdb_t0, ttmtdb_t, Taylor1.(ttmtdb_x_coeffs))

@doc raw"""
    loadjpleph()

Load JPL ephemerides (NAIF IDs, DE430 TT-TDB and ephemerides, #197 and #199 solutions for Apophis).

See also [`SPICE.furnsh`](@ref).
"""
function loadjpleph()
    furnsh(
        # NAIF IDs
        joinpath(artifact"naif0012", "naif0012.tls"),
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
    getpv(target::Int, observer::Int, et)

Return the `[x, y, z, v_x, v_y, v_z]` state vector (in units of km, km/sec) at
TDB instant `et` from SPK-formatted ephemeris file with respect to J2000 frame.

See also [`SPICE.spkgeo`](@ref).
"""
function getpv(target::Int, observer::Int, et)
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
    sun_pv(et)

Return the `[x, y, z, v_x, v_y, v_z]` state vector (in units of km, km/sec)
of the Sun at TDB instant `et` with respect to J2000 frame.

See also [`getpv`](@ref).
"""
sun_pv(et) = getpv(10, 0, constant_term(et)) # units: km, km/second

@doc raw"""
    earth_pv(et)

Return the `[x, y, z, v_x, v_y, v_z]` state vector (in units of km, km/sec)
of the Earth at TDB instant `et` with respect to J2000 frame.

See also [`getpv`](@ref).
"""
earth_pv(et) = getpv(399, 0, constant_term(et)) # units: km, km/second

@doc raw"""
    moon_pv(et)

Return the `[x, y, z, v_x, v_y, v_z]` state vector (in units of km, km/sec)
of the Moon at TDB instant `et` with respect to J2000 frame.

See also [`getpv`](@ref).
"""
moon_pv(et) = getpv(301, 0, constant_term(et)) # units: km, km/second

@doc raw"""
    apophis_pv_197(et)

Return the `[x, y, z, v_x, v_y, v_z]` state vector (in units of km, km/sec)
of Apophis at TDB instant `et` from JPL #197 solution with respect to J2000 frame.

See also [`getpv`](@ref).
"""
apophis_pv_197(et) = getpv(9904406, 0, constant_term(et)) # units: km, km/second

@doc raw"""
    apophis_pv_199(et)

Return the `[x, y, z, v_x, v_y, v_z]` state vector (in units of km, km/sec)
of Apophis at TDB instant `et` from JPL #199 solution with respect to J2000 frame.

See also [`getpv`](@ref).
"""
apophis_pv_199(et) = getpv(2099942, 0, constant_term(et)) # units: km, km/second

@doc raw"""
    tt_tdb(et)

Return the difference TT-TDB (in units of sec) at TDB instant `et` with respect to J2000
frame.

See also [`getpv`](@ref).
"""
tt_tdb(et) = getpv(1000000001, 1000000000, constant_term(et))[1] # units: seconds

@doc raw"""
    dtt_tdb(et)

Return the rate of change of TT-TDB (in units of sec/sec) at TDB instant `et` with respect
to J2000 frame.

See also [`getpv`](@ref).
"""
dtt_tdb(et) = getpv(1000000001, 1000000000, constant_term(et))[4] # units: seconds/seconds

