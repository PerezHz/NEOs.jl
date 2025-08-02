"""
    loadjpleph()

Load JPL ephemerides:
- NAIF IDs.
- DE430 TT-TDB and ephemerides.
- #197 and #199 solutions for Apophis.

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

"""
    getposvel(target::Int, observer::Int, et)

Return the cartesian state vector [km, km/sec] in the J2000.0 frame of
body `target` with respect to body `observer` at at TDB instant `et`
from SPK-formatted ephemeris file.

See also [`SPICE.spkgeo`](@ref).
"""
getposvel(target::Int, observer::Int, et) = spkgeo(target, et, "J2000", observer)[1]

# NAIF IDs:
# 0: Solar System Barycenter
# 10: Sun (heliocenter)
# 399: Earth (geocenter)
# 301: Moon
# 2099942: Apophis
# 1000000001 from body 1000000000: TT-TDB
# Here, we follow the convention from the CSPICE, library, that the ephemeris
# time is referred to the J2000.0 frame epoch:
# https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html#Terminology
# argument `et` represents "ephemeris seconds" (TDB seconds) since J2000.0 TDB epoch
# position and velocity are assumed to be returned in km, km/sec, resp., by spkgeo
# https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkgeo_c.html
# (see: Detailed ouput section)

"""
    sunposvel(et)

Return the cartesian state vector [km, km/sec] of the Sun at
TDB instant `et` with respect to the J2000.0 frame.

See also [`getposvel`](@ref).
"""
sunposvel(et) = getposvel(10, 0, cte(et))

"""
    earthposvel(et)

Return the cartesian state vector [km, km/sec] of the Earth at
TDB instant `et` with respect to the J2000.0 frame.

See also [`getposvel`](@ref).
"""
earthposvel(et) = getposvel(399, 0, cte(et))

"""
    moonposvel(et)

Return the cartesian state vector [km, km/sec] of the Moon at
TDB instant `et` with respect to the J2000.0 frame.

See also [`getposvel`](@ref).
"""
moonposvel(et) = getposvel(301, 0, cte(et))

"""
    apophisposvel197(et)

Return the cartesian state vector [km, km/sec] of Apophis at
TDB instant `et` from JPL #197 solution with respect to the
J2000.0 frame.

See also [`getposvel`](@ref).
"""
apophisposvel197(et) = getposvel(9904406, 0, cte(et))

"""
    apophisposvel199(et)

Return the cartesian state vector [km, km/sec] of Apophis at
TDB instant `et` from JPL #199 solution with respect to the
J2000.0 frame.

See also [`getposvel`](@ref).
"""
apophisposvel199(et) = getposvel(2099942, 0, cte(et))

"""
    tt_tdb(et)

Return the difference TT-TDB [sec] at TDB instant `et`
with respect to the J2000.0 frame.

See also [`getposvel`](@ref).
"""
tt_tdb(et) = getposvel(1000000001, 1000000000, cte(et))[1]

"""
    dtt_tdb(et)

Return the rate of change of TT-TDB [sec/sec] at TDB instant `et`
with respect to the J2000.0 frame.

See also [`getposvel`](@ref).
"""
dtt_tdb(et) = getposvel(1000000001, 1000000000, cte(et))[4]