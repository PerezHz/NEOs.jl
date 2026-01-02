# [Astrometry](@id Astrometry)

NEOs can fetch, read and write astrometric observations of Near-Earth Objects. Currently,
it recognizes two types of astrometry: optical and radar. Each astrometry type is distributed
in various formats, depending on the source. Hence, NEOs defines an interface to work
with different formats in a uniform way.

## Optical astrometry

An optical observation is a tuple of topocentric right ascension and declination. NEOs can
parse optical astrometry in three different formats:
- The Minor Planet Center's [80-column](https://minorplanetcenter.net/iau/info/OpticalObs.html) format,
    represented internally as `OpticalMPC80`.
- The Minor Planet Center's [ADES](https://minorplanetcenter.net/mpcops/documentation/ades/) format,
    represented internally as `OpticalADES`.
- The OrbFit's [RWO](https://neo.ssa.esa.int/objects/help#obs) format,
    represented internally as `OpticalRWO`.

To fetch the optical astrometry of a particular object, we use the following set of functions:
```julia
using NEOs

# optical = fetch_optical_[format]([designation], [source])

# 80-column
optical = fetch_optical_mpc80("2024 YR4", MPC)
optical = fetch_optical_mpc80("P12jLGB", NEOCP)
# ADES
optical = fetch_optical_ades("2024 YR4", MPC)
optical = fetch_optical_ades("P12jLGB", NEOCP)
# RWO
optical = fetch_optical_rwo("2024 YR4", NEOCC)
optical = fetch_optical_rwo("2024 YR4", NEODyS2)
```
which leverage the Minor Planet Center's [Observations API](https://minorplanetcenter.net/mpcops/documentation/observations-api/) and the NEOCC's [Automated Data Access](https://neo.ssa.esa.int/computer-access).

Also, we can read and write observations from a file in a similar fashion:
```julia
# Read
optical = read_optical_mpc80("2024YR4.txt")
optical = read_optical_ades("2024YR4.xml")
optical = read_optical_rwo("2024YR4.rwo")
# Write
write_optical_mpc80(optical, "2024YR4.txt")
write_optical_ades(optical, "2024YR4.xml")
write_optical_rwo(optical, "2024YR4.rwo")
```

Regardless of the format, we can access the most important information in an optical observation
using the same functions:
```julia
x = optical[1]
date(x)             # Date of observation [UTC]
ra(x)               # Right ascension [rad]
dec(x)              # Declination [rad]
mag(x)              # Apparent magnitude
band(x)             # Photometric band
observatory(x)      # Observing telescope
catalogue(x)        # Reference star catalogue
rms(x)              # A-priori formal RMS [arcsec].
debias(x)           # Debiasing factor [arcsec].
```

## Radar astrometry

A radar observation is a measurement of time-delay or Doppler shift. NEOs can
parse radar astrometry in two different formats:
- The JPL's [radar](https://ssd-api.jpl.nasa.gov/doc/sb_radar.html) format,
    represented internally as `RadarJPL`.
- The OrbFit's [RWO](https://neo.ssa.esa.int/objects/help#obs) format,
    represented internally as `RadarRWO`.

To fetch the radar astrometry of a particular object, we use the following set of functions:
```julia
using NEOs

# radar = fetch_radar_[format]([designation], [source])

# JPL
radar = fetch_radar_jpl("des" => "2024 YR4", JPL)
# RWO
radar = fetch_radar_rwo("2024 YR4", NEOCC)
radar = fetch_radar_rwo("2024 YR4", NEODyS2)
```
which leverage the JPL's [Small-Body Radar Astrometry API](https://ssd-api.jpl.nasa.gov/doc/sb_radar.html)
and the NEOCC's [Automated Data Access](https://neo.ssa.esa.int/computer-access).

Also, we can read and write observations from a file in a similar fashion:
```julia
# Read
radar = read_radar_jpl("2024YR4.json")
radar = read_radar_rwo("2024YR4.rwo")
# Write
write_radar_jpl(radar, "2024YR4.json")
write_radar_rwo(radar, "2024YR4.rwo")
```

Regardless of the format, we can access the most important information in a radar observation
using the same functions:
```julia
x = radar[1]
date(x)             # Date of observation [UTC]
measure(x)          # Time-delay [us] or Doppler shift [Hz]
frequency(x)        # Transmitter frequency [MHz]
observatory(x)      # Observing station
rms(x)              # A-priori formal RMS [same units as `measure(x)`]
debias(x)           # Debiasing factor [same units as `measure(x)`]
```