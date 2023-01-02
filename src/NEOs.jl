module NEOs

# __precompile__(false)

# Units 
export kmsec2auday, auday2kmsec, julian2etsecs, etsecs2julian, datetime2et, rad2arcsec, arcsec2rad, mas2rad
# JPL Ephemerides 
export loadjpleph, sun_pv, earth_pv, moon_pv, apophis_pv_197, apophis_pv_199
# Osculating 
export OsculatingElements, pv2kep
# CatalogueMPC
export unknowncat, read_catalogues_mpc, parse_catalogues_mpc, write_catalogues_mpc, update_catalogues_mpc, search_cat_code
# ObservatoryMPC 
export hascoord, read_observatories_mpc, parse_observatories_mpc, write_observatories_mpc, update_observatories_mpc,
       unknownobs, isunknown, search_obs_code
# RadecMPC
export ra, dec, read_radec_mpc, parse_radec_mpc, search_circulars_mpc, write_radec_mpc, date 
# Topocentric
export obs_pos_ECEF, obs_pv_ECI 
# Process radec 
export compute_radec, w8sveres17, radec_astrometry
# Gauss method 
export gauss_method
# Propagate 
export propagate_dense, propagate_lyap, propagate_root

export delay_doppler, ismonostatic, t2c_rotation_iau_76_80,
    process_radar_data_jpl, RadarDataJPL, RNp1BP_pN_A_J23E_J2S_ng_eph!, RNp1BP_pN_A_J23E_J2S_ng_eph_threads!,
    utcepoch, delay, delay_sigma, delay_units, doppler, doppler_sigma,
    doppler_units, freq, rcvr, xmit, bouncepoint, valsecchi_circle,
    nrms, chi2, newtonls, newtonls_6v, diffcorr,
    newtonls_Q, readmp, bopik, yarkp2adot,
    x0_JPL_s197, x0_JPL_s199

import Base: hash, ==, show, isless, isnan
import Dates: DateTime
import Statistics: mean
import ReferenceFrameRotations: orthonormalize
import PlanetaryEphemeris, SatelliteToolbox, RemoteFiles

using Distributed, JLD, JLD2, TaylorIntegration, Printf, DelimitedFiles, Test, LinearAlgebra,
      Dates, EarthOrientation, SPICE, Quadmath, LazyArtifacts, DataFrames, TaylorSeries,
      InteractiveUtils
using PlanetaryEphemeris: daysec, su, ea, α_p_sun, δ_p_sun, t2c_jpl_de430, pole_rotation, 
      au, c_au_per_day, R_sun, c_cm_per_sec, c_au_per_sec, yr, RE, TaylorInterpolant, Rx, 
      Ry, Rz, semimajoraxis, eccentricity, inclination, longascnode, argperi, timeperipass, 
      nbodyind
using Healpix: ang2pixRing, Resolution
using SatelliteToolbox: nutation_fk5, J2000toGMST, rECEFtoECI, get_ΔAT, JD_J2000, EOPData_IAU1980, 
      rECItoECI, DCM, TOD, GCRF, ITRF, rECItoECI, PEF, satsv, EOPData_IAU2000A
using StaticArrays: SVector, SArray, @SVector
using HTTP: get
using Roots: find_zeros

include("constants.jl")
include("observations/process_radar.jl")
include("propagation/propagation.jl")
include("postprocessing/least_squares.jl")

end