module NEOs

# __precompile__(false)

if !isdefined(Base, :get_extension)
    using Requires
end

import Base: hash, ==, show, isless, isnan, convert
import PlanetaryEphemeris as PE
import JLD2: writeas

using Distributed, JLD2, TaylorIntegration, Printf, DelimitedFiles, Test, LinearAlgebra,
      Dates, EarthOrientation, SPICE, Quadmath, LazyArtifacts, TaylorSeries,
      InteractiveUtils, AutoHashEquals, RemoteFiles 
using PlanetaryEphemeris: daysec, su, ea, α_p_sun, δ_p_sun, t2c_jpl_de430, pole_rotation,
      au, c_au_per_day, R_sun, c_cm_per_sec, c_au_per_sec, yr, RE, TaylorInterpolant, Rx,
      Ry, Rz, semimajoraxis, eccentricity, inclination, longascnode, argperi, timeperipass,
      nbodyind, ordpres_differentiate, numberofbodies, kmsec2auday, auday2kmsec, meanmotion,
      meananomaly
using Healpix: ang2pixRing, Resolution
using SatelliteToolbox: get_iers_eop_iau_2000A, EOPData_IAU1980, EOPData_IAU2000A, JD_J2000,
      orbsv, sv_ecef_to_eci, get_Δat, nutation_fk5, r_ecef_to_eci, T_ECIs, T_ECIs_IAU_2006,
      we, OrbitStateVector, r_ecef_to_eci, DCM
import SatelliteToolbox.sv_ecef_to_eci
using Dates: format
using HTTP: get
using IntervalRootFinding: roots, interval, Interval, mid

# Constants
export d_EM_km, d_EM_au
# CatalogueMPC
export unknowncat, isunknown, read_catalogues_mpc, parse_catalogues_mpc, write_catalogues_mpc, update_catalogues_mpc,
       search_cat_code
# ObservatoryMPC
export unknownobs, hascoord, read_observatories_mpc, parse_observatories_mpc, write_observatories_mpc,
       update_observatories_mpc, search_obs_code
# RadecMPC
export num, tmpdesig, discovery, publishnote, obstech, ra, dec, info1, mag, band, catalogue, info2, observatory,
       read_radec_mpc, parse_radec_mpc, search_circulars_mpc, write_radec_mpc
# RadarJPL
export hasdelay, hasdoppler, ismonostatic, date, delay, delay_sigma, delay_units, doppler, doppler_sigma,
       doppler_units, freq, rcvr, xmit, bouncepoint, read_radar_jpl, write_radar_jpl
# Units
export julian2etsecs, etsecs2julian, datetime2et, et_to_200X, days_to_200X, datetime_to_200X,
       datetime2days, days2datetime, rad2arcsec, arcsec2rad, mas2rad
# JPL Ephemerides
export loadjpleph, sunposvel, earthposvel, moonposvel, apophisposvel197, apophisposvel199,
       loadpeeph, bwdfwdeph
# Osculating
export pv2kep, yarkp2adot
# Topocentric
export obs_pos_ECEF, obsposvelECI, t2c_rotation_iau_76_80
# Process radec
export compute_radec, debiasing, w8sveres17, radec_astrometry, residuals
# Process radar
export compute_delay, radar_astrometry
# Gauss method
export gauss_method
# Asteroid dynamical models
export RNp1BP_pN_A_J23E_J2S_ng_eph!, RNp1BP_pN_A_J23E_J2S_ng_eph_threads!, RNp1BP_pN_A_J23E_J2S_eph_threads!
# Propagate
export propagate, propagate_lyap, propagate_root, save2jldandcheck

export valsecchi_circle, nrms, chi2, newtonls, newtonls_6v, diffcorr, newtonls_Q, bopik

include("constants.jl")
include("observations/process_radar.jl")
include("orbit_determination/gauss_method.jl")
include("propagation/propagation.jl")
include("postprocessing/least_squares.jl")

function __init__()
    @static if !isdefined(Base, :get_extension)
        @require Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c" begin
            @require DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0" include("../ext/DataFramesExt.jl")
        end
    end
end

end