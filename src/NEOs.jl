module NEOs

# __precompile__(false)

import Base: show, string, hash, isequal, ==, isless, convert, zero, iszero, isnan, in
import Tables: istable, rowaccess, rows, schema, Schema
import SatelliteToolboxTransformations: sv_ecef_to_eci, sv_ecef_to_ecef, ecef_to_geocentric
import JLD2: writeas
import PlanetaryEphemeris as PE

using Dates, InteractiveUtils, LazyArtifacts, LinearAlgebra, Printf
using Downloads: download

using AutoHashEquals, JSON, TaylorSeries, SatelliteToolboxTransformations,
        TaylorIntegration, SPICE, JLD2, Scratch
using DelimitedFiles: readdlm
using HTTP: get
using DataFrames: DataFrame, nrow, eachcol, eachrow, groupby, combine, AbstractDataFrame,
        DataFrameRow, GroupedDataFrame
using PlanetaryEphemeris: TaylorInterpCallingArgs, TaylorInterpolant, daysec, yr,
        auday2kmsec, su, ea, au, c_au_per_day, α_p_sun, δ_p_sun, pole_rotation,
        c_cm_per_sec, c_au_per_sec, t2c_jpl_de430, R_sun, RE, Rx, Ry, Rz, semimajoraxis,
        eccentricity, inclination, longascnode, argperi, timeperipass, nbodyind,
        numberofbodies, kmsec2auday, meanmotion, meananomaly, selecteph, getinterpindex
using Healpix: ang2pixRing, Resolution
using StatsBase: mean, std
using LsqFit: curve_fit, vcov
using Roots: find_zeros
using Clustering: kmeans
using HORIZONS: smb_spk
using OhMyThreads: tmap, tmap!
using TaylorIntegration: RetAlloc, _determine_parsing!, _taylorinteg!

# Constants
export d_EM_km, d_EM_au
# CatalogueMPC
export unknowncat, read_catalogues_mpc, write_catalogues_mpc, update_catalogues_mpc,
       search_cat_code
# ObservatoryMPC
export unknownobs, hascoord, read_observatories_mpc, write_observatories_mpc,
       update_observatories_mpc, search_obs_code
# RadecMPC
export ra, dec, date, observatory, read_radec_mpc, write_radec_mpc, get_radec_mpc,
       fetch_radec_mpc
# RadarJPL
export hasdelay, hasdoppler, delay, doppler, rcvr, xmit, read_radar_jpl, write_radar_jpl
# Units
export julian2etsecs, etsecs2julian, datetime2et, et_to_200X, days_to_200X, datetime_to_200X,
       datetime2days, days2datetime, rad2arcsec, arcsec2rad, mas2rad
# JPL ephemerides
export loadjpleph, sunposvel, earthposvel, moonposvel, apophisposvel197, apophisposvel199
# PE and NEOs ephemerides
export loadpeeph, bwdfwdeph
# Observer position in ECEF and ECI frames
export obsposECEF, obsposvelECI
# Optical astrometry processing
export compute_radec, select_debiasing_table, debiasing, w8sveres17, residuals, unfold,
       relax_factor, outlier
# Radar astrometry processing
export compute_delay, radar_astrometry
# Asteroid dynamical models
export RNp1BP_pN_A_J23E_J2S_ng_eph_threads!, RNp1BP_pN_A_J23E_J2S_eph_threads!, newtonian!
# Propagation
export NEOParameters, propagate, propagate_lyap, propagate_root
# B plane
export valsecchi_circle, bopik
# Least squares
export project, chi2, nms, nrms, diffcorr, newtonls, levenbergmarquardt, tryls, sigmas, snr
# Osculating
export pv2kep, yarkp2adot
# Too Short Arc
export tooshortarc
# Gauss method
export gauss_method, gaussinitcond
# Outlier rejection
export outlier_rejection
# Orbit determination
export curvature, jplcompare, uncertaintyparameter, issinglearc, isgauss, orbitdetermination

include("constants.jl")
include("observations/observations.jl")
include("propagation/propagation.jl")
include("orbit_determination/orbitdetermination.jl")
include("postprocessing/outlier_rejection.jl")
include("init.jl")

end