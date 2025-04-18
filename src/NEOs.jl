module NEOs

# __precompile__(false)

import Base: show, string, isless, convert, zero, iszero, isnan, in
import Tables: istable, rowaccess, rows, schema, Schema
import SatelliteToolboxTransformations: sv_ecef_to_eci, sv_ecef_to_ecef, ecef_to_geocentric
import JLD2: writeas
import PlanetaryEphemeris as PE

using Dates, InteractiveUtils, LazyArtifacts, LinearAlgebra, Printf, JSON,
      TaylorSeries, SatelliteToolboxTransformations, TaylorIntegration,
      SPICE, JLD2, Scratch, HTTP
using AutoHashEquals.Compat
using Base: RefValue
using Dates: epochms2datetime
using Distributions: Chisq, cdf
using Downloads: download
using DelimitedFiles: readdlm
using DataFrames: DataFrame, nrow, eachcol, eachrow, groupby, combine, AbstractDataFrame,
        DataFrameRow, GroupedDataFrame
using PlanetaryEphemeris: TaylorInterpCallingArgs, TaylorInterpolant, daysec, yr,
        auday2kmsec, su, ea, au, c_au_per_day, α_p_sun, δ_p_sun, pole_rotation,
        c_cm_per_sec, c_au_per_sec, t2c_jpl_de430, R_sun, RE, Rx, Ry, Rz, semimajoraxis,
        eccentricity, inclination, longascnode, argperi, timeperipass, nbodyind,
        numberofbodies, kmsec2auday, meanmotion, meananomaly, selecteph, getinterpindex
using Healpix: ang2pixRing, Resolution
using StatsBase: mean, std
using LinearAlgebra: inv!
using LsqFit: curve_fit, vcov, stderror
using Roots: find_zeros
using HORIZONS: smb_spk
using OhMyThreads: tmap, tmap!, @allow_boxed_captures
using TaylorIntegration: VectorCache, RetAlloc, init_cache, taylorinteg!
using Parameters: @with_kw, @unpack

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
# NEOCPObject
export fetch_objects_neocp, get_radec_neocp, fetch_radec_neocp, get_orbits_neocp
# Units
export julian2etsecs, etsecs2julian, dtutc2et, dtutc2jdtdb, et2dtutc, jdtdb2dtutc,
       et_to_200X, days_to_200X, dtutc_to_200X, dtutc2days, days2dtutc,
       rad2arcsec, arcsec2rad, mas2rad
# JPL ephemerides
export loadjpleph, sunposvel, earthposvel, moonposvel, apophisposvel197, apophisposvel199
# PE and NEOs ephemerides
export loadpeeph, bwdfwdeph
# Observer position in ECEF and ECI frames
export obsposECEF, obsposvelECI
# Tracklet
export nobs, astrometry, datediff
# Error model
export UniformWeights, Veres17, ADESWeights, NEOCCWeights, NEODyS2Weights
export Farnocchia15, Eggl20, ZeroDebiasing, NEOCCDebiasing, NEODyS2Debiasing
# Optical astrometry processing
export compute_radec, unfold, residuals, isoutlier
# Radar astrometry processing
export compute_delay, radar_astrometry
# Asteroid dynamical models
export RNp1BP_pN_A_J23E_J2S_ng_eph_threads!, RNp1BP_pN_A_J23E_J2S_eph_threads!, newtonian!
# Propagation
export propagate, propagate_lyap, propagate_root
# Osculating
export pv2kep, yarkp2adot
# Least squares
export LeastSquaresCache, Newton, DifferentialCorrections, LevenbergMarquardt,
       project, chi2, nms, nrms, leastsquares, leastsquares!, tryls
# AbstractOrbit
export PreliminaryOrbit, LeastSquaresOrbit, epoch, minmaxdates, critical_value, sigmas,
       snr, jplcompare, uncertaintyparameter
# Too Short Arc
export tsaiod
# Gauss method
export gauss_method, gaussiod
# Orbit determination
export Parameters, ODProblem, curvature, issinglearc, orbitdetermination
# B plane
export valsecchi_circle, bopik, mtp
# Outlier rejection
export outlier_rejection!
# Magnitude
export absolutemagnitude

include("constants.jl")
include("parameters.jl")
include("observations/observations.jl")
include("propagation/propagation.jl")
include("orbitdetermination/orbitdetermination.jl")
include("postprocessing/bplane.jl")
include("init.jl")
include("deprecated.jl")

end