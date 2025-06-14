module NEOs

# __precompile__(false)
import Base: show, isless, string


import Base: convert, zero, iszero, isnan, in, summary
import Tables: istable, rowaccess, rows, schema, Schema
import SatelliteToolboxTransformations: sv_ecef_to_eci, sv_ecef_to_ecef, ecef_to_geocentric
import JLD2: writeas
import PlanetaryEphemeris as PE

using AutoHashEquals, Dates, HTTP, JSON, Printf, Scratch, XML

using AstroAngles: hms2rad, rad2hms, dms2rad, rad2dms
using StaticArraysCore: SVector


using InteractiveUtils, LazyArtifacts, LinearAlgebra,
      TaylorSeries, SatelliteToolboxTransformations, TaylorIntegration,
      SPICE, JLD2
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
# Minor bodies astrometry interface
export MPC, NEOCP, NEOCC, NEODyS2, JPL
export unpacknum, packnum, unpackdesig, packdesig, date, ra, dec, observatory
export update_catalogues_mpc, search_catalogue_code, search_catalogue_value
export update_observatories_mpc, search_observatory_code, fetch_observatory_information
export fetch_optical_mpc80, read_optical_mpc80, write_optical_mpc80
export fetch_neocp_objects, read_neocp_objects, write_neocp_objects
export fetch_optical_rwo, read_optical_rwo, write_optical_rwo
export fetch_optical_ades, read_optical_ades, write_optical_ades
export fetch_radar_jpl, read_radar_jpl, write_radar_jpl
export fetch_radar_rwo, read_radar_rwo, write_radar_rwo

#=
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
export unfold, wra, wdec, μra, μdec, isoutlier, compute_radec, residuals
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
       project, chi2, nms, nrms, leastsquares, leastsquares!, tryls, outlier_rejection!
# AbstractOrbit
export GaussOrbit, MMOVOrbit, LeastSquaresOrbit, epoch, minmaxdates, critical_value,
       sigmas, snr, jplcompare, keplerian, uncertaintyparameter
# Orbit determination
export ODProblem, Parameters, mmov, gaussmethod, tsaiod, gaussiod, curvature, issinglearc,
       initialorbitdetermination, orbitdetermination
# B plane
export valsecchi_circle, bopik, mtp
# Magnitude
export absolutemagnitude
=#

include("constants.jl")
include("parameters.jl")
include("astrometry/astrometry.jl")
# include("observations/observations.jl")
# include("propagation/propagation.jl")
# include("orbitdetermination/orbitdetermination.jl")
# include("postprocessing/bplane.jl")
include("init.jl")
# include("deprecated.jl")

end