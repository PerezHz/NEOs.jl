module NEOs

# __precompile__(false)
import Base: RefValue, isless, show, string, in, zero, iszero, isnan, summary
import PlanetaryEphemeris as PE
import SatelliteToolboxTransformations: sv_ecef_to_eci, sv_ecef_to_ecef, ecef_to_geocentric
import Tables: Schema, istable, rowaccess, rows, schema

using AutoHashEquals, Dates, HTTP, InteractiveUtils, JLD2, JSON, LazyArtifacts,
      LinearAlgebra, Printf, SatelliteToolboxTransformations, Scratch, SPICE,
      TaylorIntegration, TaylorSeries, XML

using AstroAngles: hms2rad, rad2hms, dms2rad, rad2dms
using DataFrames: AbstractDataFrame, DataFrame, DataFrameRow, nrow, eachrow, eachcol,
      groupby, combine
using Dates: epochms2datetime
using DelimitedFiles: readdlm
using Distributions: Chisq, cdf
using Healpix: Resolution, ang2pixRing
using HORIZONS: smb_spk
using LinearAlgebra: inv!
using LsqFit: curve_fit, vcov, stderror
using OhMyThreads: tmap, tmap!, @allow_boxed_captures
using Parameters: @with_kw, @unpack
using PlanetaryEphemeris: TaylorInterpCallingArgs, TaylorInterpolant, au, su, ea, yr, RE,
      Rx, Ry, Rz, R_sun, α_p_sun, δ_p_sun, daysec, auday2kmsec, kmsec2auday, c_au_per_day,
      c_au_per_sec, c_cm_per_sec, semimajoraxis, eccentricity, inclination, longascnode,
      argperi, timeperipass, meanmotion, meananomaly, nbodyind, numberofbodies, selecteph,
      getinterpindex, pole_rotation, t2c_jpl_de430
using Roots: find_zeros
using StaticArraysCore: SVector
using StatsBase: mean, std
using TaylorIntegration: VectorCache, RetAlloc, init_cache, taylorinteg!

# Common
export Parameters
export d_EM_km, d_EM_au
export julian2etsecs, etsecs2julian, dtutc2et, et2dtutc, dtutc2jdtdb, jdtdb2dtutc,
       et_to_200X, days_to_200X, dtutc_to_200X, dtutc2days, days2dtutc, rad2arcsec,
       arcsec2rad, mas2rad, range2delay, rangerate2doppler
export loadjpleph, sunposvel, earthposvel, moonposvel, apophisposvel197, apophisposvel199
# Minor bodies astrometry interface
export MPC, NEOCP, NEOCC, NEODyS2, JPL
export UniformWeights, SourceWeights, Veres17
export ZeroDebiasing, SourceDebiasing, Farnocchia15, Eggl20
export numberofdays, unpacknum, packnum, unpackdesig, packdesig
export date, measure, observatory, rms, debias, ra, dec, mag, catalogue, frequency
export obsposECEF, obsposvelECI
export update_catalogues_mpc, search_catalogue_code, search_catalogue_value
export update_observatories_mpc, search_observatory_code, fetch_observatory_information
export fetch_optical_mpc80, read_optical_mpc80, write_optical_mpc80
export fetch_neocp_objects, read_neocp_objects, write_neocp_objects
export fetch_optical_rwo, read_optical_rwo, write_optical_rwo
export fetch_optical_ades, read_optical_ades, write_optical_ades
export nobs, astrometry, datediff, reduce_tracklets
export wra, wdec, dra, ddec, isoutlier, unfold, compute_radec, residuals
export fetch_radar_jpl, read_radar_jpl, write_radar_jpl
export fetch_radar_rwo, read_radar_rwo, write_radar_rwo
export compute_delay, radar_astrometry
# Propagation

#=
# Asteroid dynamical models
export RNp1BP_pN_A_J23E_J2S_ng_eph_threads!, RNp1BP_pN_A_J23E_J2S_eph_threads!, newtonian!
# Propagation
export propagate, propagate_lyap, propagate_root
# PE and NEOs ephemerides
export loadpeeph, bwdfwdeph
# Osculating
export pv2kep, yarkp2adot
# Least squares
export LeastSquaresCache, Newton, DifferentialCorrections, LevenbergMarquardt,
       project, chi2, nms, nrms, leastsquares, leastsquares!, tryls, outlier_rejection!
# AbstractOrbit
export GaussOrbit, MMOVOrbit, LeastSquaresOrbit, epoch, minmaxdates, critical_value,
       sigmas, snr, jplcompare, keplerian, uncertaintyparameter
# Orbit determination
export ODProblem, mmov, gaussmethod, tsaiod, gaussiod, issinglearc,
       initialorbitdetermination, orbitdetermination
# B plane
export valsecchi_circle, bopik, mtp
# Magnitude
export absolutemagnitude
=#

include("constants.jl")
include("units.jl")
include("jpleph.jl")
include("parameters.jl")
include("astrometry/astrometry.jl")
# include("propagation/propagation.jl")
# include("orbitdetermination/orbitdetermination.jl")
# include("postprocessing/bplane.jl")
include("init.jl")
include("deprecated.jl")

end