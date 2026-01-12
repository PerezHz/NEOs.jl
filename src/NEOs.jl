module NEOs

# __precompile__(false)
import Base: RefValue, isless, show, string, getindex, in, zero, iszero, isnan, summary,
       firstindex, lastindex, merge
import PlanetaryEphemeris as PE
import PlanetaryEphemeris: kmsec2auday, semimajoraxis, eccentricity, inclination, argperi,
       longascnode, meanmotion, meananomaly, timeperipass, eccentricanomaly, trueanomaly,
       t2c_jpl_de430
import SatelliteToolboxTransformations: sv_ecef_to_eci, sv_ecef_to_ecef, ecef_to_geocentric
import StatsBase: weights
import Tables: Schema, istable, rowaccess, rows, schema
import TaylorSeries: get_order

using AngleBetweenVectors, AutoHashEquals, Dates, HTTP, InteractiveUtils, JLD2, JSON,
      LazyArtifacts, LinearAlgebra, Printf, Scratch, SPICE, TaylorIntegration,
      TaylorSeries, XML

using AstroAngles: hms2rad, rad2hms, dms2rad, rad2dms
using DataFrames: AbstractDataFrame, DataFrame, DataFrameRow, nrow, eachrow, eachcol,
      groupby, combine
using Dates: epochms2datetime
using DelimitedFiles: readdlm
using Distributions: Chisq, Normal, Uniform, cdf, quantile
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
using QuadGK: quadgk
using Roots: Bisection, find_zero, find_zeros
using SatelliteToolboxTransformations: DCM, EARTH_ANGULAR_SPEED, EopIau1980, EopIau2000A,
      GCRF, ITRF, J2000, JD_J2000, OrbitStateVector, PEF, TIRS, T_ECIs, T_ECIs_IAU_2006,
      fetch_iers_eop, geodetic_to_ecef, get_Δat, r_ecef_to_ecef, r_ecef_to_eci,
      sv_eci_to_ecef
using SpecialFunctions: erf
using StaticArraysCore: SVector, MVector, SMatrix, MMatrix
using StatsBase: mean, std
using TaylorIntegration: VectorCache, RetAlloc, init_cache, taylorinteg!, update_cache!,
      taylorstep!, set_psol!, findroot!

# Common
export Parameters
export d_EM_km, d_EM_au, MJD2000
export julian2etsecs, etsecs2julian, dtutc2et, et2dtutc, dtutc2jdtdb, jdtdb2dtutc,
       et_to_200X, days_to_200X, dtutc_to_200X, dtutc2days, days2dtutc, rad2arcsec,
       arcsec2rad, mas2rad, range2delay, rangerate2doppler, chi2, nms, nrms
export loadjpleph, sunposvel, earthposvel, moonposvel, apophisposvel197, apophisposvel199,
       tt_tdb, dtt_tdb
# Minor bodies astrometry interface
export MPC, NEOCP, NEOCC, NEODyS2, JPL
export UniformWeights, SourceWeights, Veres17
export ZeroDebiasing, SourceDebiasing, Farnocchia15, Eggl20
export numberofdays, unpacknum, packnum, unpackdesig, packdesig
export date, measure, observatory, rms, debias, corr, ra, dec, mag, band, catalogue,
       frequency, residual, weight, weights, isoutlier, nout, notout, notoutobs
export obsposECEF, obsposvelECI
export update_catalogues_mpc, search_catalogue_code, search_catalogue_value
export update_observatories_mpc, search_observatory_code, fetch_observatory_information
export fetch_optical_mpc80, read_optical_mpc80, write_optical_mpc80
export fetch_neocp_objects, read_neocp_objects, write_neocp_objects
export fetch_optical_rwo, read_optical_rwo, write_optical_rwo
export fetch_optical_ades, read_optical_ades, write_optical_ades
export nobs, datediff, reduce_tracklets
export wra, wdec, dra, ddec, unfold, compute_radec, residuals
export fetch_radar_jpl, read_radar_jpl, write_radar_jpl
export fetch_radar_rwo, read_radar_rwo, write_radar_rwo
export compute_delay, radar_astrometry
# Propagation
export nongravs!, gravityonly!, newtonian!
export rvelea, scaled_variables, propagate, propagate_lyap, propagate_root
# Orbit determination
export ODProblem, LeastSquaresCache, Newton, DifferentialCorrections, LevenbergMarquardt,
       GaussOrbit, MMOVOrbit, LeastSquaresOrbit, AdmissibleRegion
export gm, frame, elements, iscircular, iselliptic, isparabolic, ishyperbolic, conicsection,
       yarkp2adot, cartesian2keplerian, cartesian2equinoctial, cartesian2attributable,
       keplerian2cartesian, equinoctial2cartesian, attributable2cartesian,
       keplerian2equinoctial, equinoctial2keplerian, pericenter
export curvature
export bwdfwdeph, propres, propres!
export leastsquares, leastsquares!, tryls, outlier_rejection!, project, critical_value
export variables, epoch, noptical, nradar, minmaxdates, optical, sigmas, snr, keplerian,
       equinoctial, attributable, uncertaintyparameter, absolutemagnitude, diameter, mass
export topo2bary, bary2topo, attr2bary, tsaiod
export mmov, gaussmethod, gaussiod, jtls, issinglearc, initialorbitdetermination,
       orbitdetermination
# Impact monitoring
export ImpactTarget, IMProblem, BPlane, MTP, bopik, mtp, targetplane, crosssection,
       valsecchi_circle
export LineOfVariations, VirtualAsteroid, CloseApproach, Return, lineofvariations,
       virtualasteroids, closeapproaches, showersnreturns, sigma, lbound, ubound
# export VirtualImpactor, virtualimpactors, impact_probability, impactor_table, impactenergy,
#       palermoscale, torinoscale, semiwidth, stretching

include("constants.jl")
include("units.jl")
include("jpleph.jl")
include("parameters.jl")
include("fasttaylors.jl")
include("astrometry/astrometry.jl")
include("propagation/propagation.jl")
include("orbitdetermination/orbitdetermination.jl")
include("impactmonitoring/impactmonitoring.jl")
include("init.jl")
include("deprecated.jl")

end