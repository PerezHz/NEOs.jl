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
export propagate

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

using Distributed
using TaylorIntegration
using Printf, DelimitedFiles, Test, LinearAlgebra
using Dates: DateTime, julian2datetime, datetime2julian
import PlanetaryEphemeris
using PlanetaryEphemeris: daysec, su, ea, α_p_sun, δ_p_sun,
    t2c_jpl_de430, pole_rotation, au, c_au_per_day, R_sun,
    c_cm_per_sec, c_au_per_sec, yr, RE, TaylorInterpolant, Rx, Ry, Rz,
    semimajoraxis, eccentricity, inclination, longascnode, argperi,
    timeperipass, nbodyind
using JLD
using EarthOrientation, SPICE
using Dates
using Quadmath
using Healpix: ang2pixRing, Resolution
using LazyArtifacts
using DataFrames
import SatelliteToolbox
using SatelliteToolbox: nutation_fk5, J2000toGMST, rECEFtoECI,
    get_ΔAT, JD_J2000, EOPData_IAU1980, rECItoECI, DCM,
    TOD, GCRF, ITRF, rECItoECI, PEF, satsv, EOPData_IAU2000A
import RemoteFiles
using StaticArrays: SVector, SArray, @SVector
using TaylorSeries
using InteractiveUtils
using HTTP: get
using Roots: find_zeros

@doc raw"""
    omega(lod)

Returns the angular velocity of the earth in units of rad/sec
```math
\omega = (72921151.467064 - 0.843994809\text{LOD})\times 10^{-12},
```
where LOD is the length of the day in milliseconds. 

See https://www.iers.org/IERS/EN/Science/EarthRotation/UT1LOD.html.
"""
omega(lod) = (1e-12)*(72921151.467064 - 0.843994809lod)

# isless overloads workaround to make rECEFtoECI work with TaylorSeries
import Base.isless
isless(a::Union{Taylor1,TaylorN}, b::Union{Taylor1,TaylorN}) = isless(constant_term(a), constant_term(b))
isless(a::Real, b::Union{Taylor1,TaylorN}) = isless(a, constant_term(b))
isless(a::Union{Taylor1,TaylorN}, b::Real) = isless(constant_term(a), b)

import ReferenceFrameRotations: orthonormalize
# This method extends orthonormalize to handle Taylor1 and TaylorN
# The ReferenceFrameRotations.jl package is licensed under the MIT "Expat" License
# Copyright (c) 2014-2019: Ronan Arraes Jardim Chagas.
# See https://github.com/JuliaSpace/ReferenceFrameRotations.jl
function orthonormalize(dcm::SArray{Tuple{3,3},T,2,9} where {T<:Union{Taylor1,TaylorN}})
    e₁ = dcm[:,1]
    e₂ = dcm[:,2]
    e₃ = dcm[:,3]

    en₁  = e₁/norm(e₁)
    enj₂ =   e₂ - (en₁⋅e₂)*en₁
    en₂  = enj₂ / norm(enj₂)
    enj₃ =   e₃ - (en₁⋅e₃)*en₁
    enj₃ = enj₃ - (en₂⋅enj₃)*en₂
    en₃  = enj₃ / norm(enj₃)

    SArray(  hcat(hcat(en₁, en₂), en₃)  )
end

@doc raw"""
    yarkp2adot(A2, a, e, μ_S)

Returns the average semimajor axis drift due to the Yarkovsky effect
```math
\begin{align*}
    \left\langle\dot{a}\right\rangle & = \frac{2A_2(1-e^2)}{n}\left(\frac{1 \ \text{AU}}{p}\right)^2 \\
    & = \frac{2A_2}{(1-e^2)\sqrt{a\mu_\odot}}(1 \ \text{AU})^2,
\end{align*}
```
where ``A_2`` is the Yarkovsky parameter, ``\mu_\odot = GM_\odot`` is the Sun's mass parameter,
``e`` is the eccentricity, ``n = \sqrt{\mu/a^3}`` is the mean motion, ``p = a(1-e^2)`` is the 
semilatus rectum, and ``a`` is the semimajor axis. 

See https://doi.org/10.1016/j.icarus.2013.02.004.

# Arguments 

- `A2`: Yarkovsky parameter.
- `a`: semimajor axis. 
- `e`: eccentricity. 
- `μ_S`: mass parameter of the Sun. 
"""
function yarkp2adot(A2, a, e, μ_S)
    return 2A2/(sqrt(a)*(1-e^2)*sqrt(μ_S))
end

include("constants.jl")
include("observations/process_radec.jl")
include("propagation/gauss_method.jl")
include("propagation/propagation.jl")
include("process_radar_data_jpl.jl")
include("delay_doppler.jl")
include("asteroid_dynamical_models.jl")
include("initial_conditions.jl")
include("integration_methods.jl")
include("b_plane.jl")
include("least_squares.jl")

end