module Apophis

__precompile__(false)

export propagate, observer_position, apophisdofs,
    ssdofs, delay_doppler, ismonostatic,
    mas2rad, t2c_rotation_iau_00_06, process_radar_data_jpl, RadarDataJPL,
    julian2etsecs, etsecs2julian,
    RNp1BP_pN_A_J23E_J2S_ng_eph!, RNp1BP_pN_A_J23E_J2S_ng_eph_threads!,
    propagate_distributed, parallel_run,
    utcepoch, delay, delay_sigma, delay_units, doppler, doppler_sigma,
    doppler_units, freq, rcvr, xmit, bouncepoint, valsecchi_circle,
    nrms, chi2, newtonls, newtonls_6v

using Distributed
using TaylorIntegration
using Printf, DelimitedFiles, Test, LinearAlgebra
using Dates: DateTime, julian2datetime, datetime2julian
import PlanetaryEphemeris
using PlanetaryEphemeris: daysec, su, ea, α_p_sun, δ_p_sun,
    t2c_jpl_de430, pole_rotation, au, J2000, c_au_per_day, R_sun,
    c_cm_per_sec, c_au_per_sec, yr, RE
using JLD
using EarthOrientation, SOFA, SPICE
using Dates
using Quadmath

# integration parameters
const order = 30
const abstol = 1.0E-30

# vector of G*m values
# Bodies included: Sun, Mercury, Venus, Earth, Moon, Mars, Jupiter, Saturn,
# Uranus, Neptune, Pluto,
# 1 Ceres, 4 Vesta, 2 Pallas, 10 Hygiea, 31 Euphrosyne, 704 Interamnia,
# 511 Davida, 15 Eunomia, 3 Juno, 16 Psyche, 65 Cybele, 88 Thisbe, 48 Doris,
# 52 Europa, 451 Patientia, 87 Sylvia
# and Apophis as a massless test particle
μ_ast = PlanetaryEphemeris.μ[12:27]
# index 11 corresponds to Pluto, its mass is set to zero
# const μ = vcat(PlanetaryEphemeris.μ[1:10], 0.0, μ_ast, 0.0)
const μ = vcat(PlanetaryEphemeris.μ[1:11], μ_ast, 0.0)
const N = length(μ)

# Matrix of J2 interactions included in DE430 ephemeris, according to Folkner et al., 2014
const UJ_interaction = fill(false, N)
# UJ_interaction[su] = true
UJ_interaction[ea] = true
# UJ_interaction[mo] = true

const j2_body_index = findall(x->x, UJ_interaction)

const apophisdofs = union(3N-2:3N, 6N-2:6N)
const ssdofs = setdiff(1:6N, apophisdofs)

# standard value of nominal mean angular velocity of Earth (rad/sec), ESAA 2014 Sec 7.4.3.3 p. 296
const ω = 7.2921151467e-5 # 7.292115e-5 rad/sec
# The relationship of the angular velocity of the earth omega with LOD is (https://www.iers.org/IERS/EN/Science/EarthRotation/UT1LOD.html)
# where `omega` is in units of rad/sec, and `lod` is in units of milliseconds
omega(lod) = (1e-12)*(72921151.467064 - 0.843994809lod)

const A_sun = 1.06e8 # Solar corona parameter A [cm^-3] (ESAA 2014, Table 8.5 p. 329)
const a_sun = 4.89e5 # Solar corona parameter a [cm^-3] (ESAA 2014, Table 8.5 p. 329)
const b_sun = 3.91e5 # Solar corona parameter b [cm^-3] (ESAA 2014, Table 8.5 p. 329)

const S0_sun = 63.15E6 # Sun radiated power intensity at photosphere surface, Watt/meter^2
# const m2_s3_to_au2_day3 = 1e-6daysec^3/au^2 # conversion factor from m^2/sec^3 to au^2/day^3

# vector of J2*R^2 values
const Λ2 = zeros(N)
Λ2[ea] = 1.9679542578489185e-12
# vector of J3*R^3 values
const Λ3 = zeros(N)
Λ3[ea] = -1.962633335678878e-19
const clightkms = 2.99792458E5 # speed of light, km/sec

const jd0 = 2.4547335e6

include("process_radar_data_jpl.jl")
include("topocentric.jl")
include("delay_doppler.jl")
include("radec.jl")
include("asteroid_dynamical_models.jl")
include("initial_conditions.jl")
include("integration_methods.jl")
include("propagation.jl")
include("b_plane.jl")
include("least_squares.jl")

end