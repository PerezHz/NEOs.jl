module Apophis

__precompile__(false)

export propagate, au, yr, observer_position, apophisdofs,
    ssdofs, c_au_per_day, μ, delay_doppler, ismonostatic,
    mas2rad, t2c_rotation_iau_00_06, process_radar_data_jpl, RadarDataJPL,
    RNp1BP_pN_A_J23E_J2S_ng_eph!, semimajoraxis, eccentricity, inclination,
    julian2etsecs, etsecs2julian

using Printf, DelimitedFiles, Test, LinearAlgebra
using Dates: DateTime, julian2datetime, datetime2julian
using TaylorIntegration, PlanetaryEphemeris
using JLD
using AstroTime, EarthOrientation, SOFA, SPICE

# integration parameters
const varorder = 10
const order = 30
const abstol = 1.0E-30

# vector of G*m values
# Bodies included: Sun, Sun, Mercury, Venus, Earth, Moon, Mars, Jupiter, Saturn, Uranus,
# Neptune, Pluto,
# 1 Ceres, 4 Vesta, 2 Pallas, 10 Hygiea, 31 Euphrosyne, 704 Interamnia,
# 511 Davida, 15 Eunomia, 3 Juno, 16 Psyche, 65 Cybele, 88 Thisbe, 48 Doris,
# 52 Europa, 451 Patientia, 87 Sylvia
# and Apophis as a massless test particle
const μ = vcat(PlanetaryEphemeris.μ[1:27], 0.0)
const N = length(μ)

# Matrix of J2 interactions included in DE430 ephemeris, according to Folkner et al., 2014
const UJ_interaction = fill(false, N)
UJ_interaction[su] .= true
UJ_interaction[ea] .= true
UJ_interaction[mo] .= true

const apophisdofs = union(3N-2:3N, 6N-2:6N)
const ssdofs = setdiff(1:6N, apophisdofs)

# standard value of nominal mean angular velocity of Earth (rad/day), ESAA 2014 Sec 7.4.3.3 p. 296
const ω = daysec*7.292115e-5
# The relationship of the angular velocity of the earth Omega with LOD is (https://www.iers.org/IERS/EN/Science/EarthRotation/UT1LOD.html)
# where `omega` is in units of rad/day, and `lod` is in units of milliseconds
omega(lod) = (1e-12daysec)*(72921151.467064 - 0.843994809lod)

const A_sun = 1.06e8 # Solar corona parameter A [cm^-3] (ESAA 2014, Table 8.5 p. 329)
const a_sun = 4.89e5 # Solar corona parameter a [cm^-3] (ESAA 2014, Table 8.5 p. 329)
const b_sun = 3.91e5 # Solar corona parameter b [cm^-3] (ESAA 2014, Table 8.5 p. 329)

const S0_sun = 63.15E6 # Sun radiated power intensity at photosphere surface, Watt/meter^2
# const m2_s3_to_au2_day3 = 1e-6daysec^3/au^2 # conversion factor from m^2/sec^3 to au^2/day^3

function __init__()
    @show length(methods(RNp1BP_pN_A_J23E_J2S_ng_eph!))
    @show length(methods(TaylorIntegration.jetcoeffs!))
    @show methods(RNp1BP_pN_A_J23E_J2S_ng_eph!)
    # load JPL ephemerides
    loadjpleph()
end

include("process_radar_data_jpl.jl")
include("topocentric.jl")
include("asteroid_dynamical_models.jl")
include("initial_conditions.jl")
include("delay_doppler.jl")
include("osculating.jl")
include("propagation.jl")

end