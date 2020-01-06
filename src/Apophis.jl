module Apophis

__precompile__(false)

export propagate, observer_position, apophisdofs,
    ssdofs, delay_doppler, ismonostatic,
    mas2rad, t2c_rotation_iau_00_06, process_radar_data_jpl, RadarDataJPL,
    julian2etsecs, etsecs2julian,
    RNp1BP_pN_A_J23E_J2S_ng_eph!, RNp1BP_pN_A_J23E_J2S_ng_eph_threads!,
    propagate_distributed, parallel_run,
    utcepoch, delay, delay_sigma, delay_units, doppler, doppler_sigma,
    doppler_units, freq, rcvr, xmit, bouncepoint

using Distributed
using TaylorIntegration
using Printf, DelimitedFiles, Test, LinearAlgebra
using Dates: DateTime, julian2datetime, datetime2julian
# import PlanetaryEphemeris
using PlanetaryEphemeris: daysec, su, ea, Λ2, Λ3, α_p_sun, δ_p_sun, moon_pole_ra,
    moon_pole_dec, t2c_jpl_de430, pole_rotation, au, J2000, c_au_per_day, R_sun,
    c_cm_per_sec, c_au_per_sec, yr
using JLD
using AstroTime, EarthOrientation, SOFA, SPICE
using Dates

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
# const μ = vcat(PlanetaryEphemeris.μ[1:27], 0.0)
# const μ = vcat(PlanetaryEphemeris.μ[1:10], 0.0, PlanetaryEphemeris.μ[12:27], 0.0) # set Pluto's mass = 0
# BODYN_GM values correspond to GM values retrieved from https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de431.tpc
BODY1_GM   = 2.2031780000000021E+04
BODY2_GM   = 3.2485859200000006E+05
BODY399_GM = 3.9860043543609598E+05
BODY301_GM = 4.9028000661637961E+03
BODY4_GM   = 4.2828375214000022E+04
BODY5_GM   = 1.2671276480000021E+08
BODY6_GM   = 3.7940585200000003E+07
BODY7_GM   = 5.7945486000000080E+06
BODY8_GM   = 6.8365271005800236E+06
BODY9_GM   = 0.0 #9.7700000000000068E+02 # set Pluto's mass = 0
BODY10_GM  = 1.3271244004193938E+11

BODY2000001_GM = 6.2809393000000000E+01
BODY2000002_GM = 1.3923011000000001E+01
BODY2000003_GM = 1.6224149999999999E+00
BODY2000004_GM = 1.7288008999999999E+01
BODY2000010_GM = 5.5423920000000004E+00
BODY2000015_GM = 2.0981550000000002E+00
BODY2000016_GM = 1.5300480000000001E+00
BODY2000031_GM = 2.8448720000000001E+00
BODY2000048_GM = 1.1351590000000000E+00
BODY2000052_GM = 1.1108039999999999E+00
BODY2000065_GM = 1.4264810000000001E+00
BODY2000087_GM = 9.8635300000000004E-01
BODY2000088_GM = 1.1557990000000000E+00
BODY2000451_GM = 1.0295259999999999E+00
BODY2000511_GM = 2.3312860000000000E+00
BODY2000704_GM = 2.3573170000000001E+00

# conversion of GM from km^3/s^2 to au^3/day^2
mu_units = (daysec^2)/(au^3)

μ = [BODY10_GM, BODY1_GM, BODY2_GM, BODY399_GM, BODY301_GM, BODY4_GM, BODY5_GM,
    BODY6_GM, BODY7_GM, BODY8_GM, BODY9_GM,
    BODY2000001_GM, BODY2000004_GM, BODY2000002_GM, BODY2000010_GM,
    BODY2000031_GM, BODY2000704_GM, BODY2000511_GM, BODY2000015_GM,
    BODY2000003_GM, BODY2000016_GM, BODY2000065_GM, BODY2000088_GM,
    BODY2000048_GM, BODY2000052_GM, BODY2000451_GM, BODY2000087_GM
]*mu_units

const N = length(μ)

# Matrix of J2 interactions included in DE430 ephemeris, according to Folkner et al., 2014
const UJ_interaction = fill(false, N)
# UJ_interaction[su] = true
UJ_interaction[ea] = true
# UJ_interaction[mo] = true

const j2_body_index = findall(x->x, UJ_interaction)

const apophisdofs = union(3N-2:3N, 6N-2:6N)
const ssdofs = setdiff(1:6N, apophisdofs)

# standard value of nominal mean angular velocity of Earth (rad/day), ESAA 2014 Sec 7.4.3.3 p. 296
const ω = daysec*7.292115e-5 #7.2921151467e-5
# The relationship of the angular velocity of the earth Omega with LOD is (https://www.iers.org/IERS/EN/Science/EarthRotation/UT1LOD.html)
# where `omega` is in units of rad/day, and `lod` is in units of milliseconds
omega(lod) = (1e-12daysec)*(72921151.467064 - 0.843994809lod)

const A_sun = 1.06e8 # Solar corona parameter A [cm^-3] (ESAA 2014, Table 8.5 p. 329)
const a_sun = 4.89e5 # Solar corona parameter a [cm^-3] (ESAA 2014, Table 8.5 p. 329)
const b_sun = 3.91e5 # Solar corona parameter b [cm^-3] (ESAA 2014, Table 8.5 p. 329)

const S0_sun = 63.15E6 # Sun radiated power intensity at photosphere surface, Watt/meter^2
# const m2_s3_to_au2_day3 = 1e-6daysec^3/au^2 # conversion factor from m^2/sec^3 to au^2/day^3

include("process_radar_data_jpl.jl")
include("topocentric.jl")
include("delay_doppler.jl")
include("asteroid_dynamical_models.jl")
include("initial_conditions.jl")
include("integration_methods.jl")
include("propagation.jl")

function __init__()
    @show length(methods(TaylorIntegration.jetcoeffs!))
    @show methods(RNp1BP_pN_A_J23E_J2S_ng_eph!)
    @show methods(RNp1BP_pN_A_J23E_J2S_ng_eph_threads!)
    @show methods(TaylorIntegration.jetcoeffs!)
    # load JPL ephemerides
    # loadjpleph()
end

end