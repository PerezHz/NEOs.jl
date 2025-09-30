## Internal paths

# Path to NEOs src directory
const SRC_PATH = dirname(pathof(NEOs))
# Path to scratch space
const SCRATCH_PATH = Ref{String}("")

## Minor bodies astrometry interface

# Characters for MPC base 62 encoding
const BASE_62_ENCODING = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"

# CatalogueMPC

# See https://github.com/IAU-ADES/ADES-Master/blob/83c7232d5a019c36244363cf3c7bcdeb3f10af6e/Python/bin/packUtil.py#L184
const CATALOGUE_MPC_CODES_TO_NAMES = Dict(
    ' ' => "UNK",
    'a' => "USNOA1",
    'b' => "USNOSA1",
    'c' => "USNOA2",
    'd' => "USNOSA2",
    'e' => "UCAC1",
    'f' => "Tyc1",
    'g' => "Tyc2",
    'h' => "GSC1.0",
    'i' => "GSC1.1",
    'j' => "GSC1.2",
    'k' => "GSC2.2",
    'l' => "ACT",
    'm' => "GSCACT",
    'n' => "SDSS8",
    'o' => "USNOB1",
    'p' => "PPM",
    'q' => "UCAC4",
    'r' => "UCAC2",
    's' => "USNOB2",  # USNOB2 missing on ADES web page
    't' => "PPMXL",
    'u' => "UCAC3",
    'v' => "NOMAD",
    'w' => "CMC14",
    'x' => "Hip2",
    'y' => "Hip1",
    'z' => "GSC",
    'A' => "AC",
    'B' => "SAO1984",
    'C' => "SAO",
    'D' => "AGK3",
    'E' => "FK4",
    'F' => "ACRS",
    'G' => "LickGas",
    'H' => "Ida93",
    'I' => "Perth70",
    'J' => "COSMOS",
    'K' => "Yale",
    'L' => "2MASS",
    'M' => "GSC2.3",
    'N' => "SDSS7",
    'O' => "SSTRC1",
    'P' => "MPOSC3",
    'Q' => "CMC15",
    'R' => "SSTRC4",
    'S' => "URAT1",
    'T' => "URAT2",  # URAT2 missing on ADES web page
    'U' => "Gaia1",
    'V' => "Gaia2",
    'W' => "Gaia3",
    'X' => "Gaia3E",
    'Y' => "UCAC5",
    'Z' => "ATLAS2",
    '0' => "IHW",
    '1' => "PS1_DR1",
    '2' => "PS1_DR2",
    '3' => "Gaia_Int",
    '4' => "GZ",
    '5' => "UBSC",
    '6' => "Gaia2016",
)

const CATALOGUE_MPC_NAMES_TO_CODES = Dict(value => key for (key, value) in CATALOGUE_MPC_CODES_TO_NAMES)

const CATALOGUES_MPC_FILE_URL = "https://www.minorplanetcenter.net/iau/info/astCat_photCat.json"

# ObservatoryMPC

const OBSERVATORIES_MPC_API = "https://data.minorplanetcenter.net/api/obscodes"

# OpticalMPC80

const MPC80_OPTICAL_COLUMNS = [
    1:5, 6:12, 13:13, 14:14, 15:15, 16:32, 33:44,
    45:56, 57:65, 66:70, 71:71, 72:72, 73:77, 78:80
]
const OBSERVATIONS_MPC_API = "https://data.minorplanetcenter.net/api/get-obs"

const OBSERVATIONS_NEOCP_API = "https://data.minorplanetcenter.net/api/get-obs-neocp"

# NEOCPObject

const NEOCP_OBJECT_COLUMNS = [
    1:7, 9:11, 13:25, 26:34, 35:42, 44:47,
    49:69, 79:82, 83:89, 90:94, 95:101
]

const NEOCP_OBJECTS_FILE_URL = "https://www.minorplanetcenter.net/iau/NEO/neocp.txt"

# OpticalRWO

const RWO_OPTICAL_HEADER = """
! Object   Obser ============= Date ============= ================== Right Ascension =================  ================= Declination ===================== ==== Magnitude ==== Ast Obs  Residual SEL
! Design   K T N YYYY MM DD.dddddddddd   Accuracy HH MM SS.sss  Accuracy      RMS  F     Bias    Resid sDD MM SS.ss  Accuracy      RMS  F     Bias    Resid Val  B   RMS  Resid Cat Cod       Chi A M
"""

const RWO_OPTICAL_COLUMNS = [
    1:11, 12:12, 14:14, 16:16, 18:38, 41:49, 51:62, 65:73, 74:82, 84:84, 88:93,
    97:102, 104:115, 118:126, 127:135, 137:137, 141:146, 147:155, 157:161, 162:162,
    163:170, 171: 176, 177:180, 181:183, 184:193, 195:195, 197:197
]

const NEOCC_URL = "https://neo.ssa.esa.int/"
const OBSERVATIONS_NEOCC_API = NEOCC_URL * "PSDB-portlet/download?file="

const NEODyS2_URL = "https://newton.spacedys.com/neodys/"
const OBSERVATIONS_NEODyS2_API = "https://newton.spacedys.com/~neodys2/mpcobs/"

# TimeOfDay

const EOPIAU = Union{EopIau1980, EopIau2000A}

# Earth orientation parameters (EOP) 2000
const EOP_IAU2000A::EopIau2000A = fetch_iers_eop(Val(:IAU2000A))

# RadarJPL

const JPL_TO_MPC_OBSCODES = Dict(
    "-1"  => "251", # Arecibo
    "-2"  => "254",	# Haystack, Westford
    "-9"  => "256",	# Green Bank
    "-13" => "252",	# Goldstone DSS 13, Fort Irwin
    "-14" => "253",	# Goldstone DSS 14, Fort Irwin
    "-35" => "263", # Canberra DSS 35
    "-36" => "264", # Canberra DSS 36
    "-38" => "255",	# Yevpatoriya
    "-43" => "265", # Canberra DSS 43
    "-47" => "271", # ATCA DSS 47
    "-73" => "259",	# EISCAT Tromso UHF
)

const MPC_TO_JPL_OBSCODES = Dict(value => key for (key, value) in JPL_TO_MPC_OBSCODES)

const RADAR_JPL_DATEFORMAT = "yyyy-mm-dd HH:MM:SS"

const RADAR_JPL_API = "https://ssd-api.jpl.nasa.gov/sb_radar.api"

# RadarRWO

const RWO_RADAR_HEADER = """
! Object   Obser ====== Date =======  ============ Radar range/range rate (km or km/d) =============  Station  ====  Residual
! Design   K T N YYYY MM DD hh:mm:ss         Measure  Accuracy       rms F        Bias       Resid    TRX  RCX        Chi   S
"""

const RWO_RADAR_COLUMNS = [
    1:11,  12:12, 14:14, 16:16, 18:36, 37:52, 53:62, 63:72,
    74:74, 75:86, 87:98, 103:105, 108:110, 111:122, 125:125
]

# AbstractDebiasingScheme

# MPC catalogues corresponding to debiasing tables included in:
# - https://doi.org/10.1016/j.icarus.2014.07.033
const CATALOGUE_MPC_CODES_2014 = [
    'a', 'b', 'c', 'd', 'e', 'g', 'i', 'j', 'l', 'm',
    'o', 'p', 'q', 'r', 'u', 'v', 'w', 'L', 'N'
]

# MPC catalogues corresponding to debiasing tables included in:
# - https://doi.org/10.1016/j.icarus.2019.113596
const CATALOGUE_MPC_CODES_2018 = [
    'a', 'b', 'c', 'd', 'e', 'g', 'i', 'j', 'l', 'm', 'n', 'o', 'p',
    'q', 'r', 't', 'u', 'v', 'w', 'L', 'N', 'Q', 'R', 'S', 'U', 'W'
]

# Propagation

# Abbreviation for a dense TaylorInterpolant
const DensePropagation1{T, U} = TaylorInterpolant{T, U, 1, Vector{T}, Vector{Taylor1{U}}}
const DensePropagation2{T, U} = TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}}

# Load Solar System, accelerations, newtonian potentials and TT-TDB 2000-2100 ephemeris
const sseph_artifact_path = joinpath(artifact"sseph_p100", "sseph343ast016_p100y_et.jld2")
const sseph::DensePropagation2{Float64, Float64} = JLD2.load(sseph_artifact_path, "ss16ast_eph")
const acceph::DensePropagation2{Float64, Float64} = JLD2.load(sseph_artifact_path, "acc_eph")
const poteph::DensePropagation2{Float64, Float64} = JLD2.load(sseph_artifact_path, "pot_eph")
const ttmtdb::DensePropagation1{Float64, Float64} = TaylorInterpolant(sseph.t0, sseph.t, sseph.x[:,end])
const SSEPHORDER::Int = get_order(sseph.x[1])

# Modified julian date offset
const MJD2000 = JD_J2000 - 2400000.5

# Milliseconds between rounding epoch and J2000
const EPOCHMSJ2000::Int = (DateTime(2000, 1, 1, 12) - DateTime(0)).value

# Abbreviations
const cte = constant_term

# Integration parameters
const order = 30
const abstol = 1.0E-30

# Vector of GM's (DE430 values) [au^2 / day^3]
const μ_DE430 = PE.μ
const μ_B16_DE430 = μ_DE430[12:27]     # DE430 GM's of 16 most massive asteroids
const μ_ast343_DE430 = μ_DE430[12:end] # DE430 GM's of 343 main belt asteroids included in DE430 integration
# Gravitational parameter of the Sun [au^2 / day^3]
const μ_S::Float64 = PE.GMS

# Standard value of nominal mean angular velocity of Earth (rad/sec)
# See Explanatory Supplement to the Astronomical Almanac 2014 Sec 7.4.3.3 p. 296
const ω = 7.2921151467e-5 # 7.292115e-5 rad/sec

# Solar corona parameters
# See Table 8.5 in page 329 of Explanatory Supplement to the Astronomical Almanac 2014
const A_sun = 1.06e8  # [cm^-3]
const a_sun = 4.89e5  # [cm^-3]
const b_sun = 3.91e5  # [cm^-3]

# Sun radiated power intensity at photosphere surface, Watt/meter^2
const S0_sun = 63.15E6
# Conversion factor from m^2/sec^3 to au^2/day^3
const m2_s3_to_au2_day3 = 1e-6daysec^3/au^2

# Vector of J_2*R^2 values
# J_2: second zonal harmonic coefficient
# R: radius of the body
const Λ2 = zeros(11)
Λ2[ea] = 1.9679542578489185e-12  # Earth
# Vector of J_3*R^3 values
# J_3: third zonal harmonic coefficient
# R: radius of the body
const Λ3 = zeros(11)
Λ3[ea] = -1.962633335678878e-19  # Earth

# Speed of light
const clightkms = 2.99792458E5   # km/sec
const c_km_per_day = c_au_per_day * au # km/day
const c_km_per_us = c_cm_per_sec * 1e-11 # km/microsecond
# Parameters related to speed of light, c
const c_p2 = 29979.063823897606      # c^2 = 29979.063823897606 au^2/d^2
const c_m2 = 3.3356611996764786e-5   # c^-2 = 3.3356611996764786e-5 d^2/au^2

# Earth-Moon distance in [km]
const d_EM_km = 384_400
# Earth-Moon distance in [au]
const d_EM_au = 384_400 / au

# Zeroth order obliquity of the ecliptic in degrees
# See equation (5-153) in page 5-61 of https://doi.org/10.1002/0471728470.
const ϵ0_deg = 84381.448/3_600

# Gauss gravitational constant
const k_gauss = 0.017_202_098_95
# Earth's sphere of influence radius [AU]
const R_SI = 0.010044
# Earth's physical radius [AU]
const R_EA = 4.24e-5
# Ratio between the mass of the Earth and the mass of the Sun
const μ_ES = PE.μ[ea] / PE.μ[su] # 1 / 328_900.5614

# Postprocessing

# Conversion to V band used by MPC
# See https://minorplanetcenter.net/iau/info/BandConversion.txt
const V_BAND_CORRECTION = Dict{Char, Float64}(
    ' ' => -0.8,
    'U' => -1.3,
    'B' => -0.8,
    'g' => -0.35,
    'V' =>  0,
    'r' =>  0.14,
    'R' =>  0.4,
    'C' =>  0.4,
    'W' =>  0.4,
    'i' =>  0.32,
    'z' =>  0.26,
    'I' =>  0.8,
    'J' =>  1.2,
    'w' => -0.13,
    'y' =>  0.32,
    'L' =>  0.2,
    'H' =>  1.4,
    'K' =>  1.7,
    'Y' =>  0.7,
    'G' =>  0.28,
    'v' =>  0,
    'c' => -0.05,
    'o' =>  0.33,
    'u' => +2.5
)

# Parameters of the linear H and G magnitude system for asteroids
# See https://minorplanetcenter.net/iau/ECS/MPCArchive/1985/MPC_19851227.pdf
const SLOPE_PARAMETER = 0.15
const PHASE_INTEGRAL_A1 = 3.33
const PHASE_INTEGRAL_A2 = 1.87
const PHASE_INTEGRAL_B1 = 0.63
const PHASE_INTEGRAL_B2 = 1.22

# Earth escape velocity
# See https://doi.org/10.1006/icar.2002.6910
const EARTH_ESCAPE_VELOCITY = 11.18 # km/s