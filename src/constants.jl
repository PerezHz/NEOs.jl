# Internal paths

# Path to NEOs src directory
const src_path = dirname(pathof(NEOs))
# Path to scratch space
const scratch_path = Ref{String}("")

# Load Solar System, accelerations, newtonian potentials and TT-TDB 2000-2100 ephemeris
const sseph_artifact_path = joinpath(artifact"sseph_p100", "sseph343ast016_p100y_et.jld2")
const sseph::TaylorInterpolant{Float64, Float64, 2} = JLD2.load(sseph_artifact_path, "ss16ast_eph")
const acceph::TaylorInterpolant{Float64, Float64, 2} = JLD2.load(sseph_artifact_path, "acc_eph")
const poteph::TaylorInterpolant{Float64, Float64, 2} = JLD2.load(sseph_artifact_path, "pot_eph")
const ttmtdb::TaylorInterpolant{Float64, Float64, 1} = TaylorInterpolant(sseph.t0, sseph.t, sseph.x[:,end])

# Earth orientation parameters (eop) 2000
const eop_IAU2000A::EopIau2000A = fetch_iers_eop(Val(:IAU2000A))

# Parsing

# Regular expression to parse a catalogue in MPC format
const CATALOGUE_MPC_REGEX = r"\s{2}(?P<code>\w{1})\s{4}(?P<name>.*)"
# Header of MPC catalogues file
const CATALOGUES_MPC_HEADER = "Char   Catalogue"
# Regular expression to parse an observatory in MPC format
const OBSERVATORY_MPC_REGEX = Regex(string(
    # Code regex + space (columns 1-3)
    raw"(?P<code>[A-Z\d]{3})",
    # Longitude regex (columns 4-13)
    raw"(?P<long>[\.\d\s]{10})",
    # Cosine regex + space (column 14-21)
    raw"(?P<cos>[\.\d\s]{8})",
    # Sine regex (column 22-30)
    raw"(?P<sin>[\+\-\.\d\s]{9})",
    # Name regex (columns 31-80)
    raw"(?P<name>.*)",
))
# Header of MPC observatories file 
const OBSERVATORIES_MPC_HEADER = "Code  Long.   cos      sin    Name"
# Regular expression to parse an optical measurement in MPC format
const RADEC_MPC_REGEX = Regex(string(
    # Number regex (columns 1-5)
    raw"(?P<num>.{5})",
    # Temporary designation regex (columns 6-12)
    raw"(?P<tmpdesig>.{7})",
    # Discovery asterisk regex (column 13)
    raw"(?P<discovery>[\*\s]{1})",
    # Publishable note regex (column 14)
    raw"(?P<publishnote>.{1})",
    # Observation technique regex (column 15)
    raw"(?P<obstech>[^xX]{1})",
    # Date of observation regex (columns 16-32)
    raw"(?P<date>\d{4}\s\d{2}\s\d{2}\.[\d\s]{6})",
    # Right ascension regex (columns 33-44)
    raw"(?P<α>\d{2}\s\d{2}\s\d{2}\.[\d\s]{3})",
    # Declination regex (columns 45-56)
    raw"(?P<δ>[\+|\-]{1}\d{2}\s\d{2}\s\d{2}\.[\d\s]{2})",
    # Info 1 regex (columns 57-65)
    raw"(?P<info1>.{9})",
    # Magnitude regex (columns 66-70)
    raw"(?P<mag>[\.\s\d]{5})",
    # Band regex (column 71)
    raw"(?P<band>[\w\s]{1})",
    # Catalogue regex (column 72)
    raw"(?P<catalogue>[\w\s]{1})",
    # Info 2 regex (columns 73-77)
    raw"(?P<info2>.{5})",
    # Observatory code regex (columns 78-80)
    raw"(?P<obscode>\w{3})",
    # Optional fields (in case of satellite observations)
    # Breakline regex
    raw"(?:\n)?",
    # Number regex (columns 1-5)
    raw"(?<optional>(?P<_num_>.{5})?",
    # Temporary designation regex (columns 6-12)
    raw"(?P<_tmpdesig_>.{7})?",
    # Blank space regex (column 13)
    raw"(?P<_discovery_>\s)?",
    # Publishable note regex (column 14)
    raw"(?P<_publishnote_>.{1})?",
    # s regex (column 15)
    raw"(?P<_obstech_>s)?",
    # Date of observation regex (columns 16-32)
    raw"(?P<_date_>\d{4}\s\d{2}\s\d{2}\.[\d\s]{6})",
    # Units + space regex (columns 33-34)
    raw"(?P<_units_>\d\s)",
    # X component of geocentric vector (columns 35-46)
    raw"(?P<_x_>[\-\+]{1}[\.\d\s]{11})",
    # Y component of geocentric vector (columns 47-58)
    raw"(?P<_y_>[\-\+]{1}[\.\d\s]{11})",
    # Z component of geocentric vector (columns 59-70)
    raw"(?P<_z_>[\-\+]{1}[\.\d\s]{11})",
    # Band regex (column 71)
    raw"(?P<_band_>[\w\s]{1})?",
    # Catalogue regex (column 72)
    raw"(?P<_catalogue_>[\w\s]{1})?",
    # Info 2 regex (columns 73-77)
    raw"(?P<_info2_>.{5})?",
    # Observatory code regex (columns 78-80)
    raw"(?P<_obscode_>\w{3})?)?",
))
# Regular expression to parse a radar measurement in JPL format
const RADAR_JPL_REGEX = Regex(string(
    # ID regex + tab
    raw"(?P<id>.*)\t",
    # Date regex + tab
    raw"(?P<date>.*)\t",
    # Measurement regex + tab
    raw"(?P<measurement>.*)\t",
    # Uncertainty regex + tab
    raw"(?P<uncertainty>.*)\t",
    # Units regex + tab
    raw"(?P<units>.*)\t",
    # Frequency regex + tab
    raw"(?P<freq>.*)\t",
    # Reciever regex + tab
    raw"(?P<rcvr>.*)\t",
    # Emitter regex + tab
    raw"(?P<xmit>.*)\t",
    # Bouncepoint regex + end of line
    raw"(?P<bouncepoint>.*)"
))
# Format of date in JPL radar data files
const RADAR_JPL_DATEFORMAT = "yyyy-mm-dd HH:MM:SS"

# MPC catalogues corresponding to debiasing tables included in https://doi.org/10.1016/j.icarus.2014.07.033
const mpc_catalogue_codes_2014 = ["a", "b", "c", "d", "e", "g", "i", "j", "l", "m", "o", "p", "q", "r",
                                  "u", "v", "w", "L", "N"]

# MPC catalogues corresponding to debiasing tables included in https://doi.org/10.1016/j.icarus.2019.113596
const mpc_catalogue_codes_2018 = ["a", "b", "c", "d", "e", "g", "i", "j", "l", "m", "n", "o", "p", "q",
                                  "r", "t", "u", "v", "w", "L", "N", "Q", "R", "S", "U", "W"]

# URLs

# MPC catalogues file url 
const CATALOGUES_MPC_URL = "https://minorplanetcenter.net/iau/info/CatalogueCodes.html"
# MPC observatories file url 
const OBSERVATORIES_MPC_URL = "https://minorplanetcenter.net/iau/lists/ObsCodes.html"

# MPC database search url 
const search_mpc_url = "https://www.minorplanetcenter.net/db_search/show_object?utf8=%E2%9C%93&object_id="
# MPC observations url
const obs_mpc_url = "https://www.minorplanetcenter.net/tmp2/"

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
const μ_S = PE.GMS

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
