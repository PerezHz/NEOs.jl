# Internal paths 

# Path to NEOs src directory 
const src_path = dirname(pathof(NEOs))
# Path to scratch space 
const scratch_path = Ref{String}("")

# URLs

# MPC observatories file url 
const mpc_observatories_url = "https://minorplanetcenter.net/iau/lists/ObsCodes.html"
# MPC database search url 
const search_mpc_url = "https://www.minorplanetcenter.net/db_search/show_object?utf8=%E2%9C%93&object_id="
# MPC observations url 
const obs_mpc_url = "https://www.minorplanetcenter.net/tmp/"

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