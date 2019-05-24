# JPL-DE-430/431 model for the orientation of the Earth

# From JPL DE430 documentation (Folkner et al., 2014):
# Only the long-term change of the Earth's orientation is modeled in the ephemeris
# integration. The Earth orientation model used for the DE 430 and 431 integration
# is based on the International Astronomical Union (IAU) 1976 precession model\* with
# an estimated linear correction and on a modified IAU 1980 nutation model \* including
# only terms with a period of 18.6 years.

# The mean longitude of the ascending node of the lunar orbit measured on the ecliptic
# plane from the mean equinox of date is calculated by

Ω(t) = deg2rad( (125+2/60+40.280/3600)-(1934+8/60+10.539/3600)*(t/36525)+(7.455/3600)*(t/36525)^2+(0.008/3600)*(t/36525)^3 )

# where `t` is the TDB time in Julian days from J2000.0. The nutations in longitude
# $\Delta \psi$ and obliquity $\Delta \epsilon$ are given by

Delta_psi(t) = deg2rad( (-17.206262/3600)*sin(Ω(t)) - 0.0004*sin(deg2rad(200.9+1.97129t)) )
Delta_epsilon(t) = deg2rad( (9.205348/3600)*cos(Ω(t)) + 0.0002*cos(deg2rad(200.9+1.97129t)) )

# The true pole of date unit vector $\vec p_\mathrm{d}$ is computed by rotating the
# Earth-fixed pole vector by the effect of the 18.6-year nutation term to give

function pole_date(t)
    # ϵ0 = ϵ̄(t)
    # Δϵ = Delta_epsilon(t)
    # Δψ = Delta_psi(t)
    # ϵ = ϵ0 + Δϵ
    # pdx = sin(Δψ)*sin(ϵ)
    # pdy = cos(Δψ)*sin(ϵ)*cos(ϵ0) - cos(ϵ)*sin(ϵ0)
    # pdz = cos(Δψ)*sin(ϵ)*sin(ϵ0) + cos(ϵ)*cos(ϵ0)
    Nut_mat = nutation_iau80(t)
    # inv(Nut_mat)*[0,0,1] = Nut_mat[3,:]
    return Nut_mat[3,:]
end

# where the mean obliquity $\bar \epsilon$ is given by

ϵ̄(t) = deg2rad( (84381.448/3600)-(46.815/3600)*(t/36525)-(0.00059/3600)*(t/36525)^2+(0.001813/3600)*(t/36525)^3 )

# The pole unit vector in the intertial frame $\vec p_\mathrm{E}$ is computed by
# precessing the pole of date with an estimated linear correction,

function pole_frame(t)
    p_d = pole_date(t)
    # inverse of precession matrix with linear corrections to precession angles
    pm = (Rz(Zeta(t))*Ry(-Theta(t))*Rz(zeta(t))*Rx(-phi_x(t))*Ry(-phi_y(t)))
    return pm*p_d

end

# where

phi_x(t) = deg2rad( (phi_x0 + 100*(t/36525)*Dt_phi_x)/3600 )
phi_y(t) = deg2rad( (phi_y0 + 100*(t/36525)*Dt_phi_y)/3600 )

phi_x0 = 5.6754203322893470E-03 #x-axis rotation at J2000.0 (arcseconds)
phi_y0 = -1.7022656914989530E-02 #y-axis rotation at J2000.0 (arcseconds)
Dt_phi_x = 2.7689915574483550E-04 #Negative obliquity rate correction (arcseconds/year)
Dt_phi_y = -1.2118591216559240E-03 #Precession rate correction times sine of obliquity (arcseconds/year)

# are estimated linear corrections with offsets and rates given by `phi_x0`, `phi_y0`,
# `Dt_phi_x` and `Dt_phi_y` and the precession angles are given by

# Zeta = t -> deg2rad( (2306.2181/3600)*(t/36525)+(0.30188/3600)*(t/36525)^2+(0.017998/3600)*(t/36525)^3 )
# Theta = t -> deg2rad( (2004.3109/3600)*(t/36525)-(0.42665/3600)*(t/36525)^2-(0.041833/3600)*(t/36525)^3 )
# zeta = t -> deg2rad( (2306.2181/3600)*(t/36525)+(1.09468/3600)*(t/36525)^2+(0.018203/3600)*(t/36525)^3 )

function Zeta(t)
    t_cy = t/36525
    Zeta_arcsec = 0.017998*t_cy
    Zeta_arcsec = (0.30188 + Zeta_arcsec)*t_cy
    Zeta_arcsec = (2306.2181 + Zeta_arcsec)*t_cy
    return deg2rad(Zeta_arcsec/3600)
end

function Theta(t)
    t_cy = t/36525
    Theta_arcsec = -0.041833*t_cy
    Theta_arcsec = (-0.42665 + Theta_arcsec)*t_cy
    Theta_arcsec = (2004.3109 + Theta_arcsec)*t_cy
    return deg2rad(Theta_arcsec/3600)
end

function zeta(t)
    t_cy = t/36525
    zeta_arcsec = 0.018203*t_cy
    zeta_arcsec = (1.09468 + zeta_arcsec)*t_cy
    zeta_arcsec = (2306.2181 + zeta_arcsec)*t_cy
    return deg2rad(zeta_arcsec/3600)
end

# The rotation matrices are defined by

function Rx(alpha::T) where {T<:Number}
    res = Array{T}(undef, 3, 3)
    res[1, 1] = one(alpha)
    res[2, 1] = zero(alpha)
    res[3, 1] = zero(alpha)
    res[1, 2] = zero(alpha)
    res[2, 2] = cos(alpha)
    res[3, 2] = -sin(alpha)
    res[1, 3] = zero(alpha)
    res[2, 3] = sin(alpha)
    res[3, 3] = cos(alpha)
    return res
end

# Lieske et al. (1977, 1979) introduces the convention that this matrix has
# opposite signs outside the sine terms (wrt Rx, Rz)
function Ry(alpha::T) where {T<:Number}
    res = Array{T}(undef, 3, 3)
    res[1, 1] = cos(alpha)
    res[2, 1] = zero(alpha)
    res[3, 1] = sin(alpha)
    res[1, 2] = zero(alpha)
    res[2, 2] = one(alpha)
    res[3, 2] = zero(alpha)
    res[1, 3] = -sin(alpha)
    res[2, 3] = zero(alpha)
    res[3, 3] = cos(alpha)
    return res
end

function Rz(alpha::T) where {T<:Number}
    res = Array{T}(undef, 3, 3)
    res[1, 1] = cos(alpha)
    res[2, 1] = -sin(alpha)
    res[3, 1] = zero(alpha)
    res[1, 2] = sin(alpha)
    res[2, 2] = cos(alpha)
    res[3, 2] = zero(alpha)
    res[1, 3] = zero(alpha)
    res[2, 3] = zero(alpha)
    res[3, 3] = one(alpha)
    return res
end

# myatan(x, y) = y>=zero(x)?( x>=zero(x)?atan(y/x):(atan(y/x)+pi) ):( x>=zero(x)?(atan(y/x)+2pi):(atan(y/x)+pi) )
# myatan2(x, y) = y>=zero(x)?( x>=zero(x)?atan(y/x):(atan(y/x)-pi) ):( x>=zero(x)?(atan(y/x)):(atan(y/x)+pi) )

function pole_ra(t)
    pole_frame_t = pole_frame(t)
    return atan( pole_frame_t[2]/pole_frame_t[1] ) + π
end

function pole_dec(t)
    pole_frame_t = pole_frame(t)
    return atan(  pole_frame_t[3]/sqrt(pole_frame_t[1]^2+pole_frame_t[2]^2)  )
end

function pole_radec(t)
    pole_frame_t = pole_frame(t)
    pole_ra = atan( pole_frame_t[2]/pole_frame_t[1] ) + π
    pole_dec = atan(  pole_frame_t[3]/sqrt(pole_frame_t[1]^2+pole_frame_t[2]^2)  )
    return pole_ra, pole_dec
end

# rotation from inertial frame to frame with pole at right ascension α and declination δ
function pole_rotation(α::T, δ::T) where {T <: Number}
    m = Matrix{T}(undef, 3, 3)
    # diagonal terms
    m[1,1] = sin(α)^2 + sin(δ)*(cos(α)^2)
    m[2,2] = cos(α)^2 + sin(δ)*(sin(α)^2)
    m[3,3] = sin(δ)
    # off-diagonal terms
    m[1,2] = cos(α)*sin(α)*(-1+sin(δ))
    m[3,1] = -cos(δ)*cos(α)
    m[3,2] = -cos(δ)*sin(α)
    m[2,1] = m[1,2]
    m[1,3] = -m[3,1]
    m[2,3] = -m[3,2]
    return m
end

# rotation matrix from inertial frame to Earth pole at time t (days) since J2000.0
function earth_pole_rotation(t)
    # α_ep = Earth pole at time t (Julian days) since J2000.0
    # δ_ep = Earth pole at time t (Julian days) since J2000.0
    α_ep, δ_ep = pole_radec(t)
    return pole_rotation(α_ep, δ_ep)
end

# IAU 1980 nutation matrix (ESAA 1992)
# Δψ (rad): nutation in longitude
# ϵ (rad): mean obliquity
# ϵ0 (rad): true obliquity of date
function nutation_iau80(t)
    ϵ0 = ϵ̄(t)
    Δϵ = Delta_epsilon(t)
    Δψ = Delta_psi(t)
    ϵ = ϵ0 + Δϵ
    return Rx(-ϵ)*Rz(-Δψ)*Rx(ϵ0)
end

# terrestrial-to-celestial vector transformation, JPL DE 430/431 Earth orientation model
function t2c_jpl_de430(t)
    inv_P_iau7680 = Rz(Zeta(t))*Ry(-Theta(t))*Rz(zeta(t))
    corrections = Rx(-phi_x(t))*Ry(-phi_y(t))
    inv_N_iau80 = inv(nutation_iau80(t))
    return inv_P_iau7680*corrections*inv_N_iau80
end

# celestial-to-terrestrial vector transformation, JPL DE 430/431 Earth orientation model
function c2t_jpl_de430(t)
    P_iau7680 = Rz(-zeta(t))*Ry(Theta(t))*Rz(-Zeta(t))
    corrections = Ry(phi_y(t))*Rx(phi_x(t))
    N_iau80 = nutation_iau80(t)
    return N_iau80*corrections*P_iau7680
end

# The following was taken from "Report of the IAU/IAG Working Group", Seidelmann et. al, 2006
#
# Recommended values for the direction of the north pole of rotation and the
# prime meridian of the satellites
#
# d = interval in days from the standard epoch (J2000.0)
# T = interval in Julian centuries (of 36,525 days) from the standard epoch

function WGCCRE2006_moon_E1(d)
    return deg2rad(125.045-0.0529921d)
end

function WGCCRE2006_moon_E2(d)
    return deg2rad(250.089-0.1059842d)
end

function WGCCRE2006_moon_E3(d)
    return deg2rad(260.008+13.0120009d)
end

function WGCCRE2006_moon_E4(d)
    return deg2rad(176.625+13.3407154d)
end

function WGCCRE2006_moon_E5(d)
    return deg2rad(357.529+0.9856003d)
end

function WGCCRE2006_moon_E6(d)
    return deg2rad(311.589+26.4057084d)
end

function WGCCRE2006_moon_E7(d)
    return deg2rad(134.963+13.0649930d)
end

function WGCCRE2006_moon_E8(d)
    return deg2rad(276.617+0.3287146d)
end

function WGCCRE2006_moon_E9(d)
    return deg2rad(34.226+1.7484877d)
end

function WGCCRE2006_moon_E10(d)
    return deg2rad(15.134-0.1589763d)
end

function WGCCRE2006_moon_E11(d)
    return deg2rad(119.743+0.0036096d)
end

function WGCCRE2006_moon_E12(d)
    return deg2rad(239.961+0.1643573d)
end

function WGCCRE2006_moon_E13(d)
    return deg2rad(25.053+12.9590088d)
end

function moon_pole_ra(d)

    ans = 269.9949+0.0031d/36525
    ans += -3.8787sin(WGCCRE2006_moon_E1(d))
    ans += -0.1204sin(WGCCRE2006_moon_E2(d))
    ans +=  0.0700sin(WGCCRE2006_moon_E3(d))
    ans += -0.0172sin(WGCCRE2006_moon_E4(d))
    ans +=  0.0072sin(WGCCRE2006_moon_E6(d))
    ans += -0.0052sin(WGCCRE2006_moon_E10(d))
    ans +=  0.0043sin(WGCCRE2006_moon_E13(d))
    return deg2rad(ans)
end

function moon_pole_dec(d)

    ans = 66.5392+0.013d/36525
    ans +=  1.5419cos(WGCCRE2006_moon_E1(d))
    ans +=  0.0239cos(WGCCRE2006_moon_E2(d))
    ans += -0.0278cos(WGCCRE2006_moon_E3(d))
    ans +=  0.0068cos(WGCCRE2006_moon_E4(d))
    ans += -0.0029cos(WGCCRE2006_moon_E6(d))
    ans +=  0.0009cos(WGCCRE2006_moon_E7(d))
    ans +=  0.0008cos(WGCCRE2006_moon_E10(d))
    ans += -0.0009cos(WGCCRE2006_moon_E13(d))

    return deg2rad(ans)
end

function moon_pole_w(d)

    ans = 38.3213+13.17635815d-1.4E-12d^2
    ans += 3.5610sin(WGCCRE2006_moon_E1(d))
    ans += 0.1208sin(WGCCRE2006_moon_E2(d))
    ans += -0.0642sin(WGCCRE2006_moon_E3(d))
    ans += +0.0158sin(WGCCRE2006_moon_E4(d))
    ans += +0.0252sin(WGCCRE2006_moon_E5(d))
    ans += -0.0066sin(WGCCRE2006_moon_E6(d))
    ans += -0.0047sin(WGCCRE2006_moon_E7(d))
    ans += -0.0046sin(WGCCRE2006_moon_E8(d))
    ans += +0.0028sin(WGCCRE2006_moon_E9(d))
    ans += +0.0052sin(WGCCRE2006_moon_E10(d))
    ans += +0.0040sin(WGCCRE2006_moon_E11(d))
    ans += +0.0019sin(WGCCRE2006_moon_E12(d))
    ans += -0.0044sin(WGCCRE2006_moon_E13(d))

    return deg2rad(ans)
end
