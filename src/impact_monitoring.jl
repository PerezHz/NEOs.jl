# computes Valsecchi circle associated to a mean motion resonance, following Valsecchi et al (2003)
# Returns radius R (au) and zeta-axis coordinate D (au)
#a: asteroid heliocentric semimajor axis (au)
#e: asteroid heliocentric eccentricity
#i: asteroid heliocentric inclination, ecliptic (rad)
#k/h; h heliocentric revolutions of asteroid per k heliocentric revolutions of Earth
# m_pl: planet mass normalized to Sun's mass, equal to Earth mass in solar masses by default
function valsecchi_circle(a, e, i, k, h, m_pl=3.003489614915764e-6)
    U_x = sqrt( 2 - (1/a) - a*(1-(e^2)) ) # TODO: CHECK SIGN
    U_y = sqrt( a*(1-(e^2)) )*cos(i) - 1
    U_z = sqrt( a*(1-(e^2)) )*sin(i) # TODO: CHECK SIGN
    U_norm = sqrt( 3 - (1/a) - 2*sqrt(a*(1-(e^2)))*cos(i) )
    # should be equal to asteroid heliocentric elliptic a (Carusi et al, 1990)
    #@show 1/(1-U_^2-2U_y)
    a0p = (k/h)^(2/3)
    cosθ = U_y/U_norm
    sinθ = sin(acos(cosθ)) #sqrt(1-cos_θ_^2) # TODO: CHECK SIGN
    cosθ0p = (1-(U_norm^2)-(1/a0p))/(2U_norm)
    sinθ0p = sin(acos(cosθ0p))
    # c = m/U^2 (Valsecchi et al, 2003, Sec. 2.3, first sentence below Eq. 3)
    c = m_pl/(U_norm^2)
    R0 = abs( c*sinθ0p/(cosθ0p-cosθ) )
    D0 = c*sinθ/(cosθ0p-cosθ)

    return R0, D0
end
