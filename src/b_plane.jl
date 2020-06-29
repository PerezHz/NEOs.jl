# computes Valsecchi circle associated to a mean motion resonance, following Valsecchi et al (2003)
# Returns radius R (au) and zeta-axis coordinate D (au)
#a: asteroid heliocentric semimajor axis (au)
#e: asteroid heliocentric eccentricity
#i: asteroid heliocentric inclination, ecliptic (rad)
#k/h; h heliocentric revolutions of asteroid per k heliocentric revolutions of Earth
# m_pl: planet mass normalized to Sun's mass, equal to Earth mass in solar masses by default
function valsecchi_circle(a, e, i, k, h; m_pl=3.003489614915764e-6)
    # U_x = sqrt( 2 - (1/a) - a*(1-(e^2)) ) # TODO: CHECK SIGN
    U_y = sqrt( a*(1-(e^2)) )*cos(i) - 1
    # U_z = sqrt( a*(1-(e^2)) )*sin(i) # TODO: CHECK SIGN
    U_norm = sqrt( 3 - (1/a) - 2*sqrt(a*(1-(e^2)))*cos(i) )
    # expression below should be equal to asteroid heliocentric elliptic semimamajor axis in au units
    #@show 1/(1-U_^2-2U_y)
    return valsecchi_circle(U_y, U_norm, k, h, m_pl=m_pl)
end

# computes Valsecchi circle associated to a mean motion resonance, following Valsecchi et al (2003)
# Returns radius R (au) and zeta-axis coordinate D (au)
# U_y: Y-component of unperturbed planetocentric velocity (Y-axis coincides with the direction of motion of the planet)
# U_norm: Euclidean norm of unperturbed planetocentric velocity
# both U_y, U_norm are in units such that the heliocentric velocity of the planet is 1
#k/h; h heliocentric revolutions of asteroid per k heliocentric revolutions of Earth
# m_pl: planet mass normalized to Sun's mass, equal to Earth mass in solar masses by default
# a_pl: planetary heliocentric semimajor axis in au; default value is 1
function valsecchi_circle(U_y, U_norm, k, h; m_pl=3.003489614915764e-6, a_pl=1.0)
    a0p = a_pl*(k/h)^(2/3)
    cosθ = U_y/U_norm
    # @show cosθ, U_norm, U_y
    sinθ = sin(acos(cosθ)) #sqrt(1-cosθ^2) #  # TODO: CHECK SIGN
    cosθ0p = (1-(U_norm^2)-(1/a0p))/(2U_norm)
    sinθ0p = sin(acos(cosθ0p))
    # c = m/U^2 (Valsecchi et al, 2003, Sec. 2.3, first sentence below Eq. 3)
    c = m_pl/(U_norm^2)
    R0 = abs( c*sinθ0p/(cosθ0p-cosθ) )
    D0 = c*sinθ/(cosθ0p-cosθ)
    return R0, D0
end
