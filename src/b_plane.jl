# Planetary impact cross section; "critical" B
# derived from conservation of energy and angular momentum
# in hyperbolic close encounter
#  = what impact parameter B corresponds to a grazing impact,
# given planet radius and mass, and asymptotic vel v_infty
# if actual B is equal or less to this, then impact happens
# μ_P: planetary gravitational parameter (au^3/day^2)
# R_P: planetary radius (au)
# vinf: asymptotic inbound velocity (au/day)
# output is in planet radii
function crosssection(μ_P, R_P, vinf)
    return sqrt( 1 + (2μ_P)/(R_P*(vinf)^2) )
end

# computes Öpik's coordinates of impact parameter vector B in hyperbolic planetary close encounter
# following Farnocchia et al (2019)
# xae: asteroid's geocentric position/velocity vector at closest approach in au,au/day
# xes: planet's heliocentric position/velocity vector at asteroid's closest approach in au,au/day
function bopik(xae, xes)
    μ_E = PlanetaryEphemeris.μ[ea]
    # asteroid geocentric range at 2029 closest approach (au)
    rae = sqrt(xae[1]^2 + xae[2]^2 + xae[3]^2)
    # osculating semimajor axis at closest approach (negative since hyperbolic)
    a = semimajoraxis(xae..., μ_E, 0.0)
    # asymptotic inbound velocity v_\infty (au/day)
    v_infty = sqrt(μ_E/(-a))
    # h = r × v
    hvec = cross(xae[1:3], xae[4:6])
    h = sqrt(hvec[1]^2 + hvec[2]^2 + hvec[3]^2)
    # Laplace-Runge-Lenz (eccentricity) vector
    # \vec{e} = (\vec{v} × \vec{h})/μ - \vec{r}/r
    evec = cross(xae[4:6], hvec)/μ_E - xae[1:3]/rae
    # osculating eccentricity
    e = eccentricity(xae..., μ_E, 0.0)
    # osculating inclination
    i = inclination(xae...)
    # periapsis position (unit vector)
    P_v = evec./e
    # periapsis velocity (unit vector)
    Q_v = cross(hvec, P_v)./h
    #inbound asymptote direction
    S_v = (P_v + (sqrt(e^2 - 1))Q_v)/e
    # B-vector: "vector from the planet center to the intersection between the B-plane and the asymptote"
    Bvec = cross(S_v, hvec)./v_infty
    # Earth impact cross section ("critical B")
    b_E = crosssection(μ_E, RE/au, v_infty)
    # @show b_E
    # Earth's heliocentric velocity at CA
    v_pl = xes[4:6]
    #ξ-axis is essentially the MOID (Valsecchi et al, 2003)
    ξ_v_unnormalized = cross(v_pl, S_v)
    ξ_v_norm = sqrt(ξ_v_unnormalized[1]^2 + ξ_v_unnormalized[2]^2 + ξ_v_unnormalized[3]^2)
    ξ_v = ξ_v_unnormalized./ξ_v_norm
    # @show ξ_v
    #ζ-axis: delay/advance in CA time (Valsecchi et al, 2003)
    ζ_v = -cross(S_v, ξ_v)
    B_dot_ξ = dot(Bvec, ξ_v)
    B_dot_ζ = dot(Bvec, ζ_v)

    # computation of U vector
    # res = xes[1:3] # Earth's heliocentric position, au
    # res_norm = norm(res) # Earth's heliocentric range, au
    # res_unit = res/res_norm # Earth's heliocentric radius unit vector
    # @show res_norm
    ves = v_pl*au/daysec # Earth's velocity, km/s
    ves_norm = norm(ves) # Earth's speed, km/s
    ves_unit = ves/ves_norm # Earth's velocity unit vector
    # X-Y angle (degrees)
    # @show rad2deg(acos(dot(ves_unit, res_unit)))
    # angle between Y-axis and \vec{U}
    cosθ = dot(S_v, ves_unit)
    # @show cosθ()
    v_infty_kms = v_infty*au/daysec # asteroid unperturbed speed, km/sec
    # @show v_infty_kms()
    # The norm of \vec{U} in appropriate units
    U_unit = (2pi*au)/(yr*daysec) # 1 U in km/s
    U_norm = v_infty_kms/U_unit
    # U_y
    U_y = U_norm*cosθ
    # @show U_y U_norm

    return (ξ=B_dot_ξ*au/RE, ζ=B_dot_ζ*au/RE, U=(y=U_y, norm=U_norm), b=b_E)
end

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
