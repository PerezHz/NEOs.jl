@doc raw"""
    crosssection(μ_P, R_P, vinf)

Returns the "critical" ``B``, derived from conservation of energy and angular momentum
```math
B = \sqrt{1 + \frac{2\mu_P}{R_P v_\infty^2}},
```
i.e., what impact parameter ``B`` corresponds to a grazing impact in hyperbolic close 
encounter. ``\mu_P`` is the planet's gravitational parameter, ``R_P`` is the planet's 
radius and ``v_\infty`` is the asymptotic inbound velocity. If actual ``B`` is equal or 
less to this, then impact happens. Output is in planet radii. 

See equations (13)-(14) in pages 4-5 of https://doi.org/10.1007/s10569-019-9914-4.

# Arguments 

- `μ_P`: planetary gravitational parameter (au^3/day^2).
- `R_P`: planetary radius (au).
- `vinf`: asymptotic inbound velocity (au/day).
"""
function crosssection(μ_P, R_P, vinf)
    return sqrt( 1 + (2μ_P)/(R_P*(vinf)^2) )
end

@doc raw"""
    bopik(xae, xes)

Computes Öpik's coordinates of impact parameter vector ``\mathbf{B}`` in hyperbolic planetary
close encounter. Returns a named tuple with the following fields:

- `ξ` = ``\mathbf{B}\cdot\hat{\mathbf{\xi}}/R_E`` and `ζ` = ``\mathbf{B}\cdot\hat{\mathbf{\zeta}}/R_E``, where ``\mathbf{B}`` is the impact parameter vector, ``(\hat{\mathbf{\xi}}, \hat{\mathbf{\zeta}})`` is Öpik's frame and ``R_E`` is the Earth's radius in au. See equations (37)-(38) in page 14 of https://doi.org/10.1007/s10569-019-9914-4.

- `U` is another named tuple with the following fields:
    - `y` = ``U_y``.
    - `norm` = ``||\mathbf{U}||`` where ``\mathbf{U}`` is the planetocentric velocity vector in km/s.

- `b` = ``b_E``, where ``b_E`` is the Earth impact cross section ("critical B"). See [`crosssection`](@ref).

# Arguments 

- `xae`: asteroid's geocentric position/velocity vector at closest approach in au, au/day.
- `xes`: planet's heliocentric position/velocity vector at asteroid's closest approach in au, au/day.
"""
function bopik(xae, xes)

    # Computation of Öpik's frame (\mathbf{ξ}, \mathbf{ζ})
    # See equations (37)-(38) in page 14 of https://doi.org/10.1007/s10569-019-9914-4

    # Earth's gravitational parameter
    μ_E = PE.μ[ea]
    # Asteroid geocentric range at closest approach [au]
    rae = sqrt(xae[1]^2 + xae[2]^2 + xae[3]^2)
    # Osculating semimajor axis at closest approach (negative since hyperbolic)
    a = semimajoraxis(xae..., μ_E, 0.0)
    # Asymptotic inbound velocity v_\infty (au/day)
    # See equation (1) in page 2 of https://doi.org/10.1007/s10569-019-9914-4
    v_infty = sqrt(μ_E/(-a))
    # Angular momentum per unit mass h = r × v
    hvec = cross(xae[1:3], xae[4:6])
    # Magntude of h
    h = sqrt(hvec[1]^2 + hvec[2]^2 + hvec[3]^2)
    # Laplace-Runge-Lenz (eccentricity) vector
    # \vec{e} = (\vec{v} × \vec{h})/μ - \vec{r}/r
    evec = cross(xae[4:6], hvec)/μ_E - xae[1:3]/rae
    # Osculating eccentricity
    e = eccentricity(xae..., μ_E, 0.0)
    # Osculating inclination
    i = inclination(xae...)
    # Periapsis position (unit vector) and periapsis velocity (unit vector)
    # See first sentence below equation (3) in page 2 of https://doi.org/10.1007/s10569-019-9914-4
    P_v = evec./e
    Q_v = cross(hvec, P_v)./h
    # Inbound asymptote direction
    # See equation (2) in page 2 of https://doi.org/10.1007/s10569-019-9914-4
    S_v = (P_v + (sqrt(e^2 - 1))Q_v)/e
    # Earth's heliocentric velocity at closest approach (CA)
    v_pl = xes[4:6]
    # ξ-axis is essentially the MOID (Valsecchi et al, 2003)
    # See equation (37) in page 14 of https://doi.org/10.1007/s10569-019-9914-4
    ξ_v_unnormalized = cross(v_pl, S_v)
    ξ_v_norm = sqrt(ξ_v_unnormalized[1]^2 + ξ_v_unnormalized[2]^2 + ξ_v_unnormalized[3]^2)
    ξ_v = ξ_v_unnormalized./ξ_v_norm
    # ζ-axis: delay/advance in CA time (Valsecchi et al, 2003)
    # See equation (37) in page 14 of https://doi.org/10.1007/s10569-019-9914-4
    ζ_v = -cross(S_v, ξ_v)

    # Computation of Öpik's coordinates of impact parameter vector \mathbf{B}

    # B-vector: "vector from the planet center to the intersection between the B-plane 
    # and the asymptote".
    # See equation (4) in page 3 of https://doi.org/10.1007/s10569-019-9914-4
    Bvec = cross(S_v, hvec)./v_infty
    # Impact parameter vector Öpik's coordinates
    B_dot_ξ = dot(Bvec, ξ_v)
    B_dot_ζ = dot(Bvec, ζ_v)

    # Computation of planetocentric velocity (U) vector
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
    v_infty_kms = v_infty*au/daysec # Asteroid unperturbed speed, km/sec
    # @show v_infty_kms()
    # The norm of \vec{U} in appropriate units
    U_unit = (2pi*au)/(yr*daysec) # 1 U in km/s
    U_norm = v_infty_kms/U_unit
    # U_y
    U_y = U_norm*cosθ
    # @show U_y U_norm

    # Earth impact cross section ("critical B")
    b_E = crosssection(μ_E, RE/au, v_infty)
    # @show b_E

    return (ξ=B_dot_ξ*au/RE, ζ=B_dot_ζ*au/RE, U=(y=U_y, norm=U_norm), b=b_E)
end

@doc raw"""
    valsecchi_circle(a, e, i, k, h; m_pl=3.003489614915764e-6)  

Computes Valsecchi circle associated to a mean motion resonance. Returns radius ``R`` (au) and
``\zeta``-axis coordinate ``D`` (au). This function first computes the Y-component and norm
of the planetocentric velocity vector ``\mathbf{U}``
```math
U_y = \sqrt{a(1-e^2)}\cos i - 1 \quad \text{and} \quad 
U = ||\mathbf{U}|| = \sqrt{3 - \frac{1}{a} - 2\sqrt{a(1-e^2)}\cos i},
```
and then substitutes into `valsecchi_circle(U_y, U_norm, k, h; m_pl=3.003489614915764e-6, a_pl=1.0)`.
`a `, `e` and `i` are the asteroid heliocentric semimajor axis (au), eccentricity and 
inclination (rad) respectively. 

See section 2.1 in page 1181 of https://doi.org/10.1051/0004-6361:20031039.

# Arguments

- `a`: asteroid heliocentric semimajor axis (au).
- `e`: asteroid heliocentric eccentricity.
- `i`: asteroid heliocentric inclination, ecliptic (rad).
- `k/h`: `h` heliocentric revolutions of asteroid per `k` heliocentric revolutions of Earth.
- `m_pl`: planet mass normalized to Sun's mass, equal to Earth mass in solar masses by default.
"""
function valsecchi_circle(a, e, i, k, h; m_pl=3.003489614915764e-6)
    # Components and norm of the planetocentric velocity vector
    # See section 2.1 in page 1181 of https://doi.org/10.1051/0004-6361:20031039
    # U_x = sqrt( 2 - (1/a) - a*(1-(e^2)) ) # TODO: CHECK SIGN
    U_y = sqrt( a*(1-(e^2)) )*cos(i) - 1
    # U_z = sqrt( a*(1-(e^2)) )*sin(i) # TODO: CHECK SIGN
    U_norm = sqrt( 3 - (1/a) - 2*sqrt(a*(1-(e^2)))*cos(i) )
    # Expression below should be equal to asteroid heliocentric elliptic semimamajor axis
    # in au units
    # @show 1/(1-U_^2-2U_y)
    return valsecchi_circle(U_y, U_norm, k, h, m_pl=m_pl)
end

@doc raw"""
    valsecchi_circle(U_y, U_norm, k, h; m_pl=3.003489614915764e-6, a_pl=1.0)

Computes Valsecchi circle associated to a mean motion resonance. Returns radius ``R`` (au) and
``\zeta``-axis coordinate ``D`` (au)
```math
R = \left|\frac{c\sin\theta_0'}{\cos\theta_0' - \cos\theta}\right| \quad \text{and} \quad
D = \frac{c\sin\theta}{\cos\theta_0' - \cos\theta},
```
where ``c = m/U^2`` with ``m`` the mass of the planet and ``U = ||\mathbf{U}||`` the norm of 
the planetocentric velocity vector; and ``\theta``, ``\theta_0'`` are the angles between Y-axis
and ``\mathbf{U}`` pre and post encounter respectively. 

See pages 1181, 1182 and 1187 of https://doi.org/10.1051/0004-6361:20031039.

# Arguments 
- `U_y`: Y-component of unperturbed planetocentric velocity (Y-axis coincides with the direction of motion of the planet).
- `U_norm`: Euclidean norm of unperturbed planetocentric velocity. Both `U_y`, `U_norm` are in units such that the heliocentric velocity of the planet is 1.
- `k/h`: `h` heliocentric revolutions of asteroid per `k` heliocentric revolutions of Earth.
- `m_pl`: planet mass normalized to Sun's mass, equal to Earth mass in solar masses by default.
- `a_pl`: planetary heliocentric semimajor axis in au; default value is 1.
"""
function valsecchi_circle(U_y, U_norm, k, h; m_pl=3.003489614915764e-6, a_pl=1.0)
    # Post-encounter semimajor axis
    # See page 1187 of https://doi.org/10.1051/0004-6361:20031039
    a0p = a_pl*(k/h)^(2/3)

    # θ: angle between Y-axis and planetocentric velocity vector \vec{U}

    # Trigonometric functions of pre encounter θ
    # See page 1181 of https://doi.org/10.1051/0004-6361:20031039
    cosθ = U_y/U_norm
    sinθ = sin(acos(cosθ)) # sqrt(1-cosθ^2)  # TODO: CHECK SIGN
    # Trigonometric functions of post-encounter θ
    # See page 1187 of https://doi.org/10.1051/0004-6361:20031039
    cosθ0p = (1-(U_norm^2)-(1/a0p))/(2U_norm)
    sinθ0p = sin(acos(cosθ0p))
    # c = m/U^2 
    # See first sentence below equation (3) in section 2.3, page 1182 of 
    # https://doi.org/10.1051/0004-6361:20031039
    c = m_pl/(U_norm^2)
    # Radius of the Valsecchi circle R and its ζ-axis component D
    # See page 1187 of https://doi.org/10.1051/0004-6361:20031039
    R0 = abs( c*sinθ0p/(cosθ0p-cosθ) )
    D0 = c*sinθ/(cosθ0p-cosθ)

    return R0, D0
end
