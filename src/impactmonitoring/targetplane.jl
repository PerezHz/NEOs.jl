"""
    BPlane{U} <: AbstractTargetPlane{U}

B-Plane in Öpik's coordinates for a hyperbolic planetary close encounter.

# Fields

- `ξ/ζ::U`: Öpik's coordinates [planet radii].
- `b::U`: planet's impact cross section ("critical B") [planet radii].
- `U_y/U_norm::U`: y-component and norm of the planetocentric velocity vector [km/s].
"""
struct BPlane{U} <: AbstractTargetPlane{U}
    ξ::U
    ζ::U
    b::U
    U_y::U
    U_norm::U
end

targetplane(x::BPlane) = [x.ξ, x.ζ, x.b]

@doc raw"""
    bopik(xae, xes)

Return the [`BPlane`](@ref) of a hyperbolic planetary close encounter. `xae`
is the asteroid's planetocentric cartesian state vector at close approach
[au, day], and `xes` is the planet's heliocentric cartesian state vector at
close approach [au, day].

!!! reference
    See equations (37)-(38) in page 14 of:
    - https://doi.org/10.1007/s10569-019-9914-4

# Extended help

Öpik's coordinates are given by:
```math
ξ = \mathbf{B}\cdot\hat{\mathbf{\xi}}/R_E \quad \text{and} \quad
ζ = \mathbf{B}\cdot\hat{\mathbf{\zeta}}/R_E,
```
where ``\mathbf{B}`` is the impact parameter vector, ``(\hat{\mathbf{\xi}},
\hat{\mathbf{\zeta}})`` is Öpik's frame and ``R_E`` is the planet's radius [au].

For the computation of the planet's impact cross section ("critical B"), see
[`crosssection`](@ref).
"""
function bopik(xae::AbstractVector{U}, xes::AbstractVector{U}) where {U <: Number}

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
    # See first sentence below equation (3) in page 2 of
    # https://doi.org/10.1007/s10569-019-9914-4
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
    ves = v_pl # Earth's velocity, au/day
    ves_norm = sqrt(ves[1]^2 + ves[2]^2 + ves[3]^2) # Earth's speed, au/day
    ves_unit = ves/ves_norm # Earth's velocity unit vector
    # angle between Y-axis and \vec{U}
    cosθ = dot(S_v, ves_unit)
    v_infty_kms = v_infty*au/daysec # Asteroid unperturbed speed, km/sec
    # The norm of \vec{U} in appropriate units
    U_unit = ves_norm # 1 U = v_ast/v_pl
    U_norm = v_infty/U_unit
    # U_y
    U_y = U_norm*cosθ

    # Earth impact cross section ("critical B")
    b_E = crosssection(μ_E, RE/au, v_infty)

    return BPlane{U}(B_dot_ξ*au/RE, B_dot_ζ*au/RE, b_E, U_y, U_norm)
end

"""
    MTP{U} <: AbstractTargetPlane{U}

Modified Target Plane for a planetary close encounter.

# Fields

- `X/Y::U`: target plane cartesian coordinates [planet radii].
"""
struct MTP{U} <: AbstractTargetPlane{U}
    X::U
    Y::U
end

targetplane(x::MTP) = [x.X, x.Y, one(x.X)]

"""
    mtp(xae)

Return the [`MTP`](@ref) of a planetary close encounter. `xae` is the asteroid's
planetocentric cartesian state vector at close approach [au, day].

!!! reference
    See equations (43)-(44) in page 15 of:
    - https://doi.org/10.1007/s10569-019-9914-4
"""
function mtp(xae::Vector{U}) where {U <: Number}
    # Unfold geocentric position and velocity
    r, v = xae[1:3], xae[4:6]
    # Reference direction
    V = [0, 0, -1]
    # Unit vectors
    ez = v ./ euclid3D(v)
    _ey_ = cross(V, ez)
    ey = _ey_ ./ euclid3D(_ey_)
    ex = cross(ey, ez)
    # Cartesian coordinates [Earth radii]
    X = dot3D(r, ex)*au/RE
    Y = dot3D(r, ey)*au/RE

    return MTP{U}(X, Y)
end

@doc raw"""
    crosssection(μ_P, R_P, vinf)

Return the effective cross section [planet radii] of a planet with gravitational
parameter `μ_P` [au³/day²] and physical radius `R_P` [au], considering an object
approaching with asymptotic inbound velocity `vinf` [au/day].

!!! reference
    See equations (13)-(14) in pages 4-5 of:
    - https://doi.org/10.1007/s10569-019-9914-4

# Extended help

The effective cross section ``B`` is derived from the conservation of energy and
angular momentum; it represents the impact parameter corresponding to a grazing
impact in a hyperbolic close encounter. ``B`` is given by:
```math
B = \sqrt{1 + \frac{2\mu_P}{R_P v_\infty^2}},
```
where ``\mu_P`` is the planet's gravitational parameter, ``R_P`` is the planet's
physical radius and ``v_\infty`` is the asymptotic inbound velocity. If actual ``B``
is equal or less to this, then impact happens.
"""
crosssection(μ_P, R_P, vinf) = sqrt( 1 + (2μ_P)/(R_P*(vinf)^2) )

@doc raw"""
    valsecchi_circle(a, e, i, k, h; kwargs...)

Compute the Valsecchi circle associated to the `k:h` mean motion resonance between
a planet and an asteroid with heliocentric semimajor axis `a` [au], eccentricity `e`
and (ecliptic) inclination `i` [rad]. Returns the radius [au] and the ``\zeta``-axis
coordinate [au].

# Keyword argument

- `m_pl`: planet mass normalized to Sun's mass (default: Earth mass in solar masses).

!!! reference
    See section 2.1 in page 1181 of:
    - https://doi.org/10.1051/0004-6361:20031039

# Extended help

A `k:h` resonance corresponds to `h` heliocentric revolutions of the asteroid per `k`
heliocentric revolutions of the planet.

This function first computes the Y-component and norm of the planetocentric velocity
vector ``\mathbf{U}``, which are given by:
```math
\begin{align*}
    U_y & = \sqrt{a(1-e^2)}\cos i - 1\\
    U   & = ||\mathbf{U}|| = \sqrt{3 - \frac{1}{a} - 2\sqrt{a(1-e^2)}\cos i},
\end{align*}
```
where `a `, `e` and `i` are the asteroid's heliocentric semimajor axis [au], eccentricity
and (ecliptic) inclination [rad], respectively. Then, it substitutes into
`valsecchi_circle(U_y, U_norm, k, h; m_pl)`.
"""
function valsecchi_circle(a, e, i, k, h; m_pl=3.003489614915764e-6)
    # Components and norm of the planetocentric velocity vector
    # See section 2.1 in page 1181 of https://doi.org/10.1051/0004-6361:20031039
    # U_x = sqrt( 2 - (1/a) - a*(1-(e^2)) ) # TODO: CHECK SIGN
    U_y = sqrt( a*(1-(e^2)) )*cos(i) - 1
    # U_z = sqrt( a*(1-(e^2)) )*sin(i) # TODO: CHECK SIGN
    U_norm = sqrt( 3 - (1/a) - 2*sqrt(a*(1-(e^2)))*cos(i) )
    # The following expression should be equal to asteroid heliocentric elliptic
    # semimamajor axis in au units: 1/(1-U_^2-2U_y)
    return valsecchi_circle(U_y, U_norm, k, h, m_pl=m_pl)
end

@doc raw"""
    valsecchi_circle(U_y, U_norm, k, h; kwargs...)

Compute the Valsecchi circle associated to the `k:h` mean motion resonance between
a planet and an asteroid, where `U_y` and `U_norm` are the Y-component and Euclidean
norm of the unperturbed planetocentric velocity vector. The Y-axis coincides with the
direction of motion of the planet. Both `U_y`, `U_norm` are in units such that the
heliocentric velocity of the planet is 1. Returns the radius [au] and the ``\zeta``-axis
coordinate [au].

# Keyword arguments

- `m_pl`: planet mass normalized to Sun's mass (default: Earth mass in solar masses).
- `a_pl`: planetary heliocentric semimajor axis [a] (default: `1.0`).

!!! reference
    See pages 1181, 1182 and 1187 of:
    - https://doi.org/10.1051/0004-6361:20031039

# Extended help

A `k:h` resonance corresponds to `h` heliocentric revolutions of the asteroid per `k`
heliocentric revolutions of the planet.

The radius ``R`` [au] and ``\zeta``-axis coordinate ``D`` [au] are given by:
```math
\begin{\align*}
    R & = \left|\frac{c\sin\theta_0'}{\cos\theta_0' - \cos\theta}\right| \\
    D & = \frac{c\sin\theta}{\cos\theta_0' - \cos\theta},
\end{align*}
```
where ``c = m/U^2`` with ``m`` the mass of the planet and ``U = ||\mathbf{U}||`` the norm
of the planetocentric velocity vector; and ``\theta``, ``\theta_0'`` are the angles between
Y-axis and ``\mathbf{U}`` pre and post encounter respectively.
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