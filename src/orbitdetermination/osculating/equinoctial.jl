"""
    EquinoctialElements{T, U} <: AbstractOsculatingElements{T, U}

A set of six equinoctial elements at a given epoch.

# Fields

- `gm::T`: gravitational parameter of the central body [au³/day²].
- `epoch::T`: reference epoch [MJD TDB].
- `frame::Symbol`: reference plane, either `:equatorial` or `:ecliptic`.
- `elements::SVector{6, U}`: set of six equinoctial elements.
- `covariance::SMatrix{6, 6, U, 36}`: covariance matrix.

# Extended help

The components of `elements` are given by:
- `a::U`: semimajor axis [au].
- `h/k::U`: components of the eccentricity vector.
- `λ::U`: mean longitude [deg].
- `p/q::U`: components of the ascending node vector.
"""

@auto_hash_equals struct EquinoctialElements{T, U} <: AbstractOsculatingElements{T, U}
    gm::T
    epoch::T
    frame::Symbol
    elements::SVector{6, U}
    covariance::SMatrix{6, 6, U, 36}
end

elementstype(::EquinoctialElements) = :Equinoctial
elementsnames(::EquinoctialElements) = ["a", "h", "k", "λ", "p", "q"]
elementsunits(::EquinoctialElements) = ["au", "", "", "deg", "", ""]

semimajoraxis(x::EquinoctialElements) = x.elements[1]
eccentricity(x::EquinoctialElements) = sqrt(x.elements[2]^2 + x.elements[3]^2)
meanlongitude(x::EquinoctialElements) = x.elements[4]
eccentricities(x::EquinoctialElements) = x.elements[2:3]
longascnodes(x::EquinoctialElements) = x.elements[5:6]
meanmotion(x::EquinoctialElements) = rad2deg(sqrt(gm(x) / abs(semimajoraxis(x))^3))

keplerequation(h::Number, k::Number, F::Number) = F + h * cos(F) - k * sin(F)
keplerderivative(h::Number, k::Number, F::Number) = 1 - h * sin(F) - k * cos(F)

function eccentriclongitude(h::Number, k::Number, λ::Number)
    # Initial estimate
    F0 = λ
    # Successive approximations via Newtons' method
    for _ in 0:20
        F1 = F0 - (keplerequation(h, k, F0) - λ) / keplerderivative(h, k, F0)
        F0 = F1
    end
    return rad2deg(F0)
end

"""
    keplerian2equinoctial(x, t; kwargs...)

Convert a set of keplerian elements `x` [deg, au] to a set of
equinoctial elements [deeg, au] at an epoch `t` [MJD TDB].

# Keyword arguments

- `μ::Real`: gravitational parameter of the central body
    [au³/day²] (default: `μ_S`).
"""
function keplerian2equinoctial(x::AbstractVector{U}, t::T;
                               μ::Real = μ_S) where {T <: Real, U <: Number}
    # Eccentricity
    e = x[2]
    # Inclination [rad], argument of pericentre [rad] and
    # longitude of ascending node [rad]
    i, ω, Ω = deg2rad(x[3]), deg2rad(x[4]), deg2rad(x[5])
    # Semimajor axis [au] and mean anomaly [rad]
    if 0 < e < 1
        a, M = x[1], deg2rad(x[6])
    elseif e > 1
        q, tp = x[1], x[6]
        a = q / (1 - e)
        n = sqrt(μ / abs(a)^3)
        M = n * (t - tp)
    end
    # Components of the eccentricity vector
    h = e * sin(ω + Ω)
    k = e * cos(ω + Ω)
    # Mean longitude [rad]
    λ = M + ω + Ω
    if 0 < e < 1
        λ = mod2pi(λ)
    end
    # Components of the ascending node vector
    p = tan(i/2) * sin(Ω)
    q = tan(i/2) * cos(Ω)
    # Vector of elements [deg, au]
    λ = rad2deg(λ)
    elements = SVector{6, U}(a, h, k, λ, p, q)

    return elements
end

"""
    equinoctial2keplerian(x, t; kwargs...)

Convert a set of equinoctial elements `x` [deg, au] to a set
of keplerian elements [deg, au] at an epoch `t` [MJD TDB].

# Keyword arguments

- `μ::Real`: gravitational parameter of the central body
    [au³/day²] (default: `μ_S`).
"""
function equinoctial2keplerian(x::AbstractVector{U}, t::T;
                               μ::Real = μ_S) where {T <: Real, U <: Number}
    # Equinoctial elements [au, rad]
    a, h, k, λ, p, q = x[1], x[2], x[3], deg2rad(x[4]), x[5], x[6]
    # Eccentricity
    e = sqrt(h^2 + k^2)
    # Inclination [rad]
    i = 2 * atan(sqrt(p^2 + q^2))
    # Longitude of the ascending node [rad]
    Ω = mod2pi(atan(p, q))
    # Auxiliary angle [rad]
    ζ = mod2pi(atan(h, k))
    # Argument of periapsis [rad]
    ω = mod2pi(ζ - Ω)
    # Mean anomaly [rad]
    M = λ - (ω + Ω)
    if 0 < e < 1
        M = mod2pi(M)
    end
    # Mean motion [rad]
    n = sqrt(μ / abs(a)^3)
    # Time of perihelion passage [MJD TDB]
    tp = t - M / n
    # Perihelion distance [au]
    q = a * (1 - e)
    # Vector of elements [deg, au]
    i, ω, Ω, M = rad2deg(i), rad2deg(ω), rad2deg(Ω), rad2deg(M)
    if 0 < e < 1
        elements = SVector{6, U}(a, e, i, ω, Ω, M)
    elseif e > 1
        elements = SVector{6, U}(q, e, i, ω, Ω, tp)
    end

    return elements
end

"""
    cartesian2equinoctial(x, t; kwargs...)

Convert a cartesian state vector `x` [au, au/day], referred to
an epoch `t` [MJD TDB], to equinoctial elements [deg, au].

# Keyword arguments

- `μ::Real`: gravitational parameter of the central body
    [au³/day²] (default: `μ_S`).
"""
function cartesian2equinoctial(x::AbstractVector{U}; μ::Real = μ_S) where {U <: Number}
    # Position and velocity vectors [au, au/day]
    r_vec, v_vec = x[1:3], x[4:6]
    r = euclid3D(r_vec)
    # Orbital momentum vector [au^2/day]
    h_vec = cross(r_vec, v_vec)
    h = euclid3D(h_vec)
    # Eccentricity vector
    e_vec = cross(v_vec, h_vec) / μ - r_vec / r
    # Semimajor axis
    a = 1 / (2/r - dot(v_vec, v_vec)/μ)
    # Components of the eccentricity vector
    hx, hy, hz = h_vec / h
    p = hx / (1 + hz)
    q = -hy / (1 + hz)
    # Equinoctial reference frame basis vectors
    C = 1 + p^2 + q^2
    f_vec = [1 - p^2 + q^2, 2 * p * q, -2 * p] ./ C
    g_vec = [2 * p * q, 1 + p^2 - q^2, 2 * q] ./ C
    # Components of the ascending node vector
    h = dot(e_vec, g_vec)
    k = dot(e_vec, f_vec)
    # Position in the equinoctial frame [au, au/day]
    X = dot(r_vec, f_vec)
    Y = dot(r_vec, g_vec)
    # Eccentric longitude [rad]
    b = 1 / (1 + sqrt(1 - h^2 - k^2))
    F = mod2pi(atan(
        h + ((1 − h^2 * b) * Y − h * k * b * X) / (a * sqrt(1 - h^2 - k^2)),
        k + ((1 − k^2 * b) * X − h * k * b * Y) / (a * sqrt(1 - h^2 - k^2)),
    ))
    # Mean longitude [rad]
    λ = keplerequation(h, k, F)
    # Vector of elements
    λ = rad2deg(λ)
    elements = SVector{6, U}(a, h, k, λ, p, q)

    return elements
end

function equinoctial2cartesian(x::AbstractVector{U}; μ::Real = μ_S) where {U <: Number}
    # Equinoctial elements [rad, au]
    a, h, k, λ, p, q = x[1], x[2], x[3], deg2rad(x[4]), x[5], x[6]
    # Equinoctial reference frame basis vectors
    C = 1 + p^2 + q^2
    f_vec = [1 - p^2 + q^2, 2 * p * q, -2 * p] ./ C
    g_vec = [2 * p * q, 1 + p^2 - q^2, 2 * q] ./ C
    # w_vec = [2 * p, -2 * q, 1 - p^2 - q^2] ./ C
    # Eccentric longitude [rad]
    F = deg2rad(eccentriclongitude(h, k, λ))
    # Mean motion [rad]
    n = sqrt(μ / abs(a)^3)
    # True longitude
    b = 1 / (1 + sqrt(1 - h^2 - k^2))
    sinL = ( (1 − k^2 * b) * sin(F) + h * k * b * cos(F) − h) / keplerderivative(h, k, F)
    cosL = ( (1 − h^2 * b) * cos(F) + h * k * b * sin(F) − k) / keplerderivative(h, k, F)
    # Distance to the central body [au]
    r = a * keplerderivative(h, k, F)
    # Position and velocity in the equinoctial frame [au, au/day]
    r_o = r .* [cosL, sinL, zero(F)]
    v_o = (n * a / sqrt(1 - h^2 - k^2)) .* [-(h + sinL), k + cosL, zero(F)]
    # Transform r_o and v_o from the equinoctial to the inertial frame [au, au/day]
    r_i = r_o[1] * f_vec + r_o[2] * g_vec
    v_i = v_o[1] * f_vec + v_o[2] * g_vec
    # Cartesian state vector [au, au/day]
    rv_i = SVector{6, U}(r_i[1], r_i[2], r_i[3], v_i[1], v_i[2], v_i[3])

    return rv_i
end

# Convert a set of equinoctial elements `x` [deg, au] to a cartesian state
# vector [au, au/day] at epoch `t` [MJD TDB] (default: `epoch(x)`).
function (x::EquinoctialElements)(t::Number = epoch(x))
    # Equinoctial elements [deg, au]
    eqn = collect(elements(x))
    # Mean longitude [deg]
    eqn[4] = meanlongitude(x) + meanmotion(x) * (t - epoch(x))
    # Cartesian state vector [au, au/day]
    rv = equinoctial2cartesian(eqn; μ = gm(x))

    return rv
end