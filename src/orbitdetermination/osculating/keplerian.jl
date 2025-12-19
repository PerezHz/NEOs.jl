"""
    KeplerianElements{T, U} <: AbstractOsculatingElements{T, U}

A set of six keplerian orbital elements at a given epoch.

# Fields

- `gm::T`: gravitational parameter of the central body [au³/day²].
- `epoch::T`: reference epoch [MJD TDB].
- `frame::Symbol`: reference plane, either `:equatorial` or `:ecliptic`.
- `elements::SVector{6, U}`: set of six keplerian orbital elements.
- `covariance::SMatrix{6, 6, U, 36}`: covariance matrix.

# Extended help

The components of `elements` depend on the type of orbit. For an ellipse:
- `a::U`: semimajor axis [au].
- `e::U`: eccentricity.
- `i::U`: inclination [deg].
- `ω::U`: argument of pericentre [deg].
- `Ω::U`: longitude of ascending node [deg].
- `M::U`: mean anomaly [deg].
For a hyperbola:
- `q::U`: pericenter distance [au].
- `e::U`: eccentricity.
- `i::U`: inclination [deg].
- `ω::U`: argument of pericentre [deg].
- `Ω::U`: longitude of ascending node [deg].
- `tp::U`: time of pericenter passage [MJD TDB].
"""
@auto_hash_equals struct KeplerianElements{T, U} <: AbstractOsculatingElements{T, U}
    gm::T
    epoch::T
    frame::Symbol
    elements::SVector{6, U}
    covariance::SMatrix{6, 6, U, 36}
end

elementstype(::KeplerianElements) = :Keplerian

function elementsnames(x::KeplerianElements)
    flag = iselliptic(x)
    return [flag ? "a" : "q", "e", "i", "ω", "Ω", flag ? "M" : "tp"]
end

function elementsunits(x::KeplerianElements)
    flag = iselliptic(x)
    return ["au", "", "deg", "deg", "deg", flag ? "deg" : "MJD"]
end

function semimajoraxis(x::KeplerianElements)
    iselliptic(x) && return x.elements[1]
    # Pericenter distance [au] and eccentricity
    q, e = x.elements[1], x.elements[2]
    # Semimajor axis
    a = q / (1 - e)
    return a
end

eccentricity(x::KeplerianElements) = x.elements[2]
inclination(x::KeplerianElements) = x.elements[3]
argperi(x::KeplerianElements) = x.elements[4]
longascnode(x::KeplerianElements) = x.elements[5]
meanmotion(x::KeplerianElements) = rad2deg(sqrt(gm(x) / abs(semimajoraxis(x))^3))
meananomaly(x::KeplerianElements, t::Real) = meanmotion(x) * (t - timeperipass(x))

function timeperipass(x::KeplerianElements)
    ishyperbolic(x) && return x.elements[end]
    # Semimajor axis [au] and mean anomaly [rad]
    a, M = x.elements[1], deg2rad(x.elements[end])
    # Time of pericenter passage [MJD TDB]
    tp = x.epoch - M / sqrt(gm(x) / abs(a)^3)
    return tp
end

function eccentricanomaly(x::KeplerianElements, t::Real)
    # Boolean flag
    flag = iselliptic(x)
    # Eccentricity
    e = eccentricity(x)
    # Mean anomaly [rad]
    M = deg2rad(meananomaly(x, t))
    # Curtis (2020) initial estimate
    _M_ = mod2pi(M)
    if flag
        E0 = _M_ < π ? _M_ + e/2 : _M_ - e/2
    else
        E0 = _M_
    end
    # Successive approximations via Newtons' method
    for _ in 0:10
        # TODO: implement modified Newton's method for Kepler's equation (Murray-Dermott)
        E1 = E0 - (keplerequation(e, E0, Val(flag)) - M) / keplerderivative(e, E0, Val(flag))
        E0 = E1
    end
    return rad2deg(E0)
end

keplerequation(e::Number, E::Number, ::Val{true}) = E - e * sin(E)
keplerequation(e::Number, E::Number, ::Val{false}) = e * sinh(E) - E

keplerderivative(e::Number, E::Number, ::Val{true}) = 1 - e * cos(E)
keplerderivative(e::Number, E::Number, ::Val{false}) = e * cosh(E) - 1

function trueanomaly(x::KeplerianElements, t::Real)
    # Eccentricity
    e = eccentricity(x)
    # Eccentric anomaly [rad]
    E = deg2rad(eccentricanomaly(x, t))
    # True anomaly [rad]
    if iselliptic(x)
        ν = 2 * atan( sqrt( (1+e) / (1-e) ) * tan(E / 2) )
    elseif ishyperbolic(x)
        ν = 2 * atan( sqrt( (e+1) / (e-1) ) * tanh(E / 2) )
    end
    return rad2deg(ν)
end

"""
    cartesian2keplerian(x, t; kwargs...)

Convert a cartesian state vector `x` [au, au/day], referred to an epoch `t`
[MJD TDB], to keplerian orbital elements.

# Keyword arguments

- `μ::Real`: gravitational parameter of the central body [au³/day²] (default: `μ_S`).
- `frame::Symbol`: reference plane, either `:equatorial` (default) or `:ecliptic`.
"""
function cartesian2keplerian(x::AbstractVector{U}, t::T; μ::Real = μ_S,
                             frame::Symbol = :equatorial) where {T <: Real, U <: Number}
    if frame == :ecliptic
        x = equatorial2ecliptic(x)
    end
    # Position and velocity vectors [au, au/day]
    r_vec, v_vec = x[1:3], x[4:6]
    r = euclid3D(r_vec)
    # Orbital momentum vector [au^2/day]
    h_vec = cross(r_vec, v_vec)
    h = euclid3D(h_vec)
    # Eccentricity vector
    e_vec = cross(v_vec, h_vec) / μ - r_vec / r
    e = euclid3D(e_vec)
    # Semimajor axis
    a = 1 / (2/r - dot(v_vec, v_vec)/μ)
    # Perihelion distance
    q = a * (1 - e)
    # Orbit inclination [rad]
    i = acos(h_vec[3] / h)
    # Vector n pointing towards the ascending node [au^2/day]
    n_vec = [-h_vec[2], h_vec[1], zero(T)]
    n = euclid3D(n_vec)
    # Argument of periapsis [rad]
    if e_vec[3] ≥ 0
        ω = acos(dot(n_vec, e_vec) / (n * e))
    else
        ω = 2π - acos(dot(n_vec, e_vec) / (n * e))
    end
    # Longitude of the ascending node [rad]
    if n_vec[2] ≥ 0
        Ω = acos(n_vec[1] / n)
    else
        Ω = 2π - acos(n_vec[1] / n)
    end
    # True anomaly [rad]
    if dot(r_vec, v_vec) ≥ 0
        ν = acos(dot(e_vec, r_vec) / (e * r))
    else
        ν = 2π - acos(dot(e_vec, r_vec) / (e * r))
    end
    # Eccentric anomaly [rad]
    if 0 < e < 1
        E = 2 * atan(tan(ν/2) / sqrt((1 + e) / (1 - e)))
    elseif e > 1
        E = 2 * atanh(tan(ν/2) / sqrt((e + 1) / (e - 1)))
    end
    # Mean anomaly [rad]
    if 0 < e < 1
        M = mod2pi(E - e * sin(E))
    elseif e > 1
        M = e * sinh(E) - E
    end
    # Time of perihelion passage [MJD TDB]
    tp = t - M / sqrt(μ / abs(a)^3)
    # Vector of elements
    i, ω, Ω, M = rad2deg(i), rad2deg(ω), rad2deg(Ω), rad2deg(M)
    if 0 < e < 1
        elements = SVector{6, U}(a, e, i, ω, Ω, M)
    elseif e > 1
        elements = SVector{6, U}(q, e, i, ω, Ω, tp)
    end

    return elements
end

# Convert a set of keplerian orbital elements `x` to a cartesian state
# vector [au, au/day] at epoch `t` [MJD TDB] (default: `epoch(x)`).
function (x::KeplerianElements)(t::Number = epoch(x))
    # Eccentricity, semimajor axis [au]
    e, a = eccentricity(x), semimajoraxis(x)
    # Inclination [rad], argument of pericentre [rad] and
    # longitude of ascending node [rad]
    i = deg2rad(inclination(x))
    ω = deg2rad(argperi(x))
    Ω = deg2rad(longascnode(x))
    # Eccentric anomaly [rad]
    E = deg2rad(eccentricanomaly(x, t))
    # True anomaly [rad]
    f = deg2rad(trueanomaly(x, t))
    # Distance to the central body [au]
    r = a * (1 - e^2) / (1 + e * cos(f))
    # Obtain position and velocity in the orbital frame [au, au/day]
    r_o = r .* [cos(f), sin(f), zero(f)]
    if iselliptic(x)
        v_o = (sqrt(gm(x) * a) / r) .* [-sin(E), sqrt(1 - e^2) * cos(E), zero(E)]
    elseif ishyperbolic(x)
        v_o = (sqrt(-gm(x) * a) / r) .* [-sinh(E), sqrt(e^2 - 1) * cosh(E), zero(E)]
    end
    # Transform r_o and v_o from the orbital to the inertial frame [au, au/day]
    A = Rz(-Ω) * Rx(-i) * Rz(-ω)
    r_i = A * r_o
    v_i = A * v_o
    # State vector [au, au/day]
    pv_i = vcat(r_i, v_i)

    return pv_i
end