"""
    KeplerianElements{T, U} <: AbstractOsculatingElements{T, U}

A set of six keplerian elements at a given epoch.

# Fields

- `gm::T`: gravitational parameter of the central body [au³/day²].
- `epoch::T`: reference epoch [MJD TDB].
- `frame::Symbol`: reference plane, either `:equatorial` or `:ecliptic`.
- `elements::SVector{6, U}`: set of six keplerian elements.
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
    q, e = pericenter(x), eccentricity(x)
    # Semimajor axis [au]
    a = q / (1 - e)
    return a
end

function pericenter(x::KeplerianElements)
    ishyperbolic(x) && return x.elements[1]
    # Semimajor axis [au] and eccentricity
    a, e = semimajoraxis(x), eccentricity(x)
    # Pericenter distance [au]
    q = a * (1 - e)
    return q
end

eccentricity(x::KeplerianElements) = x.elements[2]
inclination(x::KeplerianElements) = x.elements[3]
argperi(x::KeplerianElements) = x.elements[4]
longascnode(x::KeplerianElements) = x.elements[5]
meanmotion(x::KeplerianElements) = rad2deg(sqrt(gm(x) / abs(semimajoraxis(x))^3))

function meananomaly(x::KeplerianElements)
    iselliptic(x) && return x.elements[end]
    # Mean motion [deg] and time of pericenter passage [MJD TDB]
    n, tp = meanmotion(x), timeperipass(x)
    # Mean anomaly [deg]
    M = n * (epoch(x) - tp)
    return M
end

function timeperipass(x::KeplerianElements)
    ishyperbolic(x) && return x.elements[end]
    # Mean motion [deg] and mean anomaly [deg]
    n, M = meanmotion(x), meananomaly(x)
    # Time of pericenter passage [MJD TDB]
    tp = epoch(x) - M / n
    return tp
end

keplerequation(e::Number, E::Number, ::Val{true}) = E - e * sin(E)
keplerequation(e::Number, E::Number, ::Val{false}) = e * sinh(E) - E

keplerderivative(e::Number, E::Number, ::Val{true}) = 1 - e * cos(E)
keplerderivative(e::Number, E::Number, ::Val{false}) = e * cosh(E) - 1

for V in (:(true), :(false))
    @eval begin

        function eccentricanomaly(e::Number, M::Number, ::Val{$V})
            # Curtis (2020) initial estimate
            if $V
                M = mod2pi(M)
                E0 = M < π ? M + e/2 : M - e/2
            else
                E0 = M
            end
            # Successive approximations via Newtons' method
            for _ in 0:20
                # TODO: implement modified Newton's method for Kepler's equation (Murray-Dermott)
                E1 = E0 - (keplerequation(e, E0, Val($V)) - M) / keplerderivative(e, E0, Val($V))
                E0 = E1
            end
            return rad2deg(E0)
        end

        function trueanomaly(e::Number, E::Number, ::Val{$V})
            # True anomaly [rad]
            if $V
                ν = 2 * atan( sqrt( (1+e) / (1-e) ) * tan(E / 2) )
            else
                ν = 2 * atan( sqrt( (e+1) / (e-1) ) * tanh(E / 2) )
            end
            return rad2deg(ν)
        end

    end
end

"""
    cartesian2keplerian(x, t; kwargs...)

Convert a cartesian state vector `x` [au, au/day], referred to
an epoch `t` [MJD TDB], to keplerian elements [deg, au].

# Keyword arguments

- `μ::Real`: gravitational parameter of the central body
    [au³/day²] (default: `μ_S`).
"""
function cartesian2keplerian(x::AbstractVector{U}, t::T;
                             μ::Real = μ_S) where {T <: Real, U <: Number}
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
    # Time of pericenter passage [MJD TDB]
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

"""
    keplerian2cartesian(x, t; kwargs...)

Convert a set of keplerian elements `x` [deg, au], referred to an
epoch `t` [MJD TDB], to a cartesian state vector [au, au/day].

# Keyword arguments

- `μ::Real`: gravitational parameter of the central body
    [au³/day²] (default: `μ_S`).
"""
function keplerian2cartesian(x::AbstractVector{U}, t::T;
                             μ::Real = μ_S) where {T <: Real, U <: Number}
    # Eccentricity
    e = x[2]
    # Boolean flag
    flag = 0 < e < 1
    # Inclination [rad], argument of pericentre [rad] and
    # longitude of ascending node [rad]
    i, ω, Ω = deg2rad(x[3]), deg2rad(x[4]), deg2rad(x[5])
    # Semimajor axis [au] and mean anomaly [rad]
    if flag
        a, M = x[1], deg2rad(x[6])
    else
        q, tp = x[1], x[6]
        a = q / (1 - e)
        n = sqrt(μ / abs(a)^3)
        M = n * (t - tp)
    end
    # Eccentric anomaly [rad]
    E = deg2rad(eccentricanomaly(e, M, Val(flag)))
    # True anomaly [rad]
    f = deg2rad(trueanomaly(e, E, Val(flag)))
    # Distance to the central body [au]
    r = a * (1 - e^2) / (1 + e * cos(f))
    # Obtain position and velocity in the orbital frame [au, au/day]
    r_o = r .* [cos(f), sin(f), zero(f)]
    if flag
        v_o = (sqrt(μ * a) / r) .* [-sin(E), sqrt(1 - e^2) * cos(E), zero(E)]
    else
        v_o = (sqrt(-μ * a) / r) .* [-sinh(E), sqrt(e^2 - 1) * cosh(E), zero(E)]
    end
    # Transform r_o and v_o from the orbital to the inertial frame [au, au/day]
    A = Rz(-Ω) * Rx(-i) * Rz(-ω)
    r_i = A * r_o
    v_i = A * v_o
    # Cartesian state vector [au, au/day]
    rv_i = SVector{6, U}(r_i[1], r_i[2], r_i[3], v_i[1], v_i[2], v_i[3])

    return rv_i
end

# Convert a set of keplerian elements `x` [deg, au] to a cartesian state
# vector [au, au/day] at epoch `t` [MJD TDB] (default: `epoch(x)`).
function (x::KeplerianElements)(t::Number = epoch(x))
    # Keplerian elements [deg, au]
    kep = collect(elements(x))
    # Mean anomaly [deg]
    if iselliptic(x)
        kep[end] = meananomaly(x) + meanmotion(x) * (t - epoch(x))
    end
    # Cartesian state vector [au, au/day]
    rv = keplerian2cartesian(kep, t; μ = gm(x))

    return rv
end