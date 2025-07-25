"""
    AbstractOsculatingElements{T <: Real, U <: Number}

Supertype for the osculating orbital elements interface.
"""
abstract type AbstractOsculatingElements{T <: Real, U <: Number} end

"""
    OsculatingElements{T, U} <: AbstractOsculatingElements{T, U}

A set of six osculating orbital elements at a given epoch.

# Fields

- `mu::T`: gravitational parameter of the central body [au³/day²].
- `epoch::T`: reference epoch [MJD TDB].
- `frame::Symbol`: reference plane, either `:equatorial` or `:ecliptic`.
- `elements::SVector{6, U}`: set of six osculating orbital elements.
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
@auto_hash_equals struct OsculatingElements{T, U} <: AbstractOsculatingElements{T, U}
    mu::T
    epoch::T
    frame::Symbol
    elements::SVector{6, U}
    covariance::SMatrix{6, 6, U, 36}
end

numtypes(::OsculatingElements{T, U}) where {T, U} = T, U

epoch(x::OsculatingElements) = x.epoch
elements(x::OsculatingElements) = x.elements
covariance(x::OsculatingElements) = x.covariance
sigmas(x::OsculatingElements) = sqrt.(diag(covariance(x)))

iscircular(x::OsculatingElements) = iszero(x.elements[2])
iselliptic(x::OsculatingElements) = 0 < x.elements[2] < 1
isparabolic(x::OsculatingElements) = isone(x.elements[2])
ishyperbolic(x::OsculatingElements) = x.elements[2] > 1

function evaldeltas(osc::OsculatingElements{T, TaylorN{T}},
                    y::Vector{T} = zeros(T, get_numvars())) where {T <: Real}
    return OsculatingElements{T, T}(osc.mu, osc.epoch, osc.frame, osc.elements(y),
                                    osc.covariance(y))
end

# Print method for OsculatingElements
function show(io::IO, x::OsculatingElements)
    flag = iselliptic(x)
    T, U = numtypes(x)
    header = flag ? "Elliptic" : "Hyperbolic"
    mu = x.mu
    t0 = x.epoch
    d0 = julian2datetime(t0 + JD_J2000 - MJD2000)
    frame = x.frame
    e0 = cte.(x.elements)
    σ0 = sqrt.(diag(x.covariance))
    se0 = [rpad(@sprintf("%+.12E", e0[i]), 25) for i in eachindex(e0)]
    sσ0 = [rpad(@sprintf("%+.12E", σ0[i]), 25) for i in eachindex(σ0)]
    names = [flag ? "a" : "q", "e", "i", "ω", "Ω", flag ? "M" : "tp"]
    units = ["au", "", "deg", "deg", "deg", flag ? "deg" : "MJD"]

    print(io,
        "$header{$T, $U} osculating elements\n",
        repeat("-", 67), "\n",
        "Mu: $mu au³/day²\n",
        "Epoch: $t0 MJD ($d0 TDB)\n",
        "Frame: $frame\n",
        repeat("-", 67), "\n",
        "Variable    Nominal value            Uncertainty              Units\n",
        rpad(names[1], 12), se0[1], sσ0[1], units[1], "\n",
        rpad(names[2], 12), se0[2], sσ0[2], units[2], "\n",
        rpad(names[3], 12), se0[3], sσ0[3], units[3], "\n",
        rpad(names[4], 12), se0[4], sσ0[4], units[4], "\n",
        rpad(names[5], 12), se0[5], sσ0[5], units[5], "\n",
        rpad(names[6], 12), se0[6], sσ0[6], units[6], "\n",
    )
end

# Rotate state vector `xas` from equatorial plane to the ecliptic
function equatorial2ecliptic(xas::AbstractVector)
    # Rotation matrix (only positions)
    m_eq2ecl = Rx(deg2rad(ϵ0_deg))
    # Rotational matrix (positions + velocities)
    m_xv_eq2ecl = hcat(vcat(m_eq2ecl, zeros(3,3)), vcat(zeros(3,3), m_eq2ecl))
    # Rotated state vector
    return m_xv_eq2ecl * xas
end

evaluate(y::SMatrix{S1, S2, TaylorN{T}, L}, x::Vector{<:Number}) where {S1, S2,
    T <: Number, L} = [y[i, j](x) for i in axes(y, 1), j in axes(y, 2)]

(y::SMatrix{S1, S2, TaylorN{T}, L})(x::Vector{T}) where {S1, S2, T <: Number, L} =
    evaluate(y, x)

function semimajoraxis(x::OsculatingElements)
    iselliptic(x) && return x.elements[1]
    # Pericenter distance [au] and eccentricity
    q, e = x.elements[1], x.elements[2]
    # Semimajor axis
    a = q / (1 - e)
    return a
end

eccentricity(x::OsculatingElements) = x.elements[2]
inclination(x::OsculatingElements) = x.elements[3]
argperi(x::OsculatingElements) = x.elements[4]
longascnode(x::OsculatingElements) = x.elements[5]
meanmotion(x::OsculatingElements) = rad2deg(sqrt(x.mu / abs(semimajoraxis(x))^3))
meananomaly(x::OsculatingElements, t::Real) = meanmotion(x) * (t - timeperipass(x))

function timeperipass(x::OsculatingElements)
    ishyperbolic(x) && return x.elements[end]
    # Semimajor axis [au] and mean anomaly [rad]
    a, M = x.elements[1], deg2rad(x.elements[end])
    # Time of pericenter passage [MJD TDB]
    tp = x.epoch - M / sqrt(x.mu / abs(a)^3)
    return tp
end

function eccentricanomaly(x::OsculatingElements, t::Real)
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

function trueanomaly(x::OsculatingElements, t::Real)
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
    cartesian2osculating(xas, epoch; kwargs...)

Convert a cartesian state vector `xas` [au, au/day], referred to an `epoch`
[MJD TDB], to osculating orbital elements.

# Keyword arguments

- `μ::Real`: gravitational parameter of the central body [au³/day²] (default: `μ_S`).
- `frame::Symbol`: reference plane, either `:equatorial` (default) or `:ecliptic`.
- `Γ_car::AbstractMatrix`: covariance matrix of `xas`.
"""
function cartesian2osculating(xas::AbstractVector{U}, epoch::T; μ::Real = μ_S,
                              frame::Symbol = :equatorial,
                              Γ_car::AbstractMatrix{T}) where {T <: Real, U <: Number}
    if frame == :ecliptic
        xas = equatorial2ecliptic(xas)
    end
    # Position and velocity vectors [au, au/day]
    r_vec, v_vec = xas[1:3], xas[4:6]
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
        M = E - e * sin(E)
    elseif e > 1
        M = e * sinh(E) - E
    end
    # Time of perihelion passage [MJD TDB]
    tp = epoch - M / sqrt(μ / abs(a)^3)
    # Vector of elements
    i, ω, Ω, M = rad2deg(i), rad2deg(ω), rad2deg(Ω), rad2deg(M)
    if 0 < e < 1
        elements = SVector{6, U}(a, e, i, ω, Ω, M)
    elseif e > 1
        elements = SVector{6, U}(q, e, i, ω, Ω, tp)
    end
    # Covariance matrix
    Γ_osc = project(elements, Γ_car)

    return OsculatingElements{T, U}(μ, epoch, frame, elements, Γ_osc)
end

# Return the cartesian state vector of x at time t [MJD TDB]
function (x::OsculatingElements)(t::Number = epoch(x))
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
        v_o = (sqrt(x.mu * a) / r) .* [-sin(E), sqrt(1 - e^2) * cos(E), zero(E)]
    elseif ishyperbolic(x)
        v_o = (sqrt(-x.mu * a) / r) .* [-sinh(E), sqrt(e^2 - 1) * cosh(E), zero(E)]
    end
    # Transform r_o and v_o from the orbital to the inertial frame [au, au/day]
    A = Rz(-Ω) * Rx(-i) * Rz(-ω)
    r_i = A * r_o
    v_i = A * v_o
    # State vector [au, au/day]
    pv_i = vcat(r_i, v_i)

    return  pv_i
end

@doc raw"""
    yarkp2adot(A2, a, e; kwargs...)

Return the average semimajor axis drift of an orbit with Yarkovsky coefficient `A2`,
semimajor axis `a` and eccentricity `e`.

# Keyword argument

- `μ`: gravitational parameter of the central body (default: `μ_S`).

!!! reference
    See:
    - https://doi.org/10.1016/j.icarus.2013.02.004

# Extended help

The average semimajor axis drift is given by:
```math
\begin{align*}
    \left\langle\dot{a}\right\rangle
    & = \frac{2A_2(1-e^2)}{n}\left(\frac{1 \ \text{AU}}{p}\right)^2 \\
    & = \frac{2A_2}{(1-e^2)\sqrt{a\mu_\odot}}(1 \ \text{AU})^2,
\end{align*}
```
where ``A_2`` is the Yarkovsky parameter, ``\mu_\odot = GM_\odot`` is the Sun's
gravitational parameter, ``e`` is the eccentricity, ``n = \sqrt{\mu/a^3}`` is the
mean motion, ``p = a(1-e^2)`` is the semilatus rectum, and ``a`` is the semimajor axis.
"""
yarkp2adot(A2, a, e; μ = μ_S) = 2A2 / (sqrt(a) * (1 - e^2) * sqrt(μ))
