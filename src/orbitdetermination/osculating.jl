"""
    OsculatingElements{T <: Number}

Keplerian osculating orbital elements.

# Fields

- `e::T`: eccentricity.
- `q::T`: perihelion distance [au].
- `tp::T`: time of pericenter passage [julian date].
- `M::T`: mean anomaly [deg].
- `Ω::T`: longitude of ascending node [deg].
- `ω::T`: argument of pericentre [deg].
- `i::T`: inclination [deg].
- `a::T`: semimajor axis [au].
"""
@auto_hash_equals struct OsculatingElements{T <: Number}
    e::T
    q::T
    tp::T
    Ω::T
    ω::T
    i::T
    M::T
    a::T
end

OsculatingElements() = OsculatingElements(NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)

# A OsculatingElements is NaN if all its fields are NaN
isnan(x::OsculatingElements) = isnan(x.e) && isnan(x.q) && isnan(x.tp) && isnan(x.Ω) &&
    isnan(x.ω) && isnan(x.i) && isnan(x.M) && isnan(x.a)

# Print method for OsculatingElements
function show(io::IO, x::OsculatingElements)
    print(io,
        rpad("Eccentricity (e):", 36), @sprintf("%.5f", cte(x.e)), "\n",
        rpad("Pericenter distance (q):", 36), @sprintf("%.5f", cte(x.q)), " au\n",
        rpad("Time of pericenter passage (tp):", 36), isnan(x.tp) ? "NaN" :
        julian2datetime(cte(x.tp)), " JDTDB\n",
        rpad("Longitude of Ascending Node (Ω):", 36), @sprintf("%.5f", cte(x.Ω)), " deg\n",
        rpad("Argument of pericenter (ω):", 36), @sprintf("%.5f", cte(x.ω)), " deg\n",
        rpad("Inclination (i):", 36), @sprintf("%.5f", cte(x.i)), " deg\n",
        rpad("Mean anomaly (M):", 36), @sprintf("%.5f", cte(x.M)), " deg\n",
        rpad("Semimajor axis (a): ", 36), @sprintf("%.5f", cte(x.a)), " au\n"
    )
end

# Rotate state vector `xas` from equatorial plane to the ecliptic.
function equatorial2ecliptic(xas::AbstractVector{T}) where {T <: Number}
    # Rotation matrix (only positions)
    m_eq2ecl = Rx(deg2rad(ϵ0_deg))
    # Rotational matrix (positions + velocities)
    m_xv_eq2ecl = hcat(vcat(m_eq2ecl, zeros(3,3)), vcat(zeros(3,3), m_eq2ecl))
    # Rotated state vector
    return m_xv_eq2ecl * xas
end

"""
    pv2kep(xas; kwargs...)

Return the Keplerian osculating orbital elements of a cartesian state vector `xas`.

# Keyword arguments

- `μ::Real`: gravitational parameter of the central body (default: `μ_S`).
- `jd::Real`: epoch of reference [julian days] (default: `JD_J2000`).
- `frame::Symbol`: plane of reference, one of `:equatorial` (default) or `:ecliptic`.
"""
function pv2kep(xas::AbstractVector{T}; μ::Real = μ_S, jd::Real = JD_J2000,
                frame::Symbol = :equatorial) where {T <: Number}
    if frame == :ecliptic
        xas = equatorial2ecliptic(xas)
    end
    e = eccentricity(xas..., μ, 0.0)
    a = semimajoraxis(xas..., μ, 0.0)
    q = a * (1 - e)
    tp = timeperipass(jd, xas..., μ, 0.0)
    n = meanmotion(μ,a)
    M = rad2deg(meananomaly(n, jd, tp))
    Ω = rad2deg(longascnode(xas...))
    ω = rad2deg(argperi(xas..., μ, 0.0))
    i = rad2deg(inclination(xas...))
    return OsculatingElements(e, q, tp, Ω, ω, i, M, a)
end

# Return the cartesian state vector of x at time t [julian days]
function (x::OsculatingElements)(t::Number)
    # Mean motion
    n = PE.meanmotion(μ_S, x.a)
    # Mean anomaly
    M = PE.meananomaly(n, t, x.tp)
    # Eccentric anomaly
    E = PE.eccentricanomaly(x.e, M)
    # True anomaly
    f = PE.trueanomaly(x.e, E)
    # Distance to the central body
    r = x.a * (1 - x.e^2) / (1 + x.e * cos(f))
    # Obtain position and velocity in the orbital frame
    r_o = r .* [cos(f), sin(f), 0.0]
    v_o = (sqrt(μ_S*x.a)/r) .* [-sin(E), sqrt(1 - x.e^2) * cos(E), 0.0]
    # Transform r_o and v_o to the inertial frame
    ω = deg2rad(x.ω)
    i = deg2rad(x.i)
    Ω = deg2rad(x.Ω)
    # Rotation from orbital to inertial frame
    A = Rz(-Ω) * Rx(-i) * Rz(-ω)
    r_i = A * r_o
    v_i = A * v_o
    # State vector
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
