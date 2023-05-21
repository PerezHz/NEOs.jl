@doc raw"""
    OsculatingElements{T <: Number}

Osculating orbital elements of a NEO. 

# Fields 

- `e::T`: eccentricity.
- `q::T`: perihelion distance [au].
- `tp::T`: time of pericenter passage [jd]. 
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
    a::T
    # Inner constructor 
    function OsculatingElements{T}(e::T, q::T, tp::T, Ω::T, ω::T, i::T, a::T) where {T <: Number}
        return new{T}(e, q, tp, Ω, ω, i, a)
    end
end

# Outer constructors
function OsculatingElements(e::T, q::T, tp::T, Ω::T, ω::T, i::T, a::T) where {T <: Number}
    return OsculatingElements{T}(e, q, tp, Ω, ω, i, a)
end

function OsculatingElements()
    return OsculatingElements(NaN, NaN, NaN, NaN, NaN, NaN, NaN)
end

# A OsculatingElements is NaN if all its fields are NaN 
function isnan(osc::OsculatingElements{T}) where {T <: Number} 
    return isnan(osc.e) && isnan(osc.q) && isnan(osc.tp) && isnan(osc.Ω) && isnan(osc.ω) && isnan(osc.i) && isnan(osc.a)
end

# Print method for OsculatingElements
# Example: 
# Semimajor axis (a):                 0.8717319220347314 au
# Eccentricity (e):                   0.4231715487782969
# Time of pericenter passage (tp):    2021-03-19T03:45:01.293 JDTDB
# Pericenter distance (q):            0.5028397744678126 au
# Argument of pericenter (ω):         194.74654283451 deg
# Inclination (i):                    34.81327005431841 deg
# Longitude of Ascending Node (Ω):    17.855086873010706 deg
function show(io::IO, m::OsculatingElements{T}) where {T <: Number} 
    
    print(io, rpad("Semimajor axis (a): ", 36), cte(m.a), " au\n")
    print(io, rpad("Eccentricity (e): ", 36), cte(m.e), "\n")
    if isnan(m.tp)
        print(io, rpad("Time of pericenter passage (tp): ", 36), "NaN JDTDB\n")
    else 
        print(io, rpad("Time of pericenter passage (tp): ", 36), julian2datetime(cte(m.tp)), " JDTDB\n")
    end 
    print(io, rpad("Pericenter distance (q): ", 36), cte(m.q), " au\n")
    print(io, rpad("Argument of pericenter (ω): ", 36), cte(m.ω), " deg\n")
    print(io, rpad("Inclination (i): ", 36), cte(m.i), " deg\n")
    print(io, rpad("Longitude of Ascending Node (Ω): ", 36), cte(m.Ω), " deg\n")

end

@doc raw"""
    equatorial2ecliptic(xas::Vector{T}) where {T <: Number}

Rotate state vector `xas` from equatorial plane to the ecliptic. 
"""
function equatorial2ecliptic(xas::Vector{T}) where {T <: Number}
    # Rotation matrix (only positions)
    m_eq2ecl = Rx(deg2rad(ϵ0_deg))
    # Rotational matrix (positions + velocities)
    m_xv_eq2ecl = hcat(vcat(m_eq2ecl, zeros(3,3)), vcat(zeros(3,3), m_eq2ecl))
    # Rotated state vector 
    return m_xv_eq2ecl*xas
end 

@doc raw"""
    pv2kep(xas, μ = μ_S, jd = JD_J2000, frame::Symbol = :equatorial)

Compute the orbital elements of the NEO with state vector `xas`. Return a `OsculatingElements` object. 

See also [`equatorial2ecliptic`](@ref), [`eccentricity`](@ref), [`semimajoraxis`](@ref), [`timeperipass`](@ref),
[`longascnode`](@ref), [`argperi`](@ref) and [`inclination`](@ref).

# Arguments 

- `xas`: State vector of the asteroid `[x, y, z, v_x, v_y, v_z]`. 
- `μ_S`: Mass parameter of the central body (Sun).
- `jd`: Orbit epoch of reference in julian days. 
- `frame::Symbol`: plane of reference (`:equatorial` or `:ecliptic`).
"""
function pv2kep(xas, μ = μ_S, jd = JD_J2000, frame::Symbol = :equatorial)
    if frame == :ecliptic
        xas = equatorial2ecliptic(xas)
    end 
    e = eccentricity(xas..., μ, 0.0)
    a = semimajoraxis(xas..., μ, 0.0)
    q = a * (1 - e)
    tp = timeperipass(jd, xas..., μ, 0.0)
    Ω = rad2deg(longascnode(xas...))
    ω = rad2deg(argperi(xas..., μ, 0.0))
    i = rad2deg(inclination(xas...))
    return OsculatingElements(e, q, tp, Ω, ω, i, a)
end

@doc raw"""
    (osc::OsculatingElements{T})(t::T) where {T <: Number}

Return cartesian state vector of orbit `osc` at time `t`.
"""
function (osc::OsculatingElements{T})(t::T) where {T <: Number}

    # Mean motion 
    n = PE.meanmotion(μ_S, osc.a)
    # Mean anomaly 
    M = PE.meananomaly(n, t, osc.tp)
    # Eccentric anomaly
    E = PE.eccentricanomaly(osc.e, M)
    # True anomaly
    f = PE.trueanomaly(osc.e, E)
    
    # Distance to the central body 
    r = osc.a * (1 - osc.e^2) / (1 + osc.e * cos(f))
    
    # Obtain position and velocity in the orbital frame 
    r_o = r .* [cos(f), sin(f), 0.0]
    v_o = (sqrt(μ_S*osc.a)/r) .* [-sin(E), sqrt(1 - osc.e^2) * cos(E), 0.0]
    
    # Transform r_o and v_o to the inertial frame 
    ω = deg2rad(osc.ω)
    i = deg2rad(osc.i)
    Ω = deg2rad(osc.Ω)
    # Rotation from orbital to inertial frame 
    A = Rz(-Ω) * Rx(-i) * Rz(-ω)
    r_i = A * r_o
    v_i = A * v_o

    # State vector 
    pv_i = vcat(r_i, v_i)

    return  pv_i
end 

@doc raw"""
    yarkp2adot(A2, a, e, μ_S)

Return the average semimajor axis drift due to the Yarkovsky effect
```math
\begin{align*}
    \left\langle\dot{a}\right\rangle & = \frac{2A_2(1-e^2)}{n}\left(\frac{1 \ \text{AU}}{p}\right)^2 \\
    & = \frac{2A_2}{(1-e^2)\sqrt{a\mu_\odot}}(1 \ \text{AU})^2,
\end{align*}
```
where ``A_2`` is the Yarkovsky parameter, ``\mu_\odot = GM_\odot`` is the Sun's gravitational parameter,
``e`` is the eccentricity, ``n = \sqrt{\mu/a^3}`` is the mean motion, ``p = a(1-e^2)`` is the 
semilatus rectum, and ``a`` is the semimajor axis. 

See https://doi.org/10.1016/j.icarus.2013.02.004.

# Arguments 

- `A2`: Yarkovsky parameter.
- `a`: semimajor axis. 
- `e`: eccentricity. 
- `μ_S`: mass parameter of the Sun. 
"""
function yarkp2adot(A2, a, e, μ_S)
    return 2A2/(sqrt(a)*(1-e^2)*sqrt(μ_S))
end
