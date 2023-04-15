@doc raw"""
    OsculatingElements{T <: AbstractFloat}

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
@auto_hash_equals struct OsculatingElements{T <: AbstractFloat}
    e::T
    q::T
    tp::T
    Ω::T
    ω::T
    i::T
    a::T
    # Inner constructor 
    function OsculatingElements{T}(e::T, q::T, tp::T, Ω::T, ω::T, i::T, a::T) where {T <: AbstractFloat}
        return new{T}(e, q, tp, Ω, ω, i, a)
    end
end

# Outer constructors
function OsculatingElements(e::T, q::T, tp::T, Ω::T, ω::T, i::T, a::T) where {T <: AbstractFloat}
    return OsculatingElements{T}(e, q, tp, Ω, ω, i, a)
end

function OsculatingElements()
    return OsculatingElements(NaN, NaN, NaN, NaN, NaN, NaN, NaN)
end

# A OsculatingElements is NaN if all its fields are NaN 
function isnan(osc::OsculatingElements{T}) where {T <: AbstractFloat} 
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
function show(io::IO, m::OsculatingElements{T}) where {T <: AbstractFloat} 
    
    print(io, rpad("Semimajor axis (a): ", 36), m.a, " au\n")
    print(io, rpad("Eccentricity (e): ", 36), m.e, "\n")
    if isnan(m.tp)
        print(io, rpad("Time of pericenter passage (tp): ", 36), "NaN JDTDB\n")
    else 
        print(io, rpad("Time of pericenter passage (tp): ", 36), julian2datetime(m.tp), " JDTDB\n")
    end 
    print(io, rpad("Pericenter distance (q): ", 36), m.q, " au\n")
    print(io, rpad("Argument of pericenter (ω): ", 36), m.ω, " deg\n")
    print(io, rpad("Inclination (i): ", 36), m.i, " deg\n")
    print(io, rpad("Longitude of Ascending Node (Ω): ", 36), m.Ω, " deg\n")

end

@doc raw"""
    pv2kep(xas, μ_S = μ_S, jd = JD_J2000)

Compute the orbital elements of the NEO with state vector `xas`. Return a `OsculatingElements` object. 

See also [`eccentricity`](@ref), [`semimajoraxis`](@ref), [`timeperipass`](@ref),
[`longascnode`](@ref), [`argperi`](@ref) and [`inclination`](@ref).

# Arguments 

- `xas`: State vector of the asteroid `[x, y, z, v_x, v_y, v_z]`. 
- `μ_S`: Mass parameter of the Sun.
- `jd`: Julian days since J2000. 
"""
function pv2kep(xas, μ_S = μ_S, jd = JD_J2000)
    e = eccentricity(xas..., μ_S, 0.0)
    a = semimajoraxis(xas..., μ_S, 0.0)
    q = a * (1 - e)
    tp = timeperipass(jd, xas..., μ_S, 0.0)
    Ω = rad2deg(longascnode(xas...))
    ω = rad2deg(argperi(xas..., μ_S, 0.0))
    i = rad2deg(inclination(xas...))
    return OsculatingElements(e, q, tp, Ω, ω, i, a)
end

function mean(osc::Vector{OsculatingElements{T}}) where {T <: AbstractFloat}
    m = length(osc)

    e = zero(T)
    q = zero(T)
    tp = zero(T)
    Ω = zero(T)
    ω = zero(T)
    i = zero(T)
    a = zero(T)

    for j in 1:m
        e += osc[j].e
        q += osc[j].q
        tp += osc[j].tp
        Ω += osc[j].Ω
        ω += osc[j].ω
        i += osc[j].i
        a += osc[j].a
    end

    return OsculatingElements(e/m, q/m, tp/m, Ω/m, ω/m, i/m, a/m)
end

function (osc::OsculatingElements{T})(t::T) where {T <: AbstractFloat}

    # Mean motion 
    n = PlanetaryEphemeris.meanmotion(μ_S, osc.a)
    # Mean anomaly 
    M = PlanetaryEphemeris.meananomaly(n, t, osc.tp)
    # Eccentric anomaly
    E = PlanetaryEphemeris.eccentricanomaly(osc.e, M)
    # True anomaly
    f = PlanetaryEphemeris.trueanomaly(osc.e, E)
    
    # Distance to the central body 
    r = osc.a * (1 - osc.e^2) / (1 + osc.e * cos(f))
    
    # Obtain position and velocity in the orbital frame 
    r_o = r .* [cos(f), sin(f), 0.0]
    v_o = (sqrt(μ_S*osc.a)/r) .* [-sin(E), sqrt(1 - osc.e^2) * cos(E), 0.0]
    
    # Transform r_o and v_o to the inertial frame 
    ω = deg2rad(osc.ω)
    i = deg2rad(osc.i)
    Ω = deg2rad(osc.Ω)

    A = Rz(-Ω) * Rx(-i) * Rz(-ω)
    r_i = A * r_o
    v_i = A * v_o

    pv_i = vcat(r_i, v_i)

    # Barycentric state vector of the sun 
    et = julian2etsecs(t)
    sun_bar = kmsec2auday(sun_pv(et))

    return  pv_i .+ sun_bar
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