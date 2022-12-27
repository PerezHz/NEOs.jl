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
struct OsculatingElements{T <: AbstractFloat}
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

# Print method fot OsculatingElements
# Example: 
# Semimajor axis (a):                 1.022714645967277 au
# Eccentricity (e):                   0.18919435448507194
# Time of pericenter passage (tp):    2.4592910085352035e6 JDTDB
# Pericenter distance (q):            0.8071905529640956 au
# Argument of pericenter (ω):         251.59066974220704 deg
# Inclination (i):                    22.56425964296357 deg
# Longitude of Ascending Node (Ω):    233.31225377038814 deg
function show(io::IO, m::OsculatingElements{T}) where {T <: AbstractFloat} 
    
    print(io, rpad("Semimajor axis (a): ", 36), m.a, " au\n")
    print(io, rpad("Eccentricity (e): ", 36), m.e, "\n")
    print(io, rpad("Time of pericenter passage (tp): ", 36), m.tp, " JDTDB\n")
    print(io, rpad("Pericenter distance (q): ", 36), m.q, " au\n")
    print(io, rpad("Argument of pericenter (ω): ", 36), m.ω, " deg\n")
    print(io, rpad("Inclination (i): ", 36), m.i, " deg\n")
    print(io, rpad("Longitude of Ascending Node (Ω): ", 36), m.Ω, " deg\n")

end

@doc raw"""
    pv2kep(xas, μ_S = μ_S, jd = JD_J2000)

Computes the orbital elements of the NEO with state vector `xas`. Returns a `OsculatingElements` object. 

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

function atan2(y, x)
    if x > 0
        return atan(y/x)
    elseif (y ≥ 0) && (x < 0)
        return atan(y/x) + π
    elseif  (y < 0) && (x < 0)
        return atan(y/x) - π
    elseif (y > 0) && (x == 0)
        return π/2
    elseif (y < 0) && (x == 0)
        return -π/2
    elseif (y == 0) && (x == 0)
        return NaN
    end 
end

function (osc::OsculatingElements{T})(t::T) where {T <: AbstractFloat}

    f = PlanetaryEphemeris.time2truean(osc.a, osc.e, μ_S, t, osc.tp)
    
    # 4.- Get distance to the central body 
    r = osc.a * (1 - osc.e^2) / (1 + osc.e * cos(f))
    
    # 5.- Obtain position and velocity in the orbital frame 
    r_o = r .* [cos(f), sin(f), 0.0]
    #v_o = (sqrt(μ_S*osc.a)/r) .* [-sin(E), sqrt(1 - osc.e^2) * cos(E), 0.0]
    
    # 6.- Transform r_o and v_o to the inertial frame 
    ω = deg2rad(osc.ω)
    i = deg2rad(osc.i)
    Ω = deg2rad(osc.Ω)

    A = Rz(-Ω) * Rx(-i) * Rz(-ω)
    r_i = A * r_o
    #v_i = A * v_o

    return vcat(r_i)#, v_i)
end 