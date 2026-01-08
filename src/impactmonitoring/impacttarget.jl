"""
    ImpactTarget{T} <: AbstractImpactTarget{T}

A celestial body that can be hit by an orbit.

# Fields

- `gm::T`: gravitational parameter [au³/day²].
- `radius::T`: physical radius [au].
- `eph::DensePropagation2{T, T}`: barycentric cartesian ephemeris.
"""
struct ImpactTarget{T} <: AbstractImpactTarget{T}
    gm::T
    radius::T
    eph::DensePropagation2{T, T}
end

# Outer constructors
ImpactTarget(gm::T, radius::T, eph::DensePropagation2{T, T}) where {T <: Real} =
    ImpactTarget{T}(gm, radius, eph)

ImpactTarget(i::Int) = ImpactTarget(PE.μ[i], PLANET_RADII[i], _loadeph(i))
ImpactTarget(s::String) = ImpactTarget(PLANET_NAMES_TO_INDEX[s])
ImpactTarget(s::Symbol) = ImpactTarget(string(s))

gm(x::ImpactTarget) = x.gm
radius(x::ImpactTarget) = x.radius

# Print method for ImpactTarget
show(io::IO, x::ImpactTarget) = print(io, "Impact target with gm ",
    @sprintf("%.2E", gm(x)), " au³/day² and radius ", @sprintf("%.2E", radius(x)), " au")

# Evaluation in time method
(x::ImpactTarget)(t) = x.eph(t)