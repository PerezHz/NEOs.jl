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

gm(x::ImpactTarget) = x.gm
radius(x::ImpactTarget) = x.radius

# Evaluation in time method
(x::ImpactTarget)(t) = x.eph(t)