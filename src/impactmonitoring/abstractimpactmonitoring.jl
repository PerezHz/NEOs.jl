"""
    AbstractImpactMonitoring

Supertye for the impact monitoring interface.
"""
abstract type AbstractImpactMonitoring end

"""
    AbstractTargetPlane{U <: Number} <: AbstractImpactMonitoring

Supertype for the target planes interface.

For every target plane `x`, `targetplane(x)` returns a 3-element vector
containing two coordinates on the target plane and the planet's impact
cross section.
"""
abstract type AbstractTargetPlane{U <: Number} <: AbstractImpactMonitoring end

numtype(::AbstractTargetPlane{U}) where {U} = U

# Print method for AbstractTargetPlane
show(io::IO, x::AbstractTargetPlane) = print(io, typeof(x), " with coordinates ",
    cte(targetplane(x)))

"""
    AbstractLineOfVariations{T <: Real} <: AbstractImpactMonitoring

Supertype for the line of variations nterface.

Every instance `x` of `AbstractLineOfVariations` has a:
- `date(x)`: nominal date [UTC].
- `sigma(x)`: center of expansion.
- `lbound(x)`: domain lower bound.
- `ubound(x)`: domain upper bound.
"""
abstract type AbstractLineOfVariations{T <: Real} <: AbstractImpactMonitoring end

date(x::AbstractLineOfVariations) = days2dtutc(nominaltime(x))

lbound(x::AbstractLineOfVariations) = x.domain[1]
ubound(x::AbstractLineOfVariations) = x.domain[2]

in(σ::Real, x::AbstractLineOfVariations) = lbound(x) ≤ σ ≤ ubound(x)

width(x::NTuple{2, T}) where {T <: Real} = x[2] - x[1]
width(x::AbstractLineOfVariations) = ubound(x) - lbound(x)

"""
    AbstractVirtualImpactor{T <: Real} <: AbstractImpactMonitoring

Supertype for the virtual impactors interface.

Every instance `x` of `AbstractVirtualImpactor` has a:
- `date(x)`: nominal date [UTC].
- `sigma(x)`: sigma parameter, whose interpretation depends
    on the type of virtual impactor.
- `impact_probability(x)`: impact probability.
"""
abstract type AbstractVirtualImpactor{T <: Real} <: AbstractImpactMonitoring end

overlap(a::NTuple{2, T}, b::NTuple{2, T}) where {T <: Real} =
    (a[1] ≤ b[2]) && (b[1] ≤ a[2])

function show(io::IO, x::AbstractVirtualImpactor)
    d = round(date(x), Minute)
    t = Dates.format(d, "yyyy-mm-dd HH:MM")
    ip = @sprintf("%.2E", impact_probability(x))
    asterisk = isoutlov(x) ? " *" : ""
    print(io, "VI at ", t, " with probability ", ip, asterisk)
end
