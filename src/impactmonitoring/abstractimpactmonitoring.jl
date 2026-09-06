"""
    AbstractImpactMonitoring

Supertye for the impact monitoring interface.
"""
abstract type AbstractImpactMonitoring end

"""
    AbstractImpactTarget{T <: Real} <: AbstractImpactMonitoring

Supertye for the impact targets interface.
"""
abstract type AbstractImpactTarget{T <: Real} <: AbstractImpactMonitoring end

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
    AbstractIMProblem{D, T <: Real} <: AbstractImpactMonitoring

Supertye for the impact monitoring problems interface.
"""
abstract type AbstractIMProblem{D, T <: Real} <: AbstractImpactMonitoring end

"""
    AbstractLineOfVariations{T <: Real} <: AbstractImpactMonitoring

Supertype for the line of variations (LOV) interface.

Every instance `x` of `AbstractLineOfVariations` has a:
- `nominaltime(x)`: nominal time [days since J2000 TDB].
- `sigma(x)`: LOV index.
- `domain::NTuple{2, T}` field.
"""
abstract type AbstractLineOfVariations{T <: Real} <: AbstractImpactMonitoring end

date(x::AbstractLineOfVariations) = days2dtutc(nominaltime(x))

lbound(x::AbstractLineOfVariations) = x.domain[1]
ubound(x::AbstractLineOfVariations) = x.domain[2]

in(σ::Real, x::AbstractLineOfVariations) = lbound(x) ≤ σ ≤ ubound(x)

width(x::NTuple{2, <:Real}) = x[2] - x[1]
width(x::AbstractLineOfVariations) = ubound(x) - lbound(x)
midpoint(x::NTuple{2, <:Real}) = (x[1] + x[2]) / 2

lovdensity(x::Real) = exp(-x^2/2) / sqrt(2π)

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

isspurious(::AbstractVirtualImpactor) = false

# Print methods for AbstractVirtualImpactor
function show(io::IO, ::MIME"text/plain", x::AbstractVirtualImpactor)
    t = repeat(' ', 4)
    f = isspurious(x) ? "true " : "false"
    d = Dates.format(round(date(x), Minute), "yyyy-mm-dd HH:MM")
    σ = @sprintf("%+.4f", sigma(x))
    ip = @sprintf("%.2E", impact_probability(x))
    print(io,
        typeof(x), "\n",
        t, rpad("Spurious:", 21), f, "\n",
        t, rpad("Date:", 21), d, "\n",
        t, rpad("Sigma:", 21), σ, "\n",
        t, rpad("Impact probability:", 21), ip,
    )
    return nothing
end

function show(io::IO, x::AbstractVirtualImpactor)
    f = isspurious(x) ? "[Spurious] " : ""
    d = Dates.format(round(date(x), Minute), "yyyy-mm-dd HH:MM")
    σ = @sprintf("%+.4f", sigma(x))
    ip = @sprintf("%.2E", impact_probability(x))
    print(io, f, "VI t: ", d, " σ: ", σ, " ip: ", ip)
end
