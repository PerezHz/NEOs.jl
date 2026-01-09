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
midpoint(x::NTuple{2, <:Real}) = (x[1] + x[2]) / 2

function ceilodd(x::Real)
    N = ceil(Int, x)
    return N + iseven(N)
end

refine(σ::T, domain::NTuple{2, T}, N::Int, dist::Symbol = :mixed) where {T <: Real} =
    refine(Val(dist), σ, domain, N)

function refine(::Val{:uniform}, σ::T, domain::NTuple{2, T}, N::Int) where {T <: Real}
    endpoints = LinRange(domain[1], domain[2], N+1)
    domains = [(endpoints[i], endpoints[i+1]) for i in 1:N]
    σs = midpoint.(domains)
    if isodd(N)
        σs[(N ÷ 2) + 1] = σ
    end
    return σs, domains
end

function refine(::Val{:normal}, σ::T, domain::NTuple{2, T}, N::Int) where {T <: Real}
    d = Normal{T}(zero(T), one(T))
    xs = LinRange(cdf(d, domain[1]), cdf(d, domain[2]), N+1)
    endpoints = quantile.(d, xs)
    endpoints[1], endpoints[end] = domain[1], domain[2]
    domains = [(endpoints[i], endpoints[i+1]) for i in 1:N]
    σs = midpoint.(domains)
    if iszero(σ) && isodd(N)
        σs[(N ÷ 2) + 1] = σ
    end
    return σs, domains
end

function refine(::Val{:mixed}, σ::T, domain::NTuple{2, T}, N::Int) where {T <: Real}
    σs, domains = refine(Val(:normal), σ, domain, N)
    Δσmax = width(domain) / N
    mask = @. σs < σ && width(domains) > Δσmax
    idxs = findfirst(mask):findlast(mask)
    subdomain = (first(domains[first(idxs)]), last(domains[last(idxs)]))
    subσs, subdomains = refine(midpoint(subdomain), subdomain,
        ceilodd(width(subdomain) / Δσmax), :uniform)
    splice!(σs, idxs, subσs)
    splice!(domains, idxs, subdomains)
    mask = @. σs > σ && width(domains) > Δσmax
    idxs = findfirst(mask):findlast(mask)
    subdomain = (first(domains[first(idxs)]), last(domains[last(idxs)]))
    subσs, subdomains = refine(midpoint(subdomain), subdomain,
        ceilodd(width(subdomain) / Δσmax), :uniform)
    splice!(σs, idxs, subσs)
    splice!(domains, idxs, subdomains)
    return σs, domains
end

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
