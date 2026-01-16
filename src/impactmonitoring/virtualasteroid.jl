"""
    VirtualAsteroid{T} <: AbstractLineOfVariations{T}

A segment of the line of variations (LOV).

# Fields

- `epoch::T`: reference epoch [days since J2000 TDB].
- `σ::T`: LOV index.
- `domain::NTuple{2, T}`: segment of the line of variations.
- `q0::Vector{Taylor1{T}}`: barycentric initial condition [au, au/day].
"""
struct VirtualAsteroid{T} <: AbstractLineOfVariations{T}
    epoch::T
    σ::T
    domain::NTuple{2, T}
    q0::Vector{Taylor1{T}}
end

# Outer constructors
function VirtualAsteroid(epoch::T, σ::T, domain::NTuple{2, T},
                         q0::Vector{Taylor1{T}}) where {T <: Real}
    return VirtualAsteroid{T}(epoch, σ, domain, q0)
end

VirtualAsteroid(lov::LineOfVariations, σ::Real, domain::NTuple{2, <:Real}) =
    VirtualAsteroid(epoch(lov), σ, domain, lov(σ, domain))

# Print method for VirtualAsteroid
function show(io::IO, x::VirtualAsteroid)
    t = date(x)
    domain = x.domain
    print(io, "VA at ", t, " over ", domain)
end

# AbstractLineOfVariations interface
epoch(x::VirtualAsteroid) = x.epoch
nominaltime(x::VirtualAsteroid) = x.epoch
sigma(x::VirtualAsteroid) = x.σ
initialcondition(x::VirtualAsteroid) = x.q0
get_order(x::VirtualAsteroid) = get_order(first(x.q0))

"""
    virtualasteroids(lov, method; kwargs...)

Divide the line of variations `lov` into a vector of virtual asteroids
using `method`. Each of the following available methods have their own
keyword arguments:
- `:uniform`: domains of equal length. Keyword arguments:
    - `N::Int`: number of domains.
- `:normal`: domains of equal probability according to a normal density
    function. Keyword arguments:
    - `N::Int`: number of domains.
- `:DelVigna19`: optimal method for LOV sampling proposed by Del Vigna
    et al (2019). Keyword arguments:
    - `R_TP::Real`: target plane radius [au].
    - `R_P::Real`: planet radius [au].
    - `IP::Real`: generic completeness limit.
    - `Δσmax`: maximum interval length.
"""
virtualasteroids(lov::LineOfVariations, method::Symbol; kwargs...) =
    virtualasteroids(Val(method), lov; kwargs...)

function virtualasteroids(::Val{:uniform}, lov::LineOfVariations; N::Int)
    endpoints = LinRange(lbound(lov), ubound(lov), N+1)
    domains = [(endpoints[i], endpoints[i+1]) for i in 1:N]
    σs = midpoint.(domains)
    if isodd(N)
        σs[(N ÷ 2) + 1] = sigma(lov)
    end
    return VirtualAsteroid.(Ref(lov), σs, domains)
end

function virtualasteroids(::Val{:normal}, lov::LineOfVariations; N::Int)
    _, T = numtypes(lov)
    d = Normal{T}(zero(T), one(T))
    xs = LinRange(cdf(d, lbound(lov)), cdf(d, ubound(lov)), N+1)
    endpoints = quantile.(d, xs)
    endpoints[1], endpoints[end] = lbound(lov), ubound(lov)
    domains = [(endpoints[i], endpoints[i+1]) for i in 1:N]
    σs = midpoint.(domains)
    if isodd(N)
        σs[(N ÷ 2) + 1] = sigma(lov)
    end
    return VirtualAsteroid.(Ref(lov), σs, domains)
end

function virtualasteroids(::Val{:DelVigna19}, lov::LineOfVariations;
                          R_TP::Real, R_P::Real, IP::Real, Δσmax::Real)
    _, T = numtypes(lov)
    endpoints = [zero(T)]
    while endpoints[end] < ubound(lov)
        σ = endpoints[end]
        σ += min(R_TP * IP / (2R_P * lovdensity(σ)), Δσmax)
        if σ < ubound(lov)
            push!(endpoints, σ)
        else
            push!(endpoints, ubound(lov))
            break
        end
    end
    endpoints = vcat(reverse(-endpoints[2:end]), endpoints)
    N = length(endpoints) - 1
    domains = [(endpoints[i], endpoints[i+1]) for i in 1:N]
    σs = midpoint.(domains)
    if isodd(N)
        σs[(N ÷ 2) + 1] = sigma(lov)
    end
    return VirtualAsteroid.(Ref(lov), σs, domains)
end