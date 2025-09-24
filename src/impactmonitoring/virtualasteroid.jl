"""
    AbstractVirtualAsteroid{T <: Real} <: AbstractImpactMonitoring

Supertype for the virtual asteroids interface.
"""
abstract type AbstractVirtualAsteroid{T <: Real} <: AbstractImpactMonitoring end

"""
    VirtualAsteroid{T} <: AbstractVirtualAsteroid{T}

A set of [`CloseApproach`](@ref)s with subsequent domains over the line
of variations (LOV) and nominal times within a given time window.

# Fields

- `domain::NTuple{2, T}`: segment of the LOV.
- `CAs::Vector{CloseApproach{T}}`: vector of close approaches.
"""
struct VirtualAsteroid{T} <: AbstractVirtualAsteroid{T}
    domain::NTuple{2, T}
    CAs::Vector{CloseApproach{T}}
end

function show(io::IO, x::VirtualAsteroid)
    t = days2dtutc(nominaltime(x))
    N = length(x.CAs)
    domain = x.domain
    print(io, "VA at ", t, " with ", N, " CAs covering ", domain)
end

in(σ::Real, x::VirtualAsteroid) = x.domain[1] ≤ σ ≤ x.domain[2]

getindex(x::VirtualAsteroid, i::Int) = x.CAs[i]
lastindex(x::VirtualAsteroid) = lastindex(x.CAs)

get_order(x::VirtualAsteroid) = get_order(x[1])

center(x::VirtualAsteroid) = (lbound(x) + ubound(x)) / 2
lbound(x::VirtualAsteroid) = x.domain[1]
ubound(x::VirtualAsteroid) = x.domain[2]

width(x::VirtualAsteroid) = width(x.domain)

nominaltime(x::VirtualAsteroid) = mean(nominaltime, x.CAs)

closeapproaches(x::VirtualAsteroid) = x.CAs

isconvergent(x::VirtualAsteroid, ϵ::Real) = all(Base.Fix2(isconvergent, ϵ), x.CAs)

function convergence_domain(x::VirtualAsteroid, ϵ::Real)
    ds = convergence_domain.(x.CAs, ϵ)
    return (minimum(first, ds), maximum(last, ds))
end

function exponential_weights(x::CloseApproach, σ::Real, ϵ::Real)
    dσ = abs(σ - center(x)) / (convergence_radius(x, ϵ) * domain_radius(x))
    w = 1 / 10^(dσ - 1)
    return w
end

function findsigma(x::VirtualAsteroid, σ::Real, ϵ::Real)
    a, b = convergence_domain(x, ϵ)
    a ≤ σ ≤ ubound(x[1]) && return 1
    lbound(x[end]) ≤ σ ≤ b && return lastindex(x)
    for i in 2:lastindex(x)-1
        lbound(x[i]) ≤ σ ≤ ubound(x[i]) && return i
    end
    return 0
end

for f in (:(timeofca), :(distance), :(rvelea), :(concavity))
    @eval begin
        function $f(x::VirtualAsteroid, σ::Real, ϵ::Real)
            d = convergence_domain(x, ϵ)
            @assert d[1] ≤ σ ≤ d[2] "`σ` is outside the convergence domain of `x`"
            i = findsigma(x, σ, ϵ)
            # Normal evaluation of the Taylor polynomials
            if isconvergent(x[i], ϵ) ||
                (i == 1 && d[1] ≤ σ ≤ center(x[1])) ||
                (i == length(x.CAs) && center(x[end]) ≤ σ ≤ d[2])
                return $f(x[i], σ)
            end
            # Exponentially decaying average
            i = findlast( y -> y.σ < σ, x.CAs)
            j = findfirst(y -> σ < y.σ, x.CAs)
            wi, wj = exponential_weights(x[i], σ, ϵ), exponential_weights(x[j], σ, ϵ)
            return (wi * $f(x[i], σ) + wj * $f(x[j], σ)) / (wi + wj)
        end
    end
end

"""
    virtualasteroids(CAs [, Δt])

Aggregate a vector of [`CloseApproach`](@ref)s into a vector of
[`VirtualAsteroid`](@ref)s given a time window `Δt` (default: `45.0`).
"""
function virtualasteroids(x::AbstractVector{CloseApproach{T}},
                          Δt::Real = 45.0) where {T <: Real}
    VAs = Vector{VirtualAsteroid{T}}(undef, 0)
    y = deepcopy(x)
    while !isempty(y)
        CAs = splice!(y, 1:1)
        for _ in eachindex(y)
            j = findfirst(z -> ubound(z) == lbound(CAs[1]) && difft(z, CAs[1]) < Δt, y)
            isnothing(j) && break
            z = popat!(y, j)
            pushfirst!(CAs, z)
        end
        for _ in eachindex(y)
            j = findfirst(z -> ubound(CAs[end]) == lbound(z) && difft(CAs[end], z) < Δt, y)
            isnothing(j) && break
            z = popat!(y, j)
            push!(CAs, z)
        end
        domain = (minimum(lbound, CAs), maximum(ubound, CAs))
        push!(VAs, VirtualAsteroid{T}(domain, CAs))
    end

    return VAs
end

function virtualasteroids(x::Vector{Vector{CloseApproach{T}}},
                          Δt::Real = 45.0) where {T <: Real}
    CAs = reduce(vcat, x)
    sort!(CAs, by = nominaltime)
    VAs = virtualasteroids(CAs, Δt)
    return VAs
end