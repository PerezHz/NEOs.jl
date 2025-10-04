"""
    VirtualAsteroid{T} <: AbstractLineOfVariations{T}

A set of [`CloseApproach`](@ref)s with subsequent domains over the line
of variations and nominal dates within a given time window.

# Fields

- `domain::NTuple{2, T}`: segment of the line of variations.
- `CAs::Vector{CloseApproach{T}}`: vector of close approaches.
"""
@auto_hash_equals struct VirtualAsteroid{T} <: AbstractLineOfVariations{T}
    domain::NTuple{2, T}
    CAs::Vector{CloseApproach{T}}
end

function show(io::IO, x::VirtualAsteroid)
    t = date(x)
    N = length(x.CAs)
    domain = x.domain
    print(io, "VA at ", t, " with ", N, " CAs covering ", domain)
end

nominaltime(x::VirtualAsteroid) = mean(nominaltime, x.CAs)

sigma(x::VirtualAsteroid) = (lbound(x) + ubound(x)) / 2

get_order(x::VirtualAsteroid) = get_order(x[1])

getindex(x::VirtualAsteroid, i::Int) = x.CAs[i]
lastindex(x::VirtualAsteroid) = lastindex(x.CAs)

closeapproaches(x::VirtualAsteroid) = x.CAs

isconvergent(x::VirtualAsteroid, ctol::Real) = all(Base.Fix2(isconvergent, ctol), x.CAs)

function convergence_domain(x::VirtualAsteroid, ctol::Real)
    a, _ = convergence_domain(x[1], ctol)
    _, b = convergence_domain(x[end], ctol)
    return (a, b)
end

function exponential_weights(x::CloseApproach, σ::Real, ctol::Real)
    dσ = deltasigma(x, σ) / convergence_radius(x, ctol)
    w = 1 / 10^(abs(dσ) - 1)
    return w
end

function findsigma(x::VirtualAsteroid, σ::Real, ctol::Real)
    a, b = convergence_domain(x, ctol)
    a ≤ σ ≤ ubound(x[1]) && return 1
    lbound(x[end]) ≤ σ ≤ b && return lastindex(x)
    for i in 2:lastindex(x)-1
        lbound(x[i]) ≤ σ ≤ ubound(x[i]) && return i
    end
    return 0
end

nominalstate(x::VirtualAsteroid, ctol::Real) = targetplane(x, sigma(x), ctol)

for f in (:(targetplane), :(timeofca), :(distance), :(rvelea), :(concavity))
    @eval begin
        function $f(x::VirtualAsteroid, σ::Real, ctol::Real)
            d = convergence_domain(x, ctol)
            @assert d[1] ≤ σ ≤ d[2] "`σ` is outside the convergence domain of `x`"
            i = findsigma(x, σ, ctol)
            # Normal evaluation of the Taylor polynomials
            if isconvergent(x[i], ctol) ||
                (i == 1 && d[1] ≤ σ ≤ sigma(x[1])) ||
                (i == length(x.CAs) && sigma(x[end]) ≤ σ ≤ d[2])
                return $f(x[i], σ)
            end
            # Exponentially decaying average
            i = findlast( y -> sigma(y) < σ, x.CAs)
            j = findfirst(y -> σ < sigma(y), x.CAs)
            wi, wj = exponential_weights(x[i], σ, ctol), exponential_weights(x[j], σ, ctol)
            return (wi * $f(x[i], σ) + wj * $f(x[j], σ)) / (wi + wj)
        end
    end
end

issameva(x::CloseApproach, y::CloseApproach, Δt::Real) = (x.tp == y.tp) &&
    ubound(x) == lbound(y) && difft(x, y) < Δt

"""
    virtualasteroids(CAs [, Δt])

Aggregate a vector of [`CloseApproach`](@ref)s into a vector of
[`VirtualAsteroid`](@ref)s given a time window `Δt` [days] (default: `45.0`).
"""
function virtualasteroids(x::AbstractVector{CloseApproach{T}},
                          Δt::Real = 45.0) where {T <: Real}
    VAs = Vector{VirtualAsteroid{T}}(undef, 0)
    y = deepcopy(x)
    while !isempty(y)
        CAs = splice!(y, 1:1)
        for _ in eachindex(y)
            j = findfirst(z -> issameva(z, CAs[1], Δt), y)
            isnothing(j) && break
            z = popat!(y, j)
            pushfirst!(CAs, z)
        end
        for _ in eachindex(y)
            j = findfirst(z -> issameva(CAs[end], z, Δt), y)
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