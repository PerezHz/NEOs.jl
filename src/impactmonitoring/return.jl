"""
    Return{T, U <: Number} <: AbstractLineOfVariations{T}

A set of virtual asteroids with contiguous LOV indices that have a close
approach within the same time window (shower).

# Fields

- `domain::NTuple{2, T}`: segment of the line of variations.
- `CAs::Vector{CloseApproach{T, U}}`: vector of close approaches.
"""
@auto_hash_equals struct Return{T, U <: Number} <: AbstractLineOfVariations{T}
    domain::NTuple{2, T}
    CAs::Vector{CloseApproach{T, U}}
end

# Abbreviations
const ReturnT1{T} = Return{T, Taylor1{T}}

numtypes(::Return{T, U}) where {T, U} = T, U

# Print method for Return
function show(io::IO, x::Return)
    T, U = numtypes(x)
    t = date(x)
    N = length(x.CAs)
    domain = x.domain
    print(io, "Return{$T, $U} at ", t, " with ", N, " VAs covering ", domain)
end

# AbstractLineOfVariations interface
closeapproaches(x::Return) = x.CAs
sigma(x::Return) = (lbound(x) + ubound(x)) / 2
nominaltime(x::Return) = mean(nominaltime, x.CAs)

get_order(x::Return{T, U}) where {T, U <: AbstractSeries} = get_order(x[1])

getindex(x::Return, i::Int) = x.CAs[i]
getindex(x::Return, idxs::AbstractVector{Int}) = x.CAs[idxs]
firstindex(x::Return) = firstindex(x.CAs)
lastindex(x::Return) = lastindex(x.CAs)

isconvergent(x::Return, ctol::Real) = all(Base.Fix2(isconvergent, ctol), x.CAs)

function convergence_domain(x::Return, ctol::Real)
    ds = convergence_domain.(x.CAs, ctol)
    return (minimum(first, ds), maximum(last, ds))
end

function exponential_weights(x::CloseApproach, σ::Real, ctol::Real)
    dσ = deltasigma(x, σ) / convergence_radius(x, ctol)
    w = 1 / 1000^(abs(dσ) - 1)
    return w
end

for f in (:(targetplane), :(timeofca), :(distance), :(radialvelocity), :(concavity))
    @eval begin
        function $f(x::ReturnT1, σ::Real, ctol::Real)
            d = convergence_domain(x, ctol)
            @assert d[1] ≤ σ ≤ d[2] "`σ` is outside the convergence domain of `x`"
            ws = @. exponential_weights(x.CAs, σ, ctol)
            fs = @. $f(x.CAs, σ)
            return mean(fs, weights(ws))
        end
    end
end

nominalstate(x::ReturnT1, ctol::Real) = targetplane(x, sigma(x), ctol)

issamereturn(x::CloseApproach, y::CloseApproach, Δt::Real) = (x.tp == y.tp) &&
    ubound(x) == lbound(y) && difft(x, y) < Δt

"""
    showersnreturns(CAs [, Δt])

Aggregate a vector of close approaches into a vector of returns
given a time window `Δt` [days] (default: `45.0`).
"""
function showersnreturns(x::AbstractVector{CloseApproach{T, U}},
                         Δt::Real = 45.0) where {T <: Real, U <: Number}
    RTs = Vector{Return{T, U}}(undef, 0)
    y = deepcopy(x)
    while !isempty(y)
        CAs = splice!(y, 1:1)
        for _ in eachindex(y)
            j = findfirst(z -> issamereturn(z, CAs[1], Δt), y)
            isnothing(j) && break
            z = popat!(y, j)
            pushfirst!(CAs, z)
        end
        for _ in eachindex(y)
            j = findfirst(z -> issamereturn(CAs[end], z, Δt), y)
            isnothing(j) && break
            z = popat!(y, j)
            push!(CAs, z)
        end
        domain = (minimum(lbound, CAs), maximum(ubound, CAs))
        push!(RTs, Return{T, U}(domain, CAs))
    end

    return RTs
end

function showersnreturns(x::Vector{Vector{CloseApproach{T, U}}},
                         Δt::Real = 45.0) where {T <: Real, U <: Number}
    CAs = reduce(vcat, x)
    sort!(CAs, by = nominaltime)
    RTs = showersnreturns(CAs, Δt)
    return RTs
end