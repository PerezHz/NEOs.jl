struct VirtualAsteroid{T <: Real}
    domain::NTuple{2, T}
    CAs::Vector{CloseApproach{T}}
end

function show(io::IO, x::VirtualAsteroid)
    N = length(x.CAs)
    domain = x.domain
    print(io, "Virtual asteroid with ", N, " segments covering ", domain)
end

in(σ::Real, x::VirtualAsteroid) = x.domain[1] ≤ σ ≤ x.domain[2]
getindex(x::VirtualAsteroid, i::Int) = x.CAs[i]
lastindex(x::VirtualAsteroid) = lastindex(x.CAs)

lbound(x::VirtualAsteroid) = x.domain[1]
ubound(x::VirtualAsteroid) = x.domain[2]

width(x::NTuple{2, T}) where {T <: Real} = x[2] - x[1]
width(x::VirtualAsteroid) = width(x.domain)

isconvergent(x::VirtualAsteroid, ϵ::Real) = all(Base.Fix2(isconvergent, ϵ), x.CAs)

function convergence_domain(x::VirtualAsteroid, ϵ::Real)
    # Inner convergence domain
    ds = convergence_domain.(x.CAs, ϵ)
    a, b = minimum(first, ds), maximum(last, ds)
    # Decide whether to extend to the left
    dσ = max(-1.0, -convergence_radius(x[1], ϵ))
    σ = center(x[1]) + dσ * domain_radius(x[1])
    ξ, ζ, r = targetplane(x[1], σ)
    d2 = (ξ^2 + ζ^2) - r^2
    dξ, dζ, dr = targetplanederivatives(x[1], σ)
    v2 = 2 * (ξ*dξ + ζ*dζ - r*dr)
    if v2 > 0
        if isconvergent(x[1], ϵ) && (-2 ≤ dσ - d2 / v2 ≤ -1)
            a = center(x[1]) - 2 * domain_radius(x[1])
        elseif !isconvergent(x[1], ϵ) && (-1 ≤ dσ - d2 / v2 ≤ dσ)
            a = lbound(x[1])
        end
    end
    # Decide whether to extend to the right
    dσ = min(1.0, convergence_radius(x[end], ϵ))
    σ = center(x[end]) + dσ * domain_radius(x[end])
    ξ, ζ, r = targetplane(x[end], σ)
    d2 = (ξ^2 + ζ^2) - r^2
    dξ, dζ, dr = targetplanederivatives(x[end], σ)
    v2 = 2 * (ξ*dξ + ζ*dζ - r*dr)
    if v2 < 0
        if isconvergent(x[end], ϵ) && (1 ≤ dσ - d2 / v2 ≤ 2)
            b = center(x[end]) + 2 * domain_radius(x[end])
        elseif !isconvergent(x[end], ϵ) && (dσ ≤ dσ - d2 / v2 ≤ 1)
            b = ubound(x[end])
        end
    end
    return (a, b)
end

closeapproaches(x::VirtualAsteroid) = x.CAs

meantime(x::VirtualAsteroid) = mean(nominaltime, x.CAs)

function exponential_weights(a::CloseApproach, b::CloseApproach, σ::Real, ϵ::Real)
    dσa = abs(σ - center(a)) / (convergence_radius(a, ϵ) * domain_radius(a))
    dσb = abs(σ - center(b)) / (convergence_radius(b, ϵ) * domain_radius(b))
    wa = 1 / 10^(dσa - 1)
    wb = 1 / 10^(dσb - 1)
    return wa, wb
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

function time(x::VirtualAsteroid, σ::Real, ϵ::Real)
    d = convergence_domain(x, ϵ)
    @assert d[1] ≤ σ ≤ d[2] "`σ` is outside the convergence domain domain of `x`"
    i = findsigma(x, σ, ϵ)
    # Normal evaluation of the Taylor polynomials
    if isconvergent(x.CAs[i], ϵ) ||
        (i == 1 && d[1] ≤ σ ≤ x.CAs[1].σ) ||
        (i == length(x.CAs) && x.CAs[end].σ ≤ σ ≤ d[2])
        return time(x[i], σ)
    end
    # Exponentially decaying average
    i = findlast( y -> y.σ < σ, x.CAs)
    j = findfirst(y -> σ < y.σ, x.CAs)
    wi, wj = exponential_weights(x.CAs[i], x.CAs[j], σ, ϵ)
    return (wi * time(x[i], σ) + wj * time(x[j], σ)) / (wi + wj)
end

function distance(x::VirtualAsteroid, σ::Real, ϵ::Real)
    d = convergence_domain(x, ϵ)
    @assert d[1] ≤ σ ≤ d[2] "`σ` is outside the convergence domain domain of `x`"
    i = findsigma(x, σ, ϵ)
    # Normal evaluation of the Taylor polynomials
    if isconvergent(x.CAs[i], ϵ) ||
        (i == 1 && d[1] ≤ σ ≤ x.CAs[1].σ) ||
        (i == length(x.CAs) && x.CAs[end].σ ≤ σ ≤ d[2])
        return distance(x[i], σ)
    end
    # Exponentially decaying average
    i = findlast( y -> y.σ < σ, x.CAs)
    j = findfirst(y -> σ < y.σ, x.CAs)
    wi, wj = exponential_weights(x.CAs[i], x.CAs[j], σ, ϵ)
    return (wi * distance(x[i], σ) + wj * distance(x[j], σ)) / (wi + wj)
end

function rvelea(x::VirtualAsteroid, σ::Real, ϵ::Real)
    d = convergence_domain(x, ϵ)
    @assert d[1] ≤ σ ≤ d[2] "`σ` is outside the convergence domain domain of `x`"
    i = findsigma(x, σ, ϵ)
    # Normal evaluation of the Taylor polynomials
    if isconvergent(x.CAs[i], ϵ) ||
        (i == 1 && d[1] ≤ σ ≤ x.CAs[1].σ) ||
        (i == length(x.CAs) && x.CAs[end].σ ≤ σ ≤ d[2])
        return rvelea(x[i], σ)
    end
    # Exponentially decaying average
    i = findlast( y -> y.σ < σ, x.CAs)
    j = findfirst(y -> σ < y.σ, x.CAs)
    wi, wj = exponential_weights(x.CAs[i], x.CAs[j], σ, ϵ)
    return (wi * rvelea(x[i], σ) + wj * rvelea(x[j], σ)) / (wi + wj)
end

function concavity(x::VirtualAsteroid, σ::Real, ϵ::Real)
    d = convergence_domain(x, ϵ)
    @assert d[1] ≤ σ ≤ d[2] "`σ` is outside the convergence domain domain of `x`"
    i = findsigma(x, σ, ϵ)
    # Normal evaluation of the Taylor polynomials
    if isconvergent(x.CAs[i], ϵ) ||
        (i == 1 && d[1] ≤ σ ≤ x.CAs[1].σ) ||
        (i == length(x.CAs) && x.CAs[end].σ ≤ σ ≤ d[2])
        return concavity(x.CAs[i], σ)
    end
    # Exponentially decaying average
    i = findlast( y -> y.σ < σ, x.CAs)
    j = findfirst(y -> σ < y.σ, x.CAs)
    wi, wj = exponential_weights(x.CAs[i], x.CAs[j], σ, ϵ)
    return (wi * concavity(x[i], σ) + wj * concavity(x[j], σ)) / (wi + wj)
end

function distance_derivatives(x::VirtualAsteroid, σ::Real, ϵ::Real)
    d = convergence_domain(x, ϵ)
    @assert d[1] ≤ σ ≤ d[2] "`σ` is outside the convergence domain domain of `x`"
    i = findsigma(x, σ, ϵ)
    # Normal evaluation of the Taylor polynomials
    if isconvergent(x.CAs[i], ϵ) ||
        (i == 1 && d[1] ≤ σ ≤ x.CAs[1].σ) ||
        (i == length(x.CAs) && x.CAs[end].σ ≤ σ ≤ d[2])
        return distance_derivatives(x[i], σ)
    end
    # Exponentially decaying average
    i = findlast( y -> y.σ < σ, x.CAs)
    j = findfirst(y -> σ < y.σ, x.CAs)
    wi, wj = exponential_weights(x.CAs[i], x.CAs[j], σ, ϵ)
    di, ri, ci = distance_derivatives(x[i], σ)
    dj, rj, cj = distance_derivatives(x[j], σ)
    return (wi * di + wj * dj) / (wi + wj),
           (wi * ri + wj * rj) / (wi + wj),
           (wi * ci + wj * cj) / (wi + wj)
end

function virtualasteroids(x::AbstractVector{CloseApproach{T}}) where {T <: Real}
    VAs = Vector{VirtualAsteroid{T}}(undef, 0)
    y = deepcopy(x)
    while !isempty(y)
        CAs = splice!(y, 1:1)
        for _ in eachindex(y)
            j = findfirst(z -> ubound(z) == lbound(CAs[1]) && difft(z, CAs[1]) < 45, y)
            isnothing(j) && break
            z = popat!(y, j)
            pushfirst!(CAs, z)
        end
        for _ in eachindex(y)
            j = findfirst(z -> ubound(CAs[end]) == lbound(z) && difft(CAs[end], z) < 45, y)
            isnothing(j) && break
            z = popat!(y, j)
            push!(CAs, z)
        end
        domain = (minimum(lbound, CAs), maximum(ubound, CAs))
        push!(VAs, VirtualAsteroid{T}(domain, CAs))
    end

    return VAs
end