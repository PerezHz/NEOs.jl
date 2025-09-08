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

lbound(x::VirtualAsteroid) = x.domain[1]
ubound(x::VirtualAsteroid) = x.domain[2]

width(x::VirtualAsteroid) = x.domain[2] - x.domain[1]

isconvergent(x::VirtualAsteroid, ϵ::Real) = all(Base.Fix2(isconvergent, ϵ), x.CAs)

convergence_domain(x::VirtualAsteroid, ϵ::Real) =
    (convergence_domain(x.CAs[1], ϵ)[1], convergence_domain(x.CAs[end], ϵ)[2])

closeapproaches(x::VirtualAsteroid) = x.CAs

meantime(x::VirtualAsteroid) = mean(nominaltime, x.CAs)

function distance(x::VirtualAsteroid, σ::Real, ϵ::Real)
    d = convergence_domain(x, ϵ)
    @assert d[1] ≤ σ ≤ d[2] "`σ` is outside the convergence domain domain of `x`"
    i = findfirst(Base.Fix1(in, σ), x.CAs)
    # Normal evaluation of the Taylor polynomials
    if isconvergent(x.CAs[i], ϵ) ||
        (i == 1 && d[1] ≤ σ ≤ x.CAs[1].σ) ||
        (i == length(x.CAs) && x.CAs[end].σ ≤ σ ≤ d[2])
        return distance(x.CAs[i], σ)
    end
    # Exponentially decaying average
    i = findlast( y -> y.σ < σ, x.CAs)
    j = findfirst(y -> σ < y.σ, x.CAs)
    CAi, CAj = x.CAs[i], x.CAs[j]
    dσi = abs(σ - CAi.σ) / (convergence_radius(CAi, ϵ) * domain_radius(CAi))
    dσj = abs(σ - CAj.σ) / (convergence_radius(CAj, ϵ) * domain_radius(CAj))
    wi = 1 / 10^(dσi - 1)
    wj = 1 / 10^(dσj - 1)
    di = distance(CAi, σ)
    dj = distance(CAj, σ)
    return (wi * di + wj * dj) / (wi + wj)
end

function mindistance(x::VirtualAsteroid{T}, tol::T = width(x) / 1E3) where {T <: Real}
    # 1 / φ
    invphi = (sqrt(5) - 1) / 2
    # 1 / φ^2
    invphi2 = (3 - sqrt(5)) / 2
    # Interval bounds
    a, b = x.domain
    # Interval width
    h = b - a
    # Termination condition
    h <= tol && return distance(x, (a + b) / 2)
    # Required steps to achieve tolerance
    n = ceil(Int, log(tol/h) / log(invphi))
    # Initialize center points
    c = a + invphi2 * h
    d = a + invphi * h
    yc = distance(x, c)
    yd = distance(x, d)
    # Main loop
    for _ in 1:n
        if yc < yd
            b = d
            d = c
            yd = yc
            h = invphi * h
            c = a + invphi2 * h
            yc = distance(x, c)
        else
            a = c
            c = d
            yc = yd
            h = invphi * h
            d = a + invphi * h
            yd = distance(x, d)
        end
    end

    if yc < yd
        return distance(x, (a + d) / 2)
    else
        return distance(x, (c + b) / 2)
    end
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

impact_probability(a, b) = erf(a/sqrt(2), b/sqrt(2)) / 2

function impact_probability(x::VirtualAsteroid{T}, ϵ::T) where {T <: Real}
    a, b = convergence_domain(x, ϵ)
    rs = find_zeros(σ -> distance(x, σ, ϵ), a, b, no_pts = 100)
    if isempty(rs)
        return zero(T)
    elseif length(rs) == 2
        return impact_probability(rs[1], rs[2])
    elseif length(rs) == 4
        return impact_probability(rs[1], rs[2]) + impact_probability(rs[3], rs[4])
    else
        throw(ArgumentError("Cannot process $(length(rs)) roots"))
    end
end