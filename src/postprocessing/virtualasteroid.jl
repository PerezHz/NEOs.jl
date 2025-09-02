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

closeapproaches(x::VirtualAsteroid) = x.CAs

meantime(x::VirtualAsteroid) = mean(nominaltime, x.CAs)

function distance(x::VirtualAsteroid, σ::Real)
    i = findfirst(Base.Fix1(in, σ), x.CAs)
    return distance(x.CAs[i], σ)
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

function impact_domain(x::VirtualAsteroid{T}) where {T <: Real}
    a, b = x.domain
    rs = find_zeros(Base.Fix1(distance, x), a, b, no_pts = 1_000)
    if isempty(rs)
        return distance(x, (a+b)/2) <= 0 ? (a, b) : (b, a)
    elseif length(rs) == 1
        return (distance(x, a) < 0 && distance(x, b) > 0) ? (a, rs[1]) : (rs[1], b)
    elseif length(rs) == 2
        return (rs[1], rs[2])
    else
        throw(ArgumentError("Cannot process more than two roots"))
    end
end

impact_probability(a, b) = erf(a/sqrt(2), b/sqrt(2)) / 2