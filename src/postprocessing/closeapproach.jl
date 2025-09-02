"""
    CloseApproach{T <: Real}

A segment of the line of variations (LOV) at close approach over the
target plane.

# Fields

- `σ::T`: center of the Taylor expansions.
- `domain::NTuple{2, T}`: segment of the LOV.
- `t::Taylor1{T}`: time of close approach [days since J2000 TDB].
- `x/y::Taylor1{T}`: coordinates on the target plane [Earth radii].
- `z::Taylor1{T}`: impact parameter [Earth radii].
- `coords::Symbol`: coordinate system, either `:bopik` for the Öpik
    coordinates on the B-plane of `:mtp` for the modified target plane.
"""
struct CloseApproach{T <: Real}
    σ::T
    domain::NTuple{2, T}
    t::Taylor1{T}
    x::Taylor1{T}
    y::Taylor1{T}
    z::Taylor1{T}
    coords::Symbol
end

function show(io::IO, x::CloseApproach)
    coords = x.coords == :bopik ? "Öpik" : "MTP"
    σ = x.σ
    t = days2dtutc(nominaltime(x))
    r = nominalstate(x)

    print(io, coords, " σ: ", σ, " t: ", t, " coords: ", r)
end

in(σ::Real, x::CloseApproach) = x.domain[1] ≤ σ ≤ x.domain[2]

lbound(x::CloseApproach) = x.domain[1]
ubound(x::CloseApproach) = x.domain[2]
width(x::CloseApproach) = x.domain[2] - x.domain[1]

nominaltime(x::CloseApproach) = cte(x.t)
nominalstate(x::CloseApproach) = [cte(x.x), cte(x.y), cte(x.z)]

difft(a::CloseApproach, b::CloseApproach) = abs(nominaltime(a) - nominaltime(b))

domain_radius(x::CloseApproach) = max(x.domain[2] - x.σ, x.σ - x.domain[1])

function convergence_radius(x::CloseApproach{T}, ϵ::T) where {T <: Real}
    order = get_order(x.x)
    return min(
        (ϵ / norm(x.x[end], Inf))^(1/order),
        (ϵ / norm(x.y[end], Inf))^(1/order),
        (ϵ / norm(x.z[end], Inf))^(1/order)
    )
end

isconvergent(x::CloseApproach, ϵ::Real) = convergence_radius(x, ϵ) > 1

function distance(x::CloseApproach{T}, σ::T) where {T <: Real}
    @assert σ in x "`σ` is outside the domain of `x`"
    dσ = (σ - x.σ) / domain_radius(x)
    return hypot(x.x(dσ), x.y(dσ)) - x.z(dσ)
end

function mindistance(x::CloseApproach{T}, tol::T = width(x) / 1E3) where {T <: Real}
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

function CloseApproach(σ::T, domain::NTuple{2, T}, t::Taylor1{T},
                       eph::DensePropagation2{T, Taylor1{T}},
                       params::Parameters{T}) where {T <: Real}
    # Asteroid's geocentric state vector
    xae = eph(t) - params.eph_ea(t)
    # Asteroid's geocentric semimajor axis
    a = semimajoraxis(xae..., PE.μ[ea], zero(T))
    if a < 0
        # Earth's heliocentric state vector
        xes = params.eph_ea(t) - params.eph_su(t)
        # Öpik's coordinates
        B = bopik(xae, xes)
        x, y, z = B.ξ, B.ζ, B.b
        coords = :bopik
    else
        # Modified target plane
        x, y = mtp(xae)
        z = 1.0 * one(x)
        coords = :mtp
    end
    return CloseApproach{T}(σ, domain, t, x, y, z, coords)
end