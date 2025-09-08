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

function convergence_domain(x::CloseApproach{T}, ϵ::T) where {T <: Real}
    d = domain_radius(x)
    r = convergence_radius(x, ϵ)
    return (max(lbound(x), x.σ - r*d), min(ubound(x), x.σ + r*d))
end

function distance(x::CloseApproach{T}, σ::T) where {T <: Real}
    # @assert σ in x "`σ` is outside the domain of `x`"
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

using TaylorIntegration: update_cache!, taylorstep!, set_psol!, findroot!,
                         surfacecrossing

function closeapproaches(f!::D, q0::Vector{Taylor1{T}}, jd0::T, nyears::T,
                         buffer::PropagationBuffer{T, Taylor1{T}, T},
                         params::Parameters{T}, σ::T, domain::NTuple{2, T},
                         ϵ::T; eventorder::Int = 0, newtoniter::Int = 10,
                         nrabstol::T = eps(T)) where {D, T <: Real}
    # Unpack
    @unpack abstol, maxsteps = params
    @unpack cache, dparams = buffer
    # Update reference epoch
    dparams.jd0 = jd0
    # Initial and final time
    t0, tmax = zero(T), nyears * yr

    # @assert order ≥ eventorder "`eventorder` must be less than or equal to `order`"

    @unpack tv, xv, psol, xaux, t, x, dx, rv, parse_eqs = cache

    # Initial conditions
    x0 = deepcopy(q0)
    update_cache!(cache, t0, x0)
    @inbounds tv[1] = t0
    @inbounds xv[:, 1] .= deepcopy(q0)
    sign_tstep = copysign(1, tmax - t0)

    # Some auxiliary arrays for root-finding/event detection/Poincaré surface of section evaluation
    g_tupl = rvelea(dx, x, dparams, t)
    g_tupl_old = deepcopy(g_tupl)
    δt = zero(x[1])
    δt_old = zero(x[1])

    x_dx = vcat(x, dx)
    g_dg = vcat(g_tupl[2], g_tupl_old[2])
    x_dx_val = TS.evaluate(x_dx)
    g_dg_val = vcat(TS.evaluate(g_tupl[2]), TS.evaluate(g_tupl_old[2]))

    tvS = Array{Taylor1{T}}(undef, maxsteps + 1)
    xvS = similar(xv)
    gvS = similar(tvS)

    CAs = Vector{CloseApproach{T}}(undef, 0)

    # Integration
    nsteps = 1
    nevents = 1 #number of detected events
    while sign_tstep * t0 < sign_tstep * tmax
        δt_old = δt
        δt = taylorstep!(Val(parse_eqs), f!, t, x, dx, xaux, abstol, dparams, rv) # δt is positive!
        # Below, δt has the proper sign according to the direction of the integration
        δt = sign_tstep * min(δt, sign_tstep * (tmax - t0))
        evaluate!(x, δt, x0) # new initial condition
        set_psol!(Val(true), psol, nsteps, x) # Store the Taylor polynomial solution
        g_tupl = rvelea(dx, x, dparams, t)
        nevents = findroot!(t, x, dx, g_tupl_old, g_tupl, eventorder, tvS, xvS,
                            gvS, t0, δt_old, x_dx, x_dx_val, g_dg, g_dg_val,
                            nrabstol, newtoniter, nevents)
        if surfacecrossing(g_tupl_old, g_tupl, eventorder)
            # Time at close approach
            t_CA = tvS[nevents-1] + (jd0 - JD_J2000)
            # Asteroid's geocentric state vector
            xae = x(tvS[nevents-1] - t0) - params.eph_ea(t_CA)
            # Asteroid's geocentric semimajor axis
            a = semimajoraxis(xae..., PE.μ[ea], zero(T))
            if a < 0
                # Earth's heliocentric state vector
                xes = params.eph_ea(t_CA) - params.eph_su(t_CA)
                # Öpik's coordinates
                B = bopik(xae, xes)
                X, Y, Z = B.ξ, B.ζ, B.b
                coords = :bopik
            else
                # Modified target plane
                X, Y = mtp(xae)
                Z = 1.0 * one(x)
                coords = :mtp
            end
            CA = CloseApproach{T}(σ, domain, t_CA, X, Y, Z, coords)
            push!(CAs, CA)
            !isconvergent(CA, ϵ) && break
        end
        g_tupl_old = deepcopy(g_tupl)
        t0 += δt
        update_cache!(cache, t0, x0)
        nsteps += 1
        @inbounds tv[nsteps] = t0
        @inbounds xv[:, nsteps] .= deepcopy(x0)
        if nsteps > maxsteps
            @warn("""
            Maximum number of integration steps reached; exiting.
            """)
            break
        end
    end

    return CAs
end