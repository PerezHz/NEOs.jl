"""
    AbstractCloseApproach{T <: Real} <: AbstractImpactMonitoring

Supertype for the close approaches interface.
"""
abstract type AbstractCloseApproach{T <: Real} <: AbstractImpactMonitoring end

"""
    CloseApproach{T} <: AbstractCloseApproach{T}

The projection over the target plane of a segment of the line of variations
(LOV) at close approach with respect to a planet.

# Fields

- `σ::T`: center of the Taylor expansions.
- `domain::NTuple{2, T}`: segment of the LOV.
- `t::Taylor1{T}`: time of close approach [days since J2000 TDB].
- `x/y::Taylor1{T}`: coordinates on the target plane [planet radii].
- `z::Taylor1{T}`: impact parameter [planet radii].
- `tp::Type{<:AbstractTargetPlane}`: type of target plane, either
    [`BPlane`](@ref) or [`MTP`](@ref).
"""
struct CloseApproach{T} <: AbstractCloseApproach{T}
    σ::T
    domain::NTuple{2, T}
    t::Taylor1{T}
    x::Taylor1{T}
    y::Taylor1{T}
    z::Taylor1{T}
    tp::Type{<:AbstractTargetPlane}
end

function show(io::IO, x::CloseApproach)
    tp = x.tp.name.wrapper
    σ = center(x)
    t = days2dtutc(nominaltime(x))
    r = nominalstate(x)

    print(io, tp, " σ: ", σ, " t: ", t, " coords: ", r)
end

in(σ::Real, x::CloseApproach) = x.domain[1] ≤ σ ≤ x.domain[2]

get_order(x::CloseApproach) = get_order(x.t)

center(x::CloseApproach) = x.σ
lbound(x::CloseApproach) = x.domain[1]
ubound(x::CloseApproach) = x.domain[2]

width(x::NTuple{2, T}) where {T <: Real} = x[2] - x[1]
width(x::CloseApproach) = width(x.domain)

nominaltime(x::CloseApproach) = cte(x.t)
nominalstate(x::CloseApproach) = [cte(x.x), cte(x.y), cte(x.z)]

difft(x::CloseApproach, y::CloseApproach) = abs(nominaltime(x) - nominaltime(y))

domain_radius(x::CloseApproach) = max(x.domain[2] - x.σ, x.σ - x.domain[1])

function convergence_radius(x::CloseApproach, ϵ::Real)
    order = get_order(x)
    return min(
        (ϵ / norm(x.x[end], Inf))^(1/order),
        (ϵ / norm(x.y[end], Inf))^(1/order),
        (ϵ / norm(x.z[end], Inf))^(1/order)
    )
end

function convergence_domain(x::CloseApproach, ϵ::Real)
    d = domain_radius(x)
    r = convergence_radius(x, ϵ)
    return (x.σ - r*d, x.σ + r*d)
end

isconvergent(x::CloseApproach, ϵ::Real) = convergence_radius(x, ϵ) > 1

deltasigma(x::CloseApproach, σ::Real) = (σ - center(x)) / domain_radius(x)

function timeofca(x::CloseApproach, σ::Real)
    dσ = deltasigma(x, σ)
    return x.t(dσ)
end

function targetplane(x::CloseApproach, σ::Real)
    dσ = deltasigma(x, σ)
    return [x.x(dσ), x.y(dσ), x.z(dσ)]
end

function distance(x::CloseApproach, σ::Real)
    dσ = deltasigma(x, σ)
    return hypot(x.x(dσ), x.y(dσ)) - x.z(dσ)
end

# TO DO: Use Horner's rule to evaluate derivatives
function rvelea(x::CloseApproach, σ::Real)
    dσ = deltasigma(x, σ)
    ξ, ζ = x.x(dσ), x.y(dσ)
    dξ, dζ = TS.differentiate(x.x)(dσ), TS.differentiate(x.y)(dσ)
    return (ξ*dξ + ζ*dζ) / domain_radius(x)
end

function concavity(x::CloseApproach, σ::Real)
    dσ = deltasigma(x, σ)
    ξ, ζ = x.x(dσ), x.y(dσ)
    dξ, dζ = TS.differentiate(x.x)(dσ), TS.differentiate(x.y)(dσ)
    d2ξ, d2ζ = TS.differentiate(x.x, 2)(dσ), TS.differentiate(x.y, 2)(dσ)
    return (ξ*d2ξ + dξ^2 + ζ*d2ζ + dζ^2) / domain_radius(x)^2
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
        # B-Plane in Öpik's coordinates
        tp = bopik(xae, xes)
    else
        # Modified target plane
        tp = mtp(xae)
    end
    return CloseApproach{T}(σ, domain, t, targetplane(TP)..., typeof(TP))
end

"""
    closeapproaches(lov, σ, domain, nyears, params; kwargs...)

Return the vector of close approaches obtained by integrating the expansion of `lov`
at `σ`, with scaling factors such that `domain` is mapped into the interval `[-1, 1]`,
for a period of `nyears` [years]. For a list of parameters see the `Propagation`
section of [`Parameters`](@ref).

# Keyword arguments

- `lovorder::Int`: order of the expansion of `lov` (default: `min(6, get_order(lov))`).
- `ctol::Real`: convergence tolerance (default: `0.01`); the integration will stop at
    the first close approach not converging under `ctol`.
- `newtoniter::Int`: maximum Newton-Raphson iterations per detected event (default: `10`).
- `nrabstol::Real`: allowed tolerance for the Newton-Raphson process (default: `eps()`).
- `buffer::Union{Nothing, PropagationBuffer}`: pre-allocated memory (default: `nothing`).
"""
function closeapproaches(
        lov::LOV{D, T}, σ::T, domain::NTuple{2, T}, nyears::T,
        params::Parameters{T}; lovorder::Int = min(6, get_order(lov)),
        ctol::T = 0.01, newtoniter::Int = 10, nrabstol::T = eps(T),
        buffer::Union{Nothing, PropagationBuffer{T, Taylor1{T}, T}} = nothing
    ) where {D, T <: Real}
    # Unpack parameters
    @unpack abstol, maxsteps = params
    # Dynamical model
    dynamics = lov.dynamics
    # Reference epoch
    jd0 = epoch(lov)
    d0 = jd0 - JD_J2000
    # Initial conditions
    t0, tmax = zero(T), nyears * yr
    q0 = lov(σ + max(domain[2] - σ, σ - domain[1]) * Taylor1(lovorder))
    # Propagation buffer
    if isnothing(buffer)
        tlim = (d0, d0 + tmax)
        buffer = PropagationBuffer(dynamics, q0, jd0, tlim, params)
    end
    @unpack cache, dparams = buffer
    @unpack tv, xv, psol, xaux, t, x, dx, rv, parse_eqs = cache
    dparams.jd0 = epoch(lov)

    # Specialized version of the root-finding method of taylorinteg
    x0 = deepcopy(q0)
    update_cache!(cache, t0, x0)
    @inbounds tv[1] = t0
    @inbounds xv[:, 1] .= deepcopy(q0)
    sign_tstep = copysign(1, tmax - t0)

    # Some auxiliary arrays for root-finding/event detection/Poincaré
    # surface of section evaluation
    g_tupl = (false, rvelea(dx, x, dparams, t)[2])
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
    nevents = 1 # number of detected events
    while sign_tstep * t0 < sign_tstep * tmax
        δt_old = δt
        δt = taylorstep!(Val(parse_eqs), dynamics, t, x, dx, xaux, abstol, dparams, rv) # δt is positive!
        # Below, δt has the proper sign according to the direction of the integration
        δt = sign_tstep * min(δt, sign_tstep * (tmax - t0))
        TS.evaluate!(x, δt, x0) # new initial condition
        set_psol!(Val(true), psol, nsteps, x) # Store the Taylor polynomial solution
        g_tupl = rvelea(dx, x, dparams, t)
        nevents_old = nevents
        nevents = findroot!(t, x, dx, g_tupl_old, g_tupl, 0, tvS, xvS, gvS,
                            t0, δt_old, x_dx, x_dx_val, g_dg, g_dg_val, nrabstol,
                            newtoniter, nevents)
        if nevents > nevents_old
            # Time at close approach
            t_CA = d0 + tvS[nevents-1]
            # Asteroid's geocentric state vector
            xae = psol[:, nsteps-1]( t_CA - d0 - tv[nsteps-1] ) - params.eph_ea(t_CA)
            # Asteroid's geocentric semimajor axis
            a = semimajoraxis(xae..., PE.μ[ea], zero(T))
            if a < 0
                # Earth's heliocentric state vector
                xes = params.eph_ea(t_CA) - params.eph_su(t_CA)
                # BPlane in Öpik's coordinates
                tp = bopik(xae, xes)
            else
                # Modified target plane
                tp = mtp(xae)
            end
            # Close approach
            CA = CloseApproach{T}(σ, domain, t_CA, targetplane(tp)..., typeof(tp))
            push!(CAs, CA)
            !isconvergent(CA, ctol) && break
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