"""
    CloseApproach{T} <: AbstractLineOfVariations{T}

The projection over the target plane of a segment of the line of variations
at close approach with respect to a planet.

# Fields

- `σ::T`: center of the Taylor expansions.
- `domain::NTuple{2, T}`: segment of the line of variations.
- `t::Taylor1{T}`: time of close approach [days since J2000 TDB].
- `x/y::Taylor1{T}`: coordinates on the target plane [planet radii].
- `z::Taylor1{T}`: impact parameter [planet radii].
- `tp::Type{<:AbstractTargetPlane}`: type of target plane, either
    [`BPlane`](@ref) or [`MTP`](@ref).
"""
@auto_hash_equals struct CloseApproach{T} <: AbstractLineOfVariations{T}
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
    σ = sigma(x)
    t = date(x)
    r = nominalstate(x)

    print(io, tp, " σ: ", σ, " t: ", t, " coords: ", r)
end

nominaltime(x::CloseApproach) = cte(x.t)

sigma(x::CloseApproach) = x.σ

get_order(x::CloseApproach) = get_order(x.t)

nominalstate(x::CloseApproach) = [cte(x.x), cte(x.y), cte(x.z)]

difft(x::CloseApproach, y::CloseApproach) = abs(nominaltime(x) - nominaltime(y))

domain_radius(x::CloseApproach) = max(ubound(x) - sigma(x), sigma(x) - lbound(x))

convergence_radius(x::Taylor1, ctol::Real) = (ctol / norm(x[end], Inf))^(1 / get_order(x))

convergence_radius(x::AbstractVector, ctol::Real) = minimum(convergence_radius.(x, ctol))

function convergence_radius(x::CloseApproach, ctol::Real)
    return min(
        convergence_radius(x.x, ctol),
        convergence_radius(x.y, ctol),
        convergence_radius(x.z, ctol),
    )
end

function convergence_domain(x::CloseApproach, ctol::Real)
    d = domain_radius(x)
    r = convergence_radius(x, ctol)
    return (sigma(x) - r*d, sigma(x) + r*d)
end

isconvergent(x, ctol::Real) = convergence_radius(x, ctol) > 1

function refine(σ::T, domain::NTuple{2, T}, N::Int,
                dist::Symbol = :uniform) where {T <: Real}
    @assert N ≥ 1 && isodd(N) "N must be at least one and odd to have an \
        expansion around the domain center"
    if dist == :uniform
        d = Uniform{T}(domain[1], domain[2])
    elseif dist == :normal
        d = Normal{T}(zero(T), one(T))
    else
        throw(ArgumentError("dist must be either :uniform or :normal"))
    end
    xs = LinRange(cdf(d, domain[1]), cdf(d, domain[2]), N+1)
    endpoints = quantile.(d, xs)
    endpoints[1], endpoints[end] = domain[1], domain[2]
    domains = [(endpoints[i], endpoints[i+1]) for i in 1:N]
    σs = [(domains[i][1] + domains[i][2]) / 2 for i in eachindex(domains)]
    if iszero(σ) || dist == :uniform
        σs[(N ÷ 2) + 1] = σ
    end
    return σs, domains
end

deltasigma(x::CloseApproach, σ::Real) = (σ - sigma(x)) / domain_radius(x)

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

function rvelea(x::CloseApproach, σ::Real)
    dσ = deltasigma(x, σ)
    n = get_order(x)
    ξ, ζ, b = x.x[end], x.y[end], x.z[end]
    dξ, dζ, db = zero(ξ), zero(ζ), zero(b)
    for i in n-1:-1:0
        dξ = ξ + dσ * dξ
        dζ = ζ + dσ * dζ
        db = b + dσ * db
        ξ = x.x[i] + dσ * ξ
        ζ = x.y[i] + dσ * ζ
        b = x.z[i] + dσ * b
    end
    r = hypot(ξ, ζ)
    rv = ξ*dξ + ζ*dζ
    return (rv/r - db) / domain_radius(x)
end

function concavity(x::CloseApproach, σ::Real)
    dσ = deltasigma(x, σ)
    n = get_order(x)
    ξ = x.x[end-1] + dσ * x.x[end]
    ζ = x.y[end-1] + dσ * x.y[end]
    b = x.z[end-1] + dσ * x.z[end]
    dξ, dζ, db = x.x[end], x.y[end], x.z[end]
    d2ξ, d2ζ, d2b = zero(ξ), zero(ζ), zero(b)
    for i in n-2:-1:0
        d2ξ = 2dξ + dσ * d2ξ
        d2ζ = 2dζ + dσ * d2ζ
        d2b = 2db + dσ * d2b
        dξ = ξ + dσ * dξ
        dζ = ζ + dσ * dζ
        db = b + dσ * db
        ξ = x.x[i] + dσ * ξ
        ζ = x.y[i] + dσ * ζ
        b = x.z[i] + dσ * b
    end
    r = hypot(ξ, ζ)
    rv = ξ*dξ + ζ*dζ
    ra = ξ*d2ξ + dξ^2 + ζ*d2ζ + dζ^2
    return ( ra / r - rv^2 / r^3 - d2b) / domain_radius(x)^2
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
    return CloseApproach{T}(σ, domain, t, targetplane(tp)..., typeof(tp))
end

"""
    closeapproaches(lov, σ, domain, nyears, params; kwargs...)

Compute the close approaches of the initial condition obtained by expanding `lov` at
`σ`, with scaling factors such that `domain` is mapped into the interval `[-1, 1]`,
for a period of `nyears` [years]. For a list of parameters see the `Propagation`
section of [`Parameters`](@ref).

The propagation will stop before the final time if a close approach does not converge
under the tolerance `ctol`. Thus, the function returns (i) the vector of detected
close approaches, (ii) the convergence radius at the last timestep and (iii) a
bool indicating whether the integration reached the final time.

# Keyword arguments

- `lovorder::Int`: order of the expansion of `lov` (default: `min(6, get_order(lov))`).
- `ctol::Real`: convergence tolerance (default: `0.01`); the integration will stop at
    the first close approach not converging under `ctol`.
- `newtoniter::Int`: maximum Newton-Raphson iterations per detected event (default: `10`).
- `nrabstol::Real`: allowed tolerance for the Newton-Raphson process (default: `eps()`).
- `buffer::Union{Nothing, PropagationBuffer}`: pre-allocated memory (default: `nothing`).
"""
function closeapproaches(
        lov::LineOfVariations{D, T}, σ::T, domain::NTuple{2, T}, nyears::T,
        params::Parameters{T}; lovorder::Int = min(6, get_order(lov)),
        ctol::T = 0.01, newtoniter::Int = 10, nrabstol::T = eps(T),
        buffer::Union{Nothing, PropagationBuffer{T, Taylor1{T}, T}} = nothing
    ) where {D, T <: Real}
    # Unpack parameters
    @unpack abstol, maxsteps, order = params
    # Dynamical model
    dynamics = lov.dynamics
    # Reference epoch
    d0 = epoch(lov)
    jd0 = d0 + JD_J2000
    # Initial conditions
    tpre, t0, tmax = zero(T), zero(T), nyears * yr
    q0 = lov(σ + max(domain[2] - σ, σ - domain[1]) * Taylor1(lovorder))
    # Propagation buffer
    if isnothing(buffer)
        tlim = (d0, d0 + tmax)
        buffer = PropagationBuffer(dynamics, q0, jd0, tlim, params)
    end
    @unpack cache, dparams = buffer
    @unpack psol, xaux, t, x, dx, rv, parse_eqs = cache
    dparams.jd0 = jd0

    # Specialized version of the root-finding method of taylorinteg
    x0 = deepcopy(q0)
    update_cache!(cache, t0, x0)
    psol[:, 1] .= zero.(x)
    psol[:, 2] .= zero.(x)
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
    xvS = Array{Taylor1{T}}(undef, length(q0), maxsteps + 1)
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
        # New initial condition
        TS.evaluate!(x, δt, x0)
        # Store the Taylor polynomial solution
        for i in eachindex(x)
            for k in eachindex(x[i])
                TS.identity!(psol[:, 1][i], psol[:, 2][i], k)
                TS.identity!(psol[:, 2][i], x[i], k)
            end
        end
        g_tupl = rvelea(dx, x, dparams, t)
        nevents_old = nevents
        nevents = findroot!(t, x, dx, g_tupl_old, g_tupl, 0, tvS, xvS, gvS,
                            t0, δt_old, x_dx, x_dx_val, g_dg, g_dg_val, nrabstol,
                            newtoniter, nevents)
        if nevents > nevents_old
            # Time at close approach
            t_CA = d0 + tvS[nevents-1]
            # Asteroid's geocentric state vector
            xae = psol[:, 1]( t_CA - d0 - tpre ) - params.eph_ea(t_CA)
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
        # Update time
        tpre = t0
        t0 += δt
        update_cache!(cache, t0, x0)
        nsteps += 1
        if nsteps > maxsteps
            @warn("""
            Maximum number of integration steps reached; exiting.
            """)
            break
        end
    end
    flag = sign_tstep * t0 >= sign_tstep * tmax
    r = flag ? convergence_radius(x0, ctol) : convergence_radius(CAs[end], ctol)

    return CAs, r, flag
end