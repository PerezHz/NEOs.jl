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

function radialvelocity(x::CloseApproach, σ::Real)
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

function CloseApproach(IM::AbstractIMProblem{D, T}, σ::T, domain::NTuple{2, T},
                       t::Taylor1{T}, eph::DensePropagation2{T, Taylor1{T}},
                       params::Parameters{T}) where {D, T <: Real}
    # Unpack
    @unpack target = IM
    @unpack eph_su = params
    # Asteroid's planetocentric state vector
    xae = eph(t) - target(t)
    # Asteroid's planetocentric semimajor axis
    a = semimajoraxis(xae..., gm(target), zero(T))
    if a < 0
        # Planet's heliocentric state vector
        xes = target(t) - eph_su(t)
        # B-Plane in Öpik's coordinates
        tp = bopik(xae, xes, target)
    else
        # Modified target plane
        tp = mtp(xae, target)
    end
    return CloseApproach{T}(σ, domain, t, targetplane(tp)..., typeof(tp))
end

# Asteroid's radial velocity with respect to a planet
function radialvelocity(x, target, params, t)
    # Julian date (TDB) of start time
    jd0 = params.jd0
    # Days since J2000.0 = 2.451545e6
    dsj2k = t + (jd0 - JD_J2000)
    # Planet ephemeris at dsj2k
    xe = target(dsj2k)
    # Asteroid's planetocentric state vector
    xae = x[1:6] - xe
    # Planetocentric radial velocity
    return euclid3D(cte(xae[1:3])) < 0.2, dot3D(xae[1:3], xae[4:6])
end

"""
    closeapproaches(IM, lov, σ, domain, nyears, params; kwargs...)

Return the vector of close approaches under the impact monitoring problem `IM`,
corresponding to the expansion of `lov` at `σ`, with scaling factors such that
`domain` is mapped into the `[-1, 1]` interval, propagated for a period of
`nyears` [years]. For a list of parameters see the `Propagation` section of
 [`Parameters`](@ref).

The propagation will stop before the final time if it finds a close approach that
does not converge under the tolerance `ctol`.

# Keyword arguments

- `lovorder::Int`: order of the expansion of `lov` (default: `min(6, get_order(lov))`).
- `ctol::Real`: convergence tolerance (default: `0.01`); the integration will stop at
    the first close approach not converging under `ctol`.
- `newtoniter::Int`: maximum Newton-Raphson iterations per detected event (default: `10`).
- `nrabstol::Real`: allowed tolerance for the Newton-Raphson process (default: `eps()`).
- `buffer::Union{Nothing, PropagationBuffer}`: pre-allocated memory (default: `nothing`).
"""
function closeapproaches(
        IM::AbstractIMProblem{D, T}, lov::LineOfVariations{D, T}, σ::T,
        domain::NTuple{2, T}, nyears::T, params::Parameters{T};
        lovorder::Int = min(6, get_order(lov)), ctol::T = 0.01,
        newtoniter::Int = 10, nrabstol::T = eps(T),
        buffer::Union{Nothing, PropagationBuffer{T, Taylor1{T}, T}} = nothing
    ) where {D, T <: Real}
    # Unpack
    @unpack target = IM
    @unpack abstol, maxsteps, order, eph_su = params
    # Dynamical model
    dynamics = dynamicalmodel(lov)
    # Reference epoch
    d0 = epoch(lov)
    jd0 = d0 + JD_J2000
    # Initial conditions
    tpre, t0, tmax = zero(T), zero(T), nyears * yr
    q0 = lov(σ, domain, lovorder)
    # Propagation buffer
    tlim = (d0, d0 + tmax)
    if isnothing(buffer)
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
    pbuffer = EphemerisEvaluationBuffer(target.eph, tlim, order, q0)
    g_tupl = (false, radialvelocity(x, pbuffer, dparams, t)[2])
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
        g_tupl = radialvelocity(x, pbuffer, dparams, t)
        nevents_old = nevents
        nevents = findroot!(t, x, dx, g_tupl_old, g_tupl, 0, tvS, xvS, gvS,
                            t0, δt_old, x_dx, x_dx_val, g_dg, g_dg_val, nrabstol,
                            newtoniter, nevents)
        if nevents > nevents_old
            # Time at close approach
            t_CA = d0 + tvS[nevents-1]
            # Asteroid's planetocentric state vector
            xae = psol[1:6, 1]( t_CA - d0 - tpre ) - target(t_CA)
            # Asteroid's planetocentric semimajor axis
            a = semimajoraxis(xae..., gm(target), zero(T))
            if a < 0
                # Planet's heliocentric state vector
                xes = target(t_CA) - eph_su(t_CA)
                # BPlane in Öpik's coordinates
                tp = bopik(xae, xes, target)
            else
                # Modified target plane
                tp = mtp(xae, target)
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

    return CAs
end