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
    tp = typeof(x.tp)
    σ = center(x)
    t = days2dtutc(nominaltime(x))
    r = nominalstate(x)

    print(io, tp, " σ: ", σ, " t: ", t, " coords: ", r)
end

in(σ::Real, x::CloseApproach) = x.domain[1] ≤ σ ≤ x.domain[2]

center(x::CloseApproach) = x.σ
lbound(x::CloseApproach) = x.domain[1]
ubound(x::CloseApproach) = x.domain[2]
width(x::CloseApproach) = x.domain[2] - x.domain[1]

nominaltime(x::CloseApproach) = cte(x.t)
nominalstate(x::CloseApproach) = [cte(x.x), cte(x.y), cte(x.z)]

difft(a::CloseApproach, b::CloseApproach) = abs(nominaltime(a) - nominaltime(b))

domain_radius(x::CloseApproach) = max(x.domain[2] - x.σ, x.σ - x.domain[1])

function convergence_radius(x::CloseApproach, ϵ::Real)
    order = get_order(x.x)
    return min(
        (ϵ / norm(x.x[end], Inf))^(1/order),
        (ϵ / norm(x.y[end], Inf))^(1/order),
        (ϵ / norm(x.z[end], Inf))^(1/order)
    )
end

isconvergent(x::CloseApproach, ϵ::Real) = convergence_radius(x, ϵ) > 1

function convergence_domain(x::CloseApproach, ϵ::Real)
    d = domain_radius(x)
    r = convergence_radius(x, ϵ)
    return (x.σ - r*d, x.σ + r*d)
end

deltasigma(x::CloseApproach, σ::Real) = (σ - center(x)) / domain_radius(x)

function time(x::CloseApproach, σ::Real)
    dσ = deltasigma(x, σ)
    return x.t(dσ)
end

function targetplane(x::CloseApproach, σ::Real)
    dσ = deltasigma(x, σ)
    return [x.x(dσ), x.y(dσ), x.z(dσ)]
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
    closeapproaches(f!, q0, jd0, nyears, buffer, params, σ, domain, ϵ; kwargs...)

Return the vector of close approaches obtained by integrating the dynamical model `f!`
starting from an initial condition `q0` [au, au/day] at time `jd0` [julian date TDB]
for a period of `nyears` [years]. For a list of parameters see the `Propagation`
section of [`Parameters`](@ref).

The initial condition `q0` is assumed to come from expanding the line of variations
(LOV) at `σ` with scaling factors such that `domain` is transformed into the interval
[-1, 1]. The integration stops when the last close approach does not converge under the
tolerance `ϵ`.

# Keyword arguments

- `eventorder::Int`: order of the derivative of `NEOs.rvelea` whose roots are computed
    (default: `0`).
- `newtoniter::Int`: maximum Newton-Raphson iterations per detected root (default: `10`).
- `nrabstol::T`: allowed tolerance for the Newton-Raphson process (default: `eps(T)`).
"""
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
    # Find close approaches
    CAs = lov_taylorinteg!(Val(true), f!, rvelea, q0, zero(T), nyears * yr, abstol,
                           cache, dparams, params, σ, domain, ϵ; maxsteps, eventorder,
                           newtoniter, nrabstol)

    return CAs
end

# Specialized version of the root-finding method of taylorinteg
function lov_taylorinteg!(dense::Val{D}, f!, g, q0::Array{U,1}, t0::T, tmax::T,
                          abstol::T, cache::VectorCache, dparams, params, σ::T,
                          domain::NTuple{2, T}, ϵ::T; maxsteps::Int = 500,
                          eventorder::Int = 0, newtoniter::Int = 10,
                          nrabstol::T = eps(T)) where {D, T <: Real, U <: Number}

    @unpack tv, xv, psol, xaux, t, x, dx, rv, parse_eqs = cache
    @unpack jd0 = dparams

    # Initial conditions
    epoch0 = jd0 - JD_J2000
    x0 = deepcopy(q0)
    update_cache!(cache, t0, x0)
    @inbounds tv[1] = t0
    @inbounds xv[:, 1] .= deepcopy(q0)
    sign_tstep = copysign(1, tmax - t0)

    # Some auxiliary arrays for root-finding/event detection/Poincaré
    # surface of section evaluation
    g_tupl = (false, g(dx, x, dparams, t)[2])
    g_tupl_old = deepcopy(g_tupl)
    δt = zero(x[1])
    δt_old = zero(x[1])

    x_dx = vcat(x, dx)
    g_dg = vcat(g_tupl[2], g_tupl_old[2])
    x_dx_val = TS.evaluate(x_dx)
    g_dg_val = vcat(TS.evaluate(g_tupl[2]), TS.evaluate(g_tupl_old[2]))

    tvS = Array{U}(undef, maxsteps + 1)
    xvS = similar(xv)
    gvS = similar(tvS)

    CAs = Vector{CloseApproach{T}}(undef, 0)

    # Integration
    nsteps = 1
    nevents = 1 # number of detected events
    while sign_tstep * t0 < sign_tstep * tmax
        δt_old = δt
        δt = taylorstep!(Val(parse_eqs), f!, t, x, dx, xaux, abstol, dparams, rv) # δt is positive!
        # Below, δt has the proper sign according to the direction of the integration
        δt = sign_tstep * min(δt, sign_tstep * (tmax - t0))
        TS.evaluate!(x, δt, x0) # new initial condition
        set_psol!(dense, psol, nsteps, x) # Store the Taylor polynomial solution
        g_tupl = g(dx, x, dparams, t)
        nevents_old = nevents
        nevents = findroot!(t, x, dx, g_tupl_old, g_tupl, eventorder, tvS, xvS, gvS,
                            t0, δt_old, x_dx, x_dx_val, g_dg, g_dg_val, nrabstol,
                            newtoniter, nevents)
        if nevents > nevents_old
            # Time at close approach
            t_CA = epoch0 + tvS[nevents-1]
            # Asteroid's geocentric state vector
            xae = psol[:, nsteps-1]( t_CA - epoch0 - tv[nsteps-1] ) - params.eph_ea(t_CA)
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