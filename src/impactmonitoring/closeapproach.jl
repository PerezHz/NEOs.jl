"""
    CloseApproach{T, U <: Number} <: AbstractLineOfVariations{T}

The projection over the target plane of a virtual asteroid at
close approach with respect to a planet.

# Fields

- `σ::T`: LOV index.
- `domain::NTuple{2, T}`: segment of the line of variations.
- `t::U`: time of close approach [days since J2000 TDB].
- `x/y::U`: coordinates on the target plane [planet radii].
- `z::U`: impact parameter [planet radii].
- `tp::Type{<:AbstractTargetPlane}`: type of target plane, either
    [`BPlane`](@ref) or [`MTP`](@ref).
"""
@auto_hash_equals struct CloseApproach{T, U <: Number} <: AbstractLineOfVariations{T}
    σ::T
    domain::NTuple{2, T}
    t::U
    x::U
    y::U
    z::U
    tp::Type{<:AbstractTargetPlane}
end

# Abbreviations
const CloseApproachT1{T} = CloseApproach{T, Taylor1{T}}

# Outer constructors
function CloseApproach(σ::T, domain::NTuple{2, T}, t::U, x::U, y::U, z::U,
                       tp::Type{<:AbstractTargetPlane}) where {T, U}
    return CloseApproach{T, U}(σ, domain, t, x, y, z, tp)
end

CloseApproach(VA::VirtualAsteroid{T, U}, t::U, x::U, y::U, z::U, tp) where {T, U} =
    CloseApproach(sigma(VA), VA.domain, t, x, y, z, tp)

# Print method for CloseApproach
function show(io::IO, x::CloseApproach)
    tp = targetplane(x)
    σ = sigma(x)
    t = date(x)
    r = nominalstate(x)

    print(io, tp, " σ: ", σ, " t: ", t, " coords: ", r)
end

# AbstractLineOfVariations interface
sigma(x::CloseApproach) = x.σ
targetplane(x::CloseApproach) = x.tp
nominaltime(x::CloseApproach) = cte(x.t)
nominalstate(x::CloseApproach) = [cte(x.x), cte(x.y), cte(x.z)]

difft(x::CloseApproach, y::CloseApproach) = abs(nominaltime(x) - nominaltime(y))

get_order(x::CloseApproach{T, U}) where {T, U <: AbstractSeries} = get_order(x.t)

domain_radius(x::CloseApproach) = max(ubound(x) - sigma(x), sigma(x) - lbound(x))

convergence_radius(x::NumberNotSeries, ctol::Real) = ctol / zero(x)

convergence_radius(x::AbstractSeries, ctol::Real) = (ctol / norm(x[end], Inf))^(1 / get_order(x))

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

function timeofca(x::CloseApproachT1, σ::Real)
    dσ = deltasigma(x, σ)
    return x.t(dσ)
end

function targetplane(x::CloseApproachT1, σ::Real)
    dσ = deltasigma(x, σ)
    return [x.x(dσ), x.y(dσ), x.z(dσ)]
end

function distance(x::CloseApproachT1, σ::Real)
    dσ = deltasigma(x, σ)
    return hypot(x.x(dσ), x.y(dσ)) - x.z(dσ)
end

function radialvelocity(x::CloseApproachT1, σ::Real)
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

function concavity(x::CloseApproachT1, σ::Real)
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

# Root-finding function for `closeapproaches`
function closeapproach(x, teph, params, t; R_TP, R_P)
    # Julian date (TDB) of start time
    jd0 = params.jd0
    # Days since J2000.0 = 2.451545e6
    dsj2k = t + (jd0 - JD_J2000)
    # Planet ephemeris at dsj2k
    xp = teph(dsj2k)
    # Asteroid's planetocentric position and velocity
    rap = [x[i] - xp[i] for i in 1:3]
    vap = [x[i] - xp[i] for i in 4:6]
    # Planetocentric distance and radial velocity
    dap = euclid3D(rap)
    rvap = dot3D(rap, vap)
    # Root-finding tuples
    return (true, dap - R_P), (dap < R_TP, rvap)
end

"""
    closeapproaches(IM, VA, nyears, params; kwargs...)

Return the vector of close approaches, under the impact monitoring problem
`IM`, of a virtual asteroid `VA` for a period of `nyears` [years]. For a
list of parameters see the `Propagation` section of  [`Parameters`](@ref).

# Keyword arguments

- `R_TP::Real`: target plane radius in au (default: `0.2`).
- `ctol::Real`: convergence tolerance (default: `Inf`); the integration
    will stop before the final time if it finds a close approach `x`
    such that `!isconvergent(x, ctol)`.
- `newtoniter::Int`: maximum Newton-Raphson iterations per detected event
    (default: `25`).
- `nrabstol::Real`: allowed tolerance for the Newton-Raphson process
    (default: `eps()`).
- `buffer::Union{Nothing, ImpactMonitoringBuffer}`: pre-allocated memory
    (default: `nothing`).
"""
function closeapproaches(
        IM::AbstractIMProblem{D, T}, VA::VirtualAsteroid{T, U},
        nyears::T, params::Parameters{T}; R_TP::T = 0.2,
        ctol::T = T(Inf), newtoniter::Int = 25, nrabstol::T = eps(T),
        buffer::Union{Nothing, ImpactMonitoringBuffer{T, U}} = nothing
    ) where {D, T <: Real, U <: Number}
    # Unpack
    @unpack target = IM
    @unpack abstol, maxsteps, order, eph_su = params
    # Dynamical model
    dynamics = dynamicalmodel(IM)
    # Initial conditions
    t0, tmax = zero(T), nyears * yr
    q0 = initialcondition(VA)
    # Impact monitoring buffer
    if isnothing(buffer)
        buffer = ImpactMonitoringBuffer(IM, q0, nyears, params)
    end
    @unpack prop, teph, tvS, xvS, gvS = buffer
    @unpack cache, dparams = prop
    @unpack xaux, t, x, dx, rv, parse_eqs = cache
    dparams.jd0 = epoch(VA) + JD_J2000

    # Specialized version of the root-finding method of taylorinteg
    x0 = deepcopy(q0)
    update_cache!(cache, t0, x0)
    sign_tstep = copysign(1, tmax - t0)

    # Some auxiliary arrays for root-finding/event detection/Poincaré
    # surface of section evaluation
    R_P = radius(target)
    f_tupl, g_tupl = closeapproach(x, teph, dparams, t; R_TP, R_P)
    f_tupl, g_tupl = (false, f_tupl[2]), (false, g_tupl[2])
    f_tupl_old, g_tupl_old = deepcopy(f_tupl), deepcopy(g_tupl)
    δt = zero(x[1])
    δt_old = zero(x[1])

    x_dx = vcat(x, dx)
    f_df = vcat(f_tupl[2], f_tupl_old[2])
    g_dg = vcat(g_tupl[2], g_tupl_old[2])
    x_dx_val = TS.evaluate(x_dx)
    f_df_val = vcat(TS.evaluate(f_tupl[2]), TS.evaluate(f_tupl_old[2]))
    g_dg_val = vcat(TS.evaluate(g_tupl[2]), TS.evaluate(g_tupl_old[2]))

    CAs = Vector{CloseApproach{T, U}}(undef, 0)

    # Integration
    nsteps = 1
    nevents = 1 # number of detected events
    while sign_tstep * t0 < sign_tstep * tmax
        δt_old = δt
        δt = taylorstep!(Val(parse_eqs), dynamics, t, x, dx, xaux, abstol, dparams, rv) # δt is positive!
        # Below, δt has the proper sign according to the direction of the integration
        if iszero(δt)
            @warn("The step-size is zero; exiting.")
            break
        end
        δt = sign_tstep * min(δt, sign_tstep * (tmax - t0))
        # New initial condition
        TS.evaluate!(x, δt, x0)
        # Root-finding
        f_tupl, g_tupl = closeapproach(x, teph, dparams, t; R_TP, R_P)
        nevents = findroot!(t, x, dx, f_tupl_old, f_tupl, 0, tvS, xvS, gvS,
                            t0, δt_old, x_dx, x_dx_val, f_df, f_df_val, nrabstol,
                            newtoniter, nevents)
        if nevents > 1
            # Time at surface crossing
            t_CA = epoch(VA) + tvS[1]
            # Asteroid's planetocentric state vector
            xap = xvS[:, 1] - target(t_CA)
            # Planetocentric keplerian elements
            kep = KeplerianElements{T, U}(
                gm(target),
                cte(t_CA) + MJD2000,
                :equatorial,
                cartesian2keplerian(xap, t_CA + MJD2000; μ = gm(target)),
                SMatrix{6, 6, T}(fill(NaN, 6, 6))
            )
            # Time at close approach
            t_CA = timeperipass(kep) - MJD2000
            # Asteroid's planetocentric state vector
            xap = kep(t_CA + MJD2000)
            # Asteroid's planetocentric semimajor axis
            a = semimajoraxis(xap..., gm(target), zero(T))
            if a < 0
                # Planet's heliocentric state vector
                xps = target(t_CA) - eph_su(t_CA)
                # BPlane in Öpik's coordinates
                tp = bopik(xap, xps, target)
            else
                # Modified target plane
                tp = mtp(xap, target)
            end
            # Close approach
            CA = CloseApproach(VA, t_CA, targetplane(tp)..., typeof(tp))
            push!(CAs, CA)
            @warn("Integration entered the physical radius of the target; exiting.")
            break
        end
        nevents = findroot!(t, x, dx, g_tupl_old, g_tupl, 0, tvS, xvS, gvS,
                            t0, δt_old, x_dx, x_dx_val, g_dg, g_dg_val, nrabstol,
                            newtoniter, nevents)
        if nevents > 1
            # Reset events counter
            nevents = 1
            # Time at close approach
            t_CA = epoch(VA) + tvS[1]
            # Asteroid's planetocentric state vector
            xap = xvS[:, 1] - target(t_CA)
            # Asteroid's planetocentric semimajor axis
            a = semimajoraxis(xap..., gm(target), zero(T))
            if a < 0
                # Planet's heliocentric state vector
                xps = target(t_CA) - eph_su(t_CA)
                # BPlane in Öpik's coordinates
                tp = bopik(xap, xps, target)
            else
                # Modified target plane
                tp = mtp(xap, target)
            end
            # Close approach
            CA = CloseApproach(VA, t_CA, targetplane(tp)..., typeof(tp))
            push!(CAs, CA)
            !isconvergent(CA, ctol) && break
        end
        f_tupl_old, g_tupl_old = deepcopy(f_tupl), deepcopy(g_tupl)
        # Update time
        t0 += δt
        update_cache!(cache, t0, x0)
        nsteps += 1
        if nsteps > maxsteps
            @warn("Maximum number of integration steps reached; exiting.")
            break
        end
    end

    return CAs
end