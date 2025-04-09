include("asteroid_dynamical_models.jl")
include("jetcoeffs.jl")

mutable struct DynamicalParameters{T <: Real, U <: Number, V <: Number}
    sseph::TaylorInterpolant{T, T, 2, Vector{T}, Matrix{Taylor1{T}}}
    ssepht::Vector{Taylor1{U}}
    acceph::Union{Nothing, TaylorInterpolant{T, T, 2, Vector{T}, Matrix{Taylor1{T}}}}
    accepht::Union{Nothing, Vector{Taylor1{U}}}
    poteph::Union{Nothing, TaylorInterpolant{T, T, 2, Vector{T}, Matrix{Taylor1{T}}}}
    potepht::Union{Nothing, Vector{Taylor1{U}}}
    jd0::V
    UJ_interaction::Union{Nothing, Vector{Bool}}
    N::Int
    μ::Vector{T}
end

struct PropagationBuffer{T <: Real, U <: Number, V <: Number}
    cache::VectorCache{Vector{T}, Matrix{U}, Matrix{Taylor1{U}}, Vector{Taylor1{U}},
        Taylor1{T}, Vector{Taylor1{U}}, Vector{Taylor1{U}}, RetAlloc{Taylor1{U}}, Bool}
    params::DynamicalParameters{T, U, V}
end

@doc raw"""
    PropagationBuffer(dynamics::D, jd0::V, tlim::Tuple{T, T}, q0::Vector{U},
        params::NEOParameters{T}) where {D, T <: Real, U <: Number, V <: Number}

Return a `PropagationBuffer` object with pre-allocated memory for `propagate`.

## Arguments

- `dynamics::D`: dynamical model function.
- `jd0::V`: initial Julian date (TDB).
- `tlim::Tuple{T, T}`: ephemeris timespan [in days since J2000].
- `q0::Vector{U}`: vector of initial conditions.
- `params::NEOParameters{T}`: see [`NEOParameters`](@ref).
"""
function PropagationBuffer(dynamics::D, jd0::V, tlim::Tuple{T, T}, q0::Vector{U},
            params::NEOParameters{T}) where {D, T <: Real, U <: Number, V <: Number}
    # Unfold parameters
    maxsteps, order = params.maxsteps, params.order
    # Check order
    @assert order <= SSEPHORDER "order ($order) must be less or equal than SS ephemeris \
        order ($SSEPHORDER)"
    # Time limits [days since J2000]
    days_0, days_f = minmax(tlim[1], tlim[2])
    # Load Solar System ephemeris [au, au/day]
    x = Taylor1( q0[1], order )
    _sseph = convert(T, loadpeeph(sseph, days_0, days_f))
    _ssepht = [zero(x) for _ in axes(_sseph.x, 2)]
    # Number of massive bodies
    Nm1 = numberofbodies(_sseph)
    # Number of bodies, including NEA
    N = Nm1 + 1
    # Vector of G*m values
    μ = convert(Vector{T}, vcat( μ_DE430[1:11], params.μ_ast[1:Nm1-11], zero(T) ) )
    # Check: number of SS bodies (N) in ephemeris must be equal to length of GM vector (μ)
    @assert N == length(μ) "Total number of bodies in ephemeris must be equal to length \
        of GM vector μ"
    # Accelerations, Newtonian potentials and interaction matrix with flattened bodies
    if dynamics == newtonian!
        _acceph, _accepht = nothing, nothing
        _poteph, _potepht = nothing, nothing
        UJ_interaction = nothing
    else
        _acceph = convert(T, loadpeeph(acceph, days_0, days_f))
        _accepht = [zero(x) for _ in axes(_acceph.x, 2)]
        _poteph = convert(T, loadpeeph(poteph, days_0, days_f))
        _potepht = [zero(x) for _ in axes(_poteph.x, 2)]
        UJ_interaction = fill(false, N)
        # Turn on Earth interaction
        UJ_interaction[ea] = true
    end
    # Dynamical parameters for `propagate`
    _params = DynamicalParameters{T, U, V}(_sseph, _ssepht, _acceph, _accepht,
        _poteph, _potepht, jd0, UJ_interaction, N, μ)
    # TaylorIntegration cache
    cache = init_cache(Val(true), zero(T), q0, maxsteps, order, dynamics, _params)

    return PropagationBuffer{T, U, V}(cache, _params)
end

@doc raw"""
    rvelea(dx, x, params, t)

Return `true` and the asteroid's radial velocity with respect to the Earth.

## Arguments

- `dx`: asteroid's velocities.
- `x`: asteroid's degrees of freedom.
- `params`: dynamical parameters (see [`DynamicalParameters`](@ref)).
- `t`: time.
"""
function rvelea(dx, x, params, t)
    # Julian date (TDB) of start time
    jd0 = params.jd0
    # Days since J2000.0 = 2.451545e6
    dsj2k = t + (jd0 - JD_J2000)
    # Solar system ephemeris at dsj2k
    ss16asteph_t = params.ssepht
    evaleph!(ss16asteph_t, params.sseph, dsj2k)
    # Total number of bodies
    N = params.N
    # Earth's ephemeris
    xe = ss16asteph_t[nbodyind(N-1, ea)]

    return true, (x[1]-xe[1])*(x[4]-xe[4]) + (x[2]-xe[2])*(x[5]-xe[5]) +
        (x[3]-xe[3])*(x[6]-xe[6])
end

@doc raw"""
    rvelea(eph, params, t)

Return the geocentric radial velocity of `eph` at time `t`, using
the Earth's ephemerides in `params`.
"""
function rvelea(eph, params, t)
    # Geocentric state vector
    rv = eph(t) - params.eph_ea(t)
    # Derivative of geocentric distance
    return dot3D(rv[1:3], rv[4:6])
end

@doc raw"""
    scaled_variables(names::String = "δx", c::Vector{T} = fill(1e-6, 6);
        order::Int = 5) where {T <: Real}

Equivalent to `TaylorSeries.set_variables` times a scaling given by `c`.
"""
function scaled_variables(names::String = "δx", c::Vector{T} = fill(1e-6, 6);
            order::Int = 5) where {T <: Real}
    # Set TaylorN variables
    dq = set_variables(T, names, order = order, numvars = length(c))
    # Scale jet transport perturbation
    for i in eachindex(dq)
        dq[i][1][i] = c[i]
    end
    return dq
end

function issuccessfulprop(sol::TaylorInterpolant, t::T; tol::T = 10.0) where {T <: Real}
    # Zero TaylorInterpolant
    iszero(sol) && return false
    # Forward integration
    if issorted(sol.t)
        # Insufficient steps
        sol.t[end] < t && return false
        # Step that covers t
        i = searchsortedfirst(sol.t, t) - 1
    # Backward integration
    elseif issorted(sol.t, rev = true)
        # Insufficient steps
        sol.t[end] > t && return false
        # Step that covers t
        i = searchsortedfirst(sol.t, t, lt = !isless) - 1
    # This case should never happen
    else
        return false
    end
    # All coefficients are below tol
    return all( norm.(view(sol.x, 1:i, :), Inf) .< tol )
end

@doc raw"""
    propagate(dynamics::D, jd0::V, tspan::T, q0::Vector{U},
        params::NEOParameters{T}) where {T <: Real, U <: Number, V <: Number, D}

Integrate an orbit via the Taylor method. The initial Julian date `jd0` is assumed
to be in TDB time scale.

## Arguments

- `dynamics::D`: dynamical model function.
- `jd0::V`: initial Julian date (TDB).
- `tspan::T`: time span of the integration [in years].
- `q0::Vector{U}`: vector of initial conditions.
- `params::NEOParameters{T}`: see [`NEOParameters`](@ref).
"""
function propagate(dynamics::D, jd0::V, tspan::T, q0::Vector{U},
            params::NEOParameters{T}) where {T <: Real, U <: Number, V <: Number, D}
    # Pre-allocate memory
    _jd0_ = cte(cte(jd0))
    tlim = minmax(_jd0_, _jd0_ + tspan * yr) .- JD_J2000
    buffer = PropagationBuffer(dynamics, jd0, tlim, q0, params)
    # Propagate orbit
    return _propagate(dynamics, jd0, tspan, q0, buffer, params)
end

function _propagate(
        dynamics::D, jd0::V, tspan::T, q0::Vector{U},
        buffer::PropagationBuffer{T, U, V}, params::NEOParameters{T}
    ) where {T <: Real, U <: Number, V <: Number, D}
    # Unfold
    abstol, maxsteps = params.abstol, params.maxsteps
    cache, dparams = buffer.cache, buffer.params
    # Update reference epoch
    dparams.jd0 = jd0
    # Propagate orbit
    @time sol = taylorinteg!(Val(true), dynamics, q0, zero(T), tspan * yr, abstol,
        cache, dparams; maxsteps)
    # Epoch (plain)
    _jd0_ = cte(cte(jd0))
    # Output
    return TaylorInterpolant{T, U, 2}(_jd0_ - JD_J2000, sol.t, sol.p)
end

@doc raw"""
    propagate_root(dynamics::D, jd0::V, tspan::T, q0::Vector{U}, params::NEOParameters{T};
        kwargs...) where {T <: Real, U <: Number, V <: Number D}

Integrate an orbit via the Taylor method while finding the zeros of `NEOs.rvelea`.

## Arguments

- `dynamics::D`: dynamical model function.
- `jd0::V`: initial Julian date (TDB).
- `tspan::T`: time span of the integration [in years].
- `q0::Vector{U}`: vector of initial conditions.
- `params::NEOParameters{T}`: see [`NEOParameters`](@ref).

## Keyword arguments

- `eventorder::Int`: order of the derivative of `rvelea` whose roots are computed.
- `newtoniter::Int`: maximum Newton-Raphson iterations per detected root.
- `nrabstol::T`: allowed tolerance for the Newton-Raphson process.
"""
function propagate_root(dynamics::D, jd0::V, tspan::T, q0::Vector{U},
            params::NEOParameters{T}; eventorder::Int = 0, newtoniter::Int = 10,
            nrabstol::T = eps(T)) where {T <: Real, U <: Number, V <: Number, D}
    # Pre-allocate memory
    _jd0_ = cte(cte(jd0))
    tlim = minmax(_jd0_, _jd0_ + tspan * yr) .- JD_J2000
    buffer = PropagationBuffer(dynamics, jd0, tlim, q0, params)
    # Propagate orbit
    return _propagate_root(dynamics, jd0, tspan, q0, buffer, params;
        eventorder, newtoniter, nrabstol)
end

function _propagate_root(
        dynamics::D, jd0::V, tspan::T, q0::Vector{U},
        buffer::PropagationBuffer{T, U, V}, params::NEOParameters{T};
        eventorder::Int = 0, newtoniter::Int = 10, nrabstol::T = eps(T)
    ) where {T <: Real, U <: Number, V <: Number, D}
    # Unfold
    abstol, maxsteps = params.abstol, params.maxsteps
    cache, dparams = buffer.cache, buffer.params
    # Update reference epoch
    dparams.jd0 = jd0
    # Propagate orbit
    @time sol = taylorinteg!(Val(true), dynamics, rvelea, q0, zero(T), tspan * yr,
        abstol, cache, dparams; maxsteps, eventorder, newtoniter, nrabstol)
    # Epoch (plain)
    _jd0_ = cte(cte(jd0))
    # Output
    return TaylorInterpolant{T, U, 2}(_jd0_ - JD_J2000, sol.t, sol.p), sol.tevents,
        sol.xevents, sol.gresids
end

@doc raw"""
    propagate_lyap(dynamics::D, jd0::V, tspan::T, q0::Vector{U},
        params::NEOParameters{T}) where {T <: Real, U <: Number, V <: Number, D}

Compute the Lyapunov spectrum of an orbit.

# Arguments

- `dynamics::D`: dynamical model function.
- `jd0::V`: initial Julian date (TDB).
- `tspan::T`: time span of the integration [in Julian days].
- `q0::Vector{U}`: vector of initial conditions.
- `params::NEOParameters{T}`: see [`NEOParameters`](@ref).
"""
function propagate_lyap(dynamics::D, jd0::V, tspan::T, q0::Vector{U},
            params::NEOParameters{T}) where {T <: Real, U <: Number, V <: Number, D}
    # Pre-allocate memory
    _jd0_ = cte(cte(jd0))
    tlim = minmax(_jd0_, _jd0_ + tspan * yr) .- JD_J2000
    buffer = PropagationBuffer(dynamics, jd0, tlim, q0, params)
    # Unfold
    order, abstol, maxsteps = params.order, params.abstol, params.maxsteps
    dparams = buffer.params
    # Propagate orbit
    @time sol = lyap_taylorinteg(dynamics, q0, zero(T), tspan * yr, order,
          abstol, dparams; maxsteps)

    return sol.t, sol.x, sol.λ
end
