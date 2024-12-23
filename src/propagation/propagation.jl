include("asteroid_dynamical_models.jl")
include("jetcoeffs.jl")
include("parameters.jl")

mutable struct DynamicalParameters{T <: Real, V <: Number}
    sseph::TaylorInterpolant{T, T, 2, Vector{T}, Matrix{Taylor1{T}}}
    acceph::Union{Nothing, TaylorInterpolant{T, T, 2, Vector{T}, Matrix{Taylor1{T}}}}
    poteph::Union{Nothing, TaylorInterpolant{T, T, 2, Vector{T}, Matrix{Taylor1{T}}}}
    jd0::V
    UJ_interaction::Union{Nothing, Vector{Bool}}
    N::Int
    μ::Vector{T}
end

struct PropagationBuffer{T <: Real, U <: Number, V <: Number}
    t::Taylor1{T}
    x::Vector{Taylor1{U}}
    dx::Vector{Taylor1{U}}
    q0::Vector{U}
    rv::RetAlloc{Taylor1{U}}
    params::DynamicalParameters{T, V}
end

@doc raw"""
    PropagationBuffer(dynamics::D, jd0::V, tlim::Tuple{T, T}, q0::Vector{U},
                      params::NEOParameters{T}) where {T <: Real, U <: Number, V <: Number, D}

Return a `PropagationBuffer` object with pre-allocated memory for `propagate`.

# Arguments

- `dynamics::D`: dynamical model function.
- `jd0::V`: initial Julian date (TDB).
- `tlim::Tuple{T, T}`: ephemeris timespan [in days since J2000].
- `q0::Vector{U}`: vector of initial conditions.
- `params::NEOParameters{T}`: see [`NEOParameters`](@ref).
"""
function PropagationBuffer(dynamics::D, jd0::V, tlim::Tuple{T, T}, q0::Vector{U},
                           params::NEOParameters{T}) where {T <: Real, U <: Number, V <: Number, D}
    # Unfold parameters
    order = params.order
    # Check order
    @assert order <= SSEPHORDER "order ($order) must be less or equal than SS ephemeris order ($SSEPHORDER)"
    # Initialize the vector of Taylor1 expansions
    dof = length(q0)
    t = zero(T) + Taylor1( T, order )
    x = Vector{Taylor1{U}}(undef, dof)
    dx = Vector{Taylor1{U}}(undef, dof)
    @inbounds for i in eachindex(q0)
        @inbounds x[i] = Taylor1( q0[i], order )
        @inbounds dx[i] = Taylor1( zero(q0[i]), order )
    end
    # Time limits [days since J2000]
    days_0, days_f = minmax(tlim[1], tlim[2])
    # Load Solar System ephemeris [au, au/day]
    _sseph = convert(T, loadpeeph(sseph, days_0, days_f))
    # Number of massive bodies
    Nm1 = numberofbodies(_sseph)
    # Number of bodies, including NEA
    N = Nm1 + 1
    # Vector of G*m values
    μ = convert(Vector{T}, vcat( μ_DE430[1:11], params.μ_ast[1:Nm1-11], zero(T) ) )
    # Check: number of SS bodies (N) in ephemeris must be equal to length of GM vector (μ)
    @assert N == length(μ) "Total number of bodies in ephemeris must be equal to length of GM vector μ"
    # Accelerations, Newtonian potentials and interaction matrix with flattened bodies
    if dynamics == newtonian!
        _acceph = nothing
        _poteph = nothing
        UJ_interaction = nothing
    else
        _acceph = convert(T, loadpeeph(acceph, days_0, days_f))
        _poteph = convert(T, loadpeeph(poteph, days_0, days_f))
        UJ_interaction = fill(false, N)
        # Turn on Earth interaction
        UJ_interaction[ea] = true
    end
    # Dynamical parameters for `propagate`
    _params = DynamicalParameters{T, V}(_sseph, _acceph, _poteph, jd0, UJ_interaction, N, μ)
    # Determine if specialized jetcoeffs! method exists
    _, rv = _determine_parsing!(true, dynamics, t, x, dx, _params)

    return PropagationBuffer{T, U, V}(t, x, dx, q0, rv, _params)
end

@doc raw"""
    rvelea(dx, x, params, t)

Return `true` and the asteroid's radial velocity with respect to the Earth.

## Arguments

- `dx`: asteroid's velocities.
- `x`: asteroid's degrees of freedom.
- `params`: parameters (ephemeris + accelerations + newtonian N body potential + julian date of start time + matrix of extended body interactions + number of bodies + mass parameters).
- `t`: time.
"""
function rvelea(dx, x, params, t)

    jd0 = params.jd0                                  # Julian date (TDB) of start time
    dsj2k = t + (jd0 - JD_J2000)                      # Days since J2000
    ss16asteph_t = evaleph(params.sseph, dsj2k, x[1]) # Evaluate ephemeris at dsj2k
    N = params.N                                      # Total number of bodies
    xe = ss16asteph_t[nbodyind(N-1, ea)]              # Earth's ephemeris

    return true, (x[1]-xe[1])*(x[4]-xe[4]) + (x[2]-xe[2])*(x[5]-xe[5]) + (x[3]-xe[3])*(x[6]-xe[6])
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
function scaled_variables(names::String = "δx", c::Vector{T} = fill(1e-6, 6); order::Int = 5) where {T <: Real}
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

Integrate the a NEO orbit via the Taylor method. The initial Julian date `jd0` is assumed to
be in TDB time scale.

# Arguments

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

function _propagate(dynamics::D, jd0::V, tspan::T, q0::Vector{U}, buffer::PropagationBuffer{T, U, V},
                    params::NEOParameters{T}) where {T <: Real, U <: Number, V <: Number, D}
    # Unfold parameters
    order = params.order
    # Re-initialize the Taylor1 expansions
    buffer.t[:] = zero(T)
    buffer.t[1] = one(T)
    buffer.x .= Taylor1.( q0, order )
    buffer.dx .= Taylor1.( zero.(q0), order)
    buffer.params.jd0 = jd0
    # Propagate orbit
    @time sol = _taylorinteg!(Val(true), dynamics, buffer.t, buffer.x, buffer.dx, buffer.q0, zero(T),
                                       tspan * yr, params.abstol, buffer.rv, buffer.params;
                                       maxsteps = params.maxsteps)
    # Output
    if issorted(sol.t) || issorted(sol.t, rev = true)
        # Epoch (plain)
        _jd0_ = cte(cte(jd0))
        return TaylorInterpolant{T, U, 2}(_jd0_ - JD_J2000, sol.t, sol.p)
    else
        return zero(TaylorInterpolant{T, U, 2, SubArray{T, 1}, SubArray{Taylor1{U}, 2}})
    end
end

@doc raw"""
    propagate_root(dynamics::D, jd0::V, tspan::T, q0::Vector{U}, params::NEOParameters{T};
                   eventorder::Int = 0, newtoniter::Int = 10,
                   nrabstol::T = eps(T)) where {T <: Real, U <: Number, V <: Number D}

Integrate the orbit of a NEO via the Taylor method while finding the zeros of
`NEOs.rvelea`.

# Arguments

- `dynamics::D`: dynamical model function.
- `jd0::V`: initial Julian date (TDB).
- `tspan::T`: time span of the integration [in years].
- `q0::Vector{U}`: vector of initial conditions.
- `params::NEOParameters{T}`: see [`NEOParameters`](@ref).

# Keyword arguments

- `eventorder::Int`: order of the derivative of `rvelea` whose roots are computed.
- `newtoniter::Int`: maximum Newton-Raphson iterations per detected root.
- `nrabstol::T`: allowed tolerance for the Newton-Raphson process.
"""
function propagate_root(dynamics::D, jd0::V, tspan::T, q0::Vector{U}, params::NEOParameters{T};
                        eventorder::Int = 0, newtoniter::Int = 10,
                        nrabstol::T = eps(T)) where {T <: Real, U <: Number, V <: Number, D}
    # Pre-allocate memory
    _jd0_ = cte(cte(jd0))
    tlim = minmax(_jd0_, _jd0_ + tspan * yr) .- JD_J2000
    buffer = PropagationBuffer(dynamics, jd0, tlim, q0, params)
    # Propagate orbit
    return _propagate_root(dynamics, jd0, tspan, q0, buffer, params;
                           eventorder, newtoniter, nrabstol)
end

function _propagate_root(dynamics::D, jd0::V, tspan::T, q0::Vector{U}, buffer::PropagationBuffer{T, U, V},
                         params::NEOParameters{T}; eventorder::Int = 0, newtoniter::Int = 10,
                         nrabstol::T = eps(T)) where {T <: Real, U <: Number, V <: Number, D}
    # Unfold parameters
    order = params.order
    @assert order ≥ eventorder "`eventorder` must be less than or equal to `order`"
    # Re-initialize the Taylor1 expansions
    buffer.t[:] = zero(T)
    buffer.t[1] = one(T)
    buffer.x .= Taylor1.( q0, order )
    buffer.dx .= Taylor1.( zero.(q0), order)
    buffer.params.jd0 = jd0
    # Propagate orbit
    @time sol = _taylorinteg!(Val(true), dynamics, rvelea, buffer.t, buffer.x, buffer.dx,
          buffer.q0, zero(T), tspan * yr, params.abstol, buffer.rv, buffer.params;
          maxsteps = params.maxsteps, eventorder = eventorder, newtoniter = newtoniter,
          nrabstol = nrabstol)
    # Epoch (plain)
    _jd0_ = cte(cte(jd0))
    # Output
    return TaylorInterpolant{T, U, 2}(_jd0_ - JD_J2000, sol.t, sol.p), sol.tevents, sol.xevents, sol.gresids
end

@doc raw"""
    propagate_lyap(dynamics::D, jd0::V, tspan::T, q0::Vector{U},
                   params::NEOParameters{T}) where {T <: Real, U <: Number, V <: Number, D}

Compute the Lyapunov spectrum of a NEO.

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
    # Propagate orbit
    @time sol = lyap_taylorinteg(dynamics, q0, zero(T), tspan * yr, params.order,
          params.abstol, buffer.params; maxsteps = params.maxsteps, parse_eqs = params.parse_eqs)

    return sol.t, sol.x, sol.λ
end
