include("asteroid_dynamical_models.jl")
include("jetcoeffs.jl")
include("serialization.jl")
include("parameters.jl")

@doc raw"""
    rvelea(dx, x, params, t)

Return `true` and the asteroid's radial velocity with respect to the Earth.

# Arguments

- `dx`: asteroid's velocities.
- `x`: asteroid's degrees of freedom.
- `params`: parameters (ephemeris + accelerations + newtonian N body potential + julian date of start time + matrix of extended body interactions + number of bodies + mass parameters).
- `t`: time.
"""
function rvelea(dx, x, params, t)

    jd0 = params[4]                                 # Julian date of start time
    dsj2k = t + (jd0 - JD_J2000)                    # Days since J2000
    ss16asteph_t = evaleph(params[1], dsj2k, x[1])  # Evaluate ephemeris at dsj2k
    N = params[6]                                   # Total number of bodies
    xe = ss16asteph_t[nbodyind(N-1,ea)]             # Earth's ephemeris

    return true, (x[1]-xe[1])*(x[4]-xe[4]) + (x[2]-xe[2])*(x[5]-xe[5]) + (x[3]-xe[3])*(x[6]-xe[6])
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
    propagate_params(jd0::T, tspan::T, q0::Vector{U}; μ_ast::Vector = μ_ast343_DE430[1:end], order::Int = order,
                     abstol::T = abstol) where {T <: Real, U <: Number}

Return the parameters needed for `propagate`, `propagate_root` and `propagate_lyap`.

# Arguments

- `jd0::T`: initial Julian date.
- `tspan::T`: time span of the integration [in years].
- `q0::Vector{U}`: vector of initial conditions.
- `μ_ast::Vector`: vector of gravitational parameters.
- `order::Int`: order of the Taylor expansions to be used in the integration.
- `abstol::T`: absolute tolerance.
"""
function propagate_params(jd0::T, tspan::T, q0::Vector{U}; μ_ast::Vector = μ_ast343_DE430[1:end], order::Int = order,
                          abstol::T = abstol) where {T <: Real, U <: Number}

    # Time limits [days since J2000]
    days_0, days_f = minmax(jd0, jd0 + tspan*yr) .- JD_J2000

    # Load Solar System ephemeris [au, au/day]
    _sseph = convert(T, loadpeeph(sseph, days_0, days_f))

    # Load accelerations
    _acceph = convert(T, loadpeeph(acceph, days_0, days_f))

    # Load Newtonian potentials
    _poteph = convert(T, loadpeeph(poteph, days_0, days_f))

    # Number of massive bodies
    Nm1 = numberofbodies(_sseph)

    # Number of bodies, including NEA
    N = Nm1 + 1

    # Vector of G*m values
    μ = convert(Vector{T}, vcat( μ_DE430[1:11], μ_ast[1:Nm1-11], zero(T) ) )

    # Check: number of SS bodies (N) in ephemeris must be equal to length of GM vector (μ)
    @assert N == length(μ) "Total number of bodies in ephemeris must be equal to length of GM vector μ"

    # Interaction matrix with flattened bodies
    UJ_interaction = fill(false, N)

    # Turn on Earth interaction
    UJ_interaction[ea] = true

    # Initial time
    _t0 = zero(T)

    # Final time
    _tmax = _t0 + tspan*yr

    # Initial conditions
    _q0 = one(T) * q0

    # Vector of parameters for neosinteg
    _params = (_sseph, _acceph, _poteph, jd0, UJ_interaction, N, μ)

    return _q0, _t0, _tmax, _params

end

@doc raw"""
    propagate(dynamics::D, jd0::T, tspan::T, q0::Vector{U},
              params::NEOParameters{T}) where {T <: Real, U <: Number, D}

Integrate the orbit of a NEO via the Taylor method.

# Arguments

- `dynamics::D`: dynamical model function.
- `jd0::T`: initial Julian date.
- `tspan::T`: time span of the integration [in years].
- `q0::Vector{U}`: vector of initial conditions.
- `params::NEOParameters{T}`: see [`NEOParameters`](@ref).
"""
function propagate(dynamics::D, jd0::T, tspan::T, q0::Vector{U},
                   params::NEOParameters{T}) where {T <: Real, U <: Number, D}

    # Unfold
    maxsteps, μ_ast, order, abstol, parse_eqs = params.maxsteps, params.μ_ast, params.order,
                                                params.abstol, params.parse_eqs

    # Check order
    @assert order <= get_order(sseph.x[1]) "order ($(order)) must be less or equal than SS ephemeris order ($(get_order(sseph.x[1])))"

    # Parameters for taylorinteg
    _q0, _t0, _tmax, _params = propagate_params(jd0, tspan, q0; μ_ast = μ_ast, order = order, abstol = abstol)

    # Propagate orbit

    @time tv, xv, psol = taylorinteg(dynamics, _q0, _t0, _tmax, order, abstol, Val(true), _params;
                                     maxsteps = maxsteps, parse_eqs = parse_eqs)

    return TaylorInterpolant(jd0 - JD_J2000, tv .- tv[1], psol)
end

@doc raw"""
    propagate_root(dynamics::D, jd0::T, tspan::T, q0::Vector{U},
                   params::NEOParameters{T}; eventorder::Int = 0, newtoniter::Int = 10,
                   nrabstol::T = eps(T)) where {T <: Real, U <: Number, D}

Integrate the orbit of a NEO via the Taylor method while finding the zeros of
`NEOs.rvelea`.

# Arguments

- `dynamics::D`: dynamical model function.
- `jd0::T`: initial Julian date.
- `tspan::T`: time span of the integration [in years].
- `q0::Vector{U}`: vector of initial conditions.
- `params::NEOParameters{T}`: see [`NEOParameters`](@ref).

# Keyword arguments

- `eventorder::Int`: order of the derivative of `rvelea` whose roots are computed.
- `newtoniter::Int`: maximum Newton-Raphson iterations per detected root.
- `nrabstol::T`: allowed tolerance for the Newton-Raphson process.
"""
function propagate_root(dynamics::D, jd0::T, tspan::T, q0::Vector{U},
                        params::NEOParameters{T}; eventorder::Int = 0, newtoniter::Int = 10,
                        nrabstol::T = eps(T)) where {T <: Real, U <: Number, D}

    # Unfold
    maxsteps, μ_ast, order, abstol, parse_eqs = params.maxsteps, params.μ_ast, params.order,
                                                params.abstol, params.parse_eqs

    # Check order
    @assert order <= get_order(sseph.x[1]) "order ($(order)) must be less or equal than SS ephemeris order ($(get_order(sseph.x[1])))"

    # Parameters for neosinteg
    _q0, _t0, _tmax, _params = propagate_params(jd0, tspan, q0; μ_ast = μ_ast, order = order, abstol = abstol)

    # Propagate orbit
    @time tv, xv, psol, tvS, xvS, gvS = taylorinteg(dynamics, rvelea, _q0, _t0, _tmax, order, abstol, Val(true), _params;
                                                    maxsteps, parse_eqs, eventorder, newtoniter, nrabstol)

    return TaylorInterpolant(jd0 - JD_J2000, tv .- tv[1], psol), tvS, xvS, gvS

end

@doc raw"""
    propagate_lyap(dynamics::D, jd0::T, tspan::T, q0::Vector{U},
                   params::NEOParameters{T}) where {T <: Real, U <: Number}

Compute the Lyapunov spectrum of a NEO.

# Arguments

- `dynamics::D`: dynamical model function.
- `jd0::T`: initial Julian date.
- `tspan::T`: time span of the integration [in Julian days].
- `q0::Vector{U}`: vector of initial conditions.
- `params::NEOParameters{T}`: see [`NEOParameters`](@ref).
"""
function propagate_lyap(dynamics::D, jd0::T, tspan::T, q0::Vector{U},
                        params::NEOParameters{T}) where {T <: Real, U <: Number, D}

    # Unfold
    maxsteps, μ_ast, order, abstol, parse_eqs = params.maxsteps, params.μ_ast, params.order,
                                                params.abstol, params.parse_eqs

    # Check order
    @assert order <= get_order(sseph.x[1]) "order ($(order)) must be less or equal than SS ephemeris order ($(get_order(sseph.x[1])))"

    # Parameters for taylorinteg
    _q0, _t0, _tmax, _params = propagate_params(jd0, tspan, q0; μ_ast = μ_ast, order = order, abstol = abstol)

    # Propagate orbit
    @time sol = lyap_taylorinteg(dynamics, _q0, _t0, _tmax, order, abstol, _params;
                                 maxsteps = maxsteps, parse_eqs = parse_eqs)

    return sol

end