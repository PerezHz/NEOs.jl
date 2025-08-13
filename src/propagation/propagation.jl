include("dynamicalmodels.jl")
include("jetcoeffs.jl")

# Internal types used in the propagate* functions
# TO DO: ¿Merge these two structures into one?
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
    marsden_radial::NTuple{5, T}
end

"""
    PropagationBuffer{T <: Real, U <: Number, V <: Number} <: AbstractBuffer

Pre-allocated memory for the propagation functions: [`propagate`](@ref),
[`propagate_root`](@ref) and [`propagate_lyap`](@ref).

# Fields

- `cache::VectorCache{...}`: `TaylorIntegration` cache.
- `dparams::DynamicalParameters{T, U, V}`: parameters used by the dynamical model.
"""
struct PropagationBuffer{T <: Real, U <: Number, V <: Number} <: AbstractBuffer
    cache::VectorCache{Vector{T}, Matrix{U}, Matrix{Taylor1{U}}, Vector{Taylor1{U}},
        Taylor1{T}, Vector{Taylor1{U}}, Vector{Taylor1{U}}, RetAlloc{Taylor1{U}}, Bool}
    dparams::DynamicalParameters{T, U, V}
end

"""
    loadpeeph(eph, t0, tf)

Return a copy of `eph` in timerange `[t0, tf]`, where both times must have units
of TDB days since J2000. Currently, the only available options for `eph` are:
- `NEOs.sseph`: solar system ephemeris,
- `NEOs.acceph`: accelerations ephemeris,
- `NEOs.poteph`: newtonian potentials ephemeris.

!!! warning
    Running this function for the first time will download the `sseph_p100` artifact
    (885 MB) which can take several minutes.
"""
function loadpeeph(eph::TaylorInterpolant = sseph, t0::Real = sseph.t0,
                   tf::Real = sseph.t0 + sseph.t[end])
    @assert 0.0 ≤ t0 ≤ tf ≤ 36_525.0
    j0 = searchsortedlast(eph.t, t0)
    jf = searchsortedfirst(eph.t, tf)
    return TaylorInterpolant(eph.t0, eph.t[j0:jf], eph.x[j0:jf-1, :])
end

"""
    PropagationBuffer(dynamics, q0, jd0, tlim, params)

Return a `PropagationBuffer` object with pre-allocated memory for `propagate`.

# Arguments

- `dynamics`: dynamical model function.
- `q0::Vector{<:Number}`: vector of initial conditions.
- `jd0::Number`: initial Julian date (TDB).
- `tlim::NTuple{2, <:Real}`: ephemeris timespan [days since J2000].
- `params::Parameters{<:Real}`: see the `Propagation` section of [`Parameters`](@ref).
"""
function PropagationBuffer(dynamics::D, q0::Vector{U}, jd0::V, tlim::NTuple{2, T},
                           params::Parameters{T}) where {D, T <: Real, U <: Number,
                           V <: Number}
    # Unpack parameters
    @unpack order, μ_ast, maxsteps, marsden_radial = params
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
    μ = convert(Vector{T}, vcat( μ_DE430[1:11], μ_ast[1:Nm1-11], zero(T) ) )
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
        UJ_interaction = falses(N)
        # Turn on Earth interaction
        UJ_interaction[ea] = true
    end
    # Dynamical parameters for `propagate`
    dparams = DynamicalParameters{T, U, V}(_sseph, _ssepht, _acceph, _accepht,
        _poteph, _potepht, jd0, UJ_interaction, N, μ, marsden_radial)
    # TaylorIntegration cache
    cache = init_cache(Val(true), zero(T), q0, maxsteps, order, dynamics, dparams)

    return PropagationBuffer{T, U, V}(cache, dparams)
end

"""
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

# Return `true` and the asteroid's radial velocity with respect to the Earth
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
    # Geocentric radial velocity
    return true, (x[1]-xe[1])*(x[4]-xe[4]) + (x[2]-xe[2])*(x[5]-xe[5]) +
        (x[3]-xe[3])*(x[6]-xe[6])
end

"""
    scaled_variables(names::String, c::Vector{T}; [order::Int = 5])

Equivalent to:

`TaylorSeries.set_variables(T, names; order, numvars = length(c))`

times a scaling given by `c`.
"""
function scaled_variables(names::String = "δx", c::Vector{T} = fill(1e-6, 6);
                          order::Int = 5) where {T <: Real}
    # Set TaylorN variables
    dq = set_variables(T, names; order, numvars = length(c))
    # Scale jet transport perturbation
    for i in eachindex(dq)
        dq[i][1][i] = c[i]
    end
    return dq
end

# Check if an integration was successful
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

"""
    propagate(f, q0, jd0, tmax, params)

Integrate the dynamical model `f` starting from an initial condition `q0` [au, au/day]
at time `jd0` [julian date TDB] for a period of `tmax` [years]. For a list of parameters
see the `Propagation` section of [`Parameters`](@ref).

This function uses the Taylor Method implemented in `TaylorIntegration`.
"""
function propagate(f::D, q0::Vector{U}, jd0::Number, tmax::T,
                   params::Parameters{T}) where {D, T <: Real, U <: Number}
    # Pre-allocate memory
    _jd0_ = cte(cte(jd0))
    tlim = minmax(_jd0_, _jd0_ + tmax * yr) .- JD_J2000
    buffer = PropagationBuffer(f, q0, jd0, tlim, params)
    # Propagate orbit
    return _propagate(f, q0, jd0, tmax, buffer, params)
end

function _propagate(f::D, q0::Vector{U}, jd0::V, tmax::T, buffer::PropagationBuffer{T, U, V},
                    params::Parameters{T}) where {D, T <: Real, U <: Number, V <: Number}
    # Unpack
    @unpack abstol, maxsteps = params
    @unpack cache, dparams = buffer
    # Update reference epoch
    dparams.jd0 = jd0
    # Propagate orbit
    orbit = taylorinteg!(Val(true), f, q0, zero(T), tmax * yr, abstol, cache, dparams;
                         maxsteps)
    # Epoch (plain)
    _jd0_ = cte(cte(jd0))
    # Output
    return TaylorInterpolant{T, U, 2}(_jd0_ - JD_J2000, orbit.t, orbit.p)
end

"""
    propagate_root(f, q0, jd0, tmax, params; kwargs...)

Integrate the dynamical model `f` starting from an initial condition `q0` [au, au/day]
at time `jd0` [julian date TDB] for a period of `tmax` [years], while finding the zeros
of `NEOs.rvelea`. For a list of parameters see the `Propagation` section of
[`Parameters`](@ref).

This function uses the Taylor Method implemented in `TaylorIntegration`.

# Keyword arguments

- `eventorder::Int`: order of the derivative of `NEOs.rvelea` whose roots are computed
    (default: `0`).
- `newtoniter::Int`: maximum Newton-Raphson iterations per detected root (default: `10`).
- `nrabstol::T`: allowed tolerance for the Newton-Raphson process (default: `eps(T)`).
"""
function propagate_root(f::D, q0::Vector{U}, jd0::Number, tmax::T, params::Parameters{T};
                        eventorder::Int = 0, newtoniter::Int = 10,
                        nrabstol::T = eps(T)) where {D, T <: Real, U <: Number}
    # Pre-allocate memory
    _jd0_ = cte(cte(jd0))
    tlim = minmax(_jd0_, _jd0_ + tmax * yr) .- JD_J2000
    buffer = PropagationBuffer(f, q0, jd0, tlim, params)
    # Propagate orbit
    return _propagate_root(f, q0, jd0, tmax, buffer, params; eventorder,
        newtoniter, nrabstol)
end

function _propagate_root(f::D, q0::Vector{U}, jd0::V, tmax::T, buffer::PropagationBuffer{T, U, V},
                         params::Parameters{T}; eventorder::Int = 0, newtoniter::Int = 10,
                         nrabstol::T = eps(T)) where {D, T <: Real, U <: Number, V <: Number}
    # Unpack
    @unpack abstol, maxsteps = params
    @unpack cache, dparams = buffer
    # Update reference epoch
    dparams.jd0 = jd0
    # Propagate orbit
    orbit = taylorinteg!(Val(true), f, rvelea, q0, zero(T), tmax * yr, abstol, cache,
                         dparams; maxsteps, eventorder, newtoniter, nrabstol)
    # Epoch (plain)
    _jd0_ = cte(cte(jd0))
    # Output
    return TaylorInterpolant{T, U, 2}(_jd0_ - JD_J2000, orbit.t, orbit.p),
        orbit.tevents, orbit.xevents, orbit.gresids
end

"""
    propagate_lyap(f, q0, jd0, tmax, params)

Compute the Lyapunov spectrum of an orbit given by integrating the dynamical model `f`
starting from an initial condition `q0` [au, au/day] at time `jd0` [julian date TDB]
for a period of `tmax` [years]. For a list of parameters see the `Propagation` section
of [`Parameters`](@ref).

This function uses the Taylor Method implemented in `TaylorIntegration`.
"""
function propagate_lyap(f::D, q0::Vector{U}, jd0::Number, tmax::T,
                        params::Parameters{T}) where {D, T <: Real, U <: Number}
    # Pre-allocate memory
    _jd0_ = cte(cte(jd0))
    tlim = minmax(_jd0_, _jd0_ + tmax * yr) .- JD_J2000
    buffer = PropagationBuffer(f, q0, jd0, tlim, params)
    # Unpack
    @unpack order, abstol, maxsteps = params
    @unpack dparams = buffer
    # Propagate orbit
    orbit = lyap_taylorinteg(f, q0, zero(T), tmax * yr, order, abstol, dparams; maxsteps)

    return orbit.t, orbit.x, orbit.λ
end
