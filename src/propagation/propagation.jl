include("buffer.jl")
include("dynamicalmodels.jl")
include("jetcoeffs.jl")

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

# Root-finding function for geocentric close approaches search
function rvelea(dx, x, params, t)
    # Target plane radius
    R_TP = params.R_TP
    # Solar system ephemeris at dsj2k
    ss16asteph_t = params.sseph.ephU
    # Total number of bodies
    N = params.N
    # Earth's ephemeris
    xe = ss16asteph_t[nbodyind(N-1, ea)]
    # Asteroid's geocentric state vector
    xae = x[1:6] - xe
    # Geocentric radial velocity
    return euclid3D(cte(xae[1:3])) < R_TP, dot3D(xae[1:3], xae[4:6])
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

- `R_TP::Real`: target plane radius [au].
- `eventorder::Int`: order of the derivative of `NEOs.rvelea` whose roots are computed
    (default: `0`).
- `newtoniter::Int`: maximum Newton-Raphson iterations per detected root (default: `10`).
- `nrabstol::Real`: allowed tolerance for the Newton-Raphson process (default: `eps(T)`).
"""
function propagate_root(f::D, q0::Vector{U}, jd0::Number, tmax::T, params::Parameters{T};
                        R_TP::T = 0.2, eventorder::Int = 0, newtoniter::Int = 10,
                        nrabstol::T = eps(T)) where {D, T <: Real, U <: Number}
    # Pre-allocate memory
    _jd0_ = cte(cte(jd0))
    tlim = minmax(_jd0_, _jd0_ + tmax * yr) .- JD_J2000
    buffer = PropagationBuffer(f, q0, jd0, tlim, params; R_TP)
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
    @unpack order, abstol, maxsteps, parse_eqs = params
    @unpack dparams = buffer
    # Propagate orbit
    orbit = lyap_taylorinteg(f, q0, zero(T), tmax * yr, order, abstol, dparams;
                             maxsteps, parse_eqs)

    return orbit.t, orbit.x, orbit.λ
end
