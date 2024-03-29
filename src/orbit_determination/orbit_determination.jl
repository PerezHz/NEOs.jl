include("osculating.jl")
include("least_squares.jl")
include("neosolution.jl")
include("tooshortarc.jl")
include("gauss_method.jl")

@doc raw"""
    issinglearc(radec::Vector{RadecMPC{T}}, arc::Day = Day(30)) where {T <: AbstractFloat}

Check whether `radec` is a single observational arc, i.e. no two consecutive observations
are more than `arc` days apart. The function assumes `radec` is sorted.
"""
function issinglearc(radec::Vector{RadecMPC{T}}, arc::Day = Day(30)) where {T <: AbstractFloat}
    for i in 2:length(radec)
        if date(radec[i]) - date(radec[i-1]) > arc
            return false
        end
    end
    return true
end

@doc raw"""
    isgauss(sol::NEOSolution{T, T}) where {T <: AbstractFloat}

Check whether `sol` was computed via [`gaussinitcond`](@ref) (`true`) or
via [`tooshortarc`](@ref) (`false`).
"""
function isgauss(tracklets::Vector{Tracklet{T}}) where {T <: AbstractFloat}
    # Observing stations
    obs = observatory.(tracklets)
    # TSA is not well suited for satellite observatories
    any(issatellite.(obs)) && return true
    # Gauss cannot handle less than 3 tracklets
    length(tracklets) < 3 && return false
    # Time span
    Δ = numberofdays(tracklets)
    # Gauss approximation does not work with less than 1 day
    return Δ > 1
end

isgauss(sol::NEOSolution{T, T}) where {T <: AbstractFloat} = isgauss(sol.tracklets)

@doc raw"""
    adaptative_maxsteps(radec::Vector{RadecMPC{T}}) where {T <: AbstractFloat}

Empirical upper bound for the number of steps used for propagation.
"""
function adaptative_maxsteps(radec::Vector{RadecMPC{T}}) where {T <: AbstractFloat}
    #=
    # Time difference [ms]
    Δ_ms =  getfield(date(radec[end]) - date(radec[1]), :value)
    # Time difference [days]
    Δ_day = Δ_ms / 86_400_000
    # Adaptative maxsteps
    if Δ_day <= 30
        return 55 - floor(Int, 5*Δ_day/6)
    else
        return ceil(Int, (Δ_day + 360)/13)
    end
    =#
    return 500
end

@doc raw"""
    orbitdetermination(radec::Vector{RadecMPC{T}}, params::NEOParameters{T};
                       dynamics::D = newtonian!) where {T <: AbstractFloat, D}

Initial Orbit Determination (IOD) routine.

# Arguments

- `radec::Vector{RadecMPC{T}}`: vector of observations.
- `params::NEOParameters{T}`: see [`NEOParameters`](@ref).
- `dynamics::D`: dynamical model.
"""
function orbitdetermination(radec::Vector{RadecMPC{T}}, params::NEOParameters{T};
                            dynamics::D = newtonian!) where {T <: AbstractFloat, D}

    # Allocate NEOSolution for output
    sol = zero(NEOSolution{T, T})
    # Maximum number of steps
    params = NEOParameters(params; maxsteps = adaptative_maxsteps(radec))
    # Eliminate observatories without coordinates
    filter!(x -> hascoord(observatory(x)), radec)
    # Cannot handle zero observations or multiple arcs
    if iszero(length(radec)) || !issinglearc(radec)
        return sol::NEOSolution{T, T}
    end
    # Reduce tracklets by polynomial regression
    tracklets = reduce_tracklets(radec)
    # Case 1: Gauss Method
    if isgauss(tracklets)
        sol = gaussinitcond(radec, tracklets, params; dynamics)
    end
    # Case 2: Too Short Arc (TSA)
    if iszero(sol) || nrms(sol) > params.gaussQmax
        sol = tooshortarc(radec, tracklets, params; dynamics)
    end

    return sol::NEOSolution{T, T}
end

@doc raw"""
    orbitdetermination(radec::Vector{RadecMPC{T}}, sol::NEOSolution{T, T}, params::NEOParameters{T};
                       dynamics::D = newtonian!, maxiter::Int = 5) where {T <: AbstractFloat, D}

Fit an orbit to `radec` using `sol` as initial condition.

# Arguments

- `radec::Vector{RadecMPC{T}}`: vector of observations.
- `sol::NEOSolution{T, T}:` least squares orbit.
- `params::NEOParameters{T}`: see [`NEOParameters`](@ref).
- `dynamics::D`: dynamical model.
- `maxiter::Int`: maximum number of iterations.
"""
function orbitdetermination(radec::Vector{RadecMPC{T}}, sol::NEOSolution{T, T}, params::NEOParameters{T};
                            dynamics::D = newtonian!, maxiter::Int = 5) where {T <: AbstractFloat, D}
    # Reduce tracklets by polynomial regression
    tracklets = reduce_tracklets(radec)
    # Epoch [julian days]
    jd0 = sol.bwd.t0 + PE.J2000
    # Barycentric initial conditions
    q0 = sol(sol.bwd.t0)
    # Scaling factors
    scalings = abs.(q0) ./ 10^6
    # Jet transport variables
    dq = scaled_variables("δx", scalings; order = params.varorder)
    # Origin
    x0 = zeros(T, 6)
    # Allocate memory
    best_Q = nrms(sol)
    # Least squares
    for _ in 1:maxiter
        # Initial conditions
        q = q0 + dq
        # Propagation and residuals
        bwd, fwd, res = propres(radec, jd0, q, params; dynamics)
        iszero(length(res)) && break
        # Orbit fit
        fit = tryls(res, x0, params.niter)
        !fit.success && break
        # Update solution
        sol = evalfit(NEOSolution(tracklets, bwd, fwd, res, fit, scalings))
        # NRMS
        Q = nrms(res, fit)
        # Convergence condition
        if abs(best_Q - Q) < 0.1
            break
        else
            best_Q = Q
        end
        # Update initial conditions
        q0 = q(fit.x)
    end
    # Case: all solutions were unsuccesful
    if isinf(best_Q)
        return zero(NEOSolution{T, T})
    # Case: at least one solution was succesful
    else
        return sol
    end
end