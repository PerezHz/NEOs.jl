include("osculating.jl")
include("least_squares.jl")
include("neosolution.jl")
include("admissibleregion.jl")
include("tooshortarc.jl")
include("gauss_method.jl")

@doc raw"""
    issinglearc(radec::Vector{RadecMPC{T}}, arc::Day = Day(30)) where {T <: AbstractFloat}

Check whether `radec` is a single observational arc, i.e. no two consecutive observations
are more than `arc` days apart. The function assumes `radec` is sorted.
"""
function issinglearc(radec::Vector{RadecMPC{T}}, arc::Day = Day(30)) where {T <: AbstractFloat}
    return all(diff(date.(radec)) .< arc)
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
    # Reference epoch [julian days]
    jd0 = sol.bwd.t0 + PE.J2000
    # Plain barycentric initial condition
    q0 = sol(sol.bwd.t0)
    # Scaling factors
    scalings = abs.(q0) ./ 10^6
    # Jet transport variables
    dq = scaled_variables("δx", scalings; order = params.varorder)
    # Jet Transport initial condition
    q = q0 + dq
    # Jet Transport Least Squares
    return jtls(radec, tracklets, jd0, q, 1, length(tracklets), params; maxiter, dynamics)
end