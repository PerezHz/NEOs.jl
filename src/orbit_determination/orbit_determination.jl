include("osculating.jl")
include("tooshortarc.jl")
include("gauss_method.jl")

@doc raw"""
    issinglearc(radec::Vector{RadecMPC{T}}, arc::Day = Day(30)) where {T <: AbstractFloat}

Check whether `radec` is a single observational arc, i.e. no two consecutive observations
are more than `arc` days appart. The function assumes `radec` is sorted.
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
    orbitdetermination(radec::Vector{RadecMPC{T}}, params::Parameters{T};
                       kwargs...) where {T <: AbstractFloat}

Initial Orbit Determination (IOD) routine.

# Arguments

- `radec::Vector{RadecMPC{T}}`: vector of observations.
- `params::Parameters{T}`: see [`Parameters`](@ref).

# Keyword arguments

- `maxiter::Int = 200`:  maximum number of iterations.
- `max_per::T =  18.0`: maximum allowed rejection percentage.
"""
function orbitdetermination(radec::Vector{RadecMPC{T}}, params::Parameters{T};
                            maxiter::Int = 200, max_per::T = 18.0) where {T <: AbstractFloat}
    
    # Allocate memory for output
    sol = zero(NEOSolution{T, T})
    # Maximum number of steps
    params = Parameters(params; maxsteps = adaptative_maxsteps(radec))
    # Eliminate observatories without coordinates 
    filter!(x -> hascoord(observatory(x)), radec)
    # Cannot handle zero observations or multiple arcs
    if iszero(length(radec)) || !issinglearc(radec)
        return sol::NEOSolution{T, T}
    end
    # Reduce observation nights by linear regression
    nights = reduce_nights(radec)
    # Case 1: Too Short Arc (TSA)
    if length(nights) < 3 || numberofdays(radec) < 1
        sol = tooshortarc(radec, gdf, cdf, params; maxiter)
    # Case 2: Gauss Method
    else
        sol = gaussinitcond(radec, nights, params)
    end
    # Outlier rejection (if needed)
    if nrms(sol) > 1.0 && !iszero(sol)
        sol = gauss_refinement(radec, sol, params; max_per, niter, varorder)
    end

    return sol::NEOSolution{T, T}
end
    