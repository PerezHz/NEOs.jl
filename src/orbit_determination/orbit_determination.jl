include("osculating.jl")
include("tooshortarc.jl")
include("gauss_method.jl")

function issinglearc(radec::Vector{RadecMPC{T}}, arc::Day = Day(30)) where {T <: AbstractFloat}
    for i in 2:length(radec)
        if date(radec[i]) - date(radec[i-1]) > arc
            return false
        end
    end
    return true
end

@doc raw"""
    orbitdetermination(radec::Vector{RadecMPC{T}}, params::Parameters{T}; kwargs...) where {T <: AbstractFloat}

Initial Orbit Determination (IOD) routine.

# Arguments

- `radec::Vector{RadecMPC{T}}`: vector of observations.
- `params::Parameters{T}`: see [`Parameters`](@ref).

# Keyword arguments

- `maxiter::Int = 200`:  maximum number of iterations.
- `varorder::Int`: order of jet transport perturbation. 
- `max_triplets::Int`: maximum number of triplets.
- `Q_max::T`: the maximum nrms that is considered a good enough orbit.
- `niter::Int`: number of iterations for Newton's method.
- `max_per::T =  18.0`: maximum allowed rejection percentage.
"""
function orbitdetermination(radec::Vector{RadecMPC{T}}, params::Parameters{T}; maxiter::Int = 200,
                            max_triplets::Int = 10, Q_max::T = 10., niter::Int = 5,
                            varorder::Int = 5, max_per::T = 18.0) where {T <: AbstractFloat}
    
    # Allocate memory for output
    sol = zero(NEOSolution{T, T})
    # Maximum number of steps
    params = Parameters(params; maxsteps = adaptative_maxsteps(radec))
    # Eliminate observatories without coordinates 
    filter!(x -> hascoord(observatory(x)), radec)
    # Cannot handle zero observations or multiple arcs
    if iszero(length(radec)) || !issinglearc(radec)
        return sol
    end
    # Reduce observation nights by linear regression
    gdf, cdf = reduce_nights(radec)
    # Case 1: Too Short Arc (TSA)
    if length(gdf) < 3 || numberofdays(radec) < 1
        sol = tooshortarc(radec, gdf, cdf, params; maxiter)
    # Case 2: Gauss Method
    else
        sol = gaussinitcond(radec, gdf, cdf, params; max_triplets, Q_max, niter,
                            varorder)
        if nrms(sol) > 1.0 && !iszero(sol)
            sol = gauss_refinement(radec, sol, params; max_per, niter, varorder)
        end
    end

    return sol
end
    