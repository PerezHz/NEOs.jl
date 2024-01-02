include("osculating.jl")
include("tooshortarc.jl")
include("gauss_method.jl")
include("outlier_rejection.jl")

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
    istsa(sol::NEOSolution{T, T}) where {T <: AbstractFloat}
Check whether `sol` was computed via Too Short Arc (`true`) or via Gauss Method (`false`).
"""
function istsa(sol::NEOSolution{T, T}) where {T <: AbstractFloat}
    return length(sol.nights) < 3 || numberofdays(sol.nights) < 1
end

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
    orbitdetermination(radec::Vector{RadecMPC{T}}, params::Parameters{T}) where {T <: AbstractFloat}

Initial Orbit Determination (IOD) routine.

# Arguments

- `radec::Vector{RadecMPC{T}}`: vector of observations.
- `params::Parameters{T}`: see [`Parameters`](@ref).
"""
function orbitdetermination(radec::Vector{RadecMPC{T}}, params::Parameters{T}) where {T <: AbstractFloat}
    
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
        sol = tooshortarc(radec, nights, params)
    # Case 2: Gauss Method
    else
        sol = gaussinitcond(radec, nights, params)
    end
    # Outlier rejection (if needed)
    if nrms(sol) > 1 && !iszero(sol)
        sol = outlier_rejection(radec, sol, params)
    end

    return sol::NEOSolution{T, T}
end
    