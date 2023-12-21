include("b_plane.jl")
include("least_squares.jl")

@doc raw"""
    NEOSolution{T <: Real, U <: Number}

The outcome of the orbit determination process for a NEO.

# Fields

- `nights::Vector{ObservationNight{T}}`: vector of observation nights.
- `bwd/fwd::TaylorInterpolant{T, U, 2}`: backward (forward) integration.
- `t_bwd/t_fwd::Vector{U}`: time of Earth close approach.
- `x_bwd/x_fwd::Vector{U}`: state vector at Earth close approach.
- `g_bwd/g_fwd::Vector{U}`: geocentric distance at close approach.
- `res::Vector{OpticalResidual{T, U}}`: vector of optical residuals.
- `fit::OrbitFit{T}`: least squares fit.
- `scalings::Vector{T}`: jet transport scaling factors.
"""
@auto_hash_equals struct NEOSolution{T <: Real, U <: Number}
    nights::Vector{ObservationNight{T}}
    bwd::TaylorInterpolant{T, U, 2}
    t_bwd::Vector{U}
    x_bwd::Matrix{U}
    g_bwd::Vector{U}
    fwd::TaylorInterpolant{T, U, 2}
    t_fwd::Vector{U}
    x_fwd::Matrix{U}
    g_fwd::Vector{U}
    res::Vector{OpticalResidual{T, U}}
    fit::OrbitFit{T}
    scalings::Vector{T}
    # Inner constructor
    function NEOSolution{T, U}(
        nights::Vector{ObservationNight{T}},
        bwd::TaylorInterpolant{T, U, 2}, t_bwd::Vector{U}, x_bwd::Matrix{U}, g_bwd::Vector{U},
        fwd::TaylorInterpolant{T, U, 2}, t_fwd::Vector{U}, x_fwd::Matrix{U}, g_fwd::Vector{U},
        res::Vector{OpticalResidual{T, U}}, fit::OrbitFit{T}, scalings::Vector{T}
    ) where {T <: Real, U <: Number}
        @assert bwd.t0 == fwd.t0 "Backward and forward propagation initial times must match"
        new{T, U}(
            nights, 
            bwd, t_bwd, x_bwd, g_bwd, 
            fwd, t_fwd, x_fwd, g_fwd, 
            res, fit, scalings
        )
    end
end
# Outer constructors
function NEOSolution(
    nights::Vector{ObservationNight{T}},
    bwd::TaylorInterpolant{T, U, 2}, t_bwd::Vector{U}, x_bwd::Matrix{U}, g_bwd::Vector{U},
    fwd::TaylorInterpolant{T, U, 2}, t_fwd::Vector{U}, x_fwd::Matrix{U}, g_fwd::Vector{U},
    res::Vector{OpticalResidual{T, U}}, fit::OrbitFit{T}, scalings::Vector{T}
) where {T <: Real, U <: Number}
    NEOSolution{T, U}(
        nights, 
        bwd, t_bwd, x_bwd, g_bwd, 
        fwd, t_fwd, x_fwd, g_fwd, 
        res, fit, scalings
    )
end

function NEOSolution(
    nights::Vector{ObservationNight{T}},
    bwd::TaylorInterpolant{T, U, 2}, fwd::TaylorInterpolant{T, U, 2},
    res::Vector{OpticalResidual{T, U}}, fit::OrbitFit{T}, scalings::Vector{T}
) where {T <: Real, U <: Number}
    # Backward roots
    t_bwd = Vector{U}(undef, 0)
    x_bwd = Matrix{U}(undef, 0, 0)
    g_bwd = Vector{U}(undef, 0)
    # Forward roots
    t_fwd = Vector{U}(undef, 0)
    x_fwd = Matrix{U}(undef, 0, 0)
    g_fwd = Vector{U}(undef, 0)

    NEOSolution{T, U}(
        nights, 
        bwd, t_bwd, x_bwd, g_bwd, 
        fwd, t_fwd, x_fwd, g_fwd, 
        res, fit, scalings
    )
end

# Print method for NEOSolution
# Examples:
# NEO solution with 123 residuals
function show(io::IO, x::NEOSolution{T, U}) where {T <: Real, U <: Number}
    print(io, "NEO solution with ", length(x.res), " residuals")
end

# Evaluation in time method
function (sol::NEOSolution{T, U})(t::V = sol.bwd.t0) where {T <: Real, U <: Number, V <: Number}
    if t <= sol.bwd.t0
        return sol.bwd(t)
    else
        return sol.fwd(t)
    end
end

# Evaluation in fit δs method
function evalfit(sol::NEOSolution{T, TaylorN{T}}) where {T <: Real}
    # Fit δs
    δs = sol.fit.x
    # Evaluate backward integration
    new_bwd_x = map(x -> Taylor1(x.coeffs(δs)), sol.bwd.x);
    new_bwd = TaylorInterpolant{T, T, 2}(sol.bwd.t0, sol.bwd.t, new_bwd_x)
    # Evaluate backward roots
    new_t_bwd = sol.t_bwd(δs)
    new_x_bwd = sol.x_bwd(δs)
    new_g_bwd = sol.g_bwd(δs)
    # Evaluate forward integration
    new_fwd_x = map(x -> Taylor1(x.coeffs(δs)), sol.fwd.x);
    new_fwd = TaylorInterpolant{T, T, 2}(sol.fwd.t0, sol.fwd.t, new_fwd_x)
    # Evaluate forward roots
    new_t_fwd = sol.t_fwd(δs)
    new_x_fwd = sol.x_fwd(δs)
    new_g_fwd = sol.g_fwd(δs)
    # Evaluate residuals 
    new_res = sol.res(δs)

    NEOSolution{T, T}(
        sol.nights,    
        new_bwd, new_t_bwd, new_x_bwd, new_g_bwd,
        new_fwd, new_t_fwd, new_x_fwd, new_g_fwd,
        new_res, sol.fit, sol.scalings
    )
end

# Definition of zero NEOSolution
function zero(::Type{NEOSolution{T, U}}) where {T <: Real, U <: Number}
    nights = Vector{ObservationNight{T}}(undef, 0)
    bwd = TaylorInterpolant{T, U, 2}(zero(T), zeros(T, 1), Matrix{Taylor1{U}}(undef, 0, 0))
    t_bwd = Vector{U}(undef, 0)
    x_bwd = Matrix{U}(undef, 0, 0)
    g_bwd = Vector{U}(undef, 0)
    fwd = TaylorInterpolant{T, U, 2}(zero(T), zeros(T, 1), Matrix{Taylor1{U}}(undef, 0, 0))
    t_fwd = Vector{U}(undef, 0)
    x_fwd = Matrix{U}(undef, 0, 0)
    g_fwd = Vector{U}(undef, 0)
    res = Vector{OpticalResidual{T, U}}(undef, 0)
    fit = OrbitFit{T}(false, Vector{T}(undef, 0), Matrix{T}(undef, 0, 0), :newton)
    scalings = Vector{T}(undef, 0)

    NEOSolution{T, U}(
        nights,
        bwd, t_bwd, x_bwd, g_bwd, 
        fwd, t_fwd, x_fwd, g_fwd,
        res, fit, scalings
    )
end

iszero(x::NEOSolution{T, U}) where {T <: Real, U <: Number} = x == zero(NEOSolution{T, U})

# Normalized Root Mean Square Error
function nrms(sol::NEOSolution{T, T}) where {T <: Real}
    if iszero(sol)
        return T(Inf)
    else
        return nrms(sol.res)
    end
end

sigmas(sol::NEOSolution{T, T}) where {T <: Real} = sqrt.(diag(sol.fit.Γ)) .* sol.scalings