@doc raw"""
    AbstractLeastSquaresOrbit{T, U} <: AbstractOrbit{T, U}

Supertype for the least squares orbits API.
"""
abstract type AbstractLeastSquaresOrbit{T, U} <: AbstractOrbit{T, U} end

# Initialize a vector of residuals consistent with od and orbit
function init_residuals(::Type{U}, od::ODProblem{D, T},
    orbit::AbstractLeastSquaresOrbit) where {D, T <: Real, U <: Number}
    # Optical astrometry
    odradec, orbitradec = od.radec, astrometry(orbit)
    # Weights and debiasing factors
    w8s, bias = od.w8s.w8s, od.bias.bias
    # Initialize vector of residuals
    res = Vector{OpticalResidual{T, U}}(undef, length(odradec))
    for i in eachindex(odradec)
        ξ_α, ξ_δ = zero(U), zero(U)
        w_α, w_δ = w8s[i]
        μ_α, μ_δ = bias[i]
        j = findfirst(==(odradec[i]), orbitradec)
        outlier = isnothing(j) ? false : isoutlier(orbit.res[j])
        res[i] = OpticalResidual{T, U}(ξ_α, ξ_δ, w_α, w_δ, μ_α, μ_δ, outlier)
    end

    return res
end

covariance(orbit::AbstractLeastSquaresOrbit) = orbit.fit.Γ
variances(orbit::AbstractLeastSquaresOrbit) = diag(orbit.J * covariance(orbit) * orbit.J')

@doc raw"""
    LeastSquaresOrbit{D, T, U} <: AbstractLeastSquaresOrbit{T, U}

An asteroid least squares orbit.

## Fields

- `dynamics::D`: dynamical model.
- `tracklets::Vector{Tracklet{T}}`: vector of tracklets.
- `bwd/fwd::TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}}`:
    backward (forward) integration.
- `res::Vector{OpticalResidual{T, U}}`: vector of optical residuals.
- `fit::LeastSquaresFit{T}`: least squares fit.
- `J::Matrix{T}`: residuals space to barycentric coordinates jacobian.
- `qs::Matrix{T}`: history of initial conditions.
- `Qs::Vector{T}`: history of the target function.
"""
@auto_hash_equals struct LeastSquaresOrbit{D, T, U} <: AbstractLeastSquaresOrbit{T, U}
    dynamics::D
    tracklets::Vector{Tracklet{T}}
    bwd::TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}}
    fwd::TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}}
    res::Vector{OpticalResidual{T, U}}
    fit::LeastSquaresFit{T}
    J::Matrix{T}
    qs::Matrix{T}
    Qs::Vector{T}
    # Inner constructor
    function LeastSquaresOrbit{D, T, U}(dynamics::D, tracklets::Vector{Tracklet{T}},
        bwd::TaylorInterpolant{T, U, 2}, fwd::TaylorInterpolant{T, U, 2},
        res::Vector{OpticalResidual{T, U}}, fit::LeastSquaresFit{T},
        J::Matrix{T}, qs::Matrix{T}, Qs::Vector{T}) where {D, T <: Real, U <: Number}
        @assert bwd.t0 == fwd.t0 "Backward and forward integration initial \
            times must match"
        @assert nobs(tracklets) == length(res) "Number of observations must \
            match number of residuals"
        _bwd_ = TaylorInterpolant(bwd.t0, bwd.t, collect(bwd.x))
        _fwd_ = TaylorInterpolant(fwd.t0, fwd.t, collect(fwd.x))
        return new{D, T, U}(dynamics, tracklets, _bwd_, _fwd_, res, fit, J, qs, Qs)
    end
end

# Outer constructor
LeastSquaresOrbit(
    dynamics::D, tracklets::Vector{Tracklet{T}}, bwd::TaylorInterpolant{T, U, 2},
    fwd::TaylorInterpolant{T, U, 2}, res::Vector{OpticalResidual{T, U}},
    fit::LeastSquaresFit{T}, J::Matrix{T}, qs::Matrix{T}, Qs::Vector{T}
    ) where {D, T <: Real, U <: Number} = LeastSquaresOrbit{D, T, U}(dynamics,
    tracklets, bwd, fwd, res, fit, J, qs, Qs)

# Print method for LeastSquaresOrbit
show(io::IO, orbit::LeastSquaresOrbit) = print(io, "Least squares orbit with ",
    length(orbit.res), " residuals")

# Definition of zero LeastSquaresOrbit
function zero(::Type{LeastSquaresOrbit{D, T, U}}) where {D, T <: Real, U <: Number}
    dynamics = D.instance
    tracklets = Vector{Tracklet{T}}(undef, 0)
    bwd = zero(TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}})
    fwd = zero(TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}})
    res = Vector{OpticalResidual{T, U}}(undef, 0)
    fit = zero(LeastSquaresFit{T})
    J = Matrix{T}(undef, 0, 0)
    qs = Matrix{T}(undef, 0, 0)
    Qs = Vector{T}(undef, 0)
    return LeastSquaresOrbit{D, T, U}(dynamics, tracklets, bwd, fwd, res, fit, J, qs, Qs)
end

iszero(orbit::LeastSquaresOrbit{D, T, U}) where {D, T <: Real, U <: Number} =
    orbit == zero(LeastSquaresOrbit{D, T, U})

# Evaluate integrations and residuals in fit deltas
function evalfit(orbit::LeastSquaresOrbit{D, T, TaylorN{T}}) where {D, T <: Real}
    # Fit δs
    δs = orbit.fit.x
    # Evaluate integrations
    new_bwd_x = map(x -> Taylor1(x.coeffs(δs)), orbit.bwd.x)
    new_bwd = TaylorInterpolant(orbit.bwd.t0, orbit.bwd.t, new_bwd_x)
    new_fwd_x = map(x -> Taylor1(x.coeffs(δs)), orbit.fwd.x)
    new_fwd = TaylorInterpolant(orbit.fwd.t0, orbit.fwd.t, new_fwd_x)
    # Evaluate residuals
    new_res = orbit.res(δs)

    return LeastSquaresOrbit{D, T, T}(orbit.dynamics, orbit.tracklets, new_bwd,
        new_fwd, new_res, orbit.fit, orbit.J, orbit.qs, orbit.Qs)
end