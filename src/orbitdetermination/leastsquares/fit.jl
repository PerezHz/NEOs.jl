include("methods.jl")
include("targetfunctions.jl")
include("outlierrejection.jl")

@doc raw"""
    LeastSquaresFit{T} <: AbstractLeastSquares{T}

A least squares fit.

## Fields

- `success::Bool`: whether the routine converged or not.
- `x::Vector{T}`: deltas that minimize the objective function.
- `Γ::Matrix{T}`: covariance matrix.
- `routine::Type{<:AbstractLeastSquaresMethod}`: minimization routine.
"""
@auto_hash_equals struct LeastSquaresFit{T} <: AbstractLeastSquares{T}
    success::Bool
    x::Vector{T}
    Γ::Matrix{T}
    routine::Type{<:AbstractLeastSquaresMethod}
    # Inner constructor
    LeastSquaresFit(success::Bool, x::Vector{T}, Γ::Matrix{T},
        routine::Type{<:AbstractLeastSquaresMethod}) where {T <: Real} = new{T}(success, x, Γ, routine)
end

# Outer constructor
LeastSquaresFit(::Type{T}, method::Type{<:AbstractLeastSquaresMethod}) where {T <: Real} =
    LeastSquaresFit(false, Vector{T}(undef, 0), Matrix{T}(undef, 0, 0), method)

# Definition of zero LeastSquaresFit{T}
zero(::Type{LeastSquaresFit{T}}) where {T <: Real} =
    LeastSquaresFit(false, Vector{T}(undef, 0), Matrix{T}(undef, 0, 0),
    AbstractLeastSquaresMethod)

# Print method for LeastSquaresFit
# Examples:
# Succesful Newton
# Succesful Differential Corrections
function show(io::IO, fit::LeastSquaresFit)
    success_s = fit.success ? "Succesful" : "Unsuccesful"
    routine_s = String(fit.routine)
    print(io, success_s, " ", routine_s)
end

@doc raw"""
    project(y::Vector{TaylorN{T}}, fit::LeastSquaresFit{T}) where {T <: Real}

Project the covariance matrix of `fit` into `y`.
"""
function project(y::Vector{TaylorN{T}}, fit::LeastSquaresFit{T}) where {T <: Real}
    J = Matrix{T}(undef, get_numvars(), length(y))
    for i in eachindex(y)
        J[:, i] = TS.gradient(y[i])(fit.x)
    end
    return (J') * fit.Γ * J
end

# Target functions
chi2(res::AbstractVector{OpticalResidual{T, TaylorN{T}}},
    fit::LeastSquaresFit{T}) where {T <: Real} = chi2(res(fit.x))

nms(res::AbstractVector{OpticalResidual{T, TaylorN{T}}},
    fit::LeastSquaresFit{T}) where {T <: Real} = nms(res(fit.x))

nrms(res::AbstractVector{OpticalResidual{T, TaylorN{T}}},
    fit::LeastSquaresFit{T}) where {T <: Real} = nrms(res(fit.x))

# Chi-square critical value
function critical_value(res::AbstractVector{OpticalResidual{T, TaylorN{T}}},
    fit::LeastSquaresFit{T}) where {T <: Real}
    # Evaluate residuals
    _res_ = res(fit.x)
    # Chi square distribution with N degrees of freedom
    N = 2 * notout(_res_)
    χ2 = N * nms(_res_)
    d = Chisq(N)
    # Evaluate cumulative probability at χ2
    return cdf(d, χ2)
end

# NMS threshold
function nms_threshold(N::Int, significance::T) where {T <: Real}
    # Chi-square distribution
    dist = Chisq(N)
    # Find Q for which cdf(dist, Q) == significance
    return find_zeros(t -> cdf(dist, t) - significance, 0, 1e6)[1] / N
end