"""
    LeastSquaresFit{T} <: AbstractLeastSquares{T}

A least squares fit.

# Fields

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
end

# Outer constructor
LeastSquaresFit(::Type{T}, method::Type{<:AbstractLeastSquaresMethod}) where {T <: Real} =
    LeastSquaresFit(false, Vector{T}(undef, 0), Matrix{T}(undef, 0, 0), method)

# Definition of zero LeastSquaresFit{T}
zero(::Type{LeastSquaresFit{T}}) where {T <: Real} = LeastSquaresFit(false,
    Vector{T}(undef, 0), Matrix{T}(undef, 0, 0), AbstractLeastSquaresMethod)

# Print method for LeastSquaresFit
function show(io::IO, x::LeastSquaresFit)
    success_s = x.success ? "Succesful" : "Unsuccesful"
    routine_s = string(x.routine)
    print(io, success_s, " ", routine_s)
end

project(::AbstractVector{T}, Γ::AbstractMatrix{T}) where {T <: Real} =
    fill(T(NaN), size(Γ))

function project(y::AbstractVector{TaylorN{T}}, Γ::AbstractMatrix{T}) where {T <: Real}
    J = Matrix{TaylorN{T}}(undef, get_numvars(), length(y))
    for i in eachindex(y)
        J[:, i] = TS.gradient(y[i])
    end
    return (J') * Γ * J
end

function project(y::AbstractVector{TaylorN{T}}, x::AbstractVector{T},
                 Γ::AbstractMatrix{T}) where {T <: Real}
    J = Matrix{T}(undef, get_numvars(), length(y))
    for i in eachindex(y)
        J[:, i] = TS.gradient(y[i])(x)
    end
    return (J') * Γ * J
end

"""
    project(y, fit)

Project the covariance matrix of a least squares `fit` into a vector `y`.
"""
project(y::Vector{TaylorN{T}}, fit::LeastSquaresFit{T}) where {T <: Real} =
    project(y, fit.x, fit.Γ)

# Target functions
chi2(res::AbstractResidualVector, fit::LeastSquaresFit) = chi2(res(fit.x))
chi2(res::Tuple{O, R}, fit::LeastSquaresFit) where {O, R} = chi2(res[1], fit) + chi2(res[2], fit)
nms(x::AbstractResidualSet, fit::LeastSquaresFit) = chi2(x, fit) / notoutobs(x)
nrms(x::AbstractResidualSet, fit::LeastSquaresFit) = sqrt(nms(x, fit))

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
function nms_threshold(N::Int, significance::Real)
    # Chi-square distribution
    dist = Chisq(N)
    # Find Q for which cdf(dist, Q) == significance
    return find_zeros(t -> cdf(dist, t) - significance, 0, 1e6)[1] / N
end