@doc raw"""
    chi2(res::AbstractVector{OpticalResidual{T, U}} [,
        fit::LeastSquaresFit{T}]) where {T <: Real, U <: Number}

Return the chi square of `res`
```math
\chi^2 = \sum_{i=1}^m w_i \xi_i^2,
```
where ``\xi_i`` is the i-th optical residual with statistical weight ``w_i``.
If `fit` is given, evaluate the residuals in `fit.x`.
"""
chi2(res::AbstractVector{U}, w::AbstractVector{T}) where {T <: Real, U <: Number} =
    sum(w .* (res .^ 2))

chi2(x::OpticalResidual{T, U}) where {T <: Real, U <: Number} =
    x.w_α * x.ξ_α^2 + x.w_δ * x.ξ_δ^2

chi2(res::AbstractVector{OpticalResidual{T, U}}) where {T <: Real, U <: Number} =
    sum(chi2.(res))

@doc raw"""
    nms(res::AbstractVector{OpticalResidual{T, U}} [,
        fit::LeastSquaresFit{T}]) where {T <: Real, U <: Number}

Return the normalized mean square error of `res`. If `fit` is given,
evaluate the residuals in `fit.x`.

See [`chi2`](@ref).
"""
nms(res::AbstractVector{U}, w::AbstractVector{T}) where {T <: Real, U <: Number} =
    chi2(res, w) / length(res)

nms(res::AbstractVector{OpticalResidual{T, U}}) where {T <: Real, U <: Number} =
    chi2(res) / (2 * length(res))

@doc raw"""
    nrms(res::AbstractVector{OpticalResidual{T, U}} [,
        fit::LeastSquaresFit{T}]) where {T <: Real, U <: Number}

Return the normalized root mean square error of `res`. If `fit` is given,
evaluate the residuals in `fit.x`.

See also [`chi2`](@ref) and [`nms`](@ref).
"""
nrms(res::AbstractVector{U}, w::AbstractVector{T}) where {T <: Real, U <: Number} =
    sqrt(nms(res, w))

nrms(res::AbstractVector{OpticalResidual{T, U}}) where {T <: Real, U <: Number} =
    sqrt(nms(res))