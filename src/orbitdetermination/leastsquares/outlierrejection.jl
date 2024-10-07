@doc raw"""
    carpino_smoothing(n::T) where {T <: Real}

Fudge term for rejection condition in [`outlier_rejection`](@ref).

!!! reference
    See page 253 of https://doi.org/10.1016/S0019-1035(03)00051-4.
"""
carpino_smoothing(n::T) where {T <: Real} = 400*(1.2)^(-n)

@doc raw"""
    outlier_rejection!(w1::AbstractVector{T}, w2::AbstractVector{T},
        res::AbstractVector{OpticalResidual{T, TaylorN{T}}}, x::Vector{T},
        Γ::Matrix{T}; kwargs...) where {T <: Real}

Carpino et al. (2003) outlier rejection algorithm.

## Arguments

- `w1::AbstractVector{T}`: new weights.
- `w2::AbstractVector{T}`: original weights.
- `res::AbstractVector{OpticalResidual{T, TaylorN{T}}}`: vector of residuals.
- `x::Vector{T}`: evaluation deltas.
- `Γ::Matrix{T}`: covariance matrix.

## Keyword arguments

- `χ2_rec::T`: recovery threshold (default: `7.0`).
- `χ2_rej::T`: rejection threshold (default: `8.0`).
- `α::T`: scaling factor for maximum chi (default: `0.25`).
- `max_per::T`: maximum allowed drop percentage (default: `10.0`).

!!! reference
    See https://doi.org/10.1016/S0019-1035(03)00051-4.
"""
function outlier_rejection!(w1::AbstractVector{T}, w2::AbstractVector{T},
    res::AbstractVector{OpticalResidual{T, TaylorN{T}}}, x::Vector{T},
    Γ::Matrix{T}; χ2_rec::T = 7.0, χ2_rej::T = 8.0, α::T = 0.25,
    max_per::T = 10.0) where {T <: Real}
    # Number of residuals
    L = length(res)
    # Evaluate residuals
    eval_res = res(x)
    # Vector of χ2
    χ2s = Vector{T}(undef, L)
    # Maximum χ2 (over non outliers)
    χ2_max = zero(T)
    # Number of non outliers
    N_sel = 0
    # Compute χ2s
    for i in eachindex(χ2s)
        # Outlier flag
        outlier = iszero(w1[i])
        # Weights of current residual
        w_α, w_δ = w2[i], w2[i]
        # Current observation covariance matrix
        γ = diagm([inv(w_α), inv(w_δ)])
        # Current model matrix
        A = hcat(TS.gradient(res[i].ξ_α)(x), TS.gradient(res[i].ξ_δ)(x))
        # Outlier sign
        outlier_sign = 2*outlier - 1
        # Current residual covariance matrix
        γ_ξ = γ + outlier_sign * (A') * Γ * A
        # Current residual
        ξ = [eval_res[i].ξ_α, eval_res[i].ξ_δ]
        # Current chi2
        χ2s[i] = ξ' * inv(γ_ξ) * ξ
        # Update N_sel
        if !outlier
            N_sel += 1
            # Update maximum χ2
            χ2_max = max(χ2_max, χ2s[i])
        end
    end
    # Maximum allowed drops
    max_drop = ceil(Int, max_per * L / 100)
    # Number of dropped residuals
    N_drop = 0
    # Sort χ2s
    idxs = sortperm(χ2s, rev = true)
    # Rejection threshold
    χ2_rej = max(χ2_rej + carpino_smoothing(N_sel), α*χ2_max)
    # Rejection / recovery loop
    for i in idxs
        # Reject
        if χ2s[i] > χ2_rej && N_drop < max_drop
            w1[i] = zero(T)
            N_drop += 1
        # Recover
        elseif χ2s[i] < χ2_rec
            w1[i] = w2[i]
        end
    end

    return nothing
end