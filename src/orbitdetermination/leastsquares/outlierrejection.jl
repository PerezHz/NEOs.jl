struct OutlierRejectionCache{T} <: AbstractLeastSquaresCache{T}
    mask::BitVector
    eval_res::Vector{OpticalResidual{T, T}}
    χ2s::Vector{T}
    ξ::Vector{T}
    γ::Matrix{T}
    A::Matrix{T}
    γ_ξ::Matrix{T}
    L::Int
    # Inner constructor
    function OutlierRejectionCache(::Type{T}, L::Int) where {T <: Real}
        mask = BitVector(undef, L)
        eval_res = Vector{OpticalResidual{T, T}}(undef, L)
        χ2s = Vector{T}(undef, L)
        ξ = Vector{T}(undef, 2)
        γ = zeros(T, 2, 2)
        A = Matrix{T}(undef, 6, 2)
        γ_ξ = Matrix{T}(undef, 2, 2)
        return new{T}(mask, eval_res, χ2s, ξ, γ, A, γ_ξ, L)
    end
end

@doc raw"""
    outlier_rejection!(res::AbstractVector{OpticalResidual{T, TaylorN{T}}},
        x::Vector{T}, Γ::Matrix{T} [, cache::OutlierRejectionCache{T}];
        kwargs...) where {T <: Real}

Carpino et al. (2003) outlier rejection algorithm.

See also [`carpino_smoothing`](@ref).

## Arguments

- `res::AbstractVector{OpticalResidual{T, TaylorN{T}}}`: vector of residuals.
- `x::Vector{T}`: evaluation deltas.
- `Γ::Matrix{T}`: covariance matrix.
- `cache::OutlierRejectionCache{T}`: see [`OutlierRejectionCache`](@ref).

## Keyword arguments

- `χ2_rec::T`: recovery threshold (default: `7.0`).
- `χ2_rej::T`: rejection threshold (default: `8.0`).
- `fudge::T`: rejection fudge term coefficient (default: `400.0`).
- `α::T`: scaling factor for maximum chi (default: `0.25`).
- `max_per::T`: maximum allowed drop percentage (default: `10.0`).

!!! reference
    See https://doi.org/10.1016/S0019-1035(03)00051-4.
"""
function outlier_rejection!(res::AbstractVector{OpticalResidual{T, TaylorN{T}}},
    x::Vector{T}, Γ::Matrix{T}, cache::OutlierRejectionCache{T} =
    OutlierRejectionCache(T, length(res)); χ2_rec::T = 7.0, χ2_rej::T = 8.0,
    fudge::T = 400.0, α::T = 0.25, max_per::T = 10.0) where {T <: Real}
    # Number of residuals
    L = length(res)
    # Unfold
    mask, eval_res, χ2s, ξ, γ, A, γ_ξ = view(cache.mask, 1:L), view(cache.eval_res, 1:L),
        view(cache.χ2s, 1:L), cache.ξ, cache.γ, cache.A, cache.γ_ξ
    # Outliers mask
    map!(isoutlier, mask, res)
    # Evaluate residuals
    map!(r -> r(x), eval_res, res)
    # Maximum χ2 (over non outliers)
    χ2_max = zero(T)
    # Compute χ2s
    @inbounds for i in eachindex(χ2s)
        # Outlier flag
        outlier = mask[i]
        # Current observation covariance matrix
        γ[1, 1], γ[2, 2] = inv(res[i].w_α), inv(res[i].w_δ)
        # Current model matrix
        A[:, 1] = TS.gradient(res[i].ξ_α)(x)
        A[:, 2] = TS.gradient(res[i].ξ_δ)(x)
        # Outlier sign
        outlier_sign = 2 * outlier - 1
        # Current residual covariance matrix
        γ_ξ .= γ + outlier_sign * (A') * Γ * A
        # Current residual
        ξ .= eval_res[i].ξ_α, eval_res[i].ξ_δ
        # Current chi2
        χ2s[i] = abs(ξ' * inv(γ_ξ) * ξ)
        # Update maximum χ2
        if !outlier
            χ2_max = max(χ2_max, χ2s[i])
        end
    end
    # Number of dropped (outliers) and selected residuals
    N_drop = count(mask)
    N_sel = L - N_drop
    # Maximum allowed drops
    max_drop = ceil(Int, max_per * L / 100)
    # Sort χ2s
    idxs = sortperm(χ2s, rev = true)
    # Rejection threshold
    χ2_rej = max(χ2_rej + carpino_smoothing(N_sel, fudge), α*χ2_max)
    @show χ2s, χ2_rej
    # Rejection / recovery loop
    @inbounds for i in idxs
        # Reject
        if χ2s[i] > χ2_rej && N_drop < max_drop && !mask[i]
            res[i] = OpticalResidual(res[i].ξ_α, res[i].ξ_δ, res[i].w_α,
                res[i].w_δ, true)
            N_drop += 1
        # Recover
        elseif χ2s[i] < χ2_rec && mask[i]
            res[i] = OpticalResidual(res[i].ξ_α, res[i].ξ_δ, res[i].w_α,
                res[i].w_δ, false)
            N_drop -= 1
        end
    end

    return nothing
end

@doc raw"""
    carpino_smoothing(n::Int [, A::T]) where {T <: Real}

Carpino et al. (2003) rejection condition fudge term
with coefficient `A` (default: `400.0`).

See also [`outlier_rejection`](@ref).

!!! reference
    See page 253 of https://doi.org/10.1016/S0019-1035(03)00051-4.
"""
carpino_smoothing(n::Int, A::T = 400.0) where {T <: Real} = A*(1.2)^(-n)