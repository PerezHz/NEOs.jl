struct OutlierRejectionCache{T} <: AbstractLeastSquaresCache{T}
    mask::BitVector
    eval_res::Vector{OpticalResidual{T, T}}
    χ2s::Vector{T}
    ξ::MVector{2, T}
    γ::MMatrix{2, 2, T, 4}
    A::Matrix{T}
    γ_ξ::MMatrix{2, 2, T, 4}
    L::Int
end

# Constructor
function OutlierRejectionCache(::Type{T}, L::Int) where {T <: Real}
    mask = BitVector(undef, L)
    eval_res = Vector{OpticalResidual{T, T}}(undef, L)
    χ2s = Vector{T}(undef, L)
    ξ = MVector{2, T}(zero(T), zero(T))
    γ = MMatrix{2, 2, T}(I)
    A = Matrix{T}(undef, get_numvars(), 2)
    γ_ξ = MMatrix{2, 2, T}(I)
    return OutlierRejectionCache{T}(mask, eval_res, χ2s, ξ, γ, A, γ_ξ, L)
end

"""
    outlier_rejection!(res, x, Γ [, cache]; kwargs...)

Reject outliers in a vector of optical residuals `res` using the Carpino et al. (2003)
algorithm. The residuals are evaluated at `x` and have a covariance matrix `Γ`. A
pre-allocated `cache` can be passed to save memory.

See also [`carpino_smoothing`](@ref).

# Keyword arguments

- `χ2_rec::Real`: recovery threshold (default: `7.0`).
- `χ2_rej::Real`: rejection threshold (default: `8.0`).
- `fudge::Real`: rejection fudge term coefficient (default: `400.0`).
- `α::Real`: scaling factor for maximum chi (default: `0.25`).
- `max_per::Real`: maximum allowed drop percentage (default: `10.0`).

!!! reference
    See:
    - https://doi.org/10.1016/S0019-1035(03)00051-4
"""
function outlier_rejection!(res::AbstractVector{OpticalResidual{T, TaylorN{T}}},
    x::Vector{T}, Γ::Matrix{T}, cache::OutlierRejectionCache{T} =
    OutlierRejectionCache(T, length(res)); χ2_rec::T = 7.0, χ2_rej::T = 8.0,
    fudge::T = 400.0, α::T = 0.25, max_per::T = 10.0) where {T <: Real}
    # Number of residuals
    L = length(res)
    # Unfold
    mask = view(cache.mask, 1:L)
    eval_res = view(cache.eval_res, 1:L)
    χ2s = view(cache.χ2s, 1:L)
    @unpack ξ, γ, A, γ_ξ = cache
    # Outliers mask
    map!(isoutlier, mask, res)
    # Evaluate residuals
    map!(r -> r(x), eval_res, res)
    # Maximum χ2 (over non outliers)
    χ2_max = zero(T)
    # Compute χ2s
    @inbounds for i in eachindex(χ2s)
        # Skip degenerate observations
        if abs(corr(res[i])) ≥ 1
            χ2s[i] = T(Inf)
            continue
        end
        # Outlier flag
        outlier = mask[i]
        # Current observation covariance matrix
        # Note: since the residuals are already normalized, the single
        # observation covariance matrix is:
        γ[1], γ[2], γ[3], γ[4] = one(T), corr(res[i]), corr(res[i]), one(T)
        # Current model matrix
        A[:, 1] = TS.gradient(ra(res[i]))(x)
        A[:, 2] = TS.gradient(dec(res[i]))(x)
        # Outlier sign
        outlier_sign = 2 * outlier - 1
        # Current residual covariance matrix
        γ_ξ .= γ + outlier_sign * (A') * Γ * A
        # Current residual
        ξ .= ra(eval_res[i]), dec(eval_res[i])
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
    max_drop = round(Int, max_per * L / 100)
    # Sort χ2s
    idxs = sortperm(χ2s, rev = true)
    # Rejection threshold
    χ2_rej = max(χ2_rej + carpino_smoothing(N_sel, fudge), α * χ2_max)
    # Rejection / recovery loop
    @inbounds for i in idxs
        # Reject
        if χ2s[i] > χ2_rej && N_drop < max_drop && !mask[i]
            res[i] = OpticalResidual{T, TaylorN{T}}(ra(res[i]), dec(res[i]), wra(res[i]),
                wdec(res[i]), dra(res[i]), ddec(res[i]), corr(res[i]), true)
            N_drop += 1
        # Recover
        elseif χ2s[i] < χ2_rec && mask[i]
            res[i] = OpticalResidual{T, TaylorN{T}}(ra(res[i]), dec(res[i]), wra(res[i]),
                wdec(res[i]), dra(res[i]), ddec(res[i]), corr(res[i]), false)
            N_drop -= 1
        end
    end

    return nothing
end

"""
    carpino_smoothing(n::Int [, A::Real])

Carpino et al. (2003) rejection condition fudge term
with coefficient `A` (default: `400.0`).

See also [`outlier_rejection!`](@ref).

!!! reference
    See page 253 of:
    - https://doi.org/10.1016/S0019-1035(03)00051-4
"""
carpino_smoothing(n::Int, A::Real = 400.0) = A * (1.2)^(-n)