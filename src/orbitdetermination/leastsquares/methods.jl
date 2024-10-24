include("targetfunctions.jl")
include("fit.jl")

@doc raw"""
    AbstractLeastSquaresCache{T} <: AbstractLeastSquares{T}

Supertype for the least squares caches API.
"""
abstract type AbstractLeastSquaresCache{T} <: AbstractLeastSquares{T} end

struct LeastSquaresCache{T} <: AbstractLeastSquaresCache{T}
    x0::Vector{T}
    xs::Matrix{T}
    Qs::Vector{T}
    idxs::Vector{Int}
    maxiter::Int
    # Inner constructor
    function LeastSquaresCache(x0::Vector{T}, idxs::AbstractVector{Int},
        maxiter::Int = 25) where {T <: Real}
        xs = Matrix{T}(undef, length(x0), maxiter + 1)
        Qs = Vector{T}(undef, maxiter + 1)
        return new{T}(x0, xs, Qs, idxs, maxiter)
    end
end

include("outlierrejection.jl")

@doc raw"""
    AbstractLeastSquaresMethod{T} <: AbstractLeastSquares{T}

Supertype for the least squares methods API.
"""
abstract type AbstractLeastSquaresMethod{T} <: AbstractLeastSquares{T} end

# Least squares main loop
function leastsquares!(ls::LS, cache::LeastSquaresCache{T}) where {T <: Real,
    LS <: AbstractLeastSquaresMethod{T}}
    # Allocate memory for least squares fit
    fit = LeastSquaresFit(T, getid(ls))
    # Unfold
    x0, xs, Qs, idxs = cache.x0, cache.xs, cache.Qs, cache.idxs
    maxiter = cache.maxiter
    # Initial conditions
    for j in axes(xs, 2)
        xs[:, j] .= x0
    end
    Qs[1] = TS.evaluate(ls.Q, x0)
    Qs[1] < 0 && return fit
    Qs[2:end] .= T(Inf)
    # Main loop
    for i in 1:maxiter
        # Least squares step
        Δx, flag = lsstep(ls, xs[:, i])
        !flag && return fit
        # Update rule
        xs[idxs, i+1] .= xs[idxs, i] + Δx
        Qs[i+1] = TS.evaluate(ls.Q, xs[:, i+1])
        Qs[i+1] < 0 && return fit
        # TO DO: break if Q has not changed significantly in
        # 3-5 iterations
    end
    # Choose best iteration
    i = argmin(Qs)
    # Normal matrix
    C = normalmatrix(ls, xs[:, i])
    # Covariance matrix
    # TO DO: use pinv for badly conditioned normal matrices
    Γ = inv(C)
    # Update fit
    if all(diag(Γ) .> 0)
        fit = LeastSquaresFit(true, xs[:, i], Γ, Symbol(getid(ls)))
    end

    return fit
end

@doc raw"""
    leastsquares(method::LS, res::AbstractVector{OpticalResidual{T, TaylorN{T}}},
        x0::Vector{T} [, idxs::AbstractVector{Int}]; kwargs...) where {LS, T <: Real}

Use least squares method `LS` to minimize the normalized mean square residual of `res`,
starting from initial guess `x0`. If `idxs` (default: `eachindex(x0)`) is given,
perform the minimization over a subset of the parameters.

## Keyword arguments

- `maxiter::Int`: maximum number of iterations (default: `25`).
"""
function leastsquares(method::LS, res::AbstractVector{OpticalResidual{T, TaylorN{T}}},
    x0::Vector{T}, idxs::AbstractVector{Int} = eachindex(x0);
    maxiter::Int = 25) where {LS, T <: Real}
    # Consistency checks
    @assert length(idxs) ≤ length(x0) == get_numvars()
    # Initialize least squares method and cache
    ls = method(res, x0, idxs)
    cache = LeastSquaresCache(x0, idxs, maxiter)
    # Main loop
    fit = leastsquares!(ls, cache)

    return fit
end

mutable struct Newton{T} <: AbstractLeastSquaresMethod{T}
    nobs::Int
    npar::Int
    Q::TaylorN{T}
    GQ::Vector{TaylorN{T}}
    HQ::Matrix{TaylorN{T}}
    dQ::Vector{T}
    d2Q::Matrix{T}
end

@doc raw"""
    Newton(res::AbstractVector{OpticalResidual{T, TaylorN{T}}},
        x0::Vector{T}, idxs::AbstractVector{Int}) where {T <: Real}

Return a `Newton{T}` object for least squares minimization.

See also [`leastsquares`](@ref).

## Arguments

- `res::AbstractVector{OpticalResidual{T, TaylorN{T}}}`: vector of residuals.
- `x_0::Vector{T}`: initial guess.
- `idxs::AbstractVector{Int}`: subset of parameters for fit.

!!! reference
    See sections 5.2 and 5.3 of https://doi.org/10.1017/CBO9781139175371.
"""
function Newton(res::AbstractVector{OpticalResidual{T, TaylorN{T}}},
    x0::Vector{T}, idxs::AbstractVector{Int}) where {T <: Real}
    # Number of observations and degrees of freedom
    nobs, npar = 2*notout(res), length(idxs)
    # Mean square residual and its gradient
    Q = nms(res)
    GQ = TaylorSeries.gradient(Q)[idxs]
    # Hessian
    HQ = [zero(Q) for _ in 1:npar, _ in 1:npar]
    for j in eachindex(idxs)
        for i in eachindex(idxs)
            HQ[i, j] = TS.differentiate(GQ[idxs[j]], idxs[i])
        end
    end
    # Evaluate gradient and hessian
    dQ = TS.evaluate(GQ, x0)
    d2Q = TS.evaluate(HQ, x0)

    return Newton{T}(nobs, npar, Q, GQ, HQ, dQ, d2Q)
end

getid(::Newton) = "Newton"

function lsstep(ls::Newton{T}, x::Vector{T}) where {T <: Real}
    # Evaluate gradient and hessian
    TS.evaluate!(ls.GQ, x, ls.dQ)
    TS.evaluate!(ls.HQ, x, ls.d2Q)
    # Invert hessian
    lud2Q = lu!(ls.d2Q; check = false)
    !issuccess(lud2Q) && return Vector{T}(undef, 0), false
    invd2Q = inv!(lud2Q)
    # Newton update rule
    Δx = -invd2Q * ls.dQ
    # Normal matrix
    C = (ls.nobs/2) * ls.d2Q
    # Error metric
    error = (Δx') * C * Δx / ls.npar

    return Δx, error > 0
end

function normalmatrix(ls::Newton{T}, x::Vector{T}) where {T <: Real}
    TS.evaluate!(ls.HQ, x, ls.d2Q)
    return (ls.nobs/2) * ls.d2Q # C = d2Q / (2/m)
end

function update!(ls::Newton{T}, res::AbstractVector{OpticalResidual{T, TaylorN{T}}},
    x0::Vector{T}, idxs::AbstractVector{Int}) where {T <: Real}
    # Number of observations and degrees of freedom
    ls.nobs, ls.npar = 2*notout(res), length(idxs)
    # Mean square residual and its gradient
    ls.Q = nms(res)
    TS.zero!.(ls.GQ)
    ls.GQ .= TaylorSeries.gradient(ls.Q)[idxs]
    # Hessian
    TS.zero!.(ls.HQ)
    for j in eachindex(idxs)
        for i in eachindex(idxs)
            ls.HQ[i, j] = TS.differentiate(ls.GQ[idxs[j]], idxs[i])
        end
    end
    # Evaluate gradient and hessian
    TS.evaluate!(ls.GQ, x0, ls.dQ)
    TS.evaluate!(ls.HQ, x0, ls.d2Q)

    return nothing
end

mutable struct DifferentialCorrections{T} <: AbstractLeastSquaresMethod{T}
    nobs::Int
    npar::Int
    Q::TaylorN{T}
    D::Vector{TaylorN{T}}
    C::Matrix{TaylorN{T}}
    Dx::Vector{T}
    Cx::Matrix{T}
end

@doc raw"""
    DifferentialCorrections(res::AbstractVector{OpticalResidual{T, TaylorN{T}}},
        x0::Vector{T}, idxs::AbstractVector{Int}) where {T <: Real}

Return a `DifferentialCorrections{T}` object for least squares minimization.

See also [`leastsquares`](@ref).

## Arguments

- `res::AbstractVector{OpticalResidual{T, TaylorN{T}}}`: vector of residuals.
- `x_0::Vector{T}`: initial guess.
- `idxs::AbstractVector{Int}`: subset of parameters for fit.

!!! reference
    See sections 5.2 and 5.3 of https://doi.org/10.1017/CBO9781139175371.
"""
function DifferentialCorrections(res::AbstractVector{OpticalResidual{T, TaylorN{T}}},
    x0::Vector{T}, idxs::AbstractVector{Int}) where {T <: Real}
    # Number of observations and degrees of freedom
    nobs, npar = 2*notout(res), length(idxs)
    # Mean square residual and its gradient
    Q = nms(res)
    # D matrix and normal matrix C
    D, C = DCVB(res, idxs)
    # Deprecated term
    # C -> C + ξTH(res, idxs)
    # Evaluate D and C matrices
    Dx = TS.evaluate(D, x0)
    Cx = TS.evaluate(C, x0)

    return DifferentialCorrections{T}(nobs, npar, Q, D, C, Dx, Cx)
end

getid(::DifferentialCorrections) = "Differential Corrections"

function lsstep(ls::DifferentialCorrections{T}, x::Vector{T}) where {T <: Real}
    # Evaluate D and C matrices
    TS.evaluate!(ls.D, x, ls.Dx)
    TS.evaluate!(ls.C, x, ls.Cx)
    # Invert C matrix
    luCx = lu!(ls.Cx; check = false)
    !issuccess(luCx) && return Vector{T}(undef, 0), false
    invCx = inv!(luCx)
    # Differential corrections update rule
    Δx = -invCx * ls.Dx
    # Error metric
    error = (Δx') * ls.Cx * Δx / ls.npar

    return Δx, error > 0
end

function normalmatrix(ls::DifferentialCorrections{T}, x::Vector{T}) where {T <: Real}
    TS.evaluate!(ls.C, x, ls.Cx)
    return ls.Cx
end

function update!(ls::DifferentialCorrections{T},
    res::AbstractVector{OpticalResidual{T, TaylorN{T}}},
    x0::Vector{T}, idxs::AbstractVector{Int}) where {T <: Real}
    # Number of observations and degrees of freedom
    ls.nobs, ls.npar = 2*notout(res), length(idxs)
    # Mean square residual and its gradient
    ls.Q = nms(res)
    # D matrix and normal matrix C
    ls.D, ls.C = DCVB(res, idxs)
    # Deprecated term
    # C -> C + ξTH(res, idxs)
    # Evaluate D and C matrices
    TS.evaluate!(ls.D, x0, ls.Dx)
    TS.evaluate!(ls.C, x0, ls.Cx)

    return nothing
end

# D matrix and normal matrix C
# See sections 5.2 and 5.3 of https://doi.org/10.1017/CBO9781139175371
function DCVB(res::AbstractVector{OpticalResidual{T, TaylorN{T}}},
    idxs::AbstractVector{Int}) where {T <: Real}
    # Number of observations and degrees of freedom
    nobs, npar = 2*notout(res), length(idxs)
    # Allocate memory
    V = Vector{TaylorN{T}}(undef, nobs)
    D = Vector{TaylorN{T}}(undef, npar)
    B = Matrix{TaylorN{T}}(undef, nobs, npar)
    BT = Matrix{TaylorN{T}}(undef, npar, nobs)
    C = Matrix{TaylorN{T}}(undef, npar, npar)
    # Tranpose of the normalized design matrix B
    for (j, k) in enumerate(findall(!isoutlier, res))
        # Normalized j-th ra and dec residuals
        V[2j-1] = sqrt(res[k].w_α) * res[k].ξ_α
        V[2j] = sqrt(res[k].w_δ) * res[k].ξ_δ
        # Gradient of the j-th ra and dec residuals with respect to
        # to the initial conditions x_0[idxs]
        for i in eachindex(idxs)
            BT[i, 2j-1] = TS.differentiate(V[2j-1], idxs[i])
            BT[i, 2j] = TS.differentiate(V[2j], idxs[i])
        end
    end
    # D matrix: transpose(B) * W * ξ
    mul!(D, BT, V)
    # Normal matrix C
    transpose!(B, BT)
    mul!(C, BT, B)

    return D, C, V, B
end

# Deprecated term in differential corrections
# This function is not used, but it is included here for completeness
# See sections 5.2 and 5.3 of https://doi.org/10.1017/CBO9781139175371
function ξTH(res::AbstractVector{OpticalResidual{T, TaylorN{T}}},
    B::AbstractMatrix{TaylorN{T}}, idxs::AbstractVector{Int}) where {T <: Real}
    # Number of observations and degrees of freedom
    nobs, npar = 2*notout(res), length(idxs)
    # Allocate memory
    ξTHv = Array{T}(undef, npar, npar)
    A = Matrix{T}(undef, 1, nobs)
    H = Array{TaylorN{T}}(undef, nobs, npar, npar)
    # H matrix and auxiliary array
    for (i, k) in enumerate(findall(!isoutlier, res))
        A[1, 2i-1] = res[k].w_α * res[k].ξ_α
        A[1, 2i] = res[k].w_δ * res[k].ξ_δ
        for j in 1:npar
            # Gradient of the (i, j)-th element of B with respect to the initial
            # conditions x_0
            H[i, j, :] .= TS.gradient(B[i, j])[idxs]
        end
    end
    # Deprecated term ξTH
    for j in 1:npar
        for i in 1:npar
            # transpose(ξ) * W * H matrix
            ξTHv[i, j] = A * H[:, i, j]
        end
    end

    return ξTHv, H
end

mutable struct LevenbergMarquardt{T} <: AbstractLeastSquaresMethod{T}
    nobs::Int
    npar::Int
    idxs::Vector{Int}
    λ::T
    Q::TaylorN{T}
    GQ::Vector{TaylorN{T}}
    HQ::Matrix{TaylorN{T}}
    dQ::Vector{T}
    d2Q::Matrix{T}
end

@doc raw"""
    LevenbergMarquardt(res::AbstractVector{OpticalResidual{T, TaylorN{T}}},
        x0::Vector{T}, idxs::AbstractVector{Int}) where {T <: Real}

Return a `LevenbergMarquardt{T}` object for least squares minimization.

See also [`leastsquares`](@ref).

## Arguments

- `res::AbstractVector{OpticalResidual{T, TaylorN{T}}}`: vector of residuals.
- `x_0::Vector{T}`: initial guess.
- `idxs::AbstractVector{Int}`: subset of parameters for fit.

!!! reference
    See section 15.5.2 of https://numerical.recipes.
"""
function LevenbergMarquardt(res::AbstractVector{OpticalResidual{T, TaylorN{T}}},
    x0::Vector{T}, idxs::AbstractVector{Int}) where {T <: Real}
    # Number of observations and degrees of freedom
    nobs, npar = 2*notout(res), length(idxs)
    # Damping factor
    λ = T(0.001)
    # Mean square residual and its gradient
    Q = nms(res)
    GQ = TaylorSeries.gradient(Q)[idxs]
    # Hessian
    HQ = [zero(Q) for _ in 1:npar, _ in 1:npar]
    for j in eachindex(idxs)
        for i in eachindex(idxs)
            HQ[i, j] = TS.differentiate(GQ[idxs[j]], idxs[i])
        end
    end
    # Evaluate gradient and hessian
    dQ = TS.evaluate(GQ, x0)
    d2Q = TS.evaluate(HQ, x0)

    return LevenbergMarquardt{T}(nobs, npar, idxs, λ, Q, GQ, HQ, dQ, d2Q)
end

getid(::LevenbergMarquardt) = "Levenberg-Marquardt"

function lsstep(ls::LevenbergMarquardt{T}, x::Vector{T}) where {T <: Real}
    # Evaluate gradient and hessian
    TS.evaluate!(ls.GQ, x, ls.dQ)
    TS.evaluate!(ls.HQ, x, ls.d2Q)
    # Modified Hessian
    for i in axes(ls.d2Q, 1)
        ls.d2Q[i, i] *= (1 + ls.λ)
    end
    # Invert hessian
    lud2Q = lu!(ls.d2Q; check = false)
    !issuccess(lud2Q) && return Vector{T}(undef, 0), false
    invd2Q = inv!(lud2Q)
    # Levenberg-Marquardt update rule
    Δx = -invd2Q * ls.dQ
    _x_ = deepcopy(x)
    _x_[ls.idxs] .+= Δx
    # Update λ
    if 0 < ls.Q(_x_) < ls.Q(x)
        ls.λ /= 10
    else
        ls.λ *= 10
        Δx .= zero(T)
    end

    return Δx, true
end

function normalmatrix(ls::LevenbergMarquardt{T}, x::Vector{T}) where {T <: Real}
    TS.evaluate!(ls.HQ, x, ls.d2Q)
    return (ls.nobs/2) * ls.d2Q # C = d2Q / (2/m)
end

function update!(ls::LevenbergMarquardt{T},
    res::AbstractVector{OpticalResidual{T, TaylorN{T}}},
    x0::Vector{T}, idxs::AbstractVector{Int}) where {T <: Real}
    # Number of observations and degrees of freedom
    ls.nobs, ls.npar = 2*notout(res), length(idxs)
    # Damping factor
    ls.λ = T(0.001)
    # Mean square residual and its gradient
    ls.Q = nms(res)
    TS.zero!.(ls.GQ)
    ls.GQ .= TaylorSeries.gradient(ls.Q)[idxs]
    # Hessian
    TS.zero!.(ls.HQ)
    for j in eachindex(idxs)
        for i in eachindex(idxs)
            ls.HQ[i, j] = TS.differentiate(ls.GQ[idxs[j]], idxs[i])
        end
    end
    # Evaluate gradient and hessian
    TS.evaluate!(ls.GQ, x0, ls.dQ)
    TS.evaluate!(ls.HQ, x0, ls.d2Q)

    return nothing
end

_lsmethods(res::AbstractVector{OpticalResidual{T, TaylorN{T}}},
    x0::Vector{T}, idxs::AbstractVector{Int} = eachindex(x0)) where {T <: Real} =
    [Newton(res, x0, idxs), DifferentialCorrections(res, x0, idxs),
    LevenbergMarquardt(res, x0, idxs)]

@doc raw"""
    tryls(res::AbstractVector{OpticalResidual{T, TaylorN{T}}},
        x0::Vector{T} [, idxs::AbstractVector{Int}]; kwargs...) where {T <: Real}

Try the least squares minimization of the normalized mean square residual of `res`,
starting from initial guess `x0`. If `idxs` (default: `eachindex(x0)`) is given,
perform the minimization over a subset of the parameters.

See also [`LeastSquaresCache`](@ref), [`AbstractLeastSquaresMethod`](@ref),
[`_lsmethods`](@ref) and [`leastsquares!`](@ref).

## Keyword Arguments

- `maxiter::Int`: maximum number of iterations (default: `25`).
- `cache::LeastSquaresCache{T}`: least squares cache (default:
    `LeastSquaresCache(x0, idxs, maxiter)`).
- `methods::Vector{AbstractLeastSquaresMethod}`: least squares routines to try
    (default: `_lsmethods(res, x0, idxs)`).
"""
function tryls(res::AbstractVector{OpticalResidual{T, TaylorN{T}}},
    x0::Vector{T}, cache::LeastSquaresCache{T},
    methods::Vector{AbstractLeastSquaresMethod{T}}) where {T <: Real}
    # Allocate memory
    fit = zero(LeastSquaresFit{T})
    # Least squares methods in order
    for i in eachindex(methods)
        update!(methods[i], res, x0, cache.idxs)
        fit = leastsquares!(methods[i], cache)
        fit.success && break
    end

    return fit
end

function tryls(res::AbstractVector{OpticalResidual{T, TaylorN{T}}},
    x0::Vector{T}, idxs::AbstractVector{Int} = eachindex(x0);
    maxiter::Int = 25) where {T <: Real}
    cache = LeastSquaresCache(x0, idxs, maxiter)
    methods = _lsmethods(res, x0, idxs)
    return tryls(res, x0, cache, methods)
end