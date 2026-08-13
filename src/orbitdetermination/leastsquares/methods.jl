"""
    AbstractLeastSquares{T <: Real}

Supertype for the least squares interface.
"""
abstract type AbstractLeastSquares{T <: Real} end

numtype(::AbstractLeastSquares{T}) where {T} = T

"""
    AbstractLeastSquaresCache{T} <: AbstractLeastSquares{T}

Supertype for the least squares caches interface.
"""
abstract type AbstractLeastSquaresCache{T} <: AbstractLeastSquares{T} end

struct LeastSquaresCache{T} <: AbstractLeastSquaresCache{T}
    x0::Vector{T}
    xs::Matrix{T}
    Qs::Vector{T}
    idxs::Vector{Int}
    maxiter::Int
end

"""
    LeastSquaresCache(x0 [, idxs]; kwargs...)

Return a cache for [`leastsquares!`](@ref) with initial guess `x0`. If `idxs`
(default: `eachindex(x0)`) is given, perform the minimization over a subset
of the parameters.

# Keyword argument

- `maxiter::Int`: maximum number of iterations (default: `25`)
"""
function LeastSquaresCache(x0::Vector{T}, idxs::AbstractVector{Int} = eachindex(x0),
    maxiter::Int = 25) where {T <: Real}
    xs = Matrix{T}(undef, length(x0), maxiter + 1)
    Qs = Vector{T}(undef, maxiter + 1)
    return LeastSquaresCache{T}(x0, xs, Qs, idxs, maxiter)
end

"""
    AbstractLeastSquaresMethod{T} <: AbstractLeastSquares{T}

Supertype for the least squares methods interface.

Every least squares `method`:
- is a mutable subtype of `AbstractLeastSquaresMethod{T}`,
- has a `method(res, x0 [, idxs])` constructor,
- defines `getid`, `lsstep`, `normalmatrix` and `update`.
"""
abstract type AbstractLeastSquaresMethod{T} <: AbstractLeastSquares{T} end

"""
    leastsquares!(ls, cache)

Minimize the target function in `ls` via least squares, using the
containers in `cache` to save intermediate results.
"""
function leastsquares!(ls::AbstractLeastSquaresMethod{T},
                       cache::LeastSquaresCache{T}) where {T <: Real}
    # Allocate memory for least squares fit
    fit = LeastSquaresFit(T, typeof(ls))
    # Unfold
    @unpack x0, xs, Qs, idxs, maxiter = cache
    # Initial conditions
    for j in axes(xs, 2)
        xs[:, j] .= x0
    end
    Qs[1] = evaluate(ls.Q, x0)
    Qs[1] < 0 && return fit
    Qs[2:end] .= T(Inf)
    # Main loop
    for i in 1:maxiter
        # Least squares step
        Δx, flag = lsstep(ls, xs[:, i])
        !flag && return fit
        # Update rule
        xs[idxs, i+1] .= xs[idxs, i] + Δx
        Qs[i+1] = evaluate(ls.Q, xs[:, i+1])
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
        fit = LeastSquaresFit(true, xs[:, i], Γ, typeof(ls))
    end

    return fit
end

"""
    leastsquares(method, res, x0 [, idxs]; kwargs...)

Use least squares `method` to minimize the normalized mean square of a
set of O-C residuals `res`, starting from initial guess `x0`. If `idxs`
(default: `eachindex(x0)`) is given, perform the minimization over a
subset of the parameters.

# Keyword argument

- `maxiter::Int`: maximum number of iterations (default: `25`).
"""
function leastsquares(method::Type{<:AbstractLeastSquaresMethod{T}},
                      res::AbstractResidualSet{T, TaylorN{T}},
                      x0::Vector{T}, idxs::AbstractVector{Int} = eachindex(x0);
                      maxiter::Int = 25) where {T <: Real}
    # Consistency checks
    @assert length(idxs) ≤ length(x0) == get_numvars()
    # Initialize least squares method and cache
    ls = method(res, x0, idxs)
    cache = LeastSquaresCache(x0, idxs, maxiter)
    # Main loop
    fit = leastsquares!(ls, cache)

    return fit
end

"""
    Newton(res, x0 [, idxs])

Newton method for minimizing the normalized mean square of a set of O-C residuals
`res`, starting from initial guess `x0`. If `idxs` (default: `eachindex(x0)`) is
given, perform the minimization over a subset of the parameters.

See also [`leastsquares`](@ref).

!!! reference
    See sections 5.2 and 5.3 of:
    - https://doi.org/10.1017/CBO9781139175371
"""
mutable struct Newton{T} <: AbstractLeastSquaresMethod{T}
    nobs::Int
    npar::Int
    Q::TaylorN{T}
    GQ::Vector{TaylorN{T}}
    HQ::Matrix{TaylorN{T}}
    dQ::Vector{T}
    d2Q::Matrix{T}
end

getid(::Newton) = "Newton"
targetfunction(x::Newton) = x.Q

function Newton(res::AbstractResidualSet{T, TaylorN{T}}, x0::Vector{T},
                idxs::AbstractVector{Int} = eachindex(x0)) where {T <: Real}
    # Number of observations and degrees of freedom
    nobs, npar = notoutobs(res), length(idxs)
    # Mean square residual and its gradient
    Q = nms(res)
    GQ = TaylorSeries.gradient(Q)[idxs]
    # Hessian of the target function
    HQ = [zero(Q) for _ in 1:npar, _ in 1:npar]
    for j in eachindex(idxs)
        for i in eachindex(idxs)
            HQ[i, j] = TS.differentiate(GQ[idxs[j]], idxs[i])
        end
    end
    # Evaluate gradient and hessian
    dQ = evaluate(GQ, x0)
    d2Q = evaluate(HQ, x0)

    return Newton{T}(nobs, npar, Q, GQ, HQ, dQ, d2Q)
end

function lsstep(ls::Newton, x::AbstractVector)
    # Unpack
    @unpack GQ, dQ, HQ, d2Q, nobs, npar = ls
    # Evaluate gradient and hessian
    evaluate!(GQ, x, dQ)
    evaluate!(HQ, x, d2Q)
    # Invert hessian
    lud2Q = lu!(d2Q; check = false)
    !issuccess(lud2Q) && return empty(x), false
    invd2Q = inv!(lud2Q)
    # Newton update rule
    Δx = -invd2Q * dQ
    # Normal matrix
    C = (nobs/2) * d2Q
    # Error metric
    error = (Δx') * C * Δx / npar
    return Δx, error > 0
end

function normalmatrix(ls::Newton, x::AbstractVector)
    @unpack HQ, d2Q, nobs = ls
    evaluate!(HQ, x, d2Q)
    return (nobs / 2) * d2Q
end

function update!(ls::Newton{T}, res::AbstractResidualSet{T, TaylorN{T}}, x0::Vector{T},
                 idxs::AbstractVector{Int}) where {T <: Real}
    # Number of observations and degrees of freedom
    ls.nobs, ls.npar = notoutobs(res), length(idxs)
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
    evaluate!(ls.GQ, x0, ls.dQ)
    evaluate!(ls.HQ, x0, ls.d2Q)

    return nothing
end

"""
    DifferentialCorrections(res, x0 [, idxs])

Differential corrections method for minimizing the normalized mean square of a
set of O-C residuals `res`, starting from initial guess `x0`. If `idxs` (default:
`eachindex(x0)`) is given, perform the minimization over a subset of the parameters.

See also [`leastsquares`](@ref).

!!! reference
    See sections 5.2 and 5.3 of:
    - https://doi.org/10.1017/CBO9781139175371
"""
mutable struct DifferentialCorrections{T} <: AbstractLeastSquaresMethod{T}
    nobs::Int
    npar::Int
    Q::TaylorN{T}
    D::Vector{TaylorN{T}}
    C::Matrix{TaylorN{T}}
    Dx::Vector{T}
    Cx::Matrix{T}
end

getid(::DifferentialCorrections) = "Differential Corrections"
targetfunction(x::DifferentialCorrections) = x.Q

function DifferentialCorrections(res::AbstractResidualSet{T, TaylorN{T}}, x0::Vector{T},
                                 idxs::AbstractVector{Int} = eachindex(x0)) where {T <: Real}
    # Number of observations and degrees of freedom
    nobs, npar = notoutobs(res), length(idxs)
    # Mean square residual and its gradient
    Q = nms(res)
    # D matrix and normal matrix C
    D, C, _, _ = DCVB(res, idxs)
    # Deprecated term
    # C -> C + ξTH(res, V, B, idxs)
    # Evaluate D and C matrices
    Dx = evaluate(D, x0)
    Cx = evaluate(C, x0)

    return DifferentialCorrections{T}(nobs, npar, Q, D, C, Dx, Cx)
end

function lsstep(ls::DifferentialCorrections, x::AbstractVector)
    # Unpack
    @unpack D, Dx, C, Cx, npar = ls
    # Evaluate D and C matrices
    evaluate!(D, x, Dx)
    evaluate!(C, x, Cx)
    # Invert C matrix
    luCx = lu!(Cx; check = false)
    !issuccess(luCx) && return empty(x), false
    invCx = inv!(luCx)
    # Differential corrections update rule
    Δx = -invCx * Dx
    # Error metric
    error = (Δx') * Cx * Δx / npar
    return Δx, error > 0
end

function normalmatrix(ls::DifferentialCorrections, x::AbstractVector)
    @unpack C, Cx = ls
    evaluate!(C, x, Cx)
    return Cx
end

function update!(ls::DifferentialCorrections{T}, res::AbstractResidualSet{T, TaylorN{T}},
                 x0::Vector{T}, idxs::AbstractVector{Int}) where {T <: Real}
    # Number of observations and degrees of freedom
    ls.nobs, ls.npar = notoutobs(res), length(idxs)
    # Mean square residual and its gradient
    ls.Q = nms(res)
    # D matrix and normal matrix C
    ls.D, ls.C, _, _ = DCVB(res, idxs)
    # Deprecated term
    # C -> C + ξTH(res, V, B, idxs)
    # Evaluate D and C matrices
    evaluate!(ls.D, x0, ls.Dx)
    evaluate!(ls.C, x0, ls.Cx)

    return nothing
end

# D matrix and normal matrix C
# See sections 5.2 and 5.3 of https://doi.org/10.1017/CBO9781139175371
function DCVB(res::AbstractResidualSet{T, TaylorN{T}},
              idxs::AbstractVector{Int}) where {T <: Real}
    # Number of observations and degrees of freedom
    nobs, npar = notoutobs(res), length(idxs)
    # Vector of normalized residuals
    V = normalized_residuals(res)
    # Allocate memory
    D = Vector{TaylorN{T}}(undef, npar)
    B = Matrix{TaylorN{T}}(undef, nobs, npar)
    C = Matrix{TaylorN{T}}(undef, npar, npar)
    BT = Matrix{TaylorN{T}}(undef, npar, nobs)
    # Transpose of the normalized design matrix B
    for j in eachindex(V)
        # Gradient of the j-th normalized residual with respect to
        # to the initial conditions x_0[idxs]
        for i in eachindex(idxs)
            BT[i, j] = TS.differentiate(V[j], idxs[i])
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
function ξTH(res::AbstractResidualSet{T, TaylorN{T}}, V::AbstractVector{TaylorN{T}},
             B::AbstractMatrix{TaylorN{T}}, idxs::AbstractVector{Int}) where {T <: Real}
    # Number of observations and degrees of freedom
    nobs, npar = notoutobs(res), length(idxs)
    # Allocate memory
    ξTHv = Array{T}(undef, npar, npar)
    H = Array{TaylorN{T}}(undef, nobs, npar, npar)
    # H matrix
    for i in axes(B, 1)
        for j in axes(B, 2)
            # Gradient of the (i, j)-th element of B with respect to the initial
            # conditions x_0[idxs]
            for k in eachindex(idxs)
                H[i, j, k] = TS.differentiate(B[i, j], idxs[k])
            end
        end
    end
    # Deprecated term ξTH
    for j in 1:npar
        for i in 1:npar
            # transpose(ξ) * H matrix
            @views ξTHv[i, j] = dot(V, H[:, i, j])
        end
    end

    return ξTHv, H
end

"""
    LevenbergMarquardt(res, x0 [, idxs])

Levenberg-Marquardt method for minimizing the normalized mean square of a
set of O-C residuals `res`, starting from initial guess `x0`. If `idxs`
(default: `eachindex(x0)`) is given, perform the minimization over a subset
of the parameters.

See also [`leastsquares`](@ref).

!!! reference
    See section 15.5.2 of:
    - https://numerical.recipes
"""
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

getid(::LevenbergMarquardt) = "Levenberg-Marquardt"
targetfunction(x::LevenbergMarquardt) = x.Q

function LevenbergMarquardt(res::AbstractResidualSet{T, TaylorN{T}}, x0::Vector{T},
                            idxs::AbstractVector{Int} = eachindex(x0)) where {T <: Real}
    # Number of observations and degrees of freedom
    nobs, npar = notoutobs(res), length(idxs)
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
    dQ = evaluate(GQ, x0)
    d2Q = evaluate(HQ, x0)

    return LevenbergMarquardt{T}(nobs, npar, idxs, λ, Q, GQ, HQ, dQ, d2Q)
end

function lsstep(ls::LevenbergMarquardt, x::AbstractVector)
    # Unpack
    @unpack GQ, dQ, HQ, d2Q, λ, idxs, Q = ls
    # Evaluate gradient and hessian
    evaluate!(GQ, x, dQ)
    evaluate!(HQ, x, d2Q)
    # Modified Hessian
    for i in axes(d2Q, 1)
        d2Q[i, i] +=  λ * abs(d2Q[i, i])
    end
    # Invert hessian
    lud2Q = lu!(d2Q; check = false)
    !issuccess(lud2Q) && return empty(x), false
    invd2Q = inv!(lud2Q)
    # Levenberg-Marquardt update rule
    Δx = -invd2Q * dQ
    y = deepcopy(x)
    y[idxs] .+= Δx
    # Update λ
    if 0 < Q(y) < Q(x)
        ls.λ /= 10
    else
        ls.λ *= 10
        Δx .= zero(eltype(x))
    end
    return Δx, true
end

function normalmatrix(ls::LevenbergMarquardt, x::AbstractVector)
    @unpack HQ, d2Q, nobs = ls
    evaluate!(HQ, x, d2Q)
    return (nobs / 2) * d2Q
end

function update!(ls::LevenbergMarquardt{T}, res::AbstractResidualSet{T, TaylorN{T}},
                 x0::Vector{T}, idxs::AbstractVector{Int}) where {T <: Real}
    # Number of observations and degrees of freedom
    ls.nobs, ls.npar = notoutobs(res), length(idxs)
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
    evaluate!(ls.GQ, x0, ls.dQ)
    evaluate!(ls.HQ, x0, ls.d2Q)

    return nothing
end

# Default methods for tryls
function _lsmethods(
        res::AbstractResidualSet{T, TaylorN{T}}, x0::AbstractVector{T},
        idxs::AbstractVector{Int} = eachindex(x0)
    ) where {T <: Real}
    return (
        Newton(res, x0, idxs),
        DifferentialCorrections(res, x0, idxs),
        LevenbergMarquardt(res, x0, idxs)
    )
end

# Base case: we ran out of methods
_tryls(fit, res, x0, cache, methods::Tuple{}) = fit

# Recursive step
function _tryls(fit, res, x0, cache, methods::Tuple)
    method = first(methods)
    update!(method, res, x0, cache.idxs)
    newfit = leastsquares!(method, cache)
    if issuccess(newfit) && 0 < targetfunction(method, newfit) < targetfunction(method, fit)
        fit = newfit
    end
    return _tryls(fit, res, x0, cache, Base.tail(methods))
end

# Main entry point
function tryls(res::AbstractResidualSet, x0::AbstractVector,
               cache::LeastSquaresCache, methods::Tuple)
    fit = zero(LeastSquaresFit{eltype(x0)})
    return _tryls(fit, res, x0, cache, methods)
end

"""
    tryls(res, x0 [, idxs]; kwargs...)

Try the least squares minimization of the normalized mean square of a
set of O-C residuals `res`, starting from initial guess `x0`. If `idxs`
(default: `eachindex(x0)`) is given, perform the minimization over a
subset of the parameters.

See also [`LeastSquaresCache`](@ref), [`AbstractLeastSquaresMethod`](@ref),
[`_lsmethods`](@ref) and [`leastsquares!`](@ref).

# Keyword Argument

- `maxiter::Int`: maximum number of iterations (default: `25`).
"""
function tryls(
        res::AbstractResidualSet{T, TaylorN{T}}, x0::AbstractVector{T},
        idxs::AbstractVector{Int} = eachindex(x0); maxiter::Int = 25
    ) where {T <: Real}
    cache = LeastSquaresCache(x0, idxs, maxiter)
    methods = _lsmethods(res, x0, idxs)
    return tryls(res, x0, cache, methods)
end