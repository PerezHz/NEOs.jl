@doc raw"""
    LeastSquaresFit{T <: Real}

A least squares fit.

# Fields
- `success::Bool`: whether the routine converged or not.
- `x::Vector{T}`: deltas that minimize the objective function.
- `Γ::Matrix{T}`: covariance matrix.
- `routine::Symbol`: minimization routine (`:newton` or `:diffcorr`).
"""
@auto_hash_equals struct LeastSquaresFit{T <: Real}
    success::Bool
    x::Vector{T}
    Γ::Matrix{T}
    routine::Symbol
    # Inner constructor
    function LeastSquaresFit{T}(success::Bool, x::Vector{T}, Γ::Matrix{T},
                         routine::Symbol) where {T <: Real}
        return new{T}(success, x, Γ, routine)
    end
end
# Outer constructor
function LeastSquaresFit(success::Bool, x::Vector{T}, Γ::Matrix{T},
                  routine::Symbol) where {T <: Real}
    return LeastSquaresFit{T}(success, x, Γ, routine)
end

# Definition of zero LeastSquaresFit{T}
function zero(::Type{LeastSquaresFit{T}}) where {T <: Real}
    return LeastSquaresFit{T}(false, Vector{T}(undef, 0), Matrix{T}(undef, 0, 0), :zero)
end

# Print method for LeastSquaresFit
# Examples:
# Succesful Newton
# Succesful differential corrections
function show(io::IO, fit::LeastSquaresFit{T}) where {T <: Real}
    success_s = fit.success ? "Succesful" : "Unsuccesful"
    if fit.routine == :newton
        routine_s = "Newton"
    elseif fit.routine == :diffcorr
        routine_s = "Differential Corrections"
    elseif fit.routine == :lm
        routine_s = "Levenberg-Marquardt"
    else
        routine_s = "Least Squares"
    end
    print(io, success_s, " ", routine_s)
end

@doc raw"""
    carpino_smoothing(n::T) where {T <: Real}

Fudge term for rejection condition in [`outlier_rejection`](@ref).

!!! reference
    See page 253 of https://doi.org/10.1016/S0019-1035(03)00051-4.
"""
carpino_smoothing(n::T) where {T <: Real} = 400*(1.2)^(-n)

@doc raw"""
    outlier_rejection(ξs::Vector{OpticalResidual{T, TaylorN{T}}}, fit::LeastSquaresFit{T};
                      χ2_rec::T = 7.0, χ2_rej::T = 8.0, α::T = 0.25, max_per::T = 10.) where {T <: Real}

Outlier rejection algorithm.

# Arguments

- `ξs::Vector{OpticalResidual{T}}`: vector of residuals.
- `fit::LeastSquaresFit{T}`: least squares fit.
- `χ2_rec::T`: recovery condition.
- `χ2_rej::T`: rejection condition.
- `α::T`: scaling factor for maximum chi.
- `max_per::T`: maximum allowed drop percentage.

!!! reference
    See https://doi.org/10.1016/S0019-1035(03)00051-4.
"""
function outlier_rejection(res::Vector{OpticalResidual{T, TaylorN{T}}}, fit::LeastSquaresFit{T};
                           χ2_rec::T = 7., χ2_rej::T = 8., α::T = 0.25, max_per::T = 10.) where {T <: Real}

    # Number of residuals
    L = length(res)
    # Evaluate residuals
    eval_res = res(fit.x)
    # Vector of χ2
    χ2s = Vector{T}(undef, L)
    # Maximum χ2 (over non outliers)
    χ2_max = zero(T)
    # Number of non outliers
    N_sel = 0

    # Compute χ2s
    for i in eachindex(χ2s)
        # Weights of current residual
        w_α, w_δ = res[i].w_α / res[i].relax_factor, res[i].w_δ / res[i].relax_factor
        # Current observation covariance matrix
        γ = [w_α zero(T); zero(T) w_δ]
        # Current model matrix
        A = hcat(TS.gradient(res[i].ξ_α)(fit.x), TS.gradient(res[i].ξ_δ)(fit.x))
        # Outlier sign
        outlier_sign = res[i].outlier*2-1
        # Current residual covariance matrix
        γ_ξ = γ + outlier_sign*(A')*fit.Γ*A
        # Current residual
        ξ = [eval_res[i].ξ_α, eval_res[i].ξ_δ]
        # Current chi2
        χ2s[i] = ξ' * inv(γ_ξ) * ξ
        # Update N_sel
        if !res[i].outlier
            N_sel += 1
            # Update maximum χ2
            if χ2s[i] > χ2_max
                χ2_max = χ2s[i]
            end
        end
    end

    # Maximum allowed drops
    max_drop = ceil(Int, max_per * L / 100)
    # Number of dropped residuals
    N_drop = 0
    # New outliers
    new_outliers = outlier.(res)
    # Sort χ2s
    idxs = sortperm(χ2s, rev = true)
    # Rejection threshold
    χ2_rej = max(χ2_rej + carpino_smoothing(N_sel), α*χ2_max)

    for i in idxs
        if χ2s[i] >  χ2_rej && N_drop < max_drop
            new_outliers[i] = true
            N_drop += 1
        elseif χ2s[i] < χ2_rec
            new_outliers[i] = false
        end
    end

    return OpticalResidual.(ra.(res), dec.(res), weight_ra.(res), weight_dec.(res),
                            relax_factor.(res), new_outliers)
end

@doc raw"""
    project(y::Vector{TaylorN{T}}, fit::LeastSquaresFit{T}) where {T <: Real}

Project `fit`'s covariance matrix into `y`.
"""
function project(y::Vector{TaylorN{T}}, fit::LeastSquaresFit{T}) where {T <: Real}
    J = Matrix{T}(undef, get_numvars(), length(y))
    for i in eachindex(y)
        J[:, i] = TS.gradient(y[i])(fit.x)
    end
    return (J') * fit.Γ * J
end

@doc raw"""
    chi2(res::Vector{U}, w::Vector{T}) where {T <: Real, U <: Number}
    chi2(res::Vector{OpticalResidual{T, U}}) where {T <: Real, U <: Number}

Returns the chi square
```math
\chi^2 = \sum_{i=1}^m \frac{ \xi_i^2}{\sigma_i^2},
```
where ``\mathbf{w} = (1/\sigma_1^2,\ldots,1/\sigma_m^2)^T`` and ``\mathbf{\xi} = (\xi_1,\ldots,\xi_m)^T``
are the vectors of weights and residuals respectively.

# Arguments

- `res::Vector{U}/Vector{OpticalResidual{T, U}}`: vector of residuals.
- `w::Vector{T}`: vector of weights.
"""
function chi2(res::Vector{U}, w::Vector{T}) where {T <: Real, U <: Number}
    # Have as many residuals as weights
    @assert length(res) == length(w)
    # Chi square
    return sum(w .* (res.^2))
end
function chi2(res::Vector{OpticalResidual{T, U}}) where {T <: Real, U <: Number}
    _res_, _w_ = unfold(res)
    return chi2(_res_, _w_)
end

@doc raw"""
    nms(res::Vector{U}, w::Vector{T}) where {T <: Real, U <: Number}
    nms(res::Vector{OpticalResidual{T, U}}) where {T <: Real, U <: Number}
    nms(res::Vector{OpticalResidual{T, TaylorN{T}}}, fit::LeastSquaresFit{T}) where {T <: Real}

Return the normalized chi square. See [`chi2`](@ref).
"""
nms(res::Vector{U}, w::Vector{T}) where {T <: Real, U <: Number} = chi2(res, w) / length(res)

function nms(res::Vector{OpticalResidual{T, U}}) where {T <: Real, U <: Number}
    _res_, _w_ = unfold(res)
    return nms(_res_, _w_)
end

function nms(res::Vector{OpticalResidual{T, TaylorN{T}}}, fit::LeastSquaresFit{T}) where {T <: Real}
    _res_, _w_ = unfold(res)
    return nms(_res_(fit.x), _w_)
end

@doc raw"""
    nrms(res::Vector{U}, w::Vector{T}) where {T <: Real, U <: Number}
    nrms(res::Vector{OpticalResidual{T, U}}) where {T <: Real, U <: Number}
    nrms(res::Vector{OpticalResidual{T, TaylorN{T}}}, fit::LeastSquaresFit{T}) where {T <: Real}

Returns the normalized root mean square error
```math
\texttt{NRMS} = \sqrt{\frac{\chi^2}{m}},
```
where ``\chi^2`` is the chi square and ``\mathbf{\xi} = (\xi_1,\ldots,\xi_m)^T`` is the vector
of residuals.

See also [`chi2`](@ref).

# Arguments

- `res::Vector{U}/Vector{OpticalResidual{T, U}}`: Vector of residuals.
- `w::Vector{T}`: Vector of weights.
- `fit::LeastSquaresFit{T}`: least squares fit.

"""
function nrms(res::Vector{U}, w::Vector{T}) where {T <: Real, U <: Number}
    # Have as many residuals as weights
    @assert length(res) == length(w)
    # Normalized root mean square error
    return sqrt( chi2(res, w)/length(res) )
end

function nrms(res::Vector{OpticalResidual{T, U}}) where {T <: Real, U <: Number}
    _res_, _w_ = unfold(res)
    return nrms(_res_, _w_)
end

function nrms(res::Vector{OpticalResidual{T, TaylorN{T}}}, fit::LeastSquaresFit{T}) where {T <: Real}
    _res_, _w_ = unfold(res)
    return nrms(_res_(fit.x), _w_)
end

@doc raw"""
    BHC(res::Vector{TaylorN{T}}, w::Vector{T}, npar::Int) where {T <: Real}

Returns the ``\mathbf{B}``, ``\mathbf{H}`` and ``\mathbf{C}`` arrays
```math
\mathbf{B} = \frac{\partial\mathbf{\xi}}{\partial\mathbf{x}_0}(\mathbf{x}_0), \quad
\mathbf{H} = \frac{\partial^2\mathbf{\xi}}{\partial\mathbf{x}_0^2}(\mathbf{x}_0) \quad \text{and} \quad
\mathbf{C} = \mathbf{B}^T\mathbf{W}\mathbf{B},
```
where ``\mathbf{x}_0 = (x_1,\ldots,x_n)^T`` and ``\mathbf{\xi} = (\xi_1,\ldots,\xi_m)^T``
are the vectors of initial conditions and residuals respectively; and
``\mathbf{W} = \text{diag}(1/\sigma_1^2,\ldots,1/\sigma_m^2)`` is the weights matrix.

``\mathbf{B}`` is called the design matrix and is of size ``m\times n``, ``\mathbf{H}`` is a three index
array of size ``m\times n\times n`` and ``\mathbf{C}`` is called the normal matrix and is of size ``n\times n``.

# Arguments

- `res::Vector{TaylorN{T}}`: vector of residuals.
- `w::Vector{T}`: vector of weights.
- `npar::Int`: degrees of freedom ``n``.

!!! reference
    See sections 5.2 and 5.3 of https://doi.org/10.1017/CBO9781139175371.
"""
function BHC(res::Vector{TaylorN{T}}, w::Vector{T}, npar::Int) where {T <: Real}
    # Number of observations
    nobs = length(res)

    # Allocate memory for the three arrays
    B_mat = Matrix{TaylorN{T}}(undef, nobs, npar)
    H_mat = Array{TaylorN{T}}(undef, nobs, npar, npar)
    C_mat = Array{TaylorN{T}}(undef, npar, npar)

    # Design matrix B
    for i in 1:nobs
        # Gradient of the i-th residual with respect to the initial conditions x_0
        B_mat[i,:] .= TaylorSeries.gradient(res[i])
    end
    # H matrix
    for i in 1:nobs
        for j in 1:npar
                # Gradient of the (i, j)-th element of B with respect to to the initial
                # conditions x_0
                H_mat[i,j,:] .= TaylorSeries.gradient(B_mat[i,j])
        end
    end
    # Normal matrix C
    sqrtw_B = sqrt.(w) .* B_mat
    C_mat .= (sqrtw_B') * sqrtw_B

    return B_mat, H_mat, C_mat
end

@doc raw"""
    ξTH(w, res, H_mat, npar)

Returns ``\mathbf{\xi}^T\mathbf{W}\mathbf{H}``, where ``\mathbf{\xi} = (\xi_1,\ldots,\xi_m)^T``
is the vector of residuals, ``\mathbf{W} = \text{diag}(1/\sigma_1^2,\ldots,1/\sigma_m^2)`` is
the weights matrix and ``\mathbf{H}`` is the matrix of second derivatives of ``\mathbf{\xi}``
with respect to the initial conditions.

See also [`BHC`](@ref).

# Arguments

- `w`: Vector of weights.
- `res`: Vector or residuals.
- `H_mat`: matrix of second derivatives of ``\mathbf{\xi}`` with respect to the initial conditions.
- `npar`: Degrees of freedom ``n``.
"""
function ξTH(w, res, H_mat, npar)
    # Allocate memory for output
    ξTHv = Array{Float64}(undef, npar, npar)

    for j in 1:npar
        for i in 1:npar
           # transpose(ξ) * W * H matrix
           ξTHv[i,j] = (w .* res)' * (H_mat[:,i,j])
        end
    end

    return ξTHv
end

@doc raw"""
    diffcorr(res::Vector{TaylorN{T}}, w::Vector{T}, x0::Vector{T},
             niters::Int = 5) where {T <: Real}
    diffcorr(res::Vector{OpticalResidual{T, TaylorN{T}}}, x0::Vector{T},
             niters::Int = 5) where {T <: Real}

Differential corrections subroutine for least-squares fitting. Returns an
`LeastSquaresFit` with the `niters`-th
correction
```math
\mathbf{x}_{k+1} = \mathbf{x}_k - \mathbf{C}^{-1}\mathbf{D},
```
and the covariance matrix
```math
\mathbf{\Gamma} = \mathbf{C}^{-1},
```
where ``\mathbf{C} = \mathbf{B}^T\mathbf{W}\mathbf{B}`` is the normal matrix and
``\mathbf{D} = \mathbf{B}^T\mathbf{W}\mathbf{\xi}``, with ``\mathbf{\xi} = (\xi_1,\ldots,\xi_m)^T``
the vector of residuals, ``\mathbf{B}`` the design matrix and ``\mathbf{W} = \text{diag}(1/\sigma_1^2,\ldots,1/\sigma_m^2)``
the weights matrix.

See also [`BHC`](@ref).

# Arguments

- `res::Vector{TaylorN{T}/Vector{OpticalResidual{T, TaylorN{T}}}`: vector of residuals.
- `w::Vector{T}`: vector of weights.
- `x_0::Vector{T}`: first guess.
- `niters::Int`: number of iterations.

!!! reference
    See sections 5.2 and 5.3 of https://doi.org/10.1017/CBO9781139175371.
"""
function diffcorr(res::Vector{TaylorN{T}}, w::Vector{T}, x0::Vector{T},
                  niters::Int = 5, idxs::AbstractVector{Int} = eachindex(x0)) where {T <: Real}
    # Degrees of freedom
    npar = length(x0)
    # Design matrix B, H array and normal matrix C
    B_mat, H_mat, C_mat = BHC(res, w, npar)
    # D matrix: transpose(B) * W * ξ
    D_mat = B_mat' * (w .* res)
    # ξTH_mat = ξTH(w, res, H_mat, npar)
    # Vector of x
    x = Matrix{T}(undef, npar, niters + 1)
    # First guess
    for i in axes(x, 2)
        x[:, i] .= x0
    end
    # Vector of errors
    error = Vector{T}(undef, niters + 1)
    # Error of first guess
    error[1] = T(Inf)
    # Iteration
    for i in 1:niters
        # Current x
        xi = x[:, i]
        # D matrix evaluated in xi
        D = D_mat(xi)[idxs]
        # C matrix evaluated in xi
        C = C_mat(xi)[idxs, idxs] #.+ ξTH_mat(xi)
        # Update rule
        Δx = - inv(C)*D
        # New x
        x[idxs, i+1] = xi[idxs] + Δx
        # Error
        error2 = ( (Δx') * (C*Δx) ) / npar
        if error2 ≥ 0
            error[i+1] = sqrt(error2)
        # The method do not converge
        else
            return LeastSquaresFit(false, x[:, i+1], inv(C), :diffcorr)
        end
    end
    # Index with the lowest error
    i = argmin(error)
    # x with the lowest error
    x_new = x[:, i]
    # Normal C matrix evaluated in x_new
    C = C_mat(x_new)
    # Covariance matrix
    Γ = inv(C[idxs, idxs])

    if any(diag(Γ) .< 0)
        return LeastSquaresFit(false, x_new, Γ, :diffcorr)
    else
        return LeastSquaresFit(true, x_new, Γ, :diffcorr)
    end
end

function diffcorr(res::Vector{OpticalResidual{T, TaylorN{T}}}, x0::Vector{T},
                  niters::Int = 5, idxs::AbstractVector{Int} = eachindex(x0)) where {T <: Real}
    # Unfold residuals and weights
    _res_, _w_ = unfold(res)

    return diffcorr(_res_, _w_, x0, niters, idxs)
end

@doc raw"""
    newtonls(res::Vector{TaylorN{T}}, w::Vector{T}, x0::Vector{T},
             niters::Int = 5) where {T <: Real}
    newtonls(res::Vector{OpticalResidual{T, TaylorN{T}}}, x0::Vector{T},
             niters::Int = 5) where {T <: Real}

Newton method subroutine for least-squares fitting. Returns an `LeastSquaresFit`
with the `niters`-th iteration
```math
\mathbf{x}_{k+1} = \mathbf{x}_k -
\left(\frac{\partial^2 Q}{\partial\mathbf{x}_0^2}\right)^{-1}
\frac{\partial Q}{\partial\mathbf{x}_0},
```
and the covariance matrix
```math
\mathbf{\Gamma} = \mathbf{C}^{-1},
```
where ``\mathbf{C} = \frac{m}{2}\frac{\partial^2 Q}{\partial\mathbf{x}_0^2}`` is the normal
matrix, ``Q = \frac{\chi^2}{m}`` is the mean square residual, ``m`` is the number of
observations and ``\mathbf{x}_0 = (x_1,\ldots,x_n)`` is the vector of initial conditions.

See also [`chi2`](@ref).

# Arguments

- `res::Vector{TaylorN{T}/Vector{OpticalResidual{T, TaylorN{T}}}`: vector of residuals.
- `w::Vector{T}`: vector of weights.
- `x_0::Vector{T}`: first guess.
- `niters::Int`: number of iterations.

!!! reference
    See sections 5.2 and 5.3 of https://doi.org/10.1017/CBO9781139175371.
"""
function newtonls(res::Vector{TaylorN{T}}, w::Vector{T}, x0::Vector{T},
                  niters::Int = 5, idxs::AbstractVector{Int} = eachindex(x0)) where {T <: Real}
    # Number of observations
    nobs = length(res)
    # Degrees of freedom
    npar = length(x0)
    # Mean square residual
    Q = chi2(res, w)/nobs
    # Vector of x
    x = Matrix{T}(undef, npar, niters + 1)
    # First guess
    for i in axes(x, 2)
        x[:, i] .= x0
    end
    # Vector of errors
    error = Vector{T}(undef, niters + 1)
    # Error of first guess
    error[1] = T(Inf)
    # Iteration
    for i in 1:niters
        # Current x
        xi = x[:, i]
        # Gradient of Q with respect to x
        dQ = TaylorSeries.gradient(Q)(xi)[idxs]
        # Hessian of Q with respect to x
        d2Q = TaylorSeries.hessian(Q, xi)[idxs, idxs]
        # Newton update rule
        Δx = - inv(d2Q)*dQ
        # New x
        x[idxs, i+1] = xi[idxs] + Δx
        # Normal matrix
        C = d2Q/(2/nobs) # C = d2Q/(2/m)
        # Error
        error2 = ( (Δx') * (C*Δx) ) / npar
        if error2 ≥ 0
            error[i+1] = sqrt(error2)
        # The method do not converge
        else
            return LeastSquaresFit(false, x[:, i+1], inv(C), :newton)
        end
    end
    # TO DO: study Gauss method solution dependence on jt order
    # TO DO: try even varorder
    # TO DO: study optimal number of iterations

    # Index with the lowest error
    i = argmin(error)
    # x with the lowest error
    x_new = x[:, i]
    # Normal matrix
    C = TaylorSeries.hessian(Q, x_new)/(2/nobs) # C = d2Q/(2/m)
    # Covariance matrix
    Γ = inv(C[idxs, idxs])

    if any(diag(Γ) .< 0)
        return LeastSquaresFit(false, x_new, Γ, :newton)
    else
        return LeastSquaresFit(true, x_new, Γ, :newton)
    end
end

function newtonls(res::Vector{OpticalResidual{T, TaylorN{T}}}, x0::Vector{T},
                  niters::Int = 5, idxs::AbstractVector{Int} = eachindex(x0)) where {T <: Real}
    # Unfold residuals and weights
    _res_, _w_ = unfold(res)

    return newtonls(_res_, _w_, x0, niters, idxs)
end

# In-place Levenberg-Marquardt hessian
function lmhessian!(_d2Q_::AbstractMatrix{T}, d2Q::AbstractMatrix{T}, λ::T) where {T <: AbstractFloat}
    k = 1 + λ
    for j in axes(d2Q, 2)
        for i in axes(d2Q, 1)
            if i == j
                _d2Q_[i, j] = k * d2Q[i, j]
            else
                _d2Q_[i, j] = d2Q[i, j]
            end
        end
    end
    return nothing
end

@doc raw"""
    levenbergmarquardt(res::Vector{TaylorN{T}}, w::Vector{T}, x0::Vector{T},
                       niters::Int = 500) where {T <: Real}
    levenbergmarquardt(res::Vector{OpticalResidual{T, TaylorN{T}}}, x0::Vector{T},
                       niters::Int = 500) where {T <: Real}

Levenberg-Marquardt method subroutine for least-squares fitting. Returns an `LeastSquaresFit`
with the best iteration of `niters` iterations.

See also [`chi2`](@ref).

# Arguments

- `res::Vector{TaylorN{T}/Vector{OpticalResidual{T, TaylorN{T}}}`: vector of residuals.
- `w::Vector{T}`: vector of weights.
- `x_0::Vector{T}`: first guess.
- `niters::Int`: number of iterations.

!!! reference
    See section 15.5.2 of https://books.google.com.mx/books?id=1aAOdzK3FegC&lpg=PP1&pg=PP1#v=onepage&q&f=false.
"""
function levenbergmarquardt(res::Vector{TaylorN{T}}, w::Vector{T}, x0::Vector{T},
                            niters::Int = 500, idxs::AbstractVector{Int} = eachindex(x0),
                            λ::T = 0.001) where {T <: Real}
    # Number of observations
    nobs = length(res)
    # Degrees of freedom
    npar = length(x0)
    # Normalized mean square residual
    Q = chi2(res, w)/nobs
    # Vector of Qs
    Qs = fill(T(Inf), niters + 1)
    # Gradient of Q with respect to x
    dQ = TaylorSeries.gradient(Q)
    # Pre-allocate memory
    _dQ_ = zeros(T, npar)
    _d2Q_ = zeros(T, npar, npar)
    x = Matrix{T}(undef, npar, niters + 1)
    # First guess
    for i in axes(x, 2)
        x[:, i] .= x0
    end
    # Iteration
    for i in 1:niters
        # Current x
        xi = x[:, i]
        # Current Q
        Qs[i] = Q(xi)
        # Convergence condition
        i > 1 && abs(Qs[i] - Qs[i-1]) / Qs[i-1] < 0.001 && break
        # Gradient of Q with respect to x
        _dQ_ .= dQ(xi)
        # Hessian of Q with respect to x
        d2Q = TaylorSeries.hessian(Q, xi)
        # Choose λ
        for _ in 1:niters
            # Modified Hessian
            lmhessian!(_d2Q_, d2Q, λ)
            # Levenberg-Marquardt step
            Δx = - inv(_d2Q_[idxs, idxs]) * _dQ_[idxs]
            # Update point
            x[idxs, i+1] = xi[idxs] + Δx
            # Choose λ
            if 0 < Q(x[:, i+1]) < Qs[i] && isposdef(_d2Q_[idxs, idxs])
                λ /= 10
                x[idxs, i+1] = xi[idxs] + Δx
                break
            else
                λ *= 10
            end
        end
    end
    # Index with the lowest error
    i = argmin(Qs)
    # x with the lowest error
    x_new = x[:, i]
    # Hessian of Q with respect to x
    d2Q = TaylorSeries.hessian(Q, x_new)
    # Normal matrix
    if isposdef(d2Q)
        C = d2Q/(2/nobs)
    else
        lmhessian!(_d2Q_, d2Q, λ)
        C = _d2Q_/(2/nobs) # C = d2Q/(2/m)
    end
    # Covariance matrix
    Γ = inv(C[idxs, idxs])

    if any(diag(Γ) .< 0)
        return LeastSquaresFit(false, x_new, Γ, :lm)
    else
        return LeastSquaresFit(true, x_new, Γ, :lm)
    end
end

function levenbergmarquardt(res::Vector{OpticalResidual{T, TaylorN{T}}}, x0::Vector{T},
                            niters::Int = 500, idxs::AbstractVector{Int} = eachindex(x0),
                            λ::T = 0.001) where {T <: Real}
    # Unfold residuals and weights
    _res_, _w_ = unfold(res)

    return levenbergmarquardt(_res_, _w_, x0, niters, idxs, λ)
end

@doc raw"""
    tryls(res::Vector{OpticalResidual{T, TaylorN{T}}}, x0::Vector{T},
          niters::Int = 5) where {T <: Real}

Return the best least squares fit between three routines: [`newtonls`](@ref),
[`diffcorr`](@ref) and [`levenbergmarquardt`](@ref).

# Arguments

- `res::Vector{OpticalResidual{T, TaylorN{T}}}`: vector of residuals.
- `x_0::Vector{T}`: first guess.
- `niters::Int`: number of iterations.
"""
function tryls(res::Vector{OpticalResidual{T, TaylorN{T}}}, x0::Vector{T},
               niters::Int = 5, idxs::AbstractVector{Int} = eachindex(x0),
               order::Vector{Symbol} = [:newton, :diffcorr, :lm]) where {T <: Real}
    # Allocate memory
    fit = zero(LeastSquaresFit{T})
    # Least squares methods in order
    for i in eachindex(order)
        if order[i] == :newton
            fit = newtonls(res, x0, niters, idxs)
        elseif order[i] == :diffcorr
            fit = diffcorr(res, x0, niters, idxs)
        elseif order[i] == :lm
            fit = levenbergmarquardt(res, x0, niters, idxs)
        end
        fit.success && break
    end

    return fit
end

# TO DO: update / deprecate the following three functions

@doc raw"""
    newtonls_Q(Q, nobs, x0, niters=5)

Does the same as `newtonls`, but recives ``Q`` as an argument, instead of computing it.
Returns the `niters`-th iteration and the covariance matrix ``\Gamma``.

See sections 5.2 and 5.3 of https://doi.org/10.1017/CBO9781139175371.

See also [`newtonls`](@ref).

# Arguments

- `Q`: Mean square residual.
- `nobs`: Number of observations.
- `x_0`: First guess for the initial conditions.
- `niters`: Number of iterations.
"""
function newtonls_Q(Q, nobs, x0, niters=5)
    # Number of observations
    npar = length(x0)
    # First guess
    x_new = x0
    # Iteration
    for i in 1:niters
        # Gradient of Q with respect to x_0
        dQ = TaylorSeries.gradient(Q)(x_new)
        # Hessian of Q with respect to x_0
        d2Q = TaylorSeries.hessian(Q, x_new)
        # Newton update rule
        Δx = - inv(d2Q)*dQ
        x_new = x_new + Δx
        # Normal matrix
        C = d2Q/(2/nobs) # C = d2Q/(2/m)
        @show sqrt(((Δx')*(C*Δx))/npar)
    end
    # Normal matrix
    C = TaylorSeries.hessian(Q, x_new)/(2/nobs) # C = d2Q/(2/m)
    # Covariance matrix
    Γ = inv(C)

    return x_new, Γ
end

@doc raw"""
    newtonls_6v(res, w, x0, niters=5)

Specialized version of `newtonls` on 6 variables for parametrized orbit determination
with respect to ``A_2`` Yarkovsky non-gravitational coefficient. Returns the `niters`-th
iteration and the covariance matrix ``\Gamma``.

See sections 5.2 and 5.3 of https://doi.org/10.1017/CBO9781139175371.

See also [`newtonls`](@ref).

# Arguments

- `res`: Vector of residuals.
- `w`: Vector of weights.
- `x_0`: First guess for the initial conditions.
- `niters`: Number of iterations.
"""
function newtonls_6v(res, w, x0, niters=5)
    # Have as many residuals as weights
    @assert length(res) == length(w)
    # Number of observations
    nobs = length(res)
    # Degrees of freedom
    npar = 6 # length(x0)
    # Mean square residual
    Q = chi2(res, w)/nobs
    # First guess
    x_new = x0
    # Iteration
    for i in 1:niters
        # Gradient of Q with respect to x_0
        dQ = TaylorSeries.gradient(Q)(x_new)[1:6]
        # Hessian of Q with respect to x_0
        d2Q = TaylorSeries.hessian(Q, x_new)[1:6,1:6]
        # Newton update rule
        Δx = - inv(d2Q)*dQ
        x_new[1:6] = x_new[1:6] + Δx
        # Normal matrix
        C = d2Q/(2/nobs) # C = d2Q/(2/m)
        @show sqrt(((Δx')*(C*Δx))/npar)
    end
    # Normal matrix
    C = TaylorSeries.hessian(Q, x_new)/(2/nobs) # C = d2Q/(2/m)
    # Covariance matrix
    Γ = inv(C[1:6,1:6])

    return x_new, Γ
end

@doc raw"""
    newtonls_A2(res, w, x0, niters=5)

Specialized version of `newtonls` with the Newton method only over the seventh degree of
freedom, i.e., with respect to ``A_2`` Yarkovsky non-gravitational coefficient. Returns the
`niters`-th iteration, the covariance matrix ``\Gamma`` and the normal matrix ``C``.

See sections 5.2 and 5.3 of https://doi.org/10.1017/CBO9781139175371.

See also [`newtonls`](@ref).

# Arguments

- `res`: Vector of residuals.
- `w`: Vector of weights.
- `x_0`: First guess for the initial conditions.
- `niters`: Number of iterations.
"""
function newtonls_A2(res, w, x0, niters=5)
    # Have as many residuals as weights
    @assert length(res) == length(w)
    # Number of observations
    nobs = length(res)
    # Degrees of freedom
    npar = length(x0)
    # Mean square residual
    Q = chi2(res, w)/nobs
    # First guess
    x_new = x0
    # Iteration
    for i in 1:niters
        # Gradient of Q with respect to x_0
        dQ = TaylorSeries.gradient(Q)(x_new)
        # Hessian of Q with respect to x_0
        d2Q = TaylorSeries.hessian(Q, x_new)
        # Newton update rule
        Δx = - inv(d2Q)*dQ
        x_new[7] = x_new[7] + Δx[7]
        # Normal matrix
        C = d2Q/(2/nobs) # C = d2Q/(2/m)
        @show sqrt(((Δx')*(C*Δx))/npar)
    end
    # Normal matrix
    C = TaylorSeries.hessian(Q, x_new)/(2/nobs) # C = d2Q/(2/m)
    # Covariance matrix
    Γ = inv(C)

    return x_new, Γ, C
end
