@doc raw"""
    OrbitFit{T <: Real}

A least squares fit.

# Fields
- `success::Bool`: wheter the routine converged or not.
- `x::Vector{T}`: deltas that minimize the objective function.
- `Γ::Matrix{T}`: covariance matrix.
- `routine::Symbol`: minimization routine (`:newton` or `:diffcorr`).
"""
@auto_hash_equals struct OrbitFit{T <: Real}
    success::Bool
    x::Vector{T}
    Γ::Matrix{T}
    routine::Symbol
    OrbitFit{T}(success::Bool, x::Vector{T}, Γ::Matrix{T}, routine::Symbol) where {T <: Real} = new{T}(success, x, Γ, routine)
end

OrbitFit(success::Bool, x::Vector{T}, Γ::Matrix{T}, routine::Symbol) where {T <: Real} = OrbitFit{T}(success, x, Γ, routine)

# Print method for OrbitFit
# Examples:
# Succesful Newton
# Succesful differential corrections
function show(io::IO, fit::OrbitFit{T}) where {T <: Real}
    success_s = fit.success ? "Succesful" : "Unsuccesful"
    routine_s = fit.routine == :newton ? "Newton" : "differential corrections"
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
    outlier_rejection(ξs::Vector{OpticalResidual{T}}, fit::OrbitFit{T}; χ2_rec::T = 7., χ2_rej::T = 8.,
                      α::T = 0.25) where {T <: Real}

Outlier rejection algorithm.

# Arguments

- `ξs::Vector{OpticalResidual{T}}`: vector of residuals.
- `fit::OrbitFit{T}`: least squares fit.
- `χ2_rec::T`: recovery conditions.
- `χ2_rej::T`: rejection condition.
- `α::T`: scaling factor for maximum chi.

!!! reference
    See https://doi.org/10.1016/S0019-1035(03)00051-4.
"""
function outlier_rejection(ξs::Vector{OpticalResidual{T, TaylorN{T}}}, fit::OrbitFit{T}; χ2_rec::T = 7., χ2_rej::T = 8.,
                           α::T = 0.25) where {T <: Real}
    # Number of residuals
    N = length(ξs)
    # Vector of chi2s
    χ2s = Vector{T}(undef, N)

    for i in eachindex(χ2s)
        # Current observation covariance matrix
        γ = diagm([ξs[i].w_α/ξs[i].relax_factor, ξs[i].w_δ/ξs[i].relax_factor])
        # Current model matrix
        A = hcat(TS.gradient(ξs[i].ξ_α)(fit.x), TS.gradient(ξs[i].ξ_δ)(fit.x))
        # Outlier sign
        outlier_sign = ξs[i].outlier*2-1
        # Current residual covariance matrix
        γ_ξ = γ + outlier_sign*(A')*fit.Γ*A
        # Current residual
        ξ = [ξs[i].ξ_α(fit.x), ξs[i].ξ_δ(fit.x)]
        # Current chi2
        χ2s[i] = ξ' * inv(γ_ξ) * ξ
    end

    χ2_max = maximum(χ2s)
    new_ξs = Vector{OpticalResidual{T, TaylorN{T}}}(undef, N)
    N_sel = count(x -> !x.outlier, ξs)

    for i in eachindex(χ2s)
        if χ2s[i] > max(χ2_rej + carpino_smoothing(N_sel), α*χ2_max)
            new_ξs[i] = OpticalResidual(ξs[i].ξ_α, ξs[i].ξ_δ, ξs[i].w_α, ξs[i].w_δ, ξs[i].relax_factor, true)
        elseif χ2s[i] < χ2_rec
            new_ξs[i] = OpticalResidual(ξs[i].ξ_α, ξs[i].ξ_δ, ξs[i].w_α, ξs[i].w_δ, ξs[i].relax_factor, false)
        else
            new_ξs[i] = ξs[i]
        end
    end

    return new_ξs
end

@doc raw"""
    project(y::Vector{TaylorN{T}}, fit::OrbitFit{T}) where {T <: Real}

Project `fit`'s covariance matrix into `y`.
"""
function project(y::Vector{TaylorN{T}}, fit::OrbitFit{T}) where {T <: Real}
    J = Matrix{T}(undef, get_numvars(), length(y))
    for i in eachindex(y)
        J[:, i] = TS.gradient(y[i])(fit.x)
    end
    return (J') * fit.Γ * J
end

@doc raw"""
    chi2(res::Vector{U}, w::Vector{T}) where {T <: Real, U <: Number}

Returns the chi square
```math
\chi^2 = \sum_{i=1}^m \frac{ \xi_i^2}{\sigma_i^2},
```
where ``\mathbf{w} = (1/\sigma_1^2,\ldots,1/\sigma_m^2)^T`` and ``\mathbf{\xi} = (\xi_1,\ldots,\xi_m)^T`` 
are the vectors of weights and residuals respectively. 

# Arguments

- `res::Vector{U}`: Vector of residuals.
- `w::Vector{T}`: Vector of weights. 
"""
function chi2(res::Vector{U}, w::Vector{T}) where {T <: Real, U <: Number}
    # Have as many residuals as weights
    @assert length(res) == length(w)
    # Chi square
    return sum(w .* (res.^2))
end

nms(res, w) = chi2(res, w) / length(res)

@doc raw"""
    nrms(res::Vector{U}, w::Vector{T}) where {T <: Real, U <: Number}
    nrms(res::Vector{OpticalResidual{T, U}}) where {T <: Real, U <: Number}
    nrms(res::Vector{OpticalResidual{T, TaylorN{T}}}, fit::OrbitFit{T}) where {T <: Real}

Returns the normalized root mean square error
```math
\texttt{NRMS} = \sqrt{\frac{\chi^2}{m}},
```
where ``\chi^2`` is the chi square and ``\mathbf{\xi} = (\xi_1,\ldots,\xi_m)^T`` is the vector
of residuals. 

See also [`chi2`](@ref).

# Arguments

- `res`: Vector of residuals.
- `w`: Vector of weights. 
- `fit`: least squares fit.

"""
function nrms(res::Vector{U}, w::Vector{T}) where {T <: Real, U <: Number}
    # Have as many residuals as weights
    @assert length(res) == length(w)
    # Normalized root mean square error
    return sqrt( chi2(res, w)/length(res) )
end

function nrms(res::Vector{OpticalResidual{T, U}}) where {T <: Real, U <: Number}
    res_, w = unfold(res)
    return nrms(res_, w)
end

function nrms(res::Vector{OpticalResidual{T, TaylorN{T}}}, fit::OrbitFit{T}) where {T <: Real}
    res_, w = unfold(res)
    return nrms(res_(fit.x), w)
end

@doc raw"""
    BHC(res, w, npar)

Returns the ``\mathbf{B}``, ``\mathbf{H}`` and ``\mathbf{C}`` arrays
```math
\mathbf{B} = \frac{\partial\mathbf{\xi}}{\partial\mathbf{x}_0}(\mathbf{x}_0), \quad
\mathbf{H} = \frac{\partial^2\mathbf{\xi}}{\partial\mathbf{x}_0^2}(\mathbf{x}_0) \quad \text{and} \quad
\mathbf{C} = \mathbf{B}^T\mathbf{W}\mathbf{B}, 
```
where ``\mathbf{x}_0 = (x_1,\ldots,x_n)^T`` and ``\mathbf{\xi} = (\xi_1,\ldots,\xi_m)^T`` 
are the vectors of initial conditions and residuals respectively; and ``\mathbf{W} = \text{diag}(1/\sigma_1^2,\ldots,1/\sigma_m^2)``
is the weights matrix. 

``\mathbf{B}`` is called the design matrix and is of size ``m\times n``, ``\mathbf{H}`` is a three index
array of size ``m\times n\times n`` and ``\mathbf{C}`` is called the normal matrix and is of size ``n\times n``.
See sections 5.2 and 5.3 of https://doi.org/10.1017/CBO9781139175371.

# Arguments

- `res`: Vector of residuals.
- `w`: Vector of weights. 
- `npar`: Degrees of freedom ``n``. 
"""
function BHC(res, w, npar)
    # Number of observations
    nobs = length(res)

    # Allocate memory for the three arrays
    B_mat = Matrix{TaylorN{Float64}}(undef, nobs, npar)
    H_mat = Array{TaylorN{Float64}}(undef, nobs, npar, npar)
    C_mat = Array{TaylorN{Float64}}(undef, npar, npar)

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
    diffcorr(ξs::Vector{OpticalResidual{T}}, x0::Vector{T}, niters::Int = 5) where {T <: Real}

Differential corrections subroutine for least-squares fitting. Returns the `niters`-th 
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

See sections 5.2 and 5.3 of https://doi.org/10.1017/CBO9781139175371. 

See also [`BHC`](@ref).

# Arguments

- `ξs::Vector{OpticalResidual{T, TaylorN{T}}}`: Vector of residuals.
- `x_0::Vector{T}`: First guess for the initial conditions.
- `niters::Int`: Number of iterations.
"""
function diffcorr(ξs::Vector{OpticalResidual{T, TaylorN{T}}}, x0::Vector{T}, niters::Int = 5) where {T <: Real}
    # Unfold residuals and weights
    res, w = unfold(ξs)
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
    x[:, 1] = x0
    # Vector of errors 
    error = Vector{T}(undef, niters + 1)
    # Error of first guess 
    error[1] = T(Inf)
    # Iteration
    for i in 1:niters
        # Current x
        xi = x[:, i]
        # D matrix evaluated in xi
        D = D_mat(xi)
        # C matrix evaluated in xi
        C = C_mat(xi) #.+ ξTH_mat(xi)
        # Update rule
        Δx = - inv(C)*D
        # New x 
        x[:, i+1] = xi + Δx 
        # Error 
        error2 = ( (Δx') * (C*Δx) ) / npar
        if error2 ≥ 0
            error[i+1] = sqrt(error2)
        # The method do not converge
        else 
            return OrbitFit(false, x[:, i+1], inv(C), :diffcorr)
        end 
    end
    # Index with the lowest error 
    i = argmin(error)
    # x with the lowest error 
    x_new = x[:, i]
    # Normal C matrix evaluated in x_new
    C = C_mat(x_new)
    # Covariance matrix
    Γ = inv(C)
    
    return OrbitFit(true, x_new, Γ, :diffcorr)
end

@doc raw"""
    newtonls(res, w, x0, niters=5)

Newton method subroutine for least-squares fitting. Returns the `niters`-th iteration
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

See sections 5.2 and 5.3 of https://doi.org/10.1017/CBO9781139175371. 

See also [`chi2`](@ref).

# Arguments

- `ξs::Vector{OpticalResidual{T, TaylorN{T}}}`: Vector of residuals.
- `x_0::Vector{T}`: First guess for the initial conditions.
- `niters::Int`: Number of iterations.
"""
function newtonls(ξs::Vector{OpticalResidual{T, TaylorN{T}}}, x0::Vector{T}, niters::Int = 5) where {T <: Real}
    # Unfold residuals and weights
    res, w = unfold(ξs)
    # Number of observations
    nobs = length(res)
    # Degrees of freedom
    npar = length(x0)
    # Mean square residual
    Q = chi2(res, w)/nobs
    # Vector of x 
    x = Matrix{T}(undef, npar, niters + 1) 
    # First guess
    x[:, 1] = x0
    # Vector of errors 
    error = Vector{T}(undef, niters + 1)
    # Error of first guess 
    error[1] = T(Inf)
    # Iteration
    for i in 1:niters
        # Current x
        xi = x[:, i]
        # Gradient of Q with respect to x
        dQ = TaylorSeries.gradient(Q)(xi)
        # Hessian of Q with respect to x
        d2Q = TaylorSeries.hessian(Q, xi)
        # Newton update rule
        Δx = - inv(d2Q)*dQ
        # New x 
        x[:, i+1] = xi + Δx
        # Normal matrix
        C = d2Q/(2/nobs) # C = d2Q/(2/m)
        # Error 
        error2 = ( (Δx') * (C*Δx) ) / npar
        if error2 ≥ 0
            error[i+1] = sqrt(error2)
        # The method do not converge
        else 
            return OrbitFit(false, x[:, i+1], inv(C), :newton)
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
    Γ = inv(C)

    return OrbitFit(true, x_new, Γ, :newton)
end

@doc raw"""
    tryls(ξs::Vector{OpticalResidual{T}}, x0::Vector{T}, niters::Int = 5) where {T <: Real}

Returns the best least squares fit between two routines: [`newtonls`](@ref) and [`diffcorr`](@ref).

# Arguments

- `ξs::Vector{OpticalResidual{T, TaylorN{T}}}`: Vector of residuals.
- `x_0::Vector{T}`: First guess for the initial conditions.
- `niters::Int`: Number of iterations.
"""
function tryls(ξs::Vector{OpticalResidual{T, TaylorN{T}}}, x0::Vector{T}, niters::Int = 5) where {T <: Real}
    # Newton's method
    fit_1 = newtonls(ξs, x0, niters)
    # Differential corrections
    fit_2 = diffcorr(ξs, x0, niters)
    if fit_1.success && fit_2.success
        Q_1 = nrms(ξs, fit_1)
        Q_2 = nrms(ξs, fit_2)
        if Q_1 <= Q_2
            return fit_1
        else
            return fit_2
        end
    elseif fit_1.success
        return fit_1
    elseif fit_2.success
        return fit_2
    else
        return fit_1
    end

end

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
