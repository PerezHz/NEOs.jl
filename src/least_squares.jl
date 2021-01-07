#normalized root mean square error
function nrms(res, w)
    @assert length(res) == length(w)
    return sqrt( sum(res .* w .* res)/length(res) )
end

# chi square
function chi2(res, w)
    @assert length(res) == length(w)
    return sum(w .* (res.^2))
end

# B and C matrices from Milani and Gronchi (2010)
function BHC(res, w, npar)
    nobs = length(res)
    B_mat = Matrix{TaylorN{Float64}}(undef, nobs, npar)
    H_mat = Array{TaylorN{Float64}}(undef, nobs, npar, npar)
    C_mat = Array{TaylorN{Float64}}(undef, npar, npar)
    for i in 1:nobs
        B_mat[i,:] .= TaylorSeries.gradient(res[i])
    end
    for i in 1:nobs
        for j in 1:npar
                H_mat[i,j,:] .= TaylorSeries.gradient(B_mat[i,j])
        end
    end
    sqrtw_B = sqrt.(w) .* B_mat
    C_mat .= (sqrtw_B') * sqrtw_B
    return B_mat, H_mat, C_mat
end

# transpose(ξ)*H matrix, with weights
function ξTH(w, res, H_mat, npar)
    ξTHv = Array{Float64}(undef, npar, npar)
    for j in 1:npar
        for i in 1:npar
           ξTHv[i,j] = (w .* res)' * (H_mat[:,i,j])
        end
    end
    return ξTHv
end

# differential corrections subroutine for least-squares fitting
function diffcorr(res, w, x0, niters=5)
    @assert length(res) == length(w)
    nobs = length(res)
    npar = length(x0)
    B_mat, H_mat, C_mat = BHC(res, w, npar)
    D_mat = B_mat' * (w .* res)
    # ξTH_mat = ξTH(w, res, H_mat, npar)
    x_new = x0
    for i in 1:niters
        D = D_mat(x_new)
        C = C_mat(x_new) #.+ ξTH_mat(x_new)
        Δx = - inv(C)*D
        x_new = x_new .+ Δx
        @show sqrt(((Δx')*(C*Δx))/npar)
    end
    C = C_mat(x_new)
    Γ = inv(C)
    return x_new, Γ
end

# Newton method subroutine for least-squares fitting
function newtonls(res, w, x0, niters=5)
    @assert length(res) == length(w)
    nobs = length(res)
    npar = length(x0)
    Q = sum(w .* (res.^2))/nobs
    x_new = x0
    for i in 1:niters
        dQ = TaylorSeries.gradient(Q)(x_new)
        d2Q = TaylorSeries.hessian(Q, x_new)
        Δx = - inv(d2Q)*dQ
        x_new = x_new + Δx
        C = d2Q/(2/nobs) # C = d2Q/(2/m)
        @show sqrt(((Δx')*(C*Δx))/npar)
    end
    C = TaylorSeries.hessian(Q, x_new)/(2/nobs) # C = d2Q/(2/m)
    Γ = inv(C)
    return x_new, Γ
end

function newtonls_Q(Q, nobs, x0, niters=5)
    npar = length(x0)
    x_new = x0
    for i in 1:niters
        dQ = TaylorSeries.gradient(Q)(x_new)
        d2Q = TaylorSeries.hessian(Q, x_new)
        Δx = - inv(d2Q)*dQ
        x_new = x_new + Δx
        C = d2Q/(2/nobs) # C = d2Q/(2/m)
        @show sqrt(((Δx')*(C*Δx))/npar)
    end
    C = TaylorSeries.hessian(Q, x_new)/(2/nobs) # C = d2Q/(2/m)
    Γ = inv(C)
    return x_new, Γ
end

# specialized version of newtonls on 6 variables for parametrized orbit determination wrt A2 nongrav coefficient
function newtonls_6v(res, w, x0, niters=5)
    @assert length(res) == length(w)
    nobs = length(res)
    npar = 6 # length(x0)
    Q = sum(w .* (res.^2))/nobs
    x_new = x0
    for i in 1:niters
        dQ = TaylorSeries.gradient(Q)(x_new)[1:6]
        d2Q = TaylorSeries.hessian(Q, x_new)[1:6,1:6]
        Δx = - inv(d2Q)*dQ
        x_new[1:6] = x_new[1:6] + Δx
        C = d2Q/(2/nobs) # C = d2Q/(2/m)
        @show sqrt(((Δx')*(C*Δx))/npar)
    end
    C = TaylorSeries.hessian(Q, x_new)/(2/nobs) # C = d2Q/(2/m)
    Γ = inv(C[1:6,1:6])
    return x_new, Γ
end

function newtonls_A2(res, w, x0, niters=5)
    @assert length(res) == length(w)
    nobs = length(res)
    npar = length(x0)
    Q = sum(w .* (res.^2))/nobs
    x_new = x0
    for i in 1:niters
        dQ = TaylorSeries.gradient(Q)(x_new)
        d2Q = TaylorSeries.hessian(Q, x_new)
        Δx = - inv(d2Q)*dQ
        x_new[7] = x_new[7] + Δx[7]
        C = d2Q/(2/nobs) # C = d2Q/(2/m)
        @show sqrt(((Δx')*(C*Δx))/npar)
    end
    C = TaylorSeries.hessian(Q, x_new)/(2/nobs) # C = d2Q/(2/m)
    Γ = inv(C)
    return x_new, Γ, C
end
