#normalized root mean square error
function nrms(res, w)
    @assert length(res) == length(w)
    return sqrt( sum(res .* w .* res)/length(res) )
end

#normalized root mean square error, squared
function nrms2(res, w)
    @assert length(res) == length(w)
    return sum(w .* (res.^2))/length(res)
end

# B and C matrices from Milani and Gronchi (2010)
function BC(res, w, x)
    B = Matrix{TaylorN{Float64}}(undef, length(res), length(x))
    for i in 1:length(res)
        B[i,:] .= TaylorSeries.gradient(res[i])
    end
    sqrtw_B = sqrt.(w) .* B(x)
    C = sqrtw_B' * sqrtw_B
    return B, C # B is a TaylorN matrix; C is a Float64 matrix
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
        #B, C = BC(res, w, x_new)
        C = d2Q/(2/nobs) # C = (2/m)*d2Q
        @show sqrt(((Δx')*(C*Δx))/npar)
    end
    C = TaylorSeries.hessian(Q, x_new)/(2/nobs) # C = (2/m)*d2Q
    #B, C = BC(res, w, x_new)
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
        C = d2Q/(2/nobs) # C = (2/m)*d2Q
        @show sqrt(((Δx')*(C*Δx))/npar)
    end
    C = TaylorSeries.hessian(Q, x_new)/(2/nobs) # C = (2/m)*d2Q
    Γ = inv(C)
    return x_new, Γ
end
