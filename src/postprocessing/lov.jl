function covariance(od::ODProblem{D, T}, sol::NEOSolution{T, T},
    params::NEOParameters{T}) where {D, T <: Real}
    # Reference epoch [Julian days TDB]
    jd0 = epoch(sol) + PE.J2000
    # Plain initial condition
    q00 = sol(epoch(sol))
    # Scaling factors
    scalings = 3 * sigmas(sol)

    return covariance(od, jd0, q00, params; scalings)
end

function covariance(od::ODProblem{D, T}, jd0::T, q00::Vector{T}, params::NEOParameters{T};
    scalings::Vector{T} = abs.(q00) ./ 1e6) where {D, T <: Real}
    # Jet transport variables
    dq = set_variables(T, "dx"; order = 2, numvars = 6)
    # Jet transport initial condition
    q0 = q00 + scalings .* dq
    # Propagation and residuals
    _, _, res = propres(od, jd0, q0, params)
    # Covariance matrix in residuals space
    Q = nms(res)
    C = notout(res) * TS.hessian(Q)
    Γ_ξ = inv(Symmetric(C))
    # Covariance matrix in cartesian space
    J = TS.jacobian(scalings .* dq)
    Γ_c = Symmetric(J * Γ_ξ * J')

    return Γ_c
end

function lovcovariance(od::ODProblem{D, T}, jd0::T, q00::Vector{T},
    params::NEOParameters{T}; scalings::Vector{T} = abs.(q00) ./ 1e6,
    order::Int = 20) where {D, T <: Real}
    # Covariance matrix
    Γ_c = covariance(od, jd0, q00, params; scalings)
    # Greatest eigenpair
    E = eigen(Γ_c)
    k1, v1 = sqrt(E.values[end]), E.vectors[:, end]
    # Jet transport variable
    dq = set_variables(T, "dx"; order = order, numvars = 1)[1]
    # LOV initial condition
    q0 = q00 + k1 * v1 * dq
    # Propagation and residuals
    _, _, res = propres(od, jd0, q0, params)
    # Covariance matrix in residuals space
    Q = nms(res)
    C = notout(res) * TS.differentiate(Q, (2,))
    Γ_ξ = inv(C)
    # Covariance matrix in cartesian space
    J = TS.jacobian(q0 - constant_term.(q0))
    Γ_c = Symmetric(J * Γ_ξ * J')

    return Γ_c
end

function machseries(Γ::AbstractMatrix{TaylorN{T}}, order::Int) where {T <: Real}
    # Allocate memory
    λs = Vector{T}(undef, order+1)
    vs = Matrix{T}(undef, 6, order+1)
    As = Array{T}(undef, 6, 6, order+1)
    # 0th-order eigenvalues
    As[:, :, 1] = Symmetric(constant_term.(Γ))
    E0 = eigen(As[:, :, 1])
    λs[1], vs[:, 1] = E0.values[end], E0.vectors[:, end]
    # Coefficients matrix
    E = [0 vs[:, 1]'; vs[:, 1] λs[1]*I-As[:, :, 1]]
    # Main loop
    for k in 1:order
        As[:, :, k+1] = TS.differentiate.(Ref((k,)), Γ)
        y, z = zeros(T, 6), zero(T)
        for l in 0:k-1
            y += binomial(k, l) * As[:, :, k-l+1] * vs[:, l+1]
            if l >= 1
                y -= binomial(k, l) * vs[:, k-l+1] * λs[l+1]
                if l < k -1
                    z += binomial(k-1, l-1) * vs[:, k-l+1]' * vs[:, l+1]
                end
            end
        end
        Ek = E \ [-z/2; y]
        λs[k+1], vs[:, k+1] = Ek[1], Ek[2:end]
    end
    #=
    λ, v = Γ[1] * TaylorN(1, order), Γ[:, 1] * TaylorN(1, order)
    for k in 0:order
        λ.coeffs[k+1].coeffs[1] = λs[k+1]
        for j in 1:6
            v[j].coeffs[k+1].coeffs[1] = vs[j, k+1]
        end
    end
    =#
    λ, v = Taylor1(λs, order), [Taylor1(vs[i, :], order) for i in axes(vs, 1)]

    return λ, v
end

@taylorize function lov!(dq, q, lovparams, t)
    local od = lovparams[1]
    local jd0 = lovparams[2]
    local params = lovparams[3]
    local scalings = lovparams[4]
    local lovorder = lovparams[5]
    local Γ_c = lovcovariance(od, jd0, constant_term.(q), params; scalings,
        order = lovorder)
    local λ, v = machseries(Γ_c, lovorder)
    dq[1] = sqrt(λ) * v[1]
    dq[2] = sqrt(λ) * v[2]
    dq[3] = sqrt(λ) * v[3]
    dq[4] = sqrt(λ) * v[4]
    dq[5] = sqrt(λ) * v[5]
    dq[6] = sqrt(λ) * v[6]
end