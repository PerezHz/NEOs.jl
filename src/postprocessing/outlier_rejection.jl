include("bplane.jl")

@doc raw"""
    residual_norm(x::OpticalResidual{T, T}) where {T <: Real}

Return the contribution of `x` to the nrms.
"""
residual_norm(x::OpticalResidual{T, T}) where {T <: Real} = x.w_α * x.ξ_α^2 / x.relax_factor + x.w_δ * x.ξ_δ^2 / x.relax_factor

@doc raw"""
    outlier_rejection(radec::Vector{RadecMPC{T}}, sol::NEOSolution{T, T},
                      params::NEOParameters{T}; dynamics::D = newtonian!) where {T <: Real, D}

Refine an orbit, computed by [`orbitdetermination`](@ref), via propagation and/or outlier rejection.

# Arguments

- `radec::Vector{RadecMPC{T}}`: vector of observations.
- `sol::NEOSolution{T, T}`: orbit to be refined.
- `params::NEOParameters{T}`: see `Outlier Rejection Parameters` of [`NEOParameters`](@ref).
- `dynamics::D`: dynamical model.

!!! warning
    This function will change the (global) `TaylorSeries` variables.
"""
function outlier_rejection(radec::Vector{RadecMPC{T}}, sol::NEOSolution{T, T},
                           params::NEOParameters{T}; dynamics::D = newtonian!) where {T <: Real, D}
    # Julian day to start propagation
    jd0 = sol.bwd.t0 + PE.J2000
    # Initial conditions (T)
    q0 = sol(sol.bwd.t0)
    # Scaling factors
    scalings = abs.(q0) ./ 10^6
    # Jet transport perturbation
    dq = scaled_variables("δx", scalings; order = params.jtlsorder)
    # Initial conditions (jet transport)
    q = q0 .+ dq
    # Propagation and residuals
    bwd, fwd, res = propres(radec, jd0, q, params; dynamics)
    iszero(length(res)) && return zero(NEOSolution{T, T})
    # Origin
    x0 = zeros(T, 6)
    # Orbit fit
    fit = tryls(res, x0, params.newtoniter)
    # Residuals space to barycentric coordinates jacobian
    J = Matrix(TS.jacobian(dq))
    # NRMS (with 0 outliers)
    Q_0 = nrms(res, fit)

    if Q_0 < 1
        return evalfit(NEOSolution(sol.tracklets, bwd, fwd, res, fit, J))
    end

    # Number of observations
    N_radec = length(radec)
    # Maximum allowed outliers
    max_drop = ceil(Int, N_radec * params.max_per / 100)
    # Boolean mask (0: included in fit, 1: outlier)
    new_outliers = BitVector(zeros(Int, N_radec))

    # Drop loop
    for i in 1:max_drop
        # Contribution of each residual to nrms
        norms = residual_norm.(res(fit.x))
        # Iterate norms from largest to smallest
        idxs = sortperm(norms, rev = true)

        for j in idxs
            if !new_outliers[j]
                # Drop residual
                new_outliers[j] = true
                # Update residuals
                res = OpticalResidual.(ra.(res), dec.(res), weight_ra.(res), weight_dec.(res),
                                       relax_factor.(res), new_outliers)
                # Update fit
                fit = tryls(res, x0, params.newtoniter)
                break
            end
        end
    end

    # Outliers
    idxs = Vector{Int}(undef, max_drop)
    # NRMS
    Qs = Vector{T}(undef, max_drop + 1)
    # Number of outliers
    N_outliers = Vector{T}(undef, max_drop + 1)

    # Recovery loop
    for i in 1:max_drop
        # NRMS of current fit
        Qs[i] = nrms(res, fit)
        # Number of outliers in current fit
        N_outliers[i] = float(max_drop - i + 1)
        # Contribution of each residual to nrms
        norms = residual_norm.(res(fit.x))
        # Minimum norm among outliers
        j = findmin(norms[new_outliers])[2]
        # Find residual with minimum norm
        j = findall(new_outliers)[j]
        # Add j-th residual to outliers list
        idxs[i] = j
        # Recover residual
        new_outliers[j] = false
        # Update residuals
        res = OpticalResidual.(ra.(res), dec.(res), weight_ra.(res), weight_dec.(res),
                               relax_factor.(res), new_outliers)
        # Update fit
        fit = tryls(res, x0, params.newtoniter)
    end
    # Add 0 outliers fit
    Qs[end] = Q_0
    N_outliers[end] = zero(T)

    # Outlier rejection cannot reduce Q
    if all(Qs .> 1.)
        # Reset boolean mask
        new_outliers[1:end] .= false
        # Update residuals
        res = OpticalResidual.(ra.(res), dec.(res), weight_ra.(res), weight_dec.(res),
                               relax_factor.(res), new_outliers)
        # Update fit
        fit = tryls(res, x0, params.newtoniter)

        return evalfit(NEOSolution(sol.tracklets, bwd, fwd, res, fit, J))
    end

    if max_drop > 1
        # Assemble points
        points = Matrix{T}(undef, 2, max_drop + 1)
        for i in eachindex(Qs)
            points[1, i] = Qs[i]
            points[2, i] = N_outliers[i]
        end
        # K-means clustering
        cluster = kmeans(points, 2; init = [1, max_drop + 1])
        # Index of smallest cluster
        i_0 = cluster.assignments[1]
        # Find last fit of smallest cluster
        i = findfirst(x -> x != i_0, cluster.assignments) - 1
        # Update outliers indexes
        idxs = idxs[i:end]
    end

    # Reset boolean mask
    new_outliers[1:end] .= false
    # Outliers
    new_outliers[idxs] .= true
    # Update residuals
    res = OpticalResidual.(ra.(res), dec.(res), weight_ra.(res), weight_dec.(res),
                           relax_factor.(res), new_outliers)
    # Update fit
    fit = tryls(res, x0, params.newtoniter)

    return evalfit(NEOSolution(sol.tracklets, bwd, fwd, res, fit, J))
end

function outlier_rejection_carpino03(radec::Vector{RadecMPC{T}}, sol::NEOSolution{T, T},
        params::NEOParameters{T}; dynamics::D = newtonian!) where {T <: Real, D}
    # Julian day to start propagation
    jd0 = sol.bwd.t0 + PE.J2000
    # Initial conditions (T)
    q0 = sol(sol.bwd.t0)
    # Scaling factors
    scalings = abs.(q0) ./ 10^6
    # Jet transport perturbation
    dq = scaled_variables("δx", scalings; order = params.jtlsorder)
    # Initial conditions (jet transport)
    q = q0 .+ dq
    # Propagation and residuals
    bwd, fwd, res = propres(radec, jd0, q, params; dynamics)
    # Origin
    x0 = zeros(T, 6)
    # Orbit fit
    # fit = tryls(res, x0, params.newtoniter, outrej=true)
    fit, new_outliers = newtonls_outrej(res, x0, params.newtoniter, outrej=true)
    # Residuals space to barycentric coordinates jacobian
    J = Matrix(TS.jacobian(dq))
    # Update residuals
    res = OpticalResidual.(ra.(res), dec.(res), weight_ra.(res), weight_dec.(res),
                           relax_factor.(res), new_outliers)
    return evalfit(NEOSolution(sol.tracklets, bwd, fwd, res, fit, J))
end