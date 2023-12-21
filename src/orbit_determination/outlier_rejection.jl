@doc raw"""
    residual_norm(x::OpticalResidual{T, T}) where {T <: Real}

Return the contribution of `x` to the nrms.
"""
residual_norm(x::OpticalResidual{T, T}) where {T <: Real} = x.w_α * x.ξ_α^2 / x.relax_factor + x.w_δ * x.ξ_δ^2 / x.relax_factor

@doc raw"""
    gauss_refinement(radec::Vector{RadecMPC{T}}, sol::NEOSolution{T, T}, params::Parameters{T};
                     kwargs...) where {T <: AbstractFloat}

Refine an orbit computed by [`gaussinitcond`](@ref) via propagation and/or outlier rejection.

# Arguments

- `radec::Vector{RadecMPC{T}}`: vector of observations.
- `sol::NEOSolution{T, T}`: orbit to be refined.
- `params::Parameters{T}`: see [`Parameters`](@ref).

# Keyword arguments

- `max_per::T =  18.0`: maximum allowed rejection percentage.
- `niter::Int = 5`: number of iterations. 
- `varorder::Int`: order of jet transport perturbation. 
"""
function gauss_refinement(radec::Vector{RadecMPC{T}}, sol::NEOSolution{T, T}, params::Parameters{T};
                          max_per::T = 18.0, niter::Int = 5, varorder::Int = 5) where {T <: AbstractFloat}

    # Sun's ephemeris
    eph_su = selecteph(sseph, su)
    # Earth's ephemeris
    eph_ea = selecteph(sseph, ea)
    # Julian day to start propagation
    jd0 = sol.bwd.t0 + PE.J2000
    # Julian day of first (last) observation
    t0, tf = datetime2julian(date(radec[1])), datetime2julian(date(radec[end]))
    # Number of years in forward integration 
    nyears_fwd = (tf - jd0 + 2) / yr
    # Number of years in backward integration
    nyears_bwd = -(jd0 - t0 + 2) / yr
    # Dynamical function 
    dynamics = RNp1BP_pN_A_J23E_J2S_eph_threads!
    # Maximum number of steps
    params = Parameters(params; maxsteps = adaptative_maxsteps(radec))
    # Initial conditions (T)
    q00 = sol(sol.bwd.t0)
    # Jet transport perturbation
    dq = scaled_variables("δx", fill(1e-6, 6); order = varorder)
    # Initial conditions (jet transport)
    q0 = q00 .+ dq

    # Backward integration
    bwd = propagate(dynamics, jd0, nyears_bwd, q0, params)

    if bwd.t[end] > t0 - jd0
        return zero(NEOSolution{T, T})
    end

    # Forward integration
    fwd = propagate(dynamics, jd0, nyears_fwd, q0, params)

    if fwd.t[end] < tf - jd0
        return zero(NEOSolution{T, T})
    end

    # Residuals
    res = residuals(radec, params;
                    xvs = et -> auday2kmsec(eph_su(et/daysec)),
                    xve = et -> auday2kmsec(eph_ea(et/daysec)),
                    xva = et -> bwdfwdeph(et, bwd, fwd))
    # Orbit fit
    fit = tryls(res, zeros(get_numvars()), niter)

    # NRMS (with 0 outliers)
    Q_0 = nrms(res, fit)

    if Q_0 < 1
        return evalfit(NEOSolution(bwd, fwd, res, fit))
    end

    # Number of observations
    N_radec = length(radec)
    # Maximum allowed outliers
    max_drop = ceil(Int, N_radec * max_per / 100)
    # Boolean mask (0: included in fit, 1: outlier)
    new_outliers = BitVector(zeros(N_radec))

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
                fit = tryls(res, zeros(get_numvars()), niter)
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
        fit = tryls(res, zeros(get_numvars()), niter)
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
        fit = tryls(res, zeros(get_numvars()), niter)

        return evalfit(NEOSolution(bwd, fwd, res, fit)) 
    end

    if max_drop > 1
        # Assemble points
        points = Matrix{T}(undef, 2, max_drop + 1)
        for i in eachindex(Qs)
            points[1, i] = Qs[i]
            points[2, i] = N_outliers[i]
        end
        # K-means clustering
        cluster = kmeans(points, 2)
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
    fit = tryls(res, zeros(get_numvars()), niter)

    return evalfit(NEOSolution(bwd, fwd, res, fit))
end