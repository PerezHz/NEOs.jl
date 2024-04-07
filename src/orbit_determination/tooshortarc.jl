# Propagate an orbit and compute residuals
function propres(radec::Vector{RadecMPC{T}}, jd0::U, q0::Vector{V}, params::NEOParameters{T};
                 dynamics::D = newtonian!)  where {D, T <: AbstractFloat, U <: Number, V <: Number}
    # Time of first (last) observation
    t0, tf = datetime2julian(date(radec[1])), datetime2julian(date(radec[end]))
    # Epoch (plain)
    _jd0_ = cte(cte(jd0))
    # Years in backward (forward) integration
    nyears_bwd = -(_jd0_ - t0 + params.bwdoffset) / yr
    nyears_fwd = (tf - _jd0_ + params.fwdoffset) / yr
    # Backward (forward) integration
    bwd = propagate(dynamics, jd0, nyears_bwd, q0, params)
    fwd = propagate(dynamics, jd0, nyears_fwd, q0, params)
    if !issuccessfulprop(bwd, t0 - _jd0_; tol = params.coeffstol) ||
       !issuccessfulprop(fwd, tf - _jd0_; tol = params.coeffstol)
        return bwd, fwd, Vector{OpticalResidual{T, U}}(undef, 0)
    end
    # O-C residuals
    res = residuals(radec, params;
                    xvs = et -> auday2kmsec(params.eph_su(et/daysec)),
                    xve = et -> auday2kmsec(params.eph_ea(et/daysec)),
                    xva = et -> bwdfwdeph(et, bwd, fwd))

    return bwd, fwd, res
end

@doc raw"""
    adam(radec::Vector{RadecMPC{T}}, A::AdmissibleRegion{T}, ρ::T, v_ρ::T,
         params::NEOParameters{T}; scale::Symbol = :linear, η::T = 25.0,
         μ::T = 0.75, ν::T = 0.9, ϵ::T = 1e-8, Qtol::T = 0.001, order::Int = 2,
         dynamics::D = newtonian!) where {T <: AbstractFloat, D}

Adaptative moment estimation (ADAM) minimizer of normalized mean square
residual over and admissible region `A`.

!!! warning
    This function will set the (global) `TaylorSeries` variables to `dx₁ dx₂ dx₃ dx₄ dx₅ dx₆`.

!!! reference
    See Algorithm 1 of https://doi.org/10.48550/arXiv.1412.6980.
"""
function adam(radec::Vector{RadecMPC{T}}, A::AdmissibleRegion{T}, ρ::T, v_ρ::T,
              params::NEOParameters{T}; scale::Symbol = :linear, η::T = 25.0,
              μ::T = 0.75, ν::T = 0.9, ϵ::T = 1e-8, Qtol::T = 0.001, order::Int = 2,
              dynamics::D = newtonian!) where {T <: AbstractFloat, D}
    # Initial time of integration [julian days]
    jd0 = datetime2julian(A.date)
    # Attributable elements
    ae = [A.α, A.δ, A.v_α, A.v_δ, zero(T), zero(T)]
    # Scaling factors
    scalings = Vector{T}(undef, 6)
    scalings[1:4] .= abs.(ae[1:4]) ./ 1e6
    if scale == :linear
        scalings[5] = (A.ρ_domain[2] - A.ρ_domain[1]) / 1_000
    elseif scale == :log
        x = log10(ρ)
        scalings[5] = (log10(A.ρ_domain[2]) - log10(A.ρ_domain[1])) / 1_000
    end
    scalings[6] = (A.v_ρ_domain[2] - A.v_ρ_domain[1]) / 1_000
    # Jet transport variables
    dae = scaled_variables("dx", scalings; order)
    # Maximum number of iterations
    maxiter = params.maxiter
    # Allocate memory
    ρs = Vector{T}(undef, maxiter+1)
    v_ρs = Vector{T}(undef, maxiter+1)
    Qs = fill(T(Inf), maxiter+1)
    # Origin
    x0 = zeros(T, 6)
    x1 = zeros(T, 6)
    # First momentum
    m = zeros(T, 2)
    _m_ = zeros(T, 2)
    # Second momentum
    n = zeros(T, 2)
    _n_ = zeros(T, 2)
    # Gradient descent
    for t in 1:maxiter+1
        # Current position in admissible region
        ρs[t] = ρ
        v_ρs[t] = v_ρ
        # Attributable elements
        if scale == :linear
            ae[5], ae[6] = ρ, v_ρ
            AE = ae + dae
        elseif scale == :log
            ae[5], ae[6] = x, v_ρ
            AE = [ae[1] + dae[1], ae[2] + dae[2], ae[3] + dae[3],
                  ae[4] + dae[4], 10^(ae[5] + dae[5]), ae[6] + dae[6]]
        end
        # Barycentric state vector
        q = attr2bary(A, AE, params)
        # Propagation and residuals
        # TO DO: `ρ::TaylorN` is too slow for `adam` due to evaluations
        # within the dynamical model
        _, _, res = propres(radec, jd0 - ρ/c_au_per_day, q, params; dynamics)
        iszero(length(res)) && break
        # Least squares fit
        fit = tryls(res, x0, 5, 1:4)
        x1 .= fit.x
        # Current Q
        Q = nms(res)
        Qs[t] = Q(x1)
        # Convergence condition
        t > 1 && abs(Qs[t] - Qs[t-1]) / Qs[t] < Qtol && break
        # Gradient of objective function
        g_t = TaylorSeries.gradient(Q)(x1)[5:6]
        # First momentum
        m .= μ * m + (1 - μ) * g_t
        _m_ .= m / (1 - μ^t)
        # Second momentum
        n .= ν * n + (1 - ν) * g_t .^ 2
        _n_ .= n / (1 - ν^t)
        # Step
        x1[5:6] = x1[5:6] - η * _m_ ./ (sqrt.(_n_) .+ ϵ)
        # Update values
        ae .= AE(x1)
        # Projection
        ρ, v_ρ = boundary_projection(A, ae[5], ae[6])
        if scale == :log
            x = log10(ρ)
        end
    end
    # Find point with smallest Q
    t = argmin(Qs)

    return T(ρs[t]), T(v_ρs[t]), T(Qs[t])
end

@doc raw"""
    ρminmontecarlo(radec::Vector{RadecMPC{T}}, A::AdmissibleRegion{T},
                   params::NEOParameters{T}; N_samples::Int = 25) where {T <: AbstractFloat}

Monte Carlo sampling over the left boundary of `A`.
"""
function ρminmontecarlo(radec::Vector{RadecMPC{T}}, A::AdmissibleRegion{T},
                        params::NEOParameters{T}; N_samples::Int = 25, dynamics::D=newtonian!) where {T <: AbstractFloat, D}
    # Initial time of integration [julian days]
    jd0 = datetime2julian(A.date)
    # Range lower bound
    ρ = A.ρ_domain[1]
    # Sample range rate
    v_ρs = LinRange(A.v_ρ_domain[1], A.v_ρ_domain[2], N_samples)
    # Allocate memory
    Qs = fill(Inf, N_samples)
    # Monte Carlo
    for i in eachindex(Qs)
        # Barycentric initial conditions
        q = topo2bary(A, ρ, v_ρs[i])
        # Propagation & residuals
        _, _, res = propres(radec, jd0, q, params; dynamics)
        iszero(length(res)) && continue
        # NRMS
        Qs[i] = nrms(res)
    end
    # Find solution with smallest Q
    t = argmin(Qs)

    return T(ρ), T(v_ρs[t]), T(Qs[t])
end

@doc raw"""
    tsals(A::AdmissibleRegion{T}, radec::Vector{RadecMPC{T}}, tracklets::Vector{Tracklet{T}},
          i::Int, ρ::T, v_ρ::T, params::NEOParameters{T}; maxiter::Int = 5) where {T <: AbstractFloat}

Used within [`tooshortarc`](@ref) to compute an orbit from a point in an
admissible region via least squares.

!!! warning
    This function will set the (global) `TaylorSeries` variables to `dx₁ dx₂ dx₃ dx₄ dx₅ dx₆`.
"""
function tsals(A::AdmissibleRegion{T}, radec::Vector{RadecMPC{T}}, tracklets::Vector{Tracklet{T}},
               i::Int, ρ::T, v_ρ::T, params::NEOParameters{T}; maxiter::Int = 5, dynamics::D = newtonian!,
               order = 6) where {T <: AbstractFloat, D}
    # Initial time of integration [julian days]
    # (corrected for light-time)
    jd0 = datetime2julian(A.date) - ρ / c_au_per_day
    # Barycentric initial conditions
    q0 = topo2bary(A, ρ, v_ρ)
    # Scaling factors
    scalings = abs.(q0) ./ 10^5
    # Jet transport variables
    dq = scaled_variables("dx", scalings; order)
    # Origin
    x0 = zeros(T, 6)
    # Subset of radec for orbit fit
    g_0 = i
    g_f = i
    idxs = indices(tracklets[i])
    sort!(idxs)
    # Allocate memory
    best_sol = zero(NEOSolution{T, T})
    best_Q = T(Inf)
    flag = false
    # Least squares
    for _ in 1:maxiter
        # Initial conditions
        q = q0 + dq
        # Propagation & residuals
        bwd, fwd, res = propres(radec, jd0, q, params; dynamics)
        iszero(length(res)) && break
        # Orbit fit
        fit = tryls(res[idxs], x0, params.niter)
        !fit.success && break
        # Right iteration
        for k in g_f+1:length(tracklets)
            extra = indices(tracklets[k])
            fit_new = tryls(res[idxs ∪ extra], x0, params.niter)
            if fit_new.success
                fit = fit_new
                idxs = vcat(idxs, extra)
                sort!(idxs)
                g_f = k
            else
                break
            end
        end
        # Left iteration
        for k in g_0-1:-1:1
            extra = indices(tracklets[k])
            fit_new = tryls(res[idxs ∪ extra], x0, params.niter)
            if fit_new.success
                fit = fit_new
                idxs = vcat(idxs, extra)
                sort!(idxs)
                g_0 = k
            else
                break
            end
        end
        # NRMS
        Q = nrms(res, fit)
        if length(idxs) == length(radec) && abs(best_Q - Q) < 0.1
            flag = true
        end
        # Update NRMS and initial conditions
        if Q < best_Q
            best_Q = Q
            best_sol = evalfit(NEOSolution(tracklets[g_0:g_f], bwd, fwd,
                               res[idxs], fit, scalings))
            flag && break
        else
            break
        end
        # Update values
        q0 = q(fit.x)
    end
    # Case: all solutions were unsuccesful
    if isinf(best_Q)
        return zero(NEOSolution{T, T})
    # Case: at least one solution was succesful
    else
        return best_sol
    end
end

# Order in which to check tracklets in tooshortarc
function tsatrackletorder(x::Tracklet{T}, y::Tracklet{T}) where {T <: AbstractFloat}
    if x.nobs == y.nobs
        return x.date > y.date
    else
        return x.nobs > y.nobs
    end
end

# Point in admissible region -> NEOSolution{T, T}
function _tooshortarc(A::AdmissibleRegion{T}, radec::Vector{RadecMPC{T}},
                      tracklets::Vector{Tracklet{T}}, i::Int, params::NEOParameters{T};
                      scale::Symbol = :log, dynamics::D = newtonian!) where {T <: AbstractFloat, D}
    # Center
    if scale == :linear
        ρ = sum(A.ρ_domain) / 2
    elseif scale == :log
        ρ = A.ρ_domain[1]
    end
    v_ρ = sum(A.v_ρ_domain) / 2
    # ADAM minimization over admissible region
    ρ, v_ρ, Q = adam(radec, A, ρ, v_ρ, params; scale, dynamics)
    if isinf(Q)
        return zero(NEOSolution{T, T})
    else
        # 6 variables least squares
        return tsals(A, radec, tracklets, i, ρ, v_ρ, params; maxiter = 5, dynamics)
    end
end

@doc raw"""
    tooshortarc(radec::Vector{RadecMPC{T}}, tracklets::Vector{Tracklet{T}},
                params::NEOParameters{T}; dynamics::D = newtonian!) where {T <: AbstractFloat, D}

Return initial conditions by minimizing the normalized root mean square residual
over the admissible region.

# Arguments

- `radec::Vector{RadecMPC{T}}`: vector of observations.
- `tracklets::Vector{Tracklet{T}},`: vector of tracklets.
- `params::NEOParameters{T}`: see `Admissible Region Parameters` of [`NEOParameters`](@ref).
- `dynamics::D`: dynamical model.

!!! warning
    This function will set the (global) `TaylorSeries` variables to `dx₁ dx₂ dx₃ dx₄ dx₅ dx₆`.
"""
function tooshortarc(radec::Vector{RadecMPC{T}}, tracklets::Vector{Tracklet{T}},
                     params::NEOParameters{T}; dynamics::D = newtonian!) where {T <: AbstractFloat, D}

    # Allocate memory for output
    best_sol = zero(NEOSolution{T, T})
    # Sort tracklets by tsatrackletorder
    idxs = sortperm(tracklets, lt = tsatrackletorder)

    # Iterate tracklets
    for i in idxs
        # Admissible region
        A = AdmissibleRegion(tracklets[i], params)
        iszero(A) && continue
        # See Table 1 of https://doi.org/10.1051/0004-6361/201732104
        if A.ρ_domain[2] < sqrt(10)
            sol1 = _tooshortarc(A, radec, tracklets, i, params;
                                scale = :log, dynamics)
            # Break condition
            nrms(sol1) < params.tsaQmax && return sol1
            sol2 = _tooshortarc(A, radec, tracklets, i, params;
                                scale = :linear, dynamics)
            # Break condition
            nrms(sol2) < params.tsaQmax && return sol2
        else
            sol1 = _tooshortarc(A, radec, tracklets, i, params;
                                scale = :linear, dynamics)
            # Break condition
            nrms(sol1) < params.tsaQmax && return sol1
            sol2 = _tooshortarc(A, radec, tracklets, i, params;
                                scale = :log, dynamics)
            # Break condition
            nrms(sol2) < params.tsaQmax && return sol2
        end
        # Update best solution
        if nrms(sol1) < nrms(best_sol)
            best_sol = sol1
        elseif nrms(sol2) < nrms(best_sol)
            best_sol = sol2
        end
    end

    return best_sol
end