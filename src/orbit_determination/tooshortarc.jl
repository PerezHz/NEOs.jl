# Times used within propres
function _proprestimes(radec::Vector{RadecMPC{T}}, jd0::U, params::NEOParameters{T}) where {T <: Real, U <: Number}
    # Time of first (last) observation
    t0, tf = datetime2julian(date(radec[1])), datetime2julian(date(radec[end]))
    # Epoch (plain)
    _jd0_ = cte(cte(jd0))
    # Years in backward (forward) integration
    nyears_bwd = -(_jd0_ - t0 + params.bwdoffset) / yr
    nyears_fwd = (tf - _jd0_ + params.fwdoffset) / yr

    return t0, tf, _jd0_, nyears_bwd, nyears_fwd
end

# Propagate an orbit and compute residuals
function propres(radec::Vector{RadecMPC{T}}, jd0::U, q0::Vector{V}, params::NEOParameters{T};
                 buffer::Union{Nothing, PropagationBuffer{T, U, V}} = nothing,
                 dynamics::D = newtonian!)  where {D, T <: Real, U <: Number, V <: Number}
    # Times of first/last observation, epoch and years in backward/forward propagation
    t0, tf, _jd0_, nyears_bwd, nyears_fwd = _proprestimes(radec, jd0, params)
    # Propagation buffer
    if isnothing(buffer)
        tlim = (t0 - JD_J2000 - params.bwdoffset, tf - JD_J2000 + params.fwdoffset)
        buffer = PropagationBuffer(dynamics, jd0, tlim, q0, params)
    end
    # Backward (forward) integration
    bwd = _propagate(dynamics, jd0, nyears_bwd, q0, buffer, params)
    fwd = _propagate(dynamics, jd0, nyears_fwd, q0, buffer, params)
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

# In-place method of propres
function propres!(res::Vector{OpticalResidual{T, U}}, radec::Vector{RadecMPC{T}}, jd0::V, q0::Vector{U},
    params::NEOParameters{T}; buffer::Union{Nothing, PropagationBuffer{T, U, V}} = nothing,
    dynamics::D = newtonian!)  where {D, T <: Real, U <: Number, V <: Number}
    # Times of first/last observation, epoch and years in backward/forward propagation
    t0, tf, _jd0_, nyears_bwd, nyears_fwd = _proprestimes(radec, jd0, params)
    # Propagation buffer
    if isnothing(buffer)
        tlim = (t0 - JD_J2000 - params.bwdoffset, tf - JD_J2000 + params.fwdoffset)
        buffer = PropagationBuffer(dynamics, jd0, tlim, q0, params)
    end
    # Backward (forward) integration
    bwd = _propagate(dynamics, jd0, nyears_bwd, q0, buffer, params)
    fwd = _propagate(dynamics, jd0, nyears_fwd, q0, buffer, params)
    if !issuccessfulprop(bwd, t0 - _jd0_; tol = params.coeffstol) ||
       !issuccessfulprop(fwd, tf - _jd0_; tol = params.coeffstol)
        res = Vector{OpticalResidual{T, U}}(undef, 0)
        return bwd, fwd
    end
    # O-C residuals
    residuals!(res, radec, params;
        xvs = et -> auday2kmsec(params.eph_su(et/daysec)),
        xve = et -> auday2kmsec(params.eph_ea(et/daysec)),
        xva = et -> bwdfwdeph(et, bwd, fwd))

    return bwd, fwd
end

@doc raw"""
    jtls(radec::Vector{RadecMPC{T}}, tracklets::Vector{Tracklet{T}}, jd0::V, q::Vector{U},
         g0::Int, gf::Int, params::NEOParameters{T}; maxiter::Int = 5,
         dynamics::D = newtonian!) where {D, T <: Real, U <: Number, V <: Number}

Compute an orbit via Jet Transport Least Squares.

# Arguments

- `radec::Vector{RadecMPC{T}}`: vector of observations.
- `tracklets::Vector{Tracklet{T}},`: vector of tracklets.
- `jd0::V`: reference epoch [julian days].
- `q::Vector{TaylorN{T}}`: jet transport initial condition.
- `g0/gf::Int`: indices of `tracklets` to start least squares fit.
- `params::NEOParameters{T}`: see ` Least Squares Fit Parameters` of [`NEOParameters`](@ref).
- `maxiter::Int`: maximum number of iterations.
- `dynamics::D`: dynamical model.
"""
function jtls(radec::Vector{RadecMPC{T}}, tracklets::Vector{Tracklet{T}}, jd0::V, q::Vector{TaylorN{T}},
              g0::Int, gf::Int, params::NEOParameters{T}; maxiter::Int = 5,
              dynamics::D = newtonian!) where {D, T <: Real, V <: Number}
    # Plain initial condition
    q0 = constant_term.(q)
    # JT tail
    dq = q - q0
    # Vector of O-C residuals
    res = Vector{OpticalResidual{T, TaylorN{T}}}(undef, length(radec))
    # Propagation buffer
    t0, tf = datetime2days(date(radec[1])), datetime2days(date(radec[end]))
    tlim = (t0 - params.bwdoffset, tf + params.fwdoffset)
    buffer = PropagationBuffer(dynamics, jd0, tlim, q, params)
    # Origin
    x0 = zeros(T, 6)
    # Subset of radec for orbit fit
    idxs = reduce(vcat, indices.(tracklets[g0:gf]))
    sort!(idxs)
    # Residuals space to barycentric coordinates jacobian
    J = Matrix(TS.jacobian(dq))
    # Best orbit
    best_sol = zero(NEOSolution{T, T})
    best_Q = T(Inf)
    # Convergence flag
    flag = false
    # Jet transport least squares
    for _ in 1:maxiter
        # Initial conditions
        q = q0 + dq
        # Propagation & residuals
        bwd, fwd = propres!(res, radec, jd0, q, params; buffer, dynamics)
        iszero(length(res)) && break
        # Orbit fit
        fit = tryls(res[idxs], x0, params.niter)
        !fit.success && break
        # Right iteration
        for k in gf+1:length(tracklets)
            extra = indices(tracklets[k])
            fit_new = tryls(res[idxs ∪ extra], x0, params.niter)
            if fit_new.success
                fit = fit_new
                idxs = vcat(idxs, extra)
                sort!(idxs)
                gf = k
            else
                break
            end
        end
        # Left iteration
        for k in g0-1:-1:1
            extra = indices(tracklets[k])
            fit_new = tryls(res[idxs ∪ extra], x0, params.niter)
            if fit_new.success
                fit = fit_new
                idxs = vcat(idxs, extra)
                sort!(idxs)
                g0 = k
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
            J .= TS.jacobian(dq, fit.x)
            best_sol = evalfit(NEOSolution(tracklets[g0:gf], bwd, fwd,
                               res[idxs], fit, J))
            flag && break
        else
            break
        end
        # Update initial condition
        q0 .= q(fit.x)
    end

    # Case: all solutions were unsuccesful
    if isinf(best_Q)
        return zero(NEOSolution{T, T})
    # Case: at least one solution was succesful
    else
        return best_sol
    end
end

@doc raw"""
    adam(radec::Vector{RadecMPC{T}}, A::AdmissibleRegion{T}, ρ::T, v_ρ::T,
         params::NEOParameters{T}; scale::Symbol = :linear, η::T = 25.0,
         μ::T = 0.75, ν::T = 0.9, ϵ::T = 1e-8, Qtol::T = 0.001, varorder::Int = 2,
         dynamics::D = newtonian!) where {T <: AbstractFloat, D}

Adaptative moment estimation (ADAM) minimizer of normalized mean square
residual over the manifold of variations of `A`.

!!! warning
    This function will set the (global) `TaylorSeries` variables to `dx₁ dx₂ dx₃ dx₄ dx₅ dx₆`.

!!! reference
    See Algorithm 1 of https://doi.org/10.48550/arXiv.1412.6980.
"""
function adam(radec::Vector{RadecMPC{T}}, A::AdmissibleRegion{T}, ρ::T, v_ρ::T,
              params::NEOParameters{T}; scale::Symbol = :linear, η::T = 25.0,
              μ::T = 0.75, ν::T = 0.9, ϵ::T = 1e-8, Qtol::T = 0.001, varorder::Int = 2,
              dynamics::D = newtonian!) where {T <: AbstractFloat, D}
    # Initial time of integration [julian days]
    jd0 = datetime2julian(A.date)
    # Maximum number of iterations
    maxiter = params.maxiter
    # Allocate memory
    aes = Matrix{T}(undef, 6, maxiter+1)
    Qs = fill(T(Inf), maxiter+1)
    # Initial attributable elements
    aes[:, 1] .= [A.α, A.δ, A.v_α, A.v_δ, ρ, v_ρ]
    # Scaling factors
    scalings = Vector{T}(undef, 6)
    scalings[1:4] .= abs.(aes[1:4, 1]) ./ 1e6
    if scale == :linear
        scalings[5] = (A.ρ_domain[2] - A.ρ_domain[1]) / 1_000
    elseif scale == :log
        scalings[5] = (log10(A.ρ_domain[2]) - log10(A.ρ_domain[1])) / 1_000
    end
    scalings[6] = (A.v_ρ_domain[2] - A.v_ρ_domain[1]) / 1_000
    # Jet transport variables
    dae = [scalings[i] * TaylorN(i, order = varorder) for i in 1:6]
    # Propagation buffer
    t0, tf = datetime2days(date(radec[1])), datetime2days(date(radec[end]))
    tlim = (t0 - params.bwdoffset, tf + params.fwdoffset)
    buffer = PropagationBuffer(dynamics, jd0, tlim, aes[:, 1] .+ dae, params)
    # Vector of O-C residuals
    res = Vector{OpticalResidual{T, TaylorN{T}}}(undef, length(radec))
    # Origin
    x0 = zeros(T, 6)
    x1 = zeros(T, 6)
    # Gradient of objective function wrt (ρ, v_ρ)
    g_t = Vector{T}(undef, 2)
    # First momentum
    m = zeros(T, 2)
    _m_ = zeros(T, 2)
    # Second momentum
    n = zeros(T, 2)
    _n_ = zeros(T, 2)
    # Gradient descent
    for t in 1:maxiter
        # Current attributable elements (plain)
        ae = aes[:, t]
        # Attributable elements (JT)
        if scale == :linear
            AE = ae + dae
        elseif scale == :log
            AE = [ae[1] + dae[1], ae[2] + dae[2], ae[3] + dae[3],
                  ae[4] + dae[4], 10^(log10(ae[5]) + dae[5]), ae[6] + dae[6]]
        end
        # Barycentric state vector
        q = attr2bary(A, AE, params)
        # Propagation and residuals
        # TO DO: `ρ::TaylorN` is too slow for `adam` due to evaluations
        # within the dynamical model
        propres!(res, radec, jd0 - ae[5]/c_au_per_day, q, params; buffer, dynamics)
        iszero(length(res)) && break
        # Least squares fit
        fit = tryls(res, x0, 5, 1:4)
        x1 .= fit.x
        # Current Q
        Q = nms(res)
        Q(x1) < 0 && break
        Qs[t] = Q(x1)
        # Convergence condition
        t > 1 && abs(Qs[t] - Qs[t-1]) / Qs[t] < Qtol && break
        # Gradient of objective function wrt (ρ, v_ρ)
        g_t[1] = differentiate(Q, 5)(x1)
        g_t[2] = differentiate(Q, 6)(x1)
        # First momentum
        m .= μ * m + (1 - μ) * g_t
        _m_ .= m / (1 - μ^t)
        # Second momentum
        n .= ν * n + (1 - ν) * g_t .^ 2
        _n_ .= n / (1 - ν^t)
        # Step
        x1[5:6] = x1[5:6] - η * _m_ ./ (sqrt.(_n_) .+ ϵ)
        # Update attributable elements
        aes[:, t+1] .= AE(x1)
        # Projection
        aes[5:6, t+1] .= boundary_projection(A, aes[5, t+1], aes[6, t+1])
    end
    # Find attributable elements with smallest Q
    t = argmin(Qs)

    return aes[:, t], Qs[t]
end

# Order in which to check tracklets in tooshortarc
function tsatrackletorder(x::Tracklet{T}, y::Tracklet{T}) where {T <: AbstractFloat}
    if x.nobs == y.nobs
        return x.date > y.date
    else
        return x.nobs > y.nobs
    end
end

# Point in manifold of variations -> NEOSolution{T, T}
function _tooshortarc(A::AdmissibleRegion{T}, radec::Vector{RadecMPC{T}},
                      tracklets::Vector{Tracklet{T}}, i::Int, params::NEOParameters{T};
                      scale::Symbol = :log, maxiter::Int = 5, dynamics::D = newtonian!) where {T <: Real, D}
    # Center
    if scale == :linear
        ρ = sum(A.ρ_domain) / 2
    elseif scale == :log
        ρ = A.ρ_domain[1]
    end
    v_ρ = sum(A.v_ρ_domain) / 2
    # ADAM minimization over manifold of variations
    ae, Q = adam(tracklets[i].radec, A, ρ, v_ρ, params; scale, dynamics)
    # ADAM failed to converge
    if isinf(Q)
        return zero(NEOSolution{T, T})
    # Jet Transport Least Squares
    else
        # Initial time of integration [julian days]
        # (corrected for light-time)
        jd0 = datetime2julian(A.date) - ae[5] / c_au_per_day
        # Convert attributable elements to barycentric cartesian coordinates
        q0 = attr2bary(A, ae, params)
        # Scaling factors
        scalings = abs.(q0) ./ 10^5
        # Jet Transport initial condition
        q = [q0[i] + scalings[i] * TaylorN(i, order = 6) for i in 1:6]
        # Jet Transport Least Squares
        return jtls(radec, tracklets, jd0, q, i, i, params; maxiter, dynamics)
    end
end

@doc raw"""
    tooshortarc(radec::Vector{RadecMPC{T}}, tracklets::Vector{Tracklet{T}},
                params::NEOParameters{T}; dynamics::D = newtonian!) where {T <: AbstractFloat, D}

Return initial conditions by minimizing the normalized root mean square residual
over the manifold of variations.

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
    # Declare jet transport variables
    _ = scaled_variables("dx", ones(T, 6); order = 6)

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