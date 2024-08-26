@doc raw"""
    adam(radec::Vector{RadecMPC{T}}, A::AdmissibleRegion{T}, ρ::T, v_ρ::T,
         params::NEOParameters{T}; kwargs...) where {T <: Real, D}

Adaptative moment estimation (ADAM) minimizer of normalized mean square
residual of `radec` over the manifold of variations of `A`, starting
from `(ρ, v_ρ)`.

## Keyword arguments

- `scale::Symbol`: horizontal scale, either `:log` (default) or `:linear`.
- `η::T`: learning rate (default: `25.0`).
- `μ::T`: first moment (default: `0.75`).
- `ν::T`: second moment (default: `0.9`).
- `ϵ::T`: numerical stability constant (default: `1e-8`).
- `adamorder::Int`: jet transport order (default: `2`).
- `dynamics::D`: dynamical model (default: `newtonian!`).

!!! reference
    See Algorithm 1 of https://doi.org/10.48550/arXiv.1412.6980.
"""
function adam(radec::Vector{RadecMPC{T}}, A::AdmissibleRegion{T}, ρ::T, v_ρ::T,
              params::NEOParameters{T}; scale::Symbol = :linear, η::T = 25.0,
              μ::T = 0.75, ν::T = 0.9, ϵ::T = 1e-8, adamorder::Int = 2,
              dynamics::D = newtonian!) where {T <: Real, D}
    # Initial time of integration [julian days]
    jd0 = datetime2julian(A.date)
    # Maximum number of iterations
    maxiter = params.adamiter
    # Target function relative tolerance
    Qtol = params.adamQtol
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
    dae = [scalings[i] * TaylorN(i, order = adamorder) for i in 1:6]
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
function tsatrackletorder(x::Tracklet{T}, y::Tracklet{T}) where {T <: Real}
    if x.nobs == y.nobs
        return x.date > y.date
    else
        return x.nobs > y.nobs
    end
end

# Point in manifold of variations -> NEOSolution{T, T}
function _tooshortarc(A::AdmissibleRegion{T}, radec::Vector{RadecMPC{T}},
                      tracklets::Vector{Tracklet{T}}, i::Int, params::NEOParameters{T};
                      scale::Symbol = :log, dynamics::D = newtonian!) where {T <: Real, D}
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
        q = [q0[i] + scalings[i] * TaylorN(i, order = params.tsaorder) for i in 1:6]
        # Jet Transport Least Squares
        return jtls(radec, tracklets, jd0, q, i, params; dynamics)
    end
end

@doc raw"""
    tooshortarc(radec::Vector{RadecMPC{T}}, [tracklets::Vector{Tracklet{T}},]
                params::NEOParameters{T}; dynamics::D = newtonian!) where {T <: Real, D}

Compute an orbit by minimizing the normalized root mean square residual
over the manifold of variations.

# Arguments

- `radec::Vector{RadecMPC{T}}`: vector of optical astrometry.
- `tracklets::Vector{Tracklet{T}},`: vector of tracklets.
- `params::NEOParameters{T}`: see `Too Short Arc Parameters` of [`NEOParameters`](@ref).
- `dynamics::D`: dynamical model.

!!! warning
    This function will change the (global) `TaylorSeries` variables.
"""
function tooshortarc(radec::Vector{RadecMPC{T}}, params::NEOParameters{T};
                     dynamics::D = newtonian!) where {T <: Real, D}
    # Reduce tracklets by polynomial regression
    tracklets = reduce_tracklets(radec)
    # Set jet transport variables
    varorder = max(params.tsaorder, params.gaussorder)
    scaled_variables("dx", ones(T, 6); order = varorder)
    # Too Short Arc (TSA)
    return tooshortarc(radec, tracklets, params; dynamics)
end

function tooshortarc(radec::Vector{RadecMPC{T}}, tracklets::Vector{Tracklet{T}},
                     params::NEOParameters{T}; dynamics::D = newtonian!) where {T <: Real, D}

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