@doc raw"""
    adam(od::ODProblem{D, T}, i::Int, A::AdmissibleRegion{T}, ρ::T, v_ρ::T,
         params::NEOParameters{T}; kwargs...) where {D, T <: Real}

Adaptative moment estimation (ADAM) minimizer of normalized mean square
residual of `od.tracklets[i].radec` over the manifold of variations of `A`,
starting from `(ρ, v_ρ)`.

## Keyword arguments

- `scale::Symbol`: horizontal scale, either `:log` (default) or `:linear`.
- `η::T`: learning rate (default: `25.0`).
- `μ::T`: first moment (default: `0.75`).
- `ν::T`: second moment (default: `0.9`).
- `ϵ::T`: numerical stability constant (default: `1e-8`).
- `adamorder::Int`: jet transport order (default: `2`).

!!! reference
    See Algorithm 1 of https://doi.org/10.48550/arXiv.1412.6980.
"""
function adam(od::ODProblem{D, T}, i::Int, A::AdmissibleRegion{T}, ρ::T, v_ρ::T,
    params::NEOParameters{T}; scale::Symbol = :linear, η::T = 25.0,
    μ::T = 0.75, ν::T = 0.9, ϵ::T = 1e-8, adamorder::Int = 2) where {D, T <: Real}
    # Initial time of integration [julian days TDB]
    jd0 = dtutc2jdtdb(A.date)
    # Unfold maximum number of iterations and relative tolerance
    maxiter, Qtol = params.adamiter, params.adamQtol
    # Allocate memory
    aes = Matrix{T}(undef, 6, maxiter+1)
    Qs = fill(T(Inf), maxiter+1)
    # Initial attributable elements
    aes[:, 1] .= A.α, A.δ, A.v_α, A.v_δ, ρ, v_ρ
    # Scaling factors
    scalings = Vector{T}(undef, 6)
    scalings[1:4] .= abs.(aes[1:4, 1]) ./ 1e6
    if scale == :linear
        scalings[5] = (A.ρ_domain[2] - A.ρ_domain[1]) / 1_000
    elseif scale == :log
        scalings[5] = (log10(A.ρ_domain[2]) - log10(A.ρ_domain[1])) / 1_000
    end
    scalings[6] = (A.v_ρ_domain[2] - A.v_ρ_domain[1]) / 1_000
    # Jet transport variables and initial condition
    dae = [scalings[i] * TaylorN(i, order = adamorder) for i in 1:6]
    AE = aes[:, 1] .+ dae
    # Subset of radec
    idxs = indices(od.tracklets[i])
    # Propagation buffer
    buffer = PropagationBuffer(od, jd0, idxs[1], idxs[end], AE, params)
    # Vector of O-C residuals
    res = [zero(OpticalResidual{T, TaylorN{T}}) for _ in eachindex(idxs)]
    # Origin
    x0, x1 = zeros(T, 6), zeros(T, 6)
    # Least squares cache and methods
    cache = LeastSquaresCache(x0, 1:4, 5)
    methods = _lsmethods(res, x0, 1:4)
    # Gradient of objective function wrt (ρ, v_ρ)
    g_t = Vector{T}(undef, 2)
    # First momentum
    m, _m_ = zeros(T, 2), zeros(T, 2)
    # Second momentum
    n, _n_ = zeros(T, 2), zeros(T, 2)
    # Gradient descent
    for t in 1:maxiter
        # Current attributable elements (plain)
        ae = aes[:, t]
        # Attributable elements (JT)
        if scale == :linear
            AE .= ae + dae
        elseif scale == :log
            AE .= [ae[1] + dae[1], ae[2] + dae[2], ae[3] + dae[3],
                  ae[4] + dae[4], 10^(log10(ae[5]) + dae[5]), ae[6] + dae[6]]
        end
        # Barycentric state vector
        q = attr2bary(A, AE, params)
        # Propagation and residuals
        # TO DO: `ρ::TaylorN` is too slow for `adam` due to evaluations
        # within the dynamical model
        propres!(res, od, jd0 - ae[5]/c_au_per_day, q, params; buffer, idxs)
        iszero(length(res)) && break
        # Least squares fit
        fit = tryls(res, x0, cache, methods)
        !fit.success && break
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

@doc raw"""
    tsaiod(od::ODProblem{D, T}, params::NEOParameters{T};
        initcond::I = iodinitcond) where {D, I, T <: Real}

Fit a preliminary orbit to `od` via jet transport minimization
of the normalized mean square residual over the manifold of variations.

## Arguments

- `od::ODProblem{D, T}`: an orbit determination problem.
- `params::NEOParameters{T}`: see `Too Short Arc Parameters` of [`NEOParameters`](@ref).

## Keyword arguments

- `initcond::I`: naive initial conditions function; takes as input an
    `AdmissibleRegion{T}` and outputs a `Vector{Tuple{T, T, Symbol}}`,
    where each element has the form `(ρ, v_ρ, scale)`
    (default: `iodinitcond`).
"""
function tsaiod(od::ODProblem{D, T}, params::NEOParameters{T};
    initcond::I = iodinitcond) where {D, I, T <: Real}
    # Allocate memory for orbit
    sol = zero(NEOSolution{T, T})
    # Iterate tracklets
    for i in eachindex(od.tracklets)
        # ADAM requires a minimum of 2 observations
        od.tracklets[i].nobs < 2 && continue
        # Admissible region
        A = AdmissibleRegion(od.tracklets[i], params)
        # List of naive initial conditions
        I0 = initcond(A)
        # Iterate naive initial conditions
        for j in eachindex(I0)
            # ADAM minimization over manifold of variations
            ρ, v_ρ, scale = I0[j]
            ae, Q = adam(od, i, A, ρ, v_ρ, params; scale)
            # ADAM failed to converge
            isinf(Q) && continue
            # Initial time of integration [julian days]
            # (corrected for light-time)
            jd0 = dtutc2jdtdb(A.date) - ae[5] / c_au_per_day
            # Convert attributable elements to barycentric cartesian coordinates
            q0 = attr2bary(A, ae, params)
            # Scaling factors
            scalings = abs.(q0) ./ 10^5
            # Jet Transport initial condition
            q = [q0[k] + scalings[k] * TaylorN(k, order = params.tsaorder) for k in 1:6]
            # Jet Transport Least Squares
            _sol_ = jtls(od, jd0, q, od.tracklets[i:i], params, false)
            # Update solution
            sol = updatesol(sol, _sol_, od.radec)
            # Termination condition
            nrms(sol) <= params.tsaQmax && return sol
        end
    end

    return sol
end