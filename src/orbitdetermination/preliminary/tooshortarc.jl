@doc raw"""
    mmov(od, i, A, ρ, v_ρ, params; kwargs...) where {D, T <: Real}

Minimization over the MOV method for preliminary orbit determination.

## Arguments

- `od::ODProblem{D, T}`: orbit determination problem.
- `i::Int`: index of reference tracklet.
- `A::AdmissibleRegion{T}`: admissible region of reference tracklet.
- `ρ/v_ρ::T`: starting point on `A`.
- `params::Parameters{T}`: see the `Minimization over the MOV` section of
    [`Parameters`](@ref).

## Keyword arguments

- `scale::Symbol`: horizontal scale, either `:log` (default) or `:linear`.
- `η::T`: learning rate (default: `25.0`).
- `μ::T`: first moment (default: `0.75`).
- `ν::T`: second moment (default: `0.9`).
- `ϵ::T`: numerical stability constant (default: `1e-8`).
- `adamorder::Int`: jet transport order (default: `2`).

!!! reference
    See Section 4 of https://doi.org/10.1007/s10569-025-10246-2.
"""
function mmov(od::ODProblem{D, T}, i::Int, A::AdmissibleRegion{T}, ρ::T, v_ρ::T,
    params::Parameters{T}; scale::Symbol = :linear, η::T = 25.0, μ::T = 0.75,
    ν::T = 0.9, ϵ::T = 1e-8, adamorder::Int = 2) where {D, T <: Real}
    # Unpack parameters
    @unpack adamiter, adammode, adamQtol, significance = params
    # Initial time of integration [julian days TDB]
    jd0 = dtutc2jdtdb(A.date)
    # Allocate memory
    orbits = [zero(MMOVOrbit{T, T}) for _ in 1:adamiter]
    aes = Matrix{T}(undef, 6, adamiter+1)
    Qs = fill(T(Inf), adamiter+1)
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
    # Considered tracklets
    tracklets = adammode ? od.tracklets[:] : od.tracklets[i:i]
    idxs = indices(tracklets)
    # Propagation buffer
    buffer = PropagationBuffer(od, jd0, idxs[1], idxs[end], AE, params)
    # Vector of O-C residuals
    res = init_residuals(TaylorN{T}, od, idxs)
    # Origin
    x0, x1 = zeros(T, 6), zeros(T, 6)
    # Least squares cache and methods
    lscache = LeastSquaresCache(x0, 1:4, 5)
    lsmethods = _lsmethods(res, x0, 1:4)
    # Gradient of objective function wrt (ρ, v_ρ)
    g_t = Vector{T}(undef, 2)
    # First and second momentum
    m, _m_ = zeros(T, 2), zeros(T, 2)
    n, _n_ = zeros(T, 2), zeros(T, 2)
    # Detect sawtooth efect (see https://arxiv.org/abs/2410.10056#)
    Qthreshold = nms_threshold(2*length(res), significance)
    Nsawtooth = 0
    # Gradient descent
    for t in 1:adamiter
        # Current attributable elements (plain)
        ae = view(aes, :, t)
        # Attributable elements (JT)
        if scale == :linear
            AE .= ae + dae
        elseif scale == :log
            AE .= ae[1] + dae[1], ae[2] + dae[2], ae[3] + dae[3],
                ae[4] + dae[4], 10^(log10(ae[5]) + dae[5]), ae[6] + dae[6]
        end
        # Barycentric state vector
        q = attr2bary(A, AE, params)
        # Propagation and residuals
        # TO DO: `ρ::TaylorN` is too slow for `adam` due to evaluations
        # within the dynamical model
        bwd, fwd = propres!(res, od, jd0 - ae[5]/c_au_per_day, q, params; buffer, idxs)
        iszero(length(res)) && break
        # Least squares fit
        fit = tryls(res, x0, lscache, lsmethods)
        !fit.success && break
        x1 .= fit.x
        # Current Q
        Q = nms(res)
        Q(x1) < 0 && break
        Qs[t] = Q(x1)
        # Covariance matrix
        nobs = 2 * notout(res)
        C = (nobs/2) * TS.hessian(Q, x1)
        Γ = inv(C)
        # Residuals space to barycentric coordinates jacobian
        J = Matrix(TS.jacobian(q - cte.(q), x1))
        # Update orbit
        orbits[t] = evaldeltas(MMOVOrbit(tracklets, bwd, fwd, res,
            Γ, J, aes[:, 1:t], Qs[1:t]), x0)
        # Convergence conditions
        if t > 1
            (Qs[t-1] < Qthreshold < Qs[t]) && (Nsawtooth += 1)
            Nsawtooth == 2 && break
            abs(Qs[t] - Qs[t-1]) / Qs[t] < adamQtol && break
        end
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
    # Find orbit with smallest Q
    t = argmin(Qs)

    return orbits[t]
end

@doc raw"""
    mmov(od, gorbit, i, params) where {D, T <: Real}

Refine a Gauss orbit via Minimization over the MOV.

- `od::ODProblem{D, T}`: orbit determination problem.
- `gorbit::AbstractOrbit{T, T}`: a priori orbit.
- `i::Int`: index of reference tracklet.
- `params::Parameters{T}`: see the `Gauss Method` and `Minimization over the MOV`
    sections of [`Parameters`](@ref).
"""
function mmov(od::ODProblem{D, T}, orbit::O, i::Int,
    params::Parameters{T}) where {D, T <: Real, O <: AbstractOrbit{T, T}}
    # Unpack parameters
    @unpack refscale = params
    # Middle tracklet
    tracklet = od.tracklets[i]
    # Admissible region
    A = AdmissibleRegion(tracklet, params)
    # Epoch [days since J2000]
    At0 = dtutc2days(A.date)
    # Barycentric cartesian initial condition
    q0 = orbit(At0)
    # Range and range rate
    ρ, v_ρ = bary2topo(A, q0)
    # Boundary projection
    ρ, v_ρ = boundary_projection(A, ρ, v_ρ)
    # Minimization over the MOV
    orbit1 = mmov(od, i, A, ρ, v_ρ, params; scale = refscale)

    return orbit1
end

@doc raw"""
    tsaiod(od, params; kwargs...) where {D, I, T <: Real}

Compute a `LeastSquaresOrbit` via Minimization over the MOV followed by
Jet Transport Least Squares.

See also [`mmov`](@ref).

## Arguments

- `od::ODProblem{D, T}`: an orbit determination problem.
- `params::Parameters{T}`: see the `Minimization over the MOV` and `Least
    Squares` sections  of [`Parameters`](@ref).

## Keyword arguments

- `initcond::I`: naive initial conditions function; takes as input an
    `AdmissibleRegion{T}` and outputs a `Vector{Tuple{T, T, Symbol}}`,
    where each element has the form `(ρ, v_ρ, scale)` (default: `iodinitcond`).

!!! warning
    This function may change the (global) `TaylorSeries` variables.

!!! reference
    See https://doi.org/10.1007/s10569-025-10246-2.
"""
function tsaiod(od::ODProblem{D, T}, params::Parameters{T};
    initcond::I = iodinitcond) where {D, I, T <: Real}
    # Allocate memory for orbit
    orbit = zero(LeastSquaresOrbit{T, T})
    # Unpack
    @unpack tsaorder, adammode, significance, verbose = params
    @unpack radec, tracklets = od
    # Set jet transport variables
    set_od_order(params)
    # Iterate tracklets
    for i in eachindex(tracklets)
        # Minimization over the MOV requires a minimum of 2 observations
        nobs(tracklets[i]) < 2 && continue
        # Admissible region
        A = AdmissibleRegion(tracklets[i], params)
        # List of naive initial conditions
        I0 = initcond(A)
        # Iterate naive initial conditions
        for j in eachindex(I0)
            # Minimization over the MOV
            ρ, v_ρ, scale = I0[j]
            porbit = mmov(od, i, A, ρ, v_ρ, params; scale)
            # Failed to converge
            iszero(porbit) && continue
            # Jet Transport Least Squares
            _orbit_ = jtls(od, porbit, params, adammode)
            # Update orbit
            orbit = updateorbit(orbit, _orbit_, radec)
            # Termination condition
            if critical_value(orbit) < significance
                N1, N2 = length(porbit.Qs), length(orbit.Qs)
                verbose && println(
                    "* Minimization over the MOV converged in $N1 iterations to:\n\n",
                    summary(porbit), "\n",
                    "* Jet Transport Least Squares converged in $N2 iterations to: \n\n",
                    summary(orbit)
                )
                return orbit
            end
        end
        # Global MMOV should be independent of starting tracklet
        adammode && break
    end
    # Unsuccessful orbit determination
    verbose && @warn("Orbit determination did not converge within \
        the given parameters or could not fit all the astrometry")

    return orbit
end