include("osculating.jl")
include("least_squares.jl")
include("odproblem.jl")
include("neosolution.jl")
include("admissibleregion.jl")

# Times used within propres
function _proprestimes(radec::AbstractVector{RadecMPC{T}}, jd0::U,
    params::NEOParameters{T}) where {T <: Real, U <: Number}
    # Time of first (last) observation
    t0, tf = dtutc2jdtdb(date(radec[1])), dtutc2jdtdb(date(radec[end]))
    # TDB epoch (plain)
    _jd0_ = cte(cte(jd0))
    # Years in backward (forward) integration
    nyears_bwd = -(_jd0_ - t0 + params.bwdoffset) / yr
    nyears_fwd = (tf - _jd0_ + params.fwdoffset) / yr

    return t0, tf, _jd0_, nyears_bwd, nyears_fwd
end

# Propagate an orbit and compute residuals
function propres(od::ODProblem{D, T}, jd0::V, q0::Vector{U}, params::NEOParameters{T};
    buffer::Union{Nothing, PropagationBuffer{T, U, V}} = nothing,
    idxs::AbstractVector{Int} = eachindex(od.radec)) where {D, T <: Real, U <: Number,
    V <: Number}
    # Subset of radec for propagation and residuals
    radec = view(od.radec, idxs)
    # Times of first/last observation, epoch and years in backward/forward propagation
    t0, tf, _jd0_, nyears_bwd, nyears_fwd = _proprestimes(radec, jd0, params)
    # Propagation buffer
    if isnothing(buffer)
        tlim = (t0 - JD_J2000 - params.bwdoffset, tf - JD_J2000 + params.fwdoffset)
        buffer = PropagationBuffer(od.dynamics, jd0, tlim, q0, params)
    end
    # Backward (forward) integration
    bwd = _propagate(od.dynamics, jd0, nyears_bwd, q0, buffer, params)
    fwd = _propagate(od.dynamics, jd0, nyears_fwd, q0, buffer, params)
    if !issuccessfulprop(bwd, t0 - _jd0_; tol = params.coeffstol) ||
       !issuccessfulprop(fwd, tf - _jd0_; tol = params.coeffstol)
        return bwd, fwd, Vector{OpticalResidual{T, U}}(undef, 0)
    end
    # O-C residuals
    try
        w8s, bias = view(od.w8s.w8s, idxs), view(od.bias.bias, idxs)
        res = residuals(radec, w8s, bias;
            xvs = et -> auday2kmsec(params.eph_su(et/daysec)),
            xve = et -> auday2kmsec(params.eph_ea(et/daysec)),
            xva = et -> bwdfwdeph(et, bwd, fwd))
        return bwd, fwd, res
    catch
        return bwd, fwd, Vector{OpticalResidual{T, U}}(undef, 0)
    end
end

# In-place method of propres
function propres!(res::Vector{OpticalResidual{T, U}}, od::ODProblem{D, T},
    jd0::V, q0::Vector{U}, params::NEOParameters{T};
    buffer::Union{Nothing, PropagationBuffer{T, U, V}} = nothing,
    idxs::AbstractVector{Int} = eachindex(od.radec))  where {D, T <: Real, U <: Number,
    V <: Number}
    # Subset of radec for propagation and residuals
    radec = view(od.radec, idxs)
    # Times of first/last observation, epoch and years in backward/forward propagation
    t0, tf, _jd0_, nyears_bwd, nyears_fwd = _proprestimes(radec, jd0, params)
    # Propagation buffer
    if isnothing(buffer)
        tlim = (t0 - JD_J2000 - params.bwdoffset, tf - JD_J2000 + params.fwdoffset)
        buffer = PropagationBuffer(od.dynamics, jd0, tlim, q0, params)
    end
    # Backward (forward) integration
    bwd = _propagate(od.dynamics, jd0, nyears_bwd, q0, buffer, params)
    fwd = _propagate(od.dynamics, jd0, nyears_fwd, q0, buffer, params)
    if !issuccessfulprop(bwd, t0 - _jd0_; tol = params.coeffstol) ||
       !issuccessfulprop(fwd, tf - _jd0_; tol = params.coeffstol)
        empty!(res)
        return bwd, fwd
    end
    # O-C residuals
    try
        w8s, bias = view(od.w8s.w8s, idxs), view(od.bias.bias, idxs)
        residuals!(res, radec, w8s, bias;
            xvs = et -> auday2kmsec(params.eph_su(et/daysec)),
            xve = et -> auday2kmsec(params.eph_ea(et/daysec)),
            xva = et -> bwdfwdeph(et, bwd, fwd))
        return bwd, fwd
    catch
        empty!(res)
        return bwd, fwd
    end
end

# Initial subset of radec for jtls
function _initialtracklets(tracklets::Vector{Tracklet{T}}, i::Int) where {T <: Real}
    if iszero(i)
        tin = deepcopy(tracklets)
        tout = Vector{Tracklet{T}}(undef, 0)
        rin = sort!(reduce(vcat, indices.(tracklets)))
    else
        tin = [tracklets[i]]
        tout = sort(tracklets, by = x -> abs( (x.date - tracklets[i].date).value ))[2:end]
        rin = sort!(indices(tracklets[i]))
        while length(rin) < 3 && !isempty(tout)
            tracklet = popfirst!(tout)
            push!(tin, tracklet)
            sort!(tin)
            rin = vcat(rin, tracklet.indices)
            sort!(rin)
        end
    end
    return tin, tout, rin
end

# Incrementally add observations to fit

# Add as much tracklets as possible per iteration
function addradec!(::Val{true}, rin::Vector{Int}, fit::LeastSquaresFit{T},
    tin::Vector{Tracklet{T}}, tout::Vector{Tracklet{T}},
    res::Vector{OpticalResidual{T, TaylorN{T}}}, x0::Vector{T},
    params::NEOParameters{T}) where {T <: Real}
    while !isempty(tout)
        extra = indices(tout[1])
        fit_new = tryls(res[rin ∪ extra], x0, params.newtoniter)
        !fit_new.success && break
        fit = fit_new
        tracklet = popfirst!(tout)
        push!(tin, tracklet)
        sort!(tin)
        rin = vcat(rin, extra)
        sort!(rin)
    end

    return rin, fit
end

# Add at most one tracklet per iteration
function addradec!(::Val{false}, rin::Vector{Int}, fit::LeastSquaresFit{T},
    tin::Vector{Tracklet{T}}, tout::Vector{Tracklet{T}},
    res::Vector{OpticalResidual{T, TaylorN{T}}}, x0::Vector{T},
    params::NEOParameters{T}) where {T <: Real}
    if nrms(res[rin], fit) < params.tsaQmax && !isempty(tout)
        extra = indices(tout[1])
        fit_new = tryls(res[rin ∪ extra], x0, params.newtoniter)
        !fit_new.success && return rin, fit
        fit = fit_new
        tracklet = popfirst!(tout)
        push!(tin, tracklet)
        sort!(tin)
        rin = vcat(rin, extra)
        sort!(rin)
    end

    return rin, fit
end

@doc raw"""
    jtls(od::ODProblem{D, T}, jd0::V, q::Vector{TayloN{T}}, i::Int,
        params::NEOParameters{T} [, mode::Bool]) where {D, T <: Real, V <: Number}

Compute an orbit via Jet Transport Least Squares.

## Arguments

- `od::ODProblem{D, T}`: orbit determination problem.
- `jd0::V`: reference epoch [Julian days TDB].
- `q::Vector{TaylorN{T}}`: jet transport initial condition.
- `i::Int`: index of `tracklets` to start least squares fit.
- `params::NEOParameters{T}`: see `Jet Transport Least Squares Parameters`
    of [`NEOParameters`](@ref).
- `mode::Bool`: `addradec!` mode (default: `true`).
"""
function jtls(od::ODProblem{D, T}, jd0::V, q::Vector{TaylorN{T}}, i::Int,
    params::NEOParameters{T}, mode::Bool = true) where {D, T <: Real, V <: Number}
    # Plain initial condition
    q0 = constant_term.(q)
    # JT tail
    dq = q - q0
    # Vector of O-C residuals
    res = Vector{OpticalResidual{T, TaylorN{T}}}(undef, length(od.radec))
    # Propagation buffer
    t0, tf = dtutc2days(date(od.radec[1])), dtutc2days(date(od.radec[end]))
    tlim = (t0 - params.bwdoffset, tf + params.fwdoffset)
    buffer = PropagationBuffer(od.dynamics, jd0, tlim, q, params)
    # Origin
    x0 = zeros(T, 6)
    # Initial subset of radec for orbit fit
    tin, tout, rin = _initialtracklets(od.tracklets, i)
    # Residuals space to barycentric coordinates jacobian
    J = Matrix(TS.jacobian(dq))
    # Best orbit
    best_sol = zero(NEOSolution{T, T})
    best_Q = T(Inf)
    # Convergence flag
    flag = false
    # Jet transport least squares
    for _ in 1:params.jtlsiter
        # Initial conditions
        q = q0 + dq
        # Propagation & residuals
        bwd, fwd = propres!(res, od, jd0, q, params; buffer)
        iszero(length(res)) && break
        # Orbit fit
        fit = tryls(res[rin], x0, params.newtoniter)
        !fit.success && break
        # Incrementally add observations to fit
        rin, fit = addradec!(Val(mode), rin, fit, tin, tout, res, x0, params)
        # NRMS
        Q = nrms(res, fit)
        if abs(best_Q - Q) < 0.1
            flag = true
        end
        # Update NRMS and initial conditions
        best_Q = Q
        J .= TS.jacobian(dq, fit.x)
        best_sol = evalfit(NEOSolution(tin, bwd, fwd, res[rin], fit, J))
        flag && break
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

# Default naive initial conditions for iod
function iodinitcond(A::AdmissibleRegion{T}) where {T <: Real}
    v_ρ = sum(A.v_ρ_domain) / 2
    return [
        (A.ρ_domain[1], v_ρ, :log),
        (10^(sum(log10, A.ρ_domain) / 2), v_ρ, :log),
        (sum(A.ρ_domain) / 2, v_ρ, :log),
        (A.ρ_domain[2], v_ρ, :log),
    ]
end

# Update `sol` iff `_sol_` is complete and has a lower nrms
function updatesol(sol::NEOSolution{T, T}, _sol_::NEOSolution{T, T},
    radec::Vector{RadecMPC{T}}) where {T <: Real}
    if length(_sol_.res) == length(radec)
        return min(sol, _sol_)
    else
        return sol
    end
end

include("tooshortarc.jl")
include("gaussinitcond.jl")

@doc raw"""
    issinglearc(radec::Vector{RadecMPC{T}}, arc::Day = Day(30)) where {T <: Real}

Check whether `radec` is a single observational arc, i.e. no two consecutive observations
are more than `arc` days apart. The function assumes `radec` is sorted.
"""
function issinglearc(radec::Vector{RadecMPC{T}}, arc::Day = Day(30)) where {T <: Real}
    return all(diff(date.(radec)) .< arc)
end

@doc raw"""
    orbitdetermination(od::ODProblem{D, T}, params::NEOParameters{T};
        kwargs...) where {D, I, T <: Real}

Initial Orbit Determination (IOD) routine.

## Arguments

- `od::ODProblem{D, T}`: an orbit determination problem.
- `params::NEOParameters{T}`: see [`NEOParameters`](@ref).

## Keyword arguments

- `initcond::I`: naive initial conditions function; takes as input an
    `AdmissibleRegion{T}` and outputs a `Vector{Tuple{T, T, Symbol}}`,
    where each element has the form `(ρ, v_ρ, scale)`
    (default: `iodinitcond`).

!!! warning
    This function will change the (global) `TaylorSeries` variables.
"""
function orbitdetermination(od::ODProblem{D, T}, params::NEOParameters{T};
    initcond::I = iodinitcond) where {D, I, T <: Real}
    # Allocate memory for orbit
    sol = zero(NEOSolution{T, T})
    # Cannot handle observatories without coordinates
    all(x -> hascoord(observatory(x)), od.radec) || return sol
    # Cannot handle zero observations or multiple arcs
    (isempty(od.radec) || !issinglearc(od.radec)) && return sol
    # Set jet transport variables
    varorder = max(params.tsaorder, params.gaussorder, params.jtlsorder)
    scaled_variables("dx", ones(T, 6); order = varorder)
    # Gauss method
    _sol_ = gaussiod(od, params)
    # Update solution
    sol = updatesol(sol, _sol_, od.radec)
    # Termination condition
    nrms(sol) <= params.gaussQmax && return sol
    # Too short arc
    _sol_ = tsaiod(od, params; initcond)
    # Update solution
    sol = updatesol(sol, _sol_, od.radec)

    return sol
end

@doc raw"""
    orbitdetermination(od::ODProblem{D, T}, sol::NEOSolution{T, T},
        params::NEOParameters{T}) where {D, T <: Real}

Fit a least squares orbit to `od` using `sol` as an initial condition.

## Arguments

- `od::ODProblem{D, T}`: orbit determination problem.
- `sol::NEOSolution{T, T}:` preliminary orbit.
- `params::NEOParameters{T}`: see [`NEOParameters`](@ref).

!!! warning
    This function will change the (global) `TaylorSeries` variables.
"""
function orbitdetermination(od::ODProblem{D, T}, sol::NEOSolution{T, T},
    params::NEOParameters{T}) where {D, T <: Real}
    # Reference epoch [Julian days TDB]
    jd0 = sol.bwd.t0 + PE.J2000
    # Plain barycentric initial condition
    q0 = sol(sol.bwd.t0)
    # Scaling factors
    scalings = abs.(q0) ./ 10^6
    # Jet transport variables
    dq = scaled_variables("dx", scalings; order = params.jtlsorder)
    # Jet Transport initial condition
    q = q0 + dq
    # Jet Transport Least Squares
    _, i = findmin(t -> abs(dtutc2days(t.date) - sol.bwd.t0), od.tracklets)
    return jtls(od, jd0, q, i, params)
end