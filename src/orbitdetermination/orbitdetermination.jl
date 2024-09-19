include("osculating.jl")
include("least_squares.jl")
include("neosolution.jl")
include("admissibleregion.jl")

# Times used within propres
function _proprestimes(radec::Vector{RadecMPC{T}}, jd0::U, params::NEOParameters{T}) where {T <: Real, U <: Number}
    # Time of first (last) observation
    t0, tf = datetime2julian(date(radec[1])), datetime2julian(date(radec[end]))
    # TDB epoch (plain)
    _jd0_ = cte(cte(jd0))
    # Years in backward (forward) integration
    nyears_bwd = -(_jd0_ - t0 + params.bwdoffset) / yr
    nyears_fwd = (tf - _jd0_ + params.fwdoffset) / yr

    return t0, tf, _jd0_, nyears_bwd, nyears_fwd
end

# Propagate an orbit and compute residuals
function propres(radec::Vector{RadecMPC{T}}, jd0::V, q0::Vector{U}, params::NEOParameters{T};
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
    try
        res = residuals(radec, params;
            xvs = et -> auday2kmsec(params.eph_su(et/daysec)),
            xve = et -> auday2kmsec(params.eph_ea(et/daysec)),
            xva = et -> bwdfwdeph(et, bwd, fwd))
        return bwd, fwd, res
    catch
        return bwd, fwd, Vector{OpticalResidual{T, U}}(undef, 0)
    end
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
        empty!(res)
        return bwd, fwd
    end
    # O-C residuals
    try
        residuals!(res, radec, params;
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
function _initialtracklets(radec::Vector{RadecMPC{T}}, tracklets::Vector{Tracklet{T}},
    i::Int) where {T <: Real}
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
    jtls(radec::Vector{RadecMPC{T}}, tracklets::Vector{Tracklet{T}}, jd0::V,
         q::Vector{TayloN{T}}, i::Int, params::NEOParameters{T}, mode::Bool;
         kwargs...) where {D, T <: Real, V <: Number}

Compute an orbit via Jet Transport Least Squares.

## Arguments

- `radec::Vector{RadecMPC{T}}`: vector of optical astrometry.
- `tracklets::Vector{Tracklet{T}},`: vector of tracklets.
- `jd0::V`: reference epoch [Julian days TDB].
- `q::Vector{TaylorN{T}}`: jet transport initial condition.
- `i::Int`: index of `tracklets` to start least squares fit.
- `params::NEOParameters{T}`: see `Jet Transport Least Squares Parameters`
    of [`NEOParameters`](@ref).
- `mode::Bool`: `addradec!` mode (default: `true`).

## Keyword arguments

- `dynamics::D`: dynamical model (default: `newtonian!`).
"""
function jtls(radec::Vector{RadecMPC{T}}, tracklets::Vector{Tracklet{T}},
    jd0::V, q::Vector{TaylorN{T}}, i::Int, params::NEOParameters{T},
    mode::Bool = true; dynamics::D = newtonian!) where {D, T <: Real, V <: Number}
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
    # Initial subset of radec for orbit fit
    tin, tout, rin = _initialtracklets(radec, tracklets, i)
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
        bwd, fwd = propres!(res, radec, jd0, q, params; buffer, dynamics)
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
    isgauss(sol::NEOSolution{T, T}) where {T <: Real}

Check whether `sol` was computed via [`gaussinitcond`](@ref) (`true`) or
via [`tooshortarc`](@ref) (`false`).
"""
function isgauss(tracklets::Vector{Tracklet{T}}) where {T <: Real}
    # Observing stations
    obs = observatory.(tracklets)
    # TSA is not well suited for satellite observatories
    any(issatellite.(obs)) && return true
    # Gauss cannot handle less than 3 tracklets
    length(tracklets) < 3 && return false
    # Time span
    Δ = numberofdays(tracklets)
    # Gauss approximation does not work with less than 1 day
    return Δ > 1
end

isgauss(sol::NEOSolution{T, T}) where {T <: Real} = isgauss(sol.tracklets)

@doc raw"""
    orbitdetermination(radec::Vector{RadecMPC{T}}, params::NEOParameters{T};
                       dynamics::D = newtonian!) where {T <: Real, D}

Initial Orbit Determination (IOD) routine.

# Arguments

- `radec::Vector{RadecMPC{T}}`: vector of optical astrometry.
- `params::NEOParameters{T}`: see [`NEOParameters`](@ref).
- `dynamics::D`: dynamical model.

!!! warning
    This function will change the (global) `TaylorSeries` variables.
"""
function orbitdetermination(radec::Vector{RadecMPC{T}}, params::NEOParameters{T};
                            dynamics::D = newtonian!) where {T <: Real, D}
    # Allocate memory for orbit
    sol = zero(NEOSolution{T, T})
    # Eliminate observatories without coordinates
    filter!(x -> hascoord(observatory(x)), radec)
    # Cannot handle zero observations or multiple arcs
    if iszero(length(radec)) || !issinglearc(radec)
        return sol
    end
    # Reduce tracklets by polynomial regression
    tracklets = reduce_tracklets(radec)
    # Set jet transport variables
    varorder = max(params.tsaorder, params.gaussorder)
    scaled_variables("dx", ones(T, 6); order = varorder)
    # Case 1: Gauss' Method
    if isgauss(tracklets)
        sol = gaussinitcond(radec, tracklets, params; dynamics)
    end
    # Case 2: Too Short Arc (TSA)
    if iszero(sol) || nrms(sol) > params.gaussQmax
        sol = tooshortarc(radec, tracklets, params; dynamics)
    end

    return sol
end

@doc raw"""
    orbitdetermination(radec::Vector{RadecMPC{T}}, sol::NEOSolution{T, T}, params::NEOParameters{T};
                       dynamics::D = newtonian!) where {T <: Real, D}

Fit `sol` to `radec` via Jet Transport Least Squares.

# Arguments

- `radec::Vector{RadecMPC{T}}`: vector of optical astrometry.
- `sol::NEOSolution{T, T}:` preliminary orbit.
- `params::NEOParameters{T}`: see [`NEOParameters`](@ref).
- `dynamics::D`: dynamical model.

!!! warning
    This function will change the (global) `TaylorSeries` variables.
"""
function orbitdetermination(radec::Vector{RadecMPC{T}}, sol::NEOSolution{T, T}, params::NEOParameters{T};
                            dynamics::D = newtonian!) where {T <: Real, D}
    # Reduce tracklets by polynomial regression
    tracklets = reduce_tracklets(radec)
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
    _, i = findmin(tracklet -> abs(datetime2days(tracklet.date) - sol.bwd.t0), tracklets)
    return jtls(radec, tracklets, jd0, q, i, params; dynamics)
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

# Gauss method initial orbit determination
function _gaussiod(radec::Vector{RadecMPC{T}}, tracklets::Vector{Tracklet{T}},
    params::NEOParameters{T}; dynamics::D = newtonian!) where {T <: Real, D}
    # Allocate memory for orbit
    sol = zero(NEOSolution{T, T})
    # This function requires exactly 3 tracklets
    length(tracklets) != 3 && return sol
    # Jet transport scaling factors
    scalings = fill(1e-6, 6)
    # Jet transport perturbation (ra/dec)
    dq = [scalings[i] * TaylorN(i, order = params.gaussorder) for i in 1:6]
    # gauss_method input
    observatories = Vector{ObservatoryMPC{T}}(undef, 3)
    dates = Vector{DateTime}(undef, 3)
    α = Vector{T}(undef, 3)
    δ = Vector{T}(undef, 3)
    # Find best triplet of observations
    _gausstriplets2!(observatories, dates, α, δ, tracklets)
    # Julian day of middle observation
    _jd0_ = datetime2julian(dates[2])
    # Gauss method solution
    solG = gauss_method(observatories, dates, α .+ dq[1:3], δ .+ dq[4:6], params)
    # Filter non-physical (negative) rho solutions
    filter!(x -> all(cte.(x.ρ) .> 0), solG)
    # Filter Gauss solutions by heliocentric energy
    filter!(x -> cte(cte(heliocentric_energy(x.statevect))) <= 0, solG)
    # Iterate over Gauss solutions
    for i in eachindex(solG)
        # Epoch (corrected for light-time)
        jd0 = _jd0_ - cte(solG[i].ρ[2]) / c_au_per_day
        # Jet transport initial conditions
        q = solG[i].statevect .+ params.eph_su(jd0 - JD_J2000)
        # Jet Transport Least Squares
        _sol_ = jtls(radec, tracklets, jd0, q, 0, params, true; dynamics)
        # Update solution
        sol = updatesol(sol, _sol_, radec)
        # Termination condition
        nrms(sol) <= params.gaussQmax && return sol
        # ADAM help
        tracklets[2].nobs < 2 && continue
        jd0 = _adam!(q, jd0, tracklets[2], params; dynamics)
        # Jet Transport Least Squares
        _sol_ = jtls(radec, tracklets, jd0, q, 0, params, true; dynamics)
        # Update solution
        sol = updatesol(sol, _sol_, radec)
        # Termination condition
        nrms(sol) <= params.gaussQmax && return sol
    end

    return sol
end

# Too short arc initial orbit determination
function _tsaiod(radec::Vector{RadecMPC{T}}, tracklets::Vector{Tracklet{T}},
    params::NEOParameters{T}; dynamics::D1 = newtonian!,
    initcond::D2 = iodinitcond) where {T <: Real, D1, D2}
    # Allocate memory for orbit
    sol = zero(NEOSolution{T, T})
    # Iterate tracklets
    for i in eachindex(tracklets)
        # ADAM requires a minimum of 2 observations
        tracklets[i].nobs < 2 && continue
        # Admissible region
        A = AdmissibleRegion(tracklets[i], params)
        # List of naive initial conditions
        I0 = initcond(A)
        # Iterate naive initial conditions
        for j in eachindex(I0)
            # ADAM minimization over manifold of variations
            ρ, v_ρ, scale = I0[j]
            ae, Q = adam(tracklets[i].radec, A, ρ, v_ρ, params; scale, dynamics)
            # ADAM failed to converge
            isinf(Q) && continue
            # Initial time of integration [julian days]
            # (corrected for light-time)
            jd0 = datetime2julian(A.date) - ae[5] / c_au_per_day
            # Convert attributable elements to barycentric cartesian coordinates
            q0 = attr2bary(A, ae, params)
            # Scaling factors
            scalings = abs.(q0) ./ 10^5
            # Jet Transport initial condition
            q = [q0[k] + scalings[k] * TaylorN(k, order = params.tsaorder) for k in 1:6]
            # Jet Transport Least Squares
            _sol_ = jtls(radec, tracklets, jd0, q, i, params, false; dynamics)
            # Update solution
            sol = updatesol(sol, _sol_, radec)
            # Termination condition
            nrms(sol) <= params.tsaQmax && return sol
        end
    end

    return sol
end

@doc raw"""
    iod(radec::Vector{RadecMPC{T}}, params::NEOParameters{T};
        kwargs...) where {T <: Real, D1, D2}

Initial Orbit Determination (IOD) routine.

## Arguments

- `radec::Vector{RadecMPC{T}}`: vector of optical astrometry.
- `params::NEOParameters{T}`: see [`NEOParameters`](@ref).

## Keyword arguments

- `dynamics::D1`: dynamical model (default: `newtonian!`).
- `initcond::D2`: naive initial conditions function; takes as input an
    `AdmissibleRegion` and outputs a `Vector{Tuple{T, T, Symbol}}`,
    where each element has the form `(ρ, v_ρ, scale)`
    (default: `iodinitcond`).

!!! warning
    This function will change the (global) `TaylorSeries` variables.
"""
function iod(radec::Vector{RadecMPC{T}}, params::NEOParameters{T};
    dynamics::D1 = newtonian!, initcond::D2 = iodinitcond) where {T <: Real, D1, D2}
    # Allocate memory for orbit
    sol = zero(NEOSolution{T, T})
    # Eliminate observatories without coordinates
    filter!(x -> hascoord(observatory(x)), radec)
    # Cannot handle zero observations or multiple arcs
    (iszero(length(radec)) || !issinglearc(radec)) && return sol
    # Reduce tracklets by polynomial regression
    tracklets = reduce_tracklets(radec)
    # Set jet transport variables
    varorder = max(params.tsaorder, params.gaussorder, params.jtlsorder)
    scaled_variables("dx", ones(T, 6); order = varorder)
    # Gauss method
    _sol_ = _gaussiod(radec, tracklets, params; dynamics)
    # Update solution
    sol = updatesol(sol, _sol_, radec)
    # Termination condition
    nrms(sol) <= params.gaussQmax && return sol
    # Too short arc
    _sol_ = _tsaiod(radec, tracklets, params; dynamics, initcond)
    # Update solution
    sol = updatesol(sol, _sol_, radec)

    return sol
end