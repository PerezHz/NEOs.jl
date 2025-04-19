include("osculating.jl")
include("leastsquares/methods.jl")
include("leastsquares/targetfunctions.jl")
include("leastsquares/outlierrejection.jl")
include("leastsquares/fit.jl")
include("odproblem.jl")
include("propres.jl")
include("abstractorbit/abstractorbit.jl")
include("abstractorbit/preliminaryorbit.jl")
include("abstractorbit/leastsquaresorbit.jl")
include("admissibleregion.jl")
include("tooshortarc.jl")
include("gaussinitcond.jl")

# Initial subset of radec for jtls
function _initialtracklets(trksa::AbstractVector{Tracklet{T}},
    trksb::AbstractVector{Tracklet{T}}) where {T <: Real}
    # Starting tracklets
    tin = sort!(intersect(trksa, trksb))
    # Remaining tracklets
    tout = sort!(setdiff(trksa, trksb))
    # Consistency check
    @assert sort!(union(tin, tout)) == trksa
    # Sort tout by absolute time to tin
    et = mean(@. dtutc2et(date(tin)))
    dts = @. abs(dtutc2et(date(tout)) - et)
    permute!(tout, sortperm(dts))
    # Starting observations
    rin = indices(tin)
    # jtls needs at least three observations
    while length(rin) < 3 && !isempty(tout)
        tracklet = popfirst!(tout)
        push!(tin, tracklet)
        sort!(tin)
        rin = vcat(rin, indices(tracklet))
        sort!(rin)
    end

    return tin, tout, rin
end

# Incrementally add observations to fit

# Add as much tracklets as possible per iteration
function addradec!(::Val{true}, rin::Vector{Int}, fit::LeastSquaresFit{T},
    tin::Vector{Tracklet{T}}, tout::Vector{Tracklet{T}},
    res::Vector{OpticalResidual{T, TaylorN{T}}}, x0::Vector{T},
    params::Parameters{T}) where {T <: Real}
    while !isempty(tout)
        extra = indices(tout[1])
        fit_new = tryls(res[rin ∪ extra], x0; maxiter = params.lsiter)
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
    params::Parameters{T}) where {T <: Real}
    if critical_value(view(res, rin), fit) < params.significance && !isempty(tout)
        extra = indices(tout[1])
        fit_new = tryls(res[rin ∪ extra], x0; maxiter = params.lsiter)
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

# Decide whether q is suitable for jtls
function isjtlsfit(od::ODProblem{D, T}, jd0::V, q::Vector{TaylorN{T}},
    params::Parameters{T}) where {D, T <: Real, V <: Number}
    # Unpack
    @unpack eph_ea = params
    @unpack tracklets = od
    # Try plain propagation and residuals
    bwd, fwd, res = propres(od, jd0, cte.(q), params)
    isempty(res) && return false
    # Check that orbit crosses all admissible regions
    eph(t) = bwdfwdeph(t, bwd, fwd, false, false)
    mask = BitVector(undef, length(tracklets))
    for i in eachindex(tracklets)
        A = AdmissibleRegion(tracklets[i], params)
        At0 = dtutc2days(A.date)
        q0 = eph(At0)
        ρ, v_ρ = bary2topo(A, q0)
        ρmin, ρmax = A.ρ_domain
        c1 = (isnan(tracklets[i].mag) ? R_EA : ρmin) ≤ ρ ≤ ρmax
        y_range = rangerates(A, ρ, :outer)
        c2 = (length(y_range) == 2) && (y_range[1] ≤ v_ρ ≤ y_range[2])
        mask[i] = c1 && c2
    end
    all(mask) || return false
    # Check that orbit stays out of Earth's radius
    ts = find_zeros(t -> rvelea(eph, params, t), bwd.t0 + bwd.t[end] + 0.007,
        fwd.t0 + fwd.t[end] - 0.007)
    ds = map(t -> euclid3D(eph(t) - eph_ea(t)), ts)
    mask = ds .> R_EA

    return all(mask)
end

@doc raw"""
    jtls(od::ODProblem{D, T}, jd0::V, q::Vector{TayloN{T}},
        tracklets::AbstractVector{Tracklet{T}}, params::Parameters{T}
        [, mode::Bool]) where {D, T <: Real, V <: Number}

Compute an orbit via Jet Transport Least Squares.

## Arguments

- `od::ODProblem{D, T}`: orbit determination problem.
- `jd0::V`: reference epoch [Julian days TDB].
- `q::Vector{TaylorN{T}}`: jet transport initial condition.
- `tracklets::AbstractVector{Tracklet{T}}`: initial tracklets for fit.
- `params::Parameters{T}`: see `Jet Transport Least Squares Parameters`
    of [`Parameters`](@ref).
- `mode::Bool`: `addradec!` mode (default: `true`).
"""
function jtls(od::ODProblem{D, T}, jd0::V, q::Vector{TaylorN{T}},
    tracklets::AbstractVector{Tracklet{T}}, params::Parameters{T},
    mode::Bool = true) where {D, T <: Real, V <: Number}
    # Unpack
    @unpack lsiter, jtlsiter, outrej, jtlsmask, χ2_rec, χ2_rej,
        fudge, max_per = params
    # Plain initial condition
    q0 = cte.(q)
    # JT tail
    dq = q - q0
    # Propagation buffer
    buffer = PropagationBuffer(od, jd0, 1, nobs(od), q, params)
    # Number of jet transport variables
    Npar = get_numvars()
    # Origin
    x0 = zeros(T, Npar)
    # Vector of O-C residuals
    res = init_residuals(TaylorN{T}, od)
    # Least squares cache and methods
    lscache = LeastSquaresCache(x0, 1:Npar, lsiter)
    lsmethods = _lsmethods(res, x0, 1:Npar)
    # Initial subset of radec for orbit fit
    tin, tout, rin = _initialtracklets(od.tracklets, tracklets)
    # Allocate memory for orbits
    sols = [zero(LeastSquaresOrbit{T, T}) for _ in 1:jtlsiter]
    Qs = fill(T(Inf), jtlsiter)
    # Outlier rejection
    if outrej
        orcache = OutlierRejectionCache(T, nobs(od))
        outs = zeros(Int, jtlsiter)
    end
    # Jet transport least squares
    for i in 1:jtlsiter
        # Initial conditions
        q = q0 + dq
        # Decide whether q is suitable for jtls
        if jtlsmask
            isjtlsfit(od, jd0, q, params) || break
        end
        # Propagation & residuals
        bwd, fwd = propres!(res, od, jd0, q, params; buffer)
        iszero(length(res)) && break
        # Orbit fit
        fit = tryls(res[rin], x0, lscache, lsmethods)
        !fit.success && break
        # Incrementally add observations to fit
        rin, fit = addradec!(Val(mode), rin, fit, tin, tout, res, x0, params)
        # Residuals space to barycentric coordinates jacobian
        J = Matrix(TS.jacobian(dq, fit.x))
        all(diag(J * fit.Γ * J') .> 0) || break
        # Outlier rejection
        if outrej
            outlier_rejection!(view(res, rin), fit.x, fit.Γ, orcache;
                χ2_rec, χ2_rej, fudge, max_per)
        end
        # Update solution
        Qs[i] = nrms(res, fit)
        sols[i] = evalfit(LeastSquaresOrbit(tin, bwd, fwd, res[rin], fit, J))
        if outrej
            outs[i] = notout(res)
        end
        # Convergence conditions
        if i > 1
            C1 = abs(Qs[i-1] - Qs[i]) < 0.01
            C2 = outrej ? (outs[i-1] == outs[i]) : true
            C1 && C2 && break
        end
        i > 2 && issorted(view(Qs, i-2:i)) && break
        # Update initial condition
        TS.evaluate!(q, fit.x, q0)
    end
    # Find complete solutions
    mask = findall(s -> length(s.res) == length(od.radec), sols)
    # Choose best solution
    if isempty(mask)
        _, k = findmin(nrms, sols)
    else
        _, k = findmin(nrms, view(sols, mask))
        k = mask[k]
    end

    return sols[k]
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
function updatesol(sol::LeastSquaresOrbit{T, T}, _sol_::LeastSquaresOrbit{T, T},
    radec::Vector{RadecMPC{T}}) where {T <: Real}
    if nobs(_sol_) == length(radec)
        return nrms(sol) <= nrms(_sol_) ? sol : _sol_
    else
        return sol
    end
end

@doc raw"""
    issinglearc(radec::Vector{RadecMPC{T}}, arc::Day = Day(30)) where {T <: Real}

Check whether `radec` is a single observational arc, i.e. no two consecutive observations
are more than `arc` days apart. The function assumes `radec` is sorted.
"""
function issinglearc(radec::Vector{RadecMPC{T}}, arc::Day = Day(30)) where {T <: Real}
    return all(diff(date.(radec)) .< arc)
end

@doc raw"""
    orbitdetermination(od::ODProblem{D, T}, params::Parameters{T};
        kwargs...) where {D, I, T <: Real}

Initial Orbit Determination (IOD) routine.

## Arguments

- `od::ODProblem{D, T}`: an orbit determination problem.
- `params::Parameters{T}`: see [`Parameters`](@ref).

## Keyword arguments

- `initcond::I`: naive initial conditions function; takes as input an
    `AdmissibleRegion{T}` and outputs a `Vector{Tuple{T, T, Symbol}}`,
    where each element has the form `(ρ, v_ρ, scale)`
    (default: `iodinitcond`).

!!! warning
    This function will change the (global) `TaylorSeries` variables.
"""
function orbitdetermination(od::ODProblem{D, T}, params::Parameters{T};
    initcond::I = iodinitcond) where {D, I, T <: Real}
    # Allocate memory for orbit
    sol = zero(LeastSquaresOrbit{T, T})
    # Unpack
    @unpack tsaorder, gaussorder, jtlsorder, significance = params
    @unpack radec = od
    # Cannot handle observatories without coordinates
    all(x -> hascoord(observatory(x)), radec) || return sol
    # Cannot handle zero observations or multiple arcs
    (isempty(radec) || !issinglearc(radec)) && return sol
    # Set jet transport variables
    varorder = max(tsaorder, gaussorder, jtlsorder)
    scaled_variables("dx", ones(T, 6); order = varorder)
    # Gauss method
    _sol_ = gaussiod(od, params)
    # Update solution
    sol = updatesol(sol, _sol_, radec)
    # Termination condition
    critical_value(sol) < significance && return sol
    # Too short arc
    _sol_ = tsaiod(od, params; initcond)
    # Update solution
    sol = updatesol(sol, _sol_, radec)

    return sol
end

@doc raw"""
    orbitdetermination(od::ODProblem{D, T}, sol::LeastSquaresOrbit{T, T},
        params::Parameters{T}) where {D, T <: Real}

Fit a least squares orbit to `od` using `sol` as an initial condition.

## Arguments

- `od::ODProblem{D, T}`: orbit determination problem.
- `sol::LeastSquaresOrbit{T, T}:` preliminary orbit.
- `params::Parameters{T}`: see [`Parameters`](@ref).

!!! warning
    This function will change the (global) `TaylorSeries` variables.
"""
function orbitdetermination(od::ODProblem{D, T}, sol::LeastSquaresOrbit{T, T},
    params::Parameters{T}) where {D, T <: Real}
    # Unpack parameters
    @unpack jtlsorder, significance, adammode = params
    @unpack radec, tracklets = od
    # Reference epoch [TDB]
    t = epoch(sol)
    jd0 = t + PE.J2000
    # Plain barycentric initial condition
    q0 = sol(t)
    # Scaling factors
    scalings = abs.(q0) ./ 10^6
    # Jet transport variables
    dq = scaled_variables("dx", scalings; order = jtlsorder)
    # Jet Transport initial condition
    q = q0 + dq
    # Jet Transport Least Squares
    sol1 = jtls(od, jd0, q, sol.tracklets, params, true)
    # Termination condition
    (nobs(sol1) == nobs(od) && critical_value(sol1) < significance) && return sol1
    # ADAM refinement
    _, i = findmin(@. abs(t - dtutc2days(date(tracklets))))
    jd0 = _adam!(od, i, q, jd0, params)
    # Jet Transport Least Squares
    trks = adammode ? tracklets[:] : tracklets[i:i]
    sol2 = jtls(od, jd0, q, trks, params, true)
    # Update solution
    sol1 = updatesol(sol1, sol2, radec)

    return sol1
end