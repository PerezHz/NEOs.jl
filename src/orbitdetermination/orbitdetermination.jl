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
include("preliminary/admissibleregion.jl")
include("preliminary/tooshortarc.jl")
include("preliminary/gaussinitcond.jl")

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
    jtls(od, orbit, params [, mode::Bool]) where {D, T <: Real}

Jet Transport Least Squares method for orbit determination.

## Arguments

- `od::ODProblem{D, T}`: orbit determination problem.
- `orbit::AbstractOrbit{T, T}`: a priori orbit.
- `params::Parameters{T}`: see the `Least Squares` section of [`Parameters`](@ref).
- `mode::Bool`: `addradec!` mode (default: `true`).
"""
function jtls(od::ODProblem{D, T}, orbit::O, params::Parameters{T},
              mode::Bool = true) where {D, T <: Real, O <: AbstractOrbit{T, T}}
    # Unpack
    @unpack jtlsorder, lsiter, jtlsiter, outrej, jtlsmask, χ2_rec, χ2_rej,
            fudge, max_per = params
    # Reference epoch [Julian days TDB]
    jd0 = epoch(orbit) + PE.J2000
    # Number of jet transport variables
    Npar = get_numvars()
    # Jet transport initial condition
    q00 = orbit()
    scalings = all(>(0), variances(orbit)) ? sigmas(orbit) : abs.(q00) ./ 1e6
    dq = [scalings[i] * TaylorN(i, order = jtlsorder) for i in 1:Npar]
    q0 = q00 + dq
    # Initialize propagation buffer and vector of residuals
    buffer = PropagationBuffer(od, jd0, 1, nobs(od), q0, params)
    res = init_residuals(TaylorN{T}, od, orbit)
    # Origin
    x0 = zeros(T, Npar)
    # Least squares cache and methods
    lscache = LeastSquaresCache(x0, 1:Npar, lsiter)
    lsmethods = _lsmethods(res, x0, 1:Npar)
    # Initial subset of radec for orbit fit
    tin, tout, rin = _initialtracklets(od.tracklets, orbit.tracklets)
    # Allocate memory for orbits
    orbits = [zero(LeastSquaresOrbit{T, T}) for _ in 1:jtlsiter]
    Qs = fill(T(Inf), jtlsiter)
    # Outlier rejection
    if outrej
        orcache = OutlierRejectionCache(T, nobs(od))
        outs = zeros(Int, jtlsiter)
    end
    # Jet Transport Least Squares
    for i in 1:jtlsiter
        # Initial conditions
        @. q0 = q00 + dq
        # Decide whether q0 is suitable for jtls
        if jtlsmask
            isjtlsfit(od, jd0, q0, params) || break
        end
        # Propagation & residuals
        bwd, fwd = propres!(res, od, jd0, q0, params; buffer)
        iszero(length(res)) && break
        # Orbit fit
        fit = tryls(res[rin], x0, lscache, lsmethods)
        !fit.success && break
        # Incrementally add observations to fit
        rin, fit = addradec!(Val(mode), rin, fit, tin, tout, res, x0, params)
        # Residuals space to barycentric coordinates jacobian
        J = Matrix(TS.jacobian(dq, fit.x))
        all(>(0), diag(J * fit.Γ * J')) || break
        # Outlier rejection
        if outrej
            outlier_rejection!(view(res, rin), fit.x, fit.Γ, orcache;
                χ2_rec, χ2_rej, fudge, max_per)
        end
        # Update solution
        Qs[i] = nrms(res, fit)
        orbits[i] = evalfit(LeastSquaresOrbit(tin, bwd, fwd, res[rin], fit, J))
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
        TS.evaluate!(q0, fit.x, q00)
    end
    # Find complete solutions
    mask = findall(o -> length(o.res) == length(od.radec), orbits)
    # Choose best solution
    if isempty(mask)
        _, k = findmin(nrms, orbits)
    else
        _, k = findmin(nrms, view(orbits, mask))
        k = mask[k]
    end

    return orbits[k]
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

# Update `orbit` iff `_orbit_` is complete and has a lower nrms
function updateorbit(orbit::LeastSquaresOrbit{T, T}, _orbit_::LeastSquaresOrbit{T, T},
    radec::Vector{RadecMPC{T}}) where {T <: Real}
    if nobs(_orbit_) == length(radec)
        return nrms(orbit) <= nrms(_orbit_) ? orbit : _orbit_
    else
        return orbit
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
    initialorbitdetermination(od, params; kwargs...) where {D, I, T <: Real}

Jet Transport Initial Orbit Determination.

## Arguments

- `od::ODProblem{D, T}`: an orbit determination problem.
- `params::Parameters{T}`: see [`Parameters`](@ref).

## Keyword arguments

- `initcond::I`: naive initial conditions function; takes as input an
    `AdmissibleRegion{T}` and outputs a `Vector{Tuple{T, T, Symbol}}`,
    where each element has the form `(ρ, v_ρ, scale)` (default: `iodinitcond`).

!!! warning
    This function will change the (global) `TaylorSeries` variables.

!!! reference
    See https://doi.org/10.1007/s10569-025-10246-2.
"""
function initialorbitdetermination(od::ODProblem{D, T}, params::Parameters{T};
    initcond::I = iodinitcond) where {D, I, T <: Real}
    # Allocate memory for orbit
    sol = zero(LeastSquaresOrbit{T, T})
    # Unpack
    @unpack tsaorder, gaussorder, jtlsorder, significance = params
    @unpack radec = od
    # Cannot handle observatories without coordinates
    all(x -> hascoord(observatory(x)), radec) || return sol
    # Cannot handle zero observations or multiple arcs
    # (isempty(radec) || !issinglearc(radec)) && return sol
    # Set jet transport variables
    varorder = max(tsaorder, gaussorder, jtlsorder)
    scaled_variables("dx", ones(T, 6); order = varorder)
    # Gauss method
    _sol_ = gaussiod(od, params)
    # Update solution
    sol = updateorbit(sol, _sol_, radec)
    # Termination condition
    critical_value(sol) < significance && return sol
    # Too short arc
    _sol_ = tsaiod(od, params; initcond)
    # Update solution
    sol = updateorbit(sol, _sol_, radec)

    return sol
end

@doc raw"""
    orbitdetermination(od, orbit, params) where {D, T <: Real}

Refine an existing orbit via Jet Transport Orbit Determination.

## Arguments

- `od::ODProblem{D, T}`: orbit determination problem.
- `orbit::LeastSquaresOrbit{T, T}:` a priori orbit.
- `params::Parameters{T}`: see [`Parameters`](@ref).

!!! warning
    This function will change the (global) `TaylorSeries` variables.

!!! reference
    See https://doi.org/10.1007/s10569-025-10246-2.
"""
function orbitdetermination(od::ODProblem{D, T}, orbit::LeastSquaresOrbit{T, T},
    params::Parameters{T}) where {D, T <: Real}
    # Unpack parameters
    @unpack significance = params
    @unpack radec, tracklets = od
    # Jet Transport Least Squares
    orbit1 = jtls(od, orbit, params, true)
    # Termination condition
    (nobs(orbit1) == nobs(od) && critical_value(orbit1) < significance) && return orbit1
    # Refine via minimization over the MOV
    j = closest_tracklet(epoch(orbit), tracklets)
    porbit = mmov(od, orbit, j, params)
    iszero(porbit) && return orbit1
    # Jet Transport Least Squares
    orbit2 = jtls(od, porbit, params, true)
    # Update orbit
    orbit1 = updateorbit(orbit1, orbit2, radec)

    return orbit1
end