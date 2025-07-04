# Initial subset of optical astrometry for jtls
function _initialtracklets(trksa::AbstractTrackletVector{T},
                           trksb::AbstractTrackletVector{T}) where {T <: Real}
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
function addoptical!(::Val{true}, rin::Vector{Int}, fit::LeastSquaresFit{T},
                     tin::TrackletVector{T}, tout::TrackletVector{T},
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
function addoptical!(::Val{false}, rin::Vector{Int}, fit::LeastSquaresFit{T},
                     tin::TrackletVector{T}, tout::TrackletVector{T},
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
function isjtlsfit(od::AbstractODProblem{D, T}, jd0::Number, q::Vector{TaylorN{T}},
                   params::Parameters{T}) where {D, T <: Real}
    # Unpack
    @unpack eph_ea = params
    @unpack tracklets = od
    # Try plain propagation and residuals
    bwd, fwd, res = propres(od, cte.(q), jd0, params)
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

"""
    jtls(od, orbit, params [, mode::Bool])

Given an orbit determination problem `od`, return a `LeastSquaresOrbit`
computed via the Jet Transport Least Squares method, starting from a
preliminary `orbit`. For a list of parameters, see the `Least Squares`
section of [`Parameters`](@ref).

The `mode` optional argument determines how many tracklets are added
to the fit per iteration, either as much as possible (`true`, default)
or one (`false`).
"""
function jtls(od::IODProblem{D, T, O}, orbit::AbstractOrbit{D, T, T}, params::Parameters{T},
              mode::Bool = true) where {D, T <: Real, O <: AbstractOpticalVector{T}}
    # Unpack
    @unpack jtlsorder, lsiter, jtlsiter, outrej, jtlsmask, χ2_rec, χ2_rej,
            fudge, max_per = params
    # Reference epoch [Julian days TDB]
    jd0 = epoch(orbit) + PE.J2000
    # Allocate memory for orbits
    orbits = [zero(LeastSquaresOrbit{D, T, T, O, Nothing, Nothing}) for _ in 1:jtlsiter]
    q00s = Matrix{T}(undef, 6, jtlsiter+1)
    Qs = fill(T(Inf), jtlsiter+1)
    # Number of jet transport variables
    Npar = get_numvars()
    # Jet transport initial condition
    q00s[:, 1] .= orbit()
    scalings = all(>(0), variances(orbit)) ? sigmas(orbit) : abs.(q00s[:, 1]) ./ 1e6
    dq = [scalings[i] * TaylorN(i, order = jtlsorder) for i in 1:Npar]
    q0 = q00s[:, 1] + dq
    # Initialize propagation buffer and vector of residuals
    buffer = PropagationBuffer(od, q0, jd0, 1, nobs(od), params)
    res = init_optical_residuals(TaylorN{T}, od, orbit)
    # Origin
    x0 = zeros(T, Npar)
    # Least squares cache and methods
    lscache = LeastSquaresCache(x0, 1:Npar, lsiter)
    lsmethods = _lsmethods(res, x0, 1:Npar)
    # Initial subset of optical astrometry for orbit fit
    tin, tout, rin = _initialtracklets(od.tracklets, orbit.tracklets)
    # Outlier rejection
    if outrej
        orcache = OutlierRejectionCache(T, nobs(od))
        outs = zeros(Int, jtlsiter)
    end
    # Jet Transport Least Squares
    for i in 1:jtlsiter
        # Initial conditions
        @. q0 = q00s[:, i] + dq
        # Decide whether q0 is suitable for jtls
        if jtlsmask
            isjtlsfit(od, jd0, q0, params) || break
        end
        # Propagation & residuals
        bwd, fwd = propres!(res, od, q0, jd0, params; buffer)
        isempty(res) && break
        # Orbit fit
        fit = tryls(res[rin], x0, lscache, lsmethods)
        !fit.success && break
        # Incrementally add observations to fit
        rin, fit = addoptical!(Val(mode), rin, fit, tin, tout, res, x0, params)
        # Residuals space to barycentric coordinates jacobian
        jacobian = Matrix(TS.jacobian(dq, fit.x))
        all(>(0), diag(jacobian * fit.Γ * jacobian')) || break
        # Outlier rejection
        if outrej
            outlier_rejection!(view(res, rin), fit.x, fit.Γ, orcache;
                χ2_rec, χ2_rej, fudge, max_per)
        end
        # Update solution
        Qs[i] = nrms(res, fit)
        orbits[i] = evalfit(LeastSquaresOrbit(od.dynamics, od.optical[indices(tin)], tin,
            nothing, bwd, fwd, res[rin], nothing, fit, jacobian, q00s[:, 1:i], Qs[1:i]))
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
        q00s[:, i+1] .= q0(fit.x)
    end
    # Find complete solutions
    mask = findall(o -> noptical(o) == noptical(od), orbits)
    # Choose best solution
    if isempty(mask)
        _, k = findmin(nrms, orbits)
    else
        _, k = findmin(nrms, view(orbits, mask))
        k = mask[k]
    end

    return orbits[k]
end