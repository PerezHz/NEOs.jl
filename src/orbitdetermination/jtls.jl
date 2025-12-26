# Initial subset of optical astrometry for jtls
function _initialtracklets(trksa::AbstractTrackletVector{T},
                           trksb::AbstractTrackletVector{T}) where {T <: Real}
    # Starting tracklets
    trksin = sort!(intersect(trksa, trksb))
    # Remaining tracklets
    trksout = sort!(setdiff(trksa, trksb))
    # Consistency check
    @assert sort!(union(trksin, trksout)) == trksa
    # jtls needs at least three observations
    while nobs(trksin) < 3 && !isempty(trksout)
        tracklet = popfirst!(trksout)
        push!(trksin, tracklet)
        sort!(trksin)
    end
    # Starting indices
    oidxs = indices(trksin)
    # Sort trksout by absolute time to trksin
    et = mean(@. dtutc2et(date(trksin)))
    dts = @. abs(dtutc2et(date(trksout)) - et)
    permute!(trksout, sortperm(dts))

    return trksin, trksout, oidxs
end

# Initial subset of radar astrometry for jtls
function _initialradar(od::MixedODProblem, orbit::AbstractOrbit)
    # Get radar astrometry
    radara = od.radar
    radarb = hasradar(orbit) ? orbit.radar : od.radar # empty(radara)
    # Starting radar observations
    radarin = sort!(intersect(radara, radarb))
    # Remaining radar observations
    radarout = sort!(setdiff(radara, radarb))
    # Consistency check
    @assert sort!(union(radarin, radarout)) == radara
    # Starting indices
    ridxs::Vector{Int} = indexin(radarin, od.radar)
    # Sort radarout by absolute time to radarin
    et = mean(@. dtutc2et(date(radarin)))
    dts = @. abs(dtutc2et(date(radarout)) - et)
    permute!(radarout, sortperm(dts))

    return radarin, radarout, ridxs
end

# Decide whether q is suitable for jtls
function isjtlsfit(od::OpticalODProblem, q::Vector{<:Number},
                   jd0::Number, params::Parameters)
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

# Incrementally add observations to fit

# Add as much tracklets as possible per iteration
function addoptical!(::Val{true}, oidxs::Vector{Int}, fit::LeastSquaresFit{T},
                     lscache::LeastSquaresCache{T},
                     lsmethods::Vector{AbstractLeastSquaresMethod{T}},
                     trksin::TrackletVector{T}, trksout::TrackletVector{T},
                     res::Vector{OpticalResidual{T, TaylorN{T}}},
                     x0::Vector{T}, params::Parameters{T}) where {T <: Real}
    while !isempty(trksout)
        extra = indices(trksout[1])
        fit_new = tryls(view(res, oidxs ∪ extra), x0, lscache, lsmethods)
        !fit_new.success && break
        fit = fit_new
        tracklet = popfirst!(trksout)
        push!(trksin, tracklet)
        sort!(trksin)
        oidxs = vcat(oidxs, extra)
        sort!(oidxs)
    end

    return oidxs, fit
end

# Add at most one tracklet per iteration
function addoptical!(::Val{false}, oidxs::Vector{Int}, fit::LeastSquaresFit{T},
                     lscache::LeastSquaresCache{T},
                     lsmethods::Vector{AbstractLeastSquaresMethod{T}},
                     trksin::TrackletVector{T}, trksout::TrackletVector{T},
                     res::Vector{OpticalResidual{T, TaylorN{T}}},
                     x0::Vector{T}, params::Parameters{T}) where {T <: Real}
    if critical_value(view(res, oidxs), fit) < params.significance && !isempty(trksout)
        extra = indices(trksout[1])
        fit_new = tryls(view(res, oidxs ∪ extra), x0, lscache, lsmethods)
        !fit_new.success && return oidxs, fit
        fit = fit_new
        tracklet = popfirst!(trksout)
        push!(trksin, tracklet)
        sort!(trksin)
        oidxs = vcat(oidxs, extra)
        sort!(oidxs)
    end

    return oidxs, fit
end

function addobservations!(od::MixedODProblem, oidxs::Vector{Int}, ridxs::Vector{Int},
                          fit::LeastSquaresFit{T}, lscache::LeastSquaresCache{T},
                          lsmethods::Vector{AbstractLeastSquaresMethod{T}},
                          trksin::TrackletVector{T}, trksout::TrackletVector{T},
                          radarin::AbstractRadarVector{T}, radarout::AbstractRadarVector{T},
                          res::AbstractResidualSet{T, TaylorN{T}},
                          x0::Vector{T}, params::Parameters{T}) where {T <: Real}
    # Add optical astrometry
    while !isempty(trksout)
        extra = indices(trksout[1])
        fit_new = tryls((res[1][oidxs ∪ extra], res[2][ridxs]), x0, lscache, lsmethods)
        !fit_new.success && break
        fit = fit_new
        tracklet = popfirst!(trksout)
        push!(trksin, tracklet)
        sort!(trksin)
        oidxs = vcat(oidxs, extra)
        sort!(oidxs)
    end
    # Add radar astrometry
    while !isempty(radarout)
        extra = findfirst(==(radarout[1]), od.radar)
        isnothing(extra) && break
        fit_new = tryls((res[1][oidxs], res[2][ridxs ∪ extra]), x0, lscache, lsmethods)
        !fit_new.success && break
        fit = fit_new
        radar = popfirst!(radarout)
        push!(radarin, radar)
        sort!(radarin)
        ridxs = vcat(ridxs, extra)
        sort!(ridxs)
    end

    return oidxs, ridxs, fit
end

function updateorbit(orbit::AbstractOrbit{D, T, TaylorN{T}},
                     lsmethods::Vector{AbstractLeastSquaresMethod{T}},
                     params::Parameters{T}) where {D, T <: Real}
    # Unpack parameters
    @unpack jtlsproject, H_max = params
    # Projection onto the admissible region
    if jtlsproject
        # Reference epoch
        t0 = epoch(orbit)
        # Lest squares fit
        @unpack fit = orbit
        # Barycentric initial condition
        q0 = orbit()
        # Geocentric attributable elements
        attr = cartesian2attributable(q0(fit.x) - params.eph_ea(t0))
        α, δ, v_α, v_δ, ρ, v_ρ = attr
        h = H_max + 5 * log10(R_EA)
        # Admissible region
        A = AdmissibleRegion(days2dtutc(t0), α, δ, v_α, v_δ, h,
                             search_observatory_code("500"), params)
        # Projection onto the admissible region
        if !((ρ, v_ρ) in A)
            ρ, v_ρ = boundary_projection(A, ρ, v_ρ)
            # Barycentric initial condition
            q1 = attributable2cartesian(SVector{6, T}(α, δ, v_α, v_δ, ρ, v_ρ)) +
                                        params.eph_ea(t0)
            # New deltas and covariance matrix
            @. fit.x = (q1 - cte(q0)) / $diag(orbit.jacobian)
            for method in lsmethods
                if typeof(method) == fit.routine
                    fit.Γ .= inv(normalmatrix(method, fit.x))
                end
            end
        end
    end
    # Evaluate deltas
    return evalfit(orbit)
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
function jtls(
        od::OpticalODProblem{D, T, O}, orbit::AbstractOrbit,
        params::Parameters{T}, mode::Bool = true
    ) where {D, T <: Real, O <: AbstractOpticalVector{T}}
    # Unpack parameters
    @unpack jtlsorder, lsiter, jtlsiter, outrej, jtlsmask, χ2_rec, χ2_rej,
            fudge, max_per, marsden_scalings = params
    # Number of degrees of freedom
    Ndof = dof(od)
    # Set jet transport variables
    Npar = numvars(Val(od.dynamics), params)
    set_od_order(params, Npar)
    # Reference epoch [Julian days TDB]
    jd0 = epoch(orbit) + PE.J2000
    # Pre-allocate memory
    orbits = zeros(OpticalLeastSquaresOrbit{D, T, T, O}, jtlsiter)
    q00s = Matrix{T}(undef, Ndof, jtlsiter+1)
    Qs = fill(T(Inf), jtlsiter+1)
    # Jet transport initial condition
    q00s[:, 1] = initialcondition(orbit, Ndof, params)
    variables = collect(1:6)
    if Ndof > 6
        variables = vcat(variables, findall(!iszero, marsden_scalings) .+ 6)
    end
    dq = jtperturbation(orbit, variables, Ndof, jtlsorder, params)
    q0 = q00s[:, 1] + dq
    # Initialize buffer and set of residuals
    buffer = PropresBuffer(od, q0, jd0, params)
    res = init_residuals(TaylorN{T}, od, orbit)
    # Origin
    x0 = zeros(T, Npar)
    # Least squares cache and methods
    lscache = LeastSquaresCache(x0, 1:Npar, lsiter)
    lsmethods = _lsmethods(res, x0, 1:Npar)
    # Initial subset of optical astrometry for orbit fit
    trksin, trksout, oidxs = _initialtracklets(od.tracklets, orbit.tracklets)
    # Outlier rejection
    if outrej
        orcache = OutlierRejectionCache(T, noptical(od))
        outs = zeros(Int, jtlsiter)
    end
    # Jet Transport Least Squares
    for i in 1:jtlsiter
        # Initial conditions
        @. q0 = q00s[:, i] + dq
        # Decide whether q0 is suitable for jtls
        if jtlsmask
            isjtlsfit(od, q0, jd0, params) || break
        end
        # Propagation & residuals
        bwd, fwd = propres!(res, od, q0, jd0, params; buffer)
        isempty(res) && break
        # Orbit fit
        fit = tryls(view(res, oidxs), x0, lscache, lsmethods)
        !fit.success && break
        # Incrementally add observations to fit
        oidxs, fit = addoptical!(Val(mode), oidxs, fit, lscache, lsmethods,
                                 trksin, trksout, res, x0, params)
        # Residuals space to barycentric coordinates jacobian
        jacobian = Matrix(TS.jacobian(dq[variables], fit.x))
        all(>(0), diag(jacobian * fit.Γ * jacobian')) || break
        # Outlier rejection
        if outrej
            outlier_rejection!(view(res, oidxs), fit.x, fit.Γ, orcache;
                χ2_rec, χ2_rej, fudge, max_per)
        end
        # Update orbit
        orbits[i] = updateorbit(LeastSquaresOrbit(
            od.dynamics, variables, od.optical[oidxs], trksin, nothing, bwd, fwd,
            res[oidxs], nothing, fit, jacobian, q00s[variables, 1:i], Qs[1:i]
        ), lsmethods, params)
        Qs[i] = orbits[i].Qs[end] = nrms(res, orbits[i].fit)
        if outrej
            outs[i] = notout(orbits[i])
        end
        # Convergence conditions
        if i > 1
            C1 = abs(Qs[i-1] - Qs[i]) < 0.01
            C2 = outrej ? (outs[i-1] == outs[i]) : true
            C1 && C2 && break
        end
        i > 2 && issorted(view(Qs, i-2:i)) && break
        # Update initial condition
        q00s[:, i+1] .= q0(orbits[i].fit.x)
    end
    # Find complete solutions
    mask = findall(o -> nobs(o) == nobs(od), orbits)
    # Choose best solution
    if isempty(mask)
        _, k = findmin(nrms, orbits)
    else
        _, k = findmin(nrms, view(orbits, mask))
        k = mask[k]
    end

    return orbits[k]
end

function jtls(
        od::MixedODProblem{D, T, O, R}, orbit::AbstractOrbit,
        params::Parameters{T}
    ) where {D, T <: Real, O <: AbstractOpticalVector{T}, R <: AbstractRadarVector{T}}
    # Unpack parameters
    @unpack jtlsorder, lsiter, jtlsiter, outrej, jtlsmask, χ2_rec, χ2_rej,
            fudge, max_per, marsden_scalings = params
    # Number of degrees of freedom
    Ndof = dof(od)
    # Set jet transport variables
    Npar = numvars(Val(od.dynamics), params)
    set_od_order(params, Npar)
    # Reference epoch [Julian days TDB]
    jd0 = epoch(orbit) + PE.J2000
    # Pre-allocate memory
    orbits = zeros(MixedLeastSquaresOrbit{D, T, T, O, R}, jtlsiter)
    q00s = Matrix{T}(undef, Ndof, jtlsiter+1)
    Qs = fill(T(Inf), jtlsiter+1)
    # Jet transport initial condition
    q00s[:, 1] = initialcondition(orbit, Ndof, params)
    variables = collect(1:6)
    if Ndof > 6
        variables = vcat(variables, findall(!iszero, marsden_scalings) .+ 6)
    end
    dq = jtperturbation(orbit, variables, Ndof, jtlsorder, params)
    q0 = q00s[:, 1] + dq
    # Initialize buffer and set of residuals
    buffer = PropresBuffer(od, q0, jd0, params)
    res = init_residuals(TaylorN{T}, od, orbit)
    # Origin
    x0 = zeros(T, Npar)
    # Least squares cache and methods
    lscache = LeastSquaresCache(x0, 1:Npar, lsiter)
    lsmethods = _lsmethods(res, x0, 1:Npar)
    # Initial subset of optical astrometry for orbit fit
    trksin, trksout, oidxs = _initialtracklets(od.tracklets, orbit.tracklets)
    # Initial subset of radar astrometry for orbit fit
    radarin, radarout, ridxs = _initialradar(od, orbit)
    # Outlier rejection
    if outrej
        orcache = OutlierRejectionCache(T, noptical(od))
        outs = zeros(Int, jtlsiter)
    end
    # Jet Transport Least Squares
    for i in 1:jtlsiter
        # Initial conditions
        @. q0 = q00s[:, i] + dq
        # Decide whether q0 is suitable for jtls
        # if jtlsmask
        #     isjtlsfit(od, q0, jd0, params) || break
        # end
        # Propagation & residuals
        bwd, fwd = propres!(res, od, q0, jd0, params; buffer)
        any(isempty, res) && break
        # Orbit fit
        fit = tryls((res[1][oidxs], res[2][ridxs]), x0, lscache, lsmethods)
        !fit.success && break
        # Incrementally add observations to fit
        oidxs, ridxs, fit = addobservations!(od, oidxs, ridxs, fit, lscache, lsmethods,
            trksin, trksout, radarin, radarout, res, x0, params)
        # Residuals space to barycentric coordinates jacobian
        jacobian = Matrix(TS.jacobian(dq[variables], fit.x))
        all(>(0), diag(jacobian * fit.Γ * jacobian')) || break
        # Outlier rejection
        if outrej
            outlier_rejection!(view(res[1], oidxs), fit.x, fit.Γ, orcache;
                χ2_rec, χ2_rej, fudge, max_per)
        end
        # Update orbit
        orbits[i] = updateorbit(LeastSquaresOrbit(
            od.dynamics, variables, od.optical[oidxs], trksin, od.radar[ridxs], bwd, fwd,
            res[1][oidxs], res[2][ridxs], fit, jacobian, q00s[variables, 1:i], Qs[1:i]
        ), lsmethods, params)
        Qs[i] = orbits[i].Qs[end] = nrms(res, orbits[i].fit)
        if outrej
            outs[i] = notout(orbits[i])
        end
        # Convergence conditions
        if i > 1
            C1 = abs(Qs[i-1] - Qs[i]) < 0.01
            C2 = outrej ? (outs[i-1] == outs[i]) : true
            C1 && C2 && break
        end
        i > 2 && issorted(view(Qs, i-2:i)) && break
        # Update initial condition
        q00s[:, i+1] .= q0(orbits[i].fit.x)
    end
    # Find complete solutions
    mask = findall(o -> nobs(o) == nobs(od), orbits)
    # Choose best solution
    if isempty(mask)
        _, k = findmin(nrms, orbits)
    else
        _, k = findmin(nrms, view(orbits, mask))
        k = mask[k]
    end

    return orbits[k]
end