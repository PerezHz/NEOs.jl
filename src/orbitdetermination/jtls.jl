"""
    JTLSBuffer{
        D,
        T <: Real,
        R <: AbstractResidualSet{T, TaylorN{T}},
        O <: AbstractOrbitVector{D, T, T}
    } <: AbstractBuffer

Pre-allocated memory for [`jtls`](@ref).

# Fields

- `res::R`:set of astrometric residuals.
- `orbits::Vector{O}`: orbit history.
- `Qs::Vector{T}`: nrms history.
- `q00s::Matrix{T}`: initial conditions history.
- `outs::Union{Nothing, Vector{Int}}`: number of non rejected observations history.
- `prbuffer::PropresBuffer{T, TaylorN{T}, T}`: buffer for [`propres`](@ref).
- `lscache::LeastSquaresCache{T}`: buffer for [`leastsquares!`](@ref).
- `orcache::Union{Nothing, OutlierRejectionCache{T}}`:  buffer for
    [`outlier_rejection!`](@ref).
"""
struct JTLSBuffer{
        D,
        T <: Real,
        R <: AbstractResidualSet{T, TaylorN{T}},
        O <: AbstractOrbitVector{D, T, T}
    } <: AbstractBuffer
    res::R
    orbits::O
    Qs::Vector{T}
    q00s::Matrix{T}
    outs::Union{Nothing, Vector{Int}}
    prbuffer::PropresBuffer{T, TaylorN{T}, T}
    lscache::LeastSquaresCache{T}
    orcache::Union{Nothing, OutlierRejectionCache{T}}
end

# Abbreviation
const AbstractJTLSBuffer{D, T} = JTLSBuffer{D, T, R, O} where {R, O}

# Constructor
function JTLSBuffer(od::AbstractODProblem{D, T}, orbit::AbstractOrbit,
                    params::Parameters) where {D, T <: Real}
    # Unpack
    @unpack jtlsorder, jtlsiter, lsiter, outrej = params
    # Jet transport initial condition
    Ndof = dof(od)
    Npar = numvars(Val(od.dynamics), params)
    jd0 = epoch(orbit) + PE.J2000
    q0 = jtinitialcondition(od, orbit, jd0, params)
    # Memory allocation
    res = init_residuals(TaylorN{T}, od, orbit)
    O, R = Vector{opticaltype(od)}, Vector{radartype(od)}
    orbits = Vector{hasradar(od) ? MixedLeastSquaresOrbit{D, T, T, O, R} :
        OpticalLeastSquaresOrbit{D, T, T, O}}(undef, jtlsiter)
    Qs = Vector{T}(undef, jtlsiter + 1)
    q00s = Matrix{T}(undef, Ndof, jtlsiter+1)
    outs = outrej ? Vector{Int}(undef, jtlsiter) : nothing
    prbuffer = PropresBuffer(od, q0, jd0, params)
    lscache = LeastSquaresCache(zeros(T, Npar), 1:Npar, lsiter)
    orcache = OutlierRejectionCache(T, noptical(od))
    return JTLSBuffer{D, T, typeof(res), typeof(orbits)}(res, orbits, Qs, q00s, outs,
        prbuffer, lscache, orcache)
end

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
    ts = find_zeros(t -> rvelea(eph, params, t), lasttime(bwd) + 0.007,
        lasttime(fwd) - 0.007)
    ds = map(t -> euclid3D(eph(t) - eph_ea(t)), ts)
    mask = ds .> R_EA

    return all(mask)
end

# Incrementally add observations to fit

# Add as much tracklets as possible per iteration
function addoptical!(
        ::Val{true}, oidxs::Vector{Int}, fit::LeastSquaresFit{T},
        lscache::LeastSquaresCache{T}, lsmethods::Tuple,
        trksin::TrackletVector{T}, trksout::TrackletVector{T},
        res::Vector{OpticalResidual{T, TaylorN{T}}},
        x0::Vector{T}, params::Parameters{T}
    ) where {T <: Real}
    Qtol, Mtol = params.lsQtol, params.lsMtol
    while !isempty(trksout)
        extra = indices(trksout[1])
        fit_new = tryls(view(res, oidxs ∪ extra), x0, lscache, lsmethods; Qtol, Mtol)
        !issuccess(fit_new) && break
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
function addoptical!(
        ::Val{false}, oidxs::Vector{Int}, fit::LeastSquaresFit{T},
        lscache::LeastSquaresCache{T}, lsmethods::Tuple,
        trksin::TrackletVector{T}, trksout::TrackletVector{T},
        res::Vector{OpticalResidual{T, TaylorN{T}}},
        x0::Vector{T}, params::Parameters{T}
    ) where {T <: Real}
    Qtol, Mtol = params.lsQtol, params.lsMtol
    if critical_value(view(res, oidxs), fit) < params.significance && !isempty(trksout)
        extra = indices(trksout[1])
        fit_new = tryls(view(res, oidxs ∪ extra), x0, lscache, lsmethods; Qtol, Mtol)
        !issuccess(fit_new) && return oidxs, fit
        fit = fit_new
        tracklet = popfirst!(trksout)
        push!(trksin, tracklet)
        sort!(trksin)
        oidxs = vcat(oidxs, extra)
        sort!(oidxs)
    end

    return oidxs, fit
end

function addobservations!(
        od::MixedODProblem, oidxs::Vector{Int}, ridxs::Vector{Int},
        fit::LeastSquaresFit{T}, lscache::LeastSquaresCache{T}, lsmethods::Tuple,
        trksin::TrackletVector{T}, trksout::TrackletVector{T},
        radarin::AbstractRadarVector{T}, radarout::AbstractRadarVector{T},
        res::AbstractResidualSet{T, TaylorN{T}}, x0::Vector{T}, params::Parameters{T}
    ) where {T <: Real}
    Qtol, Mtol = params.lsQtol, params.lsMtol
    # Add optical astrometry
    while !isempty(trksout)
        extra = indices(trksout[1])
        fit_new = tryls((res[1][oidxs ∪ extra], res[2][ridxs]), x0, lscache, lsmethods;
                         Qtol, Mtol)
        !issuccess(fit_new) && break
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
        fit_new = tryls((res[1][oidxs], res[2][ridxs ∪ extra]), x0, lscache, lsmethods;
                        Qtol, Mtol)
        !issuccess(fit_new) && break
        fit = fit_new
        radar = popfirst!(radarout)
        push!(radarin, radar)
        sort!(radarin)
        ridxs = vcat(ridxs, extra)
        sort!(ridxs)
    end

    return oidxs, ridxs, fit
end

function updateorbit(orbit::AbstractOrbit{D, T, TaylorN{T}}, lsmethods::Tuple,
                     params::Parameters{T}) where {D, T <: Real}
    # Unpack parameters
    @unpack jtlsproject, H_max = params
    # Projection onto the admissible region
    !jtlsproject && return evalfit(orbit)
    # Reference epoch
    t0 = epoch(orbit)
    # Lest squares fit
    @unpack fit = orbit
    # Barycentric initial condition
    q0 = orbit()
    # Geocentric attributable elements
    attr = cartesian2attributable(q0(fit.x)[1:6] - params.eph_ea(t0))
    ρ, v_ρ = attr[5], attr[6]
    h = H_max + 5 * log10(R_EA)
    # Admissible region
    A = AdmissibleRegion(
        days2dtutc(t0), deg2rad(attr[1]), deg2rad(attr[2]), deg2rad(attr[3]),
        deg2rad(attr[4]), h, search_observatory_code("500"), params
    )
    iszero(A) && return evalfit(orbit)
    # Projection onto the admissible region
    if !((ρ, v_ρ) in A)
        ρ, v_ρ = boundary_projection(A, ρ, v_ρ)
        # Barycentric initial condition
        q1 = attributable2cartesian(SVector{6, T}(attr[1], attr[2], attr[3],
                attr[4], ρ, v_ρ)) + params.eph_ea(t0)
        # New deltas and covariance matrix
        q0cart = view(q0, 1:6)
        @. fit.x[1:6] = (q1 - cte(q0cart)) / scalingfactor(q0cart)
        for method in lsmethods
            if typeof(method) == fit.routine
                fit.Γ .= inv(normalmatrix(method, fit.x))
            end
        end
    end
    # Evaluate deltas
    return evalfit(orbit)
end

function convergenceconditions(i::Int, Qs::AbstractVector, outs, params::Parameters)
    @unpack outrej = params
    if i > 1
        C1 = abs(Qs[i-1] - Qs[i]) < 0.01
        C2 = outrej ? (outs[i-1] == outs[i]) : true
        C3 = i > 2 && Qs[i-2] < Qs[i-1] < Qs[i]
        return (C1 && C2) || C3
    end
    return false
end

"""
    jtls(od, orbit, params [, mode::Bool]; kwargs...)

Given an orbit determination problem `od`, return a `LeastSquaresOrbit`
computed via the Jet Transport Least Squares method, starting from a
preliminary `orbit`. For a list of parameters, see the `Least Squares`
section of [`Parameters`](@ref).

The `mode` optional argument determines how many tracklets are added
to the fit per iteration, either as much as possible (`true`, default)
or one (`false`).

# Keyword arguments

- `buffer::Union{Nothing, JTLSBuffer}`: pre-allocated memory (default: `nothing`).
"""
function jtls(
        od::OpticalODProblem{D, T, O}, orbit::AbstractOrbit,
        params::Parameters{T}, mode::Bool = true;
        buffer::Union{Nothing, AbstractJTLSBuffer{D, T}} = nothing
    ) where {D, T <: Real, O <: AbstractOpticalVector{T}}
    # Set jet transport variables
    Ndof = dof(od)
    Npar = numvars(Val(od.dynamics), params)
    set_od_order(params, Npar)
    # Buffer
    if isnothing(buffer)
        buffer = JTLSBuffer(od, orbit, params)
    end
    # Unpack
    Qtol, Mtol = params.lsQtol, params.lsMtol
    @unpack jtlsorder, jtlsiter, outrej, jtlsmask, χ2_rec, χ2_rej,
            fudge, max_per = params
    @unpack orbits, res, Qs, q00s, outs, prbuffer, lscache, orcache = buffer
    # Reference epoch [Julian days TDB]
    jd0 = meanepoch(od) + PE.J2000
    # Jet transport initial condition
    variables = fittedvariables(Ndof, params)
    q0 = jtinitialcondition(od, orbit, jd0, params)
    # Least squares methods
    x0 = zeros(T, Npar)
    lsmethods = _lsmethods(res, x0, 1:Npar)
    # Initial subset of astrometry for orbit fit
    trksin, trksout, oidxs = _initialtracklets(od.tracklets, orbit.tracklets)
    # Jet Transport Least Squares
    TS.evaluate!(q0, x0, view(q00s, :, 1))
    for i in 1:jtlsiter
        # Initial conditions
        TS.constant_term!.(q0, view(q00s, :, i))
        # Decide whether q0 is suitable for jtls
        if jtlsmask
            isjtlsfit(od, q0, jd0, params) || break
        end
        # Propagation & residuals
        bwd, fwd = propres!(res, od, q0, jd0, params; buffer = prbuffer)
        isempty(res) && break
        # Orbit fit
        fit = tryls(view(res, oidxs), x0, lscache, lsmethods; Qtol, Mtol)
        !issuccess(fit) && break
        # Incrementally add observations to fit
        oidxs, fit = addoptical!(Val(mode), oidxs, fit, lscache, lsmethods,
                                 trksin, trksout, res, x0, params)
        fit.Γ .= project(q0[variables], fit)
        # Outlier rejection
        if outrej
            mro = view(outliers(od), oidxs)
            outlier_rejection!(view(res, oidxs), fit.x, fit.Γ, orcache;
                mro, χ2_rec, χ2_rej, fudge, max_per)
        end
        # Update orbit
        orbits[i] = updateorbit(LeastSquaresOrbit(
            od.dynamics, variables, od.optical[oidxs], trksin, nothing, bwd, fwd,
            res[oidxs], nothing, fit, q00s[variables, 1:i], Qs[1:i]
        ), lsmethods, params)
        Qs[i] = orbits[i].Qs[end] = nrms(orbits[i])
        if outrej
            outs[i] = notout(orbits[i])
        end
        # Convergence conditions
        if convergenceconditions(i, Qs, outs, params)
            any(isnan, sigmas(orbits[i])) && @warn "Final covariance matrix \
                is not positive-definite"
            return orbits[i]
        end
        # Update initial condition
        TS.evaluate!(q0, fit.x, view(q00s, :, i+1))
    end
    return zero(eltype(orbits))
end

function jtls(
        od::MixedODProblem{D, T, O, R}, orbit::AbstractOrbit,
        params::Parameters{T}, mode::Bool = true;
        buffer::Union{Nothing, AbstractJTLSBuffer{D, T}} = nothing
    ) where {D, T <: Real, O <: AbstractOpticalVector{T}, R <: AbstractRadarVector{T}}
    # Set jet transport variables
    Ndof = dof(od)
    Npar = numvars(Val(od.dynamics), params)
    set_od_order(params, Npar)
    # Buffer
    if isnothing(buffer)
        buffer = JTLSBuffer(od, orbit, params)
    end
    # Unpack
    Qtol, Mtol = params.lsQtol, params.lsMtol
    @unpack jtlsorder, jtlsiter, outrej, jtlsmask, χ2_rec, χ2_rej,
            fudge, max_per = params
    @unpack orbits, res, Qs, q00s, outs, prbuffer, lscache, orcache = buffer
    # Reference epoch [Julian days TDB]
    jd0 = meanepoch(od) + PE.J2000
    # Jet transport initial condition
    variables = fittedvariables(Ndof, params)
    q0 = jtinitialcondition(od, orbit, jd0, params)
    # Least squares methods
    x0 = zeros(T, Npar)
    lsmethods = _lsmethods(res, x0, 1:Npar)
    # Initial subset of astrometry for orbit fit
    trksin, trksout, oidxs = _initialtracklets(od.tracklets, orbit.tracklets)
    radarin, radarout, ridxs = _initialradar(od, orbit)
    # Jet Transport Least Squares
    TS.evaluate!(q0, x0, view(q00s, :, 1))
    for i in 1:jtlsiter
        # Initial conditions
        TS.constant_term!.(q0, view(q00s, :, i))
        # Decide whether q0 is suitable for jtls
        # if jtlsmask
        #     isjtlsfit(od, q0, jd0, params) || break
        # end
        # Propagation & residuals
        bwd, fwd = propres!(res, od, q0, jd0, params; buffer = prbuffer)
        any(isempty, res) && break
        # Orbit fit
        fit = tryls((res[1][oidxs], res[2][ridxs]), x0, lscache, lsmethods; Qtol, Mtol)
        !issuccess(fit) && break
        # Incrementally add observations to fit
        oidxs, ridxs, fit = addobservations!(od, oidxs, ridxs, fit, lscache, lsmethods,
            trksin, trksout, radarin, radarout, res, x0, params)
        fit.Γ .= project(q0[variables], fit)
        # Outlier rejection
        if outrej
            mro = view(outliers(od), oidxs)
            outlier_rejection!(view(res[1], oidxs), fit.x, fit.Γ, orcache;
                mro, χ2_rec, χ2_rej, fudge, max_per)
        end
        # Update orbit
        orbits[i] = updateorbit(LeastSquaresOrbit(
            od.dynamics, variables, od.optical[oidxs], trksin, od.radar[ridxs], bwd, fwd,
            res[1][oidxs], res[2][ridxs], fit, q00s[variables, 1:i], Qs[1:i]
        ), lsmethods, params)
        Qs[i] = orbits[i].Qs[end] = nrms(orbits[i])
        if outrej
            outs[i] = notout(orbits[i])
        end
        # Convergence conditions
        if convergenceconditions(i, Qs, outs, params)
            any(isnan, sigmas(orbits[i])) && @warn "Final covariance matrix \
                is not positive-definite"
            return orbits[i]
        end
        # Update initial condition
        TS.evaluate!(q0, fit.x, view(q00s, :, i+1))
    end
    return zero(eltype(orbits))
end

"""
    linkage(od, porbit, params; kwargs...)

Fit an orbit to an orbit determination problem `od` by successively
scaling the weights of the observations that are in `od` but not in
a preliminary orbit `porbit`. For a list of parameters, see the
`Least Squares` section of [`Parameters`](@ref).

# Keyword arguments

- `maxiter::Int`: maximum number of iterations (default: `10`).
"""
function linkage(
        od::OpticalODProblem{D, T, O}, porbit::AbstractOrbit, params::Parameters{T};
        maxiter::Int = 10
    ) where {D, T <: Real, O <: AbstractOpticalVector{T}}
    # Unpack
    @unpack verbose = params
    # Reference epoch [Julian days TDB]
    jd0 = epoch(porbit) + PE.J2000
    # Find observations in od that are not in porbit
    trks = setdiff(od.tracklets, porbit.tracklets)
    idxs = indices(trks)
    # Original weights
    w8s = deepcopy(od.weights.weights)
    # Dynamical initial condition
    q0 = initialcondition(porbit, jd0, dof(od), params)
    # Initialize buffer and set of residuals
    buffer = PropresBuffer(od, q0, jd0, params)
    res = init_residuals(T, od, porbit)
    orbit = zero(OpticalLeastSquaresOrbit{D, T, T, O})
    # Select first scaling factor
    mags = Vector{Int}(undef, length(trks))
    propres!(res, od, q0, jd0, params; buffer)
    if isempty(res)
        verbose && @warn("Linkage did not converge within the given parameters \
            or could not fit all the astrometry")
        od.weights.weights .= w8s
        return orbit
    end
    for i in eachindex(mags)
        x = maximum(log10chi, view(res, indices(trks[i])))
        mags[i] = ceil(Int, x)
    end
    k = max(0, minimum(mags) - 1)
    # Main cycle
    for i in 1:maxiter
        # Scale down weights
        for j in idxs
            wra, wdec = w8s[j] ./ 10^k
            res[j] = OpticalResidual{T, T}(ra(res[j]), dec(res[j]), wra, wdec,
                dra(res[j]), ddec(res[j]), corr(res[j]), isoutlier(res[j]))
            od.weights.weights[j] = (wra, wdec)
        end
        # Jet transport least squares
        if i == 1
            orbit = jtls(od, porbit, params)
        else
            orbit = jtls(od, orbit, params)
        end
        # Break conditions
        if iszero(k) || nrms(orbit) > 10
            if nrms(orbit) < 10 && nobs(orbit) == nobs(od)
                verbose && println(
                    "* Linkage converged in $i iterations to: \n\n",
                    summary(orbit)
                )
            else
                verbose && @warn("Linkage did not converge within the \
                    given parameters or could not fit all the astrometry")
            end
            break
        end
        # Scale up weights
        q0 = initialcondition(orbit, jd0, dof(od), params)
        propres!(res, od, q0, jd0, params; buffer)
        if isempty(res)
            verbose && @warn("Linkage did not converge within the given parameters \
                or could not fit all the astrometry")
            break
        end
        for i in eachindex(mags)
            x = maximum(log10chi, view(res, indices(trks[i])))
            mags[i] = ceil(Int, x)
        end
        k += min(-1, minimum(mags))
        k = max(0, k)
    end
    # Return weights to original
    od.weights.weights .= w8s

    return orbit
end
