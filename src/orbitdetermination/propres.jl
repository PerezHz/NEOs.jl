"""
    PropresBuffer{T <: Real, U <: Number, V <: Number} <: AbstractBuffer

Pre-allocated memory for [`propres`](@ref).

# Fields

- `prop::PropagationBuffer{T, U, V}`: buffer for [`propagate`](@ref).
- `res::Vector{OpticalBuffer{U}}`: buffer for [`compute_radec`](@ref).
"""
struct PropresBuffer{T <: Real, U <: Number, V <: Number} <: AbstractBuffer
    prop::PropagationBuffer{T, U, V}
    res::Vector{OpticalBuffer{U}}
end

# Special PropresBuffer constructors
function PropresBuffer(
        od::AbstractODProblem{D, T}, q0::Vector{U},
        jd0::V, idxs::AbstractVector{Int}, params::Parameters{T}
    ) where {D, T <: Real, U <: Number, V <: Number}
    t0 = dtutc2days(date(od.optical[idxs[1]]))
    tf = dtutc2days(date(od.optical[idxs[end]]))
    tlim = (t0 - params.bwdoffset, tf + params.fwdoffset)
    prop = PropagationBuffer(od.dynamics, q0, jd0, tlim, params)
    res = [OpticalBuffer(q0[1]) for _ in eachindex(idxs)]
    return PropresBuffer{T, U, V}(prop, res)
end

function PropresBuffer(
        od::AbstractODProblem{D, T}, q0::Vector{U},
        jd0::V, params::Parameters
    ) where {D, T <: Real, U <: Number, V <: Number}
    t0, tf = dtutc2days.(minmaxdates(od))
    tlim = (t0 - params.bwdoffset, tf + params.fwdoffset)
    prop = PropagationBuffer(od.dynamics, q0, jd0, tlim, params)
    res = [OpticalBuffer(q0[1]) for _ in 1:noptical(od)]
    return PropresBuffer{T, U, V}(prop, res)
end

"""
    bwdfwdeph(t, bwd, fwd [, et [, kmsec]])

Evaluate an ephemerides at time `t`, where `bwd` (`fwd`) is the backward
(forward) propagation.

# Optional arguments

- `et::Bool`: whether `t` is in ephemeris seconds since J2000 (default: `true`).
- `kmsec::Bool`: whether to convert the state vector to [km, km/s] (default: `true`).
"""
function bwdfwdeph(t::Number, bwd::TaylorInterpolant, fwd::TaylorInterpolant,
                   et::Bool = true, kmsec::Bool = true)
    @assert bwd.t0 == fwd.t0 "Backward and forward initial times must match"
    _t_ = et ? t/daysec : t
    _rv_ = _t_ <= bwd.t0 ? bwd(_t_) : fwd(_t_)
    rv = kmsec ? auday2kmsec(_rv_) : _rv_
    return rv
end

# Times used within propres
function _proprestimes(x, jd0::Number, params::Parameters)
    # Time of first (last) observation
    t0, tf = dtutc2jdtdb.(minmaxdates(x))
    # TDB epoch (plain)
    _jd0_ = cte(cte(jd0))
    # Years in backward (forward) integration
    nyears_bwd = -(_jd0_ - t0 + params.bwdoffset) / yr
    nyears_fwd = (tf - _jd0_ + params.fwdoffset) / yr

    return t0, tf, _jd0_, nyears_bwd, nyears_fwd
end

"""
    propres(od, q0, jd0, params; kwargs...)

Propagate initial condition `q0` [au, au/day] referred to time `jd0` [julian date
TDB] for the time needed to cover the astrometry in `od`. In addition, compute the
O-C residuals of the resulting orbit with respect to the astrometry in `od`.

See also [`propagate`](@ref) and [`residuals`](@ref).

# Keyword arguments

- `buffer::Union{Nothing, PropresBuffer}`: pre-allocated memory (default: `nothing`).
- `idxs::AbstractVector{Int}`: indices of the observations in `od.optical` to be included
    in the computation.
"""
function propres(
        od::OpticalODProblem{D, T, O}, q0::Vector{U}, jd0::V, params::Parameters{T};
        buffer::Union{Nothing, PropresBuffer{T, U, V}} = nothing,
        idxs::AbstractVector{Int} = opticalindices(od)
    ) where {D, T <: Real, U <: Number, V <: Number, O <: AbstractOpticalVector{T}}
    # Unpack parameters
    @unpack coeffstol, eph_su, eph_ea = params
    # Subset of optical astrometry for propagation and residuals
    optical = view(od.optical, idxs)
    # Times of first/last observation, epoch and years in backward/forward propagation
    t0, tf, _jd0_, nyears_bwd, nyears_fwd = _proprestimes(optical, jd0, params)
    # Buffer
    if isnothing(buffer)
        buffer = PropresBuffer(od, q0, jd0, idxs, params)
    end
    # Backward (forward) integration
    bwd = _propagate(od.dynamics, q0, jd0, nyears_bwd, buffer.prop, params)
    fwd = _propagate(od.dynamics, q0, jd0, nyears_fwd, buffer.prop, params)
    if !issuccessfulprop(bwd, t0 - _jd0_; tol = coeffstol) ||
       !issuccessfulprop(fwd, tf - _jd0_; tol = coeffstol)
        return bwd, fwd, Vector{OpticalResidual{T, U}}()
    end
    # O-C residuals
    res = init_optical_residuals(U, od, idxs)
    try
        residuals!(res, optical, buffer.res; xvs = eph_su, xve = eph_ea,
                   xva = (bwd, fwd))
        return bwd, fwd, res
    catch
        empty!(res)
        return bwd, fwd, res
    end
end

function propres(
        od::MixedODProblem{D, T, O, R}, q0::Vector{U}, jd0::V, params::Parameters{T};
        buffer::Union{Nothing, PropresBuffer{T, U, V}} = nothing,
    ) where {D, T <: Real, U <: Number, V <: Number, O <: AbstractOpticalVector{T},
             R <: AbstractRadarVector{T}}
    # Unpack
    @unpack coeffstol, eph_su, eph_ea = params
    @unpack dynamics, optical, radar = od
    # Times of first/last observation, epoch and years in backward/forward propagation
    t0, tf, _jd0_, nyears_bwd, nyears_fwd = _proprestimes(od, jd0, params)
    # Buffer
    if isnothing(buffer)
        buffer = PropresBuffer(od, q0, jd0, params)
    end
    # Backward (forward) integration
    bwd = _propagate(dynamics, q0, jd0, nyears_bwd, buffer.prop, params)
    fwd = _propagate(dynamics, q0, jd0, nyears_fwd, buffer.prop, params)
    if !issuccessfulprop(bwd, t0 - _jd0_; tol = coeffstol) ||
       !issuccessfulprop(fwd, tf - _jd0_; tol = coeffstol)
        return bwd, fwd, (Vector{OpticalResidual{T, U}}(), Vector{RadarResidual{T, U}}())
    end
    # O-C residuals
    res = (init_optical_residuals(U, od), init_radar_residuals(U, od))
    try
        residuals!(res[1], optical, buffer.res; xvs = eph_su, xve = eph_ea,
                   xva = (bwd, fwd))
        residuals!(res[2], radar;
            xvs = et -> auday2kmsec(eph_su(et/daysec)),
            xve = et -> auday2kmsec(eph_ea(et/daysec)),
            xva = et -> bwdfwdeph(et, bwd, fwd))
        return bwd, fwd, res
    catch
        empty!(res[1])
        empty!(res[2])
        return bwd, fwd, res
    end
end

"""
    propres!(res, od, q0, jd0, params; kwargs...)

Equivalent to [`propres`](@ref), but computes the O-C residuals in-place over
a pre-allocated set `res`.
"""
function propres!(
        res::Vector{OpticalResidual{T, U}}, od::OpticalODProblem{D, T, O},
        q0::Vector{U}, jd0::V, params::Parameters{T};
        buffer::Union{Nothing, PropresBuffer{T, U, V}} = nothing,
        idxs::AbstractVector{Int} = opticalindices(od)
    )  where {D, T <: Real, U <: Number, V <: Number, O <: AbstractOpticalVector{T}}
    # Unpack parameters
    @unpack coeffstol, eph_su, eph_ea = params
    # Subset of optical astrometry for propagation and residuals
    optical = view(od.optical, idxs)
    # Times of first/last observation, epoch and years in backward/forward propagation
    t0, tf, _jd0_, nyears_bwd, nyears_fwd = _proprestimes(optical, jd0, params)
    # Buffer
    if isnothing(buffer)
        buffer = PropresBuffer(od, q0, jd0, idxs, params)
    end
    # Backward (forward) integration
    bwd = _propagate(od.dynamics, q0, jd0, nyears_bwd, buffer.prop, params)
    fwd = _propagate(od.dynamics, q0, jd0, nyears_fwd, buffer.prop, params)
    if !issuccessfulprop(bwd, t0 - _jd0_; tol = coeffstol) ||
       !issuccessfulprop(fwd, tf - _jd0_; tol = coeffstol)
        empty!(res)
        return bwd, fwd
    end
    # O-C residuals
    try
        residuals!(res, optical, buffer.res; xvs = eph_su, xve = eph_ea,
                   xva = (bwd, fwd))
        return bwd, fwd
    catch
        empty!(res)
        return bwd, fwd
    end
end

function propres!(
        res::Tuple{Vector{OpticalResidual{T, U}}, Vector{RadarResidual{T, U}}},
        od::MixedODProblem{D, T, O, R}, q0::Vector{U}, jd0::V, params::Parameters{T};
        buffer::Union{Nothing, PropresBuffer{T, U, V}} = nothing,
    )  where {D, T <: Real, U <: Number, V <: Number, O <: AbstractOpticalVector{T},
              R <: AbstractRadarVector{T}}
    # Unpack
    @unpack coeffstol, eph_su, eph_ea = params
    @unpack dynamics, optical, radar = od
    # Times of first/last observation, epoch and years in backward/forward propagation
    t0, tf, _jd0_, nyears_bwd, nyears_fwd = _proprestimes(od, jd0, params)
    # Buffer
    if isnothing(buffer)
        buffer = PropresBuffer(od, q0, jd0, params)
    end
    # Backward (forward) integration
    bwd = _propagate(dynamics, q0, jd0, nyears_bwd, buffer.prop, params)
    fwd = _propagate(dynamics, q0, jd0, nyears_fwd, buffer.prop, params)
    if !issuccessfulprop(bwd, t0 - _jd0_; tol = coeffstol) ||
       !issuccessfulprop(fwd, tf - _jd0_; tol = coeffstol)
        empty!(res[1])
        empty!(res[2])
        return bwd, fwd
    end
    # O-C residuals
    try
        residuals!(res[1], optical, buffer.res; xvs = eph_su, xve = eph_ea,
                   xva = (bwd, fwd))
        residuals!(res[2], radar;
            xvs = et -> auday2kmsec(eph_su(et/daysec)),
            xve = et -> auday2kmsec(eph_ea(et/daysec)),
            xva = et -> bwdfwdeph(et, bwd, fwd))
        return bwd, fwd
    catch
        empty!(res[1])
        empty!(res[2])
        return bwd, fwd
    end
end