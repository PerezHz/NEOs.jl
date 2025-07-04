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
function _proprestimes(optical::AbstractOpticalVector{T}, jd0::Number,
                       params::Parameters{T}) where {T <: Real}
    # Time of first (last) observation
    t0, tf = dtutc2jdtdb(date(optical[1])), dtutc2jdtdb(date(optical[end]))
    # TDB epoch (plain)
    _jd0_ = cte(cte(jd0))
    # Years in backward (forward) integration
    nyears_bwd = -(_jd0_ - t0 + params.bwdoffset) / yr
    nyears_fwd = (tf - _jd0_ + params.fwdoffset) / yr

    return t0, tf, _jd0_, nyears_bwd, nyears_fwd
end

"""
    propres(od, q0, jd0, params; kwargs...)

Propagate initial condition `q0` [au, au/day] referred to time `jd0` [julian date TDB]
for the time needed to cover the optical astrometry in `od`. In addition, compute the
O-C residuals of the resulting orbit with respect to the optical astrometry in `od`.

See also [`propagate`](@ref) and [`residuals`](@ref).

# Keyword arguments

- `buffer::Union{Nothing, PropagationBuffer}`: propagation buffer
    (default: `nothing`).
- `idxs::AbstractVector{Int}`: indices of the observations in `od.optical` to be included
    in the computation.
"""
function propres(
        od::AbstractODProblem{D, T}, q0::Vector{U}, jd0::V, params::Parameters{T};
        buffer::Union{Nothing, PropagationBuffer{T, U, V}} = nothing,
        idxs::AbstractVector{Int} = eachindex(od.optical)
    ) where {D, T <: Real, U <: Number, V <: Number}
    # Unpack parameters
    @unpack coeffstol, eph_su, eph_ea = params
    # Subset of optical astrometry for propagation and residuals
    optical = view(od.optical, idxs)
    # Times of first/last observation, epoch and years in backward/forward propagation
    t0, tf, _jd0_, nyears_bwd, nyears_fwd = _proprestimes(optical, jd0, params)
    # Propagation buffer
    if isnothing(buffer)
        buffer = PropagationBuffer(od, q0, jd0, idxs[1], idxs[end], params)
    end
    # Backward (forward) integration
    bwd = _propagate(od.dynamics, q0, jd0, nyears_bwd, buffer, params)
    fwd = _propagate(od.dynamics, q0, jd0, nyears_fwd, buffer, params)
    if !issuccessfulprop(bwd, t0 - _jd0_; tol = coeffstol) ||
       !issuccessfulprop(fwd, tf - _jd0_; tol = coeffstol)
        return bwd, fwd, Vector{OpticalResidual{T, U}}(undef, 0)
    end
    # O-C residuals
    res = init_optical_residuals(U, od, idxs)
    try
        residuals!(res, optical;
            xvs = et -> auday2kmsec(eph_su(et/daysec)),
            xve = et -> auday2kmsec(eph_ea(et/daysec)),
            xva = et -> bwdfwdeph(et, bwd, fwd))
        return bwd, fwd, res
    catch
        empty!(res)
        return bwd, fwd, res
    end
end

"""
    propres!(res, od, q0, jd0, params; kwargs...)

Equivalent to [`propres`](@ref), but computes the O-C residuals in-place over
a pre-allocated vector `res`.
"""
function propres!(
        res::Vector{OpticalResidual{T, U}}, od::AbstractODProblem{D, T},
        q0::Vector{U}, jd0::V, params::Parameters{T};
        buffer::Union{Nothing, PropagationBuffer{T, U, V}} = nothing,
        idxs::AbstractVector{Int} = eachindex(od.optical)
    )  where {D, T <: Real, U <: Number, V <: Number}
    # Unpack parameters
    @unpack coeffstol, eph_su, eph_ea = params
    # Subset of optical astrometry for propagation and residuals
    optical = view(od.optical, idxs)
    # Times of first/last observation, epoch and years in backward/forward propagation
    t0, tf, _jd0_, nyears_bwd, nyears_fwd = _proprestimes(optical, jd0, params)
    # Propagation buffer
    if isnothing(buffer)
        buffer = PropagationBuffer(od, q0, jd0, idxs[1], idxs[end], params)
    end
    # Backward (forward) integration
    bwd = _propagate(od.dynamics, q0, jd0, nyears_bwd, buffer, params)
    fwd = _propagate(od.dynamics, q0, jd0, nyears_fwd, buffer, params)
    if !issuccessfulprop(bwd, t0 - _jd0_; tol = coeffstol) ||
       !issuccessfulprop(fwd, tf - _jd0_; tol = coeffstol)
        empty!(res)
        return bwd, fwd
    end
    # O-C residuals
    try
        residuals!(res, optical;
            xvs = et -> auday2kmsec(eph_su(et/daysec)),
            xve = et -> auday2kmsec(eph_ea(et/daysec)),
            xva = et -> bwdfwdeph(et, bwd, fwd))
        return bwd, fwd
    catch
        empty!(res)
        return bwd, fwd
    end
end