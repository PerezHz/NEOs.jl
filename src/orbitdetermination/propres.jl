# Times used within propres
function _proprestimes(radec::AbstractVector{RadecMPC{T}}, jd0::U,
    params::Parameters{T}) where {T <: Real, U <: Number}
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
function propres(od::ODProblem{D, T}, jd0::V, q0::Vector{U}, params::Parameters{T};
    buffer::Union{Nothing, PropagationBuffer{T, U, V}} = nothing,
    idxs::AbstractVector{Int} = eachindex(od.radec)) where {D, T <: Real, U <: Number,
    V <: Number}
    # Unpack parameters
    @unpack coeffstol, eph_su, eph_ea = params
    # Subset of radec for propagation and residuals
    radec = view(od.radec, idxs)
    # Times of first/last observation, epoch and years in backward/forward propagation
    t0, tf, _jd0_, nyears_bwd, nyears_fwd = _proprestimes(radec, jd0, params)
    # Propagation buffer
    if isnothing(buffer)
        buffer = PropagationBuffer(od, jd0, idxs[1], idxs[end], q0, params)
    end
    # Backward (forward) integration
    bwd = _propagate(od.dynamics, jd0, nyears_bwd, q0, buffer, params)
    fwd = _propagate(od.dynamics, jd0, nyears_fwd, q0, buffer, params)
    if !issuccessfulprop(bwd, t0 - _jd0_; tol = coeffstol) ||
       !issuccessfulprop(fwd, tf - _jd0_; tol = coeffstol)
        return bwd, fwd, Vector{OpticalResidual{T, U}}(undef, 0)
    end
    # O-C residuals
    res = init_residuals(U, od, idxs)
    try
        residuals!(res, radec;
            xvs = et -> auday2kmsec(eph_su(et/daysec)),
            xve = et -> auday2kmsec(eph_ea(et/daysec)),
            xva = et -> bwdfwdeph(et, bwd, fwd))
        return bwd, fwd, res
    catch
        empty!(res)
        return bwd, fwd, res
    end
end

# In-place method of propres
function propres!(res::Vector{OpticalResidual{T, U}}, od::ODProblem{D, T},
    jd0::V, q0::Vector{U}, params::Parameters{T};
    buffer::Union{Nothing, PropagationBuffer{T, U, V}} = nothing,
    idxs::AbstractVector{Int} = eachindex(od.radec))  where {D, T <: Real, U <: Number,
    V <: Number}
    # Unpack parameters
    @unpack coeffstol, eph_su, eph_ea = params
    # Subset of radec for propagation and residuals
    radec = view(od.radec, idxs)
    # Times of first/last observation, epoch and years in backward/forward propagation
    t0, tf, _jd0_, nyears_bwd, nyears_fwd = _proprestimes(radec, jd0, params)
    # Propagation buffer
    if isnothing(buffer)
        buffer = PropagationBuffer(od, jd0, idxs[1], idxs[end], q0, params)
    end
    # Backward (forward) integration
    bwd = _propagate(od.dynamics, jd0, nyears_bwd, q0, buffer, params)
    fwd = _propagate(od.dynamics, jd0, nyears_fwd, q0, buffer, params)
    if !issuccessfulprop(bwd, t0 - _jd0_; tol = coeffstol) ||
       !issuccessfulprop(fwd, tf - _jd0_; tol = coeffstol)
        empty!(res)
        return bwd, fwd
    end
    # O-C residuals
    try
        residuals!(res, radec;
            xvs = et -> auday2kmsec(eph_su(et/daysec)),
            xve = et -> auday2kmsec(eph_ea(et/daysec)),
            xva = et -> bwdfwdeph(et, bwd, fwd))
        return bwd, fwd
    catch
        empty!(res)
        return bwd, fwd
    end
end