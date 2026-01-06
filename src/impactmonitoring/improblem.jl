"""
    IMProblem{D, T, O <: AbstractOrbit{D, T, T}} <: AbstractIMProblem{D, T}

An impact monitoring problem.

# Fields

- `orbit::O`: reference orbit.
- `target::ImpactTarget{T}`: target celestial body.
"""
struct IMProblem{D, T, O <: AbstractOrbit{D, T, T}} <: AbstractIMProblem{D, T}
    orbit::O
    target::ImpactTarget{T}
end

# Outer constructor
IMProblem(orbit::AbstractOrbit{D, T, T}, target::ImpactTarget{T}) where {D, T <: Real} =
    IMProblem{D, T, typeof(orbit)}(orbit, target)

# Print method for IMProblem
function show(io::IO, x::IMProblem)
    t = repeat(' ', 4)
    print(io,
        "Impact monitoring problem\n",
        t, rpad("Orbit:", 10), string(x.orbit), "\n",
        t, rpad("Target:", 10), string(x.target), "\n",
    )
end

# AbstractIMProblem interface
dynamicalmodel(x::IMProblem) = x.orbit.dynamics
dof(x::IMProblem) = dof(Val(dynamicalmodel(x)))

noptical(x::IMProblem) = length(x.orbit.optical)
opticalindices(x::IMProblem) = eachindex(x.orbit.optical)

function minmaxdates(x::IMProblem)
    t0, tf = minmaxdates(x.orbit.optical)
    if hasradar(x.orbit)
        _t0_, _tf_ = minmaxdates(x.orbit.radar)
        t0, tf = min(t0, _t0_), max(tf, _tf_)
    end
    return t0, tf
end

function init_optical_residuals(
        ::Type{U}, x::AbstractIMProblem{D, T}
    ) where {D, T <: Real, U <: Number}
    # Initialize vector of optical residuals
    res = Vector{OpticalResidual{T, U}}(undef, length(x.orbit.ores))
    for (i, r) in enumerate(x.orbit.ores)
        res[i] = OpticalResidual{T, U}(zero(U), zero(U), wra(r), wdec(r), dra(r), ddec(r),
                                       corr(r), isoutlier(r))
    end
    return res
end

function PropresBuffer(
        IM::AbstractIMProblem{D, T}, q0::Vector{U},
        jd0::V, idxs::AbstractVector{Int}, params::Parameters{T}
    ) where {D, T <: Real, U <: Number, V <: Number}
    t0 = dtutc2days(date(IM.orbit.optical[idxs[1]]))
    tf = dtutc2days(date(IM.orbit.optical[idxs[end]]))
    tlim = (t0 - params.bwdoffset, tf + params.fwdoffset)
    prop = PropagationBuffer(dynamicalmodel(IM), q0, jd0, tlim, params)
    res = [OpticalBuffer(q0[1]) for _ in eachindex(idxs)]
    return PropresBuffer{T, U, V}(prop, res)
end

function PropresBuffer(
        IM::AbstractIMProblem{D, T}, q0::Vector{U},
        jd0::V, params::Parameters
    ) where {D, T <: Real, U <: Number, V <: Number}
    t0, tf = dtutc2days.(minmaxdates(IM))
    tlim = (t0 - params.bwdoffset, tf + params.fwdoffset)
    prop = PropagationBuffer(dynamicalmodel(IM), q0, jd0, tlim, params)
    res = [OpticalBuffer(q0[1]) for _ in 1:noptical(IM)]
    return PropresBuffer{T, U, V}(prop, res)
end

function propres!(
        res::Vector{OpticalResidual{T, U}}, IM::AbstractIMProblem{D, T},
        q0::Vector{U}, jd0::V, params::Parameters{T};
        buffer::Union{Nothing, PropresBuffer{T, U, V}} = nothing,
        idxs::AbstractVector{Int} = opticalindices(IM)
    )  where {D, T <: Real, U <: Number, V <: Number}
    # Unpack parameters
    @unpack coeffstol, eph_su, eph_ea = params
    # Subset of optical astrometry for propagation and residuals
    optical = view(IM.orbit.optical, idxs)
    # Times of first/last observation, epoch and years in backward/forward propagation
    t0, tf, _jd0_, nyears_bwd, nyears_fwd = _proprestimes(optical, jd0, params)
    # Buffer
    if isnothing(buffer)
        buffer = PropresBuffer(IM, q0, jd0, idxs, params)
    end
    # Backward (forward) integration
    bwd = _propagate(dynamicalmodel(IM), q0, jd0, nyears_bwd, buffer.prop, params)
    fwd = _propagate(dynamicalmodel(IM), q0, jd0, nyears_fwd, buffer.prop, params)
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