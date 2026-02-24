"""
    ImpactMonitoringBuffer{T <: Real, U <: Number} <: AbstractBuffer

Pre-allocated memory for [`closeapproaches`](@ref).

# Fields

- `prop::PropagationBuffer{T, U, T}`: propagation buffer.
- `teph::EphemerisEvaluationBuffer{T, U}`: target ephemeris evaluation buffer.
- `g_constant::Vector{Taylor1{T}}`: surface crossing buffer.
- `g_dg::Vector{Taylor1{U}}`: objective function buffer.
- `g_dg_val::Vector{U}`: objective function evaluation buffer.
"""
struct ImpactMonitoringBuffer{T <: Real, U <: Number} <: AbstractBuffer
    prop::PropagationBuffer{T, U, T}
    teph::EphemerisEvaluationBuffer{T, U}
    g_constant::Vector{Taylor1{T}}
    g_dg::Vector{Taylor1{U}}
    g_dg_val::Vector{U}
end

"""
    ImpactMonitoringBuffer(IM, q0, nyears, params)

Return an `ImpactMonitoringBuffer` object with pre-allocated
memory for [`closeapproaches`](@ref).

# Arguments

- `IM::IMProblem`: impact monitoring problem.
- `q0::Vector{<:Number}`: initial condition.
- `nyears::Real`: number of years.
- `params::Parameters`: see the `Propagation` section of [`Parameters`](@ref).
"""
function ImpactMonitoringBuffer(IM::AbstractIMProblem{D, T}, q0::Vector{U}, nyears::T,
                                params::Parameters{T}) where {D, T <: Real, U <: Number}
    # Unpack
    @unpack orbit, target = IM
    @unpack order, maxsteps = params
    # Propagation buffer
    jd0 = epoch(orbit) + PE.J2000
    tlim = (epoch(orbit), epoch(orbit) + nyears * yr)
    prop = PropagationBuffer(dynamicalmodel(IM), q0, jd0, tlim, params)
    # Target ephemeris evaluation buffer
    teph = EphemerisEvaluationBuffer(target.eph, tlim, order, q0)
    # Objective function buffers
    @unpack t, x = prop.cache
    g_constant = [zero(t), zero(t)]
    g_dg = [zero(x[1]), zero(x[1])]
    g_dg_val = [zero(q0[1]), zero(q0[1])]

    return ImpactMonitoringBuffer{T, U}(prop, teph, g_constant, g_dg, g_dg_val)
end