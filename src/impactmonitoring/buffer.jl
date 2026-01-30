"""
    ImpactMonitoringBuffer{T <: Real, U <: Number} <: AbstractBuffer

Pre-allocated memory for [`closeapproaches`](@ref).

# Fields

- `prop::PropagationBuffer{T, U, T}`: propagation buffer.
- `teph::EphemerisEvaluationBuffer{T, U}`: target ephemeris evaluation buffer.
- `tvS/xvS/gvS::Array{U, _}`: root-finding arrays.
"""
struct ImpactMonitoringBuffer{T <: Real, U <: Number} <: AbstractBuffer
    prop::PropagationBuffer{T, U, T}
    teph::EphemerisEvaluationBuffer{T, U}
    tvS::Vector{U}
    xvS::Matrix{U}
    gvS::Vector{U}
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
    # Root-finding arrays
    tvS = Array{U}(undef, maxsteps + 1)
    xvS = Array{U}(undef, length(q0), maxsteps + 1)
    gvS = similar(tvS)

    return ImpactMonitoringBuffer{T, U}(prop, teph, tvS, xvS, gvS)
end