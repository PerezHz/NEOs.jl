"""
    ImpactMonitoringBuffer{T <: Real} <: AbstractBuffer

Pre-allocated memory for [`closeapproaches`](@ref).

# Fields

- `prop::PropagationBuffer{T, Taylor1{T}, T}`: propagation buffer.
- `teph::EphemerisEvaluationBuffer{T, Taylor1{T}}`: target ephemeris
    evaluation buffer.
- `tvS/xvS/gvS::Array{Taylor1{T}, _}`: root-finding arrays.
"""
struct ImpactMonitoringBuffer{T <: Real} <: AbstractBuffer
    prop::PropagationBuffer{T, Taylor1{T}, T}
    teph::EphemerisEvaluationBuffer{T, Taylor1{T}}
    tvS::Vector{Taylor1{T}}
    xvS::Matrix{Taylor1{T}}
    gvS::Vector{Taylor1{T}}
end

"""
    ImpactMonitoringBuffer(IM, nyears, vaorder, params)

Return an `ImpactMonitoringBuffer` object with pre-allocated
memory for `closeapproaches`.

# Arguments

- `IM::IMProblem`: impact monitoring problem.
- `nyears::Real`: number of years.
- `vaorder::Int`: order of Taylor expansions wrt LOV index.
- `params::Parameters`: see the `Propagation` section of [`Parameters`](@ref).
"""
function ImpactMonitoringBuffer(IM::AbstractIMProblem{D, T}, nyears::T, vaorder::Int,
                                params::Parameters{T}) where {D, T <: Real}
    # Unpack
    @unpack orbit, target = IM
    @unpack order, maxsteps = params
    # Propagation buffer
    jd0 = epoch(orbit) + PE.J2000
    q0 = orbit() .+ 1E-8 * Taylor1(vaorder)
    tlim = (epoch(orbit), epoch(orbit) + nyears * yr)
    prop = PropagationBuffer(dynamicalmodel(IM), q0, jd0, tlim, params)
    # Target ephemeris evaluation buffer
    teph = EphemerisEvaluationBuffer(target.eph, tlim, order, q0)
    # Root-finding arrays
    tvS = Array{Taylor1{T}}(undef, maxsteps + 1)
    xvS = Array{Taylor1{T}}(undef, length(q0), maxsteps + 1)
    gvS = similar(tvS)

    return ImpactMonitoringBuffer{T}(prop, teph, tvS, xvS, gvS)
end