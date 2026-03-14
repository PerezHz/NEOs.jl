"""
    RootFindingEvent{U <: Number}

An evaluation of an event function, e.g. [`closeapproach!`](@ref).

# Fields

- `flag::Bool`: a boolean flag indicating whether the event is
    considered or not.
- `func::Taylor1{U}`: current value of the event function.
"""
mutable struct RootFindingEvent{U <: Number}
    flag::Bool
    func::Taylor1{U}
end

# RootFindingEvent interface
first(x::RootFindingEvent) = x.flag
last(x::RootFindingEvent) = x.func
scalarzero(x::RootFindingEvent) = zero(constant_term(x))
constant_term(x::RootFindingEvent) = constant_term(last(x))
surfacecrossing(old::RootFindingEvent, new::RootFindingEvent, eventorder::Int) =
    surfacecrossing((first(old), last(old)), (first(new), last(new)), eventorder)

function identity!(x::RootFindingEvent, y::RootFindingEvent)
    x.flag = y.flag
    for k in eachindex(x.func)
        TS.identity!(x.func, y.func, k)
    end
    return nothing
end

"""
    RootFindingBuffer{T <: Real, U <: Number} <: AbstractBuffer

Pre-allocated memory for the root-finding in [`closeapproaches`](@ref).

# Fields

- `jd0::T`: reference epoch [JDTDB].
- `rv::RetAlloc{Taylor1{U}}`: [`closeapproach!`](@ref) buffer.
- `teph::EphemerisEvaluationBuffer{T, U}`: target ephemeris evaluation buffer.
- `g_dg_val::Vector{U}`, `g_dg::Vector{Taylor1{U}}` and `g_constant::Vector{Taylor1{T}}`:
    [`findroot`](@ref) buffer.
- `f_tupl/g_tupl/f_tupl_old/g_tupl_old::RootFindingEvent{U}`: pre-allocated events.
"""
mutable struct RootFindingBuffer{T <: Real, U <: Number} <: AbstractBuffer
    jd0::T
    rv::RetAlloc{Taylor1{U}}
    teph::EphemerisEvaluationBuffer{T, U}
    g_dg_val::Vector{U}
    g_dg::Vector{Taylor1{U}}
    g_constant::Vector{Taylor1{T}}
    f_tupl::RootFindingEvent{U}
    g_tupl::RootFindingEvent{U}
    f_tupl_old::RootFindingEvent{U}
    g_tupl_old::RootFindingEvent{U}
end

"""
    ImpactMonitoringBuffer{T <: Real, U <: Number} <: AbstractBuffer

Pre-allocated memory for [`closeapproaches`](@ref).

# Fields

- `prop::PropagationBuffer{T, U, T}`: propagation buffer.
- `root::RootFindingBuffer{T, U}`: root-finding buffer.
"""
struct ImpactMonitoringBuffer{T <: Real, U <: Number} <: AbstractBuffer
    prop::PropagationBuffer{T, U, T}
    root::RootFindingBuffer{T, U}
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
    # Root finding buffer
    @unpack t, x = prop.cache
    y, z = zero(x[1]), zero(q0[1])
    rv = RetAlloc{Taylor1{U}}(
        [zero(y) for _ in 1:12],
        [[zero(y) for _ in 1:6]],
        [Array{Taylor1{U}, 2}(undef, 0, 0)],
        [Array{Taylor1{U}, 3}(undef, 0, 0, 0)],
        [Array{Taylor1{U}, 4}(undef, 0, 0, 0, 0)]
    )
    teph = EphemerisEvaluationBuffer(target.eph, tlim, order, q0)
    g_dg_val = [zero(z), zero(z)]
    g_dg = [zero(y), zero(y)]
    g_constant = [zero(t), zero(t)]
    f_tupl = RootFindingEvent{U}(false, zero(y))
    g_tupl = RootFindingEvent{U}(false, zero(y))
    f_tupl_old = RootFindingEvent{U}(false, zero(y))
    g_tupl_old = RootFindingEvent{U}(false, zero(y))
    root = RootFindingBuffer{T, U}(jd0, rv, teph,  g_dg_val, g_dg, g_constant,
        f_tupl, g_tupl, f_tupl_old, g_tupl_old)
    # Impact monitoring buffer
    return ImpactMonitoringBuffer{T, U}(prop, root)
end