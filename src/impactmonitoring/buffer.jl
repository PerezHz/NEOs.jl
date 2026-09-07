"""
    LineOfVariationsBuffer{T <: Real} <: AbstractBuffer

Pre-allocated memory for [`lineofvariations`](@ref).

# Fields

- `t0::T`: reference epoch [TDB days since J2000].
- `sun::Vector{T}`: Sun barycentric cartesian state vector at `t0` [au, au/day].
- `scalings::Vector{T}`: covariance matrix scaling factors.
- `resTN::Vector{OpticalResidual{T, TaylorN{T}}}`: buffer for `TaylorN{T}` residuals.
- `bufferTN::PropresBuffer{T, TaylorN{T}, T}`: buffer for `TaylorN{T}` propagations.
"""
struct LineOfVariationsBuffer{T <: Real} <: AbstractBuffer
    t0::T
    sun::Vector{T}
    scalings::Vector{T}
    resTN::Vector{OpticalResidual{T, TaylorN{T}}}
    bufferTN::PropresBuffer{T, TaylorN{T}, T}
end

TaylorSeries.order(x::LineOfVariationsBuffer) =
    TaylorSeries.order(x.bufferTN.prop.cache.x[1][0])

"""
    LineOfVariationsBuffer(IM, lovorder, params)

Return a `LineOfVariationsBuffer` object with pre-allocated
memory for [`lineofvariations`](@ref).

# Arguments

- `IM::IMProblem`: impact monitoring problem.
- `lovorder::Int`: order of Taylor expansions wrt LOV index.
- `params::Parameters`: see the `Propagation` section of [`Parameters`](@ref).
"""
function LineOfVariationsBuffer(IM::AbstractIMProblem{D, T}, lovorder::Int,
                                params::Parameters{T}) where {D, T <: Real}
    # Unpack
    @unpack orbit = IM
    @unpack eph_su = params
    # Set jet transport order
    Ndof = dof(IM)
    set_od_order(T, lovorder, Ndof)
    # Refence epoch [julian date TDB]
    t0 = epoch(orbit)
    jd0 = t0 + PE.J2000
    # Sun's state vector at jd0
    sun = eph_su(t0)
    # Covariance matrix scaling factors
    scalings = fill(1E-8, 6)
    if Ndof == 9
        scalings = vcat(scalings, params.marsden_scalings...)
    end
    # Initial condition
    q00 = orbit()
    q0TN = q00 + sigmas(orbit) .* TaylorSeries.variables(T, lovorder)
    # Vectors of residuals
    resTN = init_optical_residuals(TaylorN{T}, IM)
    # Propagation and residuals buffers
    bufferTN = PropresBuffer(IM, q0TN, jd0, params)

    return LineOfVariationsBuffer{T}(t0, sun, scalings, resTN, bufferTN)
end

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
    CloseApproachesBuffer{T <: Real, U <: Number} <: AbstractBuffer

Pre-allocated memory for [`closeapproaches`](@ref).

# Fields

- `prop::PropagationBuffer{T, U, T}`: propagation buffer.
- `root::RootFindingBuffer{T, U}`: root-finding buffer.
"""
struct CloseApproachesBuffer{T <: Real, U <: Number} <: AbstractBuffer
    prop::PropagationBuffer{T, U, T}
    root::RootFindingBuffer{T, U}
end

"""
    CloseApproachesBuffer(IM, q0, nyears, params)

Return an `ImpactMonitoringBuffer` object with pre-allocated
memory for [`closeapproaches`](@ref).

# Arguments

- `IM::IMProblem`: impact monitoring problem.
- `q0::Vector{<:Number}`: initial condition.
- `nyears::Real`: number of years.
- `params::Parameters`: see the `Propagation` section of [`Parameters`](@ref).
"""
function CloseApproachesBuffer(
        IM::AbstractIMProblem{D, T}, q0::Vector{U},
        nyears::T, params::Parameters{T}
    ) where {D, T <: Real, U <: Number}
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
    # Close approaches buffer
    return CloseApproachesBuffer{T, U}(prop, root)
end

"""
    VirtualImpactorsBuffer{T <: Real} <: AbstractBuffer

Pre-allocated memory for [`verifyvirtualimpactor`](@ref).

# Fields

- `res::Vector{OpticalResidual{T, TaylorN{T}}}`: buffer for `TaylorN{T}`
    residuals.
- `prop::PropresBuffer{T, TaylorN{T}, T}`: buffer for `TaylorN{T}`
    propagations.
- `CAs::CloseApproachesBuffer{T, TaylorN{T}}`: buffer for `TaylorN{T}`
    close approaches search.
"""
struct VirtualImpactorsBuffer{T <: Real} <: AbstractBuffer
    res::Vector{OpticalResidual{T, TaylorN{T}}}
    prop::PropresBuffer{T, TaylorN{T}, T}
    CAs::CloseApproachesBuffer{T, TaylorN{T}}
end

"""
    VirtualImpactorsBuffer(IM, params)

Return a `VirtualImpactorsBuffer` object with pre-allocated
memory for [`verifyvirtualimpactor`](@ref).

# Arguments

- `IM::IMProblem`: impact monitoring problem.
- `params::Parameters`: see the `Propagation` section of [`Parameters`](@ref).
"""
function VirtualImpactorsBuffer(
        IM::AbstractIMProblem{D, T}, params::Parameters{T}
    ) where {D, T <: Real}
    # Unpack
    @unpack orbit = IM
    # Vector of residuals
    res = init_optical_residuals(TaylorN{T}, IM; iobs = true)
    # Propagation buffer
    jd0 = epoch(orbit) + PE.J2000
    q0 = orbit() + 1E-8 * sigmas(orbit) .* TaylorSeries.variables(T, 2)
    prop = PropresBuffer(IM, q0, jd0, params)
    # Close approaches buffer
    nyears = ( datetime2julian(DateTime(2100, 1, 1, 12)) - jd0 ) / yr
    CAs = CloseApproachesBuffer(IM, q0, nyears, params)
    # Virtual impactors buffer
    return VirtualImpactorsBuffer{T}(res, prop, CAs)
end