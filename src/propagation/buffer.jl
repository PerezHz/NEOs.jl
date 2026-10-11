"""
    EphemerisEvaluationBuffer{T <: Real, U <: Number} <: AbstractBuffer

Pre-allocated memory for the evaluation of solar system ephemerides within
the dynamical models.

# Fields

- `eph::DensePropagation2{T, T}`: ephemeris restricted to the propagation timespan
    and to the columns (bodies' coordinates) used by the dynamical model.
- `forward::Bool`: whether `eph.t` is sorted in increasing order.
- `ind::RefValue{Int}`: index of the segment of `eph` used in the last evaluation.
- `active::Vector{Int}`: columns of `eph` evaluated at each call; the others are
    left untouched.
- `ephT::Vector{Taylor1{T}}`: evaluation of each column of `eph` as a `Taylor1{T}`.
- `ephU::Vector{Taylor1{U}}`: evaluation of each column of `eph` embedded into a
    `Taylor1{U}`.
"""
struct EphemerisEvaluationBuffer{T <: Real, U <: Number} <: AbstractBuffer
    eph::DensePropagation2{T, T}
    forward::Bool
    ind::RefValue{Int}
    active::Vector{Int}
    ephT::Vector{Taylor1{T}}
    ephU::Vector{Taylor1{U}}
end

# Outer constructor
function EphemerisEvaluationBuffer(
        eph::DensePropagation2{T, T}, tlim::NTuple{2, T}, order::Int,
        q0::Vector{U}; cols::AbstractVector{Int} = axes(eph.p, 2),
        active::AbstractVector{Int} = eachindex(cols)
    ) where {T <: Real, U <: Number}
    # Check order
    @assert order <= SSEPHORDER "order ($order) must be less or equal than solar system \
        ephemerides order ($SSEPHORDER)"
    # Check active columns
    @assert issubset(active, eachindex(cols)) "active columns must be indices of cols"
    # Load ephemeris
    t0, tf = minmax(tlim[1], tlim[2])
    j0 = searchsortedlast(eph.t, t0)
    jf = searchsortedfirst(eph.t, tf)
    _eph_ = TaylorSolution(eph.t[j0:jf], eph.p[j0:jf-1, cols])
    # Direction of the time grid
    forward = issorted(_eph_.t)
    @assert forward || issorted(_eph_.t, rev = true) "ephemeris times must be sorted \
        in increasing or decreasing order"
    # Evaluation vectors
    ephT = [Taylor1(zeros(T, order + 1)) for _ in cols]
    zeroU = Taylor1(zero(q0[1]), order)
    ephU = [zero(zeroU) for _ in cols]

    return EphemerisEvaluationBuffer{T, U}(_eph_, forward, Ref(firstindex(_eph_.t)),
        collect(active), ephT, ephU)
end

# Return the index of the segment of `y.eph` that contains (the constant term of) `t`
# and the time elapsed since the start of that segment.
function ephemeristime(y::EphemerisEvaluationBuffer{T}, t::Taylor1{T}) where {T <: Real}
    ts = y.eph.t
    x = constant_term(t)
    tmin, tmax = minmax(ts[end], ts[1])
    @assert tmin ≤ x ≤ tmax "Evaluation time outside range of interpolation"
    if length(ts) == 1
        ind = firstindex(ts)
    # Final time: evaluate the expansion of the last segment
    elseif x == ts[end]
        ind = lastindex(ts) - 1
    else
        ind = y.ind[]
        if !insegment(ts, ind, x, y.forward)
            if insegment(ts, ind + 1, x, y.forward)
                ind += 1
            elseif insegment(ts, ind - 1, x, y.forward)
                ind -= 1
            else
                ind = y.forward ? searchsortedlast(ts, x) :
                    searchsortedlast(ts, x, rev = true)
            end
        end
        y.ind[] = ind
    end
    return ind, x - ts[ind]
end

# Whether time `x` lies within the `i`-th segment of a time grid `ts` sorted in
# increasing (`forward = true`) or decreasing (`forward = false`) order
insegment(ts::AbstractVector, i::Int, x, forward::Bool) =
    firstindex(ts) ≤ i < lastindex(ts) &&
    (forward ? ts[i] ≤ x < ts[i+1] : ts[i] ≥ x > ts[i+1])

# Evaluation method for EphemerisEvaluationBuffer
function (y::EphemerisEvaluationBuffer{T, U})(t::Taylor1{T},
                                              threads::Bool = true) where {T <: Real, U <: Number}
    # Index of the segment of y.eph that contains t and time since its start
    ind, δt = ephemeristime(y, t)
    # Evaluate the active columns of y.eph at t and embed them into Taylor1{U}
    if threads
        _evaluate_ephemeris!(Val(true), y, ind, δt)
    else
        _evaluate_ephemeris!(Val(false), y, ind, δt)
    end
    return y.ephU
end

# Serial and multi-threaded evaluation loops of an EphemerisEvaluationBuffer
function _evaluate_ephemeris!(::Val{false}, y::EphemerisEvaluationBuffer, ind::Int, δt)
    for i in y.active
        _evaluate_ephemeris_component!(y, ind, δt, i)
    end
    return nothing
end

function _evaluate_ephemeris!(::Val{true}, y::EphemerisEvaluationBuffer, ind::Int, δt)
    Threads.@threads for i in y.active
        _evaluate_ephemeris_component!(y, ind, δt, i)
    end
    return nothing
end

# Evaluate the i-th column of an EphemerisEvaluationBuffer at time `δt` since the
# start of the `ind`-th segment.
function _evaluate_ephemeris_component!(y::EphemerisEvaluationBuffer, ind::Int, δt, i::Int)
    @inbounds taylorshift!(y.ephT[i], y.eph.p[ind, i], δt)
    @inbounds taylorembed!(y.ephU[i], y.ephT[i])
    return nothing
end

# Parameters used within dynamical model functions
mutable struct DynamicalParameters{T <: Real, U <: Number, V <: Number}
    N::Int
    jd0::V
    R_TP::T
    μ::Vector{T}
    marsden_radial::NTuple{5, T}
    sseph::EphemerisEvaluationBuffer{T, U}
    UJ_interaction::Union{Nothing, Vector{Bool}}
    acceph::Union{Nothing, EphemerisEvaluationBuffer{T, U}}
    poteph::Union{Nothing, EphemerisEvaluationBuffer{T, U}}
    orientAlloc::Union{Nothing, RetAlloc{Taylor1{T}}}
    Mmatrix::Array{Taylor1{U},3}
    zeroq1::Taylor1{U}
    threads::Bool
end

"""
    PropagationBuffer{T <: Real, U <: Number, V <: Number} <: AbstractBuffer

Pre-allocated memory for the propagation functions: [`propagate`](@ref),
[`propagate_root`](@ref) and [`propagate_lyap`](@ref).

# Fields

- `cache::VectorCache{...}`: `TaylorIntegration` cache.
- `dparams::DynamicalParameters{T, U, V}`: parameters used by the dynamical model.
"""
struct PropagationBuffer{T <: Real, U <: Number, V <: Number} <: AbstractBuffer
    cache::VectorCache{Vector{T}, Matrix{U}, Matrix{Taylor1{U}}, Vector{Taylor1{U}},
        Taylor1{T}, Vector{Taylor1{U}}, Vector{Taylor1{U}}, RetAlloc{Taylor1{U}}, Bool}
    dparams::DynamicalParameters{T, U, V}
end

"""
    PropagationBuffer(dynamics, q0, jd0, tlim, params; kwargs...)

Return a `PropagationBuffer` object with pre-allocated memory for `propagate`.

# Arguments

- `dynamics`: dynamical model function.
- `q0::Vector{<:Number}`: vector of initial conditions.
- `jd0::Number`: initial Julian date (TDB).
- `tlim::NTuple{2, <:Real}`: ephemeris timespan [days since J2000].
- `params::Parameters{<:Real}`: see the `Propagation` section of [`Parameters`](@ref).

# Keyword arguments

- `R_TP::Real`: target plane radius [au].
"""
function PropagationBuffer(
        dynamics::D, q0::Vector{U}, jd0::V, tlim::NTuple{2, T},
        params::Parameters{T}; R_TP::T = 0.2
    ) where {D, T <: Real, U <: Number, V <: Number}
    # Unpack parameters
    @unpack order, maxsteps, marsden_radial, parse_eqs, threads = params
    # Number of bodies (perturbers + asteroid)
    N = numberofbodies(Val(dynamics))
    # Gravitational parameters
    μ = T.(gm(Val(dynamics)))
    # Indices of perturbers
    idxs = bodyindices(Val(dynamics))
    cols = nbodyind(SSEPHNBODIES, idxs)
    # Solar system ephemeris [au, au/day]
    _sseph_ = EphemerisEvaluationBuffer(sseph, tlim, order, q0; cols)
    # Accelerations, Newtonian potentials and interaction matrix with flattened bodies
    if dynamics in (sunearthmoon!, newtonian!)
        _acceph_, _poteph_, UJ_interaction = nothing, nothing, nothing
        _orientation_ = nothing
    else
        _acceph_ = EphemerisEvaluationBuffer(acceph, tlim, order, q0;
            cols = view(cols, 1:3(N-1)))
        _poteph_ = EphemerisEvaluationBuffer(poteph, tlim, order, q0;
            cols = idxs)
        UJ_interaction = falses(N)
        # Turn on Earth interaction
        UJ_interaction[ea] = true
        _orientation_ = allocate_c2t_jpl_de430(Taylor1(order))
    end
    M_ = Array{Taylor1{U}}(undef, 3, 3, N)
    zero_q_1 = auxzero(Taylor1(q0[1], order))
    for j = 1:3
        for i = 1:3
            M_[i, j, ea] = zero(zero_q_1)
        end
    end
    # Dynamical parameters for `propagate`
    dparams = DynamicalParameters{T, U, V}(N, jd0, R_TP, μ, marsden_radial, _sseph_,
        UJ_interaction, _acceph_, _poteph_, _orientation_, M_, zero_q_1, threads)
    # TaylorIntegration cache
    cache = init_cache(Val(true), zero(T), q0, maxsteps, order, dynamics, dparams;
                       parse_eqs)

    return PropagationBuffer{T, U, V}(cache, dparams)
end
