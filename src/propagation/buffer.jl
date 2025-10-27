# Pre-allocated memory for the evaluation of solar system
# ephemerides within propagation functions
struct EphemerisEvaluationBuffer{T <: Real, U <: Number} <: AbstractBuffer
    eph::DensePropagation2{T, T}
    aux::Vector{Taylor1{T}}
    ephT::Vector{Taylor1{T}}
    ephU::Vector{Taylor1{U}}
end

# This function assumes that order ≤ get_order(x)
reduceorder(x::Taylor1, order::Int) = Taylor1(x[0:order], order)

# Outer constructor
function EphemerisEvaluationBuffer(eph::DensePropagation2{T, T}, tlim::NTuple{2, T},
                                   order::Int, q0::Vector{U}) where {T <: Real, U <: Number}
    # Check order
    @assert order <= SSEPHORDER "order ($order) must be less or equal than solar system \
        ephemerides order ($SSEPHORDER)"
    # Load ephemeris
    t0, tf = minmax(tlim[1], tlim[2])
    j0 = searchsortedlast(eph.t, t0)
    jf = searchsortedfirst(eph.t, tf)
    t = eph.t[j0:jf]
    x = reduceorder.(view(eph.x, j0:jf-1, :), order)
    _eph_ = TaylorInterpolant(eph.t0, t, x)
    # Evaluation vectors
    zeroT = Taylor1(zeros(order+1), order)
    aux = [zero(zeroT) for _ in axes(x, 2)]
    ephT = [zero(zeroT) for _ in axes(x, 2)]
    zeroU = Taylor1(zero(q0[1]), order)
    ephU = [zero(zeroU) for _ in axes(x, 2)]

    return EphemerisEvaluationBuffer{T, U}(_eph_, aux, ephT, ephU)
end

# Evaluation methods for EphemerisEvaluationBuffer
function (y::EphemerisEvaluationBuffer{T, T})(t::Taylor1{T}) where {T <: Real}
    @unpack eph, aux, ephT, ephU = y
    # Get index of eph.x that interpolates at time t
    ind::Int, δt::Taylor1{T} = getinterpindex(eph, t)
    # Evaluate eph at t and convert the output to T
    Threads.@threads for i in eachindex(ephU)
        TS.zero!(ephT[i])
        TS.zero!(aux[i])
        TS._horner!(ephT[i], eph.x[ind, i], δt, aux[i])
        TS.zero!(ephU[i])
        for k in eachindex(ephU[i])
            TS.identity!(ephU[i], ephT[i], k)
        end
    end
    return ephU
end

function (y::EphemerisEvaluationBuffer{T, TaylorN{T}})(t::Taylor1{T}) where {T <: Real}
    @unpack eph, aux, ephT, ephU = y
    # Get index of eph.x that interpolates at time t
    ind::Int, δt::Taylor1{T} = getinterpindex(eph, t)
    # Evaluate eph at t and convert the output to TaylorN{T}
    Threads.@threads for i in eachindex(ephU)
        TS.zero!(ephT[i])
        TS.zero!(aux[i])
        TS._horner!(ephT[i], eph.x[ind, i], δt, aux[i])
        TS.zero!(ephU[i])
        for k in eachindex(ephU[i])
            ephU[i][k][0][1] = ephT[i][k]
        end
    end
    return ephU
end

function (y::EphemerisEvaluationBuffer{T, Taylor1{T}})(t::Taylor1{T}) where {T <: Real}
    @unpack eph, aux, ephT, ephU = y
    # Get index of eph.x that interpolates at time t
    ind::Int, δt::Taylor1{T} = getinterpindex(y.eph, t)
    # Evaluate eph at t and convert the output to Taylor1{T}
    Threads.@threads for i in eachindex(ephU)
        TS.zero!(ephT[i])
        TS.zero!(aux[i])
        TS._horner!(ephT[i], eph.x[ind, i], δt, aux[i])
        TS.zero!(ephU[i])
        for k in eachindex(ephU[i])
            ephU[i][k][0] = ephT[i][k]
        end
    end
    return ephU
end

# Parameters used within dynamical model functions
mutable struct DynamicalParameters{T <: Real, U <: Number, V <: Number}
    sseph::EphemerisEvaluationBuffer{T, U}
    acceph::Union{Nothing, EphemerisEvaluationBuffer{T, U}}
    poteph::Union{Nothing, EphemerisEvaluationBuffer{T, U}}
    jd0::V
    UJ_interaction::Union{Nothing, Vector{Bool}}
    N::Int
    μ::Vector{T}
    marsden_radial::NTuple{5, T}
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
    PropagationBuffer(dynamics, q0, jd0, tlim, params)

Return a `PropagationBuffer` object with pre-allocated memory for `propagate`.

# Arguments

- `dynamics`: dynamical model function.
- `q0::Vector{<:Number}`: vector of initial conditions.
- `jd0::Number`: initial Julian date (TDB).
- `tlim::NTuple{2, <:Real}`: ephemeris timespan [days since J2000].
- `params::Parameters{<:Real}`: see the `Propagation` section of [`Parameters`](@ref).
"""
function PropagationBuffer(dynamics::D, q0::Vector{U}, jd0::V, tlim::NTuple{2, T},
                           params::Parameters{T}) where {D, T <: Real, U <: Number,
                           V <: Number}
    # Unpack parameters
    @unpack order, μ_ast, maxsteps, parse_eqs, marsden_radial = params
    # Load solar system ephemeris [au, au/day]
    _sseph_ = EphemerisEvaluationBuffer(sseph, tlim, order, q0)
    # Number of massive bodies
    Nm1 = numberofbodies(sseph)
    # Number of bodies, including NEA
    N = Nm1 + 1
    # Vector of G*m values
    μ = convert(Vector{T}, vcat( μ_DE430[1:11], μ_ast[1:Nm1-11], zero(T) ) )
    # Check: number of SS bodies (N) in ephemeris must be equal to length of GM vector (μ)
    @assert N == length(μ) "Total number of bodies in ephemeris must be equal to length \
        of GM vector μ"
    # Accelerations, Newtonian potentials and interaction matrix with flattened bodies
    if dynamics == newtonian!
        _acceph_, _poteph_, UJ_interaction = nothing, nothing, nothing
    else
        _acceph_ = EphemerisEvaluationBuffer(acceph, tlim, order, q0)
        _poteph_ = EphemerisEvaluationBuffer(poteph, tlim, order, q0)
        UJ_interaction = falses(N)
        # Turn on Earth interaction
        UJ_interaction[ea] = true
    end
    # Dynamical parameters for `propagate`
    dparams = DynamicalParameters{T, U, V}(_sseph_, _acceph_, _poteph_, jd0, UJ_interaction,
                                           N, μ, marsden_radial)
    # TaylorIntegration cache
    cache = init_cache(Val(true), zero(T), q0, maxsteps, order, dynamics, dparams;
                       parse_eqs)

    return PropagationBuffer{T, U, V}(cache, dparams)
end