"""
    AbstractODProblem{D, T <: Real}

Supertye for the orbit determination problems interface.
"""
abstract type AbstractODProblem{D, T <: Real} end

"""
    ODProblem{D, T,
              O <: AbstractOpticalVector{T},
              R <: Union{Nothing, AbstractRadarVector{T}},
              WT <: AbstractWeightingScheme{T},
              DT <: AbstractDebiasingScheme{T}} <: AbstractODProblem{D, T}

An orbit determination problem.

# Fields

- `dynamics::D`: dynamical model.
- `optical::O`: vector of optical astrometry.
- `tracklets::TrackletVector{T}`: vector of optical tracklets.
- `radar::R`: vector of radar astrometry.
- `weights::WT`: optical astrometry weighting scheme.
- `debias::DT`: optical astrometry debiasing scheme.
"""
mutable struct ODProblem{D, T,
                         O <: AbstractOpticalVector{T},
                         R <: Union{Nothing, AbstractRadarVector{T}},
                         WT <: AbstractWeightingScheme{T},
                         DT <: AbstractDebiasingScheme{T}} <: AbstractODProblem{D, T}
    dynamics::D
    optical::O
    tracklets::TrackletVector{T}
    radar::R
    weights::WT
    debias::DT
end

# Abbreviations
const OpticalODProblem{D, T, O} = ODProblem{D, T, O, Nothing, WT, DT} where {WT, DT}
const MixedODProblem{D, T, O, R} = ODProblem{D, T, O, R, WT, DT} where {WT, DT}

# Constructor
"""
    ODProblem(dynamics, optical [, radar]; kwargs...)

Return an orbit determination problem given a dynamical model `dynamics`,
a vector of optical astrometry `optical` and an optional vector of radar
astrometry `radar`.

# Keyword arguments

- `weights`: optical astrometry weighting scheme (default: `Veres17`).
- `debias`: optical astrometry debiasing scheme (default: `Eggl20`).
"""
function ODProblem(dynamics::D, optical::O, radar::R = nothing; weights::Type{WT} = Veres17,
                   debias::Type{DT} = Eggl20) where {D, T <: Real, O <: AbstractOpticalVector{T},
                   R <: Union{Nothing, AbstractRadarVector{T}}, WT <: AbstractWeightingScheme,
                   DT <: AbstractDebiasingScheme}
    # Group optical astrometry in tracklets
    tracklets = reduce_tracklets(optical)
    return ODProblem{D, T, O, R, WT{T}, DT{T}}(dynamics, optical, tracklets, radar,
                                               weights(optical), debias(optical))
end

# Print method for ODProblem
function show(io::IO, x::ODProblem)
    t = repeat(' ', 4)
    print(io,
        "Orbit determination problem\n",
        t, rpad("Dynamical model:", 21), x.dynamics, "\n",
        t, rpad("Optical astrometry:", 21), noptical(x), " observations of type ",
            opticaltype(x), "\n",
        t, rpad("Radar astrometry:", 21), hasradar(x) ? "$(nradar(x)) observations \
            of type $(radartype(x))" : "None", "\n",
        t, rpad("Weighting scheme:", 21), getid(x.weights), "\n",
        t, rpad("Debiasing scheme:", 21), getid(x.debias), "\n"
    )
end

# AbstractODProblem interface
scalartype(x::AbstractODProblem{D, T}) where {D, T} = T
opticaltype(x::ODProblem) = eltype(x.optical)
radartype(x::ODProblem) = hasradar(x) ? eltype(x.radar) : Nothing

hasradar(x::ODProblem) = !isnothing(x.radar)
optical(x::ODProblem) = x.optical
radar(x::ODProblem) = x.radar

noptical(x::ODProblem) = length(x.optical)
nradar(x::ODProblem) = hasradar(x) ? length(x.radar) : 0
nobs(x::ODProblem) = noptical(x) + nradar(x)

opticalindices(x::ODProblem) = eachindex(x.optical)
radarindices(x::MixedODProblem) = eachindex(x.radar)

opticaloutliers(x::ODProblem) = fill(false, noptical(x))
radaroutliers(x::MixedODProblem) = fill(false, nradar(x))

function minmaxdates(x::ODProblem)
    t0, tf = minmaxdates(x.optical)
    if hasradar(x)
        _t0_, _tf_ = minmaxdates(x.radar)
        t0, tf = min(t0, _t0_), max(tf, _tf_)
    end
    return t0, tf
end

function update!(x::AbstractODProblem{D, T},
                 optical::AbstractOpticalVector{T}) where {D, T <: Real}
    # Update ODProblem fields
    x.optical = optical
    x.tracklets = reduce_tracklets(optical)
    update!(x.weights, optical)
    update!(x.debias, optical)
    return nothing
end

# Special PropagationBuffer constructors
function PropagationBuffer(od::AbstractODProblem{D, T}, q0::Vector{<:Number}, jd0::Number,
                           k0::Int, kf::Int, params::Parameters{T}) where {D, T <: Real}
    t0 = dtutc2days(date(od.optical[k0]))
    tf = dtutc2days(date(od.optical[kf]))
    tlim = (t0 - params.bwdoffset, tf + params.fwdoffset)
    buffer = PropagationBuffer(od.dynamics, q0, jd0, tlim, params)

    return buffer
end

function PropagationBuffer(od::ODProblem, q0::Vector{<:Number}, jd0::Number,
                           params::Parameters)
    t0, tf = dtutc2days.(minmaxdates(od))
    tlim = (t0 - params.bwdoffset, tf + params.fwdoffset)
    buffer = PropagationBuffer(od.dynamics, q0, jd0, tlim, params)
    return buffer
end

function init_optical_residuals(::Type{U}, od::ODProblem, idxs = opticalindices(od),
                                outliers = opticaloutliers(od)) where {U <: Number}
    optical = view(od.optical, idxs)
    weights = view(od.weights.w8s, idxs)
    debias = view(od.debias.bias, idxs)
    return init_optical_residuals(U, optical, weights, debias, outliers)
end

function init_radar_residuals(::Type{U}, od::MixedODProblem, idxs = radarindices(od),
                              outliers = radaroutliers(od)) where {U <: Number}
    radar = view(od.radar, idxs)
    return init_radar_residuals(U, radar, outliers)
end

function set_od_order(::Type{T}, varorder::Int, numvars::Int = 6) where {T <: Real}
    if get_order() < varorder || get_numvars() != numvars
        set_variables(T, "dx"; order = varorder, numvars = numvars)
    end
    return nothing
end

function set_od_order(params::Parameters{T}, numvars::Int = 6) where {T <: Real}
    @unpack tsaorder, gaussorder, jtlsorder = params
    varorder = max(tsaorder, gaussorder, jtlsorder)
    set_od_order(T, varorder, numvars)
    return nothing
end