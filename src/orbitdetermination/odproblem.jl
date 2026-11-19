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

Both keyword arguments accept either a type, which is then constructed from `optical`,
or an instance with the same number of observations as `optical` (e.g. a
[`CustomWeights`](@ref)).
"""
function ODProblem(
        dynamics::D, optical::O, radar::R = nothing;
        weights::Union{Type{<:AbstractWeightingScheme}, AbstractWeightingScheme} = Veres17,
        debias::Union{Type{<:AbstractDebiasingScheme}, AbstractDebiasingScheme} = Eggl20
    ) where {D, T <: Real, O <: AbstractOpticalVector{T},
             R <: Union{Nothing, AbstractRadarVector{T}}}
    # Group optical astrometry in tracklets
    tracklets = reduce_tracklets(optical)
    # Error model
    wsm, dsm = errormodel(weights, optical), errormodel(debias, optical)
    return ODProblem{D, T, O, R, typeof(wsm), typeof(dsm)}(dynamics, optical, tracklets,
                                                           radar, wsm, dsm)
end

# Note: the first argument of `errormodel` is not annotated as
# `Type{<:AbstractAstrometryErrorModel}` because the unbounded parameter of, e.g.,
# `ZeroDebiasing{T}` makes the `UnionAll` type `ZeroDebiasing` not a subtype of
# `AbstractAstrometryErrorModel`, whose parameter must be `<: Real`; anyway, the
# keyword arguments of `ODProblem` already restrict the admissible types
errormodel(x::Type, optical::AbstractOpticalVector) = x(optical)

function errormodel(x::AbstractAstrometryErrorModel{T},
                    optical::AbstractOpticalVector{T}) where {T <: Real}
    @assert nobs(x) == length(optical) "Error model and optical astrometry must have \
        the same number of observations"
    return x
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
dof(x::ODProblem) = dof(Val(x.dynamics))

hasradar(x::ODProblem) = !isnothing(x.radar)
optical(x::ODProblem) = x.optical
radar(x::ODProblem) = x.radar

noptical(x::ODProblem) = length(x.optical)
nradar(x::ODProblem) = hasradar(x) ? length(x.radar) : 0
nobs(x::ODProblem) = noptical(x) + nradar(x)

opticalindices(x::ODProblem) = eachindex(x.optical)
radarindices(x::MixedODProblem) = eachindex(x.radar)

covariances(x::ODProblem) = covariances(x.weights)
debias(x::ODProblem) = debias(x.debias)
outliers(x::ODProblem) = outliers(x.weights)

function minmaxdates(x::ODProblem)
    t0, tf = minmaxdates(x.optical)
    if hasradar(x)
        _t0_, _tf_ = minmaxdates(x.radar)
        t0, tf = min(t0, _t0_), max(tf, _tf_)
    end
    return t0, tf
end

"""
    meandate(::Union{ODProblem, AbstractOrbit})

Return the `DateTime` [UTC] corresponding to the
weighted mean optical observational epoch.

See also [`meanepoch`](@ref).

!!! reference
    See equation (15) of:
    - https://doi.org/10.1007/s10569-024-10225-z
"""
function meandate(x::ODProblem)
    T = scalartype(x)
    t = Vector{T}(undef, noptical(x))
    w = Vector{T}(undef, noptical(x))
    for (i, (obs, Γ)) in enumerate(zip(optical(x), covariances(x)))
        δ = dec(obs)
        t[i] = datetime2julian(obs) - PE.J2000
        w[i] = 1 / (Γ[1, 1] * cos(δ)^2 + Γ[2, 2])
    end
    tmean = mean(t, weights(w))
    return julian2datetime(tmean + PE.J2000)
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

function init_optical_residuals(
        ::Type{U}, od::ODProblem,
        idxs::AbstractVector{Int} = opticalindices(od),
    ) where {U <: Number}
    return init_optical_residuals(U, od.weights, od.debias, idxs)
end

function init_radar_residuals(
        ::Type{U}, od::MixedODProblem,
        idxs::AbstractVector{Int} = radarindices(od),
        outliers = falses(nradar(od))
    ) where {U <: Number}
    radar = view(od.radar, idxs)
    return init_radar_residuals(U, radar, outliers)
end

function set_od_order(::Type{T}, varorder::Int, numvars::Int = 6) where {T <: Real}
    if TaylorSeries.order() < varorder || get_numvars() != numvars
        TaylorSeries.variables!(T, "dx"; order = varorder, numvars = numvars)
    end
    return nothing
end

function set_od_order(params::Parameters{T}, numvars::Int = 6) where {T <: Real}
    @unpack mmovorder, gaussorder, jtlsorder = params
    varorder = max(mmovorder, gaussorder, jtlsorder)
    set_od_order(T, varorder, numvars)
    return nothing
end