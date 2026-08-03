"""
    AbstractWeightingScheme{T} <: AbstractAstrometryErrorModel{T}

Supertype for the optical astrometry weighting schemes interface.

Every weighting scheme `W{T}` must:
- be a mutable subtype of `AbstractWeightingScheme{T}`,
- implement a `W(::AbstractOpticalVector{T})` constructor,
- override `getid(::W)`, `weights(::W)`, `corrs(::W)` and `outliers(::W)`.
- override `update!(::W{T}, ::AbstractOpticalVector{T})`.
"""
abstract type AbstractWeightingScheme{T} <: AbstractAstrometryErrorModel{T} end

# AbstractWeightingScheme interface
nobs(x::AbstractWeightingScheme) = length(weights(x))

weights(x::AbstractWeightingScheme, i::Int) = weights(x)[i]
corrs(x::AbstractWeightingScheme, i::Int) = corrs(x)[i]
outliers(x::AbstractWeightingScheme, i::Int) = outliers(x)[i]

weights(x::AbstractWeightingScheme, i::AbstractVector{Int}) = view(weights(x), i)
corrs(x::AbstractWeightingScheme, i::AbstractVector{Int}) = view(corrs(x), i)
outliers(x::AbstractWeightingScheme, i::AbstractVector{Int}) = view(outliers(x), i)

setoutlier!(x::AbstractWeightingScheme, i::Int, y::Bool) = outliers(x)[i] = y
setoutlier!(x::AbstractWeightingScheme, i::AbstractVector{Int}, y::AbstractVector{Bool}) =
    outliers(x, i) .= y

# Print method for AbstractWeightingScheme
show(io::IO, x::AbstractWeightingScheme) = print(io, getid(x),
    " weighting scheme with ", nobs(x), " observations")

"""
    UniformWeights{T} <: AbstractWeightingScheme{T}

Uniform optical astrometry weighting scheme.
"""
mutable struct UniformWeights{T} <: AbstractWeightingScheme{T}
    weights::Vector{NTuple{2, T}}
    corrs::Vector{T}
    outliers::BitVector
end

# Constructor
function UniformWeights(optical::AbstractOpticalVector{T}) where {T <: Real}
    weights = [(one(T), one(T)) for _ in eachindex(optical)]
    corrs = [zero(T) for _ in eachindex(optical)]
    outliers = isdeprecated.(optical)
    return UniformWeights{T}(weights, corrs, outliers)
end

# Override getid, weights, corrs and outliers
getid(::UniformWeights) = "Uniform"
weights(x::UniformWeights) = x.weights
corrs(x::UniformWeights) = x.corrs
outliers(x::UniformWeights) = x.outliers

# Override update!
function update!(x::UniformWeights{T}, optical::AbstractOpticalVector{T}) where {T <: Real}
    x.weights = [(one(T), one(T)) for _ in eachindex(optical)]
    x.corrs = [zero(T) for _ in eachindex(optical)]
    x.outliers = isdeprecated.(optical)
    return nothing
end

"""
    SourceWeights{T} <: AbstractWeightingScheme{T}

Source optical astrometry weighting scheme.
"""
mutable struct SourceWeights{T} <: AbstractWeightingScheme{T}
    weights::Vector{NTuple{2, T}}
    corrs::Vector{T}
    outliers::BitVector
end

# Constructors
function SourceWeights(optical::AbstractOpticalVector{T}) where {T <: Real}
    σs = rms.(optical)
    weights = @. tuple(1 / first(σs), 1 / last(σs))
    corrs = corr.(optical)
    outliers = isoutlier.(optical)
    return SourceWeights{T}(weights, corrs, outliers)
end

# Override getid, weights, corrs and outliers
getid(::SourceWeights) = "Source"
weights(x::SourceWeights) = x.weights
corrs(x::SourceWeights) = x.corrs
outliers(x::SourceWeights) = x.outliers

# Override update!
function update!(x::SourceWeights{T}, optical::AbstractOpticalVector{T}) where {T <: Real}
    σs = rms.(optical)
    x.weights = @. tuple(1 / first(σs), 1 / last(σs))
    x.corrs = corr.(optical)
    x.outliers = isoutlier.(optical)
    return nothing
end

"""
    Veres17{T} <: AbstractWeightingScheme{T}

Veres et al. (2017) optical astrometry weighting scheme.

!!! reference
    See:
    - https://doi.org/10.1016/j.icarus.2017.05.021.
"""
mutable struct Veres17{T} <: AbstractWeightingScheme{T}
    weights::Vector{NTuple{2, T}}
    corrs::Vector{T}
    outliers::BitVector
end

# Constructor
function Veres17(optical::AbstractOpticalVector{T}) where {T <: Real}
    weights = w8sveres17(optical)
    corrs = [zero(T) for _ in eachindex(optical)]
    outliers = isdeprecated.(optical)
    return Veres17{T}(weights, corrs, outliers)
end

# Override getid, weights, corrs and outliers
getid(::Veres17) = "Veres et al. (2017)"
weights(x::Veres17) = x.weights
corrs(x::Veres17) = x.corrs
outliers(x::Veres17) = x.outliers

# Override update!
function update!(x::Veres17{T}, optical::AbstractOpticalVector{T}) where {T <: Real}
    x.weights = w8sveres17(optical)
    x.corrs = [zero(T) for _ in eachindex(optical)]
    x.outliers = isdeprecated.(optical)
    return nothing
end

# Return the statistical uncertainty of `obs` according to Veres et al. (2017).
# https://doi.org/10.1016/j.icarus.2017.05.021
function σsveres17(obs::AbstractOpticalAstrometry{T}) where {T <: Real}

    obscode = observatory(obs).code
    dt_utc_obs = date(obs)
    catalogue = NEOs.catalogue(obs).code

    # Unit weight (arcseconds)
    w = one(T)
    # Table 2: epoch-dependent astrometric residuals
    if obscode == "703"
        return Date(dt_utc_obs) < Date(2014,1,1) ? w : 0.8w
    elseif obscode ∈ ("691", "291") # Spacewatch, Kitt Peak, Arizona
        return Date(dt_utc_obs) < Date(2003,1,1) ? 0.6w : 0.5w
    elseif obscode == "644"
        return Date(dt_utc_obs) < Date(2003,9,1) ? 0.6w : 0.4w
    # Table 3: most active CCD asteroid observers
    elseif obscode ∈ ("704", "C51", "J75")
        return w
    elseif obscode == "G96"
        return 0.5w
    elseif obscode ∈ ("F51", "F52") # Pan-STARRS 1 & 2, Haleakala, Hawaii
        return 0.2w
    elseif obscode ∈ ("G45", "608")
        return 0.6w
    elseif obscode == "699"
        return 0.8w
    elseif obscode ∈ ("D29", "E12")
        return 0.75w
    # Table 4:
    elseif obscode ∈ ("645", "673", "H01")
        return 0.3w
    elseif obscode ∈ ("J04", "K92", "K93", "Q63", "Q64", "V37", "W85", "W86", "W87",
        "K91", "E10", "F65") # Tenerife + Las Cumbres
        return 0.4w
    elseif obscode ∈ ("689", "950", "W84")
        return 0.5w
    # Applies only to program code assigned to M. Micheli
    #elseif obscode ∈ ("G83", "309")
    #    if catalogue ∈ ("q", "t") # "q"=>"UCAC-4", "t"=>"PPMXL"
    #        return 0.3w
    #    elseif catalogue ∈ ("U", "V") # Gaia-DR1, Gaia-DR2
    #        return 0.2w
    #    end
    elseif obscode ∈ ("Y28",)
        if catalogue ∈ ("t", "U", "V")
            return 0.3w
        else
            return w
        end
    elseif obscode ∈ ("568",)
        if catalogue ∈ ("o", "s") # "o"=>"USNO-B1.0", "s"=>"USNO-B2.0"
            return 0.5w
        elseif catalogue ∈ ("U", "V") # Gaia DR1, DR2
            return 0.1w
        elseif catalogue ∈ ("t",) #"t"=>"PPMXL"
            return 0.2w
        else
            return w
        end
    elseif obscode ∈ ("T09", "T12", "T14") && catalogue ∈ ("U", "V") # Gaia DR1, DR2
        return 0.1w
    elseif catalogue == ""
        return 1.5w
    elseif catalogue != ""
        return w
    else
        return w
    end
end

# Return the relax factor for each element of `optical` according to
# Veres et al. (2017), which mitigates unresolved systematic errors
# in observations taken on the same night by the same observatory.
# https://doi.org/10.1016/j.icarus.2017.05.021
function rexveres17(optical::AbstractOpticalVector{T}) where {T <: Real}
    # Construct DataFrame
    df = DataFrame(observatory = observatory.(optical), timeofday = timeofday.(optical))
    # Group by observatory and timeofday
    gdf = groupby(df, [:observatory, :timeofday])
    # Number of observations per tracklet
    cdf = combine(gdf, nrow)
    # Count observations in each group
    Nv = cdf[gdf.groups, :nrow]
    # Relaxation factor
    return map(N -> N > 4.0 ? sqrt(N/4.0) : 1.0, Nv)
end

# Return the statistical weight of each element of `optical` according
# to Veres et al. (2017)
# https://doi.org/10.1016/j.icarus.2017.05.021
function w8sveres17(optical::AbstractOpticalVector{T}) where {T <: Real}
    σs = σsveres17.(optical)
    rex = rexveres17(optical)
    return @. tuple(1 / (rex * σs), 1 / (rex * σs))
end