# Optical astrometry weighting schemes API
# All weighting schemes `W{T}` must:
# 1. be mutable subtypes of `AbstractWeightingScheme{T}`,
# 2. have a `w8s::Vector{Tuple{T, T}}` field,
# 3. implement a `W(::AbstractOpticalVector{T})` constructor,
# 4. override `getid(::W)`,
# 5. override `update!(::W{T}, ::AbstractOpticalVector{T})`.

"""
    AbstractWeightingScheme{T} <: AbstractAstrometryErrorModel{T}

Supertype for the optical astrometry weighting schemes API.
"""
abstract type AbstractWeightingScheme{T} <: AbstractAstrometryErrorModel{T} end

# Print method for AbstractWeightingScheme
show(io::IO, w::AbstractWeightingScheme) = print(io, getid(w),
    " weighting scheme with ", length(w.w8s), " observations")

"""
    UniformWeights{T} <: AbstractWeightingScheme{T}

Uniform optical astrometry weighting scheme.
"""
mutable struct UniformWeights{T} <: AbstractWeightingScheme{T}
    w8s::Vector{Tuple{T, T}}
end

# Constructor
function UniformWeights(optical::AbstractOpticalVector{T}) where {T <: Real}
    w8s = [(one(T), one(T)) for _ in eachindex(optical)]
    return UniformWeights{T}(w8s)
end

# Override getid
getid(::UniformWeights) = "Uniform"

# Override update!
function update!(w::UniformWeights{T}, optical::AbstractOpticalVector{T}) where {T <: Real}
    w.w8s = [(one(T), one(T)) for _ in eachindex(optical)]
    return nothing
end

"""
    SourceWeights{T} <: AbstractWeightingScheme{T}

Source optical astrometry weighting scheme.
"""
mutable struct SourceWeights{T} <: AbstractWeightingScheme{T}
    w8s::Vector{Tuple{T, T}}
end

# Constructors
function SourceWeights(optical::AbstractOpticalVector{T}) where {T <: Real}
    σs = rms.(optical)
    w8s = @. tuple(1 / first(σs), 1 / last(σs))
    return SourceWeights{T}(w8s)
end

# Override getid
getid(::SourceWeights) = "Source"

# Override update!
function update!(w::SourceWeights{T}, optical::AbstractOpticalVector{T}) where {T <: Real}
    σs = rms.(optical)
    w.w8s = @. tuple(1 / first(σs, 1 / last(σs)))
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
    w8s::Vector{Tuple{T, T}}
end

# Constructor
function Veres17(optical::AbstractOpticalVector{T}) where {T <: Real}
    return Veres17{T}(w8sveres17(optical))
end

# Override getid
getid(::Veres17) = "Veres et al. (2017)"

# Override update!
function update!(w::Veres17{T}, optical::AbstractOpticalVector{T}) where {T <: Real}
    w.w8s = w8sveres17(optical)
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
    df = DataFrame(observatory = observatory.(optical), TimeOfDay = TimeOfDay.(optical))
    # Group by observatory and TimeOfDay
    gdf = groupby(df, [:observatory, :TimeOfDay])
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