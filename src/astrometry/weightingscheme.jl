"""
    AbstractWeightingScheme{T} <: AbstractAstrometryErrorModel{T}

Supertype for the optical astrometry weighting schemes interface.

Every weighting scheme `W{T}` must:
- be a mutable subtype of `AbstractWeightingScheme{T}`,
- implement a `W(::AbstractOpticalVector{T})` constructor,
- override `getid(::W)`, `covariances(::W)` and `outliers(::W)`,
- override `update!(::W{T}, ::AbstractOpticalVector{T})`.

`covariances(::W)` must return a `Vector{SMatrix{2, 2, T, 4}}` with the covariance
matrix of the measurement errors of each observation, in the `(α cos(δ), δ)`
coordinates and units of arcsec², while `outliers(::W)` must return a `BitVector`
flagging the observations to be excluded a priori. Observations whose covariance
matrix is not positive definite must be flagged as outliers.

See also [`radeccovariance`](@ref).
"""
abstract type AbstractWeightingScheme{T} <: AbstractAstrometryErrorModel{T} end

# AbstractWeightingScheme interface
nobs(x::AbstractWeightingScheme) = length(covariances(x))

covariances(x::AbstractWeightingScheme, i::Int) = covariances(x)[i]
outliers(x::AbstractWeightingScheme, i::Int) = outliers(x)[i]

covariances(x::AbstractWeightingScheme, i::AbstractVector{Int}) = view(covariances(x), i)
outliers(x::AbstractWeightingScheme, i::AbstractVector{Int}) = view(outliers(x), i)

setoutlier!(x::AbstractWeightingScheme, i::Int, y::Bool) = outliers(x)[i] = y
setoutlier!(x::AbstractWeightingScheme, i::AbstractVector{Int}, y::AbstractVector{Bool}) =
    outliers(x, i) .= y

# Print method for AbstractWeightingScheme
show(io::IO, x::AbstractWeightingScheme) = print(io, getid(x),
    " weighting scheme with ", nobs(x), " observations")

"""
    radeccovariance(σα, σδ [, ρ])

Return the covariance matrix [arcsec²] of a pair of `(α cos(δ), δ)` measurement
errors with standard deviations `σα`, `σδ` [arcsec] and correlation `ρ`
(default: `0`).
"""
function radeccovariance(σα::T, σδ::T, ρ::T = zero(T)) where {T <: Real}
    c = ρ * σα * σδ
    return SMatrix{2, 2, T, 4}(σα^2, c, c, σδ^2)
end

# Whether a 2x2 (covariance or weight) matrix is symmetric positive definite
# Note: a matrix with NaN entries is never positive definite
isposdefcov(Γ::SMatrix{2, 2}) = Γ[1, 2] == Γ[2, 1] && Γ[1, 1] > 0 &&
    Γ[1, 1] * Γ[2, 2] - Γ[1, 2]^2 > 0

# Flag as outliers the observations with a non positive definite covariance matrix
covoutliers(flags::AbstractVector{Bool}, Γs::AbstractVector{<:SMatrix{2, 2}}) =
    BitVector(flags) .| .!isposdefcov.(Γs)

"""
    weightmatrix(Γ)

Return the weight matrix [arcsec⁻²], i.e. the inverse, of a 2x2 covariance matrix
`Γ` [arcsec²]. If `Γ` is not positive definite, return a zero matrix.

See also [`radeccovariance`](@ref).
"""
function weightmatrix(Γ::SMatrix{2, 2, T, 4}) where {T <: Real}
    !isposdefcov(Γ) && return SMatrix{2, 2, T, 4}(zero(T), zero(T), zero(T), zero(T))
    # Note: these expressions reduce to 1/Γ[1, 1] and 1/Γ[2, 2]
    # for uncorrelated errors
    w11 = 1 / (Γ[1, 1] - Γ[1, 2]^2 / Γ[2, 2])
    w22 = 1 / (Γ[2, 2] - Γ[1, 2]^2 / Γ[1, 1])
    w12 = -Γ[1, 2] / (Γ[1, 1] * Γ[2, 2] - Γ[1, 2]^2)
    return SMatrix{2, 2, T, 4}(w11, w12, w12, w22)
end

"""
    UniformWeights{T} <: AbstractWeightingScheme{T}

Uniform optical astrometry weighting scheme: every observation has
uncorrelated 1 arcsec uncertainties in both coordinates.
"""
mutable struct UniformWeights{T} <: AbstractWeightingScheme{T}
    covariances::Vector{SMatrix{2, 2, T, 4}}
    outliers::BitVector
end

# Constructor
function UniformWeights(optical::AbstractOpticalVector{T}) where {T <: Real}
    covariances = covsuniform(optical)
    outliers = covoutliers(isdeprecated.(optical), covariances)
    return UniformWeights{T}(covariances, outliers)
end

# Override getid, covariances and outliers
getid(::UniformWeights) = "Uniform"
covariances(x::UniformWeights) = x.covariances
outliers(x::UniformWeights) = x.outliers

# Override update!
function update!(x::UniformWeights{T}, optical::AbstractOpticalVector{T}) where {T <: Real}
    x.covariances = covsuniform(optical)
    x.outliers = covoutliers(isdeprecated.(optical), x.covariances)
    return nothing
end

covsuniform(optical::AbstractOpticalVector{T}) where {T <: Real} =
    [radeccovariance(one(T), one(T)) for _ in eachindex(optical)]

"""
    SourceWeights{T} <: AbstractWeightingScheme{T}

Source optical astrometry weighting scheme: the uncertainties and
correlations are extracted from the observations' format.
"""
mutable struct SourceWeights{T} <: AbstractWeightingScheme{T}
    covariances::Vector{SMatrix{2, 2, T, 4}}
    outliers::BitVector
end

# Constructor
function SourceWeights(optical::AbstractOpticalVector{T}) where {T <: Real}
    covariances = covssource(optical)
    outliers = covoutliers(isoutlier.(optical), covariances)
    return SourceWeights{T}(covariances, outliers)
end

# Override getid, covariances and outliers
getid(::SourceWeights) = "Source"
covariances(x::SourceWeights) = x.covariances
outliers(x::SourceWeights) = x.outliers

# Override update!
function update!(x::SourceWeights{T}, optical::AbstractOpticalVector{T}) where {T <: Real}
    x.covariances = covssource(optical)
    x.outliers = covoutliers(isoutlier.(optical), x.covariances)
    return nothing
end

covssource(optical::AbstractOpticalVector{T}) where {T <: Real} =
    [radeccovariance(rms(x)..., corr(x)) for x in optical]

"""
    Veres17{T} <: AbstractWeightingScheme{T}

Veres et al. (2017) optical astrometry weighting scheme.

!!! reference
    See:
    - https://doi.org/10.1016/j.icarus.2017.05.021.
"""
mutable struct Veres17{T} <: AbstractWeightingScheme{T}
    covariances::Vector{SMatrix{2, 2, T, 4}}
    outliers::BitVector
end

# Constructor
function Veres17(optical::AbstractOpticalVector{T}) where {T <: Real}
    covariances = covsveres17(optical)
    outliers = covoutliers(isdeprecated.(optical), covariances)
    return Veres17{T}(covariances, outliers)
end

# Override getid, covariances and outliers
getid(::Veres17) = "Veres et al. (2017)"
covariances(x::Veres17) = x.covariances
outliers(x::Veres17) = x.outliers

# Override update!
function update!(x::Veres17{T}, optical::AbstractOpticalVector{T}) where {T <: Real}
    x.covariances = covsveres17(optical)
    x.outliers = covoutliers(isdeprecated.(optical), x.covariances)
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

# Return the covariance matrix [arcsec²] of each element of `optical`
# according to Veres et al. (2017)
# https://doi.org/10.1016/j.icarus.2017.05.021
function covsveres17(optical::AbstractOpticalVector{T}) where {T <: Real}
    σs = σsveres17.(optical)
    rex = rexveres17(optical)
    return [radeccovariance(r * σ, r * σ) for (r, σ) in zip(rex, σs)]
end

"""
    CustomWeights{T} <: AbstractWeightingScheme{T}

Optical astrometry weighting scheme with user-defined covariance matrices.

# Constructors

    CustomWeights(covariances [, outliers])
    CustomWeights(x::AbstractWeightingScheme [, idxs])

The first constructor takes a vector of `(α cos(δ), δ)` covariance matrices [arcsec²]
and an optional vector of outlier flags (default: none). The second one returns a copy
of the covariance matrices and outlier flags of another weighting scheme `x`,
optionally restricted to the observations with indices `idxs`. In both cases, the
observations with a non positive definite covariance matrix are flagged as outliers.

!!! warning
    `CustomWeights` cannot compute covariance matrices for a new set of observations;
    hence, `update!` throws an `ArgumentError` if the number of observations changes.

See also [`radeccovariance`](@ref).
"""
mutable struct CustomWeights{T} <: AbstractWeightingScheme{T}
    covariances::Vector{SMatrix{2, 2, T, 4}}
    outliers::BitVector
    # Inner constructor
    # Note: it replaces the default constructors, so every instance copies the
    # covariance matrices and flags the non positive definite ones as outliers
    function CustomWeights{T}(
            covariances::AbstractVector{SMatrix{2, 2, T, 4}},
            outliers::AbstractVector{Bool}
        ) where {T <: Real}
        @assert length(covariances) == length(outliers) "Number of covariance \
            matrices and outlier flags must be equal"
        Γs = collect(covariances)
        return new{T}(Γs, covoutliers(outliers, Γs))
    end
end

# Outer constructors
CustomWeights(
    covariances::AbstractVector{SMatrix{2, 2, T, 4}},
    outliers::AbstractVector{Bool} = falses(length(covariances))
) where {T <: Real} = CustomWeights{T}(covariances, outliers)

CustomWeights(x::AbstractWeightingScheme, idxs::AbstractVector{Int} = 1:nobs(x)) =
    CustomWeights(covariances(x, idxs), outliers(x, idxs))

# Equality and hash methods
# Note: `CustomWeights` is stored in orbits, whose equality (e.g. `iszero`)
# must compare the contents of the weighting scheme
Base.:(==)(x::CustomWeights, y::CustomWeights) = covariances(x) == covariances(y) &&
    outliers(x) == outliers(y)
Base.hash(x::CustomWeights, h::UInt) = hash(outliers(x),
    hash(covariances(x), hash(:CustomWeights, h)))

# Override getid, covariances and outliers
getid(::CustomWeights) = "Custom"
covariances(x::CustomWeights) = x.covariances
outliers(x::CustomWeights) = x.outliers

# Override update!
function update!(x::CustomWeights{T}, optical::AbstractOpticalVector{T}) where {T <: Real}
    length(optical) == nobs(x) || throw(ArgumentError("CustomWeights cannot compute \
        covariance matrices for a different number of observations"))
    return nothing
end
