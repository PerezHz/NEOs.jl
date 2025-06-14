# Type tree for the minor bodies astrometry interface
abstract type AbstractAstrometry end

abstract type AbstractAstrometrySource <: AbstractAstrometry end

abstract type AbstractCatalogue <: AbstractAstrometry end

abstract type AbstractObservatory{T <: Real} <: AbstractAstrometry end

abstract type AbstractOpticalAstrometry{T <: Real} <: AbstractAstrometry end

ra(x::AbstractOpticalAstrometry) = x.ra
dec(x::AbstractOpticalAstrometry) = x.dec

# Order in AbstractOpticalAstrometry is given by date
isless(a::AbstractOpticalAstrometry, b::AbstractOpticalAstrometry) = date(a) < date(b)

abstract type AbstractRadarAstrometry{T <: Real} <: AbstractAstrometry end

# Order in AbstractRadarAstrometry is given by date
isless(a::AbstractRadarAstrometry, b::AbstractRadarAstrometry) = date(a) < date(b)

#=
# Methods to convert a Vector{<:AbstractAstrometry} to a DataFrame
istable(::Type{Vector{<:AbstractAstrometry}}) = true
rowaccess(::Type{Vector{<:AbstractAstrometry}}) = true
rows(x::Vector{<:AbstractAstrometry}) = x
schema(::Vector{T}) where {T <: AbstractAstrometry} = Schema(fieldnames(T), Tuple{fieldtypes(T)...})

# Methods to convert a DataFrame to a Vector{<:AbstractAstrometry}
function Vector{T}(df::AbstractDataFrame) where {T <: AbstractAstrometry}
    @assert all(String.(fieldnames(T)) .== names(df)) "`DataFrame` column names don't match `$T` fieldnames"
    @assert all(fieldtypes(T) .== eltype.(eachcol(df))) "`DataFrame` column types don't match `$T` fieldtypes"
    obs = Vector{T}(undef, nrow(df))
    for (i, row) in zip(eachindex(obs), eachrow(df))
        obs[i] = T(values(row)...)
    end
    return obs
end

convert(::Type{Vector{T}}, df::DataFrame) where {T <: AbstractAstrometry} = Vector{T}(df)
=#