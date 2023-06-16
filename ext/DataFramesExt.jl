module DataFramesExt

using NEOs
import NEOs: RadecMPC, RadarJPL

if isdefined(Base, :get_extension)
    using DataFrames: DataFrame, nrow
    import Tables: istable, rowaccess, rows, schema, Schema
else
    using ..DataFrames: DataFrame, nrow
    import ..Tables: istable, rowaccess, rows, schema, Schema
end

# Methods to convert a Vector{RadecMPC{T}} to a DataFrame
istable(::Type{Vector{RadecMPC{T}}}) where {T <: AbstractFloat} = true
rowaccess(::Type{Vector{RadecMPC{T}}}) where {T <: AbstractFloat} = true
rows(x::Vector{RadecMPC{T}}) where {T <: AbstractFloat} = x
schema(::Vector{RadecMPC{T}}) where {T <: AbstractFloat} = Schema(fieldnames(RadecMPC{T}), Tuple{fieldtypes(RadecMPC{T})...})

# Methods to convert a Vector{RadarJPL{T}} to a DataFrame
istable(::Type{Vector{RadarJPL{T}}}) where {T <: AbstractFloat} = true
rowaccess(::Type{Vector{RadarJPL{T}}}) where {T <: AbstractFloat} = true
rows(x::Vector{RadarJPL{T}}) where {T <: AbstractFloat} = x
schema(::Vector{RadarJPL{T}}) where {T <: AbstractFloat} = Schema(fieldnames(RadarJPL{T}), Tuple{fieldtypes(RadarJPL{T})...})

# Methods to convert a DataFrame to a Vector{RadecMPC{T}} / Vector{RadarJPL{T}}
for T in (:(RadecMPC), :(RadarJPL))
    @eval begin
        function $T(df::DataFrame)
            @assert all(string.(fieldnames($T)) .== names(df))
            coltypes = eltype.(eachcol(df))
            numtypes = filter!(x -> x <: AbstractFloat, coltypes)
            @assert all(numtypes .== numtypes[1])
            S = numtypes[1]
            obs = Vector{$T{S}}(undef, nrow(df)) 
            for i in eachindex(obs)
                obs[i] = $T(values(df[i, :])...)
            end 

            return obs 
        end 
    end
end 

end # module