module DataFramesExt

using NEOs
using NEOs: AbstractAstrometry

import NEOs: RadecMPC, RadarJPL

if isdefined(Base, :get_extension)
    using DataFrames
    import DataFrames: DataFrame
else
    using ..DataFrames
    import ..DataFrames: DataFrame
end

# Method to convert a vector of observations to DataFrame
function DataFrame(obs::Vector{T}) where {T <: AbstractAstrometry}
    colnames = fieldnames(T)
    coltypes = fieldtypes(T)
    df = DataFrame([name => type[] for (name, type) in zip(colnames, coltypes)])
    for i in eachindex(obs)
        push!(df, [getfield(obs[i], colnames[j]) for j in eachindex(colnames)])
    end 

    return df 
end 

# Methods to convert a DataFrame to a vector of AbstractAstrometry
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