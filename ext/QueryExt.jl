module QueryExt

using Dates: DateTime, Date 
using NEOs: ObservatoryMPC, RadecMPC, observatory, date, laplace_interpolation
import NEOs: reduce_nights

if isdefined(Base, :get_extension)
    using Query
else
    using ..Query
end

@doc raw"""
    reduce_nights(radec::Vector{RadecMPC{T}}) where {T <: AbstractFloat}

Return one observatory, date, right ascension and declination per observation night in `radec`. The reduction is performed 
via Laplace interpolation. 
"""
function reduce_nights(radec::Vector{RadecMPC{T}}) where {T <: AbstractFloat}
    # Observations per observatory 
    radec_per_obs = radec |> 
                    @groupby(observatory(_)) |> 
                    @filter(length(_) >= 3) .|> 
                    collect
    # Initialize output vectos 
    observatories = Vector{ObservatoryMPC{T}}(undef, 0)
    dates = Vector{DateTime}(undef, 0)
    α = Vector{T}(undef, 0)
    δ = Vector{T}(undef, 0)
    # Find observation nights per observatory 
    for i in eachindex(radec_per_obs)
        # Current observatory 
        obs = observatory(radec_per_obs[i][1])
        # Filter observation night 
        nights = radec_per_obs[i] |> 
                @groupby(Date(date(_))) |> 
                @filter(length(_) >= 3) |> 
                @map(laplace_interpolation(_)) |> 
                @filter( !isnan(_[2]) && !isnan(_[3]) ) |> 
                collect 
        # Add observation night 
        observatories = vcat(observatories, fill(obs, length(nights)))
        dates = vcat(dates, first.(nights))
        α = vcat(α, getindex.(nights, 2))
        δ = vcat(δ, last.(nights))
    end 

    return observatories, dates, α, δ
end 

end # module