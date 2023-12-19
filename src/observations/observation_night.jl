@doc raw"""
    ObservationNight{T <: AbstractFloat}

A set of optical observations taken by the same observatory on the same night.

# Fields

- `radec::Vector{RadecMPC{T}}`: vector of observations.
- `observatory::ObservatoryMPC{T}`: observing station.
- `night::TimeOfDay`: night of observation.
- `date::DateTime`: mean date of observation.
- `α::T`: mean right ascension.
- `δ::T`: mean declination.
- `v_α::T`: mean right ascension velocity.
- `v_δ::T`: mean declination velocity.
- `mag::T`: mean apparent magnitude.
- `nobs::Int`: number of observations.
"""
@auto_hash_equals struct ObservationNight{T <: AbstractFloat}
    radec::Vector{RadecMPC{T}}
    observatory::ObservatoryMPC{T}
    night::TimeOfDay
    date::DateTime
    α::T
    δ::T
    v_α::T
    v_δ::T
    mag::T
    nobs::Int
    # Inner constructor
    function ObservationNight{T}(radec::Vector{RadecMPC{T}}, observatory::ObservatoryMPC{T},
                                 night::TimeOfDay, date::DateTime, α::T, δ::T, v_α::T,
                                 v_δ::T, mag::T, nobs::Int) where {T <: AbstractFloat}

        return new{T}(radec, observatory, night, date, α, δ, v_α, v_δ, mag, nobs)
    end
end

# Print method for ObservationNight{T}
# Examples:
# 3 observation night around 2023-11-18T19:59:55.392 at WFST, Lenghu
# 4 observation night around 2023-11-22T08:07:46.336 at Mt. Lemmon Survey
function show(io::IO, x::ObservationNight{T}) where {T <: AbstractFloat}
    print(io, x.nobs, " observation night around ", x.date, " at ", x.observatory.name)
end

# Order in ObservationNight is given by date
isless(a::ObservationNight{T}, b::ObservationNight{T}) where {T <: AbstractFloat} = a.date < b.date

# Outer constructor
function ObservationNight(radec::Vector{RadecMPC{T}}, df::AbstractDataFrame) where {T <: AbstractFloat}
    # Defining quantities of an ObservationNight
    observatory = df.observatory[1]
    night = df.TimeOfDay[1]
    nobs = nrow(df)

    # Only one observation
    if isone(nobs)
        date = df.date[1]
        α, δ = df.α[1], df.δ[1]
        v_α, v_δ = zero(T), zero(T)
        mag = df.mag[1]
        return ObservationNight{T}(radec, observatory, night, date, α,
                                   δ, v_α, v_δ, mag, nobs)
    end 

    # Make sure there are no repeated dates
    gdf = groupby(df, :date)
    df = combine(gdf, [:α, :δ, :mag] .=> mean, renamecols = false)

    # Julian days of observation
    df.t_julian = datetime2julian.(df.date)
    # Days of observation [relative to first observation]
    df.t_rel = df.t_julian .- df.t_julian[1]
    # Mean date [relative to first observation]
    t_mean = mean(df.t_rel)
    # Mean date [DateTime]
    date = julian2datetime(df.t_julian[1] + t_mean)

    # Points in top quarter
    N_top = count(x -> x > 3π/2, df.α)
    # Points in bottom quarter
    N_bottom = count(x -> x < π/2, df.α)
    # Discontinuity
    if !iszero(N_top) && !iszero(N_bottom)
        df.α = map(x -> x < π ? x + 2π : x, df.α)
    end

    # Linear regression
    α_p = lm(@formula(α ~ t_rel), df)
    α_coef = coef(α_p)
    δ_p = lm(@formula(δ ~ t_rel), df)
    δ_coef = coef(δ_p)

    # Evaluate polynomials at mean date 
    α = mod2pi(α_coef[1] + α_coef[2] * t_mean)
    δ = δ_coef[1] + δ_coef[2] * t_mean
    v_α = α_coef[2]
    v_δ = δ_coef[2]

    # Mean apparent magnitude
    mag = mean(filter(!isnan, df.mag))

    return ObservationNight{T}(radec, observatory, night, date, α, δ, v_α,
                               v_δ, mag, nobs)
end 

@doc raw"""
    reduce_nights(radec::Vector{RadecMPC{T}}) where {T <: AbstractFloat}

Return a `Vector{ObservationNight{T}}` where each eleement corresponds to a batch of
observations taken by the same observatory on the same night. 
The reduction is performed via linear regression. 
"""
function reduce_nights(radec::Vector{RadecMPC{T}}) where {T <: AbstractFloat}
    # Convert to DataFrame 
    df = DataFrame(radec)
    # Compute TimeOfDay
    df.TimeOfDay = TimeOfDay.(radec)
    # Group by observatory and TimeOfDay 
    gdf = groupby(df, [:observatory, :TimeOfDay])
    # Allocate memmory for nights vector
    nights = Vector{ObservationNight{T}}(undef, gdf.ngroups)
    # Reduce nights
    Threads.@threads for i in 1:gdf.ngroups
        rows = getfield(gdf[i], :rows)
        nights[i] = ObservationNight(radec[rows], gdf[i])
    end
    # Sort by date
    sort!(nights)
    
    return nights
end