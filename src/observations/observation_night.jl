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
- `idxs::Vector{Int}`: indices of original `radec` that entered the night.
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
    indices::Vector{Int}
    # Inner constructor
    function ObservationNight{T}(radec::Vector{RadecMPC{T}}, observatory::ObservatoryMPC{T},
                                 night::TimeOfDay, date::DateTime, α::T, δ::T, v_α::T, v_δ::T,
                                 mag::T, nobs::Int, indices::Vector{Int}) where {T <: AbstractFloat}

        return new{T}(radec, observatory, night, date, α, δ, v_α, v_δ, mag, nobs, indices)
    end
end

# Functions to get specific fields of a ObservationNight object
date(x::ObservationNight{T}) where {T <: AbstractFloat} = x.date
ra(x::ObservationNight{T}) where {T <: AbstractFloat} = x.α
dec(x::ObservationNight{T}) where {T <: AbstractFloat} = x.δ
observatory(x::ObservationNight{T}) where {T <: AbstractFloat} = x.observatory
vra(x::ObservationNight{T}) where {T <: AbstractFloat} = x.v_α
vdec(x::ObservationNight{T}) where {T <: AbstractFloat} = x.v_δ
mag(x::ObservationNight{T}) where {T <: AbstractFloat} = x.mag
indices(x::ObservationNight{T}) where {T <: AbstractFloat} = x.indices

# Print method for ObservationNight{T}
# Examples:
# 3 observation night around 2023-11-18T19:59:55.392 at WFST, Lenghu
# 4 observation night around 2023-11-22T08:07:46.336 at Mt. Lemmon Survey
function show(io::IO, x::ObservationNight{T}) where {T <: AbstractFloat}
    print(io, x.nobs, " observation night around ", x.date, " at ", x.observatory.name)
end

# Order in ObservationNight is given by date
isless(a::ObservationNight{T}, b::ObservationNight{T}) where {T <: AbstractFloat} = a.date < b.date

# Evaluate a polynomial with coefficients p in every element of x
polymodel(x, p) = map(y -> evalpoly(y, p), x)

# Normalized mean square residual for polynomial fit
polyerror(x) = sum(x .^ 2) / length(x)

@doc raw"""
    polyfit(x::Vector{T}, y::Vector{T}; tol::T = 1e-4) where {T <: AbstractFloat}

Fit a polynomial to points `(x, y)`. The order of the polynomial is increased
until `polyerror` is less than `tol`.
"""
function polyfit(x::Vector{T}, y::Vector{T}; tol::T = 1e-4) where {T <: AbstractFloat}
    # Avoid odd and high orders (to have a defined concavity and avoid overfit)
    for order in [1, 2, 4, 6]
        # Initial guess for coefficients
        coeffs = ones(T, order+1)
        # Polynomial fit
        fit = curve_fit(polymodel, x, y, coeffs)
        # Convergence condition
        if polyerror(fit.resid) < tol || order == 6
            return fit.param
        end
    end
end

@doc raw"""
    diffcoeffs(x::Vector{T}) where {T <: AbstractFloat}

Return the coefficients of the derivative of a polynomial with coefficients `x`.
"""
function diffcoeffs(x::Vector{T}) where {T <: AbstractFloat}
    y = Vector{T}(undef, length(x)-1)
    for i in eachindex(y)
        y[i] = i * x[i+1]
    end
    return y
end

# Outer constructor
function ObservationNight(radec::Vector{RadecMPC{T}}, df::AbstractDataFrame) where {T <: AbstractFloat}
    # Defining quantities of an ObservationNight
    observatory = df.observatory[1]
    night = df.TimeOfDay[1]
    nobs = nrow(df)
    indices = getfield(df, :rows)
    # Only one observation
    if isone(nobs)
        date = df.date[1]
        α, δ = df.α[1], df.δ[1]
        v_α, v_δ = zero(T), zero(T)
        mag = df.mag[1]
        return ObservationNight{T}(radec, observatory, night, date, α,
                                   δ, v_α, v_δ, mag, nobs, indices)
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

    # Polynomial regression
    α_coef = polyfit(df.t_rel, df.α)
    δ_coef = polyfit(df.t_rel, df.δ)

    # Evaluate polynomials at mean date 
    α = mod2pi(polymodel(t_mean, α_coef))
    δ = polymodel(t_mean, δ_coef)
    v_α = polymodel(t_mean, diffcoeffs(α_coef))
    v_δ = polymodel(t_mean, diffcoeffs(δ_coef))

    # Mean apparent magnitude
    mag = mean(filter(!isnan, df.mag))

    return ObservationNight{T}(radec, observatory, night, date, α, δ, v_α,
                               v_δ, mag, nobs, indices)
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