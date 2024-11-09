@doc raw"""
    Tracklet{T <: AbstractFloat}

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
- `idxs::Vector{Int}`: indices of original `radec` that entered the tracklet.
"""
@auto_hash_equals struct Tracklet{T <: AbstractFloat}
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
    function Tracklet{T}(radec::Vector{RadecMPC{T}}, observatory::ObservatoryMPC{T},
                                 night::TimeOfDay, date::DateTime, α::T, δ::T, v_α::T, v_δ::T,
                                 mag::T, nobs::Int, indices::Vector{Int}) where {T <: AbstractFloat}

        return new{T}(radec, observatory, night, date, α, δ, v_α, v_δ, mag, nobs, indices)
    end
end

# Functions to get specific fields of a Tracklet object
date(x::Tracklet) = x.date
ra(x::Tracklet) = x.α
dec(x::Tracklet) = x.δ
observatory(x::Tracklet) = x.observatory
vra(x::Tracklet) = x.v_α
vdec(x::Tracklet) = x.v_δ
mag(x::Tracklet) = x.mag
nobs(x::Tracklet) = x.nobs
indices(x::Tracklet) = x.indices
indices(x::AbstractVector{Tracklet{T}}) where {T <: AbstractFloat} =
    sort!(reduce(vcat, indices.(x)))

# Print method for Tracklet{T}
# Examples:
# 3 observation tracklet around 2023-11-18T19:59:55.392 at WFST, Lenghu
# 4 observation tracklet around 2023-11-22T08:07:46.336 at Mt. Lemmon Survey
function show(io::IO, x::Tracklet{T}) where {T <: AbstractFloat}
    print(io, x.nobs, " observation tracklet around ", x.date, " at ", x.observatory.name)
end

# Order in Tracklet is given by date
isless(a::Tracklet{T}, b::Tracklet{T}) where {T <: AbstractFloat} = a.date < b.date

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
function Tracklet(radec::Vector{RadecMPC{T}}, df::AbstractDataFrame) where {T <: AbstractFloat}
    # Defining quantities of a Tracklet
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
        return Tracklet{T}(radec, observatory, night, date, α,
                                   δ, v_α, v_δ, mag, nobs, indices)
    end

    # Make sure there are no repeated dates
    gdf = groupby(df, :date)
    df = combine(gdf, [:α, :δ, :mag] .=> mean, renamecols = false)

    # Julian days of observation
    df.t_julian = dtutc2jdtdb.(df.date)
    # Days of observation [relative to first observation]
    df.t_rel = df.t_julian .- df.t_julian[1]
    # Mean date [relative to first observation]
    t_mean = mean(df.t_rel)
    # Mean date [DateTime]
    date = jdtdb2dtutc(df.t_julian[1] + t_mean)

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

    return Tracklet{T}(radec, observatory, night, date, α, δ, v_α,
                               v_δ, mag, nobs, indices)
end

@doc raw"""
    reduce_tracklets(radec::Vector{RadecMPC{T}}) where {T <: AbstractFloat}

Return a `Vector{Tracklet{T}}` where each element corresponds to a batch of
observations taken by the same observatory on the same night. The reduction
is performed via polynomial regression.
"""
function reduce_tracklets(radec::Vector{RadecMPC{T}}) where {T <: AbstractFloat}
    # Convert to DataFrame
    df = DataFrame(radec)
    # Compute TimeOfDay
    df.TimeOfDay = TimeOfDay.(radec)
    # Group by observatory and TimeOfDay
    gdf = groupby(df, [:observatory, :TimeOfDay])
    # Allocate memmory for tracklets vector
    tracklets = Vector{Tracklet{T}}(undef, gdf.ngroups)
    # Reduce tracklets
    Threads.@threads for i in 1:gdf.ngroups
        rows = getfield(gdf[i], :rows)
        tracklets[i] = Tracklet(radec[rows], gdf[i])
    end
    # Sort by date
    sort!(tracklets)

    return tracklets
end

@doc raw"""
    curvature(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}

Return:
- A vector with the geodesic curvature and along track acceleration of `radec`.
- The curvature covariance matrix.

!!! reference
    See section 5.1 of https://doi.org/10.1016/j.icarus.2007.11.033
"""
function curvature(ts::AbstractVector{T}, αs::AbstractVector{T}, δs::AbstractVector{T},
    σs::AbstractVector{T}) where {T <: Real}
    @assert length(ts) == length(αs) == length(δs) == length(σs) ≥ 3 """
    At least three observations needed for significant curvature computation."""
    # Days of observation [relative to first observation]
    ts = ts .- ts[1]
    # Weights
    wt = 1 ./ σs .^ 2
    # Initial guess for coefficients
    p0 = ones(T, 3)
    # Fit a second order polynomial to ra/dec
    fit_α = curve_fit(polymodel, ts, αs, wt, p0)
    fit_δ = curve_fit(polymodel, ts, δs, wt, p0)
    # Ra/dec covariance matrix
    Γ_α = vcov(fit_α)
    Γ_δ = vcov(fit_δ)
    Γ_αδ = [Γ_α zeros(T, 3, 3); zeros(T, 3, 3) Γ_δ]
    # Mean time of observation
    t_mean = mean(ts)
    # Evaluate ra/dec and its derivatives at mean time
    α = mod2pi(polymodel(t_mean, fit_α.param))
    δ = polymodel(t_mean, fit_δ.param)
    v_α = polymodel(t_mean, diffcoeffs(fit_α.param))
    v_δ = polymodel(t_mean, diffcoeffs(fit_δ.param))
    a_α = polymodel(t_mean, diffcoeffs(diffcoeffs(fit_α.param)))
    a_δ = polymodel(t_mean, diffcoeffs(diffcoeffs(fit_δ.param)))
    # Trigonometric functions
    sin_δ, cos_δ = sincos(δ)
    # Proper motion
    η = sqrt(v_α^2 * cos_δ^2 + v_δ^2)
    # Geodesic curvature
    κ = ((a_δ*v_α - a_α*v_δ) * cos_δ + v_α * (η^2 + v_δ^2) * sin_δ) / η^3
    # Along track acceleration
    v_η = (a_α * v_α * cos_δ^2 + a_δ * v_δ - v_α^2 * v_δ * cos_δ * sin_δ ) / η
    # Angles to curvature jacobian
    J1 = -(-2 * v_α^3 * cos_δ^2 * sin_δ * a_δ + sin_δ * a_δ * v_α * v_δ^2 +
            2 * v_α^2 * cos_δ^2 * sin_δ * a_α * v_δ - sin_δ * a_α * v_δ^3 -
            v_α^5 * cos_δ^3 - 4 * v_α^3 * cos_δ * v_δ^2 + v_α^3 * cos_δ^3 * v_δ^2 -
            2 * v_α * cos_δ * v_δ^4) / η^5
    J2 = -v_α * (sin(2*δ) * (v_α^2 * a_α * cos_δ^2 + 2 * v_δ^2 * a_α - v_α * v_δ * a_δ) +
                 2 * v_α * v_δ^3 * cos(2δ) + 2 * v_α^3 * v_δ * cos_δ^4) / (2 * η^3)
    J3 = (-v_α * cos_δ^3 * (2v_α * a_δ - 3 * v_δ * a_α) + v_δ^2 * (a_δ * cos_δ -
           v_α^2 * cos_δ^2 * sin_δ + 2 * v_δ^2 * sin_δ) ) / η^5
    J4 = - cos_δ * v_δ * (-cos_δ * a_α * v_δ + v_α^3 * sin_δ * cos_δ^2 +
                           2 * v_α * sin_δ * v_δ^2 + cos_δ * a_δ * v_α) / η^3
    J5 = -(cos_δ * (v_α^2 * a_α * cos_δ^2 - 2 * v_δ^2 * a_α + 3 * v_α * v_δ * a_δ) -
                    v_α * v_δ * sin_δ * (v_α^2 * cos_δ^2 - 2 * v_δ^2)) / η^5
    J6 = - v_α * cos_δ^2 * (-a_δ * v_α + v_α^3 * cos_δ * sin_δ + a_α * v_δ) / η^3
    J7 = - v_δ * cos_δ / η^3
    J8 = v_α * cos_δ^2 / η
    J9 = v_α * cos_δ / η^3
    J10 = v_δ / η
    J = [
        0 J1 J3 J5 J7 J9;
        0 J2 J4 J6 J8 J10
    ]
    # Curvature covariance matrix
    Γ_C = J * Γ_αδ * J'

    return [κ, v_η], Γ_C
end

function curvature(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
    # Days of observation [julian days]
    ts = dtutc2jdtdb.(date.(radec))
    # Right ascension
    αs = ra.(radec)
    # Declination
    δs = dec.(radec)
    # Standard deviations [rad]
    σs = arcsec2rad.(σsveres17.(radec))

    return curvature(ts, αs, δs, σs)
end