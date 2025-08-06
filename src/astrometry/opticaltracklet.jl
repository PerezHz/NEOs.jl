"""
    OpticalTracklet{T} <: AbstractOpticalAstrometry{T}

The averague of a set of optical astrometry taken by the
same observatory on the same night.

# Fields

- `observatory::ObservatoryMPC{T}`: observing station.
- `night::TimeOfDay`: night of observation.
- `date::DateTime`: mean date of observation.
- `ra::T`: mean right ascension [rad].
- `dec::T`: mean declination [rad].
- `vra::T`: mean right ascension velocity [rad/day].
- `vdec::T`: mean declination velocity [rad/day].
- `mag::T`: mean apparent magnitude.
- `nobs::Int`: number of observations.
- `idxs::Vector{Int}`: indices of the original optical astrometry
    vector which are included in the tracklet.
"""
@auto_hash_equals fields = (date, ra, dec, observatory) struct OpticalTracklet{T} <: AbstractOpticalAstrometry{T}
    observatory::ObservatoryMPC{T}
    night::TimeOfDay
    date::DateTime
    ra::T
    dec::T
    vra::T
    vdec::T
    mag::T
    nobs::Int
    indices::Vector{Int}
end

# Abbreviations
const TrackletVector{T} = Vector{OpticalTracklet{T}} where {T}
const AbstractTrackletVector{T} = AbstractVector{OpticalTracklet{T}} where {T}

# AbstractAstrometryObservation interface
date(x::OpticalTracklet) = x.date
observatory(x::OpticalTracklet) = x.observatory
catalogue(::OpticalTracklet) = unknowncat()
rms(::OpticalTracklet{T}) where {T <: Real} = (T(NaN), T(NaN))
debias(::OpticalTracklet{T}) where {T <: Real} = (T(NaN), T(NaN))
corr(::OpticalTracklet{T}) where {T <: Real} = T(NaN)

nobs(x::OpticalTracklet) = x.nobs
nobs(x::AbstractTrackletVector) = sum(nobs, x; init = 0)

indices(x::OpticalTracklet) = x.indices
indices(x::AbstractTrackletVector) = sort!(reduce(vcat, indices.(x)))
indices(x::AbstractTrackletVector, i::AbstractVector{Int}) = indices(view(x, i))

#=
# TO DO: rename this function (e.g. as `timerange`), as the
# AbstractAstrometryObservation interface already defines a method
# of `numberofdays` for an AbstractObservationVector
function numberofdays(trks::AbstractVector{OpticalTracklet{<:Real}})
    dates = Vector{Tuple{DateTime, DateTime}}(undef, length(trks))
    for i in eachindex(dates)
        dates[i] = extrema(date, view(trks[i].optical, trks[i].indices))
    end
    t0, tf = minimum(first, dates), maximum(last, dates)
    return (tf - t0).value / 86_400_000
end
=#

# Print method for OpticalTracklet
function show(io::IO, x::OpticalTracklet)
    print(io, nobs(x), " observation tracklet around ", date(x),
          " at ", observatory(x).name)
end

# Return the milliseconds between two dates
datediff(a::DateTime, b::DateTime) = (a - b).value
datediff(a::OpticalTracklet, b::OpticalTracklet) = datediff(date(a), date(b))

# Find the tracklet whose epoch is closest to t
closest_tracklet(t::Real, y::AbstractTrackletVector) =
    findmin(@. abs(t - dtutc2days(date(y))))[2]

# Evaluate a polynomial with coefficients p in every element of x
polymodel(x, p) = map(y -> evalpoly(y, p), x)

# Normalized mean square residual for polynomial fit
polyerror(x) = sum(x .^ 2) / length(x)

# Fit a polynomial to points `(x, y)`. The order of the polynomial is increased
# until `polyerror` is less than `tol`
function polyfit(x::AbstractVector{T}, y::AbstractVector{T};
                 tol::T = 1e-4) where {T <: Real}
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

# Compute the mean of a vector, but skip NaNs
skipnanmean(x::AbstractVector) = mean(filter(!isnan, x))

# Return the coefficients of the derivative of a polynomial with coefficients `x`
function diffcoeffs(x::AbstractVector{T}) where {T <: Real}
    y = Vector{T}(undef, length(x)-1)
    for i in eachindex(y)
        y[i] = i * x[i+1]
    end
    return y
end

# Outer constructor
function OpticalTracklet(df::AbstractDataFrame)
    # Number of observations
    nobs = nrow(df)
    # Indices
    indices = getfield(df, :rows)

    # Only one observation
    if isone(nobs)
        observatory = df.observatory[1]
        night = df.TimeOfDay[1]
        date = df.date[1]
        ra, dec = df.ra[1], df.dec[1]
        vra, vdec = zero(ra), zero(dec)
        mag = df.mag[1]
        return OpticalTracklet(observatory, night, date, ra, dec,
                               vra, vdec, mag, nobs, indices)
    end

    # Make sure there are no repeated dates
    gdf = groupby(df, :date)
    df = combine(gdf, [:ra, :dec] .=> mean, :observatory, :mag => skipnanmean,
                 :TimeOfDay, renamecols = false)
    # Julian days of observation
    df.t_julian = dtutc2jdtdb.(df.date)
    # Days of observation [relative to first observation]
    df.t_rel = df.t_julian .- df.t_julian[1]
    # Mean date [relative to first observation]
    t_mean = mean(df.t_rel)
    # Mean date [DateTime]
    date = jdtdb2dtutc(df.t_julian[1] + t_mean)

    # Points in top quarter
    N_top = count(x -> x > 3π/2, df.ra)
    # Points in bottom quarter
    N_bottom = count(x -> x < π/2, df.ra)
    # Discontinuity
    if !iszero(N_top) && !iszero(N_bottom)
        df.ra = map(x -> x < π ? x + 2π : x, df.ra)
    end

    # Polynomial regression
    ra_coef = polyfit(df.t_rel, df.ra)
    dec_coef = polyfit(df.t_rel, df.dec)

    # All observations are from a non Earth-fixed observatory
    if all(istwoliner, df.observatory)
        i = argmin(@. abs(datediff(date, df.date)))
        observatory = df.observatory[i]
        night = df.TimeOfDay[i]
        date = df.date[i]
        ra, dec = df.ra[i], df.dec[i]
        vra = polymodel(df.t_rel[i], diffcoeffs(ra_coef))
        vdec = polymodel(df.t_rel[i], diffcoeffs(dec_coef))
        mag = df.mag[i]
    else
        observatory = df.observatory[1]
        night = df.TimeOfDay[1]
        ra = mod2pi(polymodel(t_mean, ra_coef))
        dec = polymodel(t_mean, dec_coef)
        vra = polymodel(t_mean, diffcoeffs(ra_coef))
        vdec = polymodel(t_mean, diffcoeffs(dec_coef))
        mag = skipnanmean(df.mag)
    end

    return OpticalTracklet(observatory, night, date, ra, dec,
                           vra, vdec, mag, nobs, indices)
end

"""
    reduce_tracklets(::AbstractOpticalVector)

Return a vector of optical tracklets where each element corresponds to a
batch of observations taken by the same observatory on the same night.
The reduction is performed via polynomial regression.
"""
function reduce_tracklets(optical::AbstractOpticalVector{T}) where {T <: Real}
    # Construct DataFrame
    df = DataFrame(date = date.(optical), ra = ra.(optical), dec = dec.(optical),
        observatory = observatory.(optical), mag = mag.(optical))
    # Compute TimeOfDay
    df.TimeOfDay = tmap(TimeOfDay, optical)
    # Group by observatory and TimeOfDay
    if hasfield(eltype(optical), :trkid)
        df.trkid = getfield.(optical, :trkid)
        gdf = groupby(df, [:trkid])
    else
        gdf = groupby(df, [:observatory, :TimeOfDay])
    end
    # Reduce tracklets
    tracklets = TrackletVector{T}(undef, gdf.ngroups)
    Threads.@threads for i in eachindex(tracklets)
        tracklets[i] = OpticalTracklet(gdf[i])
    end
    # Sort by date
    sort!(tracklets)

    return tracklets
end