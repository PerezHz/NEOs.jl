"""
    OpticalTracklet{T, O <: AbstractOpticalVector{T}} <: AbstractOpticalAstrometry{T}

A set of optical astrometry taken by the same observatory on the same night.

# Fields

- `optical::RefValue{O}`: reference to a vector of optical astrometry.
- `observatory::ObservatoryMPC{T}`: observing station.
- `night::TimeOfDay`: night of observation.
- `date::DateTime`: mean date of observation.
- `ra::T`: mean right ascension [rad].
- `dec::T`: mean declination [rad].
- `vra::T`: mean right ascension velocity [rad/day].
- `vdec::T`: mean declination velocity [rad/day].
- `mag::T`: mean apparent magnitude.
- `nobs::Int`: number of observations.
- `idxs::Vector{Int}`: indices of `optical` astrometry which are included in the tracklet.
"""
@auto_hash_equals fields = (date, ra, dec, observatory) struct OpticalTracklet{T, O <: AbstractOpticalVector{T}} <: AbstractOpticalAstrometry{T}
    optical::RefValue{O}
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

# AbstractAstrometryObservation interface
date(x::OpticalTracklet) = x.date
observatory(x::OpticalTracklet) = x.observatory
catalogue(::OpticalTracklet) = unknowncat()
rms(::OpticalTracklet{T}) where {T <: Real} = (T(NaN), T(NaN))
debias(::OpticalTracklet{T}) where {T <: Real} = (T(NaN), T(NaN))

nobs(x::OpticalTracklet) = x.nobs
nobs(x::AbstractVector{OpticalTracklet{T, O}}) where {T, O} =
    sum(nobs, x; init = 0)

indices(x::OpticalTracklet) = x.indices
indices(x::AbstractVector{OpticalTracklet{T, O}}) where {T, O} =
    sort!(reduce(vcat, indices.(x)))

astrometry(x::OpticalTracklet) = x.optical[][x.indices]
astrometry(x::AbstractVector{OpticalTracklet{T, O}}) where {T, O} =
    sort!(reduce(vcat, astrometry.(x)))

# Check if any observation in `t` has time `date`
function in(d::DateTime, t::OpticalTracklet)
    for x in view(t.optical, t.indices)
        d == date(x) && return true
    end
    return false
end

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
    print(io, x.nobs, " observation tracklet around ", x.date, " at ", x.observatory.name)
end

# Return the milliseconds between b.date and a.date
datediff(a::OpticalTracklet, b::OpticalTracklet) = (date(a) - date(b)).value

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

# Return the coefficients of the derivative of a polynomial with coefficients `x`
function diffcoeffs(x::AbstractVector{T}) where {T <: Real}
    y = Vector{T}(undef, length(x)-1)
    for i in eachindex(y)
        y[i] = i * x[i+1]
    end
    return y
end

# Outer constructor
function OpticalTracklet(optical::RefValue{<:AbstractOpticalVector{T}},
                         df::AbstractDataFrame) where {T <: Real}
    # Defining quantities of a Tracklet
    observatory = df.observatory[1]
    night = df.TimeOfDay[1]
    nobs = nrow(df)
    indices = getfield(df, :rows)
    # Only one observation
    if isone(nobs)
        date = df.date[1]
        ra, dec = df.ra[1], df.dec[1]
        vra, vdec = zero(T), zero(T)
        mag = df.mag[1]
        return OpticalTracklet(optical, observatory, night, date, ra,
            dec, vra, vdec, mag, nobs, indices)
    end
    # Make sure there are no repeated dates
    gdf = groupby(df, :date)
    df = combine(gdf, [:ra, :dec, :mag] .=> mean, renamecols = false)
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

    # Evaluate polynomials at mean date
    ra = mod2pi(polymodel(t_mean, ra_coef))
    dec = polymodel(t_mean, dec_coef)
    vra = polymodel(t_mean, diffcoeffs(ra_coef))
    vdec = polymodel(t_mean, diffcoeffs(dec_coef))

    # Mean apparent magnitude
    mag = mean(filter(!isnan, df.mag))

    return OpticalTracklet(optical, observatory, night, date, ra,
        dec, vra, vdec, mag, nobs, indices)
end

"""
    reduce_tracklets(::AbstractOpticalVector)

Return a vector of optical tracklets where each element corresponds to a
batch of observations taken by the same observatory on the same night.
The reduction is performed via polynomial regression.
"""
function reduce_tracklets(optical::O) where {T <: Real, O <: AbstractOpticalVector{T}}
    # Construct DataFrame
    df = DataFrame(date = date.(optical), ra = ra.(optical), dec = dec.(optical),
        observatory = observatory.(optical), mag = mag.(optical))
    # Compute TimeOfDay
    df.TimeOfDay = tmap(TimeOfDay, optical)
    # Group by observatory and TimeOfDay
    gdf = groupby(df, [:observatory, :TimeOfDay])
    # Reduce tracklets
    tracklets = Vector{OpticalTracklet{T, O}}(undef, gdf.ngroups)
    Threads.@threads for i in eachindex(tracklets)
        tracklets[i] = OpticalTracklet(Ref(optical), gdf[i])
    end
    # Sort by date
    sort!(tracklets)

    return tracklets
end