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
vra(x::OpticalTracklet) = x.vra
vdec(x::OpticalTracklet) = x.vdec
observatory(x::OpticalTracklet) = x.observatory
catalogue(::OpticalTracklet) = unknowncat()
band(x::OpticalTracklet) = ' '
rms(::OpticalTracklet{T}) where {T <: Real} = (T(NaN), T(NaN))
debias(::OpticalTracklet{T}) where {T <: Real} = (T(NaN), T(NaN))
corr(::OpticalTracklet{T}) where {T <: Real} = T(NaN)

trackletid(x::OpticalTracklet) = ""

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
function OpticalTracklet(x::NamedTuple)
    # Unpack
    @unpack date, ra, dec, observatory, mag, night, indices = x
    # Number of observations
    nobs = minimum(length, x)
    # Zero or one observations
    if iszero(nobs)
        throw(ArgumentError("Zero observations do not constitute a tracklet"))
    elseif isone(nobs)
        return OpticalTracklet(
            observatory[1], night[1], date[1], ra[1], dec[1],
            zero(ra[1]), zero(dec[1]), mag[1], nobs, collect(indices)
        )
    end
    # Make sure there are no repeated dates
    if !allunique(date)
        @warn "Two or or more observations have the same date"
    end
    # Observation times [JDTDB]
    jds = dtutc2jdtdb.(date)
    # Observation times [Julian days since first observation]
    ts = jds .- jds[1]
    # Mean observation date [Julian days since first observation]
    t_mean = mean(ts)
    # Mean observation date [UTC]
    d_mean = jdtdb2dtutc(jds[1] + t_mean)
    # Points in top quarter
    N_top = count(x -> x > 3π/2, ra)
    # Points in bottom quarter
    N_bottom = count(x -> x < π/2, ra)
    # Discontinuity
    if !iszero(N_top) && !iszero(N_bottom)
        for (i, y) in enumerate(ra)
            ra[i] = y < π ? y + 2π : y
        end
    end
    # Polynomial regression
    ra_coef = polyfit(ts, ra)
    dec_coef = polyfit(ts, dec)
    # All observations are from a non Earth-fixed observatory
    if all(istwoliner, observatory)
        i = argmin(@. abs(datediff(d_mean, date)))
        t_mean, d_mean = ts[i], date[i]
        α, δ = ra[i], dec[i]
        h = mag[i]
    else
        i = 1
        α, δ = mod2pi(polymodel(t_mean, ra_coef)), polymodel(t_mean, dec_coef)
        h = skipnanmean(mag)
    end
    v_α = polymodel(t_mean, diffcoeffs(ra_coef))
    v_δ = polymodel(t_mean, diffcoeffs(dec_coef))

    return OpticalTracklet(observatory[i], night[i], d_mean, α, δ, v_α, v_δ,
                           h, nobs, collect(indices))
end

"""
    reduce_tracklets(::AbstractOpticalVector)

Return a vector of optical tracklets where each element corresponds to a
batch of observations taken by the same observatory on the same night.
The reduction is performed via polynomial regression.
"""
reduce_tracklets

for O in nameof.(subtypes(AbstractOpticalAstrometry))
    @eval begin
        function reduce_tracklets(optical::AbstractVector{$O{T}};
                                  threads::Bool = true) where {T <: Real}
            # Construct DataFrame
            df = DataFrame(
                date = date.(optical), ra = ra.(optical), dec = dec.(optical),
                observatory = observatory.(optical), mag = mag.(optical),
                night = tmap(TimeOfDay, TimeOfDay, optical), trkid = trackletid.(optical),
                indices = eachindex(optical)
            )
            # Group by ...
            if hasfield($O, :trkid)
                gdf = groupby(df, [:trkid])
            else
                gdf = groupby(df, [:observatory, :night])
            end
            # Reduce tracklets
            cdf = combine(gdf, AsTable(:) => OpticalTracklet => :tracklets; threads)
            # Sort by date
            sort!(cdf, :tracklets)
            # Return vector of tracklets
            tracklets::TrackletVector{T} = cdf.:tracklets

            return tracklets
        end
    end
end