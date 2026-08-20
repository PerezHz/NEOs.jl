"""
    OpticalTracklet{T} <: AbstractOpticalAstrometry{T}

The average of a set of optical astrometry taken by the
same observatory on the same time of day.

# Fields

- `observatory::ObservatoryMPC{T}`: observing station.
- `timeofday::TimeOfDay`: observation time of day.
- `date::DateTime`: mean date of observation.
- `ra::T`: mean right ascension [rad].
- `dec::T`: mean declination [rad].
- `vra::T`: mean right ascension velocity [rad/day].
- `vdec::T`: mean declination velocity [rad/day].
- `mag::T`: mean apparent magnitude.
- `nobs::Int`: number of observations.
- `indices::Vector{Int}`: indices of the original optical
    astrometry vector which are included in the tracklet.
"""
@auto_hash_equals fields = (date, ra, dec, observatory) struct OpticalTracklet{T} <: AbstractOpticalAstrometry{T}
    observatory::ObservatoryMPC{T}
    timeofday::TimeOfDay
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
    return (tf - t0).value / daymillisec
end
=#

# Print methods for OpticalTracklet
show(io::IO, x::OpticalTracklet) = print(io, "Tracklet with ", nobs(x),
    " observation(s) around ", date(x), " at ", observatory(x).name)

function show(io::IO, ::MIME"text/plain", x::OpticalTracklet)
    t = repeat(' ', 4)
    print(io,
        nobs(x), "-observation Tracklet{", numtype(x), "}\n",
        t, rpad("Observatory: ", 21),  observatory(x).name, '\n',
        t, rpad("Date: ", 21),         date(x), '\n',
        t, rpad("Attributable: ", 21), "[",
            @sprintf("%.5f", rad2deg(ra(x))),   ", ",
            @sprintf("%.5f", rad2deg(dec(x))),  ", ",
            @sprintf("%.5f", rad2deg(vra(x))),  ", ",
            @sprintf("%.5f", rad2deg(vdec(x))), ", ",
            @sprintf("%.2f", mag(x)),
        "]",
    )
    return nothing
end

# Return the milliseconds between two dates
datediff(a::DateTime, b::DateTime) = (a - b).value
datediff(a::OpticalTracklet, b::OpticalTracklet) = datediff(date(a), date(b))

# Find the tracklet whose epoch is closest to t
closest_tracklet(t::Real, y::AbstractTrackletVector) =
    findmin(@. abs(t - dtutc2days(date(y))))[2]

# Normalized mean square residual for polynomial fit
function polynms(x::AbstractVector{T}, y::AbstractVector{T},
                 coeffs::AbstractVector{T}) where {T <: Real}
    χ = zero(T)
    @inbounds for i in eachindex(x, y)
        χ += abs2(y[i] - evalpoly(x[i], coeffs))
    end
    return χ / length(x)
end

# Fit a polynomial to points `(x, y)`. The order of the polynomial is increased
# until `polynms` is less than `tol`
function polyfit(x::AbstractVector{T}, y::AbstractVector{T};
                 tol::T = 1E-4) where {T <: Real}
    # x and y must be the same length and non-empty
    @assert length(x) == length(y) > 0
    # Vandermonde matrix
    V = Matrix{T}(undef, length(x), 7)
    for j in axes(V, 2)
        flag = isone(j)
        for i in axes(V, 1)
            if flag
                V[i, j] = one(T)
            else
                V[i, j] = V[i, j-1] * x[i]
            end
        end
    end
    # Avoid odd and high orders (to have a defined concavity and avoid overfit)
    for order in (1, 2, 4, 6)
        # Solve the linear least squares problem (V * coeffs = y) via QR factorization
        coeffs = view(V, :, 1:order+1) \ y
        # Convergence condition
        if polynms(x, y, coeffs) < tol || order == 6
            return coeffs
        end
    end
end

# Compute the mean of a vector, but skip NaNs
skipnanmean(x::AbstractVector) = mean(Iterators.filter(!isnan, x))

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
    @unpack date, ra, dec, observatory, mag, timeofday, indices = x
    # Number of observations
    nobs = minimum(length, x)
    # Zero or one observations
    if iszero(nobs)
        throw(ArgumentError("Zero observations do not constitute a tracklet"))
    elseif isone(nobs)
        return OpticalTracklet(
            observatory[1], timeofday[1], date[1], ra[1], dec[1],
            zero(ra[1]), zero(dec[1]), mag[1], nobs, collect(indices)
        )
    end
    # Make sure there are no repeated dates
    if !allunique(date)
        @warn "Two or or more observations have the same date"
    end
    # Observation times [Julian date UTC]
    jds = datetime2julian.(date)
    # Mean observation time [Julian days since first observation UTC]
    ts = jds .- jds[1]
    t_mean = mean(ts)
    # Mean observation date [UTC]
    d_mean = julian2datetime(jds[1] + t_mean)
    # Points in top/bottom quarter
    N_top = count(>(3π/2), ra)
    N_bottom = count(<(π/2), ra)
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
        α, δ = mod2pi(evalpoly(t_mean, ra_coef)), evalpoly(t_mean, dec_coef)
        h = skipnanmean(mag)
    end
    v_α = evalpoly(t_mean, diffcoeffs(ra_coef))
    v_δ = evalpoly(t_mean, diffcoeffs(dec_coef))

    return OpticalTracklet(observatory[i], timeofday[i], d_mean, α, δ, v_α, v_δ,
                           h, nobs, collect(indices))
end

"""
    reduce_tracklets(::AbstractOpticalVector)

Return a vector of optical tracklets where each element corresponds to a
batch of observations taken by the same observatory on the same time of day.
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
                timeofday = timeofday.(optical), trkid = trackletid.(optical),
                indices = eachindex(optical)
            )
            # Group by ...
            if hasfield($O, :trkid)
                gdf = groupby(df, [:trkid])
            else
                gdf = groupby(df, [:observatory, :timeofday])
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