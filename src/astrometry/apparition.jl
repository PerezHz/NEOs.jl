"""
    Apparition{
        T <: Real,
        O <: AbstractOpticalVector{T},
        R <: AbstractRadarVector{T}
    }

A set of optical and/or radar astrometry pertaining to the
same object such that the time between any two consecutive
observations is less than or equal to a given gap.

# Fields

- `optical::O`: optical astrometry vector.
- `radar::R`: radar astrometry vector.
"""
struct Apparition{
        T <: Real,
        O <: AbstractOpticalVector{T},
        R <: AbstractRadarVector{T}
    }
    optical::O
    radar::R
    # Inner constructor
    function Apparition{T, O, R}(optical::O, radar::R) where {T, O, R}
        @assert !isempty(optical) || !isempty(radar) "To construct \
            an apparition both `optical` and `radar` cannot be empty"
        return new{T, O, R}(optical, radar)
    end
end

# Outer constructors
Apparition(x::AbstractOpticalVector{T}) where {T} = Apparition(x, RadarJPL{T}[])
Apparition(x::AbstractRadarVector{T}) where {T} = Apparition(OpticalMPC80{T}[], x)
Apparition(x::AbstractOpticalVector{T}, y::AbstractRadarVector{T}) where {T <: Real} =
    Apparition{T, typeof(x), typeof(y)}(x, y)

# Abbreviations
const AbstractApparitionVector{T} = AbstractVector{Apparition{T, O, R}} where {O, R}

# Apparition interface
radartype(x::Apparition) = eltype(x.radar)
opticaltype(x::Apparition) = eltype(x.optical)
scalartype(::Apparition{T, O, R}) where {T, O, R} = T

nradar(x::Apparition) = length(x.radar)
nradar(x::AbstractApparitionVector) = sum(nradar, x)
noptical(x::Apparition) = length(x.optical)
noptical(x::AbstractApparitionVector) = sum(noptical, x)
nobs(x::Apparition) = noptical(x) + nradar(x)
nobs(x::AbstractApparitionVector) = sum(nobs, x)

radar(x::Apparition) = collect(x.radar)
radar(x::AbstractApparitionVector) = sort!(mapreduce(radar, vcat, x))
optical(x::Apparition) = collect(x.optical)
optical(x::AbstractApparitionVector) = sort!(mapreduce(optical, vcat, x))

numberofdays(x::Apparition) = numberofdays(x.optical, x.radar)
numberofdays(x::AbstractApparitionVector) = sum(numberofdays, x)

indices(x::Vector) = eachindex(x)
indices(x::SubArray) = first(x.indices)
radarindices(x::Apparition) = indices(x.radar)
radarindices(x::AbstractApparitionVector) = sort!(mapreduce(radarindices, vcat, x))
opticalindices(x::Apparition) = indices(x.optical)
opticalindices(x::AbstractApparitionVector) = sort!(mapreduce(opticalindices, vcat, x))

# Print methods for Apparition
show(io::IO, x::Apparition) = print(io, "Apparition with ", noptical(x), " optical and ",
    nradar(x), " radar observation(s)")

function show(io::IO, ::MIME"text/plain", x::Apparition)
    oflag, rflag = noptical(x) > 0, nradar(x) > 0
    lines, columns = displaysize(io)
    lines = (oflag && rflag) ? (lines - 2) ÷ 2 : lines - 1
    context = IOContext(io, :limit => true, :displaysize => (lines, columns))
    print(io, "Apparition including:\n")
    oflag && show(context, "text/plain", x.optical)
    (oflag && rflag) && print(io, "\n\n")
    rflag && show(context, "text/plain", x.radar)
    return nothing
end

"""
    apparitions(::AbstractObservationVector [, gap::Period])

Divide a sorted vector of astrometry into apparitions such that
the maximum time period between consecutive observations is less
than or equal to `gap` (default: `Day(30)`).

See also [`Apparition`](@ref).
"""
function apparitions(x::AbstractObservationVector, gap::Period = Day(30))
    @assert !isempty(x) && issorted(x) "Observation vector must be \
        non empty and sorted"
    # Initialize the vector of apparitions
    xview = view(x, 1:1)
    apps = Vector{typeof(Apparition(xview))}(undef, 0)
    # Initialize pointer and starting boundary
    i, start_i = 1, 1
    last_date = date(x[1])
    # Iterate astrometry vector
    @inbounds for i in eachindex(x)
        current_date = date(x[i])
        # Check gap trigger
        if current_date - last_date > gap
            xview = view(x, start_i:(i-1))
            push!(apps, Apparition(xview))
            # Reset boundary
            start_i = i
        end
        last_date = current_date
    end
    # Append the final remaining apparition
    xview = view(x, start_i:length(x))
    push!(apps, Apparition(xview))
    return apps
end

"""
    apparitions(::AbstractOpticalVector, ::AbstractRadarVector [, gap::Period])

Merge and divide sorted vectors of optical and radar astrometry into
apparitions such that the maximum time period between consecutive
observations is less than or equal to `gap` (default: `Day(30)`).
"""
function apparitions(optical::AbstractOpticalVector{T},
                     radar::AbstractRadarVector{T},
                     gap::Period = Day(30)) where {T <: Real}
    @assert issorted(optical) && issorted(radar) "Both optical and radar \
        observation vectors must be sorted"
    # Initialize the vector of apparitions
    Noptical, Nradar = length(optical), length(radar)
    oview, rview = view(optical, 1:0), view(radar, 1:0)
    apps = Vector{Apparition{T, typeof(oview), typeof(rview)}}(undef, 0)
    (iszero(Noptical) && iszero(Nradar)) && return apps
    # Initialize pointers and starting boundaries
    i, j = 1, 1
    start_i, start_j = 1, 1
    # Initialize the baseline date
    last_date = if Noptical > 0 && Nradar > 0
        min(date(optical[1]), date(radar[1]))
    elseif Noptical > 0
        date(optical[1])
    else
        date(radar[1])
    end
    # Iterate both astrometry vectors
    @inbounds while i ≤ Noptical || j ≤ Nradar
        # Select the earliest available observation
        if i ≤ Noptical && j ≤ Nradar
            if date(optical[i]) <= date(radar[j])
                current_date = date(optical[i])
                advance_optical = true
            else
                current_date = date(radar[j])
                advance_optical = false
            end
        elseif i ≤ Noptical
            current_date = date(optical[i])
            advance_optical = true
        else
            current_date = date(radar[j])
            advance_optical = false
        end
        # Check gap trigger
        if current_date - last_date > gap
            oview = start_i < i ? view(optical, start_i:(i-1)) : view(optical, 1:0)
            rview = start_j < j ? view(radar, start_j:(j-1)) : view(radar, 1:0)
            push!(apps, Apparition(oview, rview))
            # Reset boundaries
            start_i, start_j = i, j
        end
        last_date = current_date
        if advance_optical
            i += 1
        else
            j += 1
        end
    end
    # Append the final remaining apparition
    oview = start_i ≤ Noptical ? view(optical, start_i:Noptical) : view(optical, 1:0)
    rview = start_j ≤ Nradar ? view(rad, start_j:Nradar) : view(radar, 1:0)
    push!(apps, Apparition(oview, rview))

    return apps
end