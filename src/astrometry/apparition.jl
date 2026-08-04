"""
    Apparition{
        T <: Real,
        O <: AbstractOpticalAstrometry{T},
        V <: AbstractVector{O},
        I <: AbstractVector{Int},
        B
    }

A set of optical astrometry pertaining to the same object
such that the time between any two consecutive observations
is less than or equal to a given gap.

# Fields

- `optical::SubArray{O, 1, V, Tuple{I}, B}`: observations vector.
"""
struct Apparition{
        T <: Real,
        O <: AbstractOpticalAstrometry{T},
        V <: AbstractVector{O},
        I <: AbstractVector{Int},
        B
    }
    optical::SubArray{O, 1, V, Tuple{I}, B}
end

# Abbreviations
const AbstractApparitionVector{T} = AbstractVector{Apparition{T, O, V, I, B}} where {O, V, I, B}

# Apparition interface
indices(x::Apparition) = first(x.optical.indices)
indices(x::AbstractApparitionVector) = sort!(mapreduce(indices, vcat, x))
optical(x::Apparition) = collect(x.optical)
optical(x::AbstractApparitionVector) = sort!(mapreduce(optical, vcat, x))
numberofdays(x::Apparition) = numberofdays(x.optical)
numberofdays(x::AbstractApparitionVector) = sum(numberofdays, x)
noptical(x::Apparition) = length(x.optical)
noptical(x::AbstractApparitionVector) = sum(noptical, x)

# Print methods for Apparition
show(io::IO, x::Apparition) = print(io, "Apparition with ", noptical(x), " observation(s)")

function show(io::IO, ::MIME"text/plain", x::Apparition)
    print(io, noptical(x), "-observation Apparition{", eltype(x.optical), "}\n")
    Base.print_array(IOContext(io, :limit => true), x.optical)
    return nothing
end

"""
    apparitions(::AbstractOpticalVector; kwargs...)

Divide a vector of optical astrometry into apparitions.

See also [`Apparition`](@ref).

# Keyword arguments

- `gap::Period`: maximum time period between consecutive
    observations within an apparition (default: `Day(30)`).
"""
function apparitions(optical::AbstractOpticalVector, gap::Period = Day(30))
    sort!(optical)
    apps = [[1]]
    for i in 2:length(optical)
        if date(optical[i]) - date(optical[i-1]) > gap
            push!(apps, [i])
        else
            push!(apps[end], i)
        end
    end
    return [Apparition(view(optical, i)) for i in apps]
end