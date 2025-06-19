"""
    AbstractAstrometry

Supertye for the minor bodies astrometry interface.
"""
abstract type AbstractAstrometry end

"""
    AbstractAstrometryCatalogue <: AbstractAstrometry

Supertype for the reference star catalogues interface.
"""
abstract type AbstractAstrometryCatalogue <: AbstractAstrometry end

"""
    AbstractAstrometryObservatory{T <: Real} <: AbstractAstrometry

Supertye for the observatories interface.
"""
abstract type AbstractAstrometryObservatory{T <: Real} <: AbstractAstrometry end

"""
    AbstractAstrometryTime <: AbstractAstrometry

Supertye for the observing time interface.
"""
abstract type AbstractAstrometryTime <: AbstractAstrometry end

"""
    AbstractAstrometrySource <: AbstractAstrometry

Supertye for the astrometry sources interface.
"""
abstract type AbstractAstrometrySource <: AbstractAstrometry end

"""
    AbstractAstrometryObservation{T <: Real} <: AbstractAstrometry

Supertye for the astrometric observations interface.

Every astrometric observation `x` has a:
- `date(x)`: date of observation.
- `measure(x)`: measurement.
- `observatory(x)`: observing station.
- `rms(x)`: a-priori formal RMS.
- `debias(x)`: debiasing factor.
"""
abstract type AbstractAstrometryObservation{T <: Real} <: AbstractAstrometry end

const AbstractObservationVector{T} = AbstractVector{<:AbstractAstrometryObservation{T}} where {T}

# Order in AbstractAstrometryObservation is given by date
isless(a::AbstractAstrometryObservation, b::AbstractAstrometryObservation) = date(a) < date(b)

dtutc2et(x::AbstractAstrometryObservation) = dtutc2et(date(x))

function numberofdays(radec::AbstractObservationVector{T}) where {T <: Real}
    t0, tf = extrema(date, radec)
    return (tf - t0).value / 86_400_000
end

# Methods to convert an AbstractObservationVector to a DataFrame
istable(::Type{<:AbstractObservationVector{<:Real}}) = true
rowaccess(::Type{<:AbstractObservationVector{<:Real}}) = true
rows(x::AbstractObservationVector{<:Real}) = x
schema(::AbstractVector{O}) where {O <: AbstractAstrometryObservation{<:Real}} =
    Schema(fieldnames(O), Tuple{fieldtypes(O)...})

"""
    AbstractOpticalAstrometry{T} <: AbstractAstrometryObservation{T}

Supertye for the optical astrometry interface.

Every optical observation `x` has a:
- `date(x)`: date of observation.
- `ra(x)`: right ascension [rad].
- `dec(x)`: declination [rad].
- `mag(x)`: apparent magnitude.
- `observatory(x)`: observing station.
- `catalogue(x)`: reference star catalogue.
- `rms(x)`: a-priori formal RMS [arcsec].
- `debias(x)`: debiasing factor [arcsec].
"""
abstract type AbstractOpticalAstrometry{T} <: AbstractAstrometryObservation{T} end

const AbstractOpticalVector{T} = AbstractVector{<:AbstractOpticalAstrometry{T}} where {T}

measure(x::AbstractOpticalAstrometry) = (ra(x), dec(x))
ra(x::AbstractOpticalAstrometry) = x.ra
dec(x::AbstractOpticalAstrometry) = x.dec
mag(x::AbstractOpticalAstrometry) = x.mag

"""
    AbstractRadarAstrometry{T} <: AbstractAstrometryObservation{T}

Supertye for the radar astrometry interface.

Every radar observation `x` has a:
- `date(x)`: date of observation.
- `measure(x)`: time-delay [us] or Doppler shift [Hz].
- `frequency(x)`: transmitter frequency [MHz].
- `observatory(x)`: observing station.
- `rms(x)`: a-priori formal RMS [same units as `measure(x)`].
- `debias(x)`: debiasing factor [same units as `measure(x)`].
"""
abstract type AbstractRadarAstrometry{T} <: AbstractAstrometryObservation{T} end

const AbstractRadarVector{T} = AbstractVector{<:AbstractRadarAstrometry{T}} where {T}

"""
    AbstractAstrometryErrorModel{T <: Real} <: AbstractAstrometry

Supertye for the optical astrometry error models interface.
"""
abstract type AbstractAstrometryErrorModel{T <: Real} <: AbstractAstrometry end

"""
    AbstractAstrometryResidual{T <: Real, U <: Number} <: AbstractAstrometry

Supertye for the astrometric observed minus computed residuals interface.
"""
abstract type AbstractAstrometryResidual{T <: Real, U <: Number} <: AbstractAstrometry end

"""
    AbstractOpticalResidual{T, U} <: AbstractAstrometryResidual{T, U}

Supertye for the optical observed minus computed residuals interface.
"""
abstract type AbstractOpticalResidual{T, U} <: AbstractAstrometryResidual{T, U} end

"""
    AbstractRadarResidual{T, U} <: AbstractAstrometryResidual{T, U}

Supertye for the radar observed minus computed residuals interface.
"""
abstract type AbstractRadarResidual{T, U} <: AbstractAstrometryResidual{T, U} end