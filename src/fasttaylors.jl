"""
    AbstractBuffer

Supertype for the buffers interface.
"""
abstract type AbstractBuffer end

show(io::IO, x::AbstractBuffer) = print(io, typeof(x))

# Return a zero of the same type as `a`
auxzero(a::AbstractSeries) = zero(a)

# Return a `TaylorN` with zero coefficients of the same type as `a.coeffs`
# auxzero(a::TaylorN{Taylor1{T}}) where {T <: Number} = TaylorN(zero.(a.coeffs))

# Extract the linear scaling factor from a TaylorN
scalingfactor(x::TaylorN{T}) where {T <: Real} = x[1][findfirst(x[1])]

# In-place methods of auday2kmsec
function auday2kmsec!(y::Vector{T}) where {T <: Real}
    y[1:3] .*= au
    y[4:6] .*= au/daysec
    return nothing
end

function auday2kmsec!(y::Vector{TaylorN{T}}) where {T <: Real}
    for i in eachindex(y)
        k = i <= 3 ? au : au/daysec
        for j in eachindex(y[i])
            y[i].coeffs[j+1].coeffs .*= k
        end
    end
    return nothing
end

function auday2kmsec!(y::Vector{Taylor1{T}}) where {T <: Real}
    for i in eachindex(y)
        k = i <= 3 ? au : au/daysec
        y[i].coeffs .*= k
    end
    return nothing
end

# Warning: functions euclid3D(x) and dot3D(x) assume length(x) >= 3
euclid3D(x::AbstractVector{T}) where {T <: Real} = sqrt(dot3D(x, x))

function euclid3D!(z::S, x::AbstractVector{S}, aux1::S,
                   aux2::S, ord::Int) where {S <: AbstractSeries}
    dot3D!(aux1, x, x, z, ord)
    TS.zero!(z, ord)
    TS.zero!(aux2, ord)
    TS.sqrt!(z, aux1, aux2, ord)
    return nothing
end

function euclid3D(x::AbstractVector{S}) where {S <: AbstractSeries}
    z, aux1, aux2 = zero(x[1]), zero(x[1]), zero(x[1])
    for ord in eachindex(z)
        euclid3D!(z, x, aux1, aux2, ord)
    end
    return z
end

dot3D(x::AbstractVector{T}, y::AbstractVector{T}) where {T <: Real} =
    x[1]*y[1] + x[2]*y[2] + x[3]*y[3]

function dot3D!(z::S, x::AbstractVector{S}, y::AbstractVector{U}, aux::S,
                ord::Int) where {S <: AbstractSeries, U <: Number}
    TS.zero!(z, ord)
    @inbounds for i in 1:3
        TS.zero!(aux, ord)
        TS.mul!(aux, x[i], y[i], ord)
        TS.add!(z, z, aux, ord)
    end
    return nothing
end

function dot3D!(z::S, x::AbstractVector{T}, y::AbstractVector{S}, aux::S,
                ord::Int) where {T <: Real, S <: AbstractSeries}
    return dot3D!(z, y, x, aux, ord)
end

function dot3D(
        x::AbstractVector{S}, y::AbstractVector{U}
    ) where {S <: AbstractSeries, U <: Number}
    z, aux = zero(x[1]), zero(x[1])
    for ord in eachindex(z)
        dot3D!(z, x, y, aux, ord)
    end
    return z
end

function dot3D(
        x::AbstractVector{T}, y::AbstractVector{S}
    ) where {T <: Real, S <: AbstractSeries}
    return dot3D(y, x)
end

# Evaluate `y` at time `t` using `OhMyThreads.tmap`
# This function is a multithreaded version of
# (::TaylorInterpolant{T, U, 2})(t::TT) where {T, U, TT<:TaylorInterpCallingArgs{T,U}}
# at
# https://github.com/PerezHz/PlanetaryEphemeris.jl/src/interpolation/TaylorInterpolant.jl
# TO DO: move this to PlanetaryEphemeris.jl
function tpeeval(y::TaylorInterpolant{T, U, 2}, t::TT) where {T, U,
                 TT <: TaylorInterpCallingArgs{T, U}}
    # Get index of y.x that interpolates at time t
    ind::Int, δt::TT = getinterpindex(y, t)
    # Evaluate y.x[ind] at δt
    return tmap(x -> x(δt), TT, view(y.x, ind, :))
end

"""
    evaleph(eph::TaylorInterpolant, t::Taylor1, q)

Evaluate `eph` at time `t` with type given by `q`.
"""
evaleph(eph::TaylorInterpolant, t::Taylor1, q::Taylor1{U}) where {U} =
    map(x -> Taylor1( x.coeffs * one(q[0]) ), tpeeval(eph, t))

# evaleph(eph::TaylorInterpolant, t::Taylor1, q::TaylorN{Taylor1{T}}) where {T <: Real} =
#    one(q) * eph(t)

evaleph(eph::NTuple{2, DensePropagation2{T, U}}, et::Number) where {T, U} =
    bwdfwdeph(et, eph[1], eph[2])

evaleph(eph::AstEph, et::Number) where {AstEph} = eph(et)

# In-place methods of evaleph
function evaleph!(c::TaylorN{T}, a::Taylor1{TaylorN{T}}, dx::Number) where {T <: Real}
    TS.zero!(c)
    d = zero(c)
    @inbounds for k in reverse(eachindex(a))
        TS.zero!(d)
        for ord in eachindex(c)
            TS.mul!(d, c, dx, ord)
        end
        for ord in eachindex(c)
            TS.add!(c, d, a[k], ord)
        end
    end
    return nothing
end

function evaleph!(c::Taylor1{T}, a::Taylor1{Taylor1{T}}, dx::Number) where {T <: Real}
    TS.zero!(c)
    d = zero(c)
    @inbounds for k in reverse(eachindex(a))
        TS.zero!(d)
        for ord in eachindex(c)
            TS.mul!(d, c, dx, ord)
        end
        for ord in eachindex(c)
            TS.add!(c, d, a[k], ord)
        end
    end
    return nothing
end

function evaleph!(y::Vector{U}, et::U,
                  eph::DensePropagation2{T, T}) where {T <: Real, U <: Number}
    # Convert time from seconds to days (TDB) since J2000
    t = et / daysec
    # Get index of bwd/fwd that interpolates at time t
    ind::Int, δt::U = getinterpindex(eph, t)
    # Evaluate bwd/fwd.x[ind] at δt
    y .= evaluate(view(eph.x, ind, eachindex(y)), δt)
    # Convert state vector from [au, au/day] to [km, km/sec]
    auday2kmsec!(y)

    return nothing
end

function evaleph!(y::Vector{U}, et::U, bwd::DensePropagation2{T, U},
                 fwd::DensePropagation2{T, U}) where {T <: Real, U <: Number}
    # Convert time to TDB days since J2000
    t = et / daysec
    # Get index of bwd/fwd that interpolates at time t
    ind::Int, δt::U = t <= bwd.t0 ? getinterpindex(bwd, t) : getinterpindex(fwd, t)
    # Evaluate bwd/fwd.x[ind] at δt
    if t <= bwd.t0
        for i in eachindex(y)
            evaleph!(y[i], bwd.x[ind, i], δt)
        end
    else
        for i in eachindex(y)
            evaleph!(y[i], fwd.x[ind, i], δt)
        end
    end
    # Convert state vector from [au, au/day] to [km, km/sec]
    auday2kmsec!(y)

    return nothing
end

function evaleph!(y::Vector{U}, et::T, bwd::DensePropagation2{T, U},
                 fwd::DensePropagation2{T, U}) where {T <: Real, U <: Number}
    # Convert time to TDB days since J2000
    t = et / daysec
    # Get index of bwd/fwd that interpolates at time t
    ind::Int, δt::T = t <= bwd.t0 ? getinterpindex(bwd, t) : getinterpindex(fwd, t)
    # Evaluate bwd/fwd.x[ind] at δt
    if t <= bwd.t0
        for i in eachindex(y)
            evaleph!(y[i], bwd.x[ind, i], δt)
        end
    else
        for i in eachindex(y)
            evaleph!(y[i], fwd.x[ind, i], δt)
        end
    end
    # Convert state vector from [au, au/day] to [km, km/sec]
    auday2kmsec!(y)

    return nothing
end

function evaleph!(y::Vector{T}, et::T, bwd::DensePropagation2{T, T},
                 fwd::DensePropagation2{T, T}) where {T <: Real}
    # Convert time to TDB days since J2000
    t = et / daysec
    # Get index of bwd/fwd that interpolates at time t
    ind::Int, δt::T = t <= bwd.t0 ? getinterpindex(bwd, t) : getinterpindex(fwd, t)
    # Evaluate bwd/fwd.x[ind] at δt
    if t <= bwd.t0
        y .= evaluate(view(bwd.x, ind, eachindex(y)), δt)
    else
        y .= evaluate(view(fwd.x, ind, eachindex(y)), δt)
    end
    # Convert state vector from [au, au/day] to [km, km/sec]
    auday2kmsec!(y)

    return nothing
end