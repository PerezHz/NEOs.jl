# Return a zero of the same type as `a`
auxzero(a::AbstractSeries) = zero(a)

# Return a `TaylorN` with zero coefficients of the same type as `a.coeffs`
auxzero(a::TaylorN{Taylor1{T}}) where {T <: Number} = TaylorN(zero.(a.coeffs))

# TO DO: move this method to PlanetaryEphemeris in order to avoid type piracy
function kmsec2auday(x::SVector{6, T}) where {T <: Number}
    k = daysec / au
    return SVector{6, T}(x[1] / au, x[2] / au, x[3] / au,
                         x[4] * k, x[5] * k, x[6] * k)
end

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

function euclid3D(x::AbstractVector{TaylorN{T}}) where {T <: Real}
    z, w = zero(x[1]), zero(x[1])
    @inbounds for i in 1:3
        TS.zero!(w)
        for k in eachindex(x[i])
            TS.mul!(w, x[i], x[i], k)
            TS.add!(z, z, w, k)
        end
    end
    TS.zero!(w)
    for k in eachindex(z)
        TS.sqrt!(w, z, k)
    end
    return w
end

function euclid3D(x::AbstractVector{Taylor1{T}}) where {T <: Real}
    z, w = zero(x[1]), zero(x[1])
    @inbounds for i in 1:3
        TS.zero!(w)
        for k in eachindex(x[i])
            TS.mul!(w, x[i], x[i], k)
            TS.add!(z, z, w, k)
        end
    end
    TS.zero!(w)
    for k in eachindex(z)
        TS.sqrt!(w, z, k)
    end
    return w
end

dot3D(x::AbstractVector{T}, y::AbstractVector{T}) where {T <: Real} =
    x[1]*y[1] + x[2]*y[2] + x[3]*y[3]

function dot3D(x::Vector{TaylorN{T}}, y::Vector{U}) where {T <: Real, U <: Number}
    z, w = zero(x[1]), zero(x[1])
    @inbounds for i in 1:3
        TS.zero!(w)
        for k in eachindex(x[i])
            TS.mul!(w, x[i], y[i], k)
            TS.add!(z, z, w, k)
        end
    end
    return z
end

dot3D(x::Vector{T}, y::Vector{TaylorN{T}}) where {T <: Real} = dot3D(y, x)

function dot3D(x::Vector{Taylor1{T}}, y::Vector{U}) where {T <: Number, U <: Number}
    z, w = zero(x[1]), zero(x[1])
    @inbounds for i in 1:3
        TS.zero!(w)
        for k in eachindex(x[i])
            TS.mul!(w, x[i], y[i], k)
            TS.add!(z, z, w, k)
        end
    end
    return z
end

dot3D(x::Vector{T}, y::Vector{Taylor1{T}}) where {T <: Number} = dot3D(y, x)

# TO DO: move this methods to TaylorSeries' Static Arrays extension,
# in order to avoid type piracy
evaluate(y::SVector{N, Taylor1{T}}, x::Number) where {N, T <: Number} =
    [y[i](x) for i in eachindex(y)]

(y::SVector{N, Taylor1{T}})(x::Number) where {N, T <: Number} = evaluate(y, x)

evaluate(y::SVector{N, TaylorN{T}}, x::Vector{<:Number}) where {N, T <: Number} =
    [y[i](x) for i in eachindex(y)]

(y::SVector{N, TaylorN{T}})(x::Vector{T}) where {N, T <: Number} = evaluate(y, x)

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

evaleph(eph::TaylorInterpolant, t::Taylor1, q::TaylorN{Taylor1{T}}) where {T <: Real} =
    one(q) * eph(t)

evaleph(eph::NTuple{2, DensePropagation2{T, U}}, et::Number) where {T, U} =
    bwdfwdeph(et, eph[1], eph[2])

evaleph(eph::AstEph, et::Number) where {AstEph} = eph(et)

# In-place methods of evaleph
function evaleph!(y::Vector{Taylor1{TaylorN{T}}}, eph::TaylorInterpolant{T, T, 2},
                  t::Taylor1{T}) where {T <: Real}
    eph_t = tpeeval(eph, t)
    @inbounds for i in eachindex(eph_t)
        TS.zero!(y[i])
        for k in eachindex(y[i])
            y[i][k][0][1] = eph_t[i][k]
        end
    end
    return nothing
end

function evaleph!(y::Vector{Taylor1{T}}, eph::TaylorInterpolant{T, T, 2},
                  t::Taylor1{T}) where {T <: Real}
    eph_t = tpeeval(eph, t)
    @inbounds for i in eachindex(eph_t)
        for k in eachindex(y[i])
            TS.identity!(y[i], eph_t[i], k)
        end
    end
    return nothing
end

function evaleph!(y::Vector{Taylor1{Taylor1{T}}}, eph::TaylorInterpolant{T, T, 2},
                  t::Taylor1{T}) where {T <: Real}
    eph_t = tpeeval(eph, t)
    @inbounds for i in eachindex(eph_t)
        TS.zero!(y[i])
        for k in eachindex(y[i])
            y[i][k][0] = eph_t[i][k]
        end
    end
    return nothing
end

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
    y .= TS.evaluate(view(eph.x, ind, eachindex(y)), δt)
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
        y .= TS.evaluate(view(bwd.x, ind, eachindex(y)), δt)
    else
        y .= TS.evaluate(view(fwd.x, ind, eachindex(y)), δt)
    end
    # Convert state vector from [au, au/day] to [km, km/sec]
    auday2kmsec!(y)

    return nothing
end

seval(dsj2k::Taylor1{T}, zero_q::Taylor1{T}) where {T <: Real} =
    t2c_jpl_de430(dsj2k) .+ zero_q

seval(dsj2k::Taylor1{T}, zero_q::Taylor1{TaylorN{T}}) where {T <: Real} =
    t2c_jpl_de430(dsj2k) .+ zero_q

function seval(dsj2k::Taylor1{T}, zero_q::Taylor1{Taylor1{T}}) where {T <: Real}
    M = t2c_jpl_de430(dsj2k)
    one_q = one(zero_q.coeffs[1])
    _M_ = @. Taylor1(getfield(M, :coeffs) * one_q)
    return _M_
end