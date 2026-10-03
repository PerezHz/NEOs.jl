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
function auday2kmsec!(y::AbstractVector{T}) where {T <: Real}
    y[1:3] .*= au
    y[4:6] .*= au/daysec
    return nothing
end

function auday2kmsec!(y::AbstractVector{S}) where {S <: AbstractSeries}
    for i in eachindex(y)
        C = i <= 3 ? au : au/daysec
        for k in eachindex(y[i])
            multscalar!(y[i], C, k)
        end
    end
    return nothing
end

# In-place multiplication of the k-th order coefficient of `a` by the scalar `C`
function multscalar!(a::TaylorN{T}, C::Real, k::Int) where {T <: Real}
    a.coeffs[k+1].coeffs .*= C
    return nothing
end

function multscalar!(a::Taylor1{T}, C::Real, k::Int) where {T <: Real}
    a.coeffs[k+1] *= C
    return nothing
end

function multscalar!(a::Taylor1{TaylorN{T}}, C::Real, k::Int) where {T <: Real}
    for l in eachindex(a.coeffs[k+1])
        multscalar!(a.coeffs[k+1], C, l)
    end
    return nothing
end

# In-place copy of the k-th order coefficient of a `Taylor1{T}` into a `Taylor1{U}`;
# if `U` is a `TaylorN{T}`, the copied coefficient is a constant `TaylorN`
function taylorembed!(c::Taylor1{T}, a::Taylor1{T}, k::Int) where {T <: Real}
    c.coeffs[k+1] = a.coeffs[k+1]
    return nothing
end

function taylorembed!(c::Taylor1{TaylorN{T}}, a::Taylor1{T}, k::Int) where {T <: Real}
    TS.zero!(c.coeffs[k+1])
    c.coeffs[k+1].coeffs[1].coeffs[1] = a.coeffs[k+1]
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

# Evaluate `y` at time `t` using `OhMyThreads.tmap`.
function tpeeval(y::TaylorSolution{T, U, 2}, t::TT) where {T, U,
                 TT <: TaylorSolutionCallingArgs{T, U}}
    # Get index of y.p that interpolates at time t
    ind::Int, δt::TT = timeindex(y, t)
    # Evaluate y.p[ind] at δt
    return tmap(x -> x(δt), TT, view(y.p, ind, :))
end

"""
    evaleph(eph::TaylorSolution, t::Taylor1, q)

Evaluate `eph` at time `t` with type given by `q`.
"""
evaleph(eph::TaylorSolution, t::Taylor1, q::Taylor1{U}) where {U} =
    map(x -> Taylor1( x.coeffs * one(q[0]) ), tpeeval(eph, t))

# evaleph(eph::TaylorSolution, t::Taylor1, q::TaylorN{Taylor1{T}}) where {T <: Real} =
#    one(q) * eph(t)

evaleph(eph::NTuple{2, DensePropagation2{T, U}}, et::Number) where {T, U} =
    bwdfwdeph(et, eph[1], eph[2])

evaleph(eph::AstEph, et::Number) where {AstEph} = eph(et)

# In-place composition `c = a(x)` of a `Taylor1` polynomial `a` with an argument `x`,
# i.e. the evaluation of `a` at `x`; `aux` is an auxiliary variable, which must be
# different from `c` and `x`. All methods return `c`, except for real (immutable)
# `c`, in which case the result is returned instead; thus, the output must always be
# assigned, e.g. `y[i] = taylorcompose!(y[i], a, x, aux)`
taylorcompose!(::T, a::Taylor1{T}, x::T, ::T) where {T <: Real} = evaluate(a, x)

function taylorcompose!(c::TaylorN{T}, a::Taylor1{T}, x::TaylorN{T},
                        aux::TaylorN{T}) where {T <: Real}
    TS.zero!(c)
    TS._horner!(c, a, x, aux)
    return c
end

function taylorcompose!(c::TaylorN{T}, a::Taylor1{TaylorN{T}}, x::Number,
                        aux::TaylorN{T}) where {T <: Real}
    TS.zero!(c)
    @inbounds for k in reverse(eachindex(a))
        TS.zero!(aux)
        for ord in eachindex(c)
            TS.mul!(aux, c, x, ord)
        end
        for ord in eachindex(c)
            TS.add!(c, aux, a[k], ord)
        end
    end
    return c
end

function taylorcompose!(c::Taylor1{T}, a::Taylor1{Taylor1{T}}, x::Number,
                        aux::Taylor1{T}) where {T <: Real}
    TS.zero!(c)
    @inbounds for k in reverse(eachindex(a))
        TS.zero!(aux)
        for ord in eachindex(c)
            TS.mul!(aux, c, x, ord)
        end
        for ord in eachindex(c)
            TS.add!(c, aux, a[k], ord)
        end
    end
    return c
end

function taylorcompose!(c::Taylor1{U}, a::Taylor1{U}, x::Taylor1{U},
                        aux::Taylor1{U}) where {U <: Number}
    TS.zero!(c)
    TS._horner!(c, a, x, aux)
    return c
end

function taylorcompose!(c::Taylor1{U}, a::Taylor1{T}, x::Taylor1{U},
                        aux::Taylor1{U}) where {T <: Real, U <: AbstractSeries}
    TS.zero!(c)
    TS._horner!(c, a, x, aux)
    return c
end

# Specialized method that avoids the allocations of the generic Horner kernel
# when the coefficients of `a` are `TaylorN`s
function taylorcompose!(c::Taylor1{TaylorN{T}}, a::Taylor1{TaylorN{T}},
                        x::Taylor1{TaylorN{T}}, aux::Taylor1{TaylorN{T}}) where {T <: Real}
    TS.zero!(c)
    @inbounds for k in reverse(eachindex(a))
        # c <- c * x
        for ord in eachindex(c)
            TS.mul!(aux, c, x, ord)
        end
        for ord in eachindex(c)
            TS.identity!(c, aux, ord)
        end
        # c <- c + a[k]
        for ordQ in eachindex(c[0])
            TS.add!(c[0], c[0], a[k], ordQ)
        end
    end
    return c
end

# In-place evaluation of an ephemeris at time `et` [TDB seconds since J2000]. Only the
# first `length(y)` components are evaluated, and the state vector is converted from
# [au, au/day] to [km, km/sec]. `aux` is an auxiliary variable of the same type as the
# elements of `y`. For the asteroid ephemeris, pass the backward and forward integrations
function evaleph!(y::AbstractVector{U}, et::Number, eph::DensePropagation2,
                  aux::Number = zero(first(y))) where {U <: Number}
    # Convert time to TDB days since J2000
    t = et / daysec
    # Get index of eph.p that interpolates at time t
    ind::Int, δt = timeindex(eph, t)
    # Evaluate eph.p[ind] at δt
    for i in eachindex(y)
        y[i] = taylorcompose!(y[i], eph.p[ind, i], δt, aux)
    end
    # Convert state vector from [au, au/day] to [km, km/sec]
    auday2kmsec!(y)

    return nothing
end

function evaleph!(y::AbstractVector{U}, et::Number, bwd::DensePropagation2,
                  fwd::DensePropagation2, aux::Number = zero(first(y))) where
                  {U <: Number}
    # Convert time to TDB days since J2000
    t = et / daysec
    # Backward or forward integration
    eph = t <= firsttime(bwd) ? bwd : fwd
    # Get index of eph.p that interpolates at time t
    ind::Int, δt = timeindex(eph, t)
    # Evaluate eph.p[ind] at δt
    for i in eachindex(y)
        y[i] = taylorcompose!(y[i], eph.p[ind, i], δt, aux)
    end
    # Convert state vector from [au, au/day] to [km, km/sec]
    auday2kmsec!(y)

    return nothing
end
