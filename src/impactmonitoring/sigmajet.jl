"""
    SigmaJet{N, T <: Real} <: Real

Truncated polynomial of degree `N - 1` in the LOV index offset, with coefficients
of type `T`.

`Taylor1{SigmaJet{N, T}}` replaces `Taylor1{Taylor1{T}}` in the propagation of
virtual asteroids. Since `SigmaJet <: Real`, it uses the in-place scalar kernels
of `TaylorSeries`, and its arithmetic reproduces the floating point operations
of the nested `Taylor1` kernels (in the same order), but without heap allocations.

# Fields

- `c::NTuple{N, T}`: coefficients; `c[i]` multiplies `dσ^(i-1)`.
"""
struct SigmaJet{N, T <: Real} <: Real
    c::NTuple{N, T}
    SigmaJet{N, T}(c::NTuple{N, T}) where {N, T <: Real} = new{N, T}(c)
end

# Outer constructors
SigmaJet(c::Tuple{T, Vararg{T}}) where {T <: Real} = SigmaJet{length(c), T}(c)
SigmaJet{N, T}(c::Tuple) where {N, T <: Real} = SigmaJet{N, T}(NTuple{N, T}(c))
SigmaJet{N, T}(x::SigmaJet{N, T}) where {N, T <: Real} = x
SigmaJet{N, T}(x::SigmaJet{N}) where {N, T <: Real} = SigmaJet{N, T}(map(T, x.c))
SigmaJet{N, T}(x::Real) where {N, T <: Real} =
    SigmaJet{N, T}(ntuple(i -> i == 1 ? T(x) : zero(T), Val(N)))

function SigmaJet{N, T}(x::Taylor1) where {N, T <: Real}
    @assert TaylorSeries.order(x) == N - 1 "Taylor1 order must be $(N - 1)"
    return SigmaJet{N, T}(ntuple(i -> T(x.coeffs[i]), Val(N)))
end

SigmaJet(x::Taylor1{T}) where {T <: Real} = SigmaJet{TaylorSeries.order(x) + 1, T}(x)

# Print method for SigmaJet
Base.show(io::IO, x::SigmaJet{N, T}) where {N, T <: Real} =
    print(io, "SigmaJet{", N, ", ", T, "}", x.c)

# Basic interface
Base.zero(::Type{SigmaJet{N, T}}) where {N, T <: Real} =
    SigmaJet{N, T}(ntuple(_ -> zero(T), Val(N)))
Base.zero(x::SigmaJet) = zero(typeof(x))
Base.one(::Type{SigmaJet{N, T}}) where {N, T <: Real} = SigmaJet{N, T}(one(T))
Base.one(x::SigmaJet) = one(typeof(x))
Base.iszero(x::SigmaJet) = all(iszero, x.c)
Base.isnan(x::SigmaJet) = any(isnan, x.c)
Base.isinf(x::SigmaJet) = any(isinf, x.c)
Base.isfinite(x::SigmaJet) = all(isfinite, x.c)
Base.float(::Type{SigmaJet{N, T}}) where {N, T <: Real} = SigmaJet{N, float(T)}
Base.float(x::SigmaJet{N, T}) where {N, T <: Real} = SigmaJet{N, float(T)}(x)
Base.isequal(x::SigmaJet, y::SigmaJet) = isequal(x.c, y.c)
Base.hash(x::SigmaJet, h::UInt) = hash(x.c, hash(:SigmaJet, h))

Base.promote_rule(::Type{SigmaJet{N, T}}, ::Type{S}) where {N, T <: Real, S <: Real} =
    SigmaJet{N, promote_type(T, S)}
Base.promote_rule(::Type{SigmaJet{N, T}}, ::Type{SigmaJet{N, S}}) where
    {N, T <: Real, S <: Real} = SigmaJet{N, promote_type(T, S)}

# As for the coefficients of a nested Taylor1, the constant term is a scalar
TaylorSeries.constant_term(x::SigmaJet) = x.c[1]

# Infinity norm as in `TaylorIntegration._norminf` for a Taylor1
function LinearAlgebra.norm(x::SigmaJet{N, T}, p::Real = 2) where {N, T <: Real}
    p == Inf || return norm(x.c, p)
    h = zero(T)
    for c in x.c
        h = max(h, abs(c))
    end
    return h
end

# Index of the first (last) non-zero coefficient, starting from 0; -1 if there is none
function _sj_findfirst(a::Tuple)
    for i in eachindex(a)
        iszero(a[i]) || return i - 1
    end
    return -1
end

function _sj_findlast(a::Tuple)
    for i in reverse(eachindex(a))
        iszero(a[i]) || return i - 1
    end
    return -1
end

# Comparisons, with the semantics of `Taylor1`
Base.:(==)(x::SigmaJet, y::SigmaJet) = x.c == y.c
Base.:(==)(x::SigmaJet, b::Real) = (x.c[1] == b) & all(iszero, Base.tail(x.c))
Base.:(==)(b::Real, x::SigmaJet) = x == b
Base.:(==)(x::SigmaJet, b::AbstractFloat) = (x.c[1] == b) & all(iszero, Base.tail(x.c))
Base.:(==)(b::AbstractFloat, x::SigmaJet) = x == b

function _sj_isless(x::SigmaJet, b::Real)
    x0 = x.c[1]
    x0 != b && return isless(x0, b)
    for i in eachindex(x.c)
        d = i == 1 ? x0 - b : x.c[i] - zero(b)
        iszero(d) || return isless(x.c[i], zero(b))
    end
    return isless(x0, b)
end

function _sj_isless(b::Real, x::SigmaJet)
    x0 = x.c[1]
    x0 != b && return isless(b, x0)
    for i in eachindex(x.c)
        d = i == 1 ? b - x0 : zero(b) - x.c[i]
        iszero(d) || return isless(zero(b), x.c[i])
    end
    return isless(b, x0)
end

for S in (:Real, :AbstractFloat)
    @eval begin
        Base.isless(x::SigmaJet, b::$S) = _sj_isless(x, b)
        Base.isless(b::$S, x::SigmaJet) = _sj_isless(b, x)
        Base.:<(x::SigmaJet, b::$S) = _sj_isless(x, b)
        Base.:<(b::$S, x::SigmaJet) = _sj_isless(b, x)
    end
end

Base.isless(x::SigmaJet{N, T}, y::SigmaJet{N, T}) where {N, T <: Real} = isless(x - y, zero(T))
Base.isless(x::SigmaJet{N, T}, y::SigmaJet{N, S}) where {N, T <: Real, S <: Real} =
    isless(promote(x, y)...)
Base.:<(x::SigmaJet{N}, y::SigmaJet{N}) where {N} = isless(x, y)

# Arithmetic. Each method reproduces the floating point operations, in the same
# order, of the analogous in-place `TaylorSeries` kernel for `Taylor1{T}`
Base.:+(x::SigmaJet) = x
Base.:-(x::SigmaJet) = SigmaJet(map(-, x.c))
Base.:+(x::SigmaJet{N, T}, y::SigmaJet{N, T}) where {N, T <: Real} = SigmaJet(map(+, x.c, y.c))
Base.:-(x::SigmaJet{N, T}, y::SigmaJet{N, T}) where {N, T <: Real} = SigmaJet(map(-, x.c, y.c))
Base.:*(x::SigmaJet{N, T}, y::SigmaJet{N, T}) where {N, T <: Real} = _sj_mul(x, y)

# Different coefficient types; same types dispatch to the methods above
for f in (:+, :-, :*, :/)
    @eval Base.$f(x::SigmaJet{N, T}, y::SigmaJet{N, S}) where {N, T <: Real, S <: Real} =
        $f(promote(x, y)...)
end

# A scalar only modifies the constant term; the zero added to the other
# coefficients mimics `add!` and `subst!`
Base.:+(x::SigmaJet{N}, b::Real) where {N} =
    SigmaJet(ntuple(i -> x.c[i] + (i == 1 ? b : zero(b)), Val(N)))
Base.:+(b::Real, x::SigmaJet{N}) where {N} =
    SigmaJet(ntuple(i -> (i == 1 ? b : zero(b)) + x.c[i], Val(N)))
Base.:-(x::SigmaJet{N}, b::Real) where {N} =
    SigmaJet(ntuple(i -> x.c[i] - (i == 1 ? b : zero(b)), Val(N)))
Base.:-(b::Real, x::SigmaJet{N}) where {N} =
    SigmaJet(ntuple(i -> (i == 1 ? b : zero(b)) - x.c[i], Val(N)))
Base.:*(x::SigmaJet, b::Real) = SigmaJet(map(c -> c * b, x.c))
Base.:*(b::Real, x::SigmaJet) = SigmaJet(map(c -> c * b, x.c))
Base.:/(x::SigmaJet, b::Real) = SigmaJet(map(c -> c / b, x.c))

# Product: each coefficient is accumulated from zero (`mul!`)
@generated function _sj_mul(x::SigmaJet{N, T}, y::SigmaJet{N, T}) where {N, T <: Real}
    exs = Expr[]
    for j in 1:N
        ex = :(zero($T))
        for l in 1:j
            ex = :($ex + x.c[$l] * y.c[$(j - l + 1)])
        end
        push!(exs, ex)
    end
    return :(@inbounds SigmaJet{$N, $T}(($(exs...),)))
end

# Square (`sqr!`)
@generated function _sj_sqr(x::SigmaJet{N, T}) where {N, T <: Real}
    exs = Expr[:(x.c[1]^2)]
    for k in 1:N-1
        kodd = k % 2
        kend = (k - 2 + kodd) >> 1
        ex = :(zero($T))
        for i in 0:kend
            ex = :($ex + x.c[$(i + 1)] * x.c[$(k - i + 1)])
        end
        ex = :(2 * $ex)
        if kodd == 0
            ex = :($ex + x.c[$((k >> 1) + 1)]^2)
        end
        push!(exs, ex)
    end
    return :(@inbounds SigmaJet{$N, $T}(($(exs...),)))
end

# Division (`div!`)
function Base.:/(x::SigmaJet{N, T}, y::SigmaJet{N, T}) where {N, T <: Real}
    iszero(x) && !iszero(y) && return zero(x)
    iszero(y.c[1]) && return _sj_divfact(x, y)
    return _sj_div(x, y)
end

# Division with a non-zero constant term in the denominator
@generated function _sj_div(x::SigmaJet{N, T}, y::SigmaJet{N, T}) where {N, T <: Real}
    q = Symbol[]
    for j in 1:N
        push!(q, Symbol(:q, j))
    end
    body = Expr[:($(q[1]) = x.c[1] / y.c[1])]
    for j in 2:N
        acc = :($(q[1]) * y.c[$j])
        for i in 2:j-1
            acc = :($acc + $(q[i]) * y.c[$(j - i + 1)])
        end
        push!(body, :($(q[j]) = (x.c[$j] - $acc) / y.c[1]))
    end
    return quote
        @inbounds begin
            $(body...)
        end
        return SigmaJet{$N, $T}(($(q...),))
    end
end

# Division with factorization of the leading zeros
function _sj_divfact(x::SigmaJet{N, T}, y::SigmaJet{N, T}) where {N, T <: Real}
    a, b, c = x.c, y.c, zero(x).c
    n = N - 1
    anz, bnz = _sj_findfirst(a), _sj_findfirst(b)
    ordfact = min(anz ≥ 0 ? anz : n, bnz ≥ 0 ? bnz : n)
    @inbounds for k in 0:n
        if k == 0
            c = Base.setindex(c, a[ordfact+1] / b[ordfact+1], 1)
            continue
        end
        imin = max(0, k + ordfact - n)
        acc = c[imin+1] * b[k+ordfact-imin+1]
        for i in imin+1:k-1
            acc += c[i+1] * b[k+ordfact-i+1]
        end
        acc = k + ordfact ≤ n ? a[k+ordfact+1] - acc : -acc
        c = Base.setindex(c, acc / b[ordfact+1], k+1)
    end
    return SigmaJet{N, T}(c)
end

# Scalar over SigmaJet (`div!` with a scalar numerator)
function Base.:/(b::Real, y::SigmaJet)
    iszero(b) && !iszero(y) && return zero(_sj_rdiv(b, y))
    return _sj_rdiv(b, y)
end

@generated function _sj_rdiv(b::Real, y::SigmaJet{N, T}) where {N, T <: Real}
    q = Symbol[]
    for j in 1:N
        push!(q, Symbol(:q, j))
    end
    body = Expr[:($(q[1]) = b / y.c[1])]
    for j in 2:N
        acc = :($(q[1]) * y.c[$j])
        for i in 2:j-1
            acc = :($acc + $(q[i]) * y.c[$(j - i + 1)])
        end
        push!(body, :($(q[j]) = -$acc / y.c[1]))
    end
    return quote
        @inbounds begin
            $(body...)
        end
        return SigmaJet(($(q...),))
    end
end

Base.inv(x::SigmaJet{N, T}) where {N, T <: Real} = one(T) / x

# Square root (`sqrt!` with first non-zero coefficient `2k0`)
function _sj_sqrt(x::SigmaJet{N, T}, k0::Int) where {N, T <: Real}
    a, c = x.c, zero(x).c
    n = N - 1
    @inbounds for k in k0:n
        if k == k0
            c = Base.setindex(c, sqrt(a[2*k0+1]), k+1)
            continue
        end
        kodd = (k - k0) % 2
        kend = (k - k0 - 2 + kodd) >> 1
        imax = min(k0 + kend, n)
        imin = max(k0 + 1, k + k0 - n)
        acc = zero(T)
        if k + k0 ≤ n
            acc = a[k+k0+1]
        end
        if kodd == 0
            acc -= (c[kend+k0+2])^2
        end
        imin ≤ imax && (acc -= 2 * c[imin+1] * c[k+k0-imin+1])
        for i in imin+1:imax
            acc -= 2 * c[i+1] * c[k+k0-i+1]
        end
        c = Base.setindex(c, acc / (2*c[k0+1]), k+1)
    end
    return SigmaJet{N, T}(c)
end

# As `sqrt!` for a nested Taylor1 (first non-zero coefficient detected)
function Base.sqrt(x::SigmaJet)
    l0 = _sj_findfirst(x.c)
    l0 < 0 && return zero(x)
    isodd(l0) && throw(DomainError(x,
        """First non-vanishing SigmaJet coefficient must correspond
        to an **even power** in order to expand `sqrt` around 0."""))
    return _sj_sqrt(x, l0 >> 1)
end

# Powers, as `^` and `pow!` for a Taylor1
function Base.:^(x::SigmaJet, n::Integer)
    n == 0 && return one(x)
    n == 1 && return x
    n == 2 && return _sj_sqr(x)
    return x^float(n)
end

Base.:^(x::SigmaJet, r::Rational) = x^float(r)
Base.:^(x::SigmaJet, r::Real) = _sj_pow(x, r)
Base.:^(::SigmaJet, ::SigmaJet) = throw(ArgumentError("SigmaJet exponents are not supported"))

function _sj_pow(x::SigmaJet{N, T}, r::Real) where {N, T <: Real}
    (r == 0) && return one(x)
    (r == 1) && return x
    (r == 2) && return _sj_sqr(x)
    (r == 0.5) && return _sj_sqrt(x, 0)
    a, c = x.c, zero(x).c
    order_a = N - 1
    l0 = _sj_findfirst(a)
    l0 < 0 && return SigmaJet{N, T}(c)
    !isinteger(r*l0) && throw(DomainError(x,
        """The 0-th order SigmaJet coefficient must be non-zero
        to raise the SigmaJet to a non-integer exponent."""))
    lnull = trunc(Int, r*l0)
    lastnz = isinteger(r) && r > 0 ? _sj_findlast(a) : order_a
    @inbounds a_l0 = a[l0+1]
    @inbounds for k in 0:order_a
        kprime = k - lnull
        (kprime < 0 || lnull > order_a) && continue
        isinteger(r) && r > 0 && (k > r*lastnz) && continue
        if k == lnull
            c = Base.setindex(c, a_l0^float(r), k+1)
            continue
        end
        acc = zero(T)
        if l0 + kprime ≤ order_a
            acc = r * kprime * c[lnull+1] * a[l0+kprime+1]
        end
        ilo = max(1, l0 + kprime - order_a)
        ihi = min(k - lnull - 1, order_a - lnull)
        for i in ilo:ihi
            aaux = r*(kprime-i) - i
            acc += aaux * c[i+lnull+1] * a[l0+kprime-i+1]
        end
        c = Base.setindex(c, acc / (kprime * a_l0), k+1)
    end
    return SigmaJet{N, T}(c)
end

function Base.abs(x::SigmaJet)
    x0 = x.c[1]
    x0 > 0 && return x
    x0 < 0 && return -x
    throw(DomainError(x,
        """The 0th order SigmaJet coefficient must be non-zero
        (abs(x) is not differentiable at x=0)."""))
end

# The generic `pow!` for `Taylor1{T}` accumulates the recursion in a different
# order than `pow!` for `Taylor1{Taylor1{T}}`, which is reproduced here
function TaylorSeries.pow!(c::Taylor1{SigmaJet{N, T}}, a::Taylor1{SigmaJet{N, T}},
        aux::Taylor1{SigmaJet{N, T}}, r::S, k::Int) where {N, T <: Real, S <: Real}
    (r == 0) && return TaylorSeries.one!(c, a, k)
    (r == 1) && return TaylorSeries.identity!(c, a, k)
    (r == 2) && return TaylorSeries.sqr!(c, a, constant_term(aux), k)
    (r == 0.5) && return TaylorSeries.sqrt!(c, a, aux, k)
    c_coeffs, a_coeffs = c.coeffs, a.coeffs
    @inbounds c_coeffs[k+1] = zero(c_coeffs[k+1])
    l0 = findfirst(a)
    l0 < 0 && return nothing
    !isinteger(r*l0) && throw(DomainError(a,
        """The 0-th order Taylor1 coefficient must be non-zero
        to raise the Taylor1 polynomial to a non-integer exponent."""))
    lnull = trunc(Int, r*l0)
    kprime = k - lnull
    order_a = TaylorSeries.order(a)
    (kprime < 0 || lnull > order_a) && return nothing
    lastnz = isinteger(r) && r > 0 ? findlast(a) : order_a
    isinteger(r) && r > 0 && (k > r*lastnz) && return nothing
    @inbounds a_l0 = a_coeffs[l0+1]
    if k == lnull
        @inbounds c_coeffs[k+1] = a_l0^float(r)
        return nothing
    end
    ilo = max(0, l0 + kprime - order_a)
    ihi = min(k - lnull - 1, order_a - lnull)
    acc = zero(a_l0)
    @inbounds for i in ilo:ihi
        rr = r*(kprime-i) - i
        acc += (c_coeffs[i+lnull+1] * a_coeffs[l0+kprime-i+1]) * rr
    end
    @inbounds c_coeffs[k+1] = (acc / a_l0) / kprime
    return nothing
end

# In-place copy of the k-th coefficient of a `Taylor1{T}` as a constant SigmaJet
function taylorembed!(c::Taylor1{SigmaJet{N, T}}, a::Taylor1{T}, k::Int) where {N, T <: Real}
    @inbounds c.coeffs[k+1] = SigmaJet{N, T}(a.coeffs[k+1])
    return nothing
end

# Earth's orientation matrix (only the constant SigmaJet terms are non-zero)
function PE.t2c_jpl_de430!(M::Array{Taylor1{SigmaJet{N, T}}, 3}, ea::Int,
        dsj2k::Taylor1{T}, ::Taylor1{SigmaJet{N, T}},
        rotatBuf::RetAlloc{Taylor1{T}}) where {N, T <: Real}
    PE.c2t_jpl_de430!(dsj2k, rotatBuf)
    for j in axes(M, 2)
        for i in axes(M, 1)
            for ord in eachindex(M[i, j, ea])
                M[i, j, ea][ord] = SigmaJet{N, T}(rotatBuf.v2[12][j, i][ord])
            end
        end
    end
    return M
end

# Default initial condition for the propagation of a virtual asteroid
propagationic(q0::Vector{Taylor1{T}}) where {T <: Real} = map(SigmaJet, q0)
propagationic(q0::Vector{<:Number}) = q0

# Conversion of the initial condition of a virtual asteroid into the
# number type `V` of the propagation
convert_ic(::Type{U}, q0::Vector{U}) where {U <: Number} = q0
convert_ic(::Type{SigmaJet{N, T}}, q0::Vector{Taylor1{T}}) where {N, T <: Real} =
    map(SigmaJet{N, T}, q0)

# Conversion of a propagation number into the number type `U` of a virtual asteroid
vanumber(::Type{U}, x::U) where {U <: Number} = x
vanumber(::Type{Taylor1{T}}, x::SigmaJet{N, T}) where {N, T <: Real} =
    Taylor1(collect(x.c), N - 1)
