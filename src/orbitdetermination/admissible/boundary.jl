
# W function of an [`AdmissibleRegion`](@ref).
# See equation (8.9) of https://doi.org/10.1017/CBO9781139175371
arW(A::AdmissibleRegion, ρ::Number) = arW(A.coeffs, ρ)
arW(coeffs::AbstractVector, ρ::Number) = coeffs[3] * ρ^2 + coeffs[4] * ρ + coeffs[5]

ardW(A::AdmissibleRegion, ρ::Number) = ardW(A.coeffs, ρ)
ardW(coeffs::AbstractVector, ρ::Number) = 2 * coeffs[3] * ρ + coeffs[4]

ard2W(A::AdmissibleRegion, ρ::Number) = ard2W(A.coeffs, ρ)
ard2W(coeffs::AbstractVector, ρ::Number) = 2 * coeffs[3]

# S function of an [`AdmissibleRegion`](@ref).
# See equation (8.9) of https://doi.org/10.1017/CBO9781139175371
arS(A::AdmissibleRegion, ρ::Number) = arS(A.coeffs, ρ)
arS(coeffs::AbstractVector, ρ::Number) = ρ^2 + coeffs[6] * ρ + coeffs[1]

ardS(A::AdmissibleRegion, ρ::Number) = ardS(A.coeffs, ρ)
ardS(coeffs::AbstractVector, ρ::Number) = 2 * ρ + coeffs[6]

ard2S(A::AdmissibleRegion, ρ::Number) = ard2S(A.coeffs, ρ)
ard2S(coeffs::AbstractVector, ρ::Number) = 2

# Auxiliary function to compute the root of G(::AdmissibleRegion)
G⁻¹0(A::AdmissibleRegion) = G⁻¹0(A.coeffs)
G⁻¹0(coeffs::AbstractVector) = cbrt(2 * k_gauss^2 * μ_ES / coeffs[3])

# G function of an [`AdmissibleRegion`](@ref).
# See equation (8.13) of https://doi.org/10.1017/CBO9781139175371
arG(A::AdmissibleRegion, ρ::Number) = arG(A.coeffs, ρ)

function arG(coeffs::AbstractVector, ρ::Number)
    if ρ == G⁻¹0(coeffs)
        return zero(coeffs[3] * ρ)
    else
        return 2 * k_gauss^2 * μ_ES / ρ - coeffs[3] * ρ^2
    end
end

# Return the coefficients of `A`'s energy as a quadratic function
# of the topocentric range-rate evaluated at range `ρ`. `boundary`
# chooses between the `:outer` (default) or `:inner` boundary.
# See between equations (8.9)-(8.10) and (8.12)-(8.13)
# of https://doi.org/10.1017/CBO9781139175371
function arenergycoeffs(A::AdmissibleRegion, ρ::Number, boundary::Symbol = :outer)
    if boundary == :outer
        return _arhelenergycoeffs(A.coeffs, A.a_max, ρ)
    elseif boundary == :inner
        return _argeoenergycoeffs(A.coeffs, ρ)
    else
        throw(ArgumentError("Argument `boundary` must be either `:outer` or `:inner`"))
    end
end

function _arhelenergycoeffs(coeffs::Vector{T}, a_max::T, ρ::Number) where {T <: Real}
    a = one(T)
    b = coeffs[2]
    c = arW(coeffs, ρ) + k_gauss^2 * (1/a_max - 2/sqrt(arS(coeffs, ρ)))
    return a, b, c
end

function _arhelenergycoeffs_derivatives(coeffs::Vector{T}, a_max::T,
                                        ρ::Number) where {T <: Real}
    a = one(T)
    b = coeffs[2]
    W, dW, d2W = arW(coeffs, ρ), ardW(coeffs, ρ), ard2W(coeffs, ρ)
    S, dS, d2S = arS(coeffs, ρ), ardS(coeffs, ρ), ard2S(coeffs, ρ)
    sqrtS = sqrt(S)
    c = W + k_gauss^2 * (1/a_max - 2/sqrtS)
    dc = dW + k_gauss^2 * dS / sqrtS^3
    d2c = d2W + k_gauss^2 * (d2S / sqrtS^3 - 3 * dS^2 / (2 * sqrtS^5))
    return a, b, c, dc, d2c
end

function _argeoenergycoeffs(coeffs::Vector{T}, ρ::Number) where {T <: Real}
    a = one(T)
    b = zero(T)
    c = -arG(coeffs, ρ)
    return a, b, c
end

# Return the discriminant of `A`'s energy as a quadratic function
# of the topocentric range-rate evaluated at range `ρ`. `boundary`
# chooses between the `:outer` (default) or `:inner` boundary.
# See between equations (8.9)-(8.10) and (8.12)-(8.13)
# of https://doi.org/10.1017/CBO9781139175371
discriminant(a::Number, b::Number, c::Number) = b^2 - 4 * a * c
discriminant_derivatives(a::Number, b::Number, c::Number, dc::Number, d2c::Number) =
    discriminant(a, b, c), -4 * a * dc, -4 * a * d2c

arenergydis(A::AdmissibleRegion, ρ::Number, boundary::Symbol = :outer) =
    discriminant(arenergycoeffs(A, ρ, boundary)...)

_arhelenergydis(coeffs::Vector{T}, a_max::T, ρ::Number) where {T <: Real} =
    discriminant(_arhelenergycoeffs(coeffs, a_max, ρ)...)

_argeoenergydis(coeffs::Vector{T}, ρ::Number) where {T <: Real} =
    discriminant(_argeoenergycoeffs(coeffs, ρ)...)

# Return  a vector with the range-rates in the boundary of `A`
# for a given range `ρ`. `boundary` chooses between the `:outer`
# (default) or `:inner` boundary.
function rangerates(A::AdmissibleRegion, ρ::Number, boundary::Symbol = :outer)
    if boundary == :outer
        return _helrangerates(A.coeffs, A.a_max, ρ)
    elseif boundary == :inner
        return _georangerates(A.coeffs, ρ)
    else
        throw(ArgumentError("Argument `boundary` must be either `:outer` or `:inner`"))
    end
end

function _helrangerates(coeffs::Vector{T}, a_max::T, ρ::Number) where {T <: Real}
    a, b, c = _arhelenergycoeffs(coeffs, a_max, ρ)
    d = discriminant(a, b, c)
    # The number of solutions depends on the discriminant
    if d > 0
        return [(-b - sqrt(d))/(2a), (-b + sqrt(d))/(2a)]
    elseif d == 0
        return [-b/(2a)]
    else
        return Vector{T}(undef, 0)
    end
end

function _georangerates(coeffs::Vector{T}, ρ::Number) where {T <: Real}
    ρ0 = min(R_SI, G⁻¹0(coeffs))
    !(0 < ρ <= ρ0) && return Vector{T}(undef, 0)
    a, b, c = _argeoenergycoeffs(coeffs, ρ)
    d = discriminant(a, b, c)
    # The number of solutions depends on the discriminant
    if d > 0
        return [-sqrt(arG(coeffs, ρ)), sqrt(arG(coeffs, ρ))]
    elseif d == 0
        return zero(T) * [arG(coeffs, ρ)]
    else
        return Vector{T}(undef, 0)
    end
end

# Return a range-rate in the boundary of `A` for a given range `ρ`.
# `m = :min/:max` chooses which rate to return, while `boundary`
# chooses between the `:outer`(default) or `:inner` boundary.
function rangerate(A::AdmissibleRegion, ρ::Number, m::Symbol, boundary::Symbol = :outer)
    if boundary == :outer
        return _helrangerate(A.coeffs, A.a_max, ρ, m)
    elseif boundary == :inner
        return _georangerate(A.coeffs, ρ, m)
    else
        throw(ArgumentError("Argument `boundary` must be either `:outer` or `:inner`"))
    end
end

function _helrangerate(coeffs::Vector{T}, a_max::T, ρ::T, m::Symbol) where {T <: Real}
    a, b, c = _arhelenergycoeffs(coeffs, a_max, ρ)
    d = discriminant(a, b, c)
    @assert d > 0 "Less than two solutions, use rangerate(::AdmissibleRegion, ::Real,\
                   :outer) instead"
    # Choose min or max solution
    if m == :min
        return (-b - sqrt(d))/(2a)
    elseif m == :max
        return (-b + sqrt(d))/(2a)
    else
        throw(ArgumentError("Argument `m` must be either `:min` or `:max`"))
    end
end

function _helrangerate_derivatives(coeffs::Vector{T}, a_max::T, ρ::T,
                                   m::Symbol) where {T <: Real}
    # Choose between min or max sign
    if m == :min
        sgn = -1
    elseif m == :max
        sgn = +1
    else
        throw(ArgumentError("Argument `m` must be either `:min` or `:max`"))
    end
    # Outer boundary coefficients and its derivatives
    a, b, c, dc, d2c = _arhelenergycoeffs_derivatives(coeffs, a_max, ρ)
    # Discriminant and its derivatives
    dis, ddis, d2dis = discriminant_derivatives(a, b, c, dc, d2c)
    sqrtdis = sqrt(dis)
    # Range rate and its derivatives
    v_ρ = (-b + sgn * sqrtdis) / (2a)
    C = sgn / (4a)
    dv_ρ = C * ddis / sqrtdis
    d2v_ρ = C * (d2dis / sqrtdis - ddis^2 / (2sqrtdis^3))
    return v_ρ, dv_ρ, d2v_ρ
end

function _georangerate(coeffs::Vector{T}, ρ::T, m::Symbol) where {T <: Real}
    ρ0 = min(R_SI, G⁻¹0(coeffs))
    @assert 0 < ρ <= ρ0 "No solutions for geocentric energy outside 0 < ρ <= ρ0"
    a, b, c = _argeoenergycoeffs(coeffs, ρ)
    d = discriminant(a, b, c)
    @assert d > 0 "Less than two solutions, use rangerate(::AdmissibleRegion, ::Real,\
                   :inner) instead"
    # Choose min or max solution
    if m == :min
        return -sqrt(arG(coeffs, ρ))
    elseif m == :max
        return sqrt(arG(coeffs, ρ))
    else
        throw(ArgumentError("Argument `m` must be either `:min` or `:max`"))
    end
end

# Return the maximum possible range in the outer boundary of
# an admissible region with coefficients `coeffs` and maximum
# semimajor axis `a_max`.
function _helmaxrange(coeffs::Vector{T}, a_max::T) where {T <: Real}
    # Initial guess
    sol = find_zeros(s -> _arhelenergydis(coeffs, a_max, s), R_EA, 100.0)
    iszero(length(sol)) && return zero(T)
    ρ_max = sol[1]
    # Make sure U(ρ) ≥ 0 and there is at least one _helrangerate solution
    niter = 0
    while _arhelenergydis(coeffs, a_max, ρ_max) < 0 ||
          length(_helrangerates(coeffs, a_max, ρ_max)) == 0
        niter += 1
        ρ_max = prevfloat(ρ_max)
        niter > 1_000 && break
    end
    return ρ_max
end

# Return the maximum possible range in the inner boundary of
# an admissible region with coefficients `coeffs`.
function _geomaxrange(coeffs::Vector{T}) where {T <: Real}
    # Initial guess
    ρ_max = min(R_SI, G⁻¹0(coeffs))
    # Make sure G(ρ) ≥ 0 and there is at least one _georangerate solution
    niter = 0
    while _argeoenergydis(coeffs, ρ_max) < 0 ||
          length(_georangerates(coeffs, ρ_max)) == 0
        niter += 1
        ρ_max = prevfloat(ρ_max)
        niter > 1_000 && break
    end

    return ρ_max
end

# Parametrization of `A`'s
# - `:outer` boundary (default) with `t ∈ [0, 3]`, or
# - `:inner` boundary with `t ∈ [0, 2]`.
# `ρscale` sets the horizontal axis scale to `:linear` (default)
# or `:log`.
function arboundary(A::AdmissibleRegion, t::Number,
                    boundary::Symbol = :outer, ρscale = :linear)
    if boundary == :outer
        return _arhelboundary(A, t, ρscale)
    elseif boundary == :inner
        return _argeoboundary(A, t, ρscale)
    else
        throw(ArgumentError("Argument `boundary` must be either `:outer` or `:inner`"))
    end
end

function _arhelboundary(A::AdmissibleRegion, t::Number, ρscale::Symbol = :linear)
    # Parametrization domain
    @assert 0.0 <= t <= 3.0
    # Lower (upper) bounds
    if ρscale == :linear
        x_min, x_max = A.ρ_domain
    elseif ρscale == :log
        x_min, x_max = log10.(A.ρ_domain)
    else
        throw(ArgumentError("Argument `ρscale` must be either `:linear` or `:log`"))
    end
    y_min, y_max = A.v_ρ_domain
    # ρ = x_min
    if 0.0 <= t < 1.0
        x, v_ρ = x_min, y_min + t * (y_max - y_min)
    # Upper curve
    elseif 1.0 <= t < 2.0
        x = x_min + (t-1)*(x_max - x_min)
        _x_ = ρscale == :linear ? x : clamp(10^x, A.ρ_domain[1], A.ρ_domain[2])
        v_ρ = rangerates(A, _x_, :outer)[end]
    # Lower curve
    elseif 2.0 <= t <= 3.0
        x = x_max - (t-2)*(x_max - x_min)
        _x_ = ρscale == :linear ? x : clamp(10^x, A.ρ_domain[1], A.ρ_domain[2])
        v_ρ = rangerates(A, _x_, :outer)[1]
    end

    return [x, v_ρ]
end

function _argeoboundary(A::AdmissibleRegion, t::Number, ρscale::Symbol = :linear)
    # Parametrization domain
    @assert 0.0 <= t <= 2.0
    # Lower (upper) bounds
    ρ_max = _geomaxrange(A.coeffs)
    if ρscale == :linear
        x_min, x_max = A.ρ_domain[1], ρ_max
    elseif ρscale == :log
        x_min, x_max = log10(A.ρ_domain[1]), log10(ρ_max)
    else
        throw(ArgumentError("Argument `ρscale` must be either `:linear` or `:log`"))
    end
    # Upper curve
    if 0.0 <= t < 1.0
        x = x_min + t*(x_max - x_min)
        _x_ = ρscale == :linear ? x : clamp(10^x, A.ρ_domain[1], ρ_max)
        v_ρ = rangerates(A, _x_, :inner)[end]
    # Lower curve
    elseif 1.0 <= t <= 2.0
        x = x_max - (t-1)*(x_max - x_min)
        _x_ = ρscale == :linear ? x : clamp(10^x, A.ρ_domain[1], ρ_max)
        v_ρ = rangerates(A, _x_, :inner)[1]
    end

    return [x, v_ρ]
end

# Use golden section search to find the `m = :min/:max` range-rate in the
# boundary of `A` in the interval `[ρmin, ρmax]`. `boundary` chooses
# between the `:outer`(default) or `:inner` boundary and `tol` is the
# absolute tolerance (default: `1e-5`).
# Adapted from https://en.wikipedia.org/wiki/Golden-section_search
function argoldensearch(A::AdmissibleRegion{T}, ρmin::T, ρmax::T, m::Symbol,
                        boundary::Symbol = :outer, tol::T = 1e-5) where {T <: Real}
    # 1 / φ
    invphi = (sqrt(5) - 1) / 2
    # 1 / φ^2
    invphi2 = (3 - sqrt(5)) / 2
    # Interval bounds
    a, b = ρmin, ρmax
    # Interval width
    h = b - a
    # Termination condition
    if h <= tol
        ρ = (a + b) / 2
        return ρ, rangerate(A, ρ, m, boundary)
    end
    # Required steps to achieve tolerance
    n = ceil(Int, log(tol/h) / log(invphi))
    # Initialize center points
    c = a + invphi2 * h
    d = a + invphi * h
    yc = rangerate(A, c, m, boundary)
    yd = rangerate(A, d, m, boundary)
    # Main loop
    for _ in 1:n
        if (m == :min && yc < yd) || (m == :max && yc > yd)
            b = d
            d = c
            yd = yc
            h = invphi * h
            c = a + invphi2 * h
            yc = rangerate(A, c, m, boundary)
        else
            a = c
            c = d
            yc = yd
            h = invphi * h
            d = a + invphi * h
            yd = rangerate(A, d, m, boundary)
        end
    end

    if (m == :min && yc < yd) || (m == :max && yc > yd)
        ρ = (a + d) / 2
    else
        ρ = (c + b) / 2
    end

    return ρ, rangerate(A, ρ, m, boundary)
end