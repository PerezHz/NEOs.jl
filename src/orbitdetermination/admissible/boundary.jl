# Return the polynomial coefficients for an [`AdmissibleRegion`](@ref).
# See equation (8.8) of https://doi.org/10.1017/CBO9781139175371
function arcoeffs(::T, δ::T, v_α::T, v_δ::T, ρ::AbstractVector{T}, ρ_α::AbstractVector{T},
                  ρ_δ::AbstractVector{T}, q::AbstractVector{T}) where {T <: Number}
    coeffs = Vector{T}(undef, 6)
    coeffs[1] = dot3D(q[1:3], q[1:3])
    coeffs[2] = 2 * dot3D(q[4:6], ρ)
    coeffs[3] = v_α^2 * cos(δ)^2 + v_δ^2  # Proper motion squared
    coeffs[4] = 2 * v_α * dot3D(q[4:6], ρ_α) + 2 * v_δ * dot3D(q[4:6], ρ_δ)
    coeffs[5] = dot3D(q[4:6], q[4:6])
    coeffs[6] = 2 * dot3D(q[1:3], ρ)
    return coeffs
end

# Return the topocentric line-of-sight unit vector and its
# partial derivatives with respect to `α` and `δ`.
# See between equations (8.5) and (8.6) of https://doi.org/10.1017/CBO9781139175371
function topounitpdv(α::Number, δ::Number)
    sin_α, cos_α = sincos(α)
    sin_δ, cos_δ = sincos(δ)
    sin_α_sin_δ = sin_α * sin_δ
    sin_α_cos_δ = sin_α * cos_δ
    cos_α_sin_δ = cos_α * sin_δ
    cos_α_cos_δ = cos_α * cos_δ
    ρ = [cos_α_cos_δ, sin_α_cos_δ, sin_δ]
    ρ_α = [-sin_α_cos_δ, cos_α_cos_δ, zero(α)]
    ρ_δ = [-cos_α_sin_δ, -sin_α_sin_δ, cos_δ]
    return ρ, ρ_α, ρ_δ
end

# W function of an [`AdmissibleRegion`](@ref).
# See equation (8.9) of https://doi.org/10.1017/CBO9781139175371
arW(A::AdmissibleRegion, ρ::Number) = arW(A.coeffs, ρ)
arW(coeffs::AbstractVector, ρ::Number) = coeffs[3] * ρ^2 + coeffs[4] * ρ + coeffs[5]

ardW(A::AdmissibleRegion, ρ::Number) = ardW(A.coeffs, ρ)
ardW(coeffs::AbstractVector, ρ::Number) = 2 * coeffs[3] * ρ + coeffs[4]

ard2W(A::AdmissibleRegion, ρ::Number) = ard2W(A.coeffs, ρ)
ard2W(coeffs::AbstractVector, ::Number) = 2 * coeffs[3]

# S function of an [`AdmissibleRegion`](@ref).
# See equation (8.9) of https://doi.org/10.1017/CBO9781139175371
arS(A::AdmissibleRegion, ρ::Number) = arS(A.coeffs, ρ)
arS(coeffs::AbstractVector, ρ::Number) = ρ^2 + coeffs[6] * ρ + coeffs[1]

ardS(A::AdmissibleRegion, ρ::Number) = ardS(A.coeffs, ρ)
ardS(coeffs::AbstractVector, ρ::Number) = 2 * ρ + coeffs[6]

ard2S(A::AdmissibleRegion, ρ::Number) = ard2S(A.coeffs, ρ)
ard2S(::AbstractVector, ::Number) = 2

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

# Auxiliary function to compute the root of G(::AdmissibleRegion)
G⁻¹0(A::AdmissibleRegion) = G⁻¹0(A.coeffs)
G⁻¹0(coeffs::AbstractVector) = cbrt(2 * k_gauss^2 * μ_ES / coeffs[3])

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

function _arhelenergycoeffs(coeffs::AbstractVector, a_max::Real, ρ::Number)
    a = one(eltype(coeffs))
    b = coeffs[2]
    c = arW(coeffs, ρ) + k_gauss^2 * (1/a_max - 2/sqrt(arS(coeffs, ρ)))
    return a, b, c
end

function _arhelenergycoeffs_derivatives(coeffs::AbstractVector, a_max::Real, ρ::Number)
    a = one(eltype(coeffs))
    b = coeffs[2]
    W, dW, d2W = arW(coeffs, ρ), ardW(coeffs, ρ), ard2W(coeffs, ρ)
    S, dS, d2S = arS(coeffs, ρ), ardS(coeffs, ρ), ard2S(coeffs, ρ)
    sqrtS = sqrt(S)
    c = W + k_gauss^2 * (1/a_max - 2/sqrtS)
    dc = dW + k_gauss^2 * dS / sqrtS^3
    d2c = d2W + k_gauss^2 * (d2S / sqrtS^3 - 3 * dS^2 / (2 * sqrtS^5))
    return a, b, c, dc, d2c
end

function _argeoenergycoeffs(coeffs::AbstractVector, ρ::Number)
    a = one(eltype(coeffs))
    b = zero(eltype(coeffs))
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

_arhelenergydis(coeffs::AbstractVector, a_max::Real, ρ::Number) =
    discriminant(_arhelenergycoeffs(coeffs, a_max, ρ)...)

_argeoenergydis(coeffs::AbstractVector, ρ::Number) =
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

function _helrangerates(coeffs::AbstractVector, a_max::Real, ρ::Number)
    a, b, c = _arhelenergycoeffs(coeffs, a_max, ρ)
    d = discriminant(a, b, c)
    # The number of solutions depends on the discriminant
    if d > 0
        return [(-b - sqrt(d))/(2a), (-b + sqrt(d))/(2a)]
    elseif d == 0
        return [-b/(2a) * one(d)]
    else # d < 0
        return Vector{typeof(d)}(undef, 0)
    end
end

function _georangerates(coeffs::AbstractVector, ρ::Number)
    a, b, c = _argeoenergycoeffs(coeffs, ρ)
    d = discriminant(a, b, c)
    # The number of solutions depends on the discriminant
    if !(0 < ρ ≤ min(R_SI, G⁻¹0(coeffs))) || d < 0
        return Vector{typeof(d)}(undef, 0)
    elseif d > 0
        return [-sqrt(arG(coeffs, ρ)), sqrt(arG(coeffs, ρ))]
    else # d == 0
        return [zero(d)]
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

function _helrangerate(coeffs::AbstractVector, a_max::Real, ρ::Number, m::Symbol)
    a, b, c = _arhelenergycoeffs(coeffs, a_max, ρ)
    d = discriminant(a, b, c)
    @assert d > 0 "Less than two solutions, use rangerates(::AdmissibleRegion, \
        ::Real, :outer) instead"
    # Choose min or max solution
    if m == :min
        return (-b - sqrt(d))/(2a)
    elseif m == :max
        return (-b + sqrt(d))/(2a)
    else
        throw(ArgumentError("Argument `m` must be either `:min` or `:max`"))
    end
end

function _helrangerate_derivatives(coeffs::AbstractVector, a_max::Real,
                                   ρ::Number, m::Symbol)
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

function _georangerate(coeffs::AbstractVector, ρ::Number, m::Symbol)
    ρ0 = min(R_SI, G⁻¹0(coeffs))
    @assert 0 < ρ <= ρ0 "No solutions for geocentric energy outside 0 < ρ <= ρ0"
    a, b, c = _argeoenergycoeffs(coeffs, ρ)
    d = discriminant(a, b, c)
    @assert d > 0 "Less than two solutions, use rangerates(::AdmissibleRegion, \
        ::Real, :inner) instead"
    # Choose min or max solution
    if m == :min
        return -sqrt(arG(coeffs, ρ))
    elseif m == :max
        return sqrt(arG(coeffs, ρ))
    else
        throw(ArgumentError("Argument `m` must be either `:min` or `:max`"))
    end
end

# Return the smallest (largest) float that satisfies a
# condition given by the function `f::Bool`. `a` and `b`
# must be such that `a < b` and `f(a) = false`, `f(b) = true`
# (`f(a) = true`, `f(b) = false`).
function smallestfloat(f, a::Real, b::Real)
    c = b
    while a < prevfloat(b)
        c = (a + b) / 2
        if f(c)
            b = c
        else
            a = c
        end
    end
    return b
end

function largestfloat(f, a::Real, b::Real)
    c = a
    while nextfloat(a) < b
        c = (a + b) / 2
        if f(c)
            a = c
        else
            b = c
        end
    end
    return a
end

# Check whether a given range `ρ` is inside either of the connected
# components of the boundaries of an admissible region with coefficients
# `coeffs` and maximum semimajor axis `a_max`
_arhelin(coeffs::AbstractVector, a_max::Real, ρ::Number) =
    _arhelenergydis(coeffs, a_max, ρ) ≥ 0 && !isempty(_helrangerates(coeffs, a_max, ρ))
_argeoin(coeffs::AbstractVector, ρ::Number) =
    _argeoenergydis(coeffs, ρ) ≥ 0 && !isempty(_georangerates(coeffs, ρ))

# Return the minimum and maximum ranges of all the connected
# components in the outer boundary of an admissible region
# with coefficients `coeffs`, maximum semimajor axis `a_max`,
# apparent and absolute magnitudes `h` and `H_max`, respectively,
# and slope parameter `slope`. `ϵ` is a small numerical offset
# used to obtain an enclosing interval for `smallest(largest)float`
function _helrangedomain(coeffs::AbstractVector, a_max::Real, h::Real, H_max::Real;
                         slope::Real = 0.15, ϵ::Real = 1E-4)
    # Find the roots of the heliocentric energy discriminant
    ρs = find_zeros(ρ -> _arhelenergydis(coeffs, a_max, ρ), R_EA, HELIOPAUSE_RADIUS)
    isempty(ρs) && return ρs
    # Find the maximum range of the first component
    flag = _arhelenergydis(coeffs, a_max, ρs[1]) ≥ 0
    ρa, ρb = (ρs[1] - !flag*ϵ, ρs[1] + flag*ϵ)
    ρs[1] = largestfloat(ρ -> _arhelin(coeffs, a_max, ρ), ρa, ρb)
    # Find the minimum range of the first component
    if isnan(h)
        # Earth's sphere of influence radius / Earth's physical radius
        pushfirst!(ρs, R_SI < ρs[1] ? R_SI : R_EA)
    else
        # Tiny object boundary
        pushfirst!(ρs, body2observer(coeffs, h, H_max; slope))
    end
    # One component and possibly a second degenerated component
    if length(ρs) < 4
        return ρs
    # Two non degenerated components
    else
        # Find the minimum/maximum range of the second component
        flag = _arhelenergydis(coeffs, a_max, ρs[3]) ≥ 0
        ρa, ρb = (ρs[3] - flag*ϵ, ρs[3] + !flag*ϵ)
        ρs[3] = smallestfloat(ρ -> _arhelin(coeffs, a_max, ρ), ρa, ρb)
        flag = _arhelenergydis(coeffs, a_max, ρs[4]) ≥ 0
        ρa, ρb = (ρs[4] - !flag*ϵ, ρs[4] + flag*ϵ)
        ρs[4] = largestfloat(ρ -> _arhelin(coeffs, a_max, ρ), ρa, ρb)
        return ρs
    end
end

# Return the maximum range in the inner boundary of
# an admissible region with coefficients `coeffs`. `ϵ`
# is a small numerical offset used to obtain an enclosing
# interval for `largestfloat`
function _geomaxrange(coeffs::AbstractVector; ϵ::Real = 1E-4)
    ρmax = min(R_SI, G⁻¹0(coeffs))
    flag = _argeoenergydis(coeffs, ρmax) ≥ 0
    ρa, ρb = (ρmax - !flag*ϵ, ρmax + flag*ϵ)
    ρmax = largestfloat(Base.Fix1(_argeoin, coeffs), ρa, ρb)
    return ρmax
end

# Parametrization of the interval `[a, b]` with parameter
# `t ∈ [0, 1]`. The evaluation can be made from left to
# right (`right = true`) or the other way around (`right = false`).
numberbetween(a::Real, b::Real, right::Bool, t::Real) =
    right ? a + t * (b - a) : b - t * (b - a)

# Return the domain of the parameter that parametrizes the
# `:outer` or `:inner` boundary of an admissible region.
# Such domain depends on the following criteria:
# - for the `:outer` boundary, `t ∈ [0, tmax]` where `tmax = 3`
# if `A` has a single connected component and `tmax = 5` otherwise,
# - for the `:inner` boundary, `t ∈ [0, 2]`.
boundarydomain(A::AdmissibleRegion, ::Val{:outer}) =
    numberofcomponents(A) < 2 ? (0.0, 3.0) : (0.0, 5.0)
boundarydomain(::AdmissibleRegion, ::Val{:inner}) = (0.0, 2.0)

# Return a point in the boundary of `A` for a given parameter `t`.
# `boundary` chooses between the `:outer` (default) or `:inner`
# boundary, while `ρscale` sets the horizontal axis scale to
# `:linear` (default) of `:log`.
function arboundary(A::AdmissibleRegion, t::Number, boundary::Symbol = :outer,
                    ρscale::Symbol = :linear)
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
    tmin, tmax = boundarydomain(A, Val(:outer))
    @assert tmin <= t <= tmax
    # Number of components
    Nc = numberofcomponents(A)
    # Lower (upper) bounds
    if ρscale == :linear
        x_domain = rangedomain(A)
    elseif ρscale == :log
        x_domain = log10.(rangedomain(A))
    else
        throw(ArgumentError("Argument `ρscale` must be either `:linear` or `:log`"))
    end
    ydomain = rangeratedomain(A)
    # Tiny object boundary
    if 0.0 ≤ t < 1.0
        x, y = x_domain[1], numberbetween(ydomain[1], ydomain[2], true, t)
    # First component
    elseif 1.0 ≤ t && ifelse(Nc == 1, t ≤ 3.0, t < 3.0)
        flag = 1.0 ≤ t < 2.0
        _t_ = flag ? t - 1 : t - 2
        x = numberbetween(x_domain[1], x_domain[2], flag, _t_)
        _x_ = ρscale == :linear ? x : clamp(10^x, A.ρ_domain[1], A.ρ_domain[2])
        ys = rangerates(A, _x_, :outer)
        y = flag ? last(ys) : first(ys)
    # Second component
    elseif 3.0 ≤ t ≤ 5.0
        flag = 3.0 ≤ t < 4.0
        _t_ = flag ? t - 3 : t - 4
        x = numberbetween(x_domain[3], x_domain[4], flag, _t_)
        _x_ = ρscale == :linear ? x : clamp(10^x, A.ρ_domain[3], A.ρ_domain[4])
        ys = rangerates(A, _x_, :outer)
        y = flag ? last(ys) : first(ys)
    end
    return [x, y]
end

function _argeoboundary(A::AdmissibleRegion, t::Number, ρscale::Symbol = :linear)
    # Parametrization domain
    tmin, tmax = boundarydomain(A, Val(:inner))
    @assert tmin <= t <= tmax
    # Lower (upper) bounds
    ρmax = _geomaxrange(A.coeffs)
    if ρscale == :linear
        xmin, xmax = A.ρ_domain[1], ρmax
    elseif ρscale == :log
        xmin, xmax = log10(A.ρ_domain[1]), log10(ρmax)
    else
        throw(ArgumentError("Argument `ρscale` must be either `:linear` or `:log`"))
    end
    flag = 0.0 ≤ t < 1.0
    _t_ = flag ? t : t - 1
    x = numberbetween(xmin, xmax, flag, _t_)
    _x_ = ρscale == :linear ? x : clamp(10^x, A.ρ_domain[1], ρmax)
    ys = rangerates(A, _x_, :inner)
    y = flag ? last(ys) : first(ys)
    return [x, y]
end

# Use golden section search to find the `m = :min/:max` range-rate in the
# boundary of `A` in the interval `[ρmin, ρmax]`. `boundary` chooses
# between the `:outer`(default) or `:inner` boundary and `tol` is the
# absolute tolerance (default: `1E-5`).
# Adapted from https://en.wikipedia.org/wiki/Golden-section_search
function argoldensearch(A::AdmissibleRegion{T}, ρmin::T, ρmax::T, m::Symbol,
                        boundary::Symbol = :outer, tol::T = 1E-5) where {T <: Real}
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

# Angle between the line of sight and the opposition direction
# See paragraph below equation (8.16) of https://doi.org/10.1017/CBO9781139175371
opposition_angle(coeffs::AbstractVector) = acos(coeffs[6] / (2 * sqrt(coeffs[1])))
opposition_angle(A::AdmissibleRegion) = opposition_angle(A.coeffs)

# Approximation for the distance between the body and the observer
# See paragraph below equation (8.16) of https://doi.org/10.1017/CBO9781139175371
function body2observer(coeffs::AbstractVector, h::Number, H::Number;
                       slope::Number = 0.15)
    β = opposition_angle(coeffs)
    Φ = phase_integral(β; slope)
    return 10^((h - H + 2.5*log10(Φ))/5)
end

body2observer(x::AdmissibleRegion, H::Number) = body2observer(x.coeffs, mag(x), H;
    slope = slopeparameter(x))

# Check whether a point P is inside A's boundary
function in(P::Union{AbstractVector, Tuple{<:Real, <:Real}}, A::AdmissibleRegion)
    @assert length(P) == 2 "Points in admissible region are of dimension 2"
    ρ_domain = rangedomain(A)
    if ρ_domain[1] ≤ P[1] ≤ ρ_domain[2] ||
        (numberofcomponents(A) > 1 && ρ_domain[3] ≤ P[1] ≤ ρ_domain[4])
        ys = rangerates(A, P[1], :outer)
        if length(ys) == 1
            return P[2] == ys[1]
        else
            return ys[1] ≤ P[2] ≤ ys[2]
        end
    else
        return false
    end
end