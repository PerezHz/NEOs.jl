@doc raw"""
    AdmissibleRegion{T <: Real}

Subset of topocentric range × range-rate space defined by
the following constraints:
- heliocentric energy ≤ `k_gauss^2/(2a_max)`,
- absolute magnitude ≤ `H_max`,
- geocentric energy ≥ `0`.

# Fields

- `date::DateTime`: time of observation.
- `α::T`: right ascension.
- `δ::T`: declination.
- `v_α::T`: right ascension velocity.
- `v_δ::T`: declination velocity.
- `H_max::T`: maximum absolute magnitude.
- `a_max::T`: maximum semimajor axis.
- `ρ_unit, ρ_α, ρ_δ::Vector{T}`: topocentric unit vector and its partials.
- `q::Vector{T}`: heliocentric position of observer.
- `coeffs::Vector{T}`: polynomial coefficients.
- `ρ_domain::Vector{T}`: range domain.
- `v_ρ_domain::Vector{T}`: range-rate domain.
- `Fs::Matrix{T}`: boundary points.
- `observatory::ObservatoryMPC{T}`: observing station.

!!! reference
    See Chapter 8 of https://doi.org/10.1017/CBO9781139175371, or
    https://doi.org/10.1007/s10569-004-6593-5.
"""
@auto_hash_equals struct AdmissibleRegion{T <: Real}
    date::DateTime
    α::T
    δ::T
    v_α::T
    v_δ::T
    H_max::T
    a_max::T
    ρ_unit::Vector{T}
    ρ_α::Vector{T}
    ρ_δ::Vector{T}
    q::Vector{T}
    coeffs::Vector{T}
    ρ_domain::Vector{T}
    v_ρ_domain::Vector{T}
    Fs::Matrix{T}
    observatory::ObservatoryMPC{T}
end

# Definition of zero AdmissibleRegion{T}
function zero(::Type{AdmissibleRegion{T}}) where {T <: Real}
    return AdmissibleRegion{T}(
        DateTime(2000), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T),
        Vector{T}(undef, 0), Vector{T}(undef, 0), Vector{T}(undef, 0),
        Vector{T}(undef, 0), Vector{T}(undef, 0), Vector{T}(undef, 0),
        Vector{T}(undef, 0), Matrix{T}(undef, 0, 0), unknownobs()
    )
end

iszero(x::AdmissibleRegion{T}) where {T <: Real} = x == zero(AdmissibleRegion{T})

# Print method for AdmissibleRegion
# Examples:
# AE: [11.55523, 13.29296, -1.01625, -0.55432] t: 2019-10-23T05:23:27.384 obs: Palomar Mountain--ZTF
# AE: [358.56604, 1.25546, -2.05305, -1.91538] t: 2019-11-01T09:03:26.280 obs: Pan-STARRS 1, Haleakala
function show(io::IO, A::AdmissibleRegion{T}) where {T <: Real}
    v = string(
        @sprintf("%.5f", rad2deg(A.α)), ", ",
        @sprintf("%.5f", rad2deg(A.δ)), ", ",
        @sprintf("%.5f", rad2deg(A.v_α)), ", ",
        @sprintf("%.5f", rad2deg(A.v_δ)), "",
    )
    print(io, "AE: [", v, "]", " t: ", A.date, " obs: ", A.observatory.name)
end

# Outer constructor
function AdmissibleRegion(tracklet::Tracklet{T}, params::NEOParameters{T}) where {T <: Real}
    # Unfold
    obs, t_datetime, α, δ = observatory(tracklet), date(tracklet), ra(tracklet), dec(tracklet)
    v_α, v_δ, h = vra(tracklet), vdec(tracklet), mag(tracklet)
    H_max, a_max = params.H_max, params.a_max
    # Topocentric unit vector and partials
    ρ, ρ_α, ρ_δ = topounitpdv(α, δ)
    # Time of observation [days since J2000]
    t_days = datetime2days(t_datetime)
    # Time of observation [et seconds]
    t_et = datetime2et(t_datetime)
    # Heliocentric position of the observer
    q = params.eph_ea(t_days) + kmsec2auday(obsposvelECI(obs, t_et)) - params.eph_su(t_days)
    # Admissible region coefficients
    coeffs = arcoeffs(α, δ, v_α, v_δ, ρ, ρ_α, ρ_δ, q)
    # Maximum range (heliocentric energy constraint)
    ρ_max = _helmaxrange(coeffs, a_max)
    iszero(ρ_max) && return zero(AdmissibleRegion{T})
    # Minimum range
    if isnan(h)
        if R_SI < ρ_max
            # Earth's sphere of influence radius
            ρ_min = R_SI
        else
            # Earth's physical radius
            ρ_min = R_EA
        end
    else
        # Tiny object boundary
        ρ_min = 10^((h - H_max)/5)
    end
    ρ_min > ρ_max && return zero(AdmissibleRegion{T})
    # Range domain
    ρ_domain = [ρ_min, ρ_max]
    # Range-rate domain
    v_ρ_min, v_ρ_max = _helrangerates(coeffs, a_max, ρ_min)[1:2]
    v_ρ_domain = [v_ρ_min, v_ρ_max]
    # Range rate symmetry level
    v_ρ_mid = _helrangerates(coeffs, a_max, ρ_max)[1]
    # Boundary points
    Fs = Matrix{T}(undef, 3, 2)
    Fs[1, :] .= [ρ_min, v_ρ_min]
    Fs[2, :] .= [ρ_min, v_ρ_max]
    Fs[3, :] .= [ρ_max, v_ρ_mid]

    return AdmissibleRegion{T}(t_datetime, α, δ, v_α, v_δ, H_max, a_max,
                               ρ, ρ_α, ρ_δ, q, coeffs, ρ_domain, v_ρ_domain,
                               Fs, obs)
end

@doc raw"""
    arcoeffs(α::T, δ::T, v_α::T, v_δ::T, ρ::Vector{T},
             ρ_α::Vector{T}, ρ_δ::Vector{T}, q::Vector{T}) where {T <: Number}

Return the polynomial coefficients for an [`AdmissibleRegion`](@ref).

!!! reference
    See equation (8.8) of https://doi.org/10.1017/CBO9781139175371.
"""
function arcoeffs(α::T, δ::T, v_α::T, v_δ::T, ρ::Vector{T},
                  ρ_α::Vector{T}, ρ_δ::Vector{T}, q::Vector{T}) where {T <: Number}
    coeffs = Vector{T}(undef, 6)
    coeffs[1] = dot(q[1:3], q[1:3])
    coeffs[2] = 2 *  dot(q[4:6], ρ)
    coeffs[3] = v_α^2 * cos(δ)^2 + v_δ^2  # proper motion squared
    coeffs[4] = 2 * v_α * dot(q[4:6], ρ_α) + 2 * v_δ * dot(q[4:6], ρ_δ)
    coeffs[5] = dot(q[4:6], q[4:6])
    coeffs[6] = 2*dot(q[1:3], ρ)
    return coeffs
end

@doc raw"""
    topounitpdv(α::T, δ::T) where {T <: Number}

Return the topocentric line-of-sight unit vector and its
partial derivatives with respect to `α` and `δ`.

!!! reference
    See between equations (8.5) and (8.6) of
    https://doi.org/10.1017/CBO9781139175371.
"""
function topounitpdv(α::T, δ::T) where {T <: Number}
    sin_α, cos_α = sincos(α)
    sin_δ, cos_δ = sincos(δ)

    ρ = [cos_δ * cos_α, cos_δ * sin_α, sin_δ]
    ρ_α = [-sin_α * cos_δ, cos_α * cos_δ, zero(T)]
    ρ_δ = [-cos_α * sin_δ, -sin_α * sin_δ, cos_δ]

    return ρ, ρ_α, ρ_δ
end

@doc raw"""
    arW(A::AdmissibleRegion{T}, ρ::S) where {T <: Real, S <: Number}

W function of an [`AdmissibleRegion`](@ref).

!!! reference
    See equation (8.9) of https://doi.org/10.1017/CBO9781139175371.
"""
function arW(coeffs::Vector{T}, ρ::S) where {T <: Real, S <: Number}
    return coeffs[3] * ρ^2 + coeffs[4] * ρ + coeffs[5]
end
arW(A::AdmissibleRegion{T}, ρ::S) where {T <: Real, S <: Number}  = arW(A.coeffs, ρ)

@doc raw"""
    arS(A::AdmissibleRegion{T}, ρ::S) where {T <: Real, S <: Number}

S function of an [`AdmissibleRegion`](@ref).

!!! reference
    See equation (8.9) of https://doi.org/10.1017/CBO9781139175371.
"""
function arS(coeffs::Vector{T}, ρ::S) where {T <: Real, S <: Number}
    return ρ^2 + coeffs[6] * ρ + coeffs[1]
end
arS(A::AdmissibleRegion{T}, ρ::S) where {T <: Real, S <: Number} = arS(A.coeffs, ρ)

@doc raw"""
    arG(A::AdmissibleRegion{T}, ρ::S) where {T <: Real, S <: Number}

G function of an [`AdmissibleRegion`](@ref).

!!! reference
    See equation (8.13) of https://doi.org/10.1017/CBO9781139175371.
"""
function arG(coeffs::Vector{T}, ρ::S) where {T <: Real, S <: Number}
    return 2 * k_gauss^2 * μ_ES / ρ - coeffs[3] * ρ^2
end
arG(A::AdmissibleRegion{T}, ρ::S) where {T <: Real, S <: Number} = arG(A.coeffs, ρ)

@doc raw"""
    arenergycoeffs(A::AdmissibleRegion{T}, ρ::S, [boundary::Symbol]) where {T <: Real, S <: Number}

Return the coefficients of `A`'s energy as a quadratic function
of the topocentric range-rate evaluated at range `ρ`. `boundary`
chooses between the `:outer` (default) or `:inner` boundary.

!!! reference
    See between equations (8.9)-(8.10) and (8.12)-(8.13)
    of https://doi.org/10.1017/CBO9781139175371.
"""
function arenergycoeffs(A::AdmissibleRegion{T}, ρ::S,
                        boundary::Symbol = :outer) where {T <: Real, S <: Number}
    if boundary == :outer
        return _arhelenergycoeffs(A.coeffs, A.a_max, ρ)
    elseif boundary == :inner
        return _argeoenergycoeffs(A.coeffs, ρ)
    else
        throw(ArgumentError("Argument `boundary` must be either `:outer` or `:inner`"))
    end
end

function _arhelenergycoeffs(coeffs::Vector{T}, a_max::T, ρ::S) where {T <: Real, S <: Number}
    a = one(T)
    b = coeffs[2]
    c = arW(coeffs, ρ) + k_gauss^2 * (1/a_max - 2/sqrt(arS(coeffs, ρ)))
    return a, b, c
end

function _argeoenergycoeffs(coeffs::Vector{T}, ρ::S) where {T <: Real, S <: Number}
    a = one(T)
    b = zero(T)
    c = -arG(coeffs, ρ)
    return a, b, c
end

@doc raw"""
    arenergydis(A::AdmissibleRegion{T}, ρ::S, [boundary::Symbol]) where {T <: Real, S <: Number}

Return the discriminant of `A`'s energy as a quadratic function
of the topocentric range-rate evaluated at range `ρ`. `boundary`
chooses between the `:outer` (default) or `:inner` boundary.

!!! reference
    See between equations (8.9)-(8.10) and (8.12)-(8.13)
    of https://doi.org/10.1017/CBO9781139175371.
"""
function arenergydis(A::AdmissibleRegion{T}, ρ::S,
                     boundary::Symbol = :outer) where {T <: Real, S <: Number}
    a, b, c = arenergycoeffs(A, ρ, boundary)
    return b^2 - 4*a*c
end

function _arhelenergydis(coeffs::Vector{T}, a_max::T, ρ::S) where {T <: Real, S <: Number}
    a, b, c = _arhelenergycoeffs(coeffs, a_max, ρ)
    return b^2 - 4*a*c
end

function _argeoenergydis(coeffs::Vector{T}, ρ::S) where {T <: Real, S <: Number}
    a, b, c = _argeoenergycoeffs(coeffs, ρ)
    return b^2 - 4*a*c
end

@doc raw"""
    rangerates(A::AdmissibleRegion{T}, ρ::S, [boundary::Symbol]) where {T, S <: Real}

Return  a vector with the range-rates in the boundary of `A`
for a given range `ρ`. `boundary` chooses between the `:outer`
(default) or `:inner` boundary.
"""
function rangerates(A::AdmissibleRegion{T}, ρ::S,
                   boundary::Symbol = :outer) where {T, S <: Real}
    if boundary == :outer
        return _helrangerates(A.coeffs, A.a_max, ρ)
    elseif boundary == :inner
        return _georangerates(A.coeffs, ρ)
    else
        throw(ArgumentError("Argument `boundary` must be either `:outer` or `:inner`"))
    end
end

function _helrangerates(coeffs::Vector{T}, a_max::T, ρ::S) where {T, S <: Real}
    a, b, c = _arhelenergycoeffs(coeffs, a_max, ρ)
    d = b^2 - 4*a*c
    # The number of solutions depends on the discriminant
    if d > 0
        return [(-b - sqrt(d))/(2a), (-b + sqrt(d))/(2a)]
    elseif d == 0
        return [-b/(2a)]
    else
        return Vector{T}(undef, 0)
    end
end

function _georangerates(coeffs::Vector{T}, ρ::S) where {T, S <: Real}
    ρ0 = min(R_SI, cbrt(2 * k_gauss^2 * μ_ES / coeffs[3]))
    !(0 < ρ <= ρ0) && return Vector{T}(undef, 0)
    a, b, c = _argeoenergycoeffs(coeffs, ρ)
    d = b^2 - 4*a*c
    # The number of solutions depends on the discriminant
    if d > 0
        return [-sqrt(arG(coeffs, ρ)), sqrt(arG(coeffs, ρ))]
    elseif d == 0
        return zero(T) * arG(coeffs, ρ)
    else
        return Vector{T}(undef, 0)
    end
end

@doc raw"""
    rangerate(A::AdmissibleRegion{T}, ρ::T, m::Symbol, [boundary::Symbol]) where {T <: Real}

Return a range-rate in the boundary of `A` for a given range `ρ`.
`m = :min/:max` chooses which rate to return, while `boundary`
chooses between the `:outer`(default) or `:inner` boundary.
"""
function rangerate(A::AdmissibleRegion{T}, ρ::S, m::Symbol,
                   boundary::Symbol = :outer) where {T, S <: Real}
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
    d = b^2 - 4*a*c
    @assert d > 0 "Less than two solutions, use rangerate(::AdmissibleRegion, ::Real, :outer) instead"
    # Choose min or max solution
    if m == :min
        return (-b - sqrt(d))/(2a)
    elseif m == :max
        return (-b + sqrt(d))/(2a)
    else
        throw(ArgumentError("Argument `m` must be either `:min` or `:max`"))
    end
end

function _georangerate(coeffs::Vector{T}, ρ::T, m::Symbol) where {T <: Real}
    ρ0 = min(R_SI, cbrt(2 * k_gauss^2 * μ_ES / coeffs[3]))
    @assert 0 < ρ <= ρ0 "No solutions for geocentric energy outside 0 < ρ <= ρ0"
    a, b, c = _argeoenergycoeffs(coeffs, ρ)
    d = b^2 - 4*a*c
    @assert d > 0 "Less than two solutions, use rangerate(::AdmissibleRegion, ::Real, :inner) instead"
    # Choose min or max solution
    if m == :min
        return -sqrt(arG(coeffs, ρ))
    elseif m == :max
        return sqrt(arG(coeffs, ρ))
    else
        throw(ArgumentError("Argument `m` must be either `:min` or `:max`"))
    end
end

@doc raw"""
    argoldensearch(A::AdmissibleRegion{T}, ρmin::T, ρmax::T, m::Symbol,
                   [boundary::Symbol,] [tol::T]) where {T <: Real}

Use golden section search to find the `m = :min/:max` range-rate in the
boundary of `A` in the interval `[ρmin, ρmax]`. `boundary` chooses
between the `:outer`(default) or `:inner` boundary and `tol` is the
absolute tolerance (default: `1e-5`).

!!! reference
    Adapted from https://en.wikipedia.org/wiki/Golden-section_search.
"""
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

@doc raw"""
    _helmaxrange(coeffs::Vector{T}, a_max::T) where {T <: Real}

Return the maximum possible range in the outer boundary of
an admissible region with coefficients `coeffs` and maximum
semimajor axis `a_max`.
"""
function _helmaxrange(coeffs::Vector{T}, a_max::T) where {T <: Real}
    # Initial guess
    sol = find_zeros(s -> _arhelenergydis(coeffs, a_max, s), R_EA, 100.0)
    iszero(length(sol)) && return zero(T)
    ρ_max = sol[1]
    # Make sure U(ρ) ≥ 0 and there is at least one _helrangerate solution
    niter = 0
    while _arhelenergydis(coeffs, a_max, ρ_max) < 0 || length(_helrangerates(coeffs, a_max, ρ_max)) == 0
        niter += 1
        ρ_max = prevfloat(ρ_max)
        niter > 1_000 && break
    end

    return ρ_max
end

@doc raw"""
    _geomaxrange(coeffs::Vector{T}) where {T <: Real}

Return the maximum possible range in the inner boundary of
an admissible region with coefficients `coeffs`.
"""
function _geomaxrange(coeffs::Vector{T}) where {T <: Real}
    # Initial guess
    ρ_max = min(R_SI, cbrt(2 * k_gauss^2 * μ_ES / coeffs[3]))
    # Make sure G(ρ) ≥ 0 and there is at least one _georangerate solution
    niter = 0
    while _argeoenergydis(coeffs, ρ_max) < 0 || length(_georangerates(coeffs, ρ_max)) == 0
        niter += 1
        ρ_max = prevfloat(ρ_max)
        niter > 1_000 && break
    end

    return ρ_max
end

@doc raw"""
    arboundary(A::AdmissibleRegion{T}, t::S, [boundary::Symbol,
        ρscale::Symbol]) where {T <: Real, S <: Number}

Parametrization of `A`'s
- `:outer` boundary (default) with `t ∈ [0, 3]`, or
- `:inner` boundary with `t ∈ [0, 2]`.
`ρscale` sets the horizontal axis scale to `:linear` (default)
or `:log`.
"""
function arboundary(A::AdmissibleRegion{T}, t::S, boundary::Symbol = :outer,
        ρscale = :linear) where {T <: Real, S <: Number}
    if boundary == :outer
        return _arhelboundary(A, t, ρscale)
    elseif boundary == :inner
        return _argeoboundary(A, t, ρscale)
    else
        throw(ArgumentError("Argument `boundary` must be either `:outer` or `:inner`"))
    end
end

function _arhelboundary(A::AdmissibleRegion{T}, t::S,
        ρscale::Symbol = :linear) where {T <: Real, S <: Number}
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

function _argeoboundary(A::AdmissibleRegion{T}, t::S,
        ρscale::Symbol = :linear) where {T <: Real, S <: Number}
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

@doc raw"""
    boundary_projection(A::AdmissibleRegion{T}, ρ::T, v_ρ::T) where {T <: Real}

Project `[ρ, v_ρ]` into `A`'s outer boundary.
"""
function boundary_projection(A::AdmissibleRegion{T}, ρ::T, v_ρ::T) where {T <: Real}
    # Project range
    ρ = clamp(ρ,  A.ρ_domain[1], A.ρ_domain[2])
    # Project range-rate
    y_domain = rangerates(A, ρ, :outer)
    if iszero(length(y_domain))
        v_ρ = sum(A.v_ρ_domain) / 2
    elseif isone(length(y_domain))
        v_ρ = y_domain[1]
    else
        v_ρ = clamp(v_ρ, y_domain[1], y_domain[2])
    end

    return ρ, v_ρ
end

# Check whether P is inside A's boundary
for U in (:(AbstractVector{T}), :(Tuple{T, T}))
    @eval begin
        function in(P::$U, A::AdmissibleRegion{T}) where {T <: Real}
            @assert length(P) == 2 "Points in admissible region are of dimension 2"
            if A.ρ_domain[1] <= P[1] <= A.ρ_domain[2]
                y_range = rangerates(A, P[1], :outer)
                if length(y_range) == 1
                    return P[2] == y_range[1]
                else
                    return y_range[1] <= P[2] <= y_range[2]
                end
            else
                return false
            end
        end
    end
end

@doc raw"""
    topo2bary(A::AdmissibleRegion{T}, ρ::U, v_ρ::U) where {T <: Real, U <: Number}

Convert topocentric range/range-rate `[ρ, v_ρ]` to barycentric cartesian coordinates.
`A` fixes the line of sight.
"""
function topo2bary(A::AdmissibleRegion{T}, ρ::U, v_ρ::U) where {T <: Real, U <: Number}
    # Barycentric position
    r = A.q[1:3] + ρ * A.ρ_unit + sseph(su, datetime2days(A.date))[1:3]
    # Barycentric velocity
    v = A.q[4:6] + v_ρ * A.ρ_unit + ρ * A.v_α * A.ρ_α + ρ * A.v_δ * A.ρ_δ
        + sseph(su, datetime2days(A.date))[4:6]
    # Barycentric state vector
    return vcat(r, v)
end

@doc raw"""
    bary2topo(A::AdmissibleRegion{T}, q0::Vector{U}) where {T <: Real, U <: Number}

Convert barycentric cartesian coordinates `q0` to topocentric range/range-rate.
`A` fixes the line of sight.
"""
function bary2topo(A::AdmissibleRegion{T}, q0::Vector{U}) where {T <: Real, U <: Number}
    # Heliocentric state vector
    r = q0 - sseph(su, datetime2days(A.date))
    # Topocentric range
    ρ = euclid3D(r - A.q)
    # Topocentric range rate
    v_ρ = dot3D(r[4:6], A.ρ_unit) - dot3D(A.q[4:6], A.ρ_unit) - ρ * A.v_α * dot3D(A.ρ_α, A.ρ_unit)
          - ρ * A.v_δ * dot3D(A.ρ_δ, A.ρ_unit)

    return ρ, v_ρ
end

@doc raw"""
    attr2bary(A::AdmissibleRegion{T}, a::Vector{U},
              params::NEOParameters{T}) where {T <: Real, U <: Number}

Convert attributable elements `a` to barycentric cartesian coordinates.
"""
function attr2bary(A::AdmissibleRegion{T}, a::Vector{U},
                   params::NEOParameters{T}) where {T <: Real, U <: Number}
    # Unfold
    α, δ, v_α, v_δ, ρ, v_ρ = a
    # Light-time correction to epoch
    t = datetime2days(A.date) - ρ/c_au_per_day
    # TO DO: `et::TaylorN` is too slow for `adam` due to
    # SatelliteToolboxTransformations overloads in src/observations/topocentric.jl
    et = datetime2et(A.date) - cte(cte(ρ))/c_au_per_sec
    # Line of sight vectors
    ρ_unit, ρ_α, ρ_δ = topounitpdv(α, δ)
    # Heliocentric position of the observer
    q = params.eph_ea(t) + kmsec2auday(obsposvelECI(A.observatory, et)) - params.eph_su(t)
    # Barycentric position
    r = q[1:3] + ρ * ρ_unit + params.eph_su(t)[1:3]
    # Barycentric velocity
    v = q[4:6] + v_ρ * ρ_unit + ρ * v_α * ρ_α + ρ * v_δ * ρ_δ
        + params.eph_su(t)[4:6]
    # Barycentric state vector
    return vcat(r, v)
end