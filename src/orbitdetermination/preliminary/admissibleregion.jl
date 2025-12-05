"""
    AdmissibleRegion{T <: Real}

Subset of topocentric range × range-rate space defined by the following constraints:
- heliocentric energy ≤ `k_gauss^2/(2a_max)`,
- absolute magnitude ≤ `H_max`,
- geocentric energy ≥ `0`.

# Fields

- `date::DateTime`: time of observation.
- `ra::T`: right ascension.
- `dec::T`: declination.
- `vra::T`: right ascension velocity.
- `vdec::T`: declination velocity.
- `H_max::T`: maximum absolute magnitude.
- `a_max::T`: maximum semimajor axis.
- `ρ_unit/ρ_α/ρ_δ::Vector{T}`: topocentric unit vector and its partials.
- `q::Vector{T}`: heliocentric position of observer.
- `coeffs::Vector{T}`: polynomial coefficients.
- `ρ_domain::Vector{T}`: range domain.
- `v_ρ_domain::Vector{T}`: range-rate domain.
- `Fs::Matrix{T}`: boundary points.
- `observatory::ObservatoryMPC{T}`: observing station.

!!! reference
    See Chapter 8 of:
    - https://doi.org/10.1017/CBO9781139175371
    or
    - https://doi.org/10.1007/s10569-004-6593-5
"""
@auto_hash_equals struct AdmissibleRegion{T <: Real}
    date::DateTime
    ra::T
    dec::T
    vra::T
    vdec::T
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
zero(::Type{AdmissibleRegion{T}}) where {T <: Real} = AdmissibleRegion{T}(
    DateTime(2000), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T),
    Vector{T}(undef, 0), Vector{T}(undef, 0), Vector{T}(undef, 0),
    Vector{T}(undef, 0), Vector{T}(undef, 0), Vector{T}(undef, 0),
    Vector{T}(undef, 0), Matrix{T}(undef, 0, 0), unknownobs(T)
)

iszero(x::AdmissibleRegion{T}) where {T <: Real} = x == zero(AdmissibleRegion{T})

# Print method for AdmissibleRegion
function show(io::IO, A::AdmissibleRegion{T}) where {T <: Real}
    v = string(
        @sprintf("%.5f", rad2deg(A.ra)), ", ",
        @sprintf("%.5f", rad2deg(A.dec)), ", ",
        @sprintf("%.5f", rad2deg(A.vra)), ", ",
        @sprintf("%.5f", rad2deg(A.vdec)), "",
    )
    print(io, "AE: [", v, "]", " t: ", A.date, " obs: ", A.observatory.name)
end

"""
    AdmissibleRegion(::OpticalTracklet, ::Parameters)

Return the admissible region associated to an optical tracklet. For a list of
parameters, see the `Minimization over the MOV` section of [`Parameters`](@ref).
"""
function AdmissibleRegion(tracklet::OpticalTracklet{T},
                          params::Parameters{T}) where {T <: Real}
    # Unpack
    @unpack observatory, date, ra, dec, vra, vdec, mag = tracklet
    @unpack eph_ea, eph_su, H_max, a_max = params
    # Topocentric unit vector and partials
    ρ, ρ_α, ρ_δ = topounitpdv(ra, dec)
    # Time of observation [days (et seconds) since J2000]
    t_days, t_et = dtutc2days(date), dtutc2et(date)
    # Heliocentric position of the observer
    q = eph_ea(t_days) + kmsec2auday(obsposvelECI(observatory, t_et)) - eph_su(t_days)
    # Admissible region coefficients
    coeffs = arcoeffs(ra, dec, vra, vdec, ρ, ρ_α, ρ_δ, q)
    # Maximum range (heliocentric energy constraint)
    ρ_max = _helmaxrange(coeffs, a_max)
    iszero(ρ_max) && return zero(AdmissibleRegion{T})
    # Minimum range
    if isnan(mag)
        # Earth's sphere of influence radius / Earth's physical radius
        ρ_min = R_SI < ρ_max ? R_SI : R_EA
    else
        # Tiny object boundary
        ρ_min = 10^((mag - H_max)/5)
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

    return AdmissibleRegion{T}(date, ra, dec, vra, vdec, H_max, a_max,
        ρ, ρ_α, ρ_δ, q, coeffs, ρ_domain, v_ρ_domain, Fs, observatory)
end

# Return the polynomial coefficients for an [`AdmissibleRegion`](@ref).
# See equation (8.8) of https://doi.org/10.1017/CBO9781139175371
function arcoeffs(α::T, δ::T, v_α::T, v_δ::T, ρ::Vector{T}, ρ_α::Vector{T},
                  ρ_δ::Vector{T}, q::AbstractVector{T}) where {T <: Number}
    coeffs = Vector{T}(undef, 6)
    coeffs[1] = dot(q[1:3], q[1:3])
    coeffs[2] = 2 *  dot(q[4:6], ρ)
    coeffs[3] = v_α^2 * cos(δ)^2 + v_δ^2  # proper motion squared
    coeffs[4] = 2 * v_α * dot(q[4:6], ρ_α) + 2 * v_δ * dot(q[4:6], ρ_δ)
    coeffs[5] = dot(q[4:6], q[4:6])
    coeffs[6] = 2*dot(q[1:3], ρ)
    return coeffs
end

# Return the topocentric line-of-sight unit vector and its
# partial derivatives with respect to `α` and `δ`.
# See between equations (8.5) and (8.6) of https://doi.org/10.1017/CBO9781139175371
function topounitpdv(α::T, δ::T) where {T <: Number}
    sin_α, cos_α = sincos(α)
    sin_δ, cos_δ = sincos(δ)

    ρ = [cos_δ * cos_α, cos_δ * sin_α, sin_δ]
    ρ_α = [-sin_α * cos_δ, cos_α * cos_δ, zero(α)]
    ρ_δ = [-cos_α * sin_δ, -sin_α * sin_δ, cos_δ]

    return ρ, ρ_α, ρ_δ
end

# W function of an [`AdmissibleRegion`](@ref).
# See equation (8.9) of https://doi.org/10.1017/CBO9781139175371
arW(A::AdmissibleRegion, ρ::Number) = arW(A.coeffs, ρ)
arW(coeffs::AbstractVector, ρ::Number) = coeffs[3] * ρ^2 + coeffs[4] * ρ + coeffs[5]

# S function of an [`AdmissibleRegion`](@ref).
# See equation (8.9) of https://doi.org/10.1017/CBO9781139175371
arS(A::AdmissibleRegion, ρ::Number) = arS(A.coeffs, ρ)
arS(coeffs::AbstractVector, ρ::Number) = ρ^2 + coeffs[6] * ρ + coeffs[1]

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
function arenergydis(A::AdmissibleRegion, ρ::Number, boundary::Symbol = :outer)
    a, b, c = arenergycoeffs(A, ρ, boundary)
    return b^2 - 4*a*c
end

function _arhelenergydis(coeffs::Vector{T}, a_max::T, ρ::Number) where {T <: Real}
    a, b, c = _arhelenergycoeffs(coeffs, a_max, ρ)
    return b^2 - 4*a*c
end

function _argeoenergydis(coeffs::Vector{T}, ρ::Number) where {T <: Real}
    a, b, c = _argeoenergycoeffs(coeffs, ρ)
    return b^2 - 4*a*c
end

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

function _georangerates(coeffs::Vector{T}, ρ::Number) where {T}
    ρ0 = min(R_SI, G⁻¹0(coeffs))
    !(0 < ρ <= ρ0) && return Vector{T}(undef, 0)
    a, b, c = _argeoenergycoeffs(coeffs, ρ)
    d = b^2 - 4*a*c
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
    d = b^2 - 4*a*c
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

function _georangerate(coeffs::Vector{T}, ρ::T, m::Symbol) where {T <: Real}
    ρ0 = min(R_SI, G⁻¹0(coeffs))
    @assert 0 < ρ <= ρ0 "No solutions for geocentric energy outside 0 < ρ <= ρ0"
    a, b, c = _argeoenergycoeffs(coeffs, ρ)
    d = b^2 - 4*a*c
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

# Project `[ρ, v_ρ]` into `A`'s outer boundary.
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
for U in (:(AbstractVector{T}), :(NTuple{2, T}))
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

"""
    topo2bary(::AdmissibleRegion, ρ, v_ρ)

Convert topocentric range `ρ` and range-rate `v_ρ` to barycentric cartesian coordinates.
The admissible region fixes the line of sight.
"""
function topo2bary(A::AdmissibleRegion, ρ::Number, v_ρ::Number)
    # Barycentric position
    r = A.q[1:3] + ρ * A.ρ_unit + sseph(su, dtutc2days(A.date))[1:3]
    # Barycentric velocity
    v = A.q[4:6] + v_ρ * A.ρ_unit + ρ * A.vra * A.ρ_α + ρ * A.vdec * A.ρ_δ
        + sseph(su, dtutc2days(A.date))[4:6]
    # Barycentric state vector
    return vcat(r, v)
end

"""
    bary2topo(::AdmissibleRegion, q0)

Convert barycentric cartesian coordinates `q0` to topocentric range and range-rate.
The admissible region fixes the line of sight.
"""
function bary2topo(A::AdmissibleRegion, q0::Vector{<:Number})
    # Heliocentric state vector
    r = q0 - sseph(su, dtutc2days(A.date))
    # Topocentric range
    ρ = euclid3D(r - A.q)
    # Topocentric range rate
    v_ρ = dot3D(r[4:6], A.ρ_unit) - dot3D(A.q[4:6], A.ρ_unit) - ρ * A.vra * dot3D(A.ρ_α, A.ρ_unit)
          - ρ * A.vdec * dot3D(A.ρ_δ, A.ρ_unit)

    return ρ, v_ρ
end

"""
    attr2bary(::AdmissibleRegion, a, ::Parameters)

Convert attributable elements `a` to barycentric cartesian coordinates.
The admissible region fixes the reference epoch and the parameters provide
Sun and Earth's ephemerides.
"""
function attr2bary(A::AdmissibleRegion{T}, a::Vector{U},
                   params::Parameters{T}) where {T <: Real, U <: Number}
    # Unfold
    α, δ, v_α, v_δ, ρ, v_ρ = a
    # Admissible region reference epoch
    # Note: we concluded both t and et should not include the relativistic
    # correction -ρ/c for consistency (05/12/2025)
    t = dtutc2days(A.date) # - ρ / c_au_per_day
    # TO DO: `et::TaylorN` is too slow for `mmov` due to
    # SatelliteToolboxTransformations overloads in src/observations/topocentric.jl
    et = dtutc2et(A.date) # - cte(cte(ρ)) / c_au_per_sec
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