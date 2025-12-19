"""
    AttributablenElements{T, U} <: AbstractOsculatingElements{T, U}

A set of six attributable elements at a given epoch.

# Fields

- `gm::T`: gravitational parameter of the central body [au³/day²].
- `epoch::T`: reference epoch [MJD TDB].
- `frame::Symbol`: reference plane, either `:equatorial` or `:ecliptic`.
- `elements::SVector{6, U}`: set of six attributable elements.
- `covariance::SMatrix{6, 6, U, 36}`: covariance matrix.

# Extended help

The components of `elements` are:
- `α::U`: right ascension [deg].
- `δ::U`: declination [deg].
- `v_α::U`: right ascension velocity [deg/day].
- `v_δ::U`: declination velocity [deg/day].
- `ρ::U`: slant range [au].
- `v_ρ::U`: radial velocity [au/day].
"""
@auto_hash_equals struct AttributableElements{T, U} <: AbstractOsculatingElements{T, U}
    gm::T
    epoch::T
    frame::Symbol
    elements::SVector{6, U}
    covariance::SMatrix{6, 6, U, 36}
end

elementstype(::AttributableElements) = :Attributable
elementsnames(::AttributableElements) = ["α", "δ", "v_α", "v_δ", "ρ", "v_ρ"]
elementsunits(::AttributableElements) = ["deg", "deg", "deg/day", "deg/day", "au", "au/day"]

ra(x::AttributableElements) = x.elements[1]
dec(x::AttributableElements) = x.elements[2]
vra(x::AttributableElements) = x.elements[3]
vdec(x::AttributableElements) = x.elements[4]
range(x::AttributableElements) = x.elements[5]
rangerate(x::AttributableElements) = x.elements[6]

function eccentricity(x::AttributableElements)
    # Declination [rad]
    δ = deg2rad(dec(x))
    # Right ascension and declination velocities [rad/day]
    v_α, v_δ = deg2rad(vra(x)), deg2rad(vdec(x))
    # Range [au] and range rate [au/day]
    ρ, v_ρ = range(x), rangerate(x)
    # Squared proper motion
    η2 = v_α^2 * cos(δ)^2 + v_δ^2
    # Specific energy
    ε = 0.5 * (v_ρ^2 + ρ^2 * η2) - gm(x) / ρ
    # Squared specific angular momentum
    h2 = ρ^4 * η2
    # Eccentricity
    e = sqrt(1 + 2*ε*h2 / gm(x)^2)

    return e
end

"""
    cartesian2attributable(x; kwargs...)

Convert a cartesian state vector `x` [au, au/day] to
attributable elements [deg, au, day].

# Keyword arguments

- `frame::Symbol`: reference plane, either `:equatorial` (default) or `:ecliptic`.
"""
function cartesian2attributable(X::AbstractVector{U};
                                frame::Symbol = :equatorial) where {U <: Number}
    # If necessary, rotate state vector from equatorial to ecliptic plane
    if frame == :ecliptic
        X = equatorial2ecliptic(X)
    end
    # Cartesian coordinates [au, au/day]
    x, y, z, v_x, v_y, v_z = X
    # Range [au]
    ρ = sqrt(x^2 + y^2 + z^2)
    # Right ascension and declination [rad]
    α = mod2pi(atan(y, x))
    δ = asin(z / ρ)
    # Sines and cosines
    sin_α, cos_α = sincos(α)
    sin_δ, cos_δ = sincos(δ)
    sin_α_sin_δ = sin_α * sin_δ
    sin_α_cos_δ = sin_α * cos_δ
    cos_α_sin_δ = cos_α * sin_δ
    cos_α_cos_δ = cos_α * cos_δ
    # Radial velocity [au/day]
    v_ρ = v_x * cos_α_cos_δ + v_y * sin_α_cos_δ + v_z * sin_δ
    # Angular Velocities [rad/day]
    v_α = (-v_x * sin_α + v_y * cos_α) / (ρ * cos_δ)
    v_δ = (-v_x * cos_α_sin_δ - v_y * sin_α_sin_δ + v_z * cos_δ) / ρ
    # Vector of attributable elements
    α, δ, v_α, v_δ = rad2deg(α), rad2deg(δ), rad2deg(v_α), rad2deg(v_δ)
    elements = SVector{6, U}(α, δ, v_α, v_δ, ρ, v_ρ)

    return elements
end

"""
    attributable2cartesian(x)

Convert a set of attributable elements `x` [deg, au, day] to
a cartesian state vector [au, au/day].
"""
function attributable2cartesian(x::AbstractVector{U}) where {U <: Number}
    # Right ascension and declination [rad]
    α, δ = deg2rad(x[1]), deg2rad(x[2])
    # Right ascension and declination velocities [rad/day]
    v_α, v_δ = deg2rad(x[3]), deg2rad(x[4])
    # Range [au] and range rate [au/day]
    ρ, v_ρ = x[5], x[6]
    # Sine and cosines
    sin_α, cos_α = sincos(α)
    sin_δ, cos_δ = sincos(δ)
    sin_α_sin_δ = sin_α * sin_δ
    sin_α_cos_δ = sin_α * cos_δ
    cos_α_sin_δ = cos_α * sin_δ
    cos_α_cos_δ = cos_α * cos_δ
    # Positions [au]
    x = ρ * cos_α_cos_δ
    y = ρ * sin_α_cos_δ
    z = ρ * sin_δ
    # Velocities [au/day]
    ρ_v_α, ρ_v_δ = ρ * v_α, ρ * v_δ
    v_x = v_ρ * cos_α_cos_δ - ρ_v_α * sin_α_cos_δ - ρ_v_δ * cos_α_sin_δ
    v_y = v_ρ * sin_α_cos_δ + ρ_v_α * cos_α_cos_δ - ρ_v_δ * sin_α_sin_δ
    v_z = v_ρ * sin_δ + ρ_v_δ * cos_δ
    # Vector of cartesian coordinates
    elements = SVector{6, U}(x, y, z, v_x, v_y, v_z)

    return elements
end