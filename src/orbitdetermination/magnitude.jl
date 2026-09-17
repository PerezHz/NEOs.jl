"""
    phase_angle(d_BS, d_BO, d_OS)

Return the phase angle given the distances between
the (B)ody, the (S)un and the (O)bserver.
"""
phase_angle(d_BS::Number, d_BO::Number, d_OS::Number) =
    acos((d_BS^2 + d_BO^2 - d_OS^2) / (2 * d_BS * d_BO))

"""
    elongation_angle(d_BS, d_BO, d_OS)

Return the elongation angle given the distances between
the (B)ody, the (S)un and the (O)bserver.
"""
elongation_angle(d_BS::Number, d_BO::Number, d_OS::Number) =
    acos((d_BO^2 + d_OS^2 - d_BS^2) / (2 * d_BO * d_OS))

"""
    phase_integral(α; kwargs...)

Return the phase integral for a given phase angle `α`
according to the Bowell et al. (1989) model.

See also [`phase_angle`](@ref).

# Keyword arguments

- `slope::Number`: slope parameter (default: `0.15`).

!!! reference
    See
    - https://ui.adsabs.harvard.edu/abs/1989aste.conf..524B/abstract
"""
function phase_integral(α::Number; slope::Number = 0.15)
    sin_α, tan_α = sin(α), tan(α/2)
    Φ1L = exp(-PHASE_INTEGRAL_A1 * tan_α^PHASE_INTEGRAL_B1)
    Φ2L = exp(-PHASE_INTEGRAL_A2 * tan_α^PHASE_INTEGRAL_B2)
    Φ1S = 1 - (PHASE_INTEGRAL_C1 * sin_α) / (0.119 + 1.341 * sin_α - 0.754 * sin_α^2)
    Φ2S = 1 - (PHASE_INTEGRAL_C2 * sin_α) / (0.119 + 1.341 * sin_α - 0.754 * sin_α^2)
    W = exp(-90.56 * tan_α^2)
    Φ1 = W * Φ1S + (1 - W) * Φ1L
    Φ2 = W * Φ2S + (1 - W) * Φ2L
    return (1 - slope) * Φ1 + slope * Φ2
end

"""
    magnitudedifference(d_BS, d_BO, d_OS; kwargs...)

Return the difference between the absolute and apparent magnitudes
given the distances between the (B)ody, the (S)un and the (O)bserver.

See also [`apparentmagnitude`](@ref) and [`absolutemagnitude`](@ref).

# Keyword arguments

- `slope::Number`: slope parameter (default: `0.15`).
"""
function magnitudedifference(d_BS::Number, d_BO::Number, d_OS::Number;
                             slope::Number = 0.15)
    α = phase_angle(d_BS, d_BO, d_OS)
    Φ = phase_integral(α; slope)
    return -5 * log10(d_BS * d_BO) + 2.5 * log10(Φ)
end

"""
    apparentmagnitude(H, d_BS, d_BO, d_OS; kwargs...)

Return the apparent magnitude given the absolute magnitude `H` and
the distances between the (B)ody, the (S)un and the (O)bserver.

See also [`absolutemagnitude`](@ref) and [`magnitudedifference`](@ref).

# Keyword arguments

- `slope::Number`: slope parameter (default: `0.15`).
"""
apparentmagnitude(H::Number, d_BS::Number, d_BO::Number, d_OS::Number;
    slope::Number = 0.15) = H - magnitudedifference(d_BS, d_BO, d_OS; slope)

"""
    absolutemagnitude(h, d_BS, d_BO, d_OS; kwargs...)

Return the absolute magnitude given the apparent magnitude `h` and
the distances between the (B)ody, the (S)un and the (O)bserver.

See also [`apparentmagnitude`](@ref) and [`magnitudedifference`](@ref).

# Keyword arguments

- `slope::Number`: slope parameter (default: `0.15`).
"""
absolutemagnitude(h::Number, d_BS::Number, d_BO::Number, d_OS::Number;
    slope::Number = 0.15) = h + magnitudedifference(d_BS, d_BO, d_OS; slope)