"""
    RadarResidual{T, U} <: AbstractRadarResidual{T, U}

An astrometric radar observed minus computed residual.

# Fields

- `residual::U`: normalized time delay [us] or Doppler shift [Hz] residual.
- `weight::T`: statistical weight [same units as `residual`⁻¹].
- `debias::T`: debiasing factor [same units as `residual`].
- `outlier::Bool`: whether the residual is an outlier or not.
"""
@auto_hash_equals struct RadarResidual{T, U} <: AbstractRadarResidual{T, U}
    residual::U
    weight::T
    debias::T
    outlier::Bool
    # Inner constructor
    function RadarResidual{T, U}(residual::U, weight::T, debias::T,
                                 outlier::Bool = false) where {T <: Real, U <:  Number}
        return new{T, U}(residual, weight, debias, outlier)
    end
end

# AbstractAstrometryResidual interface
residual(x::RadarResidual) = x.residual
weight(x::RadarResidual) = x.weight
debias(x::RadarResidual) = x.debias

dof(::Type{RadarResidual{T, U}}) where {T, U} = 1
chi2(x::RadarResidual) = !isoutlier(x) * residual(x)^2

# Definition of zero RadarResidual
zero(::Type{RadarResidual{T, U}}) where {T, U} = RadarResidual{T, U}(
        zero(U), zero(T), zero(T), false)
iszero(x::RadarResidual{T, U}) where {T, U} = x == zero(RadarResidual{T, U})

# Print method for RadarResidual
function show(io::IO, x::RadarResidual)
    outlier_flag = isoutlier(x) ? " (outlier)" : ""
    print(io, "residual: ", @sprintf("%+.5f", cte(residual(x))), outlier_flag)
end

# Evaluate methods
evaluate(y::RadarResidual{T, TaylorN{T}}, x::Vector{T}) where {T <: Real} =
    RadarResidual{T, T}(y.residual(x), y.weight, y.debias, y.outlier)

(y::RadarResidual{T, TaylorN{T}})(x::Vector{T}) where {T <: Real} = evaluate(y, x)

function evaluate(y::AbstractVector{RadarResidual{T, TaylorN{T}}},
                  x::Vector{T}) where {T <: Real}
    z = Vector{RadarResidual{T, T}}(undef, length(y))
    for i in eachindex(z)
        z[i] = evaluate(y[i], x)
    end
    return z
end

(y::AbstractVector{RadarResidual{T, TaylorN{T}}})(x::Vector{T}) where {T <: Real} =
    evaluate(y, x)

"""
    unfold(::AbstractVector{RadarResidual})

Return three vectors by concatenating the non-outlier time-delay and
Doppler shift residuals, weights and debiasing factors.
"""
function unfold(y::AbstractVector{RadarResidual{T, U}}) where {T <: Real, U <: Number}
    # Number of non outliers
    L = notout(y)
    # Vector of residuals, weights and debiasing factors
    z = Vector{U}(undef, L)
    w = Vector{T}(undef, L)
    d = Vector{T}(undef, L)
    # Global counter
    k = 1
    # Fill residuals, weights and debiasing factors
    for i in eachindex(y)
        isoutlier(y[i]) && continue
        # Residual
        z[k], w[k], d[k] = residual(y[i]), weight(y[i]), debias(y[i])
        # Update global counter
        k += 1
    end

    return z, w, d
end

function normalized_residuals(y::AbstractVector{RadarResidual{T, U}}) where {T, U}
    # Number of non outliers
    L = notout(y)
    # Vector of normalized residuals
    z = Vector{U}(undef, L)
    # Global counter
    k = 1
    # Fill residuals
    for i in eachindex(y)
        isoutlier(y[i]) && continue
        # Residual
        z[k] = residual(y[i])
        # Update global counter
        k += 1
    end

    return z
end

function init_radar_residuals(
        ::Type{U}, radar::AbstractRadarVector{T},
        outliers::AbstractVector{Bool}
    ) where {T <: Real, U <: Number}
    # Check consistency between arrays
    @assert length(radar) == length(outliers)
    # Initialize vector of residuals
    res = Vector{RadarResidual{T, U}}(undef, length(radar))
    for i in eachindex(radar)
        residual = zero(U)
        weight = 1 / rms(radar[i])
        bias = debias(radar[i])
        res[i] = RadarResidual{T, U}(residual, weight, bias, outliers[i])
    end

    return res
end

"""
    residuals(radar, [, outliers]; xva, kwargs...)  where {AstEph, T <: Real}

Compute the observed minus computed residuals for a vector of radar astrometry.
Corrections due to Earth orientation, LOD and polar motion are computed by default.

See also [`RadarResidual`](@ref) and [`radar_astrometry`](@ref).

# Arguments

- `radar::AbstractRadarVector{T}`: radar astrometry.
- `outliers::AbstractVector{Bool}`: outlier flags (default:
    `falses(length(radar))`).

# Keyword arguments

- `tc::Real`: time offset wrt echo reception time, to compute Doppler
    shifts by range differences [sec]. Offsets will be rounded to the
    nearest millisecond (default: `1.0`).
- `autodiff::Bool`: whether to compute Doppler shift via automatic diff
    of [`compute_delay`](@ref) or not (default: `true`).
- `tord::Int`: order of Taylor expansions (default: `10`).
- `niter::Int`: number of light-time solution iterations (default: `10`).
- `xve::EarthEph`: Earth ephemeris (default: `earthposvel`).
- `xvs::SunEph`: Sun ephemeris (default: `sunposvel`).
- `xva::AstEph`: asteroid ephemeris.

All ephemeris must take [et seconds since J2000] and return [barycentric
position in km and velocity in km/sec].
"""
function residuals(radar::AbstractRadarVector{T},
                   outliers::AbstractVector{Bool} = falses(length(radar));
                   xva::AstEph, kwargs...) where {AstEph, T <: Real}
    # UTC time of first radar observation
    utc1 = date(radar[1])
    # TDB seconds since J2000.0 for first radar observation
    et1 = dtutc2et(utc1)
    # Asteroid ephemeris at et1
    a1_et1 = xva(et1)[1]
    # Type of asteroid ephemeris
    U = typeof(a1_et1)
    # Vector of residuals
    res = init_radar_residuals(U, radar, outliers)
    residuals!(res, radar; xva, kwargs...)

    return res
end

function residuals!(res::Vector{RadarResidual{T, U}}, radar::AbstractRadarVector{T};
                    xva::AstEph, kwargs...) where {AstEph, T <: Real, U <: Number}

    @allow_boxed_captures tmap!(res, radar, weight.(res), debias.(res),
                                isoutlier.(res)) do x, w8, bias, outlier
        # Observed time-delay or Doppler shift
        observed = measure(x)
        # Computed time-delay and Doppler shift
        delay, doppler = radar_astrometry(x; xva, kwargs...)
        computed = isdelay(x) ? delay : doppler
        # Observed minus computed residual
        return RadarResidual{T, U}(
            w8 * ( observed - computed - bias ),
            w8,
            bias,
            outlier
        )
    end

    return nothing
end
