"""
    OpticalResidual{T, U} <: AbstractOpticalResidual{T, U}

An astrometric optical observed minus computed residual.

# Fields

- `ra/dec::U`: normalized right ascension (declination) residual [arcsec].
- `wra/wdec::T`: right ascension (declination) weight [arcsec⁻¹].
- `dra/ddec::T`: right ascension (declination) debiasing factor [arcsec].
- `corr::T` correlation between `ra` and `dec`.
- `outlier::Bool`: whether the residual is an outlier or not.
"""
@auto_hash_equals struct OpticalResidual{T, U} <: AbstractOpticalResidual{T, U}
    ra::U
    dec::U
    wra::T
    wdec::T
    dra::T
    ddec::T
    corr::T
    outlier::Bool
    # Inner constructor
    function OpticalResidual{T, U}(
        ra::U, dec::U, wra::T, wdec::T, dra::T,
        ddec::T, corr::T, outlier::Bool = false
    ) where {T <: Real, U <:  Number}
        return new{T, U}(ra, dec, wra, wdec, dra, ddec, corr, outlier)
    end
end

# AbstractAstrometryResidual interface
residual(x::OpticalResidual) = (ra(x), dec(x))
weight(x::OpticalResidual) = (wra(x), wdec(x))
debias(x::OpticalResidual) = (dra(x), ddec(x))

dof(::Type{OpticalResidual{T, U}}) where {T, U} = 2

function chi2(x::OpticalResidual)
    ξ_α, ξ_δ, ρ = ra(x), dec(x), corr(x)
    if isoutlier(x) || abs(ρ) ≥ 1
        return zero(ξ_α)
    else
        return (ξ_α^2 + ξ_δ^2 - 2ρ*ξ_α*ξ_δ) / (1 - ρ^2)
    end
end

ra(x::OpticalResidual) = x.ra
dec(x::OpticalResidual) = x.dec
wra(x::OpticalResidual) = x.wra
wdec(x::OpticalResidual) = x.wdec
dra(x::OpticalResidual) = x.dra
ddec(x::OpticalResidual) = x.ddec
corr(x::OpticalResidual) = x.corr

# Definition of zero OpticalResidual
zero(::Type{OpticalResidual{T, U}}) where {T, U} = OpticalResidual{T, U}(
        zero(U), zero(U), zero(T), zero(T), zero(T), zero(T), zero(T), false)
iszero(x::OpticalResidual{T, U}) where {T, U} = x == zero(OpticalResidual{T, U})

# Print method for OpticalResidual
function show(io::IO, x::OpticalResidual)
    outlier_flag = isoutlier(x) ? " (outlier)" : ""
    print(io, "α: ", @sprintf("%+.5f", cte(ra(x))), " δ: ",
        @sprintf("%+.5f", cte(dec(x))), outlier_flag)
end

# Evaluate methods
evaluate(y::OpticalResidual{T, TaylorN{T}}, x::Vector{T}) where {T <: Real} =
    OpticalResidual{T, T}(y.ra(x), y.dec(x), y.wra, y.wdec, y.dra, y.ddec, y.corr, y.outlier)

(y::OpticalResidual{T, TaylorN{T}})(x::Vector{T}) where {T <: Real} = evaluate(y, x)

function evaluate(y::AbstractVector{OpticalResidual{T, TaylorN{T}}},
                  x::Vector{T}) where {T <: Real}
    z = Vector{OpticalResidual{T, T}}(undef, length(y))
    for i in eachindex(z)
        z[i] = evaluate(y[i], x)
    end
    return z
end

(y::AbstractVector{OpticalResidual{T, TaylorN{T}}})(x::Vector{T}) where {T <: Real} =
    evaluate(y, x)

"""
    unfold(::AbstractVector{OpticalResidual})

Return three vectors by concatenating the non-outlier right ascension
and declination residuals, weights and debiasing factors.
"""
function unfold(y::AbstractVector{OpticalResidual{T, U}}) where {T <: Real, U <: Number}
    # Number of non outliers
    L = notout(y)
    # Vector of residuals, weights and debiasing factors
    z = Vector{U}(undef, 2L)
    w = Vector{T}(undef, 2L)
    d = Vector{T}(undef, 2L)
    # Global counter
    k = 1
    # Fill residuals, weights and debiasing factors
    for i in eachindex(y)
        isoutlier(y[i]) && continue
        # Right ascension
        z[k], w[k], d[k] = ra(y[i]), wra(y[i]), dra(y[i])
        # Declination
        z[k+L], w[k+L], d[k+L] = dec(y[i]), wdec(y[i]), ddec(y[i])
        # Update global counter
        k += 1
    end

    return z, w, d
end

function normalized_residuals(y::AbstractVector{OpticalResidual{T, U}}) where {T, U}
    # Number of non outliers
    L = notout(y)
    # Vector of normalized residuals
    z = Vector{U}(undef, 2L)
    # Global counter
    k = 1
    # Fill residuals
    for i in eachindex(y)
        isoutlier(y[i]) && continue
        # Right ascension and declination
        z[k], z[k+L] = ra(y[i]), dec(y[i])
        # Update global counter
        k += 1
    end

    return z
end

# Angle difference taking into account the discontinuity in [0, 2π) -> [0, 2π)
# x and y must be in arcsec
function anglediff(x::Number, y::Number)
    # Signed difference
    Δ = x - y
    # Absolute difference
    Δ_abs = abs(Δ)
    # Reflection
    if Δ_abs > 648_000 # Half circle in arcsec
        return -sign(cte(Δ)) * (1_296_000 - Δ_abs)
    else
        return Δ
    end
end

function init_optical_residuals(
        ::Type{U}, optical::AbstractOpticalVector{T}, w8s::AbstractVector{NTuple{2, T}},
        bias::AbstractVector{NTuple{2, T}}, corrs::AbstractVector{T},
        outliers::AbstractVector{Bool}
    ) where {T <: Real, U <: Number}
    # Check consistency between arrays
    @assert length(optical) == length(w8s) == length(bias) == length(outliers)
    # Initialize vector of residuals
    res = Vector{OpticalResidual{T, U}}(undef, length(optical))
    for i in eachindex(optical)
        ra, dec = zero(U), zero(U)
        wra, wdec = w8s[i]
        dra, ddec = bias[i]
        corr = corrs[i]
        outlier = -1 < corr < 1 ? outliers[i] : true
        res[i] = OpticalResidual{T, U}(ra, dec, wra, wdec, dra, ddec, corr, outlier)
    end

    return res
end

"""
    residuals(optical, w8s, bias [, corrs [, outliers]]; kwargs...) where {T <: Real}

Compute the observed minus computed residuals for a vector of optical astrometry.
Corrections due to Earth orientation, LOD and polar motion are computed by default.

See also [`OpticalResidual`](@ref) and [`compute_radec`](@ref).

# Arguments

- `optical::AbstractOpticalVector{T}`: optical astrometry.
- `w8s::AbstractVector{NTuple{2, T}}`: statistical weights.
- `bias::AbstractVector{NTuple{2, T}}`: debiasing corrections.
- `corrs::AbstractVector{T}`: correlations (default:
    `zeros(T, length(optical))`).
- `outliers::AbstractVector{Bool}`: outlier flags (default:
    `falses(length(optical))`).

# Keyword arguments

- `niter::Int`: number of light-time solution iterations (default: `5`).
- `xve`: Earth ephemeris (default: `earthposvel`).
- `xvs`: Sun ephemeris (default: `sunposvel`).
- `xva`: asteroid ephemeris.

All ephemeris must take [et seconds since J2000] and return [barycentric position in km
and velocity in km/sec].
"""
function residuals(optical::AbstractOpticalVector{T},
                   w8s::AbstractVector{NTuple{2, T}},
                   bias::AbstractVector{NTuple{2, T}},
                   corrs::AbstractVector{T} = zeros(T, length(optical)),
                   outliers::AbstractVector{Bool} = falses(length(optical));
                   xva::AstEph, kwargs...) where {AstEph, T <: Real}
    # UTC time of first optical observation
    utc1 = date(optical[1])
    # TDB seconds since J2000.0 for first optical observation
    et1 = dtutc2et(utc1)
    # Asteroid ephemeris at et1
    a1_et1 = evaleph(xva, et1)[1]
    # Type of asteroid ephemeris
    U = typeof(a1_et1)
    # Buffer
    buffer = [OpticalBuffer(a1_et1) for _ in eachindex(optical)]
    # Vector of residuals
    res = init_optical_residuals(U, optical, w8s, bias, corrs, outliers)
    residuals!(res, optical, buffer; xva, kwargs...)

    return res
end

function residuals!(res::Vector{OpticalResidual{T, U}},
                    optical::AbstractOpticalVector{T},
                    buffer::Vector{OpticalBuffer{U}};
                    kwargs...) where {T <: Real, U <: Number}

    @allow_boxed_captures tmap!(res, optical, buffer, weight.(res), debias.(res), corr.(res),
                                isoutlier.(res)) do x, buff, w8s, bias, rho, outlier
        # Observed ra/dec [arcsec]
        obsra, obsdec = rad2arcsec.(measure(x))
        # Computed ra/dec [arcsec]
        compra, compdec = compute_radec(x, buff; kwargs...)
        # Statistical weights [arcsec⁻²]
        wra, wdec = w8s
        # Debiasing factors [arcsec]
        dra, ddec = bias
        # Observed minus computed residual ra/dec
        # Note: ra is multiplied by a metric factor cos(dec) to match the format of
        # debiasing corrections
        return OpticalResidual{T, U}(
            wra * ( anglediff(obsra, compra) * cos(dec(x)) - dra ),
            wdec * ( obsdec - compdec - ddec ),
            wra,
            wdec,
            dra,
            ddec,
            rho,
            -1 < rho < 1 ? outlier : true
        )
    end

    return nothing
end