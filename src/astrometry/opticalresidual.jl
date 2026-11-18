"""
    OpticalResidual{T, U} <: AbstractOpticalResidual{T, U}

An astrometric optical observed minus computed residual.

# Fields

- `ra/dec::U`: debiased right ascension (times the cosine of the declination) and
    declination observed minus computed residuals [arcsec].
- `weight::SMatrix{2, 2, T, 4}`: weight matrix, i.e. the inverse of the covariance
    matrix of `(ra, dec)` [arcsec⁻²].
- `dra/ddec::T`: right ascension (declination) debiasing factor [arcsec].
- `outlier::Bool`: whether the residual is an outlier or not.

!!! note
    A non-outlier residual contributes `ξᵀ W ξ` to the chi-square, where
    `ξ = (ra, dec)` and `W = weight`.

See also [`normalized_residual`](@ref) and [`weightmatrix`](@ref).
"""
@auto_hash_equals struct OpticalResidual{T, U} <: AbstractOpticalResidual{T, U}
    ra::U
    dec::U
    weight::SMatrix{2, 2, T, 4}
    dra::T
    ddec::T
    outlier::Bool
    # Inner constructor
    function OpticalResidual{T, U}(
            ra::U, dec::U, weight::SMatrix{2, 2, T, 4}, dra::T, ddec::T,
            outlier::Bool = false
        ) where {T <: Real, U <: Number}
        return new{T, U}(ra, dec, weight, dra, ddec, outlier)
    end
end

# AbstractAstrometryResidual interface
residual(x::OpticalResidual) = (ra(x), dec(x))
weight(x::OpticalResidual) = x.weight
debias(x::OpticalResidual) = (dra(x), ddec(x))

dof(::Type{OpticalResidual{T, U}}) where {T, U} = 2

function chi2(x::OpticalResidual)
    ξ_α, ξ_δ, W = ra(x), dec(x), weight(x)
    isoutlier(x) && return zero(ξ_α)
    return W[1, 1] * ξ_α^2 + 2 * W[1, 2] * ξ_α * ξ_δ + W[2, 2] * ξ_δ^2
end
chi(x::OpticalResidual) = sqrt(chi2(x))
logchi(x::OpticalResidual) = log(chi(x))
log10chi(x::OpticalResidual) = log10(chi(x))

"""
    normalized_residual(x::OpticalResidual)

Return the normalized (whitened) right ascension and declination residuals
`z = U ξ`, where `ξ = (ra(x), dec(x))` and `W = Uᵀ U` is the Cholesky
factorization of the weight matrix `W = weight(x)`. Hence, the components of
`z` are uncorrelated with unit variance and `z₁² + z₂² == chi2(x)` for a
non-outlier residual.

!!! warning
    The weight matrix must be positive definite.
"""
function normalized_residual(x::OpticalResidual{T}) where {T <: Real}
    W = weight(x)
    u11 = sqrt(W[1, 1])
    u12 = W[1, 2] / u11
    u22 = sqrt(max(zero(T), W[2, 2] - u12^2))
    return u11 * ra(x) + u12 * dec(x), u22 * dec(x)
end

ra(x::OpticalResidual) = x.ra
dec(x::OpticalResidual) = x.dec
dra(x::OpticalResidual) = x.dra
ddec(x::OpticalResidual) = x.ddec
setoutlier(x::OpticalResidual{T, U}, y::Bool) where {T, U} = OpticalResidual{T, U}(
    ra(x), dec(x), weight(x), dra(x), ddec(x), y)

# Definition of zero OpticalResidual
zero(::Type{OpticalResidual{T, U}}) where {T, U} = OpticalResidual{T, U}(
    zero(U), zero(U), SMatrix{2, 2, T, 4}(zero(T), zero(T), zero(T), zero(T)),
    zero(T), zero(T), false)
iszero(x::OpticalResidual{T, U}) where {T, U} = x == zero(OpticalResidual{T, U})

# Print method for OpticalResidual
function show(io::IO, x::OpticalResidual)
    outlier_flag = isoutlier(x) ? " (outlier)" : ""
    print(io, "α: ", @sprintf("%+.5f", cte(ra(x))), " δ: ",
        @sprintf("%+.5f", cte(dec(x))), outlier_flag)
end

# Methods to print optical residuals in the MPEC format
function mpec_residual(x::AbstractOpticalAstrometry, y::OpticalResidual)
    sd = Dates.format(date(x), "yymmdd")
    so = observatorycode(x)
    α, δ = ra(y), dec(y)
    sα = if abs(α) < 100
        @sprintf("%5.1f", abs(α)) * if abs(α) < 0.05
            "  "
        else
            α ≥ 0 ? "+ " : "- "
        end
    else
        @sprintf("%5.2f", abs(α / 3_600)) * (α ≥ 0 ? "+ " : "- ")
    end
    sδ = if abs(δ) < 100
        @sprintf("%5.1f", abs(δ)) * if abs(δ) < 0.05
            "  "
        else
            δ ≥ 0 ? "+ " : "- "
        end
    else
        @sprintf("%5.2f", abs(δ / 3_600)) * (δ ≥ 0 ? "+ " : "- ")
    end
    s = string(sd, " ", so, sα, sδ)
    if isoutlier(y)
        c = collect(s)
        c[11], c[end] = '(', ')'
        s = join(c)
    end
    return s
end

print_mpec_residuals(x::AbstractOpticalVector, y::AbstractResidualVector) =
    print_mpec_residuals(stdout, x, y)

function print_mpec_residuals(io::IO, x::AbstractOpticalVector,
                              y::AbstractResidualVector)
    @assert length(x) == length(y) "Observations and residuals vectors must have \
        the same length"
    L = length(x)
    Nrows = ceil(Int, L / 3)
    v = Vector{String}(undef, 3Nrows)
    @. v[1:L] = mpec_residual(x, y)
    @. v[L+1:end] = ""
    m = reshape(v, Nrows, 3)
    write(
        io,
        "Residuals in seconds of arc\n",
        join(join.(eachrow(m), "   "), '\n'), '\n'
    )
    return nothing
end

# Evaluate methods
evaluate(y::OpticalResidual{T, TaylorN{T}}, x::Vector{T}) where {T <: Real} =
    OpticalResidual{T, T}(y.ra(x), y.dec(x), y.weight, y.dra, y.ddec, y.outlier)

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
    normalized_residuals(::AbstractVector{OpticalResidual})

Return a vector with the normalized (whitened) right ascension residuals of the
non-outliers, followed by the normalized declination residuals.

See also [`normalized_residual`](@ref).
"""
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
        z[k], z[k+L] = normalized_residual(y[i])
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

# Initialize vector of optical residuals
function init_optical_residuals(
        ::Type{U},
        wsm::AbstractWeightingScheme{T},
        dsm::AbstractDebiasingScheme{T},
        idxs::AbstractVector{Int} = 1:nobs(wsm)
    ) where {T <: Real, U <: Number}
    @assert nobs(wsm) == nobs(dsm) "Weighting and debiasing schemes must have the \
        same number of observations"
    Γs = covariances(wsm, idxs)
    bias = debias(dsm, idxs)
    outs = outliers(wsm, idxs)
    res = Vector{OpticalResidual{T, U}}(undef, length(idxs))
    for i in eachindex(res)
        res[i] = OpticalResidual{T, U}(zero(U), zero(U), weightmatrix(Γs[i]),
                                       bias[i]..., outs[i])
    end
    return res
end

"""
    residuals(optical, wsm, dsm; kwargs...) where {T <: Real}

Compute the observed minus computed residuals for a vector of optical astrometry.
Corrections due to Earth orientation, LOD and polar motion are computed by default.

See also [`OpticalResidual`](@ref) and [`compute_radec`](@ref).

# Arguments

- `optical::AbstractOpticalVector{T}`: optical astrometry.
- `wsm::AbstractWeightingScheme{T}`: weighting scheme.
- `dsm::AbstractDebiasingScheme{T}`: debiasing scheme.

# Keyword arguments

- `niter::Int`: number of light-time solution iterations (default: `5`).
- `xve`: Earth ephemeris (default: `earthposvel`).
- `xvs`: Sun ephemeris (default: `sunposvel`).
- `xva`: asteroid ephemeris.

All ephemeris must take [et seconds since J2000] and return [barycentric position in
km and velocity in km/sec].
"""
function residuals(
        optical::AbstractOpticalVector{T},
        wsm::AbstractWeightingScheme{T},
        dsm::AbstractDebiasingScheme{T};
        xva::AstEph, kwargs...
    ) where {AstEph, T <: Real}
    # UTC time of first optical observation
    utc1 = date(optical[1])
    # TDB seconds since J2000.0 for first optical observation
    et1 = dtutc2et(utc1)
    # Asteroid ephemeris at et1
    a1_et1 = evaleph(xva, et1)[1]
    # Type of asteroid ephemeris
    U = typeof(a1_et1)
    # Buffer
    buffer = [OpticalBuffer(a1_et1, x) for x in optical]
    # Vector of residuals
    res = init_optical_residuals(U, wsm, dsm)
    residuals!(res, optical, buffer, covariances(wsm); xva, kwargs...)

    return res
end

"""
    residuals!(res, optical, buffer, Γs; kwargs...)

Overwrite the vector of optical residuals `res` with the observed minus computed
residuals of `optical`, using the vector of buffers `buffer` and the vector of
`(α cos(δ), δ)` covariance matrices `Γs` [arcsec²] to assemble the weight matrices.
The debiasing factors and outlier flags are taken from `res`. For a list of keyword
arguments see [`residuals`](@ref).
"""
function residuals!(res::AbstractVector{OpticalResidual{T, U}},
                    optical::AbstractOpticalVector{T},
                    buffer::Vector{OpticalBuffer{T, U}},
                    Γs::AbstractVector{SMatrix{2, 2, T, 4}};
                    kwargs...) where {T <: Real, U <: Number}

    # Function barrier: the closure only captures the (never reassigned) arguments,
    # so none of its captured variables is boxed
    tforeach(eachindex(res, optical, buffer, Γs)) do i
        res[i] = optical_residual(res[i], optical[i], buffer[i], Γs[i]; kwargs...)
    end

    return nothing
end

# Observed minus computed residual of optical observation `x` with covariance matrix `Γ`;
# the debiasing factors and outlier flag are taken from the previous residual `r`, while
# the weight matrix is assembled from `Γ`. Observations with a non positive definite
# covariance matrix are flagged as outliers
function optical_residual(r::OpticalResidual{T, U}, x::AbstractOpticalAstrometry{T},
                          buffer::OpticalBuffer{T, U}, Γ::SMatrix{2, 2, T, 4};
                          kwargs...) where {T <: Real, U <: Number}
    # Observed ra/dec [arcsec]
    obsra, obsdec = rad2arcsec.(measure(x))
    # Computed ra/dec [arcsec]
    compra, compdec = compute_radec(x, buffer; kwargs...)
    # Debiasing factors [arcsec]
    dra, ddec = debias(r)
    # Weight matrix [arcsec⁻²]
    W = weightmatrix(Γ)
    # Observed minus computed residual ra/dec
    # Note: ra is multiplied by a metric factor cos(dec) to match the format of
    # debiasing corrections
    return OpticalResidual{T, U}(
        anglediff(obsra, compra) * cos(dec(x)) - dra,
        obsdec - compdec - ddec,
        W,
        dra,
        ddec,
        isoutlier(r) || !isposdefcov(Γ)
    )
end
