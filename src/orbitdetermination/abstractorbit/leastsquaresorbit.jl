"""
    AbstractLeastSquaresOrbit{D, T, U} <: AbstractOrbit{D, T, U}

Supertype for the least squares orbits interface.

Every least squares orbit has:
- a dynamical function of type `D`,
- a vector of optical astrometry of type `O <: AbstractOpticalVector{T}`,
- a vector of optical tracklets of type `TrackletVector{T}`,
- a backward and a forward integration, both of type `DensePropagation2{T, U}`,
- a vector of optical residuals of type `Vector{OpticalResidual{T, U}}`,
- a least squares fit of type `LeastSquaresFit{T}`,
- a jacobian representing the transformation from the space of residuals to
    barycentric coordinates, of type `Matrix{T}`.
"""
abstract type AbstractLeastSquaresOrbit{D, T, U} <: AbstractOrbit{D, T, U} end

# AbstractOrbit interface
hasradar(x::AbstractLeastSquaresOrbit) = !isnothing(x.radar)

covariance(x::AbstractLeastSquaresOrbit) = x.fit.Γ

variances(x::AbstractLeastSquaresOrbit) = diag(x.jacobian * covariance(x) * x.jacobian')

# Initialize a vector of optical residuals consistent with od and orbit
function init_optical_residuals(
        ::Type{U}, od::AbstractODProblem{D, T},
        orbit::AbstractLeastSquaresOrbit
    ) where {D, T <: Real, U <: Number}
    # Optical astrometry
    optical1, optical2 = optical(od), optical(orbit)
    # Weights and debiasing factors
    w8s, bias = od.weights.w8s, od.debias.bias
    # Initialize vector of optical residuals
    res = Vector{OpticalResidual{T, U}}(undef, length(optical1))
    for i in eachindex(optical1)
        ra, dec = zero(U), zero(U)
        wra, wdec = w8s[i]
        dra, ddec = bias[i]
        j = findfirst(==(optical1[i]), optical2)
        outlier = isnothing(j) ? false : isoutlier(orbit.ores[j])
        res[i] = OpticalResidual{T, U}(ra, dec, wra, wdec, dra, ddec, outlier)
    end

    return res
end

"""
    LeastSquaresOrbit{D, T, U,
                      O <: AbstractOpticalVector{T},
                      R <: Union{Nothing, AbstractRadarVector{T}},
                      RR <: Union{Nothing, Vector{RadarResidual{T}}}
                      } <: AbstractLeastSquaresOrbit{D, T, U}

An asteroid least squares orbit.

# Fields

- `dynamics::D`: dynamical model.
- `optical::O`: vector of optical astrometry.
- `tracklets::TrackletVector{T}`: vector of optical tracklets.
- `radar::R`: vector of radar astrometry.
- `bwd/fwd::DensePropagation2{T, U}`: backward (forward) integration.
- `ores::Vector{OpticalResidual{T, U}}`: vector of optical residuals.
- `rres::RR`: vector of radar residuals.
- `fit::LeastSquaresFit{T}`: least squares fit.
- `jacobian::Matrix{T}`: space of residuals to barycentric coordinates jacobian.
- `qs::Matrix{T}`: history of initial conditions.
- `Qs::Vector{T}`: history of the target function.
"""
@auto_hash_equals struct LeastSquaresOrbit{D, T, U,
                                           O <: AbstractOpticalVector{T},
                                           R <: Union{Nothing, AbstractRadarVector{T}},
                                           RR <: Union{Nothing, Vector{RadarResidual{T}}}
                                           } <: AbstractLeastSquaresOrbit{D, T, U}
    dynamics::D
    optical::O
    tracklets::TrackletVector{T}
    radar::R
    bwd::DensePropagation2{T, U}
    fwd::DensePropagation2{T, U}
    ores::Vector{OpticalResidual{T, U}}
    rres::RR
    fit::LeastSquaresFit{T}
    jacobian::Matrix{T}
    qs::Matrix{T}
    Qs::Vector{T}
    # Inner constructor
    function LeastSquaresOrbit{D, T, U, O, R, RR}(
            dynamics::D, optical::O, tracklets::TrackletVector{T}, radar::R,
            bwd::TaylorInterpolant{T, U, 2}, fwd::TaylorInterpolant{T, U, 2},
            ores::Vector{OpticalResidual{T, U}}, rres::RR, fit::LeastSquaresFit{T},
            jacobian::Matrix{T}, qs::Matrix{T}, Qs::Vector{T}
        ) where {
            D, T <: Real, U <: Number, O <: AbstractOpticalVector{T},
            R <: Union{Nothing, AbstractRadarVector{T}},
            RR <: Union{Nothing, Vector{RadarResidual{T}}}
        }
        @assert bwd.t0 == fwd.t0 "Backward and forward integration initial \
            times must match"
        @assert length(optical) == length(ores) "Number of observations must \
            match number of residuals"
        _bwd_ = TaylorInterpolant(bwd.t0, bwd.t, collect(bwd.x))
        _fwd_ = TaylorInterpolant(fwd.t0, fwd.t, collect(fwd.x))
        return new{D, T, U, O, R, RR}(dynamics, optical, tracklets, radar, _bwd_,
                                      _fwd_, ores, rres, fit, jacobian, qs, Qs)
    end
end

# Outer constructor
function LeastSquaresOrbit(
        dynamics::D, optical::O, tracklets::TrackletVector{T}, radar::R,
        bwd::TaylorInterpolant{T, U, 2}, fwd::TaylorInterpolant{T, U, 2},
        ores::Vector{OpticalResidual{T, U}}, rres::RR, fit::LeastSquaresFit{T},
        jacobian::Matrix{T}, qs::Matrix{T}, Qs::Vector{T}
    ) where {
        D, T <: Real, U <: Number, O <: AbstractOpticalVector{T},
        R <: Union{Nothing, AbstractRadarVector{T}},
        RR <: Union{Nothing, Vector{RadarResidual{T}}}
    }
    return LeastSquaresOrbit{D, T, U, O, R, RR}(dynamics, optical, tracklets, radar, bwd,
                                                fwd, ores, rres, fit, jacobian, qs, Qs)
end

# Print method for LeastSquaresOrbit
show(io::IO, x::LeastSquaresOrbit) = print(io, "Least squares orbit with ", nobs(x),
                                           " residuals")

# Definition of zero LeastSquaresOrbit
function zero(::Type{LeastSquaresOrbit{D, T, U, O, R, RR}}) where {D, T, U, O, R, RR}
    dynamics = D.instance
    optical = O()
    tracklets = TrackletVector{T}()
    radar = R == Nothing ? nothing : R()
    bwd = zero(DensePropagation2{T, U})
    fwd = zero(DensePropagation2{T, U})
    ores = Vector{OpticalResidual{T, U}}(undef, 0)
    rres = RR == Nothing ? nothing : RR()
    fit = zero(LeastSquaresFit{T})
    jacobian = Matrix{T}(undef, 0, 0)
    qs = Matrix{T}(undef, 0, 0)
    Qs = Vector{T}(undef, 0)
    return LeastSquaresOrbit{D, T, U, O, R, RR}(dynamics, optical, tracklets, radar, bwd,
                                                fwd, ores, rres, fit, jacobian, qs, Qs)
end

iszero(x::LeastSquaresOrbit{D, T, U, O, R, RR}) where {D, T, U, O, R, RR} =
    x == zero(LeastSquaresOrbit{D, T, U, O, R, RR})

# Evaluate integrations and residuals in fit deltas
function evalfit(orbit::LeastSquaresOrbit{D, T, TaylorN{T}, O, R, RR}) where {D, T, O, R, RR}
    # Fit δs
    δs = orbit.fit.x
    # Evaluate integrations
    new_bwd_x = map(x -> Taylor1(x.coeffs(δs)), orbit.bwd.x)
    new_bwd = TaylorInterpolant(orbit.bwd.t0, orbit.bwd.t, new_bwd_x)
    new_fwd_x = map(x -> Taylor1(x.coeffs(δs)), orbit.fwd.x)
    new_fwd = TaylorInterpolant(orbit.fwd.t0, orbit.fwd.t, new_fwd_x)
    # Evaluate residuals
    new_ores = orbit.ores(δs)
    new_rres = hasradar(orbit) ? orbit.rres(δs) : nothing

    return LeastSquaresOrbit{D, T, T, O, R, RR}(orbit.dynamics, orbit.optical,
        orbit.tracklets, orbit.radar, new_bwd, new_fwd, new_ores, new_rres,
        orbit.fit, orbit.jacobian, orbit.qs, orbit.Qs)
end