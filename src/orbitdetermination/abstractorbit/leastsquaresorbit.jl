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
function init_optical_residuals(::Type{U}, od::ODProblem,
                                orbit::AbstractLeastSquaresOrbit) where {U <: Number}
    # Scalar type
    T = scalartype(od)
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

function init_radar_residuals(::Type{U}, od::MixedODProblem,
                              orbit::AbstractLeastSquaresOrbit) where {U <: Number}
    # Orbit does not have radar observations
    !hasradar(orbit) && return init_radar_residuals(U, od)
    # Scalar type
    T = scalartype(od)
    # Radar astrometry
    radar1, radar2 = radar(od), radar(orbit)
    # Initialize vector of radar residuals
    res = Vector{RadarResidual{T, U}}(undef, length(radar1))
    for i in eachindex(radar1)
        residual = zero(U)
        weight = 1 / rms(radar1[i])
        bias = debias(radar1[i])
        j = findfirst(==(radar1[i]), radar2)
        outlier = isnothing(j) ? false : isoutlier(orbit.rres[j])
        res[i] = RadarResidual{T, U}(residual, weight, bias, outlier)
    end

    return res
end

"""
    LeastSquaresOrbit{D, T, U,
                      O <: AbstractOpticalVector{T},
                      R <: Union{Nothing, AbstractRadarVector{T}},
                      RR <: Union{Nothing, Vector{RadarResidual{T, U}}}
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
                                           RR <: Union{Nothing, Vector{RadarResidual{T, U}}}
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
            RR <: Union{Nothing, Vector{RadarResidual{T, U}}}
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

# Abbreviations
const OpticalLeastSquaresOrbit{D, T, U, O} = LeastSquaresOrbit{D, T, U, O, Nothing, Nothing}
const MixedLeastSquaresOrbit{D, T, U, O, R} = LeastSquaresOrbit{D, T, U, O, R, Vector{RadarResidual{T, U}}}

# Outer constructor
function LeastSquaresOrbit(
        dynamics::D, optical::O, tracklets::TrackletVector{T}, radar::R,
        bwd::TaylorInterpolant{T, U, 2}, fwd::TaylorInterpolant{T, U, 2},
        ores::Vector{OpticalResidual{T, U}}, rres::RR, fit::LeastSquaresFit{T},
        jacobian::Matrix{T}, qs::Matrix{T}, Qs::Vector{T}
    ) where {
        D, T <: Real, U <: Number, O <: AbstractOpticalVector{T},
        R <: Union{Nothing, AbstractRadarVector{T}},
        RR <: Union{Nothing, Vector{RadarResidual{T, U}}}
    }
    return LeastSquaresOrbit{D, T, U, O, R, RR}(dynamics, optical, tracklets, radar, bwd,
                                                fwd, ores, rres, fit, jacobian, qs, Qs)
end

function LeastSquaresOrbit(od::OpticalODProblem{D, T, O}, q00::Vector{T}, jd0::T,
                           params::Parameters{T}) where {D, T <: Real, O}
    # Unpack
    @unpack dynamics, optical, tracklets, radar = od
    # Number of degrees of freedom
    Npar = dof(Val(od.dynamics))
    # Jet transport initial condition
    set_od_order(T, 2, Npar)
    if length(q00) < Npar
        q00 = vcat(q00, zero(T), zero(T))
    end
    scalings = fill(1e-8, 6)
    if Npar == 8
        scalings = vcat(scalings, 1e-14, 1e-15)
    end
    dq = scalings .* get_variables(T, 2)
    q0 = q00 + dq
    # Propagation and residuals
    bwd, fwd, res = propres(od, q0, jd0, params)
    # Least squares fit
    x0 = zeros(T, Npar)
    Q = nms(res)
    C = (notoutobs(res) / 2) * TS.hessian(Q, x0)
    Γ = inv(C)
    fit = LeastSquaresFit{T}(true, x0, Γ, Newton{T})
    # Residuals space to barycentric coordinates jacobian
    jacobian = Matrix(TS.jacobian(dq, fit.x))
    # History of initial conditions and target function
    qs = reshape(q00, Npar, 1)
    Qs = [nrms(res, fit)]

    return evalfit(LeastSquaresOrbit{D, T, TaylorN{T}, O, Nothing, Nothing}(dynamics,
        optical, tracklets, nothing, bwd, fwd, res, nothing, fit, jacobian, qs, Qs))
end

function LeastSquaresOrbit(od::MixedODProblem{D, T, O, R}, q00::Vector{T}, jd0::T,
                           params::Parameters{T}) where {D, T <: Real, O, R}
    # Unpack
    @unpack dynamics, optical, tracklets, radar = od
    # Number of degrees of freedom
    Npar = dof(Val(od.dynamics))
    # Jet transport initial condition
    set_od_order(T, 2, Npar)
    if length(q00) < Npar
        q00 = vcat(q00, zero(T), zero(T))
    end
    scalings = fill(1e-8, 6)
    if Npar == 8
        scalings = vcat(scalings, 1e-14, 1e-15)
    end
    dq = scalings .* get_variables(T, 2)
    q0 = q00 + dq
    # Propagation and residuals
    bwd, fwd, res = propres(od, q0, jd0, params)
    # Least squares fit
    x0 = zeros(T, Npar)
    Q = nms(res)
    C = (notoutobs(res) / 2) * TS.hessian(Q, x0)
    Γ = inv(C)
    fit = LeastSquaresFit{T}(true, x0, Γ, Newton{T})
    # Residuals space to barycentric coordinates jacobian
    jacobian = Matrix(TS.jacobian(dq, fit.x))
    # History of initial conditions and target function
    qs = reshape(q00, Npar, 1)
    Qs = [nrms(res[1], fit) + nrms(res[2], fit)]

    return evalfit(LeastSquaresOrbit{D, T, TaylorN{T}, O, R, Vector{RadarResidual{T, TaylorN{T}}}}(
        dynamics, optical, tracklets, radar, bwd, fwd, res[1], res[2], fit, jacobian, qs, Qs))
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
function evalfit(orbit::OpticalLeastSquaresOrbit{D, T, TaylorN{T}, O}) where {D, T, O}
    # Fit δs
    δs = orbit.fit.x
    # Evaluate integrations
    new_bwd_x = map(x -> Taylor1(x.coeffs(δs)), orbit.bwd.x)
    new_bwd = TaylorInterpolant(orbit.bwd.t0, orbit.bwd.t, new_bwd_x)
    new_fwd_x = map(x -> Taylor1(x.coeffs(δs)), orbit.fwd.x)
    new_fwd = TaylorInterpolant(orbit.fwd.t0, orbit.fwd.t, new_fwd_x)
    # Evaluate residuals
    new_ores = orbit.ores(δs)

    return LeastSquaresOrbit{D, T, T, O, Nothing, Nothing}(orbit.dynamics,
        orbit.optical, orbit.tracklets, nothing, new_bwd, new_fwd, new_ores, nothing,
        orbit.fit, orbit.jacobian, orbit.qs, orbit.Qs)
end

function evalfit(orbit::MixedLeastSquaresOrbit{D, T, TaylorN{T}, O, R}) where {D, T, O, R}
    # Fit δs
    δs = orbit.fit.x
    # Evaluate integrations
    new_bwd_x = map(x -> Taylor1(x.coeffs(δs)), orbit.bwd.x)
    new_bwd = TaylorInterpolant(orbit.bwd.t0, orbit.bwd.t, new_bwd_x)
    new_fwd_x = map(x -> Taylor1(x.coeffs(δs)), orbit.fwd.x)
    new_fwd = TaylorInterpolant(orbit.fwd.t0, orbit.fwd.t, new_fwd_x)
    # Evaluate residuals
    new_ores = orbit.ores(δs)
    new_rres = orbit.rres(δs)

    return LeastSquaresOrbit{D, T, T, O, R, Vector{RadarResidual{T, T}}}(orbit.dynamics,
        orbit.optical, orbit.tracklets, orbit.radar, new_bwd, new_fwd, new_ores, new_rres,
        orbit.fit, orbit.jacobian, orbit.qs, orbit.Qs)
end