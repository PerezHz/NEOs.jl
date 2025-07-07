"""
    AbstractPreliminaryOrbit{D, T, U} <: AbstractOrbit{D, T, U}

Supertype for the preliminary orbits interface.

Every preliminary orbit has:
- a dynamical function of type `D`,
- a vector of optical astrometry of type `O <: AbstractOpticalVector{T}`,
- a vector of optical tracklets of type `TrackletVector{T}`,
- a backward and a forward integration, both of type `DensePropagation2{T, U}`,
- a vector of optical residuals of type `Vector{OpticalResidual{T, U}}`,
- a covariance matrix of type `Matrix{T}`,
- a jacobian representing the transformation from the space of residuals to
    barycentric coordinates, of type `Matrix{T}`.
"""
abstract type AbstractPreliminaryOrbit{D, T, U} <: AbstractOrbit{D, T, U} end

# AbstractOrbit interface
hasradar(x::AbstractPreliminaryOrbit) = false

covariance(x::AbstractPreliminaryOrbit) = x.covariance

variances(x::AbstractPreliminaryOrbit) = diag(x.jacobian * covariance(x) * x.jacobian')

init_optical_residuals(::Type{U}, od::ODProblem,
    ::AbstractPreliminaryOrbit) where {U <: Number} = init_optical_residuals(U, od)

init_radar_residuals(::Type{U}, od::ODProblem,
    ::AbstractPreliminaryOrbit) where {U <: Number} = init_radar_residuals(U, od)

"""
    GaussOrbit{D, T, U,
               O <: AbstractOpticalVector{T}
               } <: AbstractPreliminaryOrbit{D, T, U}

A preliminary orbit computed by Gauss method.

# Fields

- `dynamics::D`: dynamical model.
- `optical::O`: vector of optical astrometry.
- `tracklets::TrackletVector{T}`: vector of optical tracklets.
- `bwd/fwd::DensePropagation2{T, U}`: backward (forward) integration.
- `ores::Vector{OpticalResidual{T, U}}`: vector of optical residuals.
- `covariance::Matrix{T}`: covariance matrix.
- `jacobian::Matrix{T}`: space of residuals to barycentric coordinates jacobian.
- `τ_1/τ_3::T`: time between first (third) and second observations.
- `ρ_vec::Matrix{U}`: topocentric line-of-sight unit vectors.
- `R_vec::Matrix{T}`: observer's heliocentric positions.
- `D_0::U`: D scalar.
- `D_mat::Matrix{U}`: D matrix.
- `a/b/c::U`: lagrange polynomial coefficients.
- `r_2::U`: lagrange polynomial root.
- `r_vec::Matrix{U}`: heliocentric state vectors of the object.
- `ρ::Vector{U}`: topocentric slant ranges.
"""
@auto_hash_equals struct GaussOrbit{D, T, U,
                                    O <: AbstractOpticalVector{T}
                                    } <: AbstractPreliminaryOrbit{D, T, U}
    dynamics::D
    optical::O
    tracklets::TrackletVector{T}
    bwd::DensePropagation2{T, U}
    fwd::DensePropagation2{T, U}
    ores::Vector{OpticalResidual{T, U}}
    covariance::Matrix{T}
    jacobian::Matrix{T}
    τ_1::T
    τ_3::T
    ρ_vec::Matrix{U}
    R_vec::Matrix{T}
    D_0::U
    D_mat::Matrix{U}
    a::U
    b::U
    c::U
    r_2::U
    r_vec::Matrix{U}
    ρ::Vector{U}
    # Inner constructor
    function GaussOrbit{D, T, U, O}(
            dynamics::D, optical::O, tracklets::TrackletVector{T},
            bwd::TaylorInterpolant{T, U, 2}, fwd::TaylorInterpolant{T, U, 2},
            ores::Vector{OpticalResidual{T, U}}, covariance::Matrix{T},
            jacobian::Matrix{T}, τ_1::T, τ_3::T, ρ_vec::Matrix{U}, R_vec::Matrix{T},
            D_0::U, D_mat::Matrix{U}, a::U, b::U, c::U, r_2::U, r_vec::Matrix{U},
            ρ::Vector{U}
        ) where {D, T <: Real, U <: Number, O <: AbstractOpticalVector{T}}
        @assert bwd.t0 == fwd.t0 "Backward and forward integration initial \
            times must match"
        @assert length(optical) == length(ores) "Number of observations must \
            match number of residuals"
        _bwd_ = TaylorInterpolant(bwd.t0, bwd.t, collect(bwd.x))
        _fwd_ = TaylorInterpolant(fwd.t0, fwd.t, collect(fwd.x))
        return new{D, T, U, O}(dynamics, optical, tracklets, _bwd_, _fwd_, ores,
            covariance, jacobian, τ_1, τ_3, ρ_vec, R_vec, D_0, D_mat, a, b, c,
            r_2, r_vec, ρ)
    end
end

# Outer constructor
function GaussOrbit(
    dynamics::D, optical::O, tracklets::TrackletVector{T},
    bwd::TaylorInterpolant{T, U, 2}, fwd::TaylorInterpolant{T, U, 2},
    ores::Vector{OpticalResidual{T, U}}, covariance::Matrix{T},
    jacobian::Matrix{T}, τ_1::T, τ_3::T, ρ_vec::Matrix{U}, R_vec::Matrix{T},
    D_0::U, D_mat::Matrix{U}, a::U, b::U, c::U, r_2::U, r_vec::Matrix{U},
    ρ::Vector{U}
) where {D, T <: Real, U <: Number, O <: AbstractOpticalVector{T}}
    return GaussOrbit{D, T, U, O}(dynamics, optical, tracklets, bwd, fwd, ores,
        covariance, jacobian, τ_1, τ_3, ρ_vec, R_vec, D_0, D_mat, a, b, c, r_2,
        r_vec, ρ)
end

# Print method for GaussOrbit
show(io::IO, x::GaussOrbit) = print(io, "Gauss preliminary orbit with ", nobs(x),
                                    " residuals")

# Definition of zero GaussOrbit
function zero(::Type{GaussOrbit{D, T, U, O}}) where {D, T, U, O}
    dynamics = D.instance
    optical = O()
    tracklets = TrackletVector{T}()
    bwd = zero(DensePropagation2{T, U})
    fwd = zero(DensePropagation2{T, U})
    ores = Vector{OpticalResidual{T, U}}(undef, 0)
    covariance = Matrix{T}(undef, 0, 0)
    jacobian = Matrix{T}(undef, 0, 0)
    τ_1, τ_3 = zero(T), zero(T)
    ρ_vec = Matrix{U}(undef, 0, 0)
    R_vec = Matrix{T}(undef, 0, 0)
    D_0 = zero(U)
    D_mat = Matrix{U}(undef, 0, 0)
    a, b, c = zero(U), zero(U), zero(U)
    r_2 = zero(U)
    r_vec = Matrix{U}(undef, 0, 0)
    ρ = Vector{U}(undef, 0)

    return GaussOrbit{D, T, U, O}(dynamics, optical, tracklets, bwd, fwd, ores,
        covariance, jacobian, τ_1, τ_3, ρ_vec, R_vec, D_0, D_mat, a, b, c, r_2,
        r_vec, ρ)
end

iszero(x::GaussOrbit{D, T, U, O}) where {D, T, U, O} = x == zero(GaussOrbit{D, T, U, O})

# Check topocentric slant ranges are positive
isphysical(x::GaussOrbit) = iszero(x) ? false : all(x.ρ .> 0)

# Check heliocentric energy is negative
function isclosed(x::GaussOrbit)
    iszero(x) && return false
    # Heliocentric state vector [au, au/day]
    r = x.r_vec[:, 2]
    # Heliocentric energy per unit mass
    kinetic = 0.5 * (r[4]^2 + r[5]^2 + r[6]^2)
    potential = k_gauss^2 / sqrt(r[1]^2 + r[2]^2 + r[3]^2)
    E  = kinetic - potential

    return E <= 0
end

# Evaluate TaylorN's in deltas
function evaldeltas(orbit::GaussOrbit{D, T, TaylorN{T}, O},
                    δs::Vector{T} = zeros(T, get_numvars())) where {D, T, O}
    # Evaluate integrations
    new_bwd_x = map(x -> Taylor1(x.coeffs(δs)), orbit.bwd.x)
    new_bwd = TaylorInterpolant(orbit.bwd.t0, orbit.bwd.t, new_bwd_x)
    new_fwd_x = map(x -> Taylor1(x.coeffs(δs)), orbit.fwd.x)
    new_fwd = TaylorInterpolant(orbit.fwd.t0, orbit.fwd.t, new_fwd_x)
    # Evaluate residuals
    new_ores = orbit.ores(δs)
    # Evaluate Gauss method objects
    new_ρ_vec = orbit.ρ_vec(δs)
    new_D_0 = orbit.D_0(δs)
    new_D_mat = orbit.D_mat(δs)
    new_a, new_b, new_c = orbit.a(δs), orbit.b(δs), orbit.c(δs)
    new_r_2 = orbit.r_2(δs)
    new_r_vec = orbit.r_vec(δs)
    new_ρ = orbit.ρ(δs)

    return GaussOrbit{D, T, T, O}(orbit.dynamics, orbit.optical, orbit.tracklets,
        new_bwd, new_fwd, new_ores, orbit.covariance, orbit.jacobian, orbit.τ_1,
        orbit.τ_3, new_ρ_vec, orbit.R_vec, new_D_0, new_D_mat, new_a, new_b, new_c,
        new_r_2, new_r_vec, new_ρ)
end

"""
    MMOVOrbit{D, T, U,
              O <: AbstractOpticalVector{T}
              } <: AbstractPreliminaryOrbit{D, T, U}

A preliminary orbit computed by the minimization over the MOV method.

# Fields

- `dynamics::D`: dynamical model.
- `optical::O`: vector of optical astrometry.
- `tracklets::TrackletVector{T}`: vector of optical tracklets.
- `bwd/fwd::DensePropagation2{T, U}`: backward (forward) integration.
- `ores::Vector{OpticalResidual{T, U}}`: vector of optical residuals.
- `covariance::Matrix{T}`: covariance matrix.
- `jacobian::Matrix{T}`: space of residuals to barycentric coordinates jacobian.
- `aes::Matrix{T}`: history of attributable elements.
- `Qs::Vector{T}`: history of the target function.
"""
@auto_hash_equals struct MMOVOrbit{D, T, U,
                                   O <: AbstractOpticalVector{T}
                                   } <: AbstractPreliminaryOrbit{D, T, U}
    dynamics::D
    optical::O
    tracklets::TrackletVector{T}
    bwd::DensePropagation2{T, U}
    fwd::DensePropagation2{T, U}
    ores::Vector{OpticalResidual{T, U}}
    covariance::Matrix{T}
    jacobian::Matrix{T}
    aes::Matrix{T}
    Qs::Vector{T}
    # Inner constructor
    function MMOVOrbit{D, T, U, O}(
            dynamics::D, optical::O, tracklets::TrackletVector{T},
            bwd::TaylorInterpolant{T, U, 2}, fwd::TaylorInterpolant{T, U, 2},
            ores::Vector{OpticalResidual{T, U}}, covariance::Matrix{T},
            jacobian::Matrix{T}, aes::Matrix{T}, Qs::Vector{T}
        ) where {D, T <: Real, U <: Number, O <: AbstractOpticalVector{T}}
        @assert bwd.t0 == fwd.t0 "Backward and forward integration initial \
            times must match"
        @assert length(optical) == length(ores) "Number of observations must \
            match number of residuals"
        # @assert size(aes, 1) == 6
        # @assert size(aes, 2) == length(Qs)
        _bwd_ = TaylorInterpolant(bwd.t0, bwd.t, collect(bwd.x))
        _fwd_ = TaylorInterpolant(fwd.t0, fwd.t, collect(fwd.x))
        return new{D, T, U, O}(dynamics, optical, tracklets, _bwd_, _fwd_,
                               ores, covariance, jacobian, aes, Qs)
    end
end

# Outer constructor
function MMOVOrbit(
        dynamics::D, optical::O, tracklets::TrackletVector{T},
        bwd::TaylorInterpolant{T, U, 2}, fwd::TaylorInterpolant{T, U, 2},
        ores::Vector{OpticalResidual{T, U}}, covariance::Matrix{T},
        jacobian::Matrix{T}, aes::Matrix{T}, Qs::Vector{T}
    ) where {D, T <: Real, U <: Number, O <: AbstractOpticalVector{T}}
    return MMOVOrbit{D, T, U, O}(dynamics, optical, tracklets, bwd, fwd,
                                 ores, covariance, jacobian, aes, Qs)
end

# Print method for MMOVOrbit
show(io::IO, x::MMOVOrbit) = print(io, "Minimization over the MOV preliminary orbit ",
                                   "with ", nobs(x), " residuals")

# Definition of zero MMOVOrbit
function zero(::Type{MMOVOrbit{D, T, U, O}}) where {D, T, U, O}
    dynamics = D.instance
    optical = O()
    tracklets = TrackletVector{T}()
    bwd = zero(DensePropagation2{T, U})
    fwd = zero(DensePropagation2{T, U})
    ores = Vector{OpticalResidual{T, U}}(undef, 0)
    covariance = Matrix{T}(undef, 0, 0)
    jacobian = Matrix{T}(undef, 0, 0)
    aes = Matrix{T}(undef, 0, 0)
    Qs = Vector{T}(undef, 0)

    return MMOVOrbit{D, T, U, O}(dynamics, optical, tracklets, bwd, fwd,
                                 ores, covariance, jacobian, aes, Qs)
end

iszero(x::MMOVOrbit{D, T, U, O}) where {D, T, U, O} = x == zero(MMOVOrbit{D, T, U, O})

# Evaluate integrations and residuals in deltas
function evaldeltas(orbit::MMOVOrbit{D, T, TaylorN{T}, O},
                    δs::Vector{T} = zeros(T, get_numvars())) where {D, T, O}
    # Evaluate integrations
    new_bwd_x = map(x -> Taylor1(x.coeffs(δs)), orbit.bwd.x)
    new_bwd = TaylorInterpolant(orbit.bwd.t0, orbit.bwd.t, new_bwd_x)
    new_fwd_x = map(x -> Taylor1(x.coeffs(δs)), orbit.fwd.x)
    new_fwd = TaylorInterpolant(orbit.fwd.t0, orbit.fwd.t, new_fwd_x)
    # Evaluate residuals
    new_ores = orbit.ores(δs)

    return MMOVOrbit{D, T, T, O}(orbit.dynamics, orbit.optical, orbit.tracklets, new_bwd,
                                 new_fwd, new_ores, orbit.covariance, orbit.jacobian,
                                 orbit.aes, orbit.Qs)
end