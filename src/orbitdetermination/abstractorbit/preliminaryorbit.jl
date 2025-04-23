@doc raw"""
    AbstractPreliminaryOrbit{T, U} <: AbstractOrbit{T, U}

Supertype for the preliminary orbits API.
"""
abstract type AbstractPreliminaryOrbit{T, U} <: AbstractOrbit{T, U} end

init_residuals(::Type{U}, od::ODProblem, ::AbstractPreliminaryOrbit) where {U <: Number} =
    init_residuals(U, od)

variances(orbit::AbstractPreliminaryOrbit) = diag(orbit.J * orbit.Γ * orbit.J')

@doc raw"""
    GaussOrbit{T, U} <: AbstractPreliminaryOrbit{T, U}

A preliminary orbit computed by Gauss method.

## Fields

- `tracklets::Vector{Tracklet{T}}`: vector of tracklets.
- `bwd/fwd::TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}}`:
    backward (forward) integration.
- `res::Vector{OpticalResidual{T, U}}`: vector of optical residuals.
- `Γ::Matrix{T}`: covariance matrix.
- `J::Matrix{T}`: residuals space to barycentric coordinates jacobian.
- `τ_1/τ_3::T`: time between first (third) and second observations.
- `ρ_vec::Matrix{U}`: topocentric line-of-sight unit vectors.
- `R_vec::Matrix{T}`: observer's heliocentric positions.
- `D_0::U`: D scalar.
- `D::Matrix{U}`: D matrix.
- `a/b/c::U`: lagrange polynomial coefficients.
- `r_2::U`: lagrange polynomial root.
- `r_vec::Matrix{U}`: heliocentric state vectors of the object.
- `ρ::Vector{U}`: topocentric slant ranges.
"""
@auto_hash_equals struct GaussOrbit{T, U} <: AbstractPreliminaryOrbit{T, U}
    tracklets::Vector{Tracklet{T}}
    bwd::TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}}
    fwd::TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}}
    res::Vector{OpticalResidual{T, U}}
    Γ::Matrix{T}
    J::Matrix{T}
    τ_1::T
    τ_3::T
    ρ_vec::Matrix{U}
    R_vec::Matrix{T}
    D_0::U
    D::Matrix{U}
    a::U
    b::U
    c::U
    r_2::U
    r_vec::Matrix{U}
    ρ::Vector{U}
    # Inner constructor
    function GaussOrbit{T, U}(tracklets::Vector{Tracklet{T}},
        bwd::TaylorInterpolant{T, U, 2}, fwd::TaylorInterpolant{T, U, 2},
        res::Vector{OpticalResidual{T, U}}, Γ::Matrix{T}, J::Matrix{T},
        τ_1::T, τ_3::T, ρ_vec::Matrix{U}, R_vec::Matrix{T}, D_0::U,
        D::Matrix{U}, a::U, b::U, c::U, r_2::U, r_vec::Matrix{U},
        ρ::Vector{U}) where {T <: Real, U <: Number}
        @assert bwd.t0 == fwd.t0 "Backward and forward integration initial \
            times must match"
        @assert nobs(tracklets) == length(res) "Number of observations must \
            match number of residuals"
        _bwd_ = TaylorInterpolant(bwd.t0, bwd.t, collect(bwd.x))
        _fwd_ = TaylorInterpolant(fwd.t0, fwd.t, collect(fwd.x))
        return new{T, U}(tracklets, _bwd_, _fwd_, res, Γ, J, τ_1, τ_3, ρ_vec,
            R_vec, D_0, D, a, b, c, r_2, r_vec, ρ)
    end
end

# Outer constructor
GaussOrbit(tracklets::Vector{Tracklet{T}}, bwd::TaylorInterpolant{T, U, 2},
    fwd::TaylorInterpolant{T, U, 2}, res::Vector{OpticalResidual{T, U}},
    Γ::Matrix{T}, J::Matrix{T}, τ_1::T, τ_3::T, ρ_vec::Matrix{U},
    R_vec::Matrix{T}, D_0::U, D::Matrix{U}, a::U, b::U, c::U, r_2::U,
    r_vec::Matrix{U}, ρ::Vector{U}) where {T <: Real, U <: Number} =
    GaussOrbit{T, U}(tracklets, bwd, fwd, res, Γ, J, τ_1, τ_3, ρ_vec,
    R_vec, D_0, D, a, b, c, r_2, r_vec, ρ)

# Print method for GaussOrbit
show(io::IO, orbit::GaussOrbit) = print(io, "Gauss preliminary orbit with ",
    length(orbit.res), " residuals")

# Definition of zero GaussOrbit
function zero(::Type{GaussOrbit{T, U}}) where {T <: Real, U <: Number}
    tracklets = Vector{Tracklet{T}}(undef, 0)
    bwd = zero(TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}})
    fwd = zero(TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}})
    res = Vector{OpticalResidual{T, U}}(undef, 0)
    Γ = Matrix{T}(undef, 0, 0)
    J = Matrix{T}(undef, 0, 0)
    τ_1, τ_3 = zero(T), zero(T)
    ρ_vec = Matrix{U}(undef, 0, 0)
    R_vec = Matrix{T}(undef, 0, 0)
    D_0 = zero(U)
    D = Matrix{U}(undef, 0, 0)
    a, b, c = zero(U), zero(U), zero(U)
    r_2 = zero(U)
    r_vec = Matrix{U}(undef, 0, 0)
    ρ = Vector{U}(undef, 0)

    return GaussOrbit{T, U}(tracklets, bwd, fwd, res, Γ, J, τ_1, τ_3, ρ_vec,
        R_vec, D_0, D, a, b, c, r_2, r_vec, ρ)
end

iszero(orbit::GaussOrbit{T, U}) where {T <: Real, U <: Number} =
    orbit == zero(GaussOrbit{T, U})

# Evaluate TaylorN's in deltas
function evaldeltas(orbit::GaussOrbit{T, TaylorN{T}}, δs::Vector{T} =
    zeros(T, get_numvars())) where {T <: Real}
    # Evaluate integrations
    new_bwd_x = map(x -> Taylor1(x.coeffs(δs)), orbit.bwd.x)
    new_bwd = TaylorInterpolant(orbit.bwd.t0, orbit.bwd.t, new_bwd_x)
    new_fwd_x = map(x -> Taylor1(x.coeffs(δs)), orbit.fwd.x)
    new_fwd = TaylorInterpolant(orbit.fwd.t0, orbit.fwd.t, new_fwd_x)
    # Evaluate residuals
    new_res = orbit.res(δs)
    # Evaluate Gauss method objects
    new_ρ_vec = orbit.ρ_vec(δs)
    new_D_0 = orbit.D_0(δs)
    new_D = orbit.D(δs)
    new_a, new_b, new_c = orbit.a(δs), orbit.b(δs), orbit.c(δs)
    new_r_2 = orbit.r_2(δs)
    new_r_vec = orbit.r_vec(δs)
    new_ρ = orbit.ρ(δs)

    return GaussOrbit{T, T}(orbit.tracklets, new_bwd, new_fwd, new_res, orbit.Γ,
        orbit.J, orbit.τ_1, orbit.τ_3, new_ρ_vec, orbit.R_vec, new_D_0, new_D,
        new_a, new_b, new_c, new_r_2, new_r_vec, new_ρ)
end

@doc raw"""
    MMOVOrbit{T, U} <: AbstractPreliminaryOrbit{T, U}

A preliminary orbit computed by the minimization over the MOV method.

## Fields

- `tracklets::Vector{Tracklet{T}}`: vector of tracklets.
- `bwd/fwd::TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}}`:
    backward (forward) integration.
- `res::Vector{OpticalResidual{T, U}}`: vector of optical residuals.
- `Γ::Matrix{T}`: covariance matrix.
- `J::Matrix{T}`: residuals space to barycentric coordinates jacobian.
- `aes::Matrix{T}`: history of attributable elements.
- `Qs::Vector{T}`: history of the target function.
"""
@auto_hash_equals struct MMOVOrbit{T, U} <: AbstractPreliminaryOrbit{T, U}
    tracklets::Vector{Tracklet{T}}
    bwd::TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}}
    fwd::TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}}
    res::Vector{OpticalResidual{T, U}}
    Γ::Matrix{T}
    J::Matrix{T}
    aes::Matrix{T}
    Qs::Vector{T}
    # Inner constructor
    function MMOVOrbit{T, U}(
            tracklets::Vector{Tracklet{T}}, bwd::TaylorInterpolant{T, U, 2},
            fwd::TaylorInterpolant{T, U, 2}, res::Vector{OpticalResidual{T, U}},
            Γ::Matrix{T}, J::Matrix{T}, aes::Matrix{T}, Qs::Vector{T}
        ) where {T <: Real, U <: Number}
        @assert bwd.t0 == fwd.t0 "Backward and forward integration initial \
            times must match"
        @assert nobs(tracklets) == length(res) "Number of observations must \
            match number of residuals"
        # @assert size(aes, 1) == 6
        # @assert size(aes, 2) == length(Qs)
        _bwd_ = TaylorInterpolant(bwd.t0, bwd.t, collect(bwd.x))
        _fwd_ = TaylorInterpolant(fwd.t0, fwd.t, collect(fwd.x))
        return new{T, U}(tracklets, _bwd_, _fwd_, res, Γ, J, aes, Qs)
    end
end

# Outer constructor
MMOVOrbit(tracklets::Vector{Tracklet{T}}, bwd::TaylorInterpolant{T, U, 2},
    fwd::TaylorInterpolant{T, U, 2}, res::Vector{OpticalResidual{T, U}},
    Γ::Matrix{T}, J::Matrix{T}, aes::Matrix{T}, Qs::Vector{T}) where {T <: Real,
    U <: Number} = MMOVOrbit{T, U}(tracklets, bwd, fwd, res, Γ, J, aes, Qs)

# Print method for MMOVOrbit
show(io::IO, orbit::MMOVOrbit) = print(io, "Minimization over the MOV preliminary orbit ",
    "with ", length(orbit.res), " residuals")

# Definition of zero MMOVOrbit
function zero(::Type{MMOVOrbit{T, U}}) where {T <: Real, U <: Number}
    tracklets = Vector{Tracklet{T}}(undef, 0)
    bwd = zero(TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}})
    fwd = zero(TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}})
    res = Vector{OpticalResidual{T, U}}(undef, 0)
    Γ = Matrix{T}(undef, 0, 0)
    J = Matrix{T}(undef, 0, 0)
    aes = Matrix{T}(undef, 0, 0)
    Qs = Vector{T}(undef, 0)

    return MMOVOrbit{T, U}(tracklets, bwd, fwd, res, Γ, J, aes, Qs)
end

iszero(orbit::MMOVOrbit{T, U}) where {T <: Real, U <: Number} =
    orbit == zero(MMOVOrbit{T, U})

# Evaluate integrations and residuals in deltas
function evaldeltas(orbit::MMOVOrbit{T, TaylorN{T}}, δs::Vector{T} =
    zeros(T, get_numvars())) where {T <: Real}
    # Evaluate integrations
    new_bwd_x = map(x -> Taylor1(x.coeffs(δs)), orbit.bwd.x)
    new_bwd = TaylorInterpolant(orbit.bwd.t0, orbit.bwd.t, new_bwd_x)
    new_fwd_x = map(x -> Taylor1(x.coeffs(δs)), orbit.fwd.x)
    new_fwd = TaylorInterpolant(orbit.fwd.t0, orbit.fwd.t, new_fwd_x)
    # Evaluate residuals
    new_res = orbit.res(δs)

    return MMOVOrbit{T, T}(orbit.tracklets, new_bwd, new_fwd, new_res,
        orbit.Γ, orbit.J, orbit.aes, orbit.Qs)
end