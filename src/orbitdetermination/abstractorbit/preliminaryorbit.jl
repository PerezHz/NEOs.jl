@doc raw"""
    PreliminaryOrbit{T, U} <: AbstractOrbit{T, U}

An asteroid preliminary orbit.

## Fields

- `tracklets::Vector{Tracklet{T}}`: vector of tracklets.
- `bwd/fwd::TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}}`:
    backward (forward) integration.
- `res::Vector{OpticalResidual{T, U}}`: vector of optical residuals.
"""
@auto_hash_equals struct PreliminaryOrbit{T, U} <: AbstractOrbit{T, U}
    tracklets::Vector{Tracklet{T}}
    bwd::TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}}
    fwd::TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}}
    res::Vector{OpticalResidual{T, U}}
    # Inner constructor
    function PreliminaryOrbit{T, U}(tracklets::Vector{Tracklet{T}},
        bwd::TaylorInterpolant{T, U, 2}, fwd::TaylorInterpolant{T, U, 2},
        res::Vector{OpticalResidual{T, U}}) where {T <: Real, U <: Number}
        @assert bwd.t0 == fwd.t0 "Backward and forward integration initial \
            times must match"
        @assert nobs(tracklets) == length(res) "Number of observations must \
            match number of residuals"
        _bwd_ = TaylorInterpolant(bwd.t0, bwd.t, collect(bwd.x))
        _fwd_ = TaylorInterpolant(fwd.t0, fwd.t, collect(fwd.x))
        return new{T, U}(tracklets, _bwd_, _fwd_, res)
    end
end

# Outer constructor
PreliminaryOrbit(tracklets::Vector{Tracklet{T}}, bwd::TaylorInterpolant{T, U, 2},
    fwd::TaylorInterpolant{T, U, 2}, res::Vector{OpticalResidual{T, U}}) where
    {T <: Real, U <: Number} = PreliminaryOrbit{T, U}(tracklets, bwd, fwd, res)

# Print method for PreliminaryOrbit
show(io::IO, orbit::PreliminaryOrbit) = print(io, "Preliminary orbit with ",
    length(orbit.res), " residuals")

# Definition of zero PreliminaryOrbit
function zero(::Type{PreliminaryOrbit{T, U}}) where {T <: Real, U <: Number}
    tracklets = Vector{Tracklet{T}}(undef, 0)
    bwd = zero(TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}})
    fwd = zero(TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}})
    res = Vector{OpticalResidual{T, U}}(undef, 0)
    return PreliminaryOrbit{T, U}(tracklets, bwd, fwd, res)
end

iszero(orbit::PreliminaryOrbit{T, U}) where {T <: Real, U <: Number} =
    orbit == zero(PreliminaryOrbit{T, U})