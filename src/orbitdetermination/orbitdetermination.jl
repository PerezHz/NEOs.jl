include("osculating.jl")
include("curvature.jl")
include("odproblem.jl")
include("propres.jl")
include("leastsquares/methods.jl")
include("leastsquares/outlierrejection.jl")
include("leastsquares/fit.jl")
include("abstractorbit/abstractorbit.jl")
include("abstractorbit/preliminaryorbit.jl")
include("abstractorbit/leastsquaresorbit.jl")
include("preliminary/admissibleregion.jl")
include("preliminary/mmov.jl")
include("preliminary/gauss.jl")
include("jtls.jl")

# Default naive initial conditions for iod
function iodinitcond(A::AdmissibleRegion)
    v_ρ = sum(A.v_ρ_domain) / 2
    return [
        (A.ρ_domain[1], v_ρ, :log),
        (10^(sum(log10, A.ρ_domain) / 2), v_ρ, :log),
        (sum(A.ρ_domain) / 2, v_ρ, :log),
        (A.ρ_domain[2], v_ρ, :log),
    ]
end

"""
    issinglearc(::AbstractOpticalVector, arc::Day)

Check whether a (sorted) vector of optical astrometry is a single observational
arc, i.e. no two consecutive observations are more than `arc` days apart.
"""
issinglearc(x::AbstractOpticalVector, arc::Day = Day(30)) = all(diff(date.(x)) .< arc)

"""
    initialorbitdetermination(od, params; kwargs...)

Given an initial orbit determination problem `od`, return a `LeastSquaresOrbit`.
For a list of parameters, see [`Parameters`](@ref).

# Keyword arguments

- `initcond::I`: naive initial conditions function; takes as input an
    `AdmissibleRegion{T}` and outputs a `Vector{Tuple{T, T, Symbol}}`,
    where each element has the form `(ρ, v_ρ, scale)` (default: `iodinitcond`).

!!! warning
    This function may change the (global) `TaylorSeries` variables.

!!! reference
    See:
    - https://doi.org/10.1007/s10569-025-10246-2
"""
function initialorbitdetermination(od::OpticalODProblem{D, T, O}, params::Parameters{T};
                                   initcond::I = iodinitcond) where {D, I, T <: Real,
                                   O <: AbstractOpticalVector{T}}
    # Allocate memory for orbit
    orbit = zero(LeastSquaresOrbit{D, T, T, O, Nothing, Nothing})
    # Unpack
    @unpack tsaorder, gaussorder, jtlsorder, significance = params
    @unpack optical = od
    # Cannot handle observatories without coordinates
    all(x -> hascoord(observatory(x)), optical) || return orbit
    # Cannot handle zero observations or multiple arcs
    # (isempty(optical) || !issinglearc(optical)) && return orbit
    # Set jet transport variables
    set_od_order(params)
    # Gauss method
    _orbit_ = gaussiod(od, params)
    # Update orbit
    orbit = updateorbit(orbit, _orbit_, optical)
    # Termination condition
    critical_value(orbit) < significance && return orbit
    # Too short arc
    _orbit_ = tsaiod(od, params; initcond)
    # Update orbit
    orbit = updateorbit(orbit, _orbit_, optical)

    return orbit
end

"""
    orbitdetermination(od, orbit, params)

Given an orbit determination problem `od`, refine an a-priori `orbit`
via the Jet Transport Least Squares method. For a list of parameters,
see [`Parameters`](@ref).

!!! warning
    This function may change the (global) `TaylorSeries` variables.

!!! reference
    See:
    - https://doi.org/10.1007/s10569-025-10246-2
"""
function orbitdetermination(od::OpticalODProblem{D, T, O}, orbit::AbstractOrbit,
                            params::Parameters{T}) where {D, T <: Real,
                            O <: AbstractOpticalVector{T}}
    # Unpack parameters
    @unpack significance, verbose = params
    @unpack tracklets, optical = od
    # Set jet transport variables
    set_od_order(params)
    # Jet Transport Least Squares
    orbit1 = jtls(od, orbit, params, true)
    # Termination condition
    if (nobs(orbit1) == nobs(od) && critical_value(orbit1) < significance)
        N2 = length(orbit1.Qs)
        verbose && println(
            "* Jet Transport Least Squares converged in $N2 iterations to: \n\n",
            summary(orbit1)
        )
        return orbit1
    end
    # Refine via minimization over the MOV
    j = closest_tracklet(epoch(orbit), tracklets)
    for scale in (:log, :linear)
        porbit = mmov(od, orbit, j, params; scale)
        iszero(porbit) && continue
        # Jet Transport Least Squares
        orbit2 = jtls(od, porbit, params, true)
        # Update orbit
        orbit1 = updateorbit(orbit1, orbit2, optical)
        # Termination condition
        if (nobs(orbit1) == nobs(od) && critical_value(orbit1) < significance)
            N1, N2 = length(porbit.Qs), length(orbit1.Qs)
            verbose && println(
                "* Refinement of LeastSquaresOrbit via MMOV converged in \
                $N1 iterations to:\n\n",
                summary(porbit), "\n",
                "* Jet Transport Least Squares converged in $N2 iterations to: \n\n",
                summary(orbit1)
            )
            return orbit1
        end
    end
    # Unsuccessful orbit determination
    verbose && @warn("Orbit determination did not converge within \
        the given parameters or could not fit all the astrometry")

    return orbit1
end

function orbitdetermination(od::MixedODProblem{D, T, O, R}, orbit::AbstractOrbit,
                            params::Parameters{T}) where {D, T <: Real,
                            O <: AbstractOpticalVector{T}, R <: AbstractRadarVector{T}}
    # Unpack parameters
    @unpack significance, verbose = params
    @unpack tracklets, optical = od
    # Set jet transport variables
    set_od_order(params)
    # Jet Transport Least Squares
    orbit1 = jtls(od, orbit, params, true)
    # Termination condition
    if (nobs(orbit1) == nobs(od) && critical_value(orbit1) < significance)
        N2 = length(orbit1.Qs)
        verbose && println(
            "* Jet Transport Least Squares converged in $N2 iterations to: \n\n",
            summary(orbit1)
        )
        return orbit1
    end
    # Unsuccessful orbit determination
    verbose && @warn("Orbit determination did not converge within \
        the given parameters or could not fit all the astrometry")

    return orbit1
end