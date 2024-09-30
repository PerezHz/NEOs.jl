@doc raw"""
    ODProblem{D, T <: Real, WT <: AbstractWeightingScheme{T},
        DT <: AbstractDebiasingScheme{T}}

An orbit determination problem.

## Fields

- `dynamics::D`: dynamical model.
- `radec::Vector{RadecMPC{T}}`: vector of optical astrometry.
- `tracklets::Vector{Tracklet{T}}`: vector of tracklets.
- `w8s::WT`: weighting scheme.
- `bias::DT`: debiasing scheme.
"""
mutable struct ODProblem{D, T <: Real, WT <: AbstractWeightingScheme{T},
    DT <: AbstractDebiasingScheme{T}}
    dynamics::D
    radec::Vector{RadecMPC{T}}
    tracklets::Vector{Tracklet{T}}
    w8s::WT
    bias::DT
    @doc raw"""
        ODProblem(dynamics::D, radec::Vector{RadecMPC{T}} [,
            w8s::WT [, bias::DT]]) where {D, T <: Real, WT, DT}

    Return an orbit determination problem.

    ## Arguments

    - `dynamics::D`: dynamical model.
    - `radec::Vector{RadecMPC{T}}`: vector of optical astrometry.
    - `w8s::WT`: weighting scheme (default: `Veres17`).
    - `bias::DT`: debiasing scheme (default: `Eggl20`).
    """
    function ODProblem(dynamics::D, radec::Vector{RadecMPC{T}},
        w8s::WT = Veres17, bias::DT = Eggl20) where {D, T <: Real, WT, DT}
        # Reduce tracklets by polynomial regression
        tracklets = reduce_tracklets(radec)
        # Build weighting scheme
        _w8s_ = w8s(radec)
        # Build debiasing scheme
        _bias_ = bias(radec)
        # Consistency check
        @assert length(radec) == sum(nobs, tracklets) == length(_w8s_.w8s) ==
            length(_bias_.bias)
        # Assemble orbit determination problem
        new{D, T, typeof(_w8s_), typeof(_bias_)}(dynamics, radec, tracklets,
            _w8s_, _bias_)
    end
end

const ODProblem{D, T} = ODProblem{D, T, WT, DT} where {D, T <: Real,
    WT <: AbstractWeightingScheme{T}, DT <: AbstractDebiasingScheme{T}}

# Print method for ODProblem
function show(io::IO, p::ODProblem)
    t = "    "
    print(io, "Orbit determination problem:\n")
    print(io, t, rpad("Dynamical model:", 21), p.dynamics, "\n")
    print(io, t, rpad("Astrometry:", 21), length(p.radec),
        " optical observations (", length(p.tracklets), " tracklets)\n")
    print(io, t, rpad("Weighting scheme:", 21), getid(p.w8s), "\n")
    print(io, t, rpad("Debiasing scheme:", 21), getid(p.bias), "\n")
end

# Override update!
function update!(p::ODProblem{D, T}, radec::Vector{RadecMPC{T}}) where {D, T <: Real}
    p.radec = radec
    p.tracklets = reduce_tracklets(radec)
    update!(p.w8s, radec)
    update!(p.bias, radec)
    return nothing
end