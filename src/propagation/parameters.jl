@doc raw"""
    NEOParameters([params::NEOParameters{T};] kwargs...) where {T <: Real}

Parameters for all orbit determination functions.

## Propagation Parameters

- `maxsteps::Int`: maximum number of steps for the integration (default: `500`).
- `μ_ast::Vector{T}`: vector of gravitational parameters
    (default: `μ_ast343_DE430[1:end]`).
- `order::Int`: order of Taylor expansions wrt time (default: 25).
- `abstol::T`: absolute tolerance used to compute propagation timestep
    (default: `1e-20`).
- `parse_eqs::Bool`: whether to use the specialized method of `jetcoeffs` or not
    (default: `true`).
- `bwdoffset/fwdoffset::T`: days to propagate beyond first (bwd) / last (fwd) observation
    (default: `0.5`).
- `coeffstol::T`: maximum size of the coefficients (default: `10.0`).

## Gauss Method Parameters

- `max_triplets::Int`: maximum number of triplets to check for a solution
    (default: `10`).
- `gaussorder::Int`: order of the jet transport perturbation (default: `5`).
- `adamhelp::Bool`: whether to refine Gauss preliminary orbit via ADAM
    (default: `false`).
- `gaussQmax::T`: nrms threshold (default: `5.0`).

## Too Short Arc Parameters

- `H_max::T`: maximum absolute magnitude (default: `34.5`).
- `a_max::T`: maximum semimajor axis (default: `100.0`).
- `adamiter::Int`: maximum number of iterations for `ADAM` optimizer (default: `200`).
- `adammode::Bool`: whether to perform ADAM iterations with all the observations
    (default: `false`).
- `adamQtol::T`: target function relative tolerance (default: `0.001`).
- `tsaorder::Int`: order of the jet transport perturbation (default: `6`).
- `tsaQmax::T`: nrms threshold (default: `1.5`).

## Jet Transport Least Squares Parameters

- `lsiter::Int`: maximum number of iterations for `leastsquares` (default: `5`).
- `jtlsiter::Int`: maximum number of iterations for `jtls` (default: `5`).
- `jtlsorder::Int`: order of the jet transport perturbation in `jtls` (default: `5`).

## Outlier Rejection Parameters

- `outrej::Bool`: whether to perform outlier rejection during least squares
    iterations (default: `false`).
- `χ2_rec::T`: recovery threshold (default: `7.0`).
- `χ2_rej::T`: rejection threshold (default: `8.0`).
- `fudge::T`: rejection fudge term coefficient (default: `400.0`).
- `max_per::T`: maximum allowed rejection percentage (default: `10.0`).
"""
@kwdef struct NEOParameters{T <: Real}
    # Propagation Parameters
    maxsteps::Int = 500
    μ_ast::Vector{T} = μ_ast343_DE430[1:end]
    order::Int = 25
    abstol::T = 1e-20
    parse_eqs::Bool = true
    bwdoffset::T = 0.5
    fwdoffset::T = 0.5
    coeffstol::T = 10.0
    # Sun (earth) ephemeris
    eph_su::TaylorInterpolant{T, T, 2, Vector{T}, Matrix{Taylor1{T}}} = _loadephsu()
    eph_ea::TaylorInterpolant{T, T, 2, Vector{T}, Matrix{Taylor1{T}}} = _loadephea()
    # Gauss' Method Parameters
    max_triplets::Int = 10
    gaussorder::Int = 5
    adamhelp::Bool = false
    gaussQmax::T = 5.0
    # Too Short Arc Parameters
    H_max::T = 34.5
    a_max::T = 100.0
    adamiter::Int = 200
    adammode::Bool = false
    adamQtol::T = 0.001
    tsaorder::Int = 6
    tsaQmax::T = 1.5
    # Jet Transport Least Squares Parameters
    lsiter::Int = 5
    jtlsiter::Int = 5
    jtlsorder::Int = 5
    # Outlier Rejection Parameters
    outrej::Bool = false
    χ2_rec::T = 7.0
    χ2_rej::T = 8.0
    fudge::T = 400.0
    max_per::T = 10.0
end

# Outer constructors
function NEOParameters(params::NEOParameters{T}; kwargs...) where {T <: Real}
    fields = fieldnames(NEOParameters{T})
    vals = Vector{Any}(undef, length(fields))
    for i in eachindex(vals)
        if fields[i] in keys(kwargs)
            vals[i] = kwargs[fields[i]]
        else
            vals[i] = getfield(params, i)
        end
    end

    return NEOParameters{T}(vals...)
end

# Print method for NEOParameters
function show(io::IO, p::NEOParameters{T}) where {T <: Real}
    avoid = [:μ_ast, :eph_su, :eph_ea]
    params = fieldnames(NEOParameters{T})
    s = Vector{String}(undef, length(params)+1)
    s[1] = "NEOParameters{$T}:\n"
    for i in eachindex(params)
        if params[i] in avoid
            s[i+1] = ""
        else
            x = string(params[i], ":")
            s[i+1] = string("    ", rpad(x, 15), getfield(p, params[i]), "\n")
        end
    end

    print(io, join(s))
end

# Load Sun (Earth) ephemeris
function _loadephsu()
    _su = selecteph(sseph, su)
    return TaylorInterpolant(_su.t0, _su.t, collect(_su.x))
end

function _loadephea()
    _ea = selecteph(sseph, ea)
    return TaylorInterpolant(_ea.t0, _ea.t, collect(_ea.x))
end