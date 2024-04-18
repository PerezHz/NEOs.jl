struct NEOParameters{T <: AbstractFloat}
    # Propagation parameters
    maxsteps::Int
    μ_ast::Vector{T}
    order::Int
    abstol::T
    parse_eqs::Bool
    bwdoffset::T
    fwdoffset::T
    coeffstol::T
    # Residuals parameters
    mpc_catalogue_codes_201X::Vector{String}
    truth::String
    resol::Resolution
    bias_matrix::Matrix{T}
    eph_su::TaylorInterpolant{T, T, 2, Vector{T}, Matrix{Taylor1{T}}}
    eph_ea::TaylorInterpolant{T, T, 2, Vector{T}, Matrix{Taylor1{T}}}
    # Least squares fit parameters
    niter::Int
    # Gauss Method parameters
    max_triplets::Int
    varorder::Int
    gaussQmax::T
    # Admissible Region parameters
    H_max::T
    a_max::T
    maxiter::Int
    tsaQmax::T
    # Outlier rejection parameters
    max_per::T
    # Inner constructor is generated by default
end

# Print method for NEOParameters
function show(io::IO, p::NEOParameters{T}) where {T <: AbstractFloat}
    params = [
        :maxsteps, :order, :abstol, :parse_eqs, :bwdoffset,
        :fwdoffset, :coeffstol, :niter, :max_triplets, :varorder,
        :gaussQmax, :H_max, :a_max, :maxiter, :tsaQmax, :max_per
    ]
    s = Vector{String}(undef, length(params)+1)
    s[1] = "NEOParameters{$T}:\n"
    for i in eachindex(params)
        x = string(params[i], ":")
        s[i+1] = string("    ", rpad(x, 15), getfield(p, params[i]), "\n")
    end

    print(io, join(s))
end

# Outer constructors

@doc raw"""
    NEOParameters([params::NEOParameters{T};] kwargs...) where {T <: AbstractFloat}

Parameters for all orbit determination functions.

# Propagation Parameters

- `maxsteps::Int`: maximum number of steps for the integration (default: `500`).
- `μ_ast::Vector{T}`: vector of gravitational parameters (default: `μ_ast343_DE430[1:end]`).
- `order::Int`: order of Taylor expansions wrt time (default: 25).
- `abstol::T`: absolute tolerance used to compute propagation timestep (default: `1e-20`).
- `parse_eqs::Bool`: whether to use the specialized method of `jetcoeffs` or not (default: `true`).
- `bwdoffset/fwdoffset::T`: days to propagate beyond first (bwd) / last (fwd) observation (default: `0.5`).
- `coeffstol::T`: maximum size of the coefficients (default: `10.0`).

# Residuals Parameters

- `debias_table::String`: debiasing scheme (default: `"2018"`). Possible values are:
    - `"2014"` corresponds to https://doi.org/10.1016/j.icarus.2014.07.033,
    - `"2018"` corresponds to https://doi.org/10.1016/j.icarus.2019.113596,
    - `"hires2018"` corresponds to https://doi.org/10.1016/j.icarus.2019.113596.

# Least Squares Fit Parameters

- `niter::Int`: number of iterations for differential corrections / Newton's method (default: `5`).

# Gauss Method Parameters

- `max_triplets::Int`: maximum number of triplets to check for a solution (default: `10`).
- `varorder::Int`: order of jet transport perturbation (default: `5`).
- `gaussQmax::T`: nrms threshold (default: `5.0`).

# Admissible Region Parameters

- `H_max::T`: maximum absolute magnitude (default: `34.5`).
- `a_max::T`: maximum semimajor axis (default: `100.0`).
- `maxiter::Int`: maximum number of iterations for admissible region `ADAM` optimizer (default: `200`).
- `tsaQmax::T`: nrms threshold (default: `1.5`).

# Outlier Rejection Parameters

- `max_per::T`: maximum allowed rejection percentage (default: `18.0`).
"""
function NEOParameters(;
    maxsteps::Int = 500, μ_ast::Vector{T} = μ_ast343_DE430[1:end], order::Int = 25,
    abstol::T = 1e-20, parse_eqs::Bool = true, bwdoffset::T = 0.5, fwdoffset::T = 0.5,
    coeffstol::T = 10.0, debias_table::String = "2018", niter::Int = 5,
    max_triplets::Int = 10, varorder::Int = 5, gaussQmax::T = 5.0, H_max::T = 34.5,
    a_max::T = 100.0, maxiter::Int = 200, tsaQmax::T = 1.5, max_per::T = 18.0
    ) where {T <: AbstractFloat}
    # Unfold debiasing matrix
    mpc_catalogue_codes_201X, truth, resol, bias_matrix = select_debiasing_table(debias_table)
    # Sun (Earth) ephemeris
    _su = selecteph(sseph, su)
    eph_su = TaylorInterpolant(_su.t0, _su.t, collect(_su.x))
    _ea = selecteph(sseph, ea)
    eph_ea = TaylorInterpolant(_ea.t0, _ea.t, collect(_ea.x))
    # Assemble NEOParameters
    return NEOParameters{T}(
        maxsteps, μ_ast, order, abstol, parse_eqs, bwdoffset, fwdoffset,
        coeffstol, mpc_catalogue_codes_201X, truth, resol, bias_matrix,
        eph_su, eph_ea, niter, max_triplets, varorder, gaussQmax, H_max,
        a_max, maxiter, tsaQmax, max_per
    )
end

function NEOParameters(params::NEOParameters{T}; kwargs...) where {T <: AbstractFloat}
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

function residuals(obs::Vector{RadecMPC{T}}, params::NEOParameters{T}; kwargs...) where {T <: AbstractFloat}
    return residuals(obs, params.mpc_catalogue_codes_201X, params.truth, params.resol,
                     params.bias_matrix; kwargs...)
end