include("b_plane.jl")
include("least_squares.jl")

@auto_hash_equals struct NEOSolution{T <: Real, U <: Number}
    bwd::TaylorInterpolant{T, U, 2}
    t_bwd::Vector{U}
    x_bwd::Vector{U}
    g_bwd::Vector{U}
    fwd::TaylorInterpolant{T, U, 2}
    t_fwd::Vector{U}
    x_fwd::Vector{U}
    g_fwd::Vector{U}
    res::Vector{OpticalResidual{T, U}}
    fit::OrbitFit{T}
    function NEOSolution{T, U}(bwd::TaylorInterpolant{T, U, 2}, t_bwd::Vector{U}, x_bwd::Vector{U}, g_bwd::Vector{U},
                               fwd::TaylorInterpolant{T, U, 2}, t_fwd::Vector{U}, x_fwd::Vector{U}, g_fwd::Vector{U},
                               res::Vector{OpticalResidual{T, U}}, fit::OrbitFit{T}) where {T <: Real, U <: Number}
        @assert bwd.t0 == fwd.t0 "Backward and forward propagation initial times must match"
        new{T, U}(bwd, t_bwd, x_bwd, g_bwd, fwd, t_fwd, x_fwd, g_fwd, res, fit)
    end
end

function NEOSolution(bwd::TaylorInterpolant{T, U, 2}, t_bwd::Vector{U}, x_bwd::Vector{U}, g_bwd::Vector{U},
                     fwd::TaylorInterpolant{T, U, 2}, t_fwd::Vector{U}, x_fwd::Vector{U}, g_fwd::Vector{U},
                     res::Vector{OpticalResidual{T, U}}, fit::OrbitFit{T}) where {T <: Real, U <: Number}
    NEOSolution{T, U}(bwd, t_bwd, x_bwd, g_bwd, fwd, t_fwd, x_fwd, g_fwd, res, fit)
end

function NEOSolution(bwd::TaylorInterpolant{T, U, 2}, fwd::TaylorInterpolant{T, U, 2},
                     res::Vector{OpticalResidual{T, U}}, fit::OrbitFit{T}) where {T <: Real, U <: Number}
    # Backward roots
    t_bwd = Vector{U}(undef, 0)
    x_bwd = Vector{U}(undef, 0)
    g_bwd = Vector{U}(undef, 0)
    # forward roots
    t_fwd = Vector{U}(undef, 0)
    x_fwd = Vector{U}(undef, 0)
    g_fwd = Vector{U}(undef, 0)

    return NEOSolution{T, U}(bwd, t_bwd, x_bwd, g_bwd, fwd, t_fwd, x_fwd, g_fwd, res, fit)
end

# Print method for NEOSolution
# Examples:
# NEO solution with 123 residuals
function show(io::IO, x::NEOSolution{T, U}) where {T <: Real, U <: Number}
    print(io, "NEO solution with ", length(x.res), " residuals")
end

# Evaluation in time method
function (sol::NEOSolution{T, U})(t::V) where {T <: Real, U <: Number, V <: Number}
    if t <= sol.bwd.t0
        return sol.bwd(t)
    else
        return sol.fwd(t)
    end
end

# Evaluation in fit δs method
function evalfit(sol::NEOSolution{T, TaylorN{T}}) where {T <: Real}
    # Fit δs
    δs = sol.fit.x
    # Evaluate backward integration
    new_bwd_x = map(x -> Taylor1(x.coeffs(δs)), sol.bwd.x);
    new_bwd = TaylorInterpolant{T, T, 2}(sol.bwd.t0, sol.bwd.t, new_bwd_x)
    # Evaluate backward roots
    new_t_bwd = sol.t_bwd(δs)
    new_x_bwd = sol.x_bwd(δs)
    new_g_bwd = sol.g_bwd(δs)
    # Evaluate forward integration
    new_fwd_x = map(x -> Taylor1(x.coeffs(δs)), sol.fwd.x);
    new_fwd = TaylorInterpolant{T, T, 2}(sol.fwd.t0, sol.fwd.t, new_fwd_x)
    # Evaluate forward roots
    new_t_fwd = sol.t_fwd(δs)
    new_x_fwd = sol.x_fwd(δs)
    new_g_fwd = sol.g_fwd(δs)
    # Evaluate residuals 
    new_res = sol.res(δs)

    return NEOSolution{T, T}(new_bwd, new_t_bwd, new_x_bwd, new_g_bwd,
                             new_fwd, new_t_fwd, new_x_fwd, new_g_fwd,
                             new_res, sol.fit)
end

function zero(::Type{NEOSolution{T, U}}) where {T <: Real, U <: Number}
    bwd = TaylorInterpolant{T, U, 2}(zero(T), zeros(T, 1), Matrix{Taylor1{U}}(undef, 0, 0))
    t_bwd = Vector{U}(undef, 0)
    x_bwd = Vector{U}(undef, 0)
    g_bwd = Vector{U}(undef, 0)
    fwd = TaylorInterpolant{T, U, 2}(zero(T), zeros(T, 1), Matrix{Taylor1{U}}(undef, 0, 0))
    t_fwd = Vector{U}(undef, 0)
    x_fwd = Vector{U}(undef, 0)
    g_fwd = Vector{U}(undef, 0)
    res = Vector{OpticalResidual{T, U}}(undef, 0)
    fit = OrbitFit{T}(false, Vector{T}(undef, 0), Matrix{T}(undef, 0, 0), :newton)
    return NEOSolution{T, U}(bwd, t_bwd, x_bwd, g_bwd, fwd, t_fwd, x_fwd, g_fwd,
                             res, fit)
end

iszero(x::NEOSolution{T, U}) where {T <: Real, U <: Number} = x == zero(NEOSolution{T, U})

function orbitdetermination(radec::Vector{RadecMPC{T}}, dynamics::D, maxsteps::Int, jd0::T, q0::Vector{U};
                            debias_table::String = "2018",  kwargs...) where {T <: Real, U <: Number, D}
    mpc_catalogue_codes_201X, truth, resol, bias_matrix = select_debiasing_table(debias_table)
    return orbitdetermination(radec, dynamics, maxsteps, jd0, q0,
                              mpc_catalogue_codes_201X, truth, resol, bias_matrix; kwargs...)
end 

function orbitdetermination(radec::Vector{RadecMPC{T}}, dynamics::D, maxsteps::Int, jd0::T, q0::Vector{U},
                            mpc_catalogue_codes_201X::Vector{String}, truth::String, resol::Resolution,
                            bias_matrix::Matrix{T}; niter::Int = 5, order::Int = order, abstol::T = abstol,
                            parse_eqs::Bool = true, μ_ast::Vector = μ_ast343_DE430[1:end]) where {T <: Real, U <: Number, D}

    # Sun's ephemeris
    eph_su = selecteph(sseph, su)
    # Earth's ephemeris
    eph_ea = selecteph(sseph, ea)
    # Julian day of first (last) observation
    t0, tf = datetime2julian(date(radec[1])), datetime2julian(date(radec[end]))
    # Number of years in forward integration 
    nyears_fwd = (tf - jd0 + 2) / yr
    # Number of years in backward integration
    nyears_bwd = -(jd0 - t0 + 2) / yr

    # Backward integration
    bwd = propagate(dynamics, maxsteps, jd0, nyears_bwd, q0; μ_ast = μ_ast, order = order,
                    abstol = abstol, parse_eqs = parse_eqs)

    if bwd.t[end] > t0 - jd0
        return zero(NEOSolution{T, T})
    end

    # Forward integration
    fwd = propagate(dynamics, maxsteps, jd0, nyears_fwd, q0; μ_ast = μ_ast, order = order,
                    abstol = abstol, parse_eqs = parse_eqs)

    if fwd.t[end] < tf - jd0
        return zero(NEOSolution{T, T})
    end

    # Residuals
    res = residuals(radec, mpc_catalogue_codes_201X, truth, resol, bias_matrix;
                    xvs = et -> auday2kmsec(eph_su(et/daysec)), xve = et -> auday2kmsec(eph_ea(et/daysec)),
                    xva = et -> bwdfwdeph(et, bwd, fwd))
    # Orbit fit
    fit = tryls(res, zeros(get_numvars()), niter)

    if nrms(res, fit) < 1
        return evalfit(NEOSolution(bwd, fwd, res, fit))
    end
    
    # Orbit fit / outlier rejection loop
    mask, _mask_ = outlier.(res), outlier.(res)
    for i in 1:10
        # Update outliers
        res = outlier_rejection(res, fit)
        # Update fit
        fit = tryls(res, zeros(get_numvars()), niter)

        _mask_ = outlier.(res)
        if all( mask .== _mask_ )
            break
        else
            mask = _mask_
        end
    end
    
    return evalfit(NEOSolution(bwd, fwd, res, fit))

end

function nrms(sol::NEOSolution{T, T}) where {T <: Real}
    if iszero(sol)
        return T(Inf)
    else
        return nrms(sol.res)
    end
end