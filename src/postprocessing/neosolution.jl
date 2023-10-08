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

@doc raw"""
    residual_norm(x::OpticalResidual{T, T}) where {T <: Real}

Return the contribution of `x` to the nrms.
"""
residual_norm(x::OpticalResidual{T, T}) where {T <: Real} = x.w_α * x.ξ_α^2 / x.relax_factor + x.w_δ * x.ξ_δ^2 / x.relax_factor

function orbitdetermination(radec::Vector{RadecMPC{T}}, sol::NEOSolution{T, T};
                            debias_table::String = "2018",  kwargs...) where {T <: Real}
    mpc_catalogue_codes_201X, truth, resol, bias_matrix = select_debiasing_table(debias_table)
    return orbitdetermination(radec, sol, mpc_catalogue_codes_201X, truth, resol, bias_matrix; kwargs...)
end 

function orbitdetermination(radec::Vector{RadecMPC{T}}, sol::NEOSolution{T, T}, mpc_catalogue_codes_201X::Vector{String},
                            truth::String, resol::Resolution, bias_matrix::Matrix{T}; max_per = 18., niter::Int = 5,
                            order::Int = order, varorder::Int = 5, abstol::T = abstol, parse_eqs::Bool = true,
                            μ_ast::Vector = μ_ast343_DE430[1:end]) where {T <: Real}

    # Sun's ephemeris
    eph_su = selecteph(sseph, su)
    # Earth's ephemeris
    eph_ea = selecteph(sseph, ea)
    # Julian day to start propagation
    jd0 = sol.bwd.t0 + PE.J2000
    # Julian day of first (last) observation
    t0, tf = datetime2julian(date(radec[1])), datetime2julian(date(radec[end]))
    # Number of years in forward integration 
    nyears_fwd = (tf - jd0 + 2) / yr
    # Number of years in backward integration
    nyears_bwd = -(jd0 - t0 + 2) / yr
    # Dynamical function 
    dynamics = RNp1BP_pN_A_J23E_J2S_eph_threads!
    # Maximum number of steps
    maxsteps = adaptative_maxsteps(radec)
    # Initial conditions (T)
    q00 = sol(sol.bwd.t0)
    # Jet transport perturbation
    dq = scaled_variables("δx", fill(1e-6, 6); order = varorder)
    # Initial conditions (jet transport)
    q0 = q00 .+ dq

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

    # NRMS (with 0 outliers)
    Q_0 = nrms(res, fit)

    if Q_0 < 1
        return evalfit(NEOSolution(bwd, fwd, res, fit))
    end

    # Number of observations
    N_radec = length(radec)
    # Maximum allowed outliers
    max_drop = ceil(Int, N_radec * max_per / 100)
    # Boolean mask (0: included in fit, 1: outlier)
    new_outliers = BitVector(zeros(N_radec))

    # Drop loop
    for i in 1:max_drop
        # Contribution of each residual to nrms
        norms = residual_norm.(res(fit.x))
        # Iterate norms from largest to smallest
        idxs = sortperm(norms, rev = true)

        for j in idxs
            if !new_outliers[j]
                # Drop residual
                new_outliers[j] = true
                # Update residuals
                res = OpticalResidual.(ra.(res), dec.(res), weight_ra.(res), weight_dec.(res), relax_factor.(res), new_outliers)
                # Update fit
                fit = tryls(res, zeros(get_numvars()), niter)
                break
            end
        end
    end

    # Outliers
    idxs = Vector{Int}(undef, max_drop)
    # NRMS
    Qs = Vector{T}(undef, max_drop + 1)
    # Number of outliers
    N_outliers = Vector{T}(undef, max_drop + 1)

    # Recovery loop
    for i in 1:max_drop
        # NRMS of current fit
        Qs[i] = nrms(res, fit)
        # Number of outliers in current fit
        N_outliers[i] = float(max_drop - i + 1)
        # Contribution of each residual to nrms
        norms = residual_norm.(res(fit.x))
        # Minimum norm among outliers
        j = findmin(norms[new_outliers])[2]
        # Find residual with minimum norm
        j = findall(new_outliers)[j]
        # Add j-th residual to outliers list
        idxs[i] = j
        # Recover residual
        new_outliers[j] = false
        # Update residuals
        res = OpticalResidual.(ra.(res), dec.(res), weight_ra.(res), weight_dec.(res), relax_factor.(res), new_outliers)
        # Update fit
        fit = tryls(res, zeros(get_numvars()), niter)
    end
    # Add 0 outliers fit 
    Qs[end] = Q_0
    N_outliers[end] = zero(T)

    # Outlier rejection cannot reduce Q
    if all(Qs .> 1.)
        # Reset boolean mask
        new_outliers[1:end] .= false
        # Update residuals
        res = OpticalResidual.(ra.(res), dec.(res), weight_ra.(res), weight_dec.(res), relax_factor.(res), new_outliers)
        # Update fit
        fit = tryls(res, zeros(get_numvars()), niter)

        return evalfit(NEOSolution(bwd, fwd, res, fit)) 
    end

    # Assemble points
    points = [[Qs[i], N_outliers[i]] for i in eachindex(Qs)]
    # Leave only first point with Q > 1
    i_1 = findfirst(x -> x[1] > 1, points)
    points = points[1:i_1]
    # Number of remaining points
    N_points = length(points)
    # Find optimal fit
    if N_points == 0
        idxs = Vector{Int}(undef, 0)
    elseif N_points <= 2
        idxs = idxs[N_points:end]
    else
        # Difference in NRMS
        dQ = Base.diff(first.(points))
        # Mean difference 
        mQ = mean(dQ)
        # Standard deviation of difference
        sQ = std(dQ)
        # Maximum difference
        dQ_max, i_max = findmax(dQ)

        if dQ_max > mQ + sQ
            idxs = idxs[i_max:end]
        else 
            # Remove points with Q < 1
            filter!(x -> x[1] < 1, points)
            # Mean point 
            avg = sum(points) ./ length(points)
            # Distance from each point to mean points
            diff = [norm([Qs[i], N_outliers[i]] .- avg) for i in 1:max_drop]
            # Find pair closest to mean point
            i = findmin(diff)[2]
            idxs = idxs[i:end]
        end
    end

    # Reset boolean mask
    new_outliers[1:end] .= false
    # Outliers
    new_outliers[idxs] .= true
    # Update residuals
    res = OpticalResidual.(ra.(res), dec.(res), weight_ra.(res), weight_dec.(res), relax_factor.(res), new_outliers)
    # Update fit
    fit = tryls(res, zeros(get_numvars()), niter)

    return evalfit(NEOSolution(bwd, fwd, res, fit))
end

function nrms(sol::NEOSolution{T, T}) where {T <: Real}
    if iszero(sol)
        return T(Inf)
    else
        return nrms(sol.res)
    end
end