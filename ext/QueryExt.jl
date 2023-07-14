module QueryExt

using Dates: DateTime, Date, DatePeriod, Day, datetime2julian
using TaylorSeries: get_numvars
using PlanetaryEphemeris: J2000, selecteph, su, ea, yr, daysec, auday2kmsec
using NEOs: ObservatoryMPC, RadecMPC, observatory, date, polynomial_interpolation, gauss_idxs, propagate,
            RNp1BP_pN_A_J23E_J2S_eph_threads!, order, abstol, sseph, scaled_variables, gauss_method, residuals,
            bwdfwdeph, newtonls, diffcorr, nrms, hascoord
import NEOs: reduce_nights, gaussinitcond

if isdefined(Base, :get_extension)
    using Query
else
    using ..Query
end

@doc raw"""
    reduce_nights(radec::Vector{RadecMPC{T}}) where {T <: AbstractFloat}

Return one observatory, date, right ascension and declination per observation night in `radec`. The reduction is performed 
via Laplace interpolation. 
"""
function reduce_nights(radec::Vector{RadecMPC{T}}) where {T <: AbstractFloat}
    # Observations per observatory 
    radec_per_obs = radec |> 
                    @filter(hascoord(observatory(_))) |>
                    @groupby(observatory(_)) |> 
                    @filter(length(_) >= 3) .|> 
                    collect
    # Initialize output
    output = Vector{Tuple{ObservatoryMPC{T}, DateTime, T, T}}(undef, 0)
    # Find observation nights per observatory 
    for i in eachindex(radec_per_obs)
        # Current observatory 
        obs = observatory(radec_per_obs[i][1])
        # Reduce nights by interpolation
        interp = radec_per_obs[i] |> 
                 @groupby(Date(date(_))) |> 
                 @filter(length(_) >= 3) |> 
                 @map(polynomial_interpolation(_)) |> 
                 @filter( !isnan(_[2]) && !isnan(_[3]) ) |> 
                 collect 
        for j in eachindex(interp)
            output = vcat(output, (obs, interp[j]...))
        end 
    end 

    sort!(output, by = x -> x[2])

    return first.(output), getindex.(output, 2), getindex.(output, 3), last.(output)
end 

@doc raw"""
    gaussinitcond(radec::Vector{RadecMPC{T}}; Δ::DatePeriod = Day(1), niter::Int = 5, maxsteps::Int = 100, varorder::Int = 5,
                  order::Int = order, abstol::T = abstol, parse_eqs::Bool = true) where {T <: AbstractFloat}

Return initial conditions via Gauss method. 

See also [`gauss_method`](@ref).

# Arguments
- `radec::Vector{RadecMPC{T}}`: vector of oobservations.
- `Δ::DatePeriod`: see [`gauss_idxs`](@ref).
- `niter::Int`: number of iterations for Newton's method.
- `maxsteps::Int`: maximum number of steps for propagation.
- `varorder::Int`: order of jet transport perturbation. 
- `order::Int`: order of Taylor polynomials w.r.t. time.
- `abstol::T`: absolute tolerance.
- `parse_eqs::Bool`: whether to use the taylorized method of `RNp1BP_pN_A_J23E_J2S_eph_threads!` or not. 

!!! warning
    This function will set the (global) `TaylorSeries` variables to `δα₁ δα₂ δα₃ δδ₁ δδ₂ δδ₃`. 
"""
function gaussinitcond(radec::Vector{RadecMPC{T}}; Δ::DatePeriod = Day(1), niter::Int = 5, maxsteps::Int = 100, varorder::Int = 5,
                       order::Int = order, abstol::T = abstol, parse_eqs::Bool = true) where {T <: AbstractFloat}

    filter!(x -> hascoord(observatory(x)), radec)

    # Sun's ephemeris
    eph_su = selecteph(sseph, su)
    # Earth's ephemeris
    eph_ea = selecteph(sseph, ea)

    # Reduce nights by interpolation 
    observatories, dates, α, δ = reduce_nights(radec)
    # Pick three observations for Gauss method 
    idxs = gauss_idxs(dates, Δ) 
    
    # [julian days]
    t0, jd0, tf = datetime2julian.(dates[idxs])
    # Number of years in forward integration 
    nyears_fwd = (tf - jd0 + 2) / yr
    # Number of years in backward integration
    nyears_bwd = -(jd0 - t0 + 2) / yr

    # Subset of radec for residuals
    sub_radec = filter(x -> dates[idxs[1]] <= date(x) <= dates[idxs[3]], radec)

    # Jet transport perturbation (ra/dec)
    dq = scaled_variables("δα₁ δα₂ δα₃ δδ₁ δδ₂ δδ₃"; order = varorder)
    # Gauss method solution 
    sol = gauss_method(observatories[idxs], dates[idxs], α[idxs] + dq[1:3], δ[idxs] .+ dq[4:6]; niter = niter)

    # Vector of errors 
    Q = Vector{T}(undef, length(sol))
    # Vector of initial conditions
    Q0 = Matrix{T}(undef, length(sol), 6)

    # Iterate over Gauss solutions 
    for i in eachindex(sol)
        # Initial conditions (jet transport)
        q0 = sol[i].statevect .+ eph_su(jd0 - J2000)

        # Propagation 
        bwd, fwd = propagate(RNp1BP_pN_A_J23E_J2S_eph_threads!, maxsteps, jd0, nyears_bwd, nyears_fwd, q0, Val(true); 
                             order = order, abstol = abstol, parse_eqs = parse_eqs)

        # O-C residuals
        res, w = residuals(sub_radec; xvs = et -> auday2kmsec(eph_su(et/daysec)), xve = et -> auday2kmsec(eph_ea(et/daysec)), 
                           xva = et -> bwdfwdeph(et, bwd, fwd))

        # Orbit fit 
        success, x_new, Γ = newtonls(res, w, zeros(get_numvars()), niter)

        # TO DO: check cases where newton converges but diffcorr no

        if success
            # NRMS of the solution 
            Q[i] = nrms(res(x_new), w)
            # Initial conditions 
            Q0[i, :] = bwd(bwd.t0)(x_new)
        else
            success, x_new, Γ = diffcorr(res, w, zeros(get_numvars()), niter)
            if success
                Q[i] = nrms(res(x_new), w)
                Q0[i, :] = bwd(bwd.t0)(x_new)
            else
                Q[i] = T(Inf)
                Q0[i, :] .= T(Inf)
            end 
        end 
    end 

    # Solution with minimum NRMS
    i = findmin(Q)[2]
    # Case: all solutions were unsuccesful
    if isinf(Q[i])
        q00 = Vector{T}(undef, 0)
    # Case: at least one solution was succesful
    else 
        q00 = Q0[i, :]
    end

    return jd0, q00
end 

end # module