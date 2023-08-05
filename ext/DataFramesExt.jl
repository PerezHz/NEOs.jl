module DataFramesExt

using Dates: Date, Period, Hour, Day, datetime2julian, julian2datetime
using TaylorSeries: get_numvars
using PlanetaryEphemeris: J2000, selecteph, su, ea, yr, daysec, auday2kmsec
using NEOs: RadecMPC, date, gauss_triplets, propagate, RNp1BP_pN_A_J23E_J2S_eph_threads!, order, abstol, sseph,
            scaled_variables, gauss_method, residuals, bwdfwdeph, newtonls, diffcorr, nrms, hascoord

import Base: convert
import NEOs: AbstractAstrometry, extrapolation, reduce_nights, gaussinitcond

if isdefined(Base, :get_extension)
    using DataFrames: AbstractDataFrame, DataFrame, nrow, eachcol, eachrow, groupby, combine
    import Tables: istable, rowaccess, rows, schema, Schema
else
    using ..DataFrames: AbstractDataFrame, DataFrame, nrow, eachcol, eachrow, groupby, combine
    import ..Tables: istable, rowaccess, rows, schema, Schema
end

# Methods to convert a Vector{<:AbstractAstrometry} to a DataFrame
istable(::Type{Vector{<:AbstractAstrometry}}) = true
rowaccess(::Type{Vector{<:AbstractAstrometry}}) = true
rows(x::Vector{<:AbstractAstrometry}) = x
schema(::Vector{T}) where {T <: AbstractAstrometry} = Schema(fieldnames(T), Tuple{fieldtypes(T)...})

# Methods to convert a DataFrame to a Vector{<:AbstractAstrometry}
function Vector{T}(df::DataFrame) where {T <: AbstractAstrometry}
    @assert all(String.(fieldnames(T)) .== names(df)) "`DataFrame` column names don't match `$T` fieldnames"
    @assert all(fieldtypes(T) .== eltype.(eachcol(df))) "`DataFrame` column types don't match `$T` fieldtypes"
    obs = Vector{T}(undef, nrow(df)) 
    for (i, row) in zip(eachindex(obs), eachrow(df))
        obs[i] = T(values(row)...)
    end 
    return obs 
end 

convert(::Type{Vector{T}}, df::DataFrame) where {T <: AbstractAstrometry} = Vector{T}(df)

@doc raw"""
    extrapolation(df::AbstractDataFrame)

Special method of [`extrapolation`](@ref) to be used by [`gaussinitcond`](@ref).
"""
function extrapolation(df::AbstractDataFrame)
    if isone(nrow(df))
        return (observatory = df.observatory[1], date = df.date[1], α = df.α[1], δ = df.δ[1])
    end 
    
    # Julian days of observation
    t_julian = datetime2julian.(df.date)
    # Days of observation [relative to first observation]
    t_rel = t_julian .- t_julian[1]
    # Mean date 
    t_mean = sum(t_rel) / length(t_rel)

    # Extrapolate
    α_p = extrapolation(t_rel, df.α)
    δ_p = extrapolation(t_rel, df.δ)
    # Evaluate polynomials at mean date 
    α_mean = α_p(t_mean)
    δ_mean = δ_p(t_mean)

    return (observatory = df.observatory[1], date = julian2datetime(t_julian[1] + t_mean), α = α_mean, δ = δ_mean)
end 

@doc raw"""
    reduce_nights(radec::Vector{RadecMPC{T}}) where {T <: AbstractFloat}

Return one observatory, date, right ascension and declination per observation night in `radec`. The reduction is performed 
via polynomial interpolation. 
"""
function reduce_nights(radec::Vector{RadecMPC{T}}) where {T <: AbstractFloat}
    # Convert to DataFrame 
    df = DataFrame(radec)
    # Eliminate observatories without coordinates 
    filter!(:observatory => hascoord, df)
    # Group by observatory and Date 
    df.Date = Date.(df.date)
    gdf = groupby(df, [:observatory, :Date])
    # Interpolate observation nights 
    cdf = combine(gdf, extrapolation, keepkeys = false)
    # Eliminate unsuccesful interpolations 
    filter!(:α => !isnan, cdf)
    filter!(:δ => !isnan, cdf)
    # Sorbt by date 
    sort!(cdf, :date)

    return cdf.observatory, cdf.date, cdf.α, cdf.δ
end 

@doc raw"""
    gaussinitcond(radec::Vector{RadecMPC{T}}; Δ::Period = Day(1), Q_max::T = 0.75, niter::Int = 5, maxsteps::Int = 100, 
                  varorder::Int = 5, order::Int = order, abstol::T = abstol, parse_eqs::Bool = true) where {T <: AbstractFloat}

Return initial conditions via Gauss method. 

See also [`gauss_method`](@ref).

# Arguments
- `radec::Vector{RadecMPC{T}}`: vector of observations.
- `Δ::Period`: see [`gauss_triplets`](@ref).
- `Q_max::T`: The maximum nrms that is considered a good enough orbit.
- `niter::Int`: number of iterations for Newton's method.
- `maxsteps::Int`: maximum number of steps for propagation.
- `varorder::Int`: order of jet transport perturbation. 
- `order::Int`: order of Taylor polynomials w.r.t. time.
- `abstol::T`: absolute tolerance.
- `parse_eqs::Bool`: whether to use the taylorized method of `RNp1BP_pN_A_J23E_J2S_eph_threads!` or not. 

!!! warning
    This function will set the (global) `TaylorSeries` variables to `δα₁ δα₂ δα₃ δδ₁ δδ₂ δδ₃`. 
"""
function gaussinitcond(radec::Vector{RadecMPC{T}}; Δ_min::Period = Hour(20), Δ_max::Period = Day(7), max_triplets::Int = 10,
                       Q_max::T = 0.75, niter::Int = 5, maxsteps::Int = 100, varorder::Int = 5, order::Int = order, 
                       abstol::T = abstol, parse_eqs::Bool = true) where {T <: AbstractFloat}

    # Sun's ephemeris
    eph_su = selecteph(sseph, su)
    # Earth's ephemeris
    eph_ea = selecteph(sseph, ea)

    # Reduce nights by interpolation 
    observatories, dates, α, δ = reduce_nights(radec)
    # Observations triplets
    triplets = gauss_triplets(dates; Δ_min = Δ_min, Δ_max = Δ_max, max_triplets = max_triplets) 

    # Initial date of integration [julian days]
    jd0 = zero(T)

    # Jet transport perturbation (ra/dec)
    dq = scaled_variables("δα₁ δα₂ δα₃ δδ₁ δδ₂ δδ₃"; order = varorder)

    # Normalized root mean square error (NRMS)
    best_Q = T(Inf)
    # Initial conditions
    best_Q0 = zeros(T, 6)

    # Global counter
    k = 1
    # Break flag
    flag = false

    # Iterate over triplets
    for j in eachindex(triplets)
        
        # Current triplet
        idxs = triplets[j]

        # [julian days]
        t0, jd0, tf = datetime2julian.(dates[idxs])
        # Number of years in forward integration 
        nyears_fwd = (tf - jd0 + 2) / yr
        # Number of years in backward integration
        nyears_bwd = -(jd0 - t0 + 2) / yr

        # Subset of radec for residuals
        sub_radec = filter(x -> dates[idxs[1]] <= date(x) <= dates[idxs[3]], radec)

        # Gauss method solution 
        sol = gauss_method(observatories[idxs], dates[idxs], α[idxs] .+ dq[1:3], δ[idxs] .+ dq[4:6]; niter = niter)

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

            # Orbit fit (Newton)
            success, x_new, Γ = newtonls(res, w, zeros(get_numvars()), niter)
            # NRMS
            Q = nrms(res(x_new), w)

            # TO DO: check cases where newton converges but diffcorr no

            if success
                # Update NRMS and initial conditions
                if Q < best_Q
                    best_Q = Q
                    best_Q0 .= bwd(bwd.t0)(x_new)
                end 
                # Break condition
                if Q <= Q_max
                    flag = true
                end
            else
                # Orbit fit (differential corrections)
                success, x_new, Γ = diffcorr(res, w, zeros(get_numvars()), niter)
                # NRMS
                Q = nrms(res(x_new), w)

                if success
                    # Update NRMS and initial conditions
                    if Q < best_Q
                        best_Q = Q
                        best_Q0 .= bwd(bwd.t0)(x_new)
                    end
                    # Break condition
                    if Q <= Q_max
                        flag = true
                    end
                end 
            end 
            k += 1
            # Break condition
            if flag
                break
            end 
        end
        if flag
            break
        end
    end 

    # Case: all solutions were unsuccesful
    if isinf(best_Q)
        return jd0::T, Vector{T}(undef, 0)
    # Case: at least one solution was succesful
    else 
        return jd0::T, best_Q0
    end
end 

end # module