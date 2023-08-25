include("osculating.jl")

@doc raw"""
    GaussSolution{T <: Real, U <: Number}

A preliminary orbit obtained from Gauss method of orbit determination.

See also [`gauss_method`](@ref).

# Fields

- `statevect::Vector{U}`: state vector at middle observation.
- `D::Matrix{U}`: D matrix.
- `R_vec::Matrix{T}`: observer's heliocentric positions.
- `ρ_vec::Matrix{U}`: slant ranges.
- `τ_1::T`: time between first and second observations.
- `τ_3::T`: time between third and second observations.
- `f_1, g_1, f_3, g_3::U`: Lagrange coefficients.
- `status::Symbol`: the status of the solution (`:empty`, `:unkown` or `:unique`)

!!! reference
    See Algorithm 5.5 in page 274 of https://doi.org/10.1016/C2016-0-02107-1.
"""
@auto_hash_equals struct GaussSolution{T <: Real, U <: Number}
    statevect::Vector{U}
    D::Matrix{U}
    R_vec::Matrix{T}
    ρ_vec::Matrix{U}
    τ_1::T
    τ_3::T
    f_1::U
    g_1::U
    f_3::U
    g_3::U
    status::Symbol
    # Internal constructor
    function GaussSolution{T, U}(statevect::Vector{U}, D::Matrix{U}, R_vec::Matrix{T}, ρ_vec::Matrix{U}, τ_1::T, τ_3::T,
                                 f_1::U, g_1::U, f_3::U, g_3::U, status::Symbol) where {T <: Real, U <: Number}
        return new{T, U}(statevect, D, R_vec, ρ_vec, τ_1, τ_3, f_1, g_1, f_3, g_3, status)
    end
end

# Outer constructor
function GaussSolution(statevect::Vector{U}, D::Matrix{U}, R_vec::Matrix{T}, ρ_vec::Matrix{U}, τ_1::T, τ_3::T,
                       f_1::U, g_1::U, f_3::U, g_3::U, status::Symbol) where {T <: Real, U <: Number}
    return GaussSolution{T, U}(statevect, D, R_vec, ρ_vec, τ_1, τ_3, f_1, g_1, f_3, g_3, status)
end

# Print method for GaussSolution
# Examples:
# unique Gauss solution (r = 1.0800950907383229 AU)
function show(io::IO, g::GaussSolution{T, U}) where {T <: Real, U <: Number}
    print(io, uppercasefirst(string(g.status)), " Gauss solution (r = ", norm(cte.(g.statevect[1:3])), " AU)")
end

@doc raw"""
    topounit(α::T, δ::T) where {T <: Number}
    topounit(obs::RadecMPC{T}) where {T <: AbstractFloat}

Return the topocentric unit vector.
"""
function topounit(α::T, δ::T) where {T <: Number}

    sin_α, cos_α = sincos(α)
    sin_δ, cos_δ = sincos(δ)

    pos = [cos_δ * cos_α, cos_δ * sin_α, sin_δ]

    return pos
end

topounit(obs::RadecMPC{T}) where {T <: AbstractFloat} = topounit(obs.α, obs.δ)

# TO DO: Should we allow to use other μ?

@doc raw"""
    f_Lagrange(τ::T, r::U) where {T <: Real, U <: Number}

Return the 1st order approximation to Lagrange's f function.
"""
function f_Lagrange(τ::T, r::U) where {T <: Real, U <: Number}
    return 1 - μ_S * (τ^2) / 2 / (r^3)
end

@doc raw"""
    g_Lagrange(τ::T, r::U) where {T <: Real, U <: Number}

Return the 1st order approximation to Lagrange's g function.
"""
function g_Lagrange(τ::T, r::U) where {T <: Real, U <: Number}
    return τ - μ_S * (τ^3) / 6 / (r^3)
end

@doc raw"""
    vectors2matrix(x::T) where {T <: AbstractVector}

Convert a vector of vectors `x` to a matrix.
"""
vectors2matrix(x::T) where {T <: AbstractVector} = permutedims(reduce(hcat, x))

@doc raw"""
    _format_Lagrange_equation(a::T, b::T, c::T) where {T <: AbstractFloat}

Format Lagrange equation as `r⁸ + a r⁶ + b r³ + c = 0`.
"""
function _format_Lagrange_equation(a::T, b::T, c::T) where {T <: AbstractFloat}
    a_sgn = a ≥ 0 ? "+" : "-"
    b_sgn = b ≥ 0 ? "+" : "-"
    c_sgn = c ≥ 0 ? "+" : "-"

    return join(["r⁸ ", a_sgn, " ", abs(a), " r⁶ ", b_sgn, " ", abs(b), " r³ ", c_sgn, " ", abs(c), " = 0"])
end

@doc raw"""
    lagrange(x::T, a::U, b::U, c::U) where {T, U <: Number}

Lagrange polynomial to be solved during Gauss method.
"""
lagrange(x::T, a::U, b::U, c::U) where {T, U <: Number} = x^8 + a*x^6 + b*x^3 + c

@doc raw"""
    lagrange_derivative(x::T, a::T, b::T) where {T <: Number}

Derivative of Lagrange polynomial to be solved during Gauss method.
"""
lagrange_derivative(x::T, a::U, b::U) where {T, U <: Number} = 8*x^7 + 6*a*x^5 + 3*b*x^2

# TO DO: Allow to control interval over which to look for solutions
# Currently we look between the radius of the Sun (0.00465047 au) and the radius of the Solar System (40 au)

@doc raw"""
    solve_lagrange(a::T, b::T, c::T; niter::Int = 5) where {T <: Real}
    solve_lagrange(a::TaylorN{T}, b::TaylorN{T}, c::TaylorN{T}; niter::Int = 5) where {T <: Real}

Solve Lagrange polynomial.

See also [`lagrange`](@ref).
"""
function solve_lagrange(a::T, b::T, c::T; niter::Int = 5) where {T <: Real}
    sol = roots(x -> lagrange(x, a, b, c), Interval(0.00465047, 40))
    return mid.(interval.(sol)), getfield.(sol, :status)
end

function solve_lagrange(a::TaylorN{T}, b::TaylorN{T}, c::TaylorN{T}; niter::Int = 5) where {T <: Real}
    # 0-th order solution
    sol_0, status = solve_lagrange(cte(a), cte(b), cte(c))
    # Discard empty solutions
    idxs = findall(x -> x != :empty, status)
    # Vector of solutions
    sol = Vector{TaylorN{T}}(undef, length(idxs))
    # Iterate non-empty solutions
    for i in eachindex(sol)
        # Newton's method
        r_0 = sol_0[idxs[i]]
        r_2 = sol_0[idxs[i]]
        for j in 1:niter
            r_2 = r_0 - lagrange(r_0, a, b, c) / lagrange_derivative(r_0, a, b)
            r_0 = r_2
        end
        sol[i] = r_2
    end

    return sol, status[idxs]
end

@doc raw"""
    heliocentric_energy(r::Vector{T}) where {T <: Number}

Return the heliocentric energy per unit mass for heliocentric state vector `r` [au, au/day].
"""
function heliocentric_energy(r::Vector{T}) where {T <: Number}
    @assert length(r) == 6 "r must have length 6"
    kinetic = 0.5 * (r[4]^2 + r[5]^2 + r[6]^2) 
    potential = k_gauss^2 / sqrt(r[1]^2 + r[2]^2 + r[3]^2)
    return kinetic - potential 
end

@doc raw"""
    gauss_method(obs::Vector{RadecMPC{T}}; niter::Int = 5) where {T <: AbstractFloat}
    gauss_method(observatories::Vector{ObservatoryMPC{T}}, dates::Vector{DateTime}, α::Vector{U}, δ::Vector{U};
                 niter::Int = 5) where {T <: Real, U <: Number}

Core Gauss method of Initial Orbit determination (IOD).

# Arguments

- `obs::Vector{RadecMPC{T}}`: three observations.
- `observatories::Vector{ObservatoryMPC{T}}`: sites of observation.
- `dates::Vector{DateTime}`: times of observation.
- `α::Vector{U}`: right ascension [rad].
- `δ::Vector{U}`: declination [rad].
- `niter::Int`: Number of iterations for Newton's method.

!!! reference
    See Algorithm 5.5 in page 274 https://doi.org/10.1016/C2016-0-02107-1.
"""
function gauss_method(obs::Vector{RadecMPC{T}}; niter::Int = 5) where {T <: AbstractFloat}

    # Make sure observations are in temporal order
    sort!(obs)

    # Sites of observation
    observatories = observatory.(obs)

    # Dates of observation
    dates = date.(obs)

    # Right ascension
    α = ra.(obs)

    # Declination
    δ = dec.(obs)

    return gauss_method(observatories, dates, α, δ; niter = niter)
end

function gauss_method(observatories::Vector{ObservatoryMPC{T}}, dates::Vector{DateTime}, α::Vector{U}, δ::Vector{U};
                      niter::Int = 5) where {T <: Real, U <: Number}

    # Check we have exactly three observations
    @assert length(observatories) == length(dates) == length(α) == length(δ) == 3 "Gauss method requires exactly three observations"
    # Check observations are in temporal order 
    @assert issorted(dates) "Observations must be in temporal order"

    # Times of observation [et]
    t_et = datetime2et.(dates)
    # Times of observation [days since J2000]
    t_days = t_et ./ daysec

    # Time intervals
    τ_1 = t_days[1] - t_days[2]
    τ_3 = t_days[3] - t_days[2]
    τ = τ_3 - τ_1

    # NEO's topocentric position unit vectors
    ρ_vec = vectors2matrix(topounit.(α, δ))
    # Geocentric state vector of the observer [au, au/day]
    g_vec = kmsec2auday.(obsposvelECI.(observatories, t_et))
    # Sun's ephemeris 
    eph_su = selecteph(sseph, su)
    # Earth's ephemeris 
    eph_ea = selecteph(sseph, ea)
    # Heliocentric state vector of the Earth [au, au/day]
    G_vec = eph_ea.(t_days) - eph_su.(t_days)
    # Observer's heliocentric positions [au, au/day]
    R_vec = vectors2matrix(G_vec .+  g_vec)[:, 1:3]

    # Cross products
    p_vec = zeros(U, 3, 3)
    p_vec[1, :] = cross(ρ_vec[2, :], ρ_vec[3, :])
    p_vec[2, :] = cross(ρ_vec[1, :], ρ_vec[3, :])
    p_vec[3, :] = cross(ρ_vec[1, :], ρ_vec[2, :])

    # Gauss scalar
    D_0 = dot(ρ_vec[1, :], p_vec[1, :])

    # Matrix of triple products
    D = zeros(U, 3, 3)
    for i in 1:3
        for j in 1:3
            D[i, j] = dot(R_vec[i, :], p_vec[j, :])
        end
    end

    # A and B scalars
    A = (-D[1, 2]*τ_3/τ + D[2, 2] + D[3, 2]*τ_1/τ) / D_0
    B = (D[1, 2]*(τ_3^2 - τ^2)*τ_3/τ + D[3, 2]*(τ^2 - τ_1^2)*τ_1/τ) / 6 / D_0

    # E and F scalars
    E = dot(R_vec[2, :], ρ_vec[2, :])
    F = dot(R_vec[2, :], R_vec[2, :])

    # Lagrange equation coefficients
    a = -(A^2 + 2*A*E + F)
    b = -2*μ_S*B*(A + E)
    c = -(μ_S^2)*(B^2)

    # Solve Lagrange equation
    sol, status = solve_lagrange(a, b, c; niter = niter)

    # Number of solutions
    n_sol = length(sol)

    if iszero(n_sol)

        @warn("""No solutions found for Lagrange equation $(_format_Lagrange_equation(cte(a), cte(b), cte(c)))""")

        return Vector{GaussSolution{T, U}}(undef, 0)

    else

        sol_gauss = Vector{GaussSolution{T, U}}(undef, n_sol)

        for i in eachindex(sol_gauss)
            # Heliocentric range
            r_2 = sol[i]
            # Slant ranges
            ρ = zeros(U, 3)

            num_1 = 6*(D[3, 1]*τ_1/τ_3 + D[2, 1]*τ/τ_3)*(r_2^3) + μ_S*D[3, 1]*(τ^2 - τ_1^2)*τ_1/τ_3
            den_1 = 6*(r_2^3) + μ_S*(τ^2 - τ_3^2)
            ρ[1] = (num_1 / den_1 - D[1, 1]) / D_0

            ρ[2] = A + μ_S*B/(r_2^3)

            num_3 = 6*(D[1, 3]*τ_3/τ_1 - D[2, 3]*τ/τ_1)*(r_2^3) + μ_S*D[1, 3]*(τ^2 - τ_3^2)*τ_3/τ_1
            den_3 = 6*(r_2^3) + μ_S*(τ^2 - τ_1^2)
            ρ[3] = (num_3 / den_3 - D[3, 3]) / D_0

            # Heliocentric position of the NEO
            r_vec = R_vec .+ ρ.*ρ_vec

            # f, g Lagrange coefficients
            f_1 = f_Lagrange(τ_1, r_2)
            f_3 = f_Lagrange(τ_3, r_2)

            g_1 = g_Lagrange(τ_1, r_2)
            g_3 = g_Lagrange(τ_3, r_2)

            # Heliocentric velocity of the NEO
            v_2_vec = (- f_3 * r_vec[1, :] + f_1 * r_vec[3, :]) / (f_1*g_3 - f_3*g_1)

            sol_gauss[i] = GaussSolution{T, U}(vcat(r_vec[2, :], v_2_vec), D, R_vec, ρ_vec, τ_1, τ_3, f_1, g_1, f_3, g_3, status[i])
        end
        # Sort solutions by heliocentric range
        return sort!(sol_gauss, by = x -> norm(x.statevect[1:3]))

    end

end

@doc raw"""
    gauss_norm(dates::Vector{DateTime})

Return a measure of how evenly distributed in time a triplet is; used within [`gauss_triplets`](@ref) to sort triplets 
for Gauss method. The function assumes `dates` is sorted.
"""
gauss_norm(dates::Vector{DateTime}) = abs( (dates[2] - dates[1]).value - (dates[3] - dates[2]).value )

@doc raw"""
    gauss_triplets(dates::Vector{DateTime}, max_triplets::Int = 10, max_iter::Int = 100)

Return a vector of `max_triplets` triplets to be used within [`gaussinitcond`](@ref) to select the best observations for Gauss method. 
The triplets are sorted by [`gauss_norm`](@ref).
"""
function gauss_triplets(dates::Vector{DateTime}, Δ_min::Period, Δ_max::Period, avoid::Vector{Vector{Int}}, max_triplets::Int)
    
    triplets = Vector{Vector{Int}}(undef, 0)
    L = length(dates)
    for i_1 in 1:L-2
        for i_2 in i_1+1:L-1
            for i_3 in i_2+1:L
                tmp = [i_1, i_2, i_3]
                if (Δ_min <= dates[i_3] - dates[i_1] <= Δ_max) && !(tmp in avoid)
                    push!(triplets, tmp)
                end
            end
        end
    end

    sort!(triplets, by = x -> gauss_norm(dates[x]))

    n = min(length(triplets), max_triplets)

    return triplets[1:n]
end

function gauss_triplets(dates::Vector{DateTime}, max_triplets::Int = 10, max_iter::Int = 100)
    Δ_min = Hour(20)
    Δ_max = Day(7)

    triplets = Vector{Vector{Int}}(undef, 0)

    niter = 0

    while length(triplets) < max_triplets && niter < max_iter
        triplets = vcat(triplets, gauss_triplets(dates, Δ_min, Δ_max, triplets, max_triplets))
        if Δ_min >= Hour(1)
            Δ_min -= Hour(1)
        end
        Δ_max += Day(1)
        niter += 1
    end

    n = min(length(triplets), max_triplets)

    return triplets[1:n]
end

function adaptative_maxsteps(radec::Vector{RadecMPC{T}}) where {T <: AbstractFloat}
    # Time difference [ms]
    Δ_ms =  getfield(date(radec[end]) - date(radec[1]), :value)
    # Time difference [days]
    Δ_day = Δ_ms / 86_400_000
    # Adaptative maxsteps
    if Δ_day <= 6
        return ceil(Int, 5*Δ_day)
    else 
        return ceil(Int, Δ_day/2 + 27)
    end
end

@doc raw"""
    gaussinitcond(radec::Vector{RadecMPC{T}}; max_triplets::Int = 10, Q_max::T = 100., niter::Int = 5,
                  varorder::Int = 5, order::Int = order, abstol::T = abstol, parse_eqs::Bool = true) where {T <: AbstractFloat}

Return initial conditions via Gauss method. 

See also [`gauss_method`](@ref).

# Arguments
- `radec::Vector{RadecMPC{T}}`: vector of observations.
- `max_triplets::Int`: maximum number of triplets.
- `Q_max::T`: The maximum nrms that is considered a good enough orbit.
- `niter::Int`: number of iterations for Newton's method.
- `varorder::Int`: order of jet transport perturbation. 
- `order::Int`: order of Taylor polynomials w.r.t. time.
- `abstol::T`: absolute tolerance.
- `parse_eqs::Bool`: whether to use the taylorized method of `RNp1BP_pN_A_J23E_J2S_eph_threads!` or not. 

!!! warning
    This function will set the (global) `TaylorSeries` variables to `δα₁ δα₂ δα₃ δδ₁ δδ₂ δδ₃`. 
"""
function gaussinitcond(radec::Vector{RadecMPC{T}}; max_triplets::Int = 10, Q_max::T = 10., niter::Int = 5,
                       varorder::Int = 5, order::Int = order, abstol::T = abstol, parse_eqs::Bool = true) where {T <: AbstractFloat}

    # Sun's ephemeris
    eph_su = selecteph(sseph, su)
    # Earth's ephemeris
    eph_ea = selecteph(sseph, ea)

    # Eliminate observatories without coordinates 
    filter!(x -> hascoord(x.observatory), radec)
    # Reduce nights by interpolation 
    observatories, dates, α, δ, gdf = reduce_nights(radec)
    # Observations triplets
    triplets = gauss_triplets(dates, max_triplets)

    # Julian day of first (last) observation
    t0, tf = datetime2julian(radec[1].date), datetime2julian(radec[end].date)
    # Julian day when to start propagation
    jd0 = zero(T)
    # Maximum number of steps
    maxsteps = adaptative_maxsteps(radec)

    # Jet transport perturbation (ra/dec)
    dq = scaled_variables("δα₁ δα₂ δα₃ δδ₁ δδ₂ δδ₃"; order = varorder)

    # Normalized root mean square error (NRMS)
    best_Q = T(Inf)
    # Initial conditions
    best_sol = zero(NEOSolution{T, T})
    # Break flag
    flag = false

    # Iterate over triplets
    for j in eachindex(triplets)
        
        # Current triplet
        triplet = triplets[j]

        # Julian day when to start propagation
        jd0 = datetime2julian(dates[triplet[2]])
        # Number of years in forward integration 
        nyears_fwd = (tf - jd0 + 2) / yr
        # Number of years in backward integration
        nyears_bwd = -(jd0 - t0 + 2) / yr

        # Gauss method solution 
        sol = gauss_method(observatories[triplet], dates[triplet], α[triplet] .+ dq[1:3], δ[triplet] .+ dq[4:6]; niter = niter)

        # Filter Gauss solutions by heliocentric energy
        filter!(x -> heliocentric_energy(x.statevect) <= 0, sol)

        # Iterate over Gauss solutions
        for i in eachindex(sol)

            # Initial conditions (jet transport)
            q0 = sol[i].statevect .+ eph_su(jd0 - PE.J2000)

            # Propagation 
            bwd = propagate(RNp1BP_pN_A_J23E_J2S_eph_threads!, maxsteps, jd0, nyears_bwd, q0; 
                            order = order, abstol = abstol, parse_eqs = parse_eqs)
            if length(bwd.t) == maxsteps+1
                continue
            end
            
            fwd = propagate(RNp1BP_pN_A_J23E_J2S_eph_threads!, maxsteps, jd0, nyears_fwd, q0; 
                            order = order, abstol = abstol, parse_eqs = parse_eqs)

            if length(fwd.t) == maxsteps + 1
                continue
            end

            # O-C residuals
            res = residuals(radec; xvs = et -> auday2kmsec(eph_su(et/daysec)), xve = et -> auday2kmsec(eph_ea(et/daysec)), 
                            xva = et -> bwdfwdeph(et, bwd, fwd))

            # Subset of radec for orbit fit
            idxs = findall(x -> triplet[1] <= x <= triplet[3], gdf.groups)
            sort!(idxs)

            # Orbit fit
            fit = tryls(res[idxs], zeros(get_numvars()), niter)

            !fit.success && continue

            # Right iteration
            for k in triplet[3]+1:length(gdf)
                extra = findall(x -> x == k, gdf.groups)
                fit_new = tryls(res[idxs ∪ extra], zeros(get_numvars()), niter)
                if fit_new.success
                    fit = fit_new
                    idxs = vcat(idxs, extra)
                    sort!(idxs)
                else 
                    break
                end
            end

            # Left iteration
            for k in triplet[1]-1:-1:1
                extra = findall(x -> x == k, gdf.groups)
                fit_new = tryls(res[idxs ∪ extra], zeros(get_numvars()), niter)
                if fit_new.success
                    fit = fit_new
                    idxs = vcat(idxs, extra)
                    sort!(idxs)
                else 
                    break
                end
            end

            # NRMS
            Q = nrms(res, fit)

            # TO DO: check cases where newton converges but diffcorr no
            # Update NRMS and initial conditions
            if Q < best_Q
                best_Q = Q
                best_sol = evalfit(NEOSolution(bwd, fwd, res, fit))
            end 
            # Break condition
            if Q <= Q_max
                flag = true
                break
            end
        end
        if flag
            break
        end
    end 

    # Case: all solutions were unsuccesful
    if isinf(best_Q)
        return zero(NEOSolution{T, T})
    # Case: at least one solution was succesful
    else 
        return best_sol
    end
end 