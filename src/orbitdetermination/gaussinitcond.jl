@doc raw"""
    GaussSolution{T <: Real, U <: Number}

A preliminary orbit obtained from Gauss method of orbit determination.

See also [`gauss_method`](@ref).

## Fields

- `statevect::Vector{U}`: state vector at middle observation.
- `ρ::Vector{U}`: topocentric ranges.
- `D::Matrix{U}`: D matrix.
- `R_vec::Matrix{T}`: observer's heliocentric positions.
- `ρ_vec::Matrix{U}`: line of sight unit vectors.
- `τ_1::T`: time between first and second observations.
- `τ_3::T`: time between third and second observations.
- `f_1, g_1, f_3, g_3::U`: Lagrange coefficients.

!!! reference
    See Algorithm 5.5 in page 274 of https://doi.org/10.1016/C2016-0-02107-1.
"""
@auto_hash_equals struct GaussSolution{T <: Real, U <: Number}
    statevect::Vector{U}
    ρ::Vector{U}
    D::Matrix{U}
    R_vec::Matrix{T}
    ρ_vec::Matrix{U}
    τ_1::T
    τ_3::T
    f_1::U
    g_1::U
    f_3::U
    g_3::U
    # Inner constructor
    function GaussSolution{T, U}(statevect::Vector{U}, ρ::Vector{U}, D::Matrix{U},
        R_vec::Matrix{T}, ρ_vec::Matrix{U}, τ_1::T, τ_3::T, f_1::U, g_1::U, f_3::U,
        g_3::U) where {T <: Real, U <: Number}
        return new{T, U}(statevect, ρ, D, R_vec, ρ_vec, τ_1, τ_3, f_1, g_1, f_3, g_3)
    end
end

# Outer constructor
GaussSolution(statevect::Vector{U}, ρ::Vector{U}, D::Matrix{U},
    R_vec::Matrix{T}, ρ_vec::Matrix{U}, τ_1::T, τ_3::T, f_1::U,
    g_1::U, f_3::U, g_3::U) where {T <: Real, U <: Number} =
    GaussSolution{T, U}(statevect, ρ, D, R_vec, ρ_vec, τ_1, τ_3, f_1, g_1, f_3, g_3)

# Print method for GaussSolution
# Examples:
# Gauss solution (r = 1.0800950907383229 AU)
function show(io::IO, g::GaussSolution{T, U}) where {T <: Real, U <: Number}
    print(io, "Gauss solution (r = ", norm(cte.(g.statevect[1:3])), " AU)")
end

@doc raw"""
    topounit(α::T, δ::T) where {T <: Number}
    topounit(obs::RadecMPC{T}) where {T <: Real}

Return the topocentric unit vector.
"""
function topounit(α::T, δ::T) where {T <: Number}

    sin_α, cos_α = sincos(α)
    sin_δ, cos_δ = sincos(δ)

    pos = [cos_δ * cos_α, cos_δ * sin_α, sin_δ]

    return pos
end

topounit(obs::RadecMPC{T}) where {T <: Real} = topounit(obs.α, obs.δ)

# TO DO: Should we allow to use other μ?

@doc raw"""
    f_Lagrange(τ::T, r::U) where {T <: Real, U <: Number}

Return the 1st order approximation to Lagrange's f function.
"""
function f_Lagrange(τ::T, r::U) where {T <: Real, U <: Number}
    r3 = r * r * r
    return 1 - μ_S * (τ^2) / 2 / r3
end

@doc raw"""
    g_Lagrange(τ::T, r::U) where {T <: Real, U <: Number}

Return the 1st order approximation to Lagrange's g function.
"""
function g_Lagrange(τ::T, r::U) where {T <: Real, U <: Number}
    r3 = r * r * r
    return τ - μ_S * (τ^3) / 6 / r3
end

@doc raw"""
    vectors2matrix(x::T) where {T <: AbstractVector}

Convert a vector of vectors `x` to a matrix.
"""
vectors2matrix(x::T) where {T <: AbstractVector} = permutedims(reduce(hcat, x))

@doc raw"""
    _format_Lagrange_equation(a::T, b::T, c::T) where {T <: Real}

Format Lagrange equation as `r⁸ + a r⁶ + b r³ + c = 0`.
"""
function _format_Lagrange_equation(a::T, b::T, c::T) where {T <: Real}
    a_sgn = a ≥ 0 ? "+" : "-"
    b_sgn = b ≥ 0 ? "+" : "-"
    c_sgn = c ≥ 0 ? "+" : "-"

    return join(["r⁸ ", a_sgn, " ", abs(a), " r⁶ ", b_sgn, " ", abs(b), " r³ ", c_sgn,
        " ", abs(c), " = 0"])
end

@doc raw"""
    lagrange(x::T, a::U, b::U, c::U) where {T, U <: Number}

Lagrange polynomial to be solved during Gauss method.
"""
function lagrange(x::T, a::U, b::U, c::U) where {T, U <: Number}
    # Evaluate via Horner's method
    x2 = x * x
    x3 = x2 * x
    return c + x3 * (b + x3 * (a + x2))
end

@doc raw"""
    lagrange_derivative(x::T, a::T, b::T) where {T <: Number}

Derivative of Lagrange polynomial to be solved during Gauss method.
"""
function lagrange_derivative(x::T, a::U, b::U) where {T, U <: Number}
    # Evaluate via Horner's method
    x2 = x * x
    x3 = x2 * x
    return x2 * (3*b + x3 * (6*a + 8*x2))
end

# TO DO: Allow to control interval over which to look for solutions
# Currently we look between the radius of the Sun (∼0.00465047 AU) and
#  the radius of the Solar System (∼40 AU)

@doc raw"""
    solve_lagrange(a::T, b::T, c::T; niter::Int = 5) where {T <: Real}
    solve_lagrange(a::TaylorN{T}, b::TaylorN{T}, c::TaylorN{T}; niter::Int = 5) where {T <: Real}

Solve Lagrange polynomial.

See also [`lagrange`](@ref).
"""
function solve_lagrange(a::T, b::T, c::T; niter::Int = 5, rmin = 0.00465047, rmax = 40.0) where {T <: Real}
    return find_zeros(x -> lagrange(x, a, b, c), rmin, rmax)
end

function solve_lagrange(a::TaylorN{T}, b::TaylorN{T}, c::TaylorN{T}; niter::Int = 5) where {T <: Real}
    # 0-th order solution
    sol_0 = solve_lagrange(cte(a), cte(b), cte(c))
    # Vector of solutions
    sol = Vector{TaylorN{T}}(undef, length(sol_0))
    # Conversion factor to TaylorN
    oneN = one(a)
    # Iterate solutions
    for i in eachindex(sol)
        # Newton's method
        r_0::TaylorN{T} = sol_0[i] * oneN
        r_2::TaylorN{T} = sol_0[i] * oneN
        for _ in 1:niter
            r_2 = r_0 - lagrange(r_0, a, b, c) / lagrange_derivative(r_0, a, b)
            r_0 = r_2
        end
        sol[i] = r_2
    end

    return sol
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
    gauss_method(obs::Vector{RadecMPC{T}}, params::NEOParameters{T}) where {T <: Real}
    gauss_method(observatories::Vector{ObservatoryMPC{T}}, dates::Vector{DateTime},
        α::Vector{U}, δ::Vector{U}, params::NEOParameters{T}) where {T <: Real, U <: Number}

Core Gauss method of Initial Orbit determination (IOD).

## Arguments

- `obs::Vector{RadecMPC{T}}`: three observations.
- `observatories::Vector{ObservatoryMPC{T}}`: sites of observation.
- `dates::Vector{DateTime}`: times of observation.
- `α::Vector{U}`: right ascension [rad].
- `δ::Vector{U}`: declination [rad].
- `params::NEOParameters{T}`: see `Gauss Method Parameters` of [`NEOParameters`](@ref).

!!! reference
    See Algorithm 5.5 in page 274 https://doi.org/10.1016/C2016-0-02107-1.
"""
function gauss_method(obs::Vector{RadecMPC{T}}, params::NEOParameters{T}) where {T <: Real}

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

    return gauss_method(observatories, dates, α, δ, params)
end

function gauss_method(observatories::Vector{ObservatoryMPC{T}}, dates::Vector{DateTime},
                      α::Vector{U}, δ::Vector{U}, params::NEOParameters{T}) where {T <: Real, U <: Number}

    # Check we have exactly three observations
    @assert length(observatories) == length(dates) == length(α) == length(δ) == 3 "Gauss method requires exactly three observations"
    # Check observations are in temporal order
    @assert issorted(dates) "Observations must be in temporal order"

    # Times of observation [et]
    t_et = dtutc2et.(dates)
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
    # Heliocentric state vector of the Earth [au, au/day]
    G_vec = params.eph_ea.(t_days) - params.eph_su.(t_days)
    # Observer's heliocentric positions [au, au/day]
    R_vec = vectors2matrix(G_vec .+  g_vec)[:, 1:3]

    # Cross products
    p_vec = zeros(U, 3, 3)
    p_vec[1, :] = @views cross(ρ_vec[2, :], ρ_vec[3, :])
    p_vec[2, :] = @views cross(ρ_vec[1, :], ρ_vec[3, :])
    p_vec[3, :] = @views cross(ρ_vec[1, :], ρ_vec[2, :])

    # Gauss scalar
    D_0 = dot(ρ_vec[1, :], p_vec[1, :])

    # Matrix of triple products
    D = zeros(U, 3, 3)
    for i in 1:3
        for j in 1:3
            D[i, j] = @views dot(R_vec[i, :], p_vec[j, :])
        end
    end

    # A and B scalars
    A = (-D[1, 2]*τ_3/τ + D[2, 2] + D[3, 2]*τ_1/τ) / D_0
    B = (D[1, 2]*(τ_3^2 - τ^2)*τ_3/τ + D[3, 2]*(τ^2 - τ_1^2)*τ_1/τ) / 6 / D_0

    # E and F scalars
    E = @views dot(R_vec[2, :], ρ_vec[2, :])
    F = @views dot(R_vec[2, :], R_vec[2, :])

    # Lagrange equation coefficients
    a = -(A*A + 2*A*E + F)
    b = -2*μ_S*B*(A + E)
    c = -(μ_S^2)*(B*B)

    # Solve Lagrange equation
    sol = solve_lagrange(a, b, c; niter = params.lsiter)

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
            # Range cubed
            r_23 = r_2 * r_2 * r_2
            # Slant ranges
            ρ = zeros(U, 3)

            num_1 = 6*(D[3, 1]*τ_1/τ_3 + D[2, 1]*τ/τ_3)*r_23 + μ_S*D[3, 1]*(τ^2 - τ_1^2)*τ_1/τ_3
            den_1 = 6*r_23 + μ_S*(τ^2 - τ_3^2)
            ρ[1] = (num_1 / den_1 - D[1, 1]) / D_0

            ρ[2] = A + μ_S*B/r_23

            num_3 = 6*(D[1, 3]*τ_3/τ_1 - D[2, 3]*τ/τ_1)*r_23 + μ_S*D[1, 3]*(τ^2 - τ_3^2)*τ_3/τ_1
            den_3 = 6*r_23 + μ_S*(τ^2 - τ_1^2)
            ρ[3] = (num_3 / den_3 - D[3, 3]) / D_0

            # Heliocentric position of the NEO
            r_vec = R_vec .+ ρ .* ρ_vec

            # f, g Lagrange coefficients
            f_1 = f_Lagrange(τ_1, r_2)
            f_3 = f_Lagrange(τ_3, r_2)

            g_1 = g_Lagrange(τ_1, r_2)
            g_3 = g_Lagrange(τ_3, r_2)

            # Heliocentric velocity of the NEO
            v_2_vec = @views (- f_3 * r_vec[1, :] + f_1 * r_vec[3, :]) / (f_1*g_3 - f_3*g_1)

            sol_gauss[i] = GaussSolution{T, U}(
                vcat(r_vec[2, :], v_2_vec), ρ, D, R_vec,
                ρ_vec, τ_1, τ_3, f_1, g_1, f_3, g_3
            )
        end
        # Sort solutions by heliocentric range
        return sort!(sol_gauss, by = x -> norm(cte.(x.statevect[1:3])))

    end

end

@doc raw"""
    numberofdays(x)

Return the timespan of `x::AbstractVector{T}` in days, where
`T` can be `DateTime`, `RadecMPC` or `Tracklet`.
"""
function numberofdays(dates::AbstractVector{DateTime})
    t0, tf = extrema(dates)
    return (tf - t0).value / 86_400_000
end

function numberofdays(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
    t0, tf = extrema(date, radec)
    return (tf - t0).value / 86_400_000
end

function numberofdays(tracklets::AbstractVector{Tracklet{T}}) where {T <: Real}
    dates = map(t -> extrema(date, t.radec), tracklets)
    t0, tf = minimum(first, dates), maximum(last, dates)
    return (tf - t0).value / 86_400_000
end

# Update `observatories`, `dates`, `α` and `δ` with the best combination of
# three observations for Gauss' method. See the line after equation (27) of
# https://doi.org/10.1016/j.icarus.2007.11.033
function _gausstriplet!(observatories::Vector{ObservatoryMPC{T}},
    dates::Vector{DateTime}, α::Vector{T}, δ::Vector{T},
    tracklets::AbstractVector{Tracklet{T}}) where {T <: Real}
    # Unfold tracklets
    a, b, c = tracklets
    # All possible triplets of indices
    CI = CartesianIndices((a.nobs, b.nobs, c.nobs))
    # Allocate memory
    triplet = [0, 0, 0]
    τ = typemax(Int)
    # Check all possible triplets of indices
    for ci in CI
        # Unfold indices
        i, j, k = ci.I
        # Minimize asymmetry
        τ1 = (b.radec[j].date - a.radec[i].date).value
        τ3 = (c.radec[k].date - b.radec[j].date).value
        _τ_ = abs(τ3 - τ1)
        # Update triplet and τ
        if _τ_ < τ
            triplet .= i, j, k
            τ = _τ_
        end
    end
    # Unfold triplet
    i, j, k = triplet
    # Update observatories, dates, α and δ
    observatories .= a.observatory, b.observatory, c.observatory
    dates .= a.radec[i].date, b.radec[j].date, c.radec[k].date
    α .= a.radec[i].α, b.radec[j].α, c.radec[k].α
    δ .= a.radec[i].δ, b.radec[j].δ, c.radec[k].δ
    # Sort triplet
    idxs = sortperm(dates)
    permute!(observatories, idxs)
    permute!(dates, idxs)
    permute!(α, idxs)
    permute!(δ, idxs)

    return nothing
end

function _gausstriplet!(observatories::Vector{ObservatoryMPC{T}},
    dates::Vector{DateTime}, α::Vector{T}, δ::Vector{T},
    radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
    # Allocate memory
    triplet = [0, 0, 0]
    τa, τb = 0, typemax(Int)
    # Check all possible triplets of indices
    L = length(radec)
    for i in 1:L-2, j in i+1:L-1, k in j+1:L
        # Maximize timespan and minimize asymmetry
        τ1 = (radec[j].date - radec[i].date).value
        τ3 = (radec[k].date - radec[j].date).value
        _τa_, _τb_ = abs(τ1) + abs(τ3), abs(τ3 - τ1)
        # Update triplet and τ
        if _τa_ > τa || (_τa_ == τa && _τb_ < τb)
            triplet .= i, j, k
            τa, τb = _τa_, _τb_
        end
    end
    # Unfold triplet
    i, j, k = triplet
    # Update observatories, dates, α and δ
    @. observatories = observatory(radec[triplet])
    @. dates = date(radec[triplet])
    @. α = ra(radec[triplet])
    @. δ = dec(radec[triplet])
    # Sort triplet
    idxs = sortperm(dates)
    permute!(observatories, idxs)
    permute!(dates, idxs)
    permute!(α, idxs)
    permute!(δ, idxs)

    return nothing
end

# Update an initial condition `q`, obtained by Gauss' method, via ADAM
# minimization over the middle tracklet's manifold of variations.
function _adam!(od::ODProblem{D, T}, i::Int, q::Vector{TaylorN{T}}, jd0::T,
    params::NEOParameters{T}) where {D, T <: Real}
    # Middle tracklet
    tracklet = od.tracklets[i]
    # Exploratory propagation
    bwd, fwd, _ = propres(od, jd0, q(), params; idxs = indices(tracklet))
    # Admissible region
    A = AdmissibleRegion(tracklet, params)
    # Epoch [days since J2000]
    At0 = dtutc2days(A.date)
    # Barycentric cartesian initial condition
    if At0 <= jd0 - JD_J2000
        q0 = bwd(At0)
    else
        q0 = fwd(At0)
    end
    # Range and range rate
    ρ, v_ρ = bary2topo(A, q0)
    # Boundary projection
    ρ, v_ρ = boundary_projection(A, ρ, v_ρ)
    # ADAM
    aes, Qs = adam(od, i, A, ρ, v_ρ, params; scale = :log)
    ae, Q = aes[:, end], Qs[end]
    # Epoch [julian days] (corrected for light-time)
    jd0 = dtutc2jdtdb(A.date) - ae[5] / c_au_per_day
    # Convert attributable elements to barycentric cartesian coordinates
    q0 = attr2bary(A, ae, params)
    # Jet Transport initial condition
    q .= [q0[j] + (abs(q0[j]) / 10^5) * TaylorN(j, order = params.tsaorder) for j in 1:6]

    return jd0
end

@doc raw"""
    gaussiod(od::ODProblem{D, T}, params::NEOParameters{T}) where {D, T <: Real}

Fit a preliminary orbit to `od` via jet transport Gauss method.

See also [`gauss_method`](@ref).

## Arguments

- `od::ODProblem{D, T}`: an orbit determination problem.
- `params::NEOParameters{T}`: see `Gauss' Method Parameters` of [`NEOParameters`](@ref).
"""
function gaussiod(od::ODProblem{D, T}, params::NEOParameters{T}) where {D, T <: Real}
    # Allocate memory for orbit
    sol = zero(NEOSolution{T, T})
    # Unfold parameters
    safegauss, varorder, significance = params.safegauss, params.gaussorder,
        params.significance
    # This function requires exactly 3 tracklets
    (safegauss && length(od.tracklets) != 3) && return sol
    # Jet transport scaling factors
    scalings = fill(1e-6, 6)
    # Jet transport perturbation (ra/dec)
    dq = [scalings[i] * TaylorN(i, order = varorder) for i in 1:6]
    # gauss_method input
    observatories = Vector{ObservatoryMPC{T}}(undef, 3)
    dates = Vector{DateTime}(undef, 3)
    α = Vector{T}(undef, 3)
    δ = Vector{T}(undef, 3)
    # Find best triplet of observations
    if safegauss
        _gausstriplet!(observatories, dates, α, δ, od.tracklets)
    else
        _gausstriplet!(observatories, dates, α, δ, od.radec)
    end
    # Julian day of middle observation
    _jd0_ = dtutc2jdtdb(dates[2])
    # Gauss method solution
    solG = gauss_method(observatories, dates, α .+ dq[1:3], δ .+ dq[4:6], params)
    # Filter non-physical (negative) rho solutions
    filter!(x -> all(cte.(x.ρ) .> 0), solG)
    # Filter Gauss solutions by heliocentric energy
    filter!(x -> cte(cte(heliocentric_energy(x.statevect))) <= 0, solG)
    # Iterate over Gauss solutions
    for i in eachindex(solG)
        # Epoch (corrected for light-time)
        jd0 = _jd0_ - cte(solG[i].ρ[2]) / c_au_per_day
        # Jet transport initial conditions
        q = solG[i].statevect .+ params.eph_su(jd0 - JD_J2000)
        # Jet Transport Least Squares
        _sol_ = jtls(od, jd0, q, od.tracklets, params, true)
        # Update solution
        sol = updatesol(sol, _sol_, od.radec)
        # Termination condition
        critical_value(sol) < significance && return sol
        # ADAM help
        j = safegauss ? 2 : findfirst(@. dates[2] in od.tracklets)
        nobs(od.tracklets[j]) < 2 && continue
        jd0 = _adam!(od, j, q, jd0, params)
        # Jet Transport Least Squares
        _sol_ = jtls(od, jd0, q, od.tracklets, params, true)
        # Update solution
        sol = updatesol(sol, _sol_, od.radec)
        # Termination condition
        critical_value(sol) < significance && return sol
    end

    return sol
end