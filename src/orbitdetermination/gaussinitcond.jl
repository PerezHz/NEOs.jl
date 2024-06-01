@doc raw"""
    GaussSolution{T <: Real, U <: Number}

A preliminary orbit obtained from Gauss method of orbit determination.

See also [`gauss_method`](@ref).

# Fields

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
    function GaussSolution{T, U}(
        statevect::Vector{U}, ρ::Vector{U}, D::Matrix{U}, R_vec::Matrix{T},
        ρ_vec::Matrix{U}, τ_1::T, τ_3::T, f_1::U, g_1::U, f_3::U, g_3::U) where {T <: Real, U <: Number}
        return new{T, U}(statevect, ρ, D, R_vec, ρ_vec, τ_1, τ_3, f_1, g_1, f_3, g_3)
    end
end

# Outer constructor
function GaussSolution(
    statevect::Vector{U}, ρ::Vector{U}, D::Matrix{U}, R_vec::Matrix{T},
    ρ_vec::Matrix{U}, τ_1::T, τ_3::T, f_1::U, g_1::U, f_3::U, g_3::U) where {T <: Real, U <: Number}
    return GaussSolution{T, U}(statevect, ρ, D, R_vec, ρ_vec, τ_1, τ_3, f_1, g_1, f_3, g_3)
end

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

    return join(["r⁸ ", a_sgn, " ", abs(a), " r⁶ ", b_sgn, " ", abs(b), " r³ ", c_sgn, " ", abs(c), " = 0"])
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

# Arguments

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
    sol = solve_lagrange(a, b, c; niter = params.newtoniter)

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
    numberofdays(dates::Vector{DateTime})
    numberofdays(dates::Vector{RadecMPC{T}}) where {T <: Real}
    numberofdays(dates::Vector{Tracklet{T}}) where {T <: Real}

Return the timespan of `dates` in days. The function assumes `dates` is sorted.
"""
numberofdays(dates::Vector{DateTime}) = (dates[end] - dates[1]).value / 86_400_000

function numberofdays(dates::Vector{RadecMPC{T}}) where {T <: Real}
    return (dates[end].date - dates[1].date).value / 86_400_000
end

function numberofdays(dates::Vector{Tracklet{T}}) where {T <: Real}
    return (dates[end].radec[end].date - dates[1].radec[1].date).value / 86_400_000
end

function _gausstriplets1(tracklets::Vector{Tracklet{T}}, maxtriplets::Int) where {T <: Real}
    # Number of tracklets
    L = length(tracklets)
    # All possible triplets of indices
    CI = CartesianIndices((1:L-2, 2:L-1, 3:L))
    # Allocate memory
    triplets = zeros(Int, 3, maxtriplets)
    τ = fill(T(Inf), maxtriplets)
    # Check all possible triplets of indices
    for ci in CI
        # Indices cannot repeat
        !allunique(ci.I) && continue
        # Unfold indices
        i, j, k = ci.I
        # Timespan must be at least one day
        Δ = (tracklets[k].date - tracklets[i].date).value
        Δ < 86_400_000 && continue
        # Absolute difference between τ1 and τ3
        τ1 = (tracklets[j].date - tracklets[i].date).value
        τ3 = (tracklets[k].date - tracklets[j].date).value
        _τ_ = abs(τ3 - τ1)
        # Current max τ
        τmax, n = findmax(τ)
        # Update triplets and τ
        if _τ_ < τmax
            triplets[:, n] .= i, j, k
            τ[n] = _τ_
        end
    end
    # Remove Infs and sort
    idxs = findall(!isinf, τ)
    perm = sortperm(view(τ, idxs))

    return triplets[:, perm], τ[perm]
end

function _gausstriplets2!(observatories::Vector{ObservatoryMPC{T}}, dates::Vector{DateTime}, α::Vector{T},
                          δ::Vector{T}, tracklets::AbstractVector{Tracklet{T}}) where {T <: Real}
    # Unfold tracklets
    a, b, c = tracklets
    # All possible triplets of indices
    CI = CartesianIndices((a.nobs, b.nobs, c.nobs))
    # Allocate memory
    triplet = [0, 0, 0]
    τ = Inf
    # Check all possible triplets of indices
    for ci in CI
        # Unfold indices
        i, j, k = ci.I
        # Absolute difference between τ1 and τ3
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

    return nothing
end

function _adam!(q::Vector{TaylorN{T}}, jd0::T, tracklet::Tracklet, params::NEOParameters{T};
                dynamics::D = newtonian!) where {T <: Real, D}
    # Exploratory propagation
    bwd, fwd, _ = propres(tracklet.radec, jd0, q(), params)
    # Admissible region
    A = AdmissibleRegion(tracklet, params)
    # Epoch [days since J2000]
    At0 = datetime2days(A.date)
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
    ae, _ = adam(tracklet.radec, A, ρ, v_ρ, params; scale = :log, dynamics = dynamics)
    # Epoch [julian days] (corrected for light-time)
    jd0 = datetime2julian(A.date) - ae[5] / c_au_per_day
    # Convert attributable elements to barycentric cartesian coordinates
    q0 = attr2bary(A, ae, params)
    # Jet Transport initial condition
    q .= [q0[i] + (q0[i] / 10^5) * TaylorN(i, order = params.tsaorder) for i in 1:6]

    return jd0
end

@doc raw"""
    gaussinitcond(radec::Vector{RadecMPC{T}}, tracklets::Vector{Tracklet{T}},
                  params::NEOParameters{T}; dynamics::D = newtonian!) where {T <: Real, D}

Compute an orbit via Jet Transport Gauss Method.

See also [`gauss_method`](@ref).

# Arguments

- `radec::Vector{RadecMPC{T}}`: vector of optical astrometry.
- `tracklets::Vector{Tracklet{T}},`: vector of tracklets.
- `params::NEOParameters{T}`: see `Gauss Method Parameters` of [`NEOParameters`](@ref).
- `dynamics::D`: dynamical model.
"""
function gaussinitcond(radec::Vector{RadecMPC{T}}, tracklets::Vector{Tracklet{T}},
                       params::NEOParameters{T}; dynamics::D = newtonian!) where {T <: Real, D}
    # gauss_method input
    observatories = Vector{ObservatoryMPC{T}}(undef, 3)
    dates = Vector{DateTime}(undef, 3)
    α = Vector{T}(undef, 3)
    δ = Vector{T}(undef, 3)
    # Observations triplets
    triplets, _ = _gausstriplets1(tracklets, params.max_triplets)
    # Jet transport scaling factors
    scalings = fill(1e-6, 6)
    # Jet transport perturbation (ra/dec)
    dq = [scalings[i] * TaylorN(i, order = params.gaussorder) for i in 1:6]
    # Best orbit
    best_sol = zero(NEOSolution{T, T})
    best_Q = T(Inf)
    # Convergence flag
    flag = false
    # Iterate over triplets
    for triplet in eachcol(triplets)
        # Find best triplet of observations
        _gausstriplets2!(observatories, dates, α, δ, view(tracklets, triplet))
        # Julian day of middle observation
        _jd0_ = datetime2julian(dates[2])
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
            # ADAM (if requested by the user)
            if params.adamhelp
                jd0 = _adam!(q, jd0, tracklets[triplet[2]], params; dynamics)
            end
            # Jet transport least squares
            sol = jtls(radec, tracklets, jd0, q, triplet[2], params; dynamics)
            # NRMS
            Q = nrms(sol)
            # Update best orbit
            if Q < best_Q
                best_sol = sol
                best_Q = Q
                # Break condition
                if Q <= params.gaussQmax
                    flag = true
                    break
                end
            end
        end
        flag && break
    end

    return best_sol
end