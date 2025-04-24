# Topocentric (line-of-sight) unit vector
function topounit(α::T, δ::T) where {T <: Number}

    sin_α, cos_α = sincos(α)
    sin_δ, cos_δ = sincos(δ)

    pos = [cos_δ * cos_α, cos_δ * sin_α, sin_δ]

    return pos
end

topounit(obs::RadecMPC{T}) where {T <: Real} = topounit(obs.α, obs.δ)

# 1st order approximation to Lagrange's f and g functions
# TO DO: Should we allow to use other μ?
function f_Lagrange(τ::T, r::U) where {T <: Real, U <: Number}
    r3 = r * r * r
    return 1 - μ_S * (τ^2) / 2 / r3
end

function g_Lagrange(τ::T, r::U) where {T <: Real, U <: Number}
    r3 = r * r * r
    return τ - μ_S * (τ^3) / 6 / r3
end

# Format Lagrange equation as `r⁸ + a r⁶ + b r³ + c = 0`
function _format_Lagrange_equation(a::T, b::T, c::T) where {T <: Real}
    a_sgn = a ≥ 0 ? "+" : "-"
    b_sgn = b ≥ 0 ? "+" : "-"
    c_sgn = c ≥ 0 ? "+" : "-"

    return join(["r⁸ ", a_sgn, " ", abs(a), " r⁶ ", b_sgn, " ", abs(b), " r³ ",
        c_sgn, " ", abs(c), " = 0"])
end

# Lagrange polynomial to be solved during Gauss method and its derivative
function lagrange(x::T, a::U, b::U, c::U) where {T, U <: Number}
    # Evaluate via Horner's method
    x2 = x * x
    x3 = x2 * x
    return c + x3 * (b + x3 * (a + x2))
end

function lagrange_derivative(x::T, a::U, b::U) where {T, U <: Number}
    # Evaluate via Horner's method
    x2 = x * x
    x3 = x2 * x
    return x2 * (3*b + x3 * (6*a + 8*x2))
end

# Solve Lagrange polynomial
# TO DO: Allow to control interval over which to look for solutions
# Currently we look between the radius of the Sun (∼0.00465047 AU) and
# the radius of the Solar System (∼40 AU)
solve_lagrange(a::T, b::T, c::T; niter::Int = 5, rmin = 0.00465047,
    rmax = 40.0) where {T <: Real} = find_zeros(x -> lagrange(x, a, b, c), rmin, rmax)

function solve_lagrange(a::TaylorN{T}, b::TaylorN{T}, c::TaylorN{T};
    niter::Int = 5) where {T <: Real}
    # 0-th order solution
    sol0 = solve_lagrange(cte(a), cte(b), cte(c))
    # Vector of solutions
    sol = Vector{TaylorN{T}}(undef, length(sol0))
    # Conversion factor to TaylorN
    oneN = one(a)
    # Iterate solutions
    for i in eachindex(sol)
        # Newton's method
        r0::TaylorN{T} = sol0[i] * oneN
        r2::TaylorN{T} = sol0[i] * oneN
        for _ in 1:niter
            r2 = r0 - lagrange(r0, a, b, c) / lagrange_derivative(r0, a, b)
            r0 = r2
        end
        sol[i] = r2
    end

    return sol
end

# Check topocentric slant ranges are positive
isphysical(orbit::GaussOrbit) = all(orbit.ρ .> 0)
# Check heliocentric energy is negative
function isclosed(orbit::GaussOrbit)
    # Heliocentric state vector [au, au/day]
    r = orbit.r_vec[:, 2]
    # Heliocentric energy per unit mass
    kinetic = 0.5 * (r[4]^2 + r[5]^2 + r[6]^2)
    potential = k_gauss^2 / sqrt(r[1]^2 + r[2]^2 + r[3]^2)
    E  = kinetic - potential

    return E <= 0
end

# Find the tracklet whose epoch is closest to t
closest_tracklet(t::T, tracklets::Vector{Tracklet{T}}) where {T <: Real} =
    findmin(@. abs(t - dtutc2days(date(tracklets))))[2]

# Find the best combination of three observations for Gauss' method.
# See the line after equation (27) of https://doi.org/10.1016/j.icarus.2007.11.033
function gausstriplet(tracklets::AbstractVector{Tracklet{T}}) where {T <: Real}
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
    observatories = [a.observatory, b.observatory, c.observatory]
    dates = [a.radec[i].date, b.radec[j].date, c.radec[k].date]
    α = [a.radec[i].α, b.radec[j].α, c.radec[k].α]
    δ = [a.radec[i].δ, b.radec[j].δ, c.radec[k].δ]
    # Sort triplet
    idxs = sortperm(dates)
    permute!(observatories, idxs)
    permute!(dates, idxs)
    permute!(α, idxs)
    permute!(δ, idxs)

    return observatories, dates, α, δ
end

function gausstriplet(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
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
    observatories = observatory.(radec[triplet])
    dates = date.(radec[triplet])
    α = ra.(radec[triplet])
    δ = dec.(radec[triplet])
    # Sort triplet
    idxs = sortperm(dates)
    permute!(observatories, idxs)
    permute!(dates, idxs)
    permute!(α, idxs)
    permute!(δ, idxs)

    return observatories, dates, α, δ
end

@doc raw"""
    gaussmethod(od, params) where {D, T <: Real}

Gauss method for preliminary orbit determination.

## Arguments

- `od::ODProblem{D, T}`: orbit determination problem.
- `params::Parameters{T}`: see the `Gauss Method` section of [`Parameters`](@ref).

!!! reference
    See Section 3 of https://doi.org/10.1007/s10569-025-10246-2.
"""
function gaussmethod(od::ODProblem{D, T}, params::Parameters{T}) where {D, T <: Real}
    # Unpack
    @unpack safegauss, gaussorder, eph_su = params
    @unpack tracklets, radec = od
    # Find best triplet of observations
    observatories, dates, α, δ = safegauss ? gausstriplet(tracklets) : gausstriplet(radec)
    # Julian day of middle observation
    _jd0_ = dtutc2jdtdb(dates[2])
    # Scaling factors
    scalings = fill(π / (180 * 3_600), 6)
    # Jet transport perturbation (ra/dec)
    dq = [scalings[i] * TaylorN(i, order = gaussorder) for i in 1:6]
    # Gauss method solutions
    τ_1, τ_3, ρ_vec, R_vec, D_0, _D_, a, b, c, r_2s, r_vec, ρ =
        gaussmethod(observatories, dates, α .+ dq[1:3], δ .+ dq[4:6], params)
    # Preliminary orbits
    orbits = Vector{GaussOrbit{T, T}}(undef, length(r_2s))
    for i in eachindex(orbits)
        # Epoch (corrected for light-time)
        jd0 = _jd0_ - cte(ρ[2, i]) / c_au_per_day
        # Jet transport initial condition
        q0 = r_vec[:, 2, i] + eph_su(jd0 - JD_J2000)
        # Propagation and residuals
        bwd, fwd, res = propres(od, jd0, q0, params)
        if isempty(res)
            orbits[i] = zero(GaussOrbit{T, T})
            continue
        end
        # Current Q
        Q = nms(res)
        # Covariance matrix
        nobs = 2 * notout(res)
        C = (nobs/2) * TS.hessian(Q)
        Γ = inv(C)
        # Residuals space to barycentric coordinates jacobian
        J = Matrix(TS.jacobian(q0 - cte.(q0)))
        # Update orbit
        orbits[i] = evaldeltas(GaussOrbit(tracklets, bwd, fwd, res, Γ, J, τ_1, τ_3,
            ρ_vec, R_vec, D_0, _D_, a, b, c, r_2s[i], r_vec[:, :, i], ρ[:, i]))
    end
    # Sort orbits by nms
    sort!(orbits, by = nms)

    return orbits
end

function gaussmethod(observatories::Vector{ObservatoryMPC{T}}, dates::Vector{DateTime},
    α::Vector{U}, δ::Vector{U}, params::Parameters{T}) where {T <: Real, U <: Number}
    # Check we have exactly three observations
    @assert length(observatories) == length(dates) == length(α) == length(δ) == 3
        "Gauss method requires exactly three observations"
    # Check observations are in temporal order
    @assert issorted(dates) "Observations must be in temporal order"
    # Times of observation [et, days since J2000]
    t_et = dtutc2et.(dates)
    t_days = t_et ./ daysec
    # Time intervals
    τ_1 = t_days[1] - t_days[2]
    τ_3 = t_days[3] - t_days[2]
    τ = τ_3 - τ_1
    # Topocentric (line-of-sight) unit vectors
    ρ_vec = reduce(hcat, topounit.(α, δ))
    # Geocentric state vector of the observer [au, au/day]
    g_vec = @. kmsec2auday(obsposvelECI(observatories, t_et))
    # Heliocentric state vector of the Earth [au, au/day]
    G_vec = @. params.eph_ea(t_days) - params.eph_su(t_days)
    # Observer's heliocentric positions [au, au/day]
    R_vec = reduce(hcat, G_vec .+  g_vec)[1:3, :]
    # Cross products
    p_vec = zeros(U, 3, 3)
    p_vec[:, 1] = @views cross(ρ_vec[:, 2], ρ_vec[:, 3])
    p_vec[:, 2] = @views cross(ρ_vec[:, 1], ρ_vec[:, 3])
    p_vec[:, 3] = @views cross(ρ_vec[:, 1], ρ_vec[:, 2])
    # Gauss scalar
    D_0 = dot(ρ_vec[:, 1], p_vec[:, 1])
    # Matrix of triple products
    D = zeros(U, 3, 3)
    for j in 1:3, i in 1:3
        D[i, j] = @views dot(R_vec[:, i], p_vec[:, j])
    end
    # A and B scalars
    A = (-D[1, 2] * τ_3 / τ + D[2, 2] + D[3, 2] * τ_1 / τ) / D_0
    B = (D[1, 2] * (τ_3^2 - τ^2) * τ_3 / τ + D[3, 2] * (τ^2 - τ_1^2) * τ_1 / τ) / 6 / D_0
    # E and F scalars
    E = @views dot(R_vec[:, 2], ρ_vec[:, 2])
    F = @views dot(R_vec[:, 2], R_vec[:, 2])
    # Lagrange equation coefficients
    a = -(A * A + 2 * A * E + F)
    b = -2 * μ_S * B * (A + E)
    c = -(μ_S^2) * (B * B)
    # Solve Lagrange equation
    r_2s = solve_lagrange(a, b, c; niter = params.lsiter)
    # Heliocentric state vectors of the object
    r_vec = Array{U}(undef, 6, 3, length(r_2s))
    # Slant ranges
    ρ = Matrix{U}(undef, 3, length(r_2s))
    # Lagrange polynomial has no solutions
    if isempty(r_2s)
        @warn("No solutions found for Lagrange equation \
            $(_format_Lagrange_equation(cte(a), cte(b), cte(c)))")
    end
    # Lagrange polynomial has at least one solution
    for (i, r_2) in enumerate(r_2s)
        # Range cubed
        r_23 = r_2 * r_2 * r_2
        # Slant ranges
        num_1 = 6 * (D[3, 1] * τ_1 / τ_3 + D[2, 1] * τ / τ_3) * r_23 +
            μ_S * D[3, 1] * (τ^2 - τ_1^2) * τ_1 / τ_3
        den_1 = 6 * r_23 + μ_S * (τ^2 - τ_3^2)
        ρ[1, i] = (num_1 / den_1 - D[1, 1]) / D_0
        ρ[2, i] = A + μ_S * B / r_23
        num_3 = 6 * (D[1, 3] * τ_3 / τ_1 - D[2, 3] * τ / τ_1) * r_23 +
            μ_S * D[1, 3] * (τ^2 - τ_3^2) * τ_3 / τ_1
        den_3 = 6 * r_23 + μ_S * (τ^2 - τ_1^2)
        ρ[3, i] = (num_3 / den_3 - D[3, 3]) / D_0
        # Heliocentric position of the object
        @. r_vec[1:3, :, i] = R_vec + ρ[:, i]' * ρ_vec
        # f and g Lagrange coefficients
        f_1, f_3 = f_Lagrange(τ_1, r_2), f_Lagrange(τ_3, r_2)
        g_1, g_3 = g_Lagrange(τ_1, r_2), g_Lagrange(τ_3, r_2)
        # Heliocentric velocity of the object
        @. r_vec[4:6, 1, i] = NaN * one(f_1)
        @. r_vec[4:6, 2, i] = (-f_3 * r_vec[1:3, 1, i] + f_1 * r_vec[1:3, 3, i]) /
            (f_1 * g_3 - f_3 * g_1)
        @. r_vec[4:6, 3, i] = NaN * one(f_1)
    end

    return τ_1, τ_3, ρ_vec, R_vec, D_0, D, a, b, c, r_2s, r_vec, ρ
end

@doc raw"""
    gaussiod(od, params) where {D, T <: Real}

Compute a `LeastSquaresOrbit` via Jet Transport Gauss Method followed by
Jet Transport Least Squares.

See also [`gaussmethod`](@ref).

## Arguments

- `od::ODProblem{D, T}`: an orbit determination problem.
- `params::Parameters{T}`: see `Gauss Method` and `Least Squares` sections
    of [`Parameters`](@ref).

!!! warning
    This function may change the (global) `TaylorSeries` variables.

!!! reference
    See https://doi.org/10.1007/s10569-025-10246-2.
"""
function gaussiod(od::ODProblem{D, T}, params::Parameters{T}) where {D, T <: Real}
    # Allocate memory for orbit
    orbit = zero(LeastSquaresOrbit{T, T})
    # Unpack
    @unpack safegauss, significance, verbose = params
    @unpack tracklets, radec = od
    # This function requires exactly 3 tracklets
    (safegauss && length(tracklets) != 3) && return orbit
    # Set jet transport variables
    set_od_order(params)
    # Preliminary orbits
    porbits = gaussmethod(od, params)
    # Filter preliminary orbits
    filter!(!iszero, porbits)
    filter!(isphysical, porbits)
    filter!(isclosed, porbits)
    # Iterate over remaining preliminary orbits
    for i in eachindex(porbits)
        # Jet Transport Least Squares
        _orbit_ = jtls(od, porbits[i], params, true)
        # Update orbit
        orbit = updateorbit(orbit, _orbit_, radec)
        # Termination condition
        if critical_value(orbit) < significance
            verbose && println(
                "* Gauss method converged to:\n\n",
                summary(porbits[i]), "\n",
                "* Jet Transport Least Squares converged to: \n\n",
                summary(orbit)
            )
            return orbit
        end
        # Refine via Minimization over the MOV
        j = safegauss ? 2 : closest_tracklet(epoch(porbits[i]), tracklets)
        nobs(tracklets[j]) < 2 && continue
        porbit = mmov(od, porbits[i], j, params)
        # MMOV failed to converge
        iszero(porbit) && continue
        # Jet Transport Least Squares
        _orbit_ = jtls(od, porbit, params, true)
        # Update orbit
        orbit = updateorbit(orbit, _orbit_, radec)
        # Termination condition
        if critical_value(orbit) < significance
            Niter = length(porbit.Qs)
            verbose && println(
                "* Refinement of GaussOrbit via MMOV converged in \
                $Niter iterations to:\n\n",
                summary(porbit), "\n",
                "* Jet Transport Least Squares converged to: \n\n",
                summary(orbit)
            )
            return orbit
        end
    end
    # Unsuccessful orbit determination
    verbose && @warn("Orbit determination did not converge within \
        the given parameters or could not fit all the astrometry")

    return orbit
end