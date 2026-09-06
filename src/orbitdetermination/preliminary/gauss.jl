# Topocentric (line-of-sight) unit vector
topounit(x::AbstractOpticalAstrometry)= topounit(ra(x), dec(x))

function topounit(α::T, δ::T) where {T <: Number}
    sin_α, cos_α = sincos(α)
    sin_δ, cos_δ = sincos(δ)
    pos = [cos_δ * cos_α, cos_δ * sin_α, sin_δ]
    return pos
end

# 1st order approximation to Lagrange's f and g functions
# TO DO: Should we allow to use other μ?
function f_Lagrange(τ::Real, r::Number)
    r3 = r * r * r
    return 1 - μ_S * (τ^2) / (2r3)
end

function g_Lagrange(τ::Real, r::Number)
    r3 = r * r * r
    return τ - μ_S * (τ^3) / (6r3)
end

# Format Lagrange equation as `r⁸ + a r⁶ + b r³ + c = 0`
function _format_Lagrange_equation(a::T, b::T, c::T) where {T <: Real}
    a_sgn = a ≥ 0 ? "+" : "-"
    b_sgn = b ≥ 0 ? "+" : "-"
    c_sgn = c ≥ 0 ? "+" : "-"
    return string("r⁸ ", a_sgn, " ", abs(a), " r⁶ ", b_sgn, " ", abs(b),
                  " r³ ", c_sgn, " ", abs(c), " = 0")
end

# Lagrange polynomial to be solved during Gauss method and its derivative
function lagrange(x::Number, a::U, b::U, c::U) where {U <: Number}
    # Evaluate via Horner's method
    x2 = x * x
    x3 = x2 * x
    return c + x3 * (b + x3 * (a + x2))
end

function lagrange_derivative(x::Number, a::U, b::U) where {U <: Number}
    # Evaluate via Horner's method
    x2 = x * x
    x3 = x2 * x
    return x2 * (3*b + x3 * (6*a + 8*x2))
end

# Solve Lagrange polynomial
# TO DO: Allow to control interval over which to look for solutions
# Currently we look between the radius of the Sun (∼0.00465047 AU) and
# the radius of the Solar System (∼40 AU)
function solve_lagrange(a::T, b::T, c::T; niter::Int = 5, rmin::Real = 0.00465047,
                        rmax::Real = 40.0) where {T <: Real}
    return find_zeros(x -> lagrange(x, a, b, c), rmin, rmax)
end

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

"""
    gaussmetric(a, b, c; kwargs...)

Heuristic metric to quantify how fit a triplet `(a, b, c)` of optical
astrometry is for Gauss' method.

# Keyword arguments

- `τmin::Real`: minimum desired total time span in days (default: `1`).
- `τmax::Real`: maximum desired total time span in days (default: `30`).

!!! note
    This function assumes that `date(a) < date(b) < date(c)`.
"""
function gaussmetric(a::AbstractOpticalAstrometry, b::AbstractOpticalAstrometry,
                     c::AbstractOpticalAstrometry; τmin::Real = 1, τmax::Real = 30)
    τ1 = daysbetween(a, b)
    τ3 = daysbetween(b, c)
    τ = τ1 + τ3
    x = τ / sqrt(τmin * τmax)
    return ((τ3 - τ1) / τ)^2 + 0.5 * (x + 1/x)
end

function _gausstriplets(X::AbstractOpticalVector, nbest::Int, τmin::Real,
                        τmax::Real, cutoff::Real)
    L = length(X)
    triplets = [(0, 0, 0, Inf) for _ in 1:min(nbest, binomial(L, 3))]
    for i in 1:L-2
        xi = X[i]
        for j in i+1:L-1
            xj = X[j]
            daysbetween(xi, xj) > cutoff && break
            for k in j+1:L
                xk = X[k]
                daysbetween(xi, xk) > cutoff && break
                m = gaussmetric(xi, xj, xk; τmin, τmax)
                if m < triplets[end][4]
                    triplet = (i, j, k, m)
                    l = searchsortedfirst(triplets, triplet, by = last)
                    triplets[l+1:end] .= triplets[l:end-1]
                    triplets[l] = triplet
                end
            end
        end
    end
    return triplets
end

isgausstriplet(x::Tuple{Int, Int, Int, T}) where {T <: Real} = !isinf(x[4])

"""
    gausstriplets(X, nbest; kwargs...)

Find the `nbest` triplets with the smallest value of [`gaussmetric`](@ref)
in a vector of optical astrometry `X`. Returns a vector of tuples, where
each element contains three indices followed by the value of the metric.

# Keyword arguments

- `τmin::Real`: minimum desired total time span in days (default: `1`).
- `τmax::Real`: maximum desired total time span in days (default: `30`).

!!! note
    This function assumes that `X` is sorted according to `date`.
"""
function gausstriplets(X::AbstractOpticalVector, nbest::Int = 10;
                       τmin::Real = 1, τmax::Real = 30)
    # First pass: use τmax as cutoff to aggressively prune the search space
    triplets = _gausstriplets(X, nbest, τmin, τmax, τmax)
    # Fallback: rerun the search without the τmax early-exit constraint
    if !isempty(triplets) && !isgausstriplet(triplets[1])
        triplets = _gausstriplets(X, nbest, τmin, τmax, Inf)
    end
    # Clean up: remove any uninitialized slots
    filter!(isgausstriplet, triplets)
    return triplets
end

"""
    gaussmethod(od, params)

Given an orbit determination problem `od`, return a vector of preliminary
orbits computed by Gauss method. For a list of parameters, see the `Gauss
Method` section of [`Parameters`](@ref).

!!! reference
    See Section 3 of:
    - https://doi.org/10.1007/s10569-025-10246-2
"""
function gaussmethod(od::OpticalODProblem{D, T, O}, params::Parameters{T}) where {D, T, O}
    # Unpack
    @unpack gaussorder, safegauss, gaussntrip, eph_su = params
    @unpack dynamics, optical, tracklets = od
    variables = collect(1:6)
    dq = [arcsec2rad(1) * TaylorN(i, order = gaussorder) for i in 1:6]
    # Find the triplets with the smallest gaussmetric
    safegauss = safegauss && length(tracklets) ≥ 3
    τmin, τmax = params.gaussτmin, params.gaussτmax
    triplets = gausstriplets(safegauss ? tracklets : optical, gaussntrip; τmin, τmax)
    # Iterate triplets
    orbits = Vector{GaussOrbit{D, T, T, O}}(undef, 0)
    sizehint!(orbits, 3 * length(triplets))
    for triplet in triplets
        # Current triplet
        trks = view(safegauss ? tracklets : optical, [triplet[1], triplet[2], triplet[3]])
        observatories, dates, α, δ = observatory.(trks), date.(trks), ra.(trks), dec.(trks)
        # Julian day of middle observation
        _jd0_ = dtutc2jdtdb(dates[2])
        # Gauss method solutions
        τ_1, τ_3, ρ_vec, R_vec, D_0, D_mat, a, b, c, r_2s, r_vec, ρ =
            gaussmethod(observatories, dates, α .+ dq[1:3], δ .+ dq[4:6], params)
        # Iterate solutions
        for i in eachindex(r_2s)
            # Epoch (corrected for light-time)
            jd0 = _jd0_ - cte(ρ[2, i]) / c_au_per_day
            # Jet transport initial condition
            q0 = r_vec[:, 2, i] + eph_su(jd0 - JD_J2000)
            # Propagation and residuals
            bwd, fwd, res = propres(od, q0, jd0, params)
            isempty(res) && continue
            # Current Q
            Q = nms(res)
            # Covariance matrix
            x0 = zeros(T, 6)
            nobs = 2 * notout(res)
            C = (nobs/2) * TS.hessian(Q)
            Γ = project(q0, x0, inv(C))
            # Update orbit
            push!(orbits,
                evaldeltas(GaussOrbit(
                    dynamics, variables, optical, tracklets, bwd, fwd, res, Γ, τ_1, τ_3,
                    ρ_vec, R_vec, D_0, D_mat, a, b, c, r_2s[i], r_vec[:, :, i], ρ[:, i]
                ))
            )
        end
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
    # Times of observation [Julian date UTC, days since J2000 TDB]
    jds_utc = datetime2julian.(dates)
    t_days = dtutc2et.(dates) ./ daysec
    # Time intervals
    τ_1 = t_days[1] - t_days[2]
    τ_3 = t_days[3] - t_days[2]
    τ = τ_3 - τ_1
    # Topocentric (line-of-sight) unit vectors
    ρ_vec = reduce(hcat, topounit.(α, δ))
    # Geocentric state vector of the observer [au, au/day]
    g_vec = @. kmsec2auday(obsposvelECI(observatories, jds_utc))
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

"""
    gaussiod(od, params)

Given an orbit determination problem `od`, compute a `LeastSquaresOrbit`
via Gauss method followed by Jet Transport Least Squares. For a list of
parameters, see the `Gauss Method` and `Least Squares` sections of
[`Parameters`](@ref).

See also [`gaussmethod`](@ref).

!!! warning
    This function may change the (global) `TaylorSeries` variables.

!!! reference
    See section 3 of:
    - https://doi.org/10.1007/s10569-025-10246-2
"""
function gaussiod(od::OpticalODProblem{D, T, O}, params::Parameters{T}) where {D,
                  T <: Real, O <: AbstractOpticalVector{T}}
    # Unpack
    @unpack significance, verbose = params
    @unpack optical, tracklets = od
    # Set jet transport variables
    Npar = numvars(Val(od.dynamics), params)
    @assert Npar == 6 "Gauss method is not defined for dynamical models with \
        more than 6 degrees of freedom"
    set_od_order(params, Npar)
    # Pre-allocate orbit
    orbit = zero(LeastSquaresOrbit{D, T, T, O, Nothing, Nothing})
    # This function requires at least 3 observations
    length(optical) < 3 && return orbit
    # Preliminary orbits
    porbits = gaussmethod(od, params)
    # Filter preliminary orbits
    filter!(!iszero, porbits)
    filter!(isphysical, porbits)
    # filter!(isclosed, porbits)
    # Iterate over remaining preliminary orbits
    for i in eachindex(porbits)
        # Jet Transport Least Squares
        _orbit_ = jtls(od, porbits[i], params, true)
        # Update orbit
        orbit = updateorbit(orbit, _orbit_, optical)
        # Termination condition
        if critical_value(orbit) < significance
            N2 = length(orbit.Qs)
            verbose && println(
                "* Gauss method converged to:\n\n",
                summary(porbits[i]), "\n",
                "* Jet Transport Least Squares converged in $N2 iterations to: \n\n",
                summary(orbit)
            )
            return orbit
        end
        # Refine via Minimization over the MOV
        j = closest_tracklet(epoch(porbits[i]), tracklets)
        nobs(tracklets[j]) < 2 && continue
        for scale in (:log, :linear)
            # Subset of optical to be included in the calculation
            porbit = mmov(od, porbits[i], j, params; scale)
            # MMOV failed to converge
            iszero(porbit) && continue
            # Jet Transport Least Squares
            _orbit_ = jtls(od, porbit, params, true)
            # Update orbit
            orbit = updateorbit(orbit, _orbit_, optical)
            # Termination condition
            if critical_value(orbit) < significance
                N1, N2 = length(porbit.Qs), length(orbit.Qs)
                verbose && println(
                    "* Refinement of GaussOrbit via MMOV converged in $N1 iterations to:\n\n",
                    summary(porbit), "\n",
                    "* Jet Transport Least Squares converged in $N2 iterations to: \n\n",
                    summary(orbit)
                )
                return orbit
            end
        end
    end
    # Unsuccessful orbit determination
    verbose && @warn("Orbit determination did not converge within \
        the given parameters or could not fit all the astrometry")

    return orbit
end