@doc raw"""
    AdmissibleRegion{T <: AbstractFloat}

Subset of `topocentric distance` × `topocentric velocity` space defined by
some dynamical constrainits on a too short arc.

# Fields

- `date::DateTime`: time of observation.
- `α::T`: right ascension.
- `δ::T`: declination.
- `v_α::T`: right ascension velocity.
- `v_δ::T`: declination velocity.
- `H_max::T`: maximum absolute magnitude.
- `a_max::T`: maximum semimajor axis.
- `ρ_unit, ρ_α, ρ_δ::Vector{T}`: topocentric unit vector and its partials.
- `q::Vector{T}`: heliocentric position of observer.
- `coeffs::Vector{T}`: polynomial coefficients.
- `ρ_domain::Vector{T}`: range domain.
- `v_ρ_domain::Vector{T}`: range-rate domain.
- `Fs::Matrix{T}`: boundary points.
- `observatory::ObservatoryMPC{T}`: observing station.

!!! reference
    See Chapter 8 of https://doi.org/10.1017/CBO9781139175371.
"""
@auto_hash_equals struct AdmissibleRegion{T <: AbstractFloat}
    date::DateTime
    α::T
    δ::T
    v_α::T
    v_δ::T
    H_max::T
    a_max::T
    ρ_unit::Vector{T}
    ρ_α::Vector{T}
    ρ_δ::Vector{T}
    q::Vector{T}
    coeffs::Vector{T}
    ρ_domain::Vector{T}
    v_ρ_domain::Vector{T}
    Fs::Matrix{T}
    observatory::ObservatoryMPC{T}
end

# Definition of zero AdmissibleRegion{T}
function zero(::Type{AdmissibleRegion{T}}) where {T <: AbstractFloat}
    return AdmissibleRegion{T}(
        DateTime(2000), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T),
        Vector{T}(undef, 0), Vector{T}(undef, 0), Vector{T}(undef, 0),
        Vector{T}(undef, 0), Vector{T}(undef, 0), Vector{T}(undef, 0),
        Vector{T}(undef, 0), Matrix{T}(undef, 0, 0), unknownobs()
    )
end

iszero(x::AdmissibleRegion{T}) where {T <: AbstractFloat} = x == zero(AdmissibleRegion{T})

# Outer constructor
function AdmissibleRegion(tracklet::Tracklet{T}, params::NEOParameters{T}) where {T <: AbstractFloat}
    # Unfold
    obs, t_datetime, α, δ = observatory(tracklet), date(tracklet), ra(tracklet), dec(tracklet)
    v_α, v_δ, h = vra(tracklet), vdec(tracklet), mag(tracklet)
    H_max, a_max = params.H_max, params.a_max
    # Topocentric unit vector and partials
    ρ, ρ_α, ρ_δ = topounitpdv(α, δ)
    # Time of observation [days since J2000]
    t_days = datetime2days(t_datetime)
    # Time of observation [et seconds]
    t_et = datetime2et(t_datetime)
    # Heliocentric position of the observer
    q = params.eph_ea(t_days) + kmsec2auday(obsposvelECI(obs, t_et)) - params.eph_su(t_days)
    # Admissible region coefficients
    coeffs = admsreg_coeffs(α, δ, v_α, v_δ, ρ, ρ_α, ρ_δ, q)
    # Maximum range (heliocentric constraint)
    ρ_max = max_range(coeffs, a_max)
    iszero(ρ_max) && return zero(AdmissibleRegion{T})
    # Minimum range
    if isnan(h)
        if R_SI < ρ_max
            # Earth's sphere of influence radius
            ρ_min = R_SI
        else
            # Earth's physical radius
            ρ_min = R_EA
        end
    else
        # Tiny object boundary
        ρ_min = 10^((h - H_max)/5)
    end
    ρ_min > ρ_max && return zero(AdmissibleRegion{T})
    # Range domain
    ρ_domain = [ρ_min, ρ_max]
    # Range rate domain
    v_ρ_min, v_ρ_max = range_rate(coeffs, ρ_min)[1:2]
    v_ρ_domain = [v_ρ_min, v_ρ_max]
    # Range rate symmetry level
    v_ρ_mid = range_rate(coeffs, ρ_max)[1]
    # Boundary points
    Fs = Matrix{T}(undef, 3, 2)
    Fs[1, :] .= [ρ_min, v_ρ_min]
    Fs[2, :] .= [ρ_min, v_ρ_max]
    Fs[3, :] .= [ρ_max, v_ρ_mid]

    return AdmissibleRegion{T}(t_datetime, α, δ, v_α, v_δ, H_max, a_max,
                               ρ, ρ_α, ρ_δ, q, coeffs, ρ_domain, v_ρ_domain,
                               Fs, obs)
end

@doc raw"""
    admsreg_coeffs(α::T, δ::T, v_α::T, v_δ::T, ρ::Vector{T},
                   ρ_α::Vector{T}, ρ_δ::Vector{T}, q::Vector{T}) where {T <: Number}

Return the polynomial coefficients for an [`AdmissibleRegion`](@ref).

!!! reference
    See equation (8.8) of https://doi.org/10.1017/CBO9781139175371.
"""
function admsreg_coeffs(α::T, δ::T, v_α::T, v_δ::T, ρ::Vector{T},
                        ρ_α::Vector{T}, ρ_δ::Vector{T}, q::Vector{T}) where {T <: Number}
    coeffs = Vector{T}(undef, 6)
    coeffs[1] = dot(q[1:3], q[1:3])
    coeffs[2] = 2 *  dot(q[4:6], ρ)
    coeffs[3] = v_α^2 * cos(δ)^2 + v_δ^2
    coeffs[4] = 2 * v_α * dot(q[4:6], ρ_α) + 2 * v_δ * dot(q[4:6], ρ_δ)
    coeffs[5] = dot(q[4:6], q[4:6])
    coeffs[6] = 2*dot(q[1:3], ρ)
    return coeffs
end

@doc raw"""
    topounitpdv(α::T, δ::T) where {T <: Number}

Return the topocentric unit vector and its partial derivatives with
respect to `α` and `δ`.
"""
function topounitpdv(α::T, δ::T) where {T <: Number}
    sin_α, cos_α = sincos(α)
    sin_δ, cos_δ = sincos(δ)

    ρ = [cos_δ * cos_α, cos_δ * sin_α, sin_δ]
    ρ_α = [-sin_α * cos_δ, cos_α * cos_δ, zero(T)]
    ρ_δ = [-cos_α * sin_δ, -sin_α * sin_δ, cos_δ]

    return ρ, ρ_α, ρ_δ
end

@doc raw"""
    admsreg_W(A::AdmissibleRegion{T}, ρ::S) where {T <: AbstractFloat, S <: Number}

W function of an [`AdmissibleRegion`](@ref).

!!! reference
    See equation (8.9) of https://doi.org/10.1017/CBO9781139175371.
"""
function admsreg_W(coeffs::Vector{T}, ρ::S) where {T <: AbstractFloat, S <: Number}
    return coeffs[3] * ρ^2 + coeffs[4] * ρ + coeffs[5]
end
admsreg_W(A::AdmissibleRegion{T}, ρ::S) where {T <: AbstractFloat, S <: Number}  = admsreg_W(A.coeffs, ρ)

@doc raw"""
    admsreg_S(A::AdmissibleRegion{T}, ρ::S) where {T <: AbstractFloat, S <: Number}

S function of an [`AdmissibleRegion`](@ref).

!!! reference
    See equation (8.9) of https://doi.org/10.1017/CBO9781139175371.
"""
function admsreg_S(coeffs::Vector{T}, ρ::S) where {T <: AbstractFloat, S <: Number}
    return ρ^2 + coeffs[6] * ρ + coeffs[1]
end
admsreg_S(A::AdmissibleRegion{T}, ρ::S) where {T <: AbstractFloat, S <: Number} = admsreg_S(A.coeffs, ρ)

@doc raw"""
    admsreg_U(A::AdmissibleRegion{T}, ρ::S) where {T <: AbstractFloat, S <: Number}

U function of an [`AdmissibleRegion`](@ref).

!!! reference
    See second equation after (8.9) of https://doi.org/10.1017/CBO9781139175371.
"""
function admsreg_U(coeffs::Vector{T}, a_max::T, ρ::S) where {T <: AbstractFloat, S <: Number}
    return coeffs[2]^2/4 - admsreg_W(coeffs, ρ) + 2*k_gauss^2/sqrt(admsreg_S(coeffs, ρ))
           - k_gauss^2/a_max
end
admsreg_U(A::AdmissibleRegion{T}, ρ::S) where {T <: AbstractFloat, S <: Number} = admsreg_U(A.coeffs, A.a_max, ρ)

@doc raw"""
    admsreg_V(A::AdmissibleRegion{T}, ρ::S) where {T <: AbstractFloat, S <: Number}

V function of an [`AdmissibleRegion`](@ref).

!!! reference
    See first equation after (8.9) of https://doi.org/10.1017/CBO9781139175371.
"""
function admsreg_V(coeffs::Vector{T}, a_max::T, ρ::S, v_ρ::S) where {T <: AbstractFloat, S <: Number}
    return v_ρ^2 + coeffs[2] * v_ρ + admsreg_W(coeffs, ρ) - 2*k_gauss^2/sqrt(admsreg_S(coeffs, ρ))
           + k_gauss^2/a_max
end
admsreg_V(A::AdmissibleRegion{T}, ρ::S, v_ρ::S) where {T <: AbstractFloat, S <: Number} = admsreg_V(A.coeffs, A.a_max, ρ, v_ρ)

@doc raw"""
    range_rate(A::AdmissibleRegion{T}, ρ::S) where {T, S <: AbstractFloat}

Return the two possible range rates in the boundary of `A` for a given range `ρ`.
"""
function range_rate(coeffs::Vector{T}, a_max::T, ρ::S) where {T, S <: AbstractFloat}
    return find_zeros(s -> admsreg_V(coeffs, a_max, ρ, s), -10.0, 10.0)
end
range_rate(A::AdmissibleRegion{T}, ρ::S) where {T, S <: AbstractFloat} = range_rate(A.coeffs, A.a_max, ρ)

@doc raw"""
    max_range(coeffs::Vector{T}, a_max::T) where {T <: AbstractFloat}

Return the maximum possible range in the boundary of an admissible region
with coefficients `coeffs` and maximum semimajor axis `a_max`.
"""
function max_range(coeffs::Vector{T}, a_max::T) where {T <: AbstractFloat}
    # Initial guess
    sol = find_zeros(s -> admsreg_U(coeffs, a_max, s), R_EA, 100.0)
    iszero(length(sol)) && return zero(T)
    ρ_max = sol[1]
    # Make sure U(ρ) ≥ 0 and there is at least one range_rate solution
    niter = 0
    while admsreg_U(coeffs, a_max, ρ_max) < 0 || length(range_rate(coeffs, a_max, ρ_max)) == 0
        niter += 1
        ρ_max = prevfloat(ρ_max)
        niter > 1_000 && break
    end

    return ρ_max
end

@doc raw"""
    boundary(A::AdmissibleRegion{T}, t::S) where {T <: AbstractFloat, S <: Number}

Parametrization of `A` boundary with `t ∈ [0, 3]`.
"""
function boundary(A::AdmissibleRegion{T}, t::S) where {T <: AbstractFloat, S <: Number}
    # Parametrization domain
    @assert 0.0 <= t <= 3.0
    # Lower (upper) bounds
    x_min, x_max = A.ρ_domain
    y_min, y_max = A.v_ρ_domain
    # ρ = x_min
    if 0.0 <= t <= 1.0
        return [x_min, y_min + t * (y_max - y_min)]
    else
        # Upper curve
        if 1.0 <= t <= 2.0
            ρ = x_min + (t-1)*(x_max - x_min)
            v_ρ = range_rate(A, ρ)[end]
        # Lower curve
        elseif 2.0 <= t <= 3.0
            ρ = x_max - (t-2)*(x_max - x_min)
            v_ρ = range_rate(A, ρ)[1]
        end
        return [ρ, v_ρ]
    end
end

@doc raw"""
    boundary_projection(A::AdmissibleRegion{T}, ρ::T, v_ρ::T) where {T <: AbstractFloat}

Project `[ρ, v_ρ]` into `A`'s boundary.
"""
function boundary_projection(A::AdmissibleRegion{T}, ρ::T, v_ρ::T) where {T <: AbstractFloat}
    # Project range
    ρ = clamp(ρ,  A.ρ_domain[1], A.ρ_domain[2])
    # Project range-rate
    y_domain = range_rate(A, ρ)
    if iszero(length(y_domain))
        v_ρ = sum(A.v_ρ_domain) / 2
    elseif isone(length(y_domain))
        v_ρ = y_domain[1]
    else
        v_ρ = clamp(v_ρ, y_domain[1], y_domain[2])
    end

    return ρ, v_ρ
end

# Check whether P is inside A's boundary
function in(P::Vector{T}, A::AdmissibleRegion{T}) where {T <: AbstractFloat}
    @assert length(P) == 2 "Points in admissible region are of dimension 2"
    if A.ρ_domain[1] <= P[1] <= A.ρ_domain[2]
        y_range = range_rate(A, P[1])
        if length(y_range) == 1
            return P[2] == y_range[1]
        else
            return y_range[1] <= P[2] <= y_range[2]
        end
    else
        return false
    end
end

@doc raw"""
    topo2bary(A::AdmissibleRegion{T}, ρ::U, v_ρ::U) where {T <: AbstractFloat, U <: Number}

Convert topocentric range/range-rate `[ρ, v_ρ]` to barycentric cartesian coordinates.
`A` fixes the line of sight.
"""
function topo2bary(A::AdmissibleRegion{T}, ρ::U, v_ρ::U) where {T <: AbstractFloat, U <: Number}
    # Barycentric position
    r = A.q[1:3] + ρ * A.ρ_unit + sseph(su, datetime2days(A.date))[1:3]
    # Barycentric velocity
    v = A.q[4:6] + v_ρ * A.ρ_unit + ρ * A.v_α * A.ρ_α + ρ * A.v_δ * A.ρ_δ
        + sseph(su, datetime2days(A.date))[4:6]
    # Barycentric state vector
    return vcat(r, v)
end

@doc raw"""
    bary2topo(A::AdmissibleRegion{T}, q0::Vector{U}) where {T <: AbstractFloat, U <: Number}

Convert barycentric cartesian coordinates `q0` to topocentric range/range-rate.
`A` fixes the line of sight.
"""
function bary2topo(A::AdmissibleRegion{T}, q0::Vector{U}) where {T <: AbstractFloat, U <: Number}
    # Heliocentric state vector
    r = q0 - sseph(su, datetime2days(A.date))
    # Topocentric range
    ρ = euclid3D(r - A.q)
    # Topocentric range rate
    v_ρ = dot3D(r[4:6], A.ρ_unit) - dot3D(A.q[4:6], A.ρ_unit) - ρ * A.v_α * dot3D(A.ρ_α, A.ρ_unit)
          - ρ * A.v_δ * dot3D(A.ρ_δ, A.ρ_unit)

    return ρ, v_ρ
end

# Propagate an orbit and compute residuals
function propres(radec::Vector{RadecMPC{T}}, jd0::T, q0::Vector{U},
                 params::NEOParameters{T}; dynamics::D=newtonian!)  where {D, T <: AbstractFloat, U <: Number}
    # Time of first (last) observation
    t0, tf = datetime2julian(date(radec[1])), datetime2julian(date(radec[end]))
    # Years in backward (forward) integration
    nyears_bwd = -(jd0 - t0 + params.bwdoffset) / yr
    nyears_fwd = (tf - jd0 + params.fwdoffset) / yr
    # Backward (forward) integration
    bwd = propagate(dynamics, jd0, nyears_bwd, q0, params)
    fwd = propagate(dynamics, jd0, nyears_fwd, q0, params)
    if !issuccessfulprop(bwd, t0 - jd0; tol = params.coeffstol) ||
       !issuccessfulprop(fwd, tf - jd0; tol = params.coeffstol)
        return bwd, fwd, Vector{OpticalResidual{T, U}}(undef, 0)
    end
    # O-C residuals
    res = residuals(radec, params;
                    xvs = et -> auday2kmsec(params.eph_su(et/daysec)),
                    xve = et -> auday2kmsec(params.eph_ea(et/daysec)),
                    xva = et -> bwdfwdeph(et, bwd, fwd))

    return bwd, fwd, res
end

@doc raw"""
    adam(radec::Vector{RadecMPC{T}}, A::AdmissibleRegion{T}, ρ::T, v_ρ::T,
         params::NEOParameters{T}; scale::Symbol = :linear, η::T = 25.0,
         μ::T = 0.75, ν::T = 0.9, ϵ::T = 1e-8, Qtol::T = 0.001,
         dynamics::D = newtonian!) where {T <: AbstractFloat, D}

Adaptative moment estimation (ADAM) minimizer of normalized mean square
residual over and admissible region `A`.

!!! warning
    This function will set the (global) `TaylorSeries` variables to `dx dy`.

!!! reference
    See Algorithm 1 of https://doi.org/10.48550/arXiv.1412.6980.
"""
function adam(radec::Vector{RadecMPC{T}}, A::AdmissibleRegion{T}, ρ::T, v_ρ::T,
              params::NEOParameters{T}; scale::Symbol = :linear, η::T = 25.0,
              μ::T = 0.75, ν::T = 0.9, ϵ::T = 1e-8, Qtol::T = 0.001,
              dynamics::D = newtonian!) where {T <: AbstractFloat, D}
    # Initial time of integration [julian days]
    jd0 = datetime2julian(A.date)
    # Scaling factors
    if scale == :linear
        scalings = [
            A.ρ_domain[2] - A.ρ_domain[1],
            A.v_ρ_domain[2] - A.v_ρ_domain[1]
        ] / 1_000
    elseif scale == :log
        x = log10(ρ)
        scalings = [
            log10(A.ρ_domain[2]) - log10(A.ρ_domain[1]),
            A.v_ρ_domain[2] - A.v_ρ_domain[1]
        ] / 1_000
    end
    # Jet transport variables
    dq = scaled_variables("dx dy", scalings, order = 1)
    # Maximum number of iterations
    maxiter = params.maxiter
    # Allocate memory
    ρs = Vector{T}(undef, maxiter+1)
    v_ρs = Vector{T}(undef, maxiter+1)
    Qs = fill(T(Inf), maxiter+1)
    # Origin
    x0 = zeros(T, 2)
    # First momentum
    m = zeros(T, 2)
    _m_ = zeros(T, 2)
    # Second momentum
    n = zeros(T, 2)
    _n_ = zeros(T, 2)
    # Gradient descent
    for t in 1:maxiter+1
        # Current position in admissible region
        ρs[t] = ρ
        v_ρs[t] = v_ρ
        # Current barycentric state vector
        if scale == :linear
            q = topo2bary(A, ρ + dq[1], v_ρ + dq[2])
        elseif scale == :log
            q = topo2bary(A, 10^(x + dq[1]), v_ρ + dq[2])
        end
        # Propagation and residuals
        _, _, res = propres(radec, jd0 - ρ/c_au_per_day, q, params; dynamics)
        iszero(length(res)) && break
        # Current Q
        Q = nms(res)
        Qs[t] = Q(x0)
        # Convergence condition
        t > 1 && abs(Qs[t] - Qs[t-1]) / Qs[t] < Qtol && break
        # Gradient of objective function
        g_t = TaylorSeries.gradient(Q)(x0)
        # First momentum
        m .= μ * m + (1 - μ) * g_t
        _m_ .= m / (1 - μ^t)
        # Second momentum
        n .= ν * n + (1 - ν) * g_t .^ 2
        _n_ .= n / (1 - ν^t)
        # Step
        x1 = x0 - η * _m_ ./ (sqrt.(_n_) .+ ϵ)
        # Update values
        ρ, v_ρ = bary2topo(A, q(x1))
        # Projection
        ρ, v_ρ = boundary_projection(A, ρ, v_ρ)
        if scale == :log
            x = log10(ρ)
        end
    end
    # Find point with smallest Q
    t = argmin(Qs)

    return T(ρs[t]), T(v_ρs[t]), T(Qs[t])
end

@doc raw"""
    ρminmontecarlo(radec::Vector{RadecMPC{T}}, A::AdmissibleRegion{T},
                   params::NEOParameters{T}; N_samples::Int = 25) where {T <: AbstractFloat}

Monte Carlo sampling over the left boundary of `A`.
"""
function ρminmontecarlo(radec::Vector{RadecMPC{T}}, A::AdmissibleRegion{T},
                        params::NEOParameters{T}; N_samples::Int = 25, dynamics::D=newtonian!) where {T <: AbstractFloat, D}
    # Initial time of integration [julian days]
    jd0 = datetime2julian(A.date)
    # Range lower bound
    ρ = A.ρ_domain[1]
    # Sample range rate
    v_ρs = LinRange(A.v_ρ_domain[1], A.v_ρ_domain[2], N_samples)
    # Allocate memory
    Qs = fill(Inf, N_samples)
    # Monte Carlo
    for i in eachindex(Qs)
        # Barycentric initial conditions
        q = topo2bary(A, ρ, v_ρs[i])
        # Propagation & residuals
        _, _, res = propres(radec, jd0, q, params; dynamics)
        iszero(length(res)) && continue
        # NRMS
        Qs[i] = nrms(res)
    end
    # Find solution with smallest Q
    t = argmin(Qs)

    return T(ρ), T(v_ρs[t]), T(Qs[t])
end

@doc raw"""
    tsals(A::AdmissibleRegion{T}, radec::Vector{RadecMPC{T}}, tracklets::Vector{Tracklet{T}},
          i::Int, ρ::T, v_ρ::T, params::NEOParameters{T}; maxiter::Int = 5) where {T <: AbstractFloat}

Used within [`tooshortarc`](@ref) to compute an orbit from a point in an
admissible region via least squares.

!!! warning
    This function will set the (global) `TaylorSeries` variables to `dx₁ dx₂ dx₃ dx₄ dx₅ dx₆`.
"""
function tsals(A::AdmissibleRegion{T}, radec::Vector{RadecMPC{T}}, tracklets::Vector{Tracklet{T}},
               i::Int, ρ::T, v_ρ::T, params::NEOParameters{T}; maxiter::Int = 5, dynamics::D=newtonian!) where {T <: AbstractFloat, D}
    # Initial time of integration [julian days]
    # (corrected for light-time)
    jd0 = datetime2julian(A.date) - ρ / c_au_per_day
    # Barycentric initial conditions
    q0 = topo2bary(A, ρ, v_ρ)
    # Scaling factors
    scalings = abs.(q0) ./ 10^5
    # Jet transport variables
    dq = scaled_variables("dx", scalings, order = 6)
    # Origin
    x0 = zeros(T, 6)
    # Subset of radec for orbit fit
    g_0 = i
    g_f = i
    idxs = indices(tracklets[i])
    sort!(idxs)
    # Allocate memory
    best_sol = zero(NEOSolution{T, T})
    best_Q = T(Inf)
    flag = false
    # Least squares
    for _ in 1:maxiter
        # Initial conditions
        q = q0 + dq
        # Propagation & residuals
        bwd, fwd, res = propres(radec, jd0, q, params; dynamics)
        iszero(length(res)) && break
        # Orbit fit
        fit = tryls(res[idxs], x0, params.niter)
        !fit.success && break
        # Right iteration
        for k in g_f+1:length(tracklets)
            extra = indices(tracklets[k])
            fit_new = tryls(res[idxs ∪ extra], x0, params.niter)
            if fit_new.success
                fit = fit_new
                idxs = vcat(idxs, extra)
                sort!(idxs)
                g_f = k
            else
                break
            end
        end
        # Left iteration
        for k in g_0-1:-1:1
            extra = indices(tracklets[k])
            fit_new = tryls(res[idxs ∪ extra], x0, params.niter)
            if fit_new.success
                fit = fit_new
                idxs = vcat(idxs, extra)
                sort!(idxs)
                g_0 = k
            else
                break
            end
        end
        # NRMS
        Q = nrms(res, fit)
        if length(idxs) == length(radec) && abs(best_Q - Q) < 0.1
            flag = true
        end
        # Update NRMS and initial conditions
        if Q < best_Q
            best_Q = Q
            best_sol = evalfit(NEOSolution(tracklets[g_0:g_f], bwd, fwd,
                               res[idxs], fit, scalings))
            flag && break
        else
            break
        end
        # Update values
        q0 = q(fit.x)
    end
    # Case: all solutions were unsuccesful
    if isinf(best_Q)
        return zero(NEOSolution{T, T})
    # Case: at least one solution was succesful
    else
        return best_sol
    end
end

# Order in which to check tracklets in tooshortarc
function tsatrackletorder(x::Tracklet{T}, y::Tracklet{T}) where {T <: AbstractFloat}
    if x.nobs == y.nobs
        return x.date > y.date
    else
        return x.nobs > y.nobs
    end
end

@doc raw"""
    tooshortarc(radec::Vector{RadecMPC{T}}, tracklets::Vector{Tracklet{T}},
                params::NEOParameters{T}; dynamics::D = newtonian!) where {T <: AbstractFloat, D}

Return initial conditions by minimizing the normalized root mean square residual
over the admissible region.

# Arguments

- `radec::Vector{RadecMPC{T}}`: vector of observations.
- `tracklets::Vector{Tracklet{T}},`: vector of tracklets.
- `params::NEOParameters{T}`: see `Admissible Region Parameters` of [`NEOParameters`](@ref).
- `dynamics::D`: dynamical model.

!!! warning
    This function will set the (global) `TaylorSeries` variables to `dx₁ dx₂ dx₃ dx₄ dx₅ dx₆`.
"""
function tooshortarc(radec::Vector{RadecMPC{T}}, tracklets::Vector{Tracklet{T}},
                     params::NEOParameters{T}; dynamics::D = newtonian!) where {T <: AbstractFloat, D}

    # Allocate memory for output
    best_sol = zero(NEOSolution{T, T})
    # Sort tracklets by tsatrackletorder
    idxs = sortperm(tracklets, lt = tsatrackletorder)

    # Iterate tracklets
    for i in idxs
        # Admissible region
        A = AdmissibleRegion(tracklets[i], params)
        iszero(A) && continue
        # Center
        ρ = sum(A.ρ_domain) / 2
        v_ρ = sum(A.v_ρ_domain) / 2
        # ADAM minimization over admissible region
        ρ, v_ρ, Q = adam(radec, A, ρ, v_ρ, params; scale = :linear, dynamics = dynamics)
        if !isinf(Q)
            # 6 variables least squares
            sol = tsals(A, radec, tracklets, i, ρ, v_ρ, params; maxiter = 5, dynamics)
            # Update best solution
            if nrms(sol) < nrms(best_sol)
                best_sol = sol
                # Break condition
                nrms(sol) < params.tsaQmax && break
            end
        end
        # Left boundary
        ρ = A.ρ_domain[1]
        v_ρ = sum(A.v_ρ_domain) / 2
        # ADAM minimization over admissible region
        ρ, v_ρ, Q = adam(radec, A, ρ, v_ρ, params; scale = :log, dynamics = dynamics)
        if !isinf(Q)
            # 6 variables least squares
            sol = tsals(A, radec, tracklets, i, ρ, v_ρ, params; maxiter = 5, dynamics)
            # Update best solution
            if nrms(sol) < nrms(best_sol)
                best_sol = sol
                # Break condition
                nrms(sol) < params.tsaQmax && break
            end
        end
    end

    return best_sol
end