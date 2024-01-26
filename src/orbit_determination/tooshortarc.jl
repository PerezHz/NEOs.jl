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
        DateTime(2000), zero(T), zero(T), zero(T), zero(T), Vector{T}(undef, 0),
        Vector{T}(undef, 0), Vector{T}(undef, 0), Vector{T}(undef, 0), Vector{T}(undef, 0),
        Vector{T}(undef, 0), Vector{T}(undef, 0), Matrix{T}(undef, 0, 0), unknownobs()
    )
end

iszero(x::AdmissibleRegion{T}) where {T <: AbstractFloat} = x == zero(AdmissibleRegion{T})

# Outer constructor
function AdmissibleRegion(tracklet::Tracklet{T}) where {T <: AbstractFloat}
    # Unfold
    obs, t_datetime, α, δ = observatory(tracklet), date(tracklet), ra(tracklet), dec(tracklet)
    v_α, v_δ, h = vra(tracklet), vdec(tracklet), mag(tracklet)
    # Topocentric unit vector and partials
    ρ, ρ_α, ρ_δ = topounitpdv(α, δ)
    # Time of observation [days since J2000]
    t_days = datetime2days(t_datetime)
    # Time of observation [et seconds]
    t_et = datetime2et(t_datetime)
    # Sun (Earth) ephemeris
    eph_su = selecteph(sseph, su)
    eph_ea = selecteph(sseph, ea)
    # Heliocentric position of the observer
    q = eph_ea(t_days) + kmsec2auday(obsposvelECI(obs, t_et)) - eph_su(t_days)
    # Admissible region coefficients
    coeffs = admsreg_coeffs(α, δ, v_α, v_δ, ρ, ρ_α, ρ_δ, q)
    # Tiny object boundary
    H_max = 32.0                # Maximum allowed absolute magnitude
    if isnan(h)
        ρ_min = R_SI
    else
        ρ_min = 10^((h - H_max)/5)
    end
    # Maximum range (heliocentric constraint)
    ρ_max = max_range(coeffs, ρ_min)
    iszero(ρ_max) && return zero(AdmissibleRegion{Float64})
    # Range domain
    ρ_domain = [ρ_min, ρ_max]
    # Range rate domain
    v_ρ_min, v_ρ_max = range_rate(coeffs, ρ_min)[1:2]
    v_ρ_domain = [v_ρ_min, v_ρ_max]
    # Range rate symmetry level
    v_ρ_mid = range_rate(coeffs, ρ_max)[1]
    # Boundary points
    Fs = Matrix{Float64}(undef, 3, 2)
    Fs[1, :] .= [ρ_min, v_ρ_min]
    Fs[2, :] .= [ρ_min, v_ρ_max]
    Fs[3, :] .= [ρ_max, v_ρ_mid]
    
    return AdmissibleRegion{Float64}(t_datetime, α, δ, v_α, v_δ, ρ, ρ_α, ρ_δ, q,
                                     coeffs, ρ_domain, v_ρ_domain, Fs, obs)
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
function admsreg_U(coeffs::Vector{T}, ρ::S) where {T <: AbstractFloat, S <: Number} 
    return coeffs[2]^2/4 - admsreg_W(coeffs, ρ) + 2*k_gauss^2/sqrt(admsreg_S(coeffs, ρ))
end
admsreg_U(A::AdmissibleRegion{T}, ρ::S) where {T <: AbstractFloat, S <: Number} = admsreg_U(A.coeffs, ρ)

@doc raw"""
    admsreg_V(A::AdmissibleRegion{T}, ρ::S) where {T <: AbstractFloat, S <: Number} 

V function of an [`AdmissibleRegion`](@ref). 

!!! reference
    See first equation after (8.9) of https://doi.org/10.1017/CBO9781139175371.
"""
function admsreg_V(coeffs::Vector{T}, ρ::S, v_ρ::S) where {T <: AbstractFloat, S <: Number} 
    return v_ρ^2 + coeffs[2] * v_ρ + admsreg_W(coeffs, ρ) - 2*k_gauss^2/sqrt(admsreg_S(coeffs, ρ))
end
admsreg_V(A::AdmissibleRegion{T}, ρ::S, v_ρ::S) where {T <: AbstractFloat, S <: Number} = admsreg_V(A.coeffs, ρ, v_ρ)

@doc raw"""
    range_rate(A::AdmissibleRegion{T}, ρ::S) where {T, S <: AbstractFloat}

Return the two possible range rates in the boundary of `A` for a given range `ρ`.
"""
function range_rate(coeffs::Vector{T}, ρ::S) where {T, S <: AbstractFloat}
    return find_zeros(s -> admsreg_V(coeffs, ρ, s), -10.0, 10.0)
end
range_rate(A::AdmissibleRegion{T}, ρ::S) where {T, S <: AbstractFloat} = range_rate(A.coeffs, ρ)

@doc raw"""
    max_range(coeffs::Vector{T}, ρ_min::T) where {T <: AbstractFloat}

Return the maximum possible range in the boundary of an admissible region
with coefficients `coeffs` and minimum allowed range `ρ_min`.
"""
function max_range(coeffs::Vector{T}, ρ_min::T) where {T <: AbstractFloat}
    # Initial guess
    sol = find_zeros(s -> admsreg_U(coeffs, s), ρ_min, 10.0)
    iszero(length(sol)) && return zero(T)
    ρ_max = sol[1]
    # Make sure U(ρ) ≥ 0 and there is at least one range_rate solution
    niter = 0
    while admsreg_U(coeffs, ρ_max) < 0 || length(range_rate(coeffs, ρ_max)) == 0
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
    if A.ρ_domain[1] <= P[1] <= A.ρ_domain[2] && A.v_ρ_domain[1] <= P[2] <= A.v_ρ_domain[2]
        y_range = range_rate(A, P[1])
        if length(y_range) == 1
            return P[2] == y_range
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
    ρ = norm((r - A.q)[1:3])
    # Topocentric range rate
    v_ρ = dot(r[4:6], A.ρ_unit) - dot(A.q[4:6], A.ρ_unit) - ρ * A.v_α * dot(A.ρ_α, A.ρ_unit) 
          - ρ * A.v_δ * dot(A.ρ_δ, A.ρ_unit)

    return ρ, v_ρ
end

# Propagate an orbit and compute residuals
function propres(radec::Vector{RadecMPC{T}}, jd0::T, q0::Vector{U},
                 params::NEOParameters{T})  where {T <: AbstractFloat, U <: Number}
    # Time of first (last) observation
    t0, tf = datetime2julian(date(radec[1])), datetime2julian(date(radec[end]))
    # Years in backward (forward) integration
    nyears_bwd = -(jd0 - t0 + params.bwdoffset) / yr
    nyears_fwd = (tf - jd0 + params.fwdoffset) / yr
    # Backward (forward) integration
    bwd = propagate(RNp1BP_pN_A_J23E_J2S_eph_threads!, jd0, nyears_bwd, q0, params)
    fwd = propagate(RNp1BP_pN_A_J23E_J2S_eph_threads!, jd0, nyears_fwd, q0, params)
    if !issuccessfulprop(bwd, t0 - jd0; tol = params.coeffstol) || 
       !issuccessfulprop(fwd, tf - jd0; tol = params.coeffstol)
        return bwd, fwd, Vector{OpticalResidual{T, U}}(undef, 0)
    end
    # Sun (Earth) ephemeris
    eph_su = selecteph(sseph, su)
    eph_ea = selecteph(sseph, ea)
    # O-C residuals
    res = residuals(radec, params;
                    xvs = et -> auday2kmsec(eph_su(et/daysec)),
                    xve = et -> auday2kmsec(eph_ea(et/daysec)),
                    xva = et -> bwdfwdeph(et, bwd, fwd))
    
    return bwd, fwd, res
end

@doc raw"""
    adam(radec::Vector{RadecMPC{T}}, A::AdmissibleRegion{T}, ρ::T, v_ρ::T, 
         params::NEOParameters{T}; η::T = 25.0, μ::T = 0.75, ν::T = 0.9,
         ϵ::T = 1e-8, Qtol::T = 0.001) where {T <: AbstractFloat}

Adaptative moment estimation (ADAM) minimizer of normalized mean square
residual over and admissible region `A`.

!!! warning
    This function will set the (global) `TaylorSeries` variables to `δρ δvρ`. 

!!! reference
    See Algorithm 1 of https://doi.org/10.48550/arXiv.1412.6980.
"""
function adam(radec::Vector{RadecMPC{T}}, A::AdmissibleRegion{T}, ρ::T, v_ρ::T, 
              params::NEOParameters{T}; η::T = 25.0, μ::T = 0.75, ν::T = 0.9,
              ϵ::T = 1e-8, Qtol::T = 0.001) where {T <: AbstractFloat}
    # Initial time of integration [julian days]
    jd0 = datetime2julian(A.date)
    # Scaling factors
    scalings = [A.ρ_domain[2] - A.ρ_domain[1], A.v_ρ_domain[2] - A.v_ρ_domain[1]] / 1_000
    # Jet transport variables
    dq = scaled_variables("dρ dvρ", scalings, order = 1)
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
        q = topo2bary(A, ρ + dq[1], v_ρ + dq[2])
        # Propagation and residuals
        _, _, res = propres(radec, jd0, q, params)
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
    end
    # Find point with smallest Q
    t = argmin(Qs)
    # Return path
    return view(ρs, 1:t), view(v_ρs, 1:t), view(Qs, 1:t)
end

@doc raw"""
    tsals(radec::Vector{RadecMPC{T}}, tracklets::Vector{Tracklet{T}},
          jd0::T, q0::Vector{T}, params::NEOParameters{T}; maxiter::Int = 5,
          Qtol::T = 0.1) where {T <: AbstractFloat}

Used within [`tooshortarc`](@ref) to compute an orbit from a point in an
admissible region via least squares.

!!! warning
    This function will set the (global) `TaylorSeries` variables to `δx₁ δx₂ δx₃ δx₄ δx₅ δx₆`. 
"""
function tsals(radec::Vector{RadecMPC{T}}, tracklets::Vector{Tracklet{T}},
               jd0::T, q0::Vector{T}, params::NEOParameters{T}; maxiter::Int = 5,
               Qtol::T = 0.1) where {T <: AbstractFloat}
    # Scaling factors
    scalings = abs.(q0) ./ 10^5
    # Jet transport variables
    dq = scaled_variables("dx", scalings, order = 6)
    # Allocate memory
    sols = zeros(NEOSolution{T, T}, maxiter+1)
    Qs = fill(T(Inf), maxiter+1)
    # Origin
    x0 = zeros(T, 6)
    # Least squares
    for t in 1:maxiter+1
        # Initial conditions
        q = q0 + dq
        # Propagation & residuals
        bwd, fwd, res = propres(radec, jd0, q, params)
        iszero(length(res)) && break
        # Least squares fit
        fit = tryls(res, x0, params.niter)
        # Case: unsuccessful fit
        if !fit.success || any(diag(fit.Γ) .< 0)
            break
        end
        # Current solution
        sols[t] = evalfit(NEOSolution(tracklets, bwd, fwd, res, fit, scalings))
        Qs[t] = nrms(sols[t])
        # Convergence conditions
        if t > 1
            Qs[t] > Qs[t-1] && break
            abs(Qs[t] - Qs[t-1]) < Qtol && break
        end
        # Update values
        q0 = q(fit.x)
    end
    # Find solution with smallest Q
    t = argmin(Qs)
    # Return solution
    return sols[t]
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
                params::NEOParameters{T}) where {T <: AbstractFloat}

Return initial conditions by minimizing the normalized root mean square residual
over the admissible region.

# Arguments

- `radec::Vector{RadecMPC{T}}`: vector of observations.
- `tracklets::Vector{Tracklet{T}},`: vector of tracklets.
- `params::NEOParameters{T}`: see `Admissible Region Parameters` of [`NEOParameters`](@ref).

!!! warning
    This function will set the (global) `TaylorSeries` variables to `δx₁ δx₂ δx₃ δx₄ δx₅ δx₆`. 
"""
function tooshortarc(radec::Vector{RadecMPC{T}}, tracklets::Vector{Tracklet{T}},
                     params::NEOParameters{T}) where {T <: AbstractFloat}

    # Allocate memory for output
    best_sol = zero(NEOSolution{T, T})
    # Sort tracklets by tsatrackletorder
    idxs = sortperm(tracklets, lt = tsatrackletorder)

    # Iterate tracklets
    for i in idxs
        # Admissible region
        A = AdmissibleRegion(tracklets[i])
        iszero(A) && continue
        # Center
        ρ = sum(A.ρ_domain) / 2
        v_ρ = sum(A.v_ρ_domain) / 2
        # Minimization over admissible region
        ρs, v_ρs, _ = adam(radec, A, ρ, v_ρ, params)
        # Barycentric initial conditions
        q0 = topo2bary(A, ρs[end], v_ρs[end])
        # Initial time of integration [julian days]
        jd0 = datetime2julian(A.date)
        # 6 variables least squares
        sol = tsals(radec, tracklets, jd0, q0, params; maxiter = 5)
        # Update best solution
        if nrms(sol) < nrms(best_sol)
            best_sol = sol
            # Break condition
            nrms(sol) < 1.5 && break
        end
        # Heel anomaly
        ρ = A.ρ_domain[1]
        v_ρs = LinRange(A.v_ρ_domain[1], A.v_ρ_domain[2], 25)
        Qs = fill(Inf, 25)
        for i in eachindex(Qs)
            # Barycentric initial conditions
            q = topo2bary(A, ρ, v_ρs[i])
            # Propagation & residuals
            _, _, res = propres(radec, jd0, q, params)
            iszero(length(res)) && continue
            # NRMS
            Qs[i] = nrms(res) 
        end
        # Find solution with smallest Q
        t = argmin(Qs)
        isinf(Qs[t]) && continue
        # Barycentric initial conditions
        q0 = topo2bary(A, ρ, v_ρs[t])
        # 6 variables least squares
        sol = tsals(radec, tracklets, jd0, q0, params; maxiter = 5)
        # Update best solution
        if nrms(sol) < nrms(best_sol)
            best_sol = sol
            # Break condition
            nrms(sol) < 1.5 && break
        end
    end

    return best_sol
end