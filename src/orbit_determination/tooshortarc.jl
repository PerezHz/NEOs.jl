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
- `ρ_unit, ρ_α, _δ::Vector{T}`: topocentric unit vector and its partials.
- `q::Vector{T}`: heliocentric position of observer.
- `coeffs::Vector{T}`: polynomial coefficients.
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
    observatory::ObservatoryMPC{T}
end

# Outer constructor
function AdmissibleRegion(cdf::DataFrameRow)
    # Unfold
    obs, t_datetime, α, δ, v_α, v_δ = cdf
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

    return AdmissibleRegion{Float64}(t_datetime, α, δ, v_α, v_δ, ρ, ρ_α, ρ_δ, q,
                                     coeffs, obs)
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
    admsreg_coeffs(α::T, δ::T, v_α::T, v_δ::T, ρ::Vector{T},
                   ρ_α::Vector{T}, ρ_δ::Vector{T}, q::Vector{T}) where {T <: Number}

Return the polynomial coefficients for an [`AdmissibleRegioin`](@ref). 

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
    solve_quadratic(a::T, b::T, c::T) where {T <: AbstractFloat}

Return real solutions of cuadratic equation `ax² + bx + c = 0`.
"""
function solve_quadratic(a::T, b::T, c::T) where {T <: AbstractFloat}
    # Discriminant
    Δ = b^2 - 4*a*c
    # Cases
    if Δ > 0
        x_1 = (-b + sqrt(Δ))/(2*a)
        x_2 = (-b - sqrt(Δ))/(2*a)
    elseif iszero(Δ)
        x_1 = -b / (2*a)
        x_2 = T(NaN)
    else
        x_1 = T(NaN)
        x_2 = T(NaN)
    end
    
    return x_1, x_2
end

@doc raw"""
    admsreg_W(A::AdmissibleRegion{T}, ρ::S) where {T <: AbstractFloat, S <: Number} 

W function of an [`AdmissibleRegioin`](@ref). 

!!! reference
    See equation (8.9) of https://doi.org/10.1017/CBO9781139175371.
"""
function admsreg_W(A::AdmissibleRegion{T}, ρ::S) where {T <: AbstractFloat, S <: Number} 
    return A.coeffs[3] * ρ^2 + A.coeffs[4] * ρ + A.coeffs[5]
end

@doc raw"""
    admsreg_S(A::AdmissibleRegion{T}, ρ::S) where {T <: AbstractFloat, S <: Number} 

S function of an [`AdmissibleRegioin`](@ref). 

!!! reference
    See equation (8.9) of https://doi.org/10.1017/CBO9781139175371.
"""
function admsreg_S(A::AdmissibleRegion{T}, ρ::S) where {T <: AbstractFloat, S <: Number} 
    return ρ^2 + A.coeffs[6] * ρ + A.coeffs[1]
end

@doc raw"""
    range_rate(A::AdmissibleRegion{T}, ρ::S) where {T, S <: AbstractFloat}

Return the two possible range rates in the boundary of `A` for a given range `ρ`.
"""
function range_rate(A::AdmissibleRegion{T}, ρ::S) where {T, S <: AbstractFloat}
    c = admsreg_W(A, ρ) - 2*k_gauss^2/sqrt(admsreg_S(A, ρ))
    sol = solve_quadratic(one(T), A.coeffs[2], c)
    return minmax(sol[1], sol[2])
end

@doc raw"""
    max_range(A::AdmissibleRegion{T}) where {T <: AbstractFloat}

Return the maximum range for the left connected component of `A`.
"""
function max_range(A::AdmissibleRegion{T}) where {T <: AbstractFloat}
    # Find max_range    
    ρ = find_zeros(s -> A.coeffs[2]^2/4 - admsreg_W(A, s) + 2*k_gauss^2/sqrt(admsreg_S(A, s)),
                   R_SI, 10.0)[1]
    # Make sure there is no sqrt(negative)
    niter = 0
    while admsreg_S(A, ρ) < 0 && niter < 100
        niter += 1
        ρ -= eps()
    end
    
    return ρ
end

@doc raw"""
    boundary(A::AdmissibleRegion{T}, t::S) where {T <: AbstractFloat, S <: Number}

Parametrization of `A` boundary with `t ∈ [0, 3]`.
"""
function boundary(A::AdmissibleRegion{T}, t::S) where {T <: AbstractFloat, S <: Number}
    # Parametrization domain
    @assert 0.0 <= t <= 3.0
    # Lower (upper) bounds
    y_min, y_max = range_rate(A, R_SI)
    # ρ =  R_SI
    if 0.0 <= t <= 1.0
        return [R_SI, y_min + t * (y_max - y_min)]
    else 
        x_max = max_range(A)
        # Upper curve
        if 1.0 <= t <= 2.0
            ρ = R_SI + (t-1)*(x_max - R_SI)
            _, v_ρ = range_rate(A, ρ)
        # Lower curve
        elseif 2.0 <= t <= 3.0
            ρ = x_max - (t-2)*(x_max - R_SI)
            v_ρ, _ = range_rate(A, ρ)
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
    ρ = clamp(ρ,  R_SI, max_range(A))
    # Project range-rate
    y_min, y_max = range_rate(A, ρ)
    v_ρ = clamp(v_ρ, y_min, y_max)
    
    return ρ, v_ρ
end

# Check whether P is inside A's boundary
function in(P::Vector{T}, A::AdmissibleRegion{T}) where {T <: AbstractFloat}
    @assert length(P) == 2 "Points in admissible region are of dimension 2"
    if R_SI <= P[1] <= max_range(A) && boundary(A, 0.0)[2] <= P[2] <= boundary(A, 1.0)[2]
        y_min, y_max = range_rate(A, P[1])
        return y_min <= P[2] <= y_max
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
                 params::Parameters{T})  where {T <: AbstractFloat, U <: Number}
    # Time of first (last) observation
    t0, tf = datetime2julian(date(radec[1])), datetime2julian(date(radec[end]))
    # Years in backward (forward) integration
    nyears_bwd = -(jd0 - t0 + 0.5) / yr
    nyears_fwd = (tf - jd0 + 0.5) / yr
    # Backward (forward) integration
    bwd = propagate(RNp1BP_pN_A_J23E_J2S_eph_threads!, jd0, nyears_bwd, q0, params)
    fwd = propagate(RNp1BP_pN_A_J23E_J2S_eph_threads!, jd0, nyears_fwd, q0, params)
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
    gradient_descent(radec::Vector{RadecMPC{T}}, A::AdmissibleRegion{T}, ρ::T, v_ρ::T, 
                     params::Parameters{T}; maxiter::Int = 200) where {T <: AbstractFloat}

Gradient descent minimizer of root mean square error over `A`.
"""
function gradient_descent(radec::Vector{RadecMPC{T}}, A::AdmissibleRegion{T}, ρ::T, v_ρ::T, 
                          params::Parameters{T}; maxiter::Int = 200) where {T <: AbstractFloat}
    # Origin
    x0 = zeros(T, 2)
    # Scaling factors
    scalings = [1e-3, 1e-5]
    # Jet transport variables
    dq = scaled_variables("dρ dvρ", scalings, order = 1)
    # Allocate memory
    ρs = Vector{T}(undef, maxiter+1)
    v_ρs = Vector{T}(undef, maxiter+1)
    Qs = fill(T(Inf), maxiter+1)
    # Initial time of integration [julian days]
    jd0 = datetime2julian(A.date)
    # Initial conditions
    ρs[1] = ρ
    v_ρs[1] = v_ρ  
    q0 = topo2bary(A, ρ + dq[1], v_ρ + dq[2])
    _, _, res = propres(radec, jd0, q0, params)
    Q = nms(res)
    Qs[1] = Q(x0)
    # Number of iterations
    niter = 0
    # Gradient descent
    while niter < maxiter
        # Update number of iterations
        niter += 1
        # Gradient of objective function
        dQ = TaylorSeries.gradient(Q)(x0)
        # Step
        x1 = x0 - normalize(dQ)
        # Update values
        ρ, v_ρ = bary2topo(A, q0(x1))
        # Projection
        ρ, v_ρ = boundary_projection(A, ρ, v_ρ)
        # New initial conditions
        q0 = topo2bary(A, ρ + dq[1], v_ρ + dq[2])
        # Update objective function
        _, _, res = propres(radec, jd0, q0, params)
        Q = nms(res)
        # Convergence condition
        if Q(x0) > Qs[niter]
            break
        end
        # Save current variables
        ρs[niter+1] = ρ
        v_ρs[niter+1] = v_ρ
        Qs[niter+1] = Q(x0)
    end

    niter = argmin(Qs)

    return ρs[1:niter], v_ρs[1:niter], Qs[1:niter]
end

@doc raw"""
    momentum_descent(radec::Vector{RadecMPC{T}}, A::AdmissibleRegion{T}, ρ::T, v_ρ::T, 
                     params::Parameters{T}; maxiter::Int = 200, α::T = 10.0,
                     β::T = 0.75, Qtol::T = 0.01) where {T <: AbstractFloat}

Momentum gradient descent minimizer of root mean square error over `A`.
"""
function momentum_descent(radec::Vector{RadecMPC{T}}, A::AdmissibleRegion{T}, ρ::T, v_ρ::T, 
                          params::Parameters{T}; maxiter::Int = 200, α::T = 10.0,
                          β::T = 0.75, Qtol::T = 0.01) where {T <: AbstractFloat}
    # Origin
    x0 = zeros(T, 2)
    # Scaling factors
    scalings = [1e-3, 1e-5]
    # Jet transport variables
    dq = scaled_variables("dρ dvρ", scalings, order = 1)
    # Allocate memory
    ρs = Vector{T}(undef, maxiter+1)
    v_ρs = Vector{T}(undef, maxiter+1)
    Qs = fill(T(Inf), maxiter+1)
    # Initial time of integration [julian days]
    jd0 = datetime2julian(A.date)
    # Initial conditions
    ρs[1] = ρ
    v_ρs[1] = v_ρ  
    q0 = topo2bary(A, ρ + dq[1], v_ρ + dq[2])
    _, _, res = propres(radec, jd0, q0, params)
    Q = nms(res)
    Qs[1] = Q(x0)
    m = zeros(T, 2)
    # Number of iterations
    niter = 0
    # Gradient descent
    while niter < maxiter
        # Update number of iterations
        niter += 1
        # Gradient of objective function
        dQ = TaylorSeries.gradient(Q)(x0)
        # Momentum
        m = β * m + (1 - β) * dQ
        # Step
        x1 = x0 - α * m
        # Update values
        ρ, v_ρ = bary2topo(A, q0(x1))
        # Projection
        ρ, v_ρ = boundary_projection(A, ρ, v_ρ)
        # New initial conditions
        q0 = topo2bary(A, ρ + dq[1], v_ρ + dq[2])
        # Update obbjective function
        _, _, res = propres(radec, jd0, q0, params)
        Q = nms(res)
        # Convergence condition
        if abs(Q(x0) - Qs[niter]) < Qtol
            break
        end
        # Save current variables
        ρs[niter+1] = ρ
        v_ρs[niter+1] = v_ρ
        Qs[niter+1] = Q(x0)
    end

    niter = argmin(Qs)

    return ρs[1:niter], v_ρs[1:niter], Qs[1:niter]
end

@doc raw"""
    rmsprop(radec::Vector{RadecMPC{T}}, A::AdmissibleRegion{T}, ρ::T, v_ρ::T, 
            params::Parameters{T}; maxiter::Int = 200, α::T = 10.0, β::T = 0.9,
            ϵ::T = 1e-8, Qtol::T = 0.01) where {T <: AbstractFloat}

RMSPROP minimizer of root mean square error over `A`.
"""
function rmsprop(radec::Vector{RadecMPC{T}}, A::AdmissibleRegion{T}, ρ::T, v_ρ::T, 
                 params::Parameters{T}; maxiter::Int = 200, α::T = 10.0, β::T = 0.9,
                 ϵ::T = 1e-8, Qtol::T = 0.01) where {T <: AbstractFloat}
    # Origin
    x0 = zeros(T, 2)
    # Scaling factors
    scalings = [1e-3, 1e-5]
    # Jet transport variables
    dq = scaled_variables("dρ dvρ", scalings, order = 1)
    # Allocate memory
    ρs = Vector{T}(undef, maxiter+1)
    v_ρs = Vector{T}(undef, maxiter+1)
    Qs = fill(T(Inf), maxiter+1)
    # Initial time of integration [julian days]
    jd0 = datetime2julian(A.date)
    # Initial conditions
    ρs[1] = ρ
    v_ρs[1] = v_ρ  
    q0 = topo2bary(A, ρ + dq[1], v_ρ + dq[2])
    Q, _ = nms6v(radec, q0, jd0, params)
    Qs[1] = Q(x0)
    v = zeros(T, 2)
    # Number of iterations
    niter = 0
    # Gradient descent
    while niter < maxiter
        # Update number of iterations
        niter += 1
        # Gradient of objective function
        dQ = TaylorSeries.gradient(Q)(x0)
        # Sum of square of past gradients
        v = β * v + (1 - β) * dQ .^ 2
        # Step
        x1 = x0 - α * dQ ./ sqrt.(v .+ ϵ)
        # Update values
        ρ, v_ρ = bary2topo(A, q0(x1))
        # Projection
        ρ, v_ρ = boundary_projection(A, ρ, v_ρ)
        # New initial conditions
        q0 = topo2bary(A, ρ + dq[1], v_ρ + dq[2])
        # Update obbjective function
        _, _, res = propres(radec, jd0, q0, params)
        Q = nms(res)
        # Convergence condition
        if abs(Q(x0) - Qs[niter]) < Qtol
            break
        end
        # Save current variables
        ρs[niter+1] = ρ
        v_ρs[niter+1] = v_ρ
        Qs[niter+1] = Q(x0)
    end

    niter = argmin(Qs)

    return ρs[1:niter], v_ρs[1:niter], Qs[1:niter]
end

@doc raw"""
    adam(radec::Vector{RadecMPC{T}}, A::AdmissibleRegion{T}, ρ::T, v_ρ::T, 
         params::Parameters{T}; maxiter::Int = 200, α::T = 10.0, β_1::T = 0.75,
         β_2::T = 0.85, ϵ::T = 1e-8, Qtol::T =  0.01) where {T <: AbstractFloat}

ADAM minimizer of root mean square error over `A`.
"""
function adam(radec::Vector{RadecMPC{T}}, A::AdmissibleRegion{T}, ρ::T, v_ρ::T, 
              params::Parameters{T}; maxiter::Int = 200, α::T = 25.0, β_1::T = 0.5,
              β_2::T = 0.85, ϵ::T = 1e-8, Qtol::T = 0.01) where {T <: AbstractFloat}
    # Origin
    x0 = zeros(T, 2)
    # Scaling factors
    scalings = [max_range(A) - R_SI, boundary(A, 1)[2]- boundary(A, 0)[2]] / 1_000
    # Jet transport variables
    dq = scaled_variables("dρ dvρ", scalings, order = 1)
    # Allocate memory
    ρs = Vector{T}(undef, maxiter+1)
    v_ρs = Vector{T}(undef, maxiter+1)
    Qs = fill(T(Inf), maxiter+1)
    # Initial time of integration [julian days]
    jd0 = datetime2julian(A.date)
    # Initial conditions
    ρs[1] = ρ
    v_ρs[1] = v_ρ  
    q0 = topo2bary(A, ρ + dq[1], v_ρ + dq[2])
    _, _, res = propres(radec, jd0, q0, params)
    Q = nms(res)
    Qs[1] = Q(x0)
    m = zeros(T, 2)
    v = zeros(T, 2)
    # Number of iterations
    niter = 0
    # Gradient descent
    while niter < maxiter
        # Update number of iterations
        niter += 1
        # Gradient of objective function
        dQ = TaylorSeries.gradient(Q)(x0)
        # Momentum
        m = β_1 * m + (1 - β_1) * dQ
        _m_ = m / (1 - β_1^niter)
        # Sum of square of past gradients
        v = β_2 * v + (1 - β_2) * dQ .^ 2
        _v_ = v / (1 - β_2^niter)
        # Step
        x1 = x0 - α * _m_ ./ (sqrt.(_v_) .+ ϵ)
        # Update values
        ρ, v_ρ = bary2topo(A, q0(x1))
        # Projection
        ρ, v_ρ = boundary_projection(A, ρ, v_ρ)
        # New initial conditions
        q0 = topo2bary(A, ρ + dq[1], v_ρ + dq[2])
        # Update obbjective function
        _, _, res = propres(radec, jd0, q0, params)
        Q = nms(res)
        # Convergence condition
        if abs(Q(x0) - Qs[niter]) < Qtol
            break
        end
        # Save current variables
        ρs[niter+1] = ρ
        v_ρs[niter+1] = v_ρ
        Qs[niter+1] = Q(x0)
    end

    niter = argmin(Qs)

    return ρs[1:niter], v_ρs[1:niter], Qs[1:niter]
end

# Special method of tryls for tooshortarc
function tryls(radec::Vector{RadecMPC{T}}, jd0::T, q0::Vector{T}, params::Parameters{T};
               maxiter::Int = 5) where {T <: AbstractFloat}
    # Allocate memory
    sols = zeros(NEOSolution{T, T}, maxiter)
    # Origin
    x0 = zeros(T, 6)
    # Scaling factors
    scalings = abs.(q0) ./ 10^5
    # Jet transport variables
    dq = scaled_variables("dx", scalings, order = 6)
    # Number of iterations
    niter = 1
    # Minimization loop
    while niter <= maxiter
        # Initial conditions
        q = q0 + dq
        # Propagation & residuals
        bwd, fwd, res = propres(radec, jd0, q, params)
        # Least squares fit
        fit = tryls(res, x0, 5)
        # Case: unsuccessful fit
        if !fit.success || any(diag(fit.Γ) .< 0)
            break
        end
        # Current solution
        sols[niter] = evalfit(NEOSolution(bwd, fwd, res, fit))
        # Convergence condition
        if niter > 1 && nrms(sols[niter-1]) - nrms(sols[niter]) < 0.1
            break
        end
        # Update values
        q0 = q(fit.x)
        # Update number of iterations
        niter += 1
    end

    return sols[niter]

end

@doc raw"""
    tooshortarc(radec::Vector{RadecMPC{T}}, gdf::GroupedDataFrame, cdf::DataFrame,
                params::Parameters{T}; kwargs...) where {T <: AbstractFloat}

Return initial conditions by minimizing the normalized root mean square error over 
the admissible region.

# Arguments

- `radec::Vector{RadecMPC{T}}`: vector of observations.
- `gdf::GroupedDataFrame`, `cdf::DataFrame`: output of [`reduce_nights`](@ref).
- `params::Parameters{T}`: see [`Parameters`](@ref).

# Keyword arguments

- `maxiter::Int = 200`:  maximum number of iterations.
"""
function tooshortarc(radec::Vector{RadecMPC{T}}, gdf::GroupedDataFrame, cdf::DataFrame,
                     params::Parameters{T}; maxiter::Int = 200) where {T <: AbstractFloat}
    
    # Reverse temporal order
    reverse!(cdf)
    # Sort nights by number of observations 
    idxs = sortperm(cdf.nobs, rev = true)
    # Allocate memory for output
    sol = zero(NEOSolution{T, T})

    for i in idxs
        # Admissible region
        A = AdmissibleRegion(cdf[i, :])
        # Center
        ρ = (R_SI + max_range(A)) / 2
        v_ρ = sum(range_rate(A, R_SI)) / 2
        # Minimization over admissible region
        ρs, v_ρs, Qs = adam(radec, A, ρ, v_ρ, params; maxiter)
        # Barycentric initial conditions
        q0 = topo2bary(A, ρs[end], v_ρs[end])
        # 6 variables least squares
        sol = tryls(radec, datetime2julian(A.date), q0, params; maxiter = 5)

        if nrms(sol) < 1.5
            return sol
        end
    end
    
    return sol
end