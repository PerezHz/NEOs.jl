@auto_hash_equals struct GaussSolution{T <: Real, U <: Number}
    r_2_vec::Vector{U}
    v_2_vec::Vector{U}
    D::Matrix{U}
    R_vec::Matrix{T}
    ρ_vec::Matrix{U}
    τ_1::T 
    τ_3::T
    f_1::U
    g_1::U
    f_3::U 
    g_3::U 
    success::Bool
    # Internal constructor 
    function GaussSolution{T, U}(r_2_vec::Vector{U}, v_2_vec::Vector{U}, D::Matrix{U}, R_vec::Matrix{T}, ρ_vec::Matrix{U},
                                 τ_1::T, τ_3::T, f_1::U, g_1::U, f_3::U, g_3::U, success::Bool) where {T <: Real, U <: Number}
        return new{T, U}(r_2_vec, v_2_vec, D, R_vec, ρ_vec, τ_1, τ_3, f_1, g_1, f_3, g_3, success)
    end 
end 
# Outer constructor 
function GaussSolution(r_2_vec::Matrix{U}, v_2_vec::Vector{U}, D::Matrix{U}, R_vec::Matrix{T}, ρ_vec::Matrix{U}, 
                       τ_1::T, τ_3::T, f_1::U, g_1::U, f_3::U, g_3::U, success::Bool) where {T <: Real, U <: Number}
    return GaussSolution(r_2_vec, v_2_vec, D, R_vec, ρ_vec, τ_1, τ_3, f_1, g_1, f_3, g_3, success)
end 

function isless(a::GaussSolution{T, U}, b::GaussSolution{T, U}) where {T <: Real, U <: Number} 
    return norm(cte.(a.r_2_vec)) < norm(cte.(b.r_2_vec))
end 

# Print method for GaussSolution
# Examples: 
# 
function show(io::IO, g::GaussSolution{T, U}) where {T <: Real, U <: Number} 
    print(io, "Gauss solution (r = ", norm(cte.(g.r_2_vec)), ")")
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

@doc raw"""
    f_Lagrange(r, τ)

Return the 1st order approximation to Lagrange's f function. 
"""
function f_Lagrange(τ::T, r::U) where {T <: Real, U <: Number}
    return 1 - μ_S * (τ^2) / 2 / (r^3)
end

@doc raw"""
    g_Lagrange(r, τ)  
    
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

function _format_Lagrange_equation(a, b, c)
    a_sgn = a ≥ 0 ? "+" : "-"
    b_sgn = b ≥ 0 ? "+" : "-"
    c_sgn = c ≥ 0 ? "+" : "-"

    return join(["r⁸ ", a_sgn, " ", abs(a), " r⁶ ", b_sgn, " ", abs(b), " r³ ", c_sgn, " ", abs(c), " = 0"])
end 

function _format_solutions(sol)
    
    s = Vector{String}(undef, length(sol))
    for i in eachindex(sol)
        s[i] = join([string(sol[i].status), " root at ", string(sol[i].interval), ", "])
    end 

    return join(["{ ", s..., " }"])
end 

lagrange(x, a, b, c) = x^8 + a*x^6 + b*x^3 + c
lagrange_derivative(x, a, b) = 8*x^7 + 6*a*x^5 + 3*b*x^2

function solve_lagrange(a::T, b::T, c::T) where {T <: Real}

    sol = roots(x -> lagrange(x, a, b, c), Interval(0.00465047, 40))
    return mid.(interval.(sol)) 

end 

function solve_lagrange(a::TaylorN{T}, b::TaylorN{T}, c::TaylorN{T}) where {T <: Real}

    sol_0_ = roots(x -> lagrange(x, cte(a), cte(b), cte(c)), Interval(0.00465047, 40))
    sol_0 = mid.(interval.(sol_0_))

    sol = Vector{TaylorN{T}}(undef, length(sol_0))

    for i in eachindex(sol_0)
        r_0 = sol_0[i] 
        r_2 = sol_0[i] 
        for i in 1:5
            r_2 = r_0 - lagrange(r_0, a, b, c) / lagrange_derivative(r_0, a, b)
            r_0 = r_2
        end 
        sol[i] = r_2
    end

    return sol 
end 

@doc raw"""
    gauss_method(obs::Vector{RadecMPC{T}}; root_idx::Int = 1) where {T <: AbstractFloat}
    gauss_method(obs::Vector{RadecMPC{T}}; root_idx::Int = 1) where {T <: AbstractFloat}

Core Gauss method of Initial Orbit determination (IOD). See Algorithm 5.5 in page 274 https://doi.org/10.1016/C2016-0-02107-1.

# Arguments 

- `obs::Vector{RadecMPC{T}}`: three observations.
- `root_idx::Int = 1`: index of Lagrange equation solution in case of multiple roots. 
"""
function gauss_method(obs::Vector{RadecMPC{T}}; xve::Function = et -> kmsec2auday(getpv(399, 10, et)), niter::Int = 10) where {T <: AbstractFloat}
    
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

    return gauss_method(observatories, dates, α, δ; xve = xve, niter = niter)
end 

function gauss_method(observatories::Vector{ObservatoryMPC{T}}, dates::Vector{DateTime}, α::Vector{U}, δ::Vector{U};
                      xve::Function = et -> kmsec2auday(getpv(399, 10, et)), niter::Int = 10) where {T <: Real, U <: Number}

    # Check we have exactly three observations 
    @assert length(observatories) == length(dates) == length(α) == length(δ) == 3 "Gauss method requires exactly three observations"

    # Julian days of observation
    t_julian = datetime2julian.(dates)

    # Time intervals 
    τ_1 = t_julian[1] - t_julian[2]
    τ_3 = t_julian[3] - t_julian[2]
    τ = τ_3 - τ_1

    # NEO's topocentric position unit vectors 
    ρ_vec = vectors2matrix(topounit.(α, δ))

    # Times of observation [et]
    t_et = datetime2et.(dates)

    # Geocentric state vector of the observer [au, au/day] 
    g_vec = kmsec2auday.(obs_pv_ECI.(observatories, t_et))

    # Heliocentric state vector of the Earth [au, au/day]
    G_vec = xve.(t_et)

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
    sol = solve_lagrange(a, b, c)

    # Number of solutions 
    n_sol = length(sol)

    if n_sol == 0
        
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

            sol_gauss[i] = GaussSolution{T, U}(r_vec[2, :], v_2_vec, D, R_vec, ρ_vec, τ_1, τ_3, f_1, g_1, f_3, g_3, true)
        end 

        return sort!(sol_gauss)

    end 
    
end

@doc raw"""
    stumpC(z::T) where {T <: Number}

Evaluate Stumpff function C(z). 
"""
function stumpC(z::T) where {T <: Number}
    if z > 0
        c = (1 - cos(sqrt(z))) / z
    elseif z < 0
        c = (cosh(sqrt(-z)) - 1) / (-z)
    else 
        c = one(T) / 2
    end

    return c
end

@doc raw"""
    stumpS(z::T) where {T <: Number}

Evaluate Stumpff function S(z). 
"""
function stumpS(z::T) where {T <: Number}
    if z > 0
        s = (sqrt(z) - sin(sqrt(z))) / sqrt(z)^3
    elseif z < 0
        s = (sinh(sqrt(-z)) - sqrt(-z)) / sqrt(-z)^3
    else 
        s = one(T) / 6
    end

    return s
end

@doc raw"""
    _f_Lagrange(χ, z, r)

Return the current value of Lagrange's f function. 
"""
function _f_Lagrange(χ, z, r)
    return 1 - (χ^2) * stumpC(z) / r
end

@doc raw"""
    _g_Lagrange(τ, χ, z)

Return the current value of Lagrange's g function. 
"""
function _g_Lagrange(τ, χ, z)
    return τ - (χ^3) * stumpS(z) / sqrt(μ_S)
end

@doc raw"""
    univkepler(τ::T, r_2_vec::Vector{T}, v_2_vec::Vector{T}; kep_iters::Int = 10, atol::T = 3e-14) where {T <: AbstractFloat}

Solve the universal Kepler equation for the universal anomaly.

# Arguments 

- `τ::T`: time interval. 
- `r_2_vec::Vector{T}`: heliocentric position vector.
- `v_2_vec::Vector{T}`: heliocentric velocity vector.
- `kep_iters::Int`: number of iterations for Newton's method. 
- `atol::T`: absolute tolerance. 
"""
function univkepler(τ::T, r_2_vec::Vector{T}, v_2_vec::Vector{T}; kep_iters::Int = 10, atol::T = 3e-14) where {T <: AbstractFloat}
    
    # Slant range 
    r_2 = norm(r_2_vec)
    # Velocity squared 
    v2_2 = dot(v_2_vec, v_2_vec)
    # Radial velocity
    v_r = dot(r_2_vec, v_2_vec) / r_2
    # Reciprocal of semimajor axis 
    α = (2 / r_2) - (v2_2 / μ_S)

    # Initial estimate for χ 
    χ = sqrt(μ_S)*abs(α)*τ

    # Newton iteration 
    n = 0
    ratio = 1.

    while (abs(ratio) > atol) && (n <= kep_iters)

        # One more iteration 
        n += 1

        χ2 = χ^2
        z = α*χ2

        a = r_2 * v_r / sqrt(μ_S)
        b = 1 - α * r_2

        C = stumpC(z)
        S = stumpS(z)

        if isinf(C) || isinf(S)
            return NaN
        end 

        # Function 
        F = a * χ2 * C + b * (χ^3) * S + r_2 * χ - sqrt(μ_S) * τ
        # Derivative 
        dFdx = a * χ * (1 - z * S) + b * χ2 * C + r_2
        # Newton quotient 
        ratio = F / dFdx
        # Newton update rule 
        χ = χ - ratio

    end

    return χ
end

@doc raw"""
    gauss_method_refinement(τ_1::T, τ_3::T, r_2_vec::Vector{T}, v_2_vec::Vector{T}, D::Matrix{T}, R_vec::Matrix{T}, ρ_vec::Matrix{T}, 
                            f_1::T, g_1::T, f_3::T, g_3::T; kep_iters::Int = 10, atol::T = 3e-14) where {T <: AbstractFloat}

Iterative improvement of the orbit determined by core Gauss method. 

# Arguments 

- `τ_1/τ_3::T`: time intervals. 
- `r_2_vec::Vector{T}`: heliocentric position vector.
- `v_2_vec::Vector{T}`: heliocentric velocity vector.
- `D::Matrix{T}`: matrix of Gauss scalars. 
- `R_vec::Matrix{T}`: matrix of observer heliocentric position vectors. 
- `ρ_vec::Matrix{T}`: matrix of NEO topocentric position unit vectors. 
- `f_1/g_1/f_3/g_3::T`: initial values of Lagrange's f and g function. 
- `kep_iters::Int`: number of iterations for Newton's method. 
- `atol::T`: absolute tolerance. 
"""
function gauss_method_refinement(τ_1::T, τ_3::T, r_2_vec::Vector{T}, v_2_vec::Vector{T}, D::Matrix{T}, R_vec::Matrix{T}, 
                                 ρ_vec::Matrix{T}, f_1::T, g_1::T, f_3::T, g_3::T; kep_iters::Int = 10, atol::T = 3e-14) where {T <: AbstractFloat}
    
    # Refinement sucess 
    success = true

    # Solve universal Kepler problem 
    χ_1 = univkepler(τ_1, r_2_vec, v_2_vec; kep_iters = kep_iters, atol = atol)
    χ_3 = univkepler(τ_3, r_2_vec, v_2_vec; kep_iters = kep_iters, atol = atol)

    if isnan(χ_1) || isnan(χ_3)
        success = false 
        return fill(NaN, 3, 3), v_2_vec, fill(NaN, 3, 3), f_1, g_1, f_3, g_3, success
    end

    # Slant range 
    r_2 = norm(r_2_vec)
    # Velocity squared
    v2_2 = dot(v_2_vec, v_2_vec)
    # Reciprocal of semimajor axis 
    α = (2 / r_2) - (v2_2 / μ_S)

    # New values of Lagrange's f, g functions 
    z_1 = α * χ_1^2 
    f_1_new = (f_1 + _f_Lagrange(χ_1, z_1, r_2)) / 2
    g_1_new = (g_1 + _g_Lagrange(τ_1, χ_1, z_1)) / 2

    z_3 = α * χ_3^2 
    f_3_new = (f_3 + _f_Lagrange(χ_3, z_3, r_2)) / 2
    g_3_new = (g_3 + _g_Lagrange(τ_3, χ_3, z_3)) / 2

    # Denominator 
    den = f_1_new * g_3_new - f_3_new * g_1_new

    if isinf(den)
        success = false
        return fill(NaN, 3, 3), v_2_vec, fill(NaN, 3, 3), f_1, g_1, f_3, g_3, success
    end 

    # c coefficients
    c_1 = g_3_new / den
    c_3 = - g_1_new / den

    # Gauss scalar
    D_0 = dot(ρ_vec[1, :], cross(ρ_vec[2, :], ρ_vec[3, :]))

    # Updated slant ranges 
    ρ = zeros(3)
    ρ[1] = (-D[1, 1] + D[2, 1]/c_1 - D[3, 1]*c_3/c_1) / D_0
    ρ[2] = (-c_1*D[1, 2] + D[2, 2] - c_3*D[3, 2]) / D_0
    ρ[3] = (-D[1, 3]*c_1/c_3 + D[2, 3]/c_3 - D[3, 3]) / D_0

    # Updated NEO's heliocentric position
    r_vec = R_vec .+ ρ.*ρ_vec
    
    # Updated NEO's barycentric velocity
    v_2_vec = (-f_3_new * r_vec[1, :] + f_1_new * r_vec[3, :]) / den

    return r_vec, v_2_vec, ρ, f_1_new, g_1_new, f_3_new, g_3_new, success

end

@doc raw"""
    gauss_method_iterator(obs::Vector{RadecMPC{T}}; kep_iters::Int = 10, atol = 3e-14, ref_iters::Int = 10, 
                          root_idx::Int = 1) where {T <: AbstractFloat}

Iterate over consecutive triplets in `obs` and apply core Gauss method to each one. 

# Arguments 

- `obs::Vector{RadecMPC{T}}`: observations. 
- `kep_iters::Int`: number of iterations for Newton's method. 
- `atol::T`: absolute tolerance. 
- `ref_iters::Int`:: number of refinement iterations.
- `root_idx::Int`: index of Lagrange equation solution in case of multiple roots. 
"""
function gauss_method_iterator(obs::Vector{RadecMPC{T}}; kep_iters::Int = 10, atol = 3e-14, ref_iters::Int = 10, 
                               root_idx::Int = 1) where {T <: AbstractFloat}

    # Core Gauss method 
    r_vec, v_2_vec, D, R_vec, ρ_vec, τ_1, τ_3, f_1, g_1, f_3, g_3, ρ, t, success = gauss_method_core(obs; root_idx = root_idx)

    # Copy core valus as backup
    r_vec0 = copy(r_vec)
    v_2_vec0 = copy(v_2_vec)
    ρ0 = copy(ρ)

    if !success
        @warn "Unsuccessful Core Gauss method"

        return r_vec0, v_2_vec0, ρ0, t, success
    end 

    # Gauss refinment 
    for i in 1:ref_iters

        r_vec, v_2_vec, ρ, f_1, g_1, f_3, g_3, success = 
            gauss_method_refinement(τ_1, τ_3, r_vec[2, :], v_2_vec, D, R_vec, ρ_vec, f_1, g_1, f_3, g_3; kep_iters = kep_iters, atol = atol)

        if !success
            @warn "Unsuccessful refinement"
            return r_vec0, v_2_vec0, ρ0, t, success
        end

    end 

    return r_vec, v_2_vec, ρ, t, success
end 