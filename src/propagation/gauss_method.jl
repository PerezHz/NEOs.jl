@doc raw"""
    neo_pos_topo(obs::RadecMPC{T}) where {T <: AbstractFloat}

Returns the NEO's topocentric position unit vector.
"""
function neo_pos_topo(obs::RadecMPC{T}) where {T <: AbstractFloat}
    
    sin_α, cos_α = sincos(obs.α)
    sin_δ, cos_δ = sincos(obs.δ)

    pos = [cos_δ * cos_α, cos_δ * sin_α, sin_δ]

    return pos
end

@doc raw"""
    obs_pv_hel(obs::RadecMPC{T}) where {T <: AbstractFloat}

Returns the observer's  heliocentric `[x, y, z, v_x, v_y, v_z]` state vector in units of au, au/day. 
"""
function obs_pv_hel(obs::RadecMPC{T}) where {T <: AbstractFloat}

    # Geocentric state vector of the observer [km, km/s] 
    g_vec = obs_pv_ECI(obs)

    # ET secs 
    et = datetime2et(obs)
    # Barycentric state vector of the Sun [km, km/s]
    S_vec = sun_pv(et)
    # Barycentric state vector of the Earth [km, km/s]
    E_vec = earth_pv(et)

    # Heliocentric state vector of the Earth [km, km/s]
    G_vec = E_vec - S_vec

    # Heliocentric state vector of the observer [au, au/day]
    return kmsec2auday(G_vec + g_vec)
end

@doc raw"""
    f_Lagrange(r, τ)

Returns the 1st order approximation to Lagrange's f function. 
"""
function f_Lagrange(r, τ)
    return 1 - μ_S * (τ^2) / 2 / (r^3)
end

@doc raw"""
    g_Lagrange(r, τ)  
    
Returns the 1st order approximation to Lagrange's g function.
"""
function g_Lagrange(r, τ)
    return τ - μ_S * (τ^3) / 6 / (r^3)
end

function vectors2matrix(x::T) where {T <: AbstractVector}
    return permutedims(reduce(hcat, x))
end

@doc raw"""
    gauss_method_core(obs::Vector{RadecMPC{T}}; root_idx::Int = 1) where {T <: AbstractFloat}

Core Gauss method of Initial Orbit determination (IOD). See Algorithm 5.5 in page 274 https://doi.org/10.1016/C2016-0-02107-1.

# Arguments 

- `obs::Vector{RadecMPC{T}}`: three observations.
- `root_idx::Int = 1`: index of Lagrange equation solution in case of multiple roots. 
"""
function gauss_method_core(obs::Vector{RadecMPC{T}}; root_idx::Int = 1) where {T <: AbstractFloat}

    # Check we have exactly three observations 
    m = length(obs)
    @assert m == 3 "Core Gauss method requires exactly three observations, got $m"

    # Sucess of gauss core method 
    success = true 

    # Make sure observations are in temporal order 
    sort!(obs)

    # Julian dates of observation 
    t = datetime2julian.(date.(obs))

    # Time intervals 
    τ_1 = t[1] - t[2]
    τ_3 = t[3] - t[2]
    τ = τ_3 - τ_1

    # NEO's topocentric position unit vectors 
    ρ_vec = vectors2matrix(neo_pos_topo.(obs))

    # Observer's heliocentric positions 
    R_vec = vectors2matrix(obs_pv_hel.(obs))[:, 1:3]

    # Cross products 
    p_vec = zeros(3, 3)
    p_vec[1, :] = cross(ρ_vec[2, :], ρ_vec[3, :])
    p_vec[2, :] = cross(ρ_vec[1, :], ρ_vec[3, :])
    p_vec[3, :] = cross(ρ_vec[1, :], ρ_vec[2, :]) 

    # Gauss scalar 
    D_0 = dot(ρ_vec[1, :], p_vec[1, :])

    # Matrix of triple products 
    D = zeros(3, 3)
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

    # Find roots of Lagrange equation 
    sol = find_zeros(x -> x^8 + a*x^6 + b*x^3 + c, 0, 1_000)
    n_sol = length(sol)

    if n_sol == 0

        success = false 

        @warn("""No solutions found for Lagrange equation r^8 + $a r^6 + $b r^3 + $c = 0; 
        Cannot procede""")

        return fill(NaN, 3, 3), fill(NaN, 3), D, R_vec, ρ_vec, τ_1, τ_3, NaN, NaN, NaN, NaN, fill(NaN, 3), t, success

    elseif n_sol == 1
        
        r_2 = sol[1]

    else 

        r_2 = sol[root_idx]

        @warn("""More than one solution $sol found for Lagrange equation r^8 + $a r^6 + $b r^3 + $c = 0;
        Selecting root with index $root_idx""")

    end 

    # Slant ranges 
    ρ = zeros(3)

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
    f_1 = f_Lagrange(r_2, τ_1)
    f_3 = f_Lagrange(r_2, τ_3)

    g_1 = g_Lagrange(r_2, τ_1)
    g_3 = g_Lagrange(r_2, τ_3)

    # Heliocentric velocity of the NEO 
    v_2_vec = (- f_3 * r_vec[1, :] + f_1 * r_vec[3, :]) / (f_1*g_3 - f_3*g_1)

    return r_vec, v_2_vec, D, R_vec, ρ_vec, τ_1, τ_3, f_1, g_1, f_3, g_3, ρ, t, success
end

@doc raw"""
    stumpC(z::T) where {T <: Real}

Evaluates Stumpff function C(z). 
"""
function stumpC(z::T) where {T <: Real}
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
    stumpS(z::T) where {T <: Real}

Evaluates Stumpff function S(z). 
"""
function stumpS(z::T) where {T <: Real}
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

Returns the current value of Lagrange's f function. 
"""
function _f_Lagrange(χ, z, r)
    return 1 - (χ^2) * stumpC(z) / r
end

@doc raw"""
    _g_Lagrange(τ, χ, z)

Returns the current value of Lagrange's g function. 
"""
function _g_Lagrange(τ, χ, z)
    return τ - (χ^3) * stumpS(z) / sqrt(μ_S)
end

@doc raw"""
    univkepler(τ::T, r_2_vec::Vector{T}, v_2_vec::Vector{T}; kep_iters::Int = 10, atol::T = 3e-14) where {T <: AbstractFloat}

Solves the universal Kepler equation for the universal anomaly.

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

function gauss_method(obs::Vector{RadecMPC{T}}; kep_iters::Int = 10, atol = 3e-14, ref_iters::Int = 10, 
                      root_idx::Int = 1) where {T <: AbstractFloat}
    
    # Number of observations 
    m = length(obs)

    # Vector of osculating orbita elements 
    osc = Vector{OsculatingElements{T}}(undef, m-2)
    
    # Iterate over consecutive triplets 
    for j in 1:m-2
        # Current triplet
        obs_ = obs[[j, j+1, j+2]]
        # Gauss method 
        r_vec, v_2_vec, ρ, t, success = gauss_method_iterator(obs_; kep_iters = kep_iters, atol = atol, ref_iters = ref_iters,
                                                              root_idx = root_idx)

        if !success
            osc[j] = OsculatingElements()
        end

        # Correct times for light propagation                                                              
        for k in eachindex(t)
            t[k] = t[k] - ρ[k] / c_au_per_day
        end
        # Try computing orbital elements 
        try
            osc[j] = pv2kep(vcat(r_vec[2, :], v_2_vec), μ_S, t[2])
        catch
            @warn "Cannot compute orbital elements for state vector $(vcat(r_vec[2, :], v_2_vec)) at epoch $(t[2])"
            osc[j] = OsculatingElements()
        end

    end

    # Eliminate NaN orbital elements 
    filter!(!isnan, osc)

    if length(osc) == 0
        @warn("Cannot compute orbital elements for any triplet.")
        return OsculatingElements()
    end 

    # Average orbital elements 
    return mean(osc)

end