@doc raw"""
    topocentric(α::T, δ::T) where {T <: Number}
    topocentric(obs::RadecMPC{T}) where {T <: AbstractFloat}

Returns the object's topocentric position unit vector. Both angles must be in rad.
"""
function topocentric(α::T, δ::T) where {T <: Number}

    sin_α, cos_α = sincos(α)
    sin_δ, cos_δ = sincos(δ)

    return [cos_δ*cos_α, cos_δ*sin_α, sin_δ]
end
topocentric(obs::RadecMPC{T}) where {T <: AbstractFloat} = topocentric(obs.α, obs.δ)

@doc raw"""
    barycentric(obs::RadecMPC{T}) where {T <: AbstractFloat}

Returns the observer's  barycentric `[x, y, z, v_x, v_y, v_z]` geometric "state" vector. 
"""
function barycentric(obs::RadecMPC{T}) where {T <: AbstractFloat}
    # Geocentric state vector of the observer  
    g_vec = kmsec2auday(geocentric(obs))
    # Barycentric state vector of the Earth
    G_vec = kmsec2auday(earth_pv(datetime2et(obs)))

    return G_vec + g_vec
end

function f_Lagrange(μ, r_2, τ)
    return 1 - μ*(τ^2)/(2*r_2^3)
end

function g_Lagrange(μ, r_2, τ)
    return τ - μ*(τ^3)/(6*r_2^3)
end

@doc raw"""
    gauss_method_core(obs::Vector{RadecMPC{T}}, μ::T = PlanetaryEphemeris.GMS; root_idx::Int = 1) where {T <: AbstractFloat}

Gauss method of Initial Orbit determination (IOD).

# Arguments 

- `obs::Vector{RadecMPC{T}}`: three observations.
- `μ::T = PE.GMS`: gravitational parameter of center of attraction.
- `root_idx::Int = 1`: index of Lagrange equation root. 
"""
function gauss_method_core(obs::Vector{RadecMPC{T}}, μ::T = PlanetaryEphemeris.GMS; root_idx::Int = 1) where {T <: AbstractFloat}
   
    # Check we have exactly three observations 
    m = length(obs)
    @assert m == 3 "Core Gauss method requires three observations (got $m)"
    # Make sure observations are in temporal order 
    sort!(obs)

    println("***Gauss method of orbit determination***")

    # Julian dates of observation 
    t_1 = datetime2julian(obs[1].date)
    t_2 = datetime2julian(obs[2].date)
    t_3 = datetime2julian(obs[3].date)

    # Object topocentric position unit vectors 
    ρ_1_vec = topocentric(obs[1])
    ρ_2_vec = topocentric(obs[2])
    ρ_3_vec = topocentric(obs[3])

    # Observer barycentric positions 
    R_1_vec = barycentric(obs[1])[1:3]
    R_2_vec = barycentric(obs[2])[1:3]
    R_3_vec = barycentric(obs[3])[1:3]
    R = zeros(3, 3)
    R[1, :] = R_1_vec
    R[2, :] = R_2_vec
    R[3, :] = R_3_vec
    
    # Time intervals 
    τ_1 = t_1 - t_2
    τ_3 = t_3 - t_2
    τ = τ_3 - τ_1

    # Cross products 
    p = zeros(3, 3)
    p[1, :] = cross(ρ_2_vec, ρ_3_vec)
    p[2, :] = cross(ρ_1_vec, ρ_3_vec)
    p[3, :] = cross(ρ_1_vec, ρ_2_vec) 

    # Gauss scalar 
    D_0 = dot(ρ_1_vec, p[1, :])

    # Matrix of triple products 
    D = zeros(3, 3)
    for i in 1:3
        for j in 1:3
            D[i, j] = dot(R[i, :], p[j, :]) 
        end
    end 

    # A and B scalars 
    A = (-D[1, 2]*τ_3/τ + D[2, 2] + D[3, 2]*τ_1/τ) / D_0
    B = (D[1, 2]*(τ_3^2 - τ^2)*τ_3/τ + D[3, 2]*(τ^2 - τ_1^2)*τ_1/τ) / 6 / D_0

    # E and F scalars 
    E = dot(R_2_vec, ρ_2_vec)
    F = dot(R_2_vec, R_2_vec)

    # Lagrange equation coefficients
    a = -(A^2 + 2*A*E + F)
    b = -2*μ*B*(A + E)
    c = -(μ^2)*(B^2)

    # Find roots of Lagrange equation 
    sol = find_zeros(x -> x^8 + a*x^6 + b*x^3 + c, 0, 1_000)
    n_sol = length(sol)

    if n_sol == 0
        r_2 = 1.
        @warn("""No solutions found for Lagrange equation r^8 + $a r^6 + $b r^3 + $c = 0; 
        Setting default to 1.""")
    elseif n_sol == 1
        r_2 = sol[1]
        println("Lagrange equation root: ", r_2)
    else 
        r_2 = sol[root_idx]
        @warn("""More than one solution $sol found for Lagrange equation r^8 + $a r^6 + $b r^3 + $c = 0;
        Selecting root with index $root_idx""")
    end 

    # Slant ranges 
    num_1 = 6*(D[3, 1]*τ_1/τ_3 + D[2, 1]*τ/τ_3)*(r_2^3) + μ*D[3, 1]*(τ^2 - τ_1^2)*τ_1/τ_3
    den_1 = 6*(r_2^3) + μ*(τ^2 - τ_3^2)
    ρ_1 = (num_1 / den_1 - D[1, 1]) / D_0

    ρ_2 = A + μ*B/(r_2^3)

    num_3 = 6*(D[1, 3]*τ_3/τ_1 - D[2, 3]*τ/τ_1)*(r_2^3) + μ*D[1, 3]*(τ^2 - τ_3^2)*τ_3/τ_1
    den_3 = 6*(r_2^3) + μ*(τ^2 - τ_1^2)
    ρ_3 = (num_3 / den_3 - D[3, 3]) / D_0

    # Barycentric position of the body 
    r_1_vec = R_1_vec + ρ_1 * ρ_1_vec
    r_2_vec = R_2_vec + ρ_2 * ρ_2_vec
    r_3_vec = R_3_vec + ρ_3 * ρ_3_vec

    # f, g Lagrange coefficients
    f_1 = f_Lagrange(μ, r_2, τ_1)
    f_3 = f_Lagrange(μ, r_2, τ_3)

    g_1 = g_Lagrange(μ, r_2, τ_1)
    g_3 = g_Lagrange(μ, r_2, τ_3)

    # Barycentric velocity of the body 
    v_2_vec = (-f_3*r_1_vec + f_1*r_3_vec)/(f_1*g_3 - f_3*g_1)

    return r_1_vec, r_2_vec, r_3_vec, v_2_vec, D, ρ_1_vec, ρ_2_vec, ρ_3_vec, τ_1, τ_3, f_1, g_1, f_3, g_3, ρ_1, ρ_2, ρ_3 
end