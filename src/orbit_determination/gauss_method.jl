include("osculating.jl")

@doc raw"""
    GaussSolution{T <: Real, U <: Number}

A preliminary orbit obtained from Gauss method of orbit determination. See Algorithm 5.5 in page 274 of
https://doi.org/10.1016/C2016-0-02107-1.

See also [`gauss_method`](@ref).

# Fields

- `statevect::Vector{U}`: state vector at middle observation.
- `D::Matrix{U}`: D matrix.
- `R_vec::Matrix{T}`: observer's heliocentric positions.
- `ρ_vec::Matrix{U}`: slant ranges.
- `τ_1::T`: time between first and second observations.
- `τ_3::T`: time between third and second observations.
- `f_1, g_1, f_3, g_3::U`: Lagrange coefficients.
- `status::Symbol`: the status of the solution (`:empty`, `:unkown` or `:unique`)
"""
@auto_hash_equals struct GaussSolution{T <: Real, U <: Number}
    statevect::Vector{U}
    D::Matrix{U}
    R_vec::Matrix{T}
    ρ_vec::Matrix{U}
    τ_1::T
    τ_3::T
    f_1::U
    g_1::U
    f_3::U
    g_3::U
    status::Symbol
    # Internal constructor
    function GaussSolution{T, U}(statevect::Vector{U}, D::Matrix{U}, R_vec::Matrix{T}, ρ_vec::Matrix{U}, τ_1::T, τ_3::T,
                                 f_1::U, g_1::U, f_3::U, g_3::U, status::Symbol) where {T <: Real, U <: Number}
        return new{T, U}(statevect, D, R_vec, ρ_vec, τ_1, τ_3, f_1, g_1, f_3, g_3, status)
    end
end

# Outer constructor
function GaussSolution(statevect::Vector{U}, D::Matrix{U}, R_vec::Matrix{T}, ρ_vec::Matrix{U}, τ_1::T, τ_3::T,
                       f_1::U, g_1::U, f_3::U, g_3::U, status::Symbol) where {T <: Real, U <: Number}
    return GaussSolution{T, U}(statevect, D, R_vec, ρ_vec, τ_1, τ_3, f_1, g_1, f_3, g_3, status)
end

# Print method for GaussSolution
# Examples:
# unique Gauss solution (r = 1.0800950907383229)
function show(io::IO, g::GaussSolution{T, U}) where {T <: Real, U <: Number}
    print(io, g.status, " Gauss solution (r = ", norm(cte.(g.statevect[1:3])), ")")
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

# TO DO: Should we allow to use other μ?

@doc raw"""
    f_Lagrange(τ::T, r::U) where {T <: Real, U <: Number}

Return the 1st order approximation to Lagrange's f function.
"""
function f_Lagrange(τ::T, r::U) where {T <: Real, U <: Number}
    return 1 - μ_S * (τ^2) / 2 / (r^3)
end

@doc raw"""
    g_Lagrange(τ::T, r::U) where {T <: Real, U <: Number}

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

@doc raw"""
    _format_Lagrange_equation(a::T, b::T, c::T) where {T <: AbstractFloat}

Format Lagrange equation as `r⁸ + a r⁶ + b r³ + c = 0`.
"""
function _format_Lagrange_equation(a::T, b::T, c::T) where {T <: AbstractFloat}
    a_sgn = a ≥ 0 ? "+" : "-"
    b_sgn = b ≥ 0 ? "+" : "-"
    c_sgn = c ≥ 0 ? "+" : "-"

    return join(["r⁸ ", a_sgn, " ", abs(a), " r⁶ ", b_sgn, " ", abs(b), " r³ ", c_sgn, " ", abs(c), " = 0"])
end

@doc raw"""
    lagrange(x::T, a::U, b::U, c::U) where {T, U <: Number}

Lagrange polynomial to be solved during Gauss method.
"""
lagrange(x::T, a::U, b::U, c::U) where {T, U <: Number} = x^8 + a*x^6 + b*x^3 + c

@doc raw"""
    lagrange_derivative(x::T, a::T, b::T) where {T <: Number}

Derivative of Lagrange polynomial to be solved during Gauss method.
"""
lagrange_derivative(x::T, a::U, b::U) where {T, U <: Number} = 8*x^7 + 6*a*x^5 + 3*b*x^2

# TO DO: Allow to control interval over which to look for solutions
# Currently we look between the radius of the Sun ((0.00465047 au) and the radius of the Solar System (40 au)

@doc raw"""
    solve_lagrange(a::T, b::T, c::T; niter::Int = niter) where {T <: Real}
    solve_lagrange(a::TaylorN{T}, b::TaylorN{T}, c::TaylorN{T}; niter::Int = niter) where {T <: Real}

Solve Lagrange polynomial.

See also [`lagrange`](@ref).
"""
function solve_lagrange(a::T, b::T, c::T; niter::Int = 5) where {T <: Real}
    sol = roots(x -> lagrange(x, a, b, c), Interval(0.00465047, 40))
    return mid.(interval.(sol)), getfield.(sol, :status)
end

function solve_lagrange(a::TaylorN{T}, b::TaylorN{T}, c::TaylorN{T}; niter::Int = 5) where {T <: Real}
    # 0-th order solution
    sol_0, status = solve_lagrange(cte(a), cte(b), cte(c))
    # Discard empty solutions
    idxs = findall(x -> x != :empty, status)
    # Vector of solutions
    sol = Vector{TaylorN{T}}(undef, length(idxs))
    # Iterate non-empty solutions
    for i in eachindex(sol)
        # Newton's method
        r_0 = sol_0[idxs[i]]
        r_2 = sol_0[idxs[i]]
        for j in 1:niter
            r_2 = r_0 - lagrange(r_0, a, b, c) / lagrange_derivative(r_0, a, b)
            r_0 = r_2
        end
        sol[i] = r_2
    end

    return sol, status[idxs]
end

@doc raw"""
    gauss_method(obs::Vector{RadecMPC{T}}; xve::EarthEph = et -> kmsec2auday(getposvel(399, 10, et)),
                 niter::Int = 10) where {T <: AbstractFloat, EarthEph}
    gauss_method(observatories::Vector{ObservatoryMPC{T}}, dates::Vector{DateTime}, α::Vector{U}, δ::Vector{U};
                 xve::EarthEph = et -> kmsec2auday(getposvel(399, 10, et)), niter::Int = 5) where {T <: Real, U <: Number, EarthEph}

Core Gauss method of Initial Orbit determination (IOD). See Algorithm 5.5 in page 274 https://doi.org/10.1016/C2016-0-02107-1.

# Arguments

- `obs::Vector{RadecMPC{T}}`: three observations.
- `observatories::Vector{ObservatoryMPC{T}}`: sites of observation.
- `dates::Vector{DateTime}`: times of observation.
- `α::Vector{U}`: right ascension.
- `δ::Vector{U}`: declination.
- `xve::EarthEph`: Earth's ephemeris [et -> au, au/day].
- `niter::Int`: Number of iterations for Newton's method.
"""
function gauss_method(obs::Vector{RadecMPC{T}}; xve::EarthEph = et -> kmsec2auday(getposvel(399, 10, et)), niter::Int = 5) where {T <: AbstractFloat, EarthEph}

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
                      xve::EarthEph = et -> kmsec2auday(getposvel(399, 10, et)), niter::Int = 5) where {T <: Real, U <: Number, EarthEph}

    # Check we have exactly three observations
    @assert length(observatories) == length(dates) == length(α) == length(δ) == 3 "Gauss method requires exactly three observations"
    # Check observations are in temporal order 
    @assert issorted(dates)

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
    g_vec = kmsec2auday.(obsposvelECI.(observatories, t_et))

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
    sol, status = solve_lagrange(a, b, c; niter = niter)

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

            sol_gauss[i] = GaussSolution{T, U}(vcat(r_vec[2, :], v_2_vec), D, R_vec, ρ_vec, τ_1, τ_3, f_1, g_1, f_3, g_3, status[i])
        end

        return sort!(sol_gauss, by = x -> norm(x.statevect[1:3]))

    end

end

function vandermonde(x::Vector{T}, order::Int) where {T <: Real}
    # Number of points 
    L = length(x)
    # Initialize Vandermonde matrix 
    V = Matrix{T}(undef, L, order+1)
    # Fill first column
    V[:, 1] .= one(T)
    # Fill matrix 
    for j in 1:order
        for i in 1:L
            V[i, j+1] = x[i]^j
        end 
    end 
    return V
end 

function laplace_interpolation(x::Vector{T}, y::Vector{T}) where {T <: Real}
    # Check we have as many x as y 
    @assert length(x) == length(y)
    # Polynomial order
    order = length(x) - 1
    # Vandermonde matrix 
    V = vandermonde(x, order)
    # Solve the system of equations 
    if iszero(det(V))
        coeffs = fill(T(NaN), order+1)
    else 
        coeffs = V \ y
    end 
    # Return polynomial 
    return Taylor1{T}(coeffs, order)
end 

function laplace_interpolation(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
    # Julian days of observation
    t_julian = datetime2julian.(date.(radec))
    # Days of observation [relative to first observation]
    t_rel = t_julian .- t_julian[1]
    # Mean date 
    t_mean = sum(t_rel) / length(t_rel)

    # Interpolating polynomials 
    α_p = laplace_interpolation(t_rel, ra.(radec))
    δ_p = laplace_interpolation(t_rel, dec.(radec))
    # Evaluate polynomials at mean date 
    α_mean = α_p(t_mean)
    δ_mean = δ_p(t_mean)

    return julian2datetime(t_julian[1] + t_mean), α_mean, δ_mean
end 

# Empty method of reduce_nights, overloaded by QueryExt
function reduce_nights end 