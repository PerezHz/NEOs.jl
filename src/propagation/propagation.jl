include("asteroid_dynamical_models.jl")
include("jetcoeffs.jl")
include("integration_methods.jl")
include("serialization.jl")

@doc raw"""
    rvelea(dx, x, eph, params, t)

Return `true` and the asteroid's radial velocity with respect to the Earth.

# Arguments

- `dx`: asteroid's velocities.
- `x`: asteroid's degrees of freedom.
- `eph`: ephemeris.
- `params`: parameters (ephemeris + accelerations + newtonian N body potential + julian date of start time + matrix of extended body interactions + number of bodies + mass parameters).
- `t`: time.
"""
function rvelea(dx, x, eph, params, t)

    N = params[6]                                  # Number of bodies in the ephemeris
    xe = params[1][nbodyind(N-1,ea)]               # Earth's ephemeris

    return true, (x[1]-xe[1])*(x[4]-xe[4]) + (x[2]-xe[2])*(x[5]-xe[5]) + (x[3]-xe[3])*(x[6]-xe[6])
end

@doc raw"""
    loadeph(ss16asteph_::TaylorInterpolant, μ::Vector{T}) where {T <: Real}

Return the ephemeris in `ss16asteph_` with times converted from seconds to days,
the point-mass newtonian accelerations and the newtonian N body potential.

# Arguments

- `ss16asteph_`: solar system ephemeris.
- `μ::Vector{T}`: vector of mass parameters.
"""
function loadeph(ss16asteph_::TaylorInterpolant, μ::Vector{T}) where {T <: Real}

    # Read Solar System ephemeris (Sun + 8 planets + Moon + Pluto + 16 main belt asteroids)
    ss16asteph_t0 = T(ss16asteph_.t0 / daysec)    # Start time [days]
    ss16asteph_t = T.(ss16asteph_.t ./ daysec)      # Vector of times [days]
    ephord = ss16asteph_.x[1].order               # Order of the Taylor polynomials
    ss16asteph_x = map(x -> x(Taylor1(T, ephord)*daysec), ss16asteph_.x)          # Vector of Taylor polynomials [au, au/day]
    ss16asteph = TaylorInterpolant(ss16asteph_t0, ss16asteph_t, ss16asteph_x)  # TaylorInterpolant ephemeris  [au, au/day]

    # Compute point-mass Newtonian accelerations from ephemeris: all bodies except asteroid
    # accelerations of "everybody else" are needed when evaluating asteroid post-Newtonian acceleration
    # Number of bodies that contibute to the asteroid's acceleration
    Nm1 = numberofbodies(ss16asteph_x)
    # Initialize a TaylorInterpolant for the point-mass Newtonian accelerations
    acc_eph = TaylorInterpolant(ss16asteph.t0, ss16asteph.t, Matrix{eltype(ss16asteph.x)}(undef, length(ss16asteph.t)-1, 3Nm1))
    # Initialize a TaylorInterpolant for the newtonian N body potential
    newtonianNb_Potential = TaylorInterpolant(ss16asteph.t0, ss16asteph.t, Matrix{eltype(ss16asteph.x)}(undef, length(ss16asteph.t)-1, Nm1))
    # Fill TaylorInterpolant.x with zero polynomials
    fill!(acc_eph.x, zero(ss16asteph.x[1]))
    fill!(newtonianNb_Potential.x, zero(ss16asteph.x[1]))

    # Iterator over all bodies except asteroid
    _1_to_Nm1 = Base.OneTo(Nm1)
    for j in 1:Nm1
        for i in 1:Nm1
            if i == j
                #
            else
                # Difference between two positions (\mathbf{r}_i - \mathbf{r}_j)
                X_ij = ss16asteph.x[:, 3i-2] .- ss16asteph.x[:, 3j-2]  # X-axis component
                Y_ij = ss16asteph.x[:, 3i-1] .- ss16asteph.x[:, 3j-1]  # Y-axis component
                Z_ij = ss16asteph.x[:, 3i  ] .- ss16asteph.x[:, 3j  ]  # Z-axis component
                # Distance between two bodies squared ||\mathbf{r}_i - \mathbf{r}_j||^2
                r_p2_ij = ( (X_ij.^2) .+ (Y_ij.^2) ) .+ (Z_ij.^2)
                # Distance between two bodies ||\mathbf{r}_i - \mathbf{r}_j||
                r_ij = sqrt.(r_p2_ij)
                # Newtonian potential
                newtonianNb_Potential.x[:, j] .+= (μ[i]./r_ij)
            end
        end

        # Fill acelerations by differentiating velocities
        acc_eph.x[:, 3j-2] .= ordpres_differentiate.(ss16asteph.x[:, 3(Nm1+j)-2])  # X-axis component
        acc_eph.x[:, 3j-1] .= ordpres_differentiate.(ss16asteph.x[:, 3(Nm1+j)-1])  # Y-axis component
        acc_eph.x[:, 3j  ] .= ordpres_differentiate.(ss16asteph.x[:, 3(Nm1+j)  ])  # Z-axis component
    end

    return ss16asteph, acc_eph, newtonianNb_Potential
end

@doc raw"""
    save2jldandcheck(objname, sol)

Save `sol` in a file `objname_jt.jld2`.

See also [`__save2jldandcheck`](@ref).
"""
function save2jldandcheck(objname, sol)
    # Name of the file
    outfilename = string(objname, "_jt.jld2")
    # Save sol in outfilename
    return __save2jldandcheck(outfilename, sol)
end

@doc raw"""
    __save2jldandcheck(outfilename, sol)

Savs `sol` in `outfilename`.
"""
function __save2jldandcheck(outfilename, sol)
    println("Saving solution to file: $outfilename")
    # Open file
    JLD2.jldopen(outfilename, "w") do file
        # Loop over solution variables
        for ind in eachindex(sol)
            # Name of the variable
            varname = string(ind)
            println("Saving variable: ", varname)
            # Write the varaible
            write(file, varname, sol[ind])
        end
    end
    # Check that saved solution is equal to the original
    println("Checking that all variables were saved correctly...")
    # Loop over solution variables
    for ind in eachindex(sol)
        # Name of the variable
        varname = string(ind)
        # Read varname from files and assign recovered variable to recovered_sol_i
        recovered_sol_i = JLD2.load(outfilename, varname)
        # Check that varname was recovered succesfully
        @show recovered_sol_i == sol[ind]
    end
    println("Saved solution")
    return outfilename
end

@doc raw"""
    taylor_minimum(pol::Taylor1{T}, x0::T; niters::Int=10) where {T<:Real}

Return the minimum of the Taylor polynomial `pol` computed via Newton's method. `x0` is the
initial guess and `niters` is the number of iterations.
"""
function taylor_minimum(pol::Taylor1{T}, x0::T; niters::Int=10) where {T<:Real}
    # First derivative
    dpol = PE.ordpres_differentiate(pol)
    # Second derivative
    dpol2 = PE.ordpres_differentiate(dpol)
    # Initial guess
    xnewton::T = x0
    #@show xnewton
    # Newton iteration
    for i in 1:niters
        # Newton update rule
        xnewton -= dpol(xnewton)/dpol2(xnewton)
        #@show xnewton, dpol(xnewton)
    end

    return xnewton
end

@doc raw"""
    taylor_roots(pol::Taylor1{T}, x0::T; niters::Int=10) where {T<:Real}

Return the root of the Taylor polynomial `pol` computed via Newton's method. `x0` is the
initial guess and `niters` is the number of iterations.
"""
function taylor_roots(pol::Taylor1{T}, x0::T; niters::Int=10) where {T<:Real}
    # First derivative
    dpol = PE.ordpres_differentiate(pol)
    # Initial guess
    xnewton::T = x0
    #@show xnewton
    # Newton iteration
    for i in 1:niters
        # Newton update rule
        xnewton -= pol(xnewton)/dpol(xnewton)
        #@show xnewton, pol(xnewton)
    end
    return xnewton
end

@doc raw"""
    scaled_variables(c::Vector{T} = fill(1e-6, 6), names::String = "δx"; order::Int = 5) where {T <: Real}

Equivalent to `TaylorSeries.set_variables` times a scaling given by `c`.
"""
function scaled_variables(names::String = "δx", c::Vector{T} = fill(1e-6, 6); order::Int = 5) where {T <: Real}
    # Set TaylorN variables
    dq = set_variables(T, names, order = order, numvars = length(c))
    # Scale jet transport perturbation
    for i in eachindex(dq)
        dq[i][1][i] = c[i]
    end
    return dq
end

@doc raw"""
    propagate_params(jd0::T, ss16asteph_et::TaylorInterpolant, q0::Vector{U}; μ_ast::Vector = μ_ast343_DE430[1:end],
                     order::Int = order, abstol::T = abstol) where {T <: Real, U <: Number}

Return the parameters needed for `propagate`, `propagate_lyap` and `propagate_root`.

# Arguments

- `jd0::T`: initial Julian date.
- `ss16asteph_et::TaylorInterpolant`: solar system ephemeris.
- `q0::Vector{U}`: vector of initial conditions.
- `μ_ast::Vector`: vector of gravitational parameters.
- `order::Int`: order of the Taylor expansions to be used in the integration.
- `abstol::T`: absolute tolerance.
"""
function propagate_params(jd0::T, ss16asteph_et::TaylorInterpolant, q0::Vector{U}; μ_ast::Vector = μ_ast343_DE430[1:end],
                          order::Int = order, abstol::T = abstol) where {T <: Real, U <: Number}

    # Number of massive bodies
    Nm1 = numberofbodies(ss16asteph_et.x)

    # Number of bodies, including NEA
    N = Nm1 + 1

    # Vector of G*m values
    μ = convert(Vector{T}, vcat( μ_DE430[1:11], μ_ast[1:Nm1-11], zero(T) ) )

    # Check: number of SS bodies (N) in ephemeris must be equal to length of GM vector (μ)
    @assert N == length(μ) "Total number of bodies in ephemeris must be equal to length of GM vector μ"

    # Process ephemeris (convert from seconds to days)
    # Compute Newtonian accelerations and potentials (used in post-Newtonian accelerations)
    ss16asteph_auday, acc_eph, newtonianNb_Potential = loadeph(ss16asteph_et, μ)

    # Interaction matrix with flattened bodies
    UJ_interaction = fill(false, N)

    # Turn on Earth interaction
    UJ_interaction[ea] = true

    # Initial time
    _t0 = zero(T)

    # Time as Taylor variable
    t = _t0 + Taylor1( T, order )

    # Days since J2000
    dsj2k = t + (jd0 - JD_J2000)

    # Initial conditions
    _q0 = one(T) * q0

    # Auxiliary variable (to know the type of evaleph)
    aux_q0 = Taylor1(eltype(_q0), order)

    # Ephemeris at dsj2k
    ss16asteph_t = evaleph(ss16asteph_auday, dsj2k, aux_q0)

    # Accelerations at dsj2k
    acceph_t = evaleph(acc_eph, dsj2k, aux_q0)

    # Newtonian potentials at dsj2k
    newtonianNb_Potential_t = evaleph(newtonianNb_Potential, dsj2k, aux_q0)

    # Ephemeris vector
    _eph = (ss16asteph_auday, acc_eph, newtonianNb_Potential)

    # Vector of parameters for neosinteg
    _params = (ss16asteph_t, acceph_t, newtonianNb_Potential_t, jd0, UJ_interaction, N, μ)

    return _q0, _t0, _eph, _params

end

@doc raw"""
    propagate(dynamics::D, maxsteps::Int, jd0::T, tspan::T, ss16asteph_et::TaylorInterpolant,
              q0::Vector{U}, ::$V_dense; μ_ast::Vector = μ_ast343_DE430[1:end], order::Int = order,
              abstol::T = abstol, parse_eqs::Bool = true) where {T <: Real, U <: Number, D}

Integrate the orbit of a NEO via the Taylor method.

# Arguments

- `dynamics::D`: dynamical model function.
- `maxsteps::Int`: maximum number of steps for the integration.
- `jd0::T`: initial Julian date.
- `tspan::T`: time span of the integration [in Julian days].
- `ss16asteph_et::TaylorInterpolant`: solar system ephemeris.
- `q0::Vector{U}`: vector of initial conditions.
- `Val(true/false)`: whether to output the Taylor polynomials generated at each time step (`true`) or not.
- `μ_ast::Vector`: vector of gravitational parameters.
- `order::Int=order`: order of the Taylor expansions to be used in the integration.
- `abstol::T`: absolute tolerance.
- `parse_eqs::Bool`: whether to use the specialized method of `jetcoeffs` (`true`) or not.
""" propagate

@doc raw"""
    propagate_root(dynamics::D, maxsteps::Int, jd0::T, tspan::T, ss16asteph_et::TaylorInterpolant,
                   q0::Vector{U}, ::$V_dense; parse_eqs::Bool = true, eventorder::Int = 0, newtoniter::Int = 10,
                   nrabstol::T = eps(T), μ_ast::Vector = μ_ast343_DE430[1:end], order::Int = order,
                   abstol::T = abstol) where {T <: Real, U <: Number, D}

Integrate the orbit of a NEO via the Taylor method while finding the zeros of `rvelea`.

# Arguments

- `dynamics::D`: dynamical model function.
- `maxsteps::Int`: maximum number of steps for the integration.
- `jd0::T`: initial Julian date.
- `tspan::T`: time span of the integration [in Julian days].
- `ss16asteph_et::TaylorInterpolant`: solar system ephemeris.
- `q0::Vector{U}`: vector of initial conditions.
- `Val(true/false)`: whether to output the Taylor polynomials generated at each time step (`true`) or not.
- `parse_eqs::Bool`: whether to use the specialized method of `jetcoeffs` (`true`) or not.
- `eventorder::Int`: order of the derivative of `rvelea` whose roots are computed.
- `newtoniter::Int`: maximum Newton-Raphson iterations per detected root.
- `nrabstol::T`: allowed tolerance for the Newton-Raphson process.
- `μ_ast::Vector`: vector of gravitational parameters.
- `order::Int`: order of the Taylor expansions to be used in the integration.
- `abstol::T`: absolute tolerance.
""" propagate_root

for V_dense in V_true_false

    @eval begin

        function propagate(dynamics::D, maxsteps::Int, jd0::T, tspan::T, sseph::TaylorInterpolant, q0::Vector{U},
                           ::$V_dense; μ_ast::Vector = μ_ast343_DE430[1:end], order::Int = order, abstol::T = abstol,
                           parse_eqs::Bool = true) where {T <: Real, U <: Number, D}

            # Parameters for apophisinteg
            _q0, _t0, _eph, _params = propagate_params(jd0, sseph, q0; μ_ast = μ_ast, order = order, abstol = abstol)

            # Final time of integration (days)
            _tmax = _t0 + tspan*yr

            # Propagate orbit

            # Dense output (save Taylor polynomials in each step)
            @time sol = neosinteg(dynamics, _q0, _t0, _tmax, order, abstol, $V_dense(), _eph, _params;
                                  maxsteps = maxsteps, parse_eqs = parse_eqs)

            return sol

        end

        function propagate(dynamics::D, maxsteps::Int, jd0::T1, tspan::T2, sseph::TaylorInterpolant, q0::Vector{U},
                           ::$V_dense; μ_ast::Vector = μ_ast343_DE430[1:end], order::Int = order, abstol::T3 = abstol,
                           parse_eqs::Bool = true) where {T1, T2, T3 <: Real, U <: Number, D}

            _jd0, _tspan, _abstol = promote(jd0, tspan, abstol)

            return propagate(dynamics, maxsteps, _jd0, _tspan, sseph, q0, $V_dense(); μ_ast = μ_ast, order = order,
                             abstol = _abstol, parse_eqs = parse_eqs)
        end

        function propagate(dynamics::D, maxsteps::Int, jd0::T, nyears_bwd::T, nyears_fwd::T, sseph::TaylorInterpolant,
                           q0::Vector{U}, ::$V_dense; μ_ast::Vector = μ_ast343_DE430[1:end], order::Int = order, abstol::T = abstol,
                           parse_eqs::Bool = true) where {T <: Real, U <: Number, D}

            # Backward integration
            bwd = propagate(dynamics, maxsteps, jd0, nyears_bwd, sseph, q0, $V_dense(); μ_ast = μ_ast, order = order,
                            abstol = abstol, parse_eqs = parse_eqs)
            # Forward integration
            fwd = propagate(dynamics, maxsteps, jd0, nyears_fwd, sseph, q0, $V_dense(); μ_ast = μ_ast, order = order,
                            abstol = abstol, parse_eqs = parse_eqs)

            return bwd, fwd

        end

        function propagate_root(dynamics::D, maxsteps::Int, jd0::T, tspan::T, sseph::TaylorInterpolant, q0::Vector{U},
                                ::$V_dense; parse_eqs::Bool = true, eventorder::Int = 0, newtoniter::Int = 10, nrabstol::T = eps(T),
                                μ_ast::Vector = μ_ast343_DE430[1:end], order::Int = order, abstol::T = abstol) where {T <: Real, U <: Number, D}

            # Parameters for apophisinteg
            _q0, _t0, _eph, _params = propagate_params(jd0, sseph, q0; μ_ast = μ_ast, order = order, abstol = abstol)

            # Final time of integration (days)
            _tmax = _t0 + tspan*yr

            # Propagate orbit
            @time sol = neosinteg(dynamics, rvelea, _q0, _t0, _tmax, order, abstol, $V_dense(), _eph, _params;
                                  maxsteps = maxsteps, parse_eqs = parse_eqs,  eventorder = eventorder, newtoniter = newtoniter,
                                  nrabstol = nrabstol)

            return sol

        end

        function propagate_root(dynamics::D, maxsteps::Int, jd0::T1, tspan::T2, sseph::TaylorInterpolant, q0::Vector{U},
                                ::$V_dense; parse_eqs::Bool = true, eventorder::Int = 0, newtoniter::Int = 10, nrabstol::T3 = eps(T),
                                μ_ast::Vector = μ_ast343_DE430[1:end], order::Int = order, abstol::T4 = abstol) where {T1, T2, T3, T4 <: Real, U <: Number, D}

            _jd0, _tspan, _nrabstol, _abstol = promote(jd0, tspan, nrabstol, abstol)

            return propagate_root(dynamics, maxsteps, _jd0, _tspan, sseph, q0, $V_dense(); parse_eqs = parse_eqs,
                                  eventorder = eventorder, newtoniter = newtoniter, nrabstol = _nrabstol, μ_ast = μ_ast,
                                  order = order, abstol = _abstol)
        end

        function propagate_root(dynamics::D, maxsteps::Int, jd0::T, nyears_bwd::T, nyears_fwd::T, sseph::TaylorInterpolant,
                                q0::Vector{U}, ::$V_dense; parse_eqs::Bool = true, eventorder::Int = 0, newtoniter::Int = 10,
                                nrabstol::T = eps(T), μ_ast::Vector = μ_ast343_DE430[1:end], order::Int = order, abstol::T = abstol) where {T <: Real, U <: Number, D}

            # Backward integration
            bwd, tvS_bwd, xvS_bwd, gvS_bwd = propagate_root(dynamics, maxsteps, jd0, nyears_bwd, sseph, q0, $V_dense(); parse_eqs = parse_eqs,
                                                            eventorder = eventorder, newtoniter = newtoniter, nrabstol = nrabstol,
                                                            μ_ast = μ_ast, order = order, abstol = abstol)
            # Forward integration
            fwd, tvS_fwd, xvS_fwd, gvS_fwd = propagate_root(dynamics, maxsteps, jd0, nyears_fwd, sseph, q0, $V_dense(); parse_eqs = parse_eqs,
                                                            eventorder = eventorder, newtoniter = newtoniter, nrabstol = nrabstol,
                                                            μ_ast = μ_ast, order = order, abstol = abstol)

            return bwd, tvS_bwd, xvS_bwd, gvS_bwd, fwd, tvS_fwd, xvS_fwd, gvS_fwd

        end

    end
end

@doc raw"""
    propagate_lyap(dynamics::D, maxsteps::Int, jd0::T, tspan::T, ss16asteph_et::TaylorInterpolant,
                   q0::Vector{U}; μ_ast::Vector = μ_ast343_DE430[1:end], order::Int = order,
                   abstol::T = abstol, parse_eqs::Bool = true) where {T <: Real, U <: Number}

Compute the Lyapunov spectrum of a NEO.

# Arguments

- `dynamics::D`: dynamical model function.
- `maxsteps::Int`: maximum number of steps for the integration.
- `jd0::T`: initial Julian date.
- `tspan::T`: time span of the integration [in Julian days].
- `ss16asteph_et::TaylorInterpolant`: solar system ephemeris.
- `q0::Vector{U}`: vector of initial conditions.
- `μ_ast::Vector`: vector of gravitational parameters.
- `order::Int=order`: order of the Taylor expansions to be used in the integration.
- `abstol::T`: absolute tolerance.
- `parse_eqs::Bool`: whether to use the specialized method of `jetcoeffs` (`true`) or not.
""" propagate_lyap
function propagate_lyap(dynamics::D, maxsteps::Int, jd0::T, tspan::T, ss16asteph_et::TaylorInterpolant,
                        q0::Vector{U}; μ_ast::Vector = μ_ast343_DE430[1:end], order::Int = order,
                        abstol::T = abstol, parse_eqs::Bool = true) where {T <: Real, U <: Number, D}

    # Parameters for apophisinteg
    _q0, _t0, _eph, _params = propagate_params(jd0, ss16asteph_et, q0; μ_ast = μ_ast, order = order, abstol = abstol)

    # Final time of integration (days)
    _tmax = _t0 + tspan*yr

    # Propagate orbit
    @time sol = lyap_neosinteg(dynamics, _q0, _t0, _tmax, order, abstol, _eph, _params;
                               maxsteps = maxsteps, parse_eqs = parse_eqs)

    return sol

end