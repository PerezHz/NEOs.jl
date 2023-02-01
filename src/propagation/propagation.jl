include("gauss_method.jl")
include("asteroid_dynamical_models.jl")
include("jetcoeffs.jl")
include("integration_methods.jl")

@doc raw"""
    rvelea(dx, x, params, t)

Return `true` and the numerator of the asteroid's radial velocity with respect to the Earth.

# Arguments 

- `dx`: asteroid's velocities. 
- `x`: asteroid's degrees of freedom. 
- `params`: parameters (ephemeris + accelerations + newtonian N body potential + julian date of start time + matrix of extended body interactions + number of bodies + mass parameters). 
- `t`: time. 
"""
function rvelea(dx, x, params, t)
    jd0 = params[4]                                # Julian date of start time
    dsj2k = t+(jd0-JD_J2000)                       # Days since J2000.0 = 2.451545e6
    ss16asteph_t = evaleph(params[1], dsj2k, x[1]) # params[2](t)*one(q[1]) # ss16asteph(t)
    N = params[6]                                  # Number of bodies in the ephemeris
    xe = ss16asteph_t[nbodyind(N-1,ea)]            # Earth's ephemeris 
    return true, (x[1]-xe[1])*(x[4]-xe[4]) + (x[2]-xe[2])*(x[5]-xe[5]) + (x[3]-xe[3])*(x[6]-xe[6])
end

@doc raw"""
    loadeph(ss16asteph_::TaylorInterpolant, μ::Vector)

Return the ephemeris in `ss16asteph_` with times converted from seconds to days, 
the point-mass newtonian accelerations and the newtonian N body potential. 

# Arguments 

- `ss16asteph_`: solar system ephemeris. 
- `μ`: vector of mass parameters. 
"""
function loadeph(ss16asteph_::TaylorInterpolant, μ::Vector)
    # Read Solar System ephemeris (Sun+8 planets+Moon+Pluto+16 main belt asteroids)
    ss16asteph_t0 = (ss16asteph_.t0 ./ daysec)    # Start time (days)
    ss16asteph_t = (ss16asteph_.t ./ daysec)      # Vector of times (days)   
    ephord = ss16asteph_.x[1].order               # Order of the Taylor polynomials     
    ss16asteph_x = map(x->x(Taylor1(ephord)*daysec), ss16asteph_.x)            # Vector of Taylor polynomials 
    ss16asteph = TaylorInterpolant(ss16asteph_t0, ss16asteph_t, ss16asteph_x)  # TaylorInterpolant ephemeris
    # Compute point-mass Newtonian accelerations from ephemeris: all bodies except asteroid
    # accelerations of "everybody else" are needed when evaluating asteroid post-Newtonian acceleration
    # Number of bodies that contibute to the asteroid's acceleration
    Nm1 = (size(ss16asteph_x)[2]-13) ÷ 6 
    # Initialize a TaylorInterpolant for the point-mass Newtonian accelerations
    acc_eph = TaylorInterpolant(ss16asteph.t0, ss16asteph.t, Matrix{eltype(ss16asteph.x)}(undef, length(ss16asteph.t)-1, 3Nm1))
    # Initialize a TaylorInterpolant for the newtonian N body potential 
    newtonianNb_Potential = TaylorInterpolant(ss16asteph.t0, ss16asteph.t, Matrix{eltype(ss16asteph.x)}(undef, length(ss16asteph.t)-1, Nm1))
    # Fill TaylorInterpolant.x with zero polynomials
    fill!(acc_eph.x, zero(ss16asteph.x[1]))
    fill!(newtonianNb_Potential.x, zero(ss16asteph.x[1]))
    # Iterator over all bodies except asteroid
    _1_to_Nm1 = Base.OneTo(Nm1) 
    for j in _1_to_Nm1
        for i in _1_to_Nm1
            # i == j && continue
            if i == j
                # 
            else
                # Difference between two positions (\mathbf{r}_i - \mathbf{r}_j)
                X_ij = ss16asteph.x[:,3i-2] .- ss16asteph.x[:,3j-2]  # X-axis component
                Y_ij = ss16asteph.x[:,3i-1] .- ss16asteph.x[:,3j-1]  # Y-axis component
                Z_ij = ss16asteph.x[:,3i  ] .- ss16asteph.x[:,3j  ]  # Z-axis component
                # Distance between two bodies squared ||\mathbf{r}_i - \mathbf{r}_j||^2
                r_p2_ij = ( (X_ij.^2) .+ (Y_ij.^2) ) .+ (Z_ij.^2)
                # Distance between two bodies ||\mathbf{r}_i - \mathbf{r}_j||
                r_ij = sqrt.(r_p2_ij)
                # Newtonian potential
                newtonianNb_Potential.x[:,j] .+= (μ[i]./r_ij)
            end
        end
        # Fill acelerations by differentiating velocities 
        acc_eph.x[:,3j-2] .= PlanetaryEphemeris.ordpres_differentiate.(ss16asteph.x[:,3(Nm1+j)-2])  # X-axis component
        acc_eph.x[:,3j-1] .= PlanetaryEphemeris.ordpres_differentiate.(ss16asteph.x[:,3(Nm1+j)-1])  # Y-axis component
        acc_eph.x[:,3j  ] .= PlanetaryEphemeris.ordpres_differentiate.(ss16asteph.x[:,3(Nm1+j)  ])  # Z-axis component
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
    dpol = PlanetaryEphemeris.ordpres_differentiate(pol)
    # Second derivative
    dpol2 = PlanetaryEphemeris.ordpres_differentiate(dpol)
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
    dpol = PlanetaryEphemeris.ordpres_differentiate(pol)
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
    scaling(a::Taylor1{Taylor1{T}}, c::T) where {T<:Real}

Scale `a` by a factor `c`. 
"""
function scaling(a::Taylor1{Taylor1{T}}, c::T) where {T<:Real}
    x = c*Taylor1( Taylor1(a.order).coeffs*one(a[0]) )
    return a(x)
end

@doc raw"""
    propagate_dense(objname::String, dynamics::Function, maxsteps::Int, jd0::T, tspan::T, ss16asteph_et::TaylorInterpolant,
                    q0::Vector{U}, Val(true/false); output::Bool = true, μ_ast::Vector = μ_ast343_DE430[1:end], 
                    order::Int = order, abstol::T = abstol, parse_eqs::Bool = true) where {T <: Real, U <: Number}

Integrate the orbit of a NEO via the Taylor method. 

# Arguments 

- `objname::String`: name of the object. 
- `dynamics::Function`: dynamical model function.
- `maxsteps::Int`: maximum number of steps for the integration.
- `jd0::T`: initial Julian date.
- `tspan::T`: time span of the integration [in Julian days]. 
- `ss16asteph_et::TaylorInterpolant`: solar system ephemeris.
- `q0::Vector{T}`: vector of initial conditions.
- `Val(true/false)`: Whether to use quadruple precision or not.
- `output::Bool`: whether to write the output to a file (`true`) or not.
- `μ_ast::Vector`: vector of gravitational parameters. 
- `order::Int=order`: order of the Taylor expansions to be used in the integration. 
- `abstol::T`: absolute tolerance.
- `parse_eqs::Bool`: whether to use the specialized method of `jetcoeffs` (`true`) or not. 
""" propagate_dense

@doc raw"""
    propagate_lyap(objname::String, dynamics::Function, maxsteps::Int, jd0::T, tspan::T, ss16asteph_et::TaylorInterpolant,
                   q0::Vector{U}, Val(true/false); output::Bool = true, μ_ast::Vector = μ_ast343_DE430[1:end],
                   order::Int = order, abstol::T = abstol, parse_eqs::Bool = true) where {T <: Real, U <: Number}

Compute the Lyapunov spectrum of a NEO. 

# Arguments 

- `objname::String`: name of the object. 
- `dynamics::Function`: dynamical model function.
- `maxsteps::Int`: maximum number of steps for the integration.
- `jd0::T`: initial Julian date.
- `tspan::T`: time span of the integration [in Julian days]. 
- `ss16asteph_et::TaylorInterpolant`: solar system ephemeris.
- `q0::Vector{T}`: vector of initial conditions.
- `Val(true/false)`: Whether to use quadruple precision or not.
- `output::Bool`: whether to write the output to a file (`true`) or not.
- `μ_ast::Vector`: vector of gravitational parameters. 
- `order::Int=order`: order of the Taylor expansions to be used in the integration. 
- `abstol::T`: absolute tolerance.
- `parse_eqs::Bool`: whether to use the specialized method of `jetcoeffs` (`true`) or not. 
""" propagate_lyap

@doc raw"""
    propagate_root(objname::String, dynamics::Function, maxsteps::Int, jd0::T, tspan::T, ss16asteph_et::TaylorInterpolant,
                   q0::Vector{U}, Val(true/false); output::Bool = true, eventorder::Int = 0, newtoniter::Int = 10, 
                   nrabstol::T = eps(T), μ_ast::Vector = μ_ast343_DE430[1:end], order::Int = order, abstol::T = abstol, 
                   parse_eqs::Bool = true) where {T <: Real, U <: Number}

Integrate the orbit of a NEO via the Taylor method while finding the zeros of `rvelea`.  

# Arguments 

- `objname::String`: name of the object. 
- `dynamics::Function`: dynamical model function.
- `maxsteps::Int`: maximum number of steps for the integration.
- `jd0::T`: initial Julian date.
- `tspan::T`: time span of the integration [in Julian days]. 
- `ss16asteph_et::TaylorInterpolant`: solar system ephemeris.
- `q0::Vector{T}`: vector of initial conditions.
- `Val(true/false)`: Whether to use quadruple precision or not.
- `output::Bool`: whether to write the output to a file (`true`) or not.
- `eventorder::Int`: order of the derivative of `rvelea` whose roots are computed.
- `newtoniter::Int`: maximum Newton-Raphson iterations per detected root.
- `nrabstol::T`: allowed tolerance for the Newton-Raphson process.
- `μ_ast::Vector`: vector of gravitational parameters. 
- `order::Int=order`: order of the Taylor expansions to be used in the integration. 
- `abstol::T`: absolute tolerance.
- `parse_eqs::Bool`: whether to use the specialized method of `jetcoeffs` (`true`) or not. 
""" propagate_root 

for V_quadmath in V_true_false
    @eval begin
        function propagate_params(jd0::T, ss16asteph_et::TaylorInterpolant, q0::Vector{U}, ::$V_quadmath; 
                                  μ_ast::Vector = μ_ast343_DE430[1:end], abstol::T=abstol) where {T <: Real, U <: Number}
 
            # Number of massive bodies
            Nm1 = (size(ss16asteph_et.x)[2]-13) ÷ 6 
        
            # Number of bodies, including NEA
            N = Nm1 + 1 
        
            # Vector of G*m values
            μ = vcat(μ_DE430[1:11], μ_ast[1:Nm1-11], zero(μ_DE430[1]))

            # Check: number of SS bodies (N) in ephemeris must be equal to length of GM vector (μ)
            @assert N == length(μ) "Total number of bodies in ephemeris must be equal to length of GM vector μ"
        
            # Process ephemeris (switch from km, km/s units to au,au/day)
            # Compute Newtonian accelerations and potentials (used in post-Newtonian accelerations)
            ss16asteph_auday, acc_eph, newtonianNb_Potential = loadeph(ss16asteph_et, μ)
        
            # Interaction matrix with flattened bodies
            UJ_interaction = fill(false, N)
        
            # Turn on Earth interaction 
            UJ_interaction[ea] = true
        
            # Vector of parameters for apophisinteg 
            params = (ss16asteph_auday, acc_eph, newtonianNb_Potential, jd0, UJ_interaction, N, μ)

            # Use quadruple precision
            if $V_quadmath == Val{true}
                _q0 = one(Float128)*q0
                _t0 = zero(Float128)
                _abstol = Float128(abstol)
                _ss16asteph = TaylorInterpolant(
                    Float128(ss16asteph_auday.t0), 
                    Float128.(ss16asteph_auday.t), 
                    map(x->Taylor1(Float128.(x.coeffs)), ss16asteph_auday.x)
                )
                _acc_eph = TaylorInterpolant(Float128(acc_eph.t0), Float128.(acc_eph.t), map(x->Taylor1(Float128.(x.coeffs)), acc_eph.x))
                _newtonianNb_Potential = TaylorInterpolant(Float128(newtonianNb_Potential.t0), Float128.(newtonianNb_Potential.t), map(x->Taylor1(Float128.(x.coeffs)), newtonianNb_Potential.x))
                _params = (_ss16asteph, _acc_eph, _newtonianNb_Potential, Float128(jd0), UJ_interaction, N, μ)

                return _q0, _t0, _abstol, _params
            # Use double precision
            else
                
                t0 = zero(Float64)
                
                return q0, t0, abstol, params 
            end
        
        end 

        function propagate_dense(objname::String, dynamics::Function, maxsteps::Int, jd0::T, tspan::T, ss16asteph_et::TaylorInterpolant,
                                 q0::Vector{U}, ::$V_quadmath; output::Bool = true, μ_ast::Vector = μ_ast343_DE430[1:end], 
                                 order::Int = order, abstol::T = abstol, parse_eqs::Bool = true) where {T <: Real, U <: Number}

            # Parameters for apophisinteg 
            if $V_quadmath == Val{true}
                _q0, _t0, _abstol, _params = propagate_params(jd0, ss16asteph_et, q0, Val(true); μ_ast = μ_ast, abstol = abstol)
            else 
                _q0, _t0, _abstol, _params = propagate_params(jd0, ss16asteph_et, q0, Val(false); μ_ast = μ_ast, abstol = abstol)
            end

            println("Initial time of integration: ", julian2datetime(jd0))
            # Final time of integration (days)
            _tmax = _t0 + tspan*yr 
            println("Final time of integration: ", julian2datetime(jd0 + _tmax))

            # Propagate orbit

            # Dense output (save Taylor polynomials in each step)
            @time interp = neosinteg(dynamics, _q0, _t0, _tmax, order, _abstol, Val(true), _params; maxsteps = maxsteps, parse_eqs = parse_eqs)
            
            # Days since J2000 until initial integration time
            asteph_t0 = T(jd0 - JD_J2000)
            # Vector of times 
            asteph_t = T.(interp.t[:])
            # Matrix of dense polynomials
            asteph_x = convert(Array{Taylor1{U}}, interp.x[:,:])
            # TaylorInterpolant 
            asteph = TaylorInterpolant(asteph_t0, asteph_t, asteph_x)
            # Solution 
            sol = (asteph=asteph,)

            # Write solution and predicted values of observations (if requested) to .jld files
            if output
                # Name of the file 
                _ = save2jldandcheck(objname, sol)
            end

            return sol 

        end

        function propagate_lyap(objname::String, dynamics::Function, maxsteps::Int, jd0::T, tspan::T, ss16asteph_et::TaylorInterpolant,
                                q0::Vector{U}, ::$V_quadmath; output::Bool = true, μ_ast::Vector = μ_ast343_DE430[1:end],
                                order::Int = order, abstol::T = abstol, parse_eqs::Bool = true) where {T <: Real, U <: Number}

            # Parameters for apophisinteg 
            if $V_quadmath == Val{true}
                _q0, _t0, _abstol, _params = propagate_params(jd0, ss16asteph_et, q0, Val(true); μ_ast = μ_ast, abstol = abstol)
            else 
                _q0, _t0, _abstol, _params = propagate_params(jd0, ss16asteph_et, q0, Val(false); μ_ast = μ_ast, abstol = abstol)
            end

            println("Initial time of integration: ", julian2datetime(jd0))
            # Final time of integration (days)
            _tmax = _t0 + tspan*yr 
            println("Final time of integration: ", julian2datetime(jd0 + _tmax))

            # Propagate orbit
            @time sol_objs = lyap_neosinteg(dynamics, _q0, _t0, _tmax, order, _abstol, Val(false), _params; maxsteps = maxsteps, parse_eqs = parse_eqs)

            # Solution 
            sol = (
                tv = convert(Array{U}, sol_objs[1][:]),
                xv = convert(Array{U}, sol_objs[2][:,:]),
                λv = convert(Array{U}, sol_objs[3][:,:])
            )

            # Write solution and predicted values of observations (if requested) to .jld files
            if output
                # Name of the file 
                _ = save2jldandcheck(objname, sol)
            end

            return sol

        end

        function propagate_root(objname::String, dynamics::Function, maxsteps::Int, jd0::T, tspan::T, ss16asteph_et::TaylorInterpolant,
                                q0::Vector{U}, ::$V_quadmath; output::Bool = true, parse_eqs::Bool = true, eventorder::Int = 0, 
                                newtoniter::Int = 10, nrabstol::T = eps(T), μ_ast::Vector = μ_ast343_DE430[1:end], order::Int = order, 
                                abstol::T = abstol) where {T <: Real, U <: Number}

            # Parameters for apophisinteg 
            if $V_quadmath == Val{true}
                _q0, _t0, _abstol, _params = propagate_params(jd0, ss16asteph_et, q0, Val(true); μ_ast = μ_ast, abstol = abstol)
            else 
                _q0, _t0, _abstol, _params = propagate_params(jd0, ss16asteph_et, q0, Val(false); μ_ast = μ_ast, abstol = abstol)
            end

            println("Initial time of integration: ", julian2datetime(jd0))
            # Final time of integration (days)
            _tmax = _t0 + tspan*yr 
            println("Final time of integration: ", julian2datetime(jd0 + _tmax))

            # Propagate orbit
            @time sol_objs = neosinteg(dynamics, rvelea, _q0, _t0, _tmax, order, _abstol, Val(true), _params; maxsteps = maxsteps, 
                                       parse_eqs = parse_eqs,  eventorder = eventorder, newtoniter = newtoniter, nrabstol = nrabstol)

            # Days since J2000 until initial integration time
            asteph_t0 = (_params[4] - JD_J2000) 
            # Vector of times 
            asteph_t = sol_objs[1].t[:]
            # Matrix of polynomials
            asteph_x = sol_objs[1].x[:,:]
            # TaylorInterpolant
            asteph = TaylorInterpolant(asteph_t0, asteph_t, asteph_x)
            # Solution 
            sol = (                    
                asteph = asteph,                              # TaylorInterpolant
                tv = asteph_t0 .+ asteph_t,                   # Vector of times 
                xv = asteph_x(),                              # Polynomials evaluated at t = 0
                tvS1 = convert(Array{U}, sol_objs[2][:]),     # Coordinates when rvelea has a zero 
                xvS1 = convert(Array{U}, sol_objs[3][:,:]),   
                gvS1 = convert(Array{U}, sol_objs[4][:])      # Zeros of rvelea
            )

            # Write solution and predicted values of observations (if requested) to .jld files
            if output
                # Name of the file 
                _ = save2jldandcheck(objname, sol)
            end

            return sol

        end

    end
end 