include("gauss_method.jl")
include("initial_conditions.jl")
include("asteroid_dynamical_models.jl")
include("integration_methods.jl")

@doc raw"""
    rvelea(dx, x, params, t)

Returns `true` and the numerator of the asteroid's radial velocity with respect to the Earth.

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

Returns the ephemeris in `ss16asteph_` with times converted from seconds to days, 
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

Saves `sol` in a file `objname_jt.jld`. 

See also [`__save2jldandcheck`](@ref). 
"""
function save2jldandcheck(objname, sol)
    # Name of the file 
    outfilename = string(objname, "_jt.jld")
    # Save sol in outfilename
    return __save2jldandcheck(outfilename, sol)
end

@doc raw"""
    __save2jldandcheck(outfilename, sol)

Saves `sol` in `outfilename`. 
"""
function __save2jldandcheck(outfilename, sol)
    println("Saving solution to file: $outfilename")
    # Open file 
    jldopen(outfilename, "w") do file
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
        recovered_sol_i = JLD.load(outfilename, varname)
        # Check that varname was recovered succesfully
        @show recovered_sol_i == sol[ind]
    end
    println("Saved solution")
    return outfilename
end

@doc raw"""
    taylor_minimum(pol::Taylor1{T}, x0::T; niters::Int=10) where {T<:Real}

Returns the minimum of the Taylor polynomial `pol` computed via Newton's method. `x0` is the
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

Returns the root of the Taylor polynomial `pol` computed via Newton's method. `x0` is the 
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

Scales `a` by a factor `c`. 
"""
function scaling(a::Taylor1{Taylor1{T}}, c::T) where {T<:Real}
    x = c*Taylor1( Taylor1(a.order).coeffs*one(a[0]) )
    return a(x)
end

@doc raw"""
    propagate(objname::String, dynamics::Function, maxsteps::Int, jd0::T, tspan::T, 
              ephfile::String; kwargs...) where {T<:Real}

Integrates the orbit of an asteroid via the Taylor method. 

# Arguments 

- `objname::String`: name of the object. 
- `dynamics::Function`: dynamical model function.
- `maxsteps::Int`: maximum number of steps for the integration.
- `jd0::T`: initial Julian date.
- `tspan::T`: time span of the integration (in Julian days). 
- `ephfile::String`: file with solar system ephemeris. 
"""
function propagate(objname::String, dynamics::Function, maxsteps::Int, jd0::T,
        tspan::T, ephfile::String; kwargs...) where {T<:Real}
    # Load ephemeris
    ss16asteph_et = JLD.load(ephfile, "ss16ast_eph")
    # Propagate 
    propagate(objname::String, dynamics::Function, maxsteps::Int, jd0::T, tspan::T, ss16asteph_et; kwargs...)
end

const V_true = :(Val{true})
const V_false = :(Val{false})
const V_true_false = (V_true, V_false)

for V_quadmath in V_true_false
    @eval begin
        function propagate_params(ss16asteph_et::TaylorInterpolant, ::$V_quadmath; q0::Vector=initialcond(), 
                                  μ_ast::Vector = μ_ast343_DE430[1:end], abstol::T=abstol) where {T <: Real}
 
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

        function propagate_dense_sol(jd0::T, q0::Vector, interp, ::$V_quadmath) where {T <: Real}

            if $V_quadmath == Val{true}
                # Days since J2000 until initial integration time
                asteph_t0 = Float64(jd0 - JD_J2000)
                # Vector of times 
                asteph_t = Float64.(interp.t[:])
                # Matrix of dense polynomials
                asteph_x = convert(Array{Taylor1{eltype(q0)}}, interp.x[:,:])
                
            # Use double precision
            else
                # Days since J2000 until initial integration time
                asteph_t0 = (jd0 - JD_J2000) 
                # Vector of times
                asteph_t = interp.t[:]
                # Matrix of dense polynomials
                asteph_x = interp.x[:,:]
                
            end
            # TaylorInterpolant 
            asteph = TaylorInterpolant(asteph_t0, asteph_t, asteph_x)
            # Solution 
            sol = (asteph=asteph,)

            return sol 

        end 

    end
end 

for (V_quadmath, V_dense, V_lyap) in Iterators.product(V_true_false, V_true_false, V_true_false)
    @eval begin
        function propagate(objname::String, dynamics::Function, maxsteps::Int, jd0::T, tspan::T, ss16asteph_et::TaylorInterpolant,
                           ::$V_quadmath, ::$V_dense, ::$V_lyap; output::Bool = true, newtoniter::Int=10, q0::Vector=initialcond(), 
                           radarobsfile::String="", opticalobsfile::String="", debias_table::String="2018", μ_ast::Vector=μ_ast343_DE430[1:end],
                           order::Int=order, abstol::T=abstol, tord::Int=10, niter::Int=5) where {T <: Real}
        
            # Parameters for apophisinteg 
            if $V_quadmath == Val{true}
                _q0, _t0, _abstol, _params = propagate_params(ss16asteph_et, Val(true); q0 = q0, μ_ast = μ_ast, abstol = abstol)
            else 
                _q0, _t0, _abstol, _params = propagate_params(ss16asteph_et, Val(false); q0 = q0, μ_ast = μ_ast, abstol = abstol)
            end
        
            # Final time of integration
            _tmax = _t0 + tspan*yr 
            println("Final time of integration: ", _tmax)
            
            # Propagate orbit
        
            # Dense output (save Taylor polynomials in each step)
            if $V_dense == Val{true}

                # Propagation 
                @time interp = apophisinteg(dynamics, _q0, _t0, _tmax, order, _abstol, _params; maxsteps = maxsteps, dense = true)

                # Use quadruple precision 
                if $V_quadmath == Val{true}
                    sol = propagate_dense_sol(jd0, q0, interp, Val(true))
                # Use double precision
                else 
                    sol = propagate_dense_sol(jd0, q0, interp, Val(false))
                end 

            # Lyapunov spectrum
            elseif $V_lyap == Val{true}

                # Propagation 
                @time sol_objs = lyap_apophisinteg(dynamics, _q0, _t0, _tmax, order, _abstol, _params; maxsteps = maxsteps)
                # Solution 
                sol = (
                    tv=convert(Array{eltype(q0)}, sol_objs[1][:]),
                    xv=convert(Array{eltype(q0)}, sol_objs[2][:,:]),
                    λv=convert(Array{eltype(q0)}, sol_objs[3][:,:])
                )
            # Not dense or lyap
            else

                # Propagation 
                @time sol_objs = apophisinteg(dynamics, rvelea, _q0, _t0, _tmax, order, _abstol, _params; maxsteps = maxsteps, 
                                              newtoniter = newtoniter, dense = true)
                # Days since J2000 until initial integration time
                asteph_t0 = (_params[4]-JD_J2000) 
                # Vector of times 
                asteph_t = sol_objs[1].t[:]
                # Matrix of polynomials
                asteph_x = sol_objs[1].x[:,:]
                # TaylorInterpolant
                asteph = TaylorInterpolant(asteph_t0, asteph_t, asteph_x)
                # Solution 
                sol = (                    
                    asteph = asteph,                                       # TaylorInterpolant
                    tv = asteph_t0 .+ asteph_t,                            # Vector of times 
                    xv = asteph_x(),                                       # Polynomials evaluated at t = 0
                    tvS1 = convert(Array{eltype(q0)}, sol_objs[2][:]),     # 
                    xvS1 = convert(Array{eltype(q0)}, sol_objs[3][:,:]),   # 
                    gvS1 = convert(Array{eltype(q0)}, sol_objs[4][:])      # 
                )

            end
        
            # Write solution and predicted values of observations (if requested) to .jld files
            if output
                # Name of the file 
                outfilename = save2jldandcheck(objname, sol)
                # If requested by user, calculate computed (i.e., predicted) values of observations
                if $V_lyap == Val{false} # Don't compute observation ephemeris when "in Lyapunov spectra mode"
                    furnsh(
                        joinpath(artifact"naif0012", "naif0012.tls"),     # Load leapseconds kernel
                        joinpath(artifact"de430", "de430_1850-2150.bsp"), # at least one SPK file must be loaded to read .tls file
                    )
                end
            end
        
            return sol

        end
    end
end 

@doc raw"""
    propagate(objname::String, dynamics::Function, maxsteps::Int, jd0::T, tspan::T, 
              ss16asteph_et::TaylorInterpolant; output::Bool=true, newtoniter::Int=10,
              dense::Bool=false, q0::Vector=initialcond(), radarobsfile::String="",
              opticalobsfile::String="", quadmath::Bool=false, debias_table::String="2018", 
              μ_ast::Vector=μ_ast343_DE430[1:end], lyap::Bool=false, order::Int=order, 
              abstol::T=abstol, tord::Int=10, niter::Int=5) where {T<:Real}

Integrates the orbit of an asteroid via the Taylor method. 

# Arguments 
- `objname::String`: name of the object. 
- `dynamics::Function`: dynamical model function.
- `maxsteps::Int`: maximum number of steps for the integration.
- `jd0::T`: initial Julian date.
- `tspan::T`: time span of the integration (in Julian days). 
- `ss16asteph_et::TaylorInterpolant`: solar system ephemeris. 
- `output::Bool`: whether to write the output to a file (`true`) or not.
- `newtoniter::Int`: number of iterations for root-finding integration. 
- `dense::Bool`: whether to save the Taylor polynomials at each step (`true`) or not.
- `q0::Vector`: vector of initial conditions.
- `radarobsfile::String`: file with radar observations.
- `opticalobsfile::String`: file with optical observations. 
- `quadmath::Bool`: whether to use quadruple precision (`true`) or not. 
- `debias_table::String`: debias table for optical observations. 
- `μ_ast::Vector`: vector of mass parameters. 
- `lyap::Bool`: wheter to compute the Lyapunov spectrum (`true`) or not. 
- `order::Int=order`: order of the Taylor expansions to be used in the integration. 
- `abstol::T`: absolute tolerance.
- `tord::Int`: order of Taylor expansions for computing time delays and Doppler shifts. 
- `niter::Int`: number of iterations for computing radar and optical observations. 
""" propagate

@doc raw"""
    compute_radar_obs(outfilename::String, radarobsfile::String, asteph::TaylorInterpolant, 
                      ss16asteph::TaylorInterpolant; tc::Real=1.0, autodiff::Bool=true, 
                      tord::Int=10, niter::Int=5)

Computes time delays and Doppler shifts and saves the result to a file. 

# Arguments 

- `outfilename::String`: file where to save radar observations. 
- `radarobsfile::String`: file where to retrieve radar observations. 
- `asteph::TaylorInterpolant`: asteroid's ephemeris. 
- `ss16asteph::TaylorInterpolant`: solar system ephemeris. 
- `tc::Real`: time offset wrt echo reception time, to compute Doppler shifts by range differences (seconds).
- `autodiff::Bool`: wheter to use the automatic differentiation method of [`delay`](@ref) or not. 
- `tord::Int`: order of Taylor expansions. 
- `niter::Int`: number of light-time solution iterations. 
"""
function compute_radar_obs(outfilename::String, radarobsfile::String, asteph::TaylorInterpolant,
                           ss16asteph::TaylorInterpolant; tc::Real=1.0, autodiff::Bool=true, 
                           tord::Int=10, niter::Int=5)
    # Check that radarobsfile is a file 
    @assert isfile(radarobsfile) "Cannot open file: $radarobsfile"
    # Number of massive bodies 
    Nm1 = (size(ss16asteph.x)[2]-13) ÷ 6
    # Number of bodies, including NEA
    N = Nm1 + 1
    # TODO: check that first and last observation times are within interpolation interval
    # asteroid/small-body
    # Change t, x, v units, resp., from days, au, au/day to sec, km, km/sec
    asteph_ord = asteph.x[1].order
    asteph_t0 = asteph.t0*daysec
    asteph_t = asteph.t*daysec
    asteph_r = au*map(x->x(Taylor1(asteph_ord)/daysec), asteph.x[:,1:3])
    asteph_v = (au/daysec)*map(x->x(Taylor1(asteph_ord)/daysec), asteph.x[:,4:6])
    asteph_x = hcat(asteph_r, asteph_v)
    asteph_et = TaylorInterpolant(asteph_t0, asteph_t, asteph_x)
    # Sun (su=1)
    # Change x, v units, resp., from au, au/day to km, km/sec
    sseph_t0 = ss16asteph.t0
    sseph_t = ss16asteph.t
    sseph_x = ss16asteph.x
    sun_r = au*sseph_x[:,nbodyind(Nm1,su)[1:3]]
    sun_v = (au/daysec)*sseph_x[:,nbodyind(Nm1,su)[4:6]]
    sun_x = hcat(sun_r, sun_v)
    sun_et = TaylorInterpolant(sseph_t0, sseph_t, sun_x)
    # Earth (ea=4)
    # Change x, v units, resp., from au, au/day to km, km/sec
    earth_r = au*sseph_x[:,nbodyind(Nm1,ea)[1:3]]
    earth_v = (au/daysec)*sseph_x[:,nbodyind(Nm1,ea)[4:6]]
    earth_x = hcat(earth_r, earth_v)
    earth_et = TaylorInterpolant(sseph_t0, sseph_t, earth_x)
    # Construct DataFrame with delay/doppler data from JPL radar obs file, including delay/doppler ephemeris (i.e., predicted values)
    deldop_table_jld = delay_doppler(radarobsfile, niter, xve=earth_et, xvs=sun_et, xva=asteph_et, tc=tc, tord=tord, autodiff=autodiff)
    # Save data to file
    println("Saving data to file: $outfilename")
    jldopen(outfilename, "w") do file
        addrequire(file, DataFrames)         # Require DataFrames 
        addrequire(file, TaylorSeries)       # Require TaylorSeries 
        # Write variables to jld file
        JLD.write(file, "deldop_table", deldop_table_jld)
    end
    return nothing
end