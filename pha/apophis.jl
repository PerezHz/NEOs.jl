using ArgParse, NEOs, PlanetaryEphemeris, Dates, TaylorIntegration, JLD2

# Load JPL ephemeris 
loadjpleph()

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "apophis.jl"  
    # Desciption (for help screen)
    s.description = "Propagates Apophis orbit via jet transport" 

    @add_arg_table! s begin
        "--jd0"
            help = "Initial date"
            arg_type = DateTime
            default = DateTime(2020, 12, 17)
        "--varorder"
            help = "Order of the jet transport perturbation" 
            arg_type = Int
            default = 5
        "--maxsteps"
            help = "Maximum number of steps during integration"
            arg_type = Int
            default = 10_000 
        "--nyears_bwd"
            help = "Years in backward integration"
            arg_type = Float64
            default = -18.0
        "--nyears_fwd"
            help = "Years in forward integration"
            arg_type = Float64
            default = 9.0 
        "--dense"
            help = "Whether to output the Taylor polynomial solutions obtained at each time step or not"
            arg_type = Bool
            default = true
        "--order"
            help = "Order of Taylor polynomials expansions during integration"
            arg_type = Int
            default = 25
        "--abstol"
            help = "Absolute tolerance"
            arg_type = Float64
            default = 1.0E-20
        "--parse_eqs"
            help = "Whether to use the taylorized method of jetcoeffs or not"
            arg_type = Bool
            default = true 
        "--ss_eph_file"
            help = "Path to local Solar System ephemeris file"
            arg_type = String 
            default = "./sseph343ast016_p31y_et.jld2"
    end

    s.epilog = """
        examples:\n
        \n
        # Multi-threaded\n
        julia-1.6 -t 4 --project apophis.jl --maxsteps 100 --nyears_bwd -0.02 --nyears_fwd 0.02 --parse_eqs true\n
        \n
        # Single-threaded\n
        julia-1.6 --project apophis.jl --maxsteps 100 --nyears_bwd -0.02 --nyears_fwd 0.02 --parse_eqs true\n
        \n
    """

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    # Number of threads 
    N_threads = Threads.nthreads()
    println("Number of threads: ", N_threads)
    # Dynamical function 
    if N_threads == 1
        dynamics = RNp1BP_pN_A_J23E_J2S_ng_eph!
    else 
        dynamics = RNp1BP_pN_A_J23E_J2S_ng_eph_threads!
    end 
    println("Dynamical function: ", dynamics)
    # Solar system ephemeris 
    println("Loading Solar System ephemeris...")
    ss16asteph_et = JLD2.load(parsed_args["ss_eph_file"], "ss16ast_eph")
    println("Solar System ephemeris loaded")

    # Perturbation to nominal initial condition (Taylor1 jet transport)
    dq = set_variables("δx", order = parsed_args["varorder"], numvars = 8)
    for i in 1:6
        dq[i][1][i] = 1e-8
    end
    dq[7][1][7] = 1e-14
    dq[8][1][8] = 1e-13

    # Initial date 
    jd0_datetime = parsed_args["jd0"]
    jd0 = datetime2julian(jd0_datetime)

    # Initial conditions from Apophis JPL solution #197
    q00 = kmsec2auday(apophis_pv_197(datetime2et(jd0)))
    q0 = vcat(q00, 0.0, 0.0) .+ dq

    objname = "Apophis"
    maxsteps = parsed_args["maxsteps"]
    nyears_bwd = parsed_args["nyears_bwd"]
    nyears_fwd = parsed_args["nyears_fwd"]
    dense = parsed_args["dense"]
    order = parsed_args["order"]
    abstol = parsed_args["abstol"]
    parse_eqs = parsed_args["parse_eqs"]

    println("*** Integrator warmup")
    sol = NEOs.propagate(dynamics, 1, jd0, nyears_fwd, ss16asteph_et, q0, Val(dense); 
                         order = order, abstol = abstol, parse_eqs = parse_eqs)
    println("*** Finished warmup")
    
    println("*** Full jet transport integration")
    tmax = nyears_bwd*yr 
    println("Initial time of integration: ", string(jd0_datetime))
    println("Final time of integration: ", julian2datetime(jd0 + tmax))

    sol = NEOs.propagate(dynamics, maxsteps, jd0, nyears_bwd, ss16asteph_et, q0, Val(dense); 
                         order = order, abstol = abstol, parse_eqs = parse_eqs)
    save2jldandcheck(objname*"_bwd", (asteph = sol,))

    tmax = nyears_fwd*yr 
    println("Initial time of integration: ", string(jd0_datetime))
    println("Final time of integration: ", julian2datetime(jd0 + tmax))

    sol = NEOs.propagate(dynamics, maxsteps, jd0, nyears_fwd, ss16asteph_et, q0, Val(dense), 
                         order = order, abstol = abstol, parse_eqs = parse_eqs)
    save2jldandcheck(objname*"_fwd", (asteph = sol,))
    println("*** Finished asteroid ephemeris integration")
    
    nothing 
end

main()
