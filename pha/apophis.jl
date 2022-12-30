# Multi-threaded:
# julia -t <number-of-threads> --project=@. pha/apoophis.jl --help
# Single-threaded:
# julia --project=@. pha/apoophis.jl --help  

using ArgParse, NEOs, PlanetaryEphemeris, Dates, TaylorIntegration, JLD

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "apophis.jl"  
    # Desciption (for help screen)
    s.description = "Propagates Apophis orbit via jet transport" 

    @add_arg_table! s begin
        "--varorder"
            help = "Order corresponding to the jet transport perturbation" 
            arg_type = Int
            default = 5
        "--objname"
            help = "Name of the object"
            arg_type = String
            default = "Apophis"
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
            default = false
        "--quadmath"
            help = "Whether to use quadruple precision or not"
            arg_type = Bool
            default = false
        "--lyap"
            help = "Whether to compute Lyapunov exponents or not"
            arg_type = Bool
            default = false
        "--order"
            help = "Order of Taylor polynomials expansions during integration"
            arg_type = Int
            default = 25
        "--abstol"
            help = "Absolute tolerance"
            arg_type = Float64
            default = 1.0E-20
        "--ss_eph_file"
            help = "Path to local Solar System ephemeris file"
            arg_type = String 
            default = "./sseph343ast016_p31y_et.jld"
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    
    println("Number of threads: ", Threads.nthreads())

    # Number of TaylorN variables
    nv = 8 
    # Dynamical function 
    dynamics = RNp1BP_pN_A_J23E_J2S_ng_eph_threads!

    ss16asteph_et = JLD.load(parsed_args["ss_eph_file"], "ss16ast_eph")

    # TaylorN variables setup
    if parsed_args["lyap"]
        # Setup TaylorN variables with order = 1, numvars = nv
        TNvars = set_variables("δx", order = 1, numvars = nv)
        # dq corresponding to solution OR7 for Apophis
        dq_OR7 = [-8.053250543083672e-7, -6.4934453239292154e-9, -3.552581604396334e-8, 2.382431039885935e-9, -1.3384789262277344e-8, -4.6746457798167725e-9, -2.892614243659006e-14, 0.0]
        dq = dq_OR7
    else
        # dq: perturbation to nominal initial condition (Taylor1 jet transport)
        dq = set_variables("δx", order = parsed_args["varorder"], numvars = nv)
        for i in 1:6
            dq[i][1][i] = 1e-8
        end
        if get_numvars() == 8
            dq[7][1][7] = 1e-14
            dq[8][1][8] = 1e-13
        end
    end

    # Initial conditions from Apophis JPL solution #197 at the Dec-17-2020.0 (TDB) epoch
    jd0 = datetime2julian(DateTime(2020,12,17)) # JDTDB = 2459200.5
    q00 = [-0.18034747703273316, 0.9406910666200128, 0.3457360259054398, -0.016265942170279046, 4.392889725556651e-5, -0.00039519931615139716] ### JPL solution #197 at 2020Dec17.0 (TDB)
    q0 = vcat(q00, 0.0, 0.0) .+ dq

    println("*** Integrator warmup")
    NEOs.propagate(parsed_args["objname"], dynamics, 1, jd0, parsed_args["nyears_fwd"], ss16asteph_et, Val(parsed_args["quadmath"]),
                   Val(parsed_args["dense"]), Val(parsed_args["lyap"]), output = false, q0 = q0, order = parsed_args["order"], 
                   abstol = parsed_args["abstol"])
    println("*** Finished warmup")

    println("*** Full jet transport integration")
    NEOs.propagate(parsed_args["objname"]*"_bwd", dynamics, parsed_args["maxsteps"], jd0, parsed_args["nyears_bwd"], ss16asteph_et, 
                   Val(parsed_args["quadmath"]), Val(parsed_args["dense"]), Val(parsed_args["lyap"]), q0 = q0, 
                   order = parsed_args["order"], abstol = parsed_args["abstol"])
    NEOs.propagate(parsed_args["objname"]*"_fwd", dynamics, parsed_args["maxsteps"], jd0, parsed_args["nyears_fwd"], ss16asteph_et, 
                   Val(parsed_args["quadmath"]), Val(parsed_args["dense"]), Val(parsed_args["lyap"]), q0 = q0, 
                   order = parsed_args["order"], abstol = parsed_args["abstol"])
    println("*** Finished asteroid ephemeris integration")

end

main()
