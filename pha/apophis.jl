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
            default = -17.0
        "--nyears_fwd"
            help = "Years in forward integration"
            arg_type = Float64
            default = 9.0
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
    end

    s.epilog = """
        examples:\n
        \n
        # Multi-threaded\n
        julia -t 4 --project apophis.jl --maxsteps 100 --nyears_bwd -0.02 --nyears_fwd 0.02 --parse_eqs true\n
        \n
        # Single-threaded\n
        julia --project apophis.jl --maxsteps 100 --nyears_bwd -0.02 --nyears_fwd 0.02 --parse_eqs true\n
        \n
    """

    return parse_args(s)
end

function print_header(header::String, level::Int = 1)
    L = length(header)
    if level == 1
        c = "="
    else
        c = "-"
    end
    println(repeat(c, L))
    println(header)
    println(repeat(c, L))
end

function main(dynamics::D, maxsteps::Int, jd0_datetime::DateTime, nyears_bwd::T, nyears_fwd::T, order::Int, varorder::Int,
              abstol::T, parse_eqs::Bool) where {T <: Real, D}

    # Initial conditions from Apophis JPL solution #197
    q00 = kmsec2auday(apophisposvel197(datetime2et(jd0_datetime)))

    # Perturbation to nominal initial condition (Taylor1 jet transport)
    # vcat(fill(1e-8, 6), 1e-14, 1e-13) are the scaling factors for jet transport perturbation,
    # these are needed to ensure expansion coefficients remain small.
    # The magnitudes correspond to the typical order of magnitude of errors in
    # position/velocity (1e-8), Yarkovsky (1e-13) and radiation pressure (1e-14)
    if varorder == 0
        dq = zeros(8)
    else
        dq = NEOs.scaled_variables("δx", vcat(fill(1e-8, 6), 1e-14), order = varorder)
    end

    q0 = vcat(q00, 0.0, 0.0) .+ vcat(dq, zero(dq[1]))

    # Initial date (in julian days)
    jd0 = datetime2julian(jd0_datetime)

    print_header("Integrator warmup", 2)
    _ = NEOs.propagate(dynamics, 1, jd0, nyears_fwd, q0, Val(true);
                       order = order, abstol = abstol, parse_eqs = parse_eqs)
    println()

    print_header("Main integration", 2)
    tmax = nyears_bwd*yr
    println("• Initial time of integration: ", string(jd0_datetime))
    println("• Final time of integration: ", julian2datetime(jd0 + tmax))

    sol_bwd = NEOs.propagate(dynamics, maxsteps, jd0, nyears_bwd, q0, Val(true);
                             order = order, abstol = abstol, parse_eqs = parse_eqs)
    PE.save2jld2andcheck("Apophis_bwd.jld2", (asteph = sol_bwd,))

    tmax = nyears_fwd*yr
    println("• Initial time of integration: ", string(jd0_datetime))
    println("• Final time of integration: ", julian2datetime(jd0 + tmax))

    sol_fwd = NEOs.propagate(dynamics, maxsteps, jd0, nyears_fwd, q0, Val(true);
                             order = order, abstol = abstol, parse_eqs = parse_eqs)
    PE.save2jld2andcheck("Apophis_fwd.jld2", (asteph = sol_fwd,))

    println()

    # NEO
    # Change t, x, v units, resp., from days, au, au/day to sec, km, km/sec
    xva_bwd(et) = auday2kmsec(sol_bwd(et/daysec)[1:6])
    xva_fwd(et) = auday2kmsec(sol_fwd(et/daysec)[1:6])
    xva(et) = bwdfwdeph(et, sol_bwd, sol_fwd)
    # Earth
    # Change x, v units, resp., from au, au/day to km, km/sec
    eph_ea = selecteph(NEOs.sseph, ea)
    xve(et) = auday2kmsec(eph_ea(et/daysec))
    # Sun
    # Change x, v units, resp., from au, au/day to km, km/sec
    eph_su = selecteph(NEOs.sseph, su)
    xvs(et) = auday2kmsec(eph_su(et/daysec))


    radec_2004_2020 = read_radec_mpc("/Users/Jorge/projects/NEOs/data/99942_2004_2020.dat")
    radec_2020_2021 = read_radec_mpc("/Users/Jorge/projects/NEOs/data/99942_2020_2021.dat")
    radec = vcat(radec_2004_2020,radec_2020_2021)

    deldop_2005_2013 = read_radar_jpl("/Users/Jorge/projects/NEOs/data/99942_RADAR_2005_2013.dat")
    deldop_2021 = read_radar_jpl("/Users/Jorge/projects/NEOs/data/99942_RADAR_2021.dat")
    # TODO: make radar astrometry residuals work with functions!!!
    # deldop = vcat(deldop_2005_2013,deldop_2021)

    # Compute optical residuals
    print_header("Compute residuals", 2)
    res_radec, w_radec = residuals(radec, xvs = xvs, xve = xve, xva = xva)

    # Compute radar residuals
    #### TODO: make radar astrometry residuals work with functions!!!
    #### res_deldop, w_deldop = residuals(deldop, xvs = xvs, xve = xve, xva = xva_bwd)
    res_del_bwd, w_del_bwd, res_dop_bwd, w_dop_bwd = residuals(deldop_2005_2013, xvs = xvs, xve = xve, xva = xva_bwd, niter=5)
    res_del_fwd, w_del_fwd, res_dop_fwd, w_dop_fwd = residuals(deldop_2021, xvs = xvs, xve = xve, xva = xva_fwd, niter=5)

    PE.save2jld2andcheck("resw_radec.jld2", (;res_radec,w_radec,res_del_bwd,w_del_bwd,res_dop_bwd,w_dop_bwd,res_del_fwd,w_del_fwd,res_dop_fwd,w_dop_fwd))

    @show NEOs.cte(res) NEOs.cte(w)

    nothing

end

function main()

    # Parse arguments from commandline
    parsed_args = parse_commandline()

    print_header("Asteroid Apophis")
    println()
    print_header("General parameters", 2)

    # Number of threads
    N_threads = Threads.nthreads()
    println("• Number of threads: ", N_threads)

    # Dynamical function
    dynamics = RNp1BP_pN_A_J23E_J2S_ng_eph_threads!

    println("• Dynamical function: ", dynamics)

    # Maximum number of steps
    maxsteps = parsed_args["maxsteps"]
    println("• Maximum number of steps: ", maxsteps)

    # Initial date
    jd0_datetime = parsed_args["jd0"]

    # Number of years in backward integration
    nyears_bwd = parsed_args["nyears_bwd"]

    # Number of years in forward integration
    nyears_fwd = parsed_args["nyears_fwd"]

    # Order of Taylor polynomials
    order = parsed_args["order"]
    println("• Order of Taylor polynomials: ", order)

    # Order of jet transport perturbation
    varorder = parsed_args["varorder"]
    println("• Order of jet transport perturbation: ", varorder)

    # Absolute tolerance
    abstol = parsed_args["abstol"]
    println("• Absolute tolerance: ", abstol)

    # Whether to use @taylorize or not
    parse_eqs = parsed_args["parse_eqs"]
    println("• Use @taylorize: ", parse_eqs, "\n")

    main(dynamics, maxsteps, jd0_datetime,  nyears_bwd, nyears_fwd, order, varorder, abstol, parse_eqs)
end

main()
