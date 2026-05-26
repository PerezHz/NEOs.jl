using Distributed, ArgParse

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "2024YR4IM.jl"
    # Desciption (for help screen)
    s.description = "Progressive impact monitoring for asteroid 2024 YR4"

    s.epilog = """
        Example:\n
        \n
        julia -p 10 -t 5 --project 2024YR4IM.jl -i 2024YR4Orbits.jld2 -o 2024YR4IM.jld2\n
        \n
    """

    @add_arg_table! s begin
        "--input", "-i"
            help = "input orbit determination file"
            arg_type = String
        "--output", "-o"
            help = "output impact monitoring file"
            arg_type = String
        "--smax"
            help = "maximum (absolute) LOV index"
            arg_type = Float64
            default = 3.0
        "--lovorder"
            help = "order of expansions wrt LOV index"
            arg_type = Int
            default = 4
        "--lovtol"
            help = "LOV integration absolute tolerance"
            arg_type = Float64
            default = 1E-8
        "--lovsteps"
            help = "LOV integration maximum number of steps"
            arg_type = Int
            default = 10_000
        "--target"
            help = "impact target"
            arg_type = String
            default = "earth"
        "--generic"
            help = "generic completeness limit"
            arg_type = Float64
            default = 5E-6
        "--dsmax"
            help = "maximum interval length"
            arg_type = Float64
            default = 0.05
        "--vaorder"
            help = "order of expansions wrt the LOV index"
            arg_type = Int
            default = 4
        "--ctol"
            help = "convergence tolerance"
            arg_type = Float64
            default = 1.0
    end

    return parse_args(s)
end

@everywhere using Printf

@everywhere begin
    using NEOs, PlanetaryEphemeris, TaylorSeries, Dates, JLD2, ThreadPinning
    using NEOs: AbstractOrbit, ImpactMonitoringBuffer, radius, nominaltime

    computationtime(x::DateTime, y::DateTime) = @sprintf("%.2f", (y - x).value / 60_000)

    printitle(io::IO, s::AbstractString, d::AbstractString) = println(io, d ^ length(s),
        '\n', s, '\n', d ^ length(s))
    printitle(s::AbstractString, d::AbstractString) = printitle(stdout, s, d)

    function impactmonitoring(i::Int, orbit::AbstractOrbit)
        # Orbit initial time
        orbit_initial_time = now()
        printitle(io, "Orbit $i/$Norbits", "-")
        println(io, "• Run started at ", orbit_initial_time)

        # Impact monitoring problem
        IM = IMProblem(orbit, ImpactTarget(target))

        # Line of variations
        lov = lineofvariations(IM, params; coord = :cartesian, σmax, lovorder,
                               lovtol, lovsteps)

        # Virtual asteroids
        R_TP, R_P = 0.2, radius(IM.target)
        VAs = virtualasteroids(lov, :DelVigna19, vaorder; R_TP, R_P, IP, Δσmax)
        println(io, "• Sampling consists of ", length(VAs), " virtual asteroids")

        # Close approaches
        jd0 = epoch(lov) + PE.J2000
        nyears = ( datetime2julian(DateTime(2100, 1, 1, 12)) - jd0 ) / yr
        q0 = orbit() .+ 1E-8 * Taylor1(vaorder)
        buffer = ImpactMonitoringBuffer(IM, q0, nyears, params)
        CAs = Vector{Vector{CloseApproach{Float64, Taylor1{Float64}}}}(undef, length(VAs))
        for i in eachindex(CAs)
            CAs[i] = closeapproaches(IM, VAs[i], nyears, params; R_TP, ctol, buffer)
        end

        # Returns
        RTs = showersnreturns(CAs)

        # Virtual impactors
        no_pts, dmax = 100, 10.0
        VIs = virtualimpactors(IM, lov, RTs, params; ctol, no_pts, dmax)
        # Sort by time of impact
        sort!(VIs, by = nominaltime)
        # Eliminate marginal and zero probability virtual impactors
        # filter!(x -> impact_probability(x) > 0, VIs)
        # Print impactor table
        println(io, '\n', replace(summary(VIs), r"^"m => '\t'))

        # Orbit final time
        orbit_final_time = now()
        println(io, "• Run finished at ", orbit_final_time)
        orbit_computation_time = computationtime(orbit_initial_time, orbit_final_time)
        println(io, "• Orbit computation time was ", orbit_computation_time, " min\n")

        return lov, CAs, VIs, String(take!(io))
    end
end

function main()
    # Parse arguments from commandline
    parsed_args = parse_commandline()

    # Print header
    printitle("Progressive impact monitoring for asteroid 2024 YR4", "=")

    # Number of workers and threads
    println("• Detected ", nworkers(), " worker(s) with ", Threads.nthreads(),
            " thread(s) each")

    # Pin threads to sockets
    distributed_pinthreads(:sockets)

    # Input orbit determination file
    input::String = parsed_args["input"]
    println("• Input orbit determination file: ", input)

    # Output impact monitoring file
    output::String = parsed_args["output"]
    println("• Output impact monitoring file: ", output)

    # Maximum (absolute) LOV index
    σmax::Float64 = parsed_args["smax"]
    println("• Maximum (absolute) LOV index: ", σmax)

    # Order of expansions wrt LOV index
    lovorder::Int = parsed_args["lovorder"]
    println("• Order of expansions wrt LOV index: ", lovorder)

    # LOV integration absolute tolerance
    lovtol::Float64 = parsed_args["lovtol"]
    println("• LOV integration absolute tolerance: ", lovtol)

    # LOV integration maximum number of steps
    lovsteps::Int = parsed_args["lovsteps"]
    println("• LOV integration maximum number of steps: ", lovsteps)

    # Impact target
    target::String = parsed_args["target"]
    println("• Impact target: ", target)

    # Generic completeness limit
    IP::Float64 = parsed_args["generic"]
    println("• Generic completeness limit: ", IP)

    # Maximum interval length
    Δσmax::Float64 = parsed_args["dsmax"]
    println("• Maximum interval length: ", Δσmax)

    # Order of expansions wrt the LOV index
    vaorder::Int = parsed_args["vaorder"]
    println("• Order of expansions wrt the LOV index: ", vaorder)

    # Convergence tolerance
    ctol::Float64 = parsed_args["ctol"]
    println("• Convergence tolerance: ", ctol)

    # Load orbits
    orbits = JLD2.load(input, "orbits")

    # Define global constants
    @everywhere begin
        const IP = $IP
        const σmax = $σmax
        const ctol = $ctol
        const Δσmax = $Δσmax
        const io = IOBuffer()
        const lovtol = $lovtol
        const target = $target
        const vaorder = $vaorder
        const lovsteps = $lovsteps
        const lovorder = $lovorder
        const Norbits = length($orbits)
        const params = Parameters(
            maxsteps = 10_000, order = 15, abstol = 1E-12, parse_eqs = true,
            coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
            gaussorder = 2, safegauss = false, refscale = :log,
            tsaorder = 2, adamiter = 500, adamQtol = 1E-5,
            jtlsorder = 2, jtlsmask = false, jtlsiter = 20, lsiter = 20,
            jtlsproject = false, significance = 0.99, verbose = true,
            outrej = false, χ2_rec = 7.0, χ2_rej = 8.0, fudge = 100.0,
            max_per = 33.3
        )
    end

    # Global initial time
    global_initial_time = now()
    println("• Run started at ", global_initial_time)

    # Impact monitoring
    IM = pmap(impactmonitoring, eachindex(orbits), orbits)

    # Save output
    lovss, CAss, VIss, IOss = first.(IM), getindex.(IM, 2), getindex.(IM, 3), last.(IM)
    join(stdout, IOss)
    jldsave(output; lovss, CAss, VIss)
    println("• Output saved to: ", output)

    # Final time
    global_final_time = now()
    println("• Run started ", global_initial_time, " and finished ", global_final_time)
    global_computation_time = computationtime(global_initial_time, global_final_time)
    println("• Total computation time was: ", global_computation_time, " min")

    return nothing
end

main()