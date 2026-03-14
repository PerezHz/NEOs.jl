using ArgParse, NEOs, PlanetaryEphemeris, TaylorSeries, Dates, JLD2, Printf, BenchmarkTools
using NEOs: ImpactMonitoringBuffer, radius

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "impactmonitoring.jl"
    # Desciption (for help screen)
    s.description = "NEOs impact monitoring benchmark (2022 NX1)"

    s.epilog = """
        Example:\n
        \n
        julia -t 5 --project impactmonitoring.jl -i 2022NX1.jld2\n
        \n
    """

    @add_arg_table! s begin
        "--input", "-i"
            help = "input orbit determination file"
            arg_type = String
        "--output", "-o"
            help = "output benchmark file"
            arg_type = String
        "--samples", "-s"
            help = "number of samples to take"
            arg_type = Int
            default = 10
        "--seconds", "-t"
            help = "number of seconds budgeted for the benchmarking process"
            arg_type = Int
            default = 250
    end

    return parse_args(s)
end

computationtime(x::DateTime, y::DateTime) = @sprintf("%.2f", (y - x).value / 60_000)

printitle(s::AbstractString, d::AbstractString) = println(d ^ length(s),
        '\n', s, '\n', d ^ length(s))

function main()
    # Parse arguments from commandline
    parsed_args = parse_commandline()

    # Print header
    printitle("NEOs impact monitoring benchmark (2022 NX1)", "=")

    # Number of workers and threads
    println("• Detected 1 worker with ", Threads.nthreads(), " thread(s)")

    # Input orbit determination file
    input::String = parsed_args["input"]
    println("• Input orbit determination file: ", input)

    # Output benchmark file
    output::String = parsed_args["output"]
    println("• Output benchmark file: ", output)

    # Number of samples to take
    samples::Int = parsed_args["samples"]
    println("• Number of samples to take: ", samples)

    # Number of seconds budgeted for the benchmarking process
    seconds::Int = parsed_args["seconds"]
    println("• Number of seconds budgeted for the benchmarking process: ", seconds)

    # Global initial time
    global_initial_time = now()
    println("• Run started at ", global_initial_time)

    # Load orbit, line of variations and parameters
    orbit = JLD2.load(input, "orbit")
    lov = JLD2.load(input, "lov")
    params = Parameters(JLD2.load(input, "params"); order = 15, abstol = 1E-12)

    # Impact monitoring problem
    IM = IMProblem(orbit, ImpactTarget(:earth))

    # Virtual asteroids
    R_TP, R_P, IP, Δσmax, vaorder = 0.2, radius(IM.target), 5E-6, 0.05, 6
    VAs = virtualasteroids(lov, :DelVigna19, vaorder; R_TP, R_P, IP, Δσmax)
    VA = VAs[1]

    # Close approaches
    ctol = 1.0
    jd0 = epoch(lov) + PE.J2000
    nyears = ( datetime2julian(DateTime(2099, 12, 31)) - jd0 ) / yr
    q0 = orbit() .+ 1E-8 * Taylor1(vaorder)
    buffer = ImpactMonitoringBuffer(IM, q0, nyears, params)
    # Compilation run
    closeapproaches(IM, VA, nyears, params; R_TP, ctol, buffer)
    # Benchmark
    b = @benchmark closeapproaches($IM, $VA, $nyears, $params; R_TP = $R_TP, ctol = $ctol, buffer = $buffer) samples = samples seconds = seconds
    show(stdout, "text/plain", b)

    # Save results
    BenchmarkTools.save(output, b)
    println("• Output saved to: ", output)

    # Final time
    global_final_time = now()
    println("• Run started ", global_initial_time, " and finished ", global_final_time)
    global_computation_time = computationtime(global_initial_time, global_final_time)
    println("• Total computation time was: ", global_computation_time, " min")

    return nothing
end

main()