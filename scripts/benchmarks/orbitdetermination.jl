using ArgParse, Dates, NEOs, BenchmarkTools, Printf

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "orbitdetermination.jl"
    # Desciption (for help screen)
    s.description = "NEOs orbit determination benchmark (2014 AA)"

    s.epilog = """
        Example:\n
        \n
        julia -t 5 --project orbitdetermination.jl -i 2014AA.txt\n
        \n
    """

    @add_arg_table! s begin
        "--input", "-i"
            help = "input observations file"
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
            default = 100
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
    printitle("NEOs orbit determination benchmark (2014 AA)", "=")

    # Number of workers and threads
    println("• Detected 1 worker with ", Threads.nthreads(), " thread(s)")

    # Input observations file
    input::String = parsed_args["input"]
    println("• Input observations file: ", input)

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

    # Load optical astrometry
    optical = read_optical_mpc80(input)
    # Parameters
    params = Parameters(
        coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
        tsaorder = 2, adamiter = 500, adamQtol = 1e-5,
        jtlsorder = 2, jtlsiter = 200, lsiter = 1,
        significance = 0.99, outrej = false, verbose = false
    )
    # Orbit determination problem
    od = ODProblem(newtonian!, optical)
    # Compilation run
    initialorbitdetermination(od, params)
    # Benchmark
    b = @benchmark initialorbitdetermination($od, $params) samples = samples seconds = seconds
    show(stdout, "text/plain", b)
    println()

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