using ArgParse, NEOs, PlanetaryEphemeris, TaylorSeries
using Dates, BenchmarkTools, Printf

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "propagation.jl"
    # Desciption (for help screen)
    s.description = "NEOs propagation benchmark (Apophis)"

    s.epilog = """
        Example:\n
        \n
        julia -t 5 --project propagation.jl -o benchmark.json\n
        \n
    """

    @add_arg_table! s begin
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
            default = 700
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
    printitle("NEOs propagation benchmark (Apophis)", "=")

    # Number of workers and threads
    println("• Detected 1 worker with ", Threads.nthreads(), " thread(s)")

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

    # Parameters
    params = Parameters(
        order = 15, abstol = 1E-12, parse_eqs = true
    )
    # Time span
    jd0 = PE.J2000
    nyears = 25.0
    # Initial condition
    set_variables(Float64, "dx"; order = 2, numvars = 6)
    q00 = [-1.045062875223473E+00, -1.294565996082367E-01, -7.496573820257184E-02,
           +4.232752764404431E-03, -1.412783025595556E-02, -5.148374014688117E-03]
           # -2.901766637153165E-14, 5.E-13, 0.0
    q0 = q00 + get_variables(Float64, 2)
    # Compilation run
    params = Parameters(params; maxsteps = 1)
    NEOs.propagate(gravityonly!, q0, jd0, nyears, params)
    # Benchmark
    params = Parameters(params; maxsteps = 10_000)
    b = @benchmark NEOs.propagate($gravityonly!, $q0, $jd0, $nyears, $params) samples = samples seconds = seconds
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
