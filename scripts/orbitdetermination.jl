using Distributed, ArgParse

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "orbitdetermination.jl"
    # Desciption (for help screen)
    s.description = """
    Determines the orbit of a list of NEOs via jet transport
    using Julia's interface for distributed computing."""

    @add_arg_table! s begin
        "--names", "-n"
            help = "File containing the names of the NEOs to be processed"
            arg_type = String
            default = "names.txt"
        "--output", "-o"
            help = "Output directory"
            arg_type = String
            default = pwd()
    end

    s.epilog = """
        example:\n
        \n
        # Run orbitdetermination.jl with 10 workers and 5 threads each\n
        julia -p 10 -t 5 --project scripts/orbitdetermination.jl -n names.txt -o orbits/\n
        \n
    """

    return parse_args(s)
end

@everywhere begin
    using NEOs, Dates, JLD2

    # Initial orbit determination routine
    function iod(neo::String, outdir::String)
        # Output file
        filename = joinpath(outdir, replace(neo, " " => "_") * ".jld2")

        try
            # Download optical astrometry
            radec = fetch_radec_mpc("designation" => neo)
            # Parameters
            params = NEOParameters(coeffstol = Inf, bwdoffset = 0.5, fwdoffset = 0.5)
            # Start of computation
            init_time = now()
            # Orbit determination
            sol = orbitdetermination(radec, params)
            # Time of computation
            Δ = (now() - init_time).value
            # Save orbit
            jldsave(filename; sol = sol, Δ = Δ)
        catch error
            # Save error
            jldsave(filename; error = error)
        end

        nothing
    end
end

function main()
    # Parse arguments from commandline
    parsed_args = parse_commandline()

    println("Orbit determination for NEOs via jet transport")
    println()

    # NEOs to be processed
    namesfile = parsed_args["names"]
    neos = readlines(namesfile)
    println(length(neos), " NEOs to be processed with ", nworkers(), " workers (",
            Threads.nthreads(), " threads each)")
    println()
    # Output directory
    outdir = parsed_args["output"]

    # Distributed orbit determination
    pmap(neo -> iod(neo, outdir), neos)

    nothing
end

main()
