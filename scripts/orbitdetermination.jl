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
        "--input", "-i"
            help = "Input names file"
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
        julia -p 10 -t 5 --project scripts/orbitdetermination.jl -i names.txt -o orbits/\n
        \n
    """

    return parse_args(s)
end

@everywhere begin
    using NEOs, Dates, JLD2
    using NEOs: AdmissibleRegion, reduce_tracklets

    # Default naive initial conditions for iod
    function initcond(A::AdmissibleRegion{T}) where {T <: Real}
        v_ρ = sum(A.v_ρ_domain) / 2
        return [
            (A.ρ_domain[1], v_ρ, :log),
            (10^(sum(log10, A.ρ_domain) / 2), v_ρ, :log),
            (sum(A.ρ_domain) / 2, v_ρ, :log),
            (A.ρ_domain[2], v_ρ, :log),
        ]
    end

    # Initial orbit determination routine
    function iod(neo::String, filename::String)
        # Download optical astrometry
        radec = fetch_radec_mpc("designation" => neo)
        length(radec) < 3 && return false
        # Parameters
        params = NEOParameters(coeffstol = Inf, bwdoffset = 0.007,
            fwdoffset = 0.007, adamiter = 500, adamQtol = 1e-5,
            jtlsiter = 20, newtoniter = 10)
        dynamics = newtonian!
        gauss = true
        # Select at most three tracklets
        tracklets = reduce_tracklets(radec)
        if length(tracklets) > 3
            tracklets = tracklets[1:3]
            radec = reduce(vcat, getfield.(tracklets, :radec))
            sort!(radec)
        end
        # Start of computation
        init_time = now()
        # Initial orbit determination
        sol = NEOs.iod(radec, params; gauss, dynamics, initcond)
        # Time of computation
        Δ = (now() - init_time).value
        # Unsucessful orbit determination
        length(sol.res) != length(radec) && return false
        # Save orbit
        jldsave(filename; sol = sol, Δ = Δ)

        return true
    end
end

function main()
    # Parse arguments from commandline
    parsed_args = parse_commandline()
    # Input names file
    input::String = parsed_args["input"]
    # Output directory
    output::String = parsed_args["output"]
    # Print header
    println("Initial orbit determination for NEOs via jet transport")
    println("• Input names file: ", input)
    println("• Output directory: ", output)

    # Parse NEOs' designations
    neos = readlines(input)
    println("• ", length(neos), " NEOs to be processed with ", nworkers(),
        " workers (", Threads.nthreads(), " threads each)")
    # Output files
    filenames = map(neos) do neo
        return joinpath(output, replace(neo, " " => "") * ".jld2")
    end

    # Distributed orbit determination
    mask = pmap(iod, neos, filenames; on_error = ex -> false)
    println("• ", count(mask), " / ", length(neos), " successful NEOs")

    return nothing
end

main()