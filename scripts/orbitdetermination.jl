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
    using NEOs: Tracklet, reduce_tracklets, isgauss

    # Initial orbit determination routines

    # First filter
    function iod1(neo::String, outdir::String)
        # Output file
        filename = joinpath(outdir, replace(neo, " " => "_") * ".jld2")
        # Download optical astrometry
        radec = fetch_radec_mpc("designation" => neo)
        if length(radec) < 3
            jldsave(filename; tracklets = Vector{Tracklet}(undef, 0))
            return false
        end
        # Parameters
        params = NEOParameters(coeffstol = 10.0, bwdoffset = 0.007,
                 fwdoffset = 0.007, jtlsiter = 10, adamhelp = false)
        # Select at most three tracklets
        tracklets = reduce_tracklets(radec)
        if length(tracklets) > 3
            tracklets = tracklets[1:3]
            radec = reduce(vcat, getfield.(tracklets, :radec))
            sort!(radec)
        end

        try
            # Start of computation
            init_time = now()
            # Orbit determination
            sol = orbitdetermination(radec, params)
            # Time of computation
            Δ = (now() - init_time).value
            # Save orbit
            if length(sol.res) != length(radec)
                jldsave(filename; tracklets = tracklets, Δ = Δ)
                return false
            else
                jldsave(filename; sol = sol, Δ = Δ)
                return true
            end
        catch
            # An error ocurred
            jldsave(filename; tracklets = tracklets)
            return false
        end
    end

    # Second filter
    function iod2(neo::String, outdir::String)
        # Output from last filter
        filename = joinpath(outdir, replace(neo, " " => "_") * ".jld2")
        dict = JLD2.load(filename)
        # Previous filter already computed an orbit
        "sol" in keys(dict) && return true
        # Check tracklets are non empty
        tracklets = dict["tracklets"]
        isempty(tracklets) && return false
        # Computation time so far
        if "Δ" in keys(dict)
            Δ = dict["Δ"]
        else
            Δ = 0
        end
        # Optical astrometry
        radec = reduce(vcat, getfield.(tracklets, :radec))
        sort!(radec)
        # Parameters
        params = NEOParameters(coeffstol = Inf, bwdoffset = 0.5,
                 fwdoffset = 0.5, jtlsiter = 10, adamhelp = true)

        try
            # Start of computation
            init_time = now()
            # Orbit determination
            sol = orbitdetermination(radec, params)
            # Time of computation
            Δ += (now() - init_time).value
            # Save orbit
            if length(sol.res) != length(radec)
                jldsave(filename; tracklets = tracklets, Δ = Δ)
                return false
            else
                jldsave(filename; sol = sol, Δ = Δ)
                return true
            end
        catch
            # An error ocurred
            jldsave(filename; tracklets = tracklets)
            return false
        end
    end

    # Third filter
    function iod3(neo::String, outdir::String)
        # Output from last filter
        filename = joinpath(outdir, replace(neo, " " => "_") * ".jld2")
        dict = JLD2.load(filename)
        # Previous filter already computed an orbit
        "sol" in keys(dict) && return true
        # Check tracklets are non empty
        tracklets = dict["tracklets"]
        isempty(tracklets) && return false
        # Computation time so far
        if "Δ" in keys(dict)
            Δ = dict["Δ"]
        else
            Δ = 0
        end
        # Optical astrometry
        radec = reduce(vcat, getfield.(tracklets, :radec))
        sort!(radec)
        # Parameters
        params = NEOParameters(coeffstol = Inf, bwdoffset = 0.5,
                 fwdoffset = 0.5, jtlsiter = 10, adamhelp = true)

        try
            # Start of computation
            init_time = now()
            # Orbit determination
            if isgauss(tracklets)
                sol = tooshortarc(radec, tracklets, params)
            else
                sol = gaussinitcond(radec, tracklets, params)
            end
            # Time of computation
            Δ += (now() - init_time).value
            # Save orbit
            if length(sol.res) != length(radec)
                jldsave(filename; tracklets = tracklets, Δ = Δ)
                return false
            else
                jldsave(filename; sol = sol, Δ = Δ)
                return true
            end
        catch
            # An error ocurred
            jldsave(filename; tracklets = tracklets)
            return false
        end
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

    # First filter
    mask = pmap(neo -> iod1(neo, outdir), neos)
    all(mask) && return nothing
    # Second filter
    mask = pmap(neo -> iod2(neo, outdir), neos)
    all(mask) && return nothing
    # Third filter
    mask = pmap(neo -> iod3(neo, outdir), neos)
    return nothing
end

main()