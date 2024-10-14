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
            help = "input names file"
            arg_type = String
            default = "names.txt"
        "--output", "-o"
            help = "output directory"
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
    using NEOs: AdmissibleRegion, RadecMPC, reduce_tracklets, numberofdays,
        issatellite

    function radecfilter(radec::Vector{RadecMPC{T}}) where {T <: Real}
        # Eliminate observations before oficial discovery
        firstobs = findfirst(r -> !isempty(r.discovery), radec)
        isnothing(firstobs) && return false, radec
        radec = radec[firstobs:end]
        # Filter out incompatible observations
        filter!(radec) do r
            hascoord(r.observatory) && !issatellite(r.observatory) &&
            date(r) >= Date(2000)
        end
        length(radec) < 3 && return false, radec
        # Find the first set of 3 tracklets with a < 15 days timespan
        tracklets = reduce_tracklets(radec)
        for i in 1:length(tracklets)-2
            numberofdays(tracklets[i:i+2]) > 15.0 && continue
            tracklets = tracklets[i:i+2]
            radec = reduce(vcat, getfield.(tracklets, :radec))
            sort!(radec)
            break
        end
        return numberofdays(radec) <= 15.0, radec
    end

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
        radec = fetch_radec_mpc(neo)
        # Get at most 3 tracklets for orbit determination
        flag, radec = radecfilter(radec)
        !flag && return false
        # Orbit determination problem
        od = ODProblem(newtonian!, radec)
        # Parameters
        params = NEOParameters(
            coeffstol = Inf, bwdoffset = 0.042, fwdoffset = 0.042, # Propagation
            gaussorder = 6, gaussQmax = 2.0,                       # Gauss method
            adamiter = 500, adamQtol = 1e-5, tsaQmax = 2.0,        # ADAM
            jtlsiter = 20, lsiter = 10,                            # Least squares
            outrej = true, χ2_rec = 7.0, χ2_rej = 8.0,             # Outlier rejection
            fudge = 10.0, max_per = 34.0
        )
        # Start of computation
        init_time = now()
        # Initial orbit determination
        sol = orbitdetermination(od, params; initcond)
        # Termination condition
        if length(sol.res) == length(radec) && nrms(sol) < 2.0
            # Time of computation
            Δ = (now() - init_time).value
            # Save orbit
            jldsave(filename; sol = sol, Δ = Δ)
            # Sucess flag
            return true
        end
        # Parameters
        params = NEOParameters(params;
            coeffstol = 10.0, adamiter = 200, adamQtol = 0.01,
            lsiter = 5
        )
        # Initial orbit determination
        _sol_ = orbitdetermination(od, params; initcond)
        # Time of computation
        Δ = (now() - init_time).value
        # Choose best orbit
        if length(sol.res) == length(_sol_.res) == length(radec)
            sol = min(sol, _sol_)
        elseif length(sol.res) == length(radec)
            sol = sol
        elseif length(_sol_.res) == length(radec)
            sol = _sol_
        else
            sol = zero(NEOSolution{Float64, Float64})
        end
        # Save orbit
        iszero(sol) && return false
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