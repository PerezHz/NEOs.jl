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
    using NEOs: AdmissibleRegion, RadecMPC, NEOSolution, reduce_tracklets,
        numberofdays, issatellite

    function radecfilter(radec::Vector{RadecMPC{T}}) where {T <: Real}
        # Eliminate observations before oficial discovery
        firstobs = findfirst(r -> !isempty(r.discovery), radec)
        isnothing(firstobs) && return false, radec
        radec = radec[firstobs:end]
        # Filter out incompatible observations
        filter!(radec) do r
            hascoord(r.observatory) && !issatellite(r.observatory) &&
            date(r) > DateTime(2000, 1, 1, 12)
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

    # Naive initial conditions for iod
    function initcond1(A::AdmissibleRegion{T}) where {T <: Real}
        v_ρ = sum(A.v_ρ_domain) / 2
        return [
            (A.ρ_domain[1], v_ρ, :log),
            (10^(sum(log10, A.ρ_domain) / 2), v_ρ, :log),
            (sum(A.ρ_domain) / 2, v_ρ, :log),
            (A.ρ_domain[2], v_ρ, :log),
        ]
    end

    function initcond2(A::AdmissibleRegion{T}) where {T <: Real}
        v_ρ = sum(A.v_ρ_domain) / 2
        return [
            (A.ρ_domain[1], v_ρ, :linear),
            (10^(sum(log10, A.ρ_domain) / 2), v_ρ, :linear),
            (sum(A.ρ_domain) / 2, v_ρ, :linear),
            (A.ρ_domain[2], v_ρ, :linear),
        ]
    end

    function ioditer()
        # Parameters
        params = Vector{NEOParameters{Float64}}(undef, 2)
        params[1] = NEOParameters(
            coeffstol = Inf, bwdoffset = 0.042, fwdoffset = 0.042, # Propagation
            adamiter = 500, adamQtol = 1e-5,                       # ADAM
            jtlsiter = 20, lsiter = 10, significance = 0.99,       # Least squares
            outrej = true, χ2_rec = 7.0, χ2_rej = 8.0,             # Outlier rejection
            fudge = 400.0, max_per = 20.0
        )
        params[2] = NEOParameters(params[1];
            coeffstol = 10.0, adamiter = 200, adamQtol = 0.01,
            lsiter = 5
        )
        # Naive initial conditions
        initcond = [initcond1, initcond2]
        # Initial orbit determination iterator
        return enumerate(Iterators.product(params, initcond))
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
        # Pre-allocate solutions
        sols = [zero(NEOSolution{Float64, Float64}) for _ in 1:4]
        # Start of computation
        init_time = now()
        # Initial orbit determination cycle
        for (i, j) in ioditer()
            # Unfold
            params, initcond = j
            # Initial orbit determination
            sols[i] = orbitdetermination(od, params; initcond)
            # Termination condition
            length(sols[i].res) == length(radec) &&
            critical_value(sols[i]) < params.significance && break
        end
        # Time of computation
        Δ = (now() - init_time).value
        # Select complete solutions
        mask = findall(s -> length(s.res) == length(radec), sols)
        isempty(mask) && return false
        # Choose best solution
        _, k = findmin(nrms, view(sols, mask))
        # Save orbit
        k = mask[k]
        jldsave(filename; sol = sols[k], Δ = Δ, k = k)
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