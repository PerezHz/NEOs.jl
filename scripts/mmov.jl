using Distributed, ArgParse
using Dates, TaylorSeries, PlanetaryEphemeris, JLD2
@everywhere using NEOs
using NEOs: AdmissibleRegion, reduce_tracklets, arboundary

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "mmov.jl"
    # Desciption (for help screen)
    s.description = """
    Evaluate `NEOs.mmov` over `N` equally spaced points
    in the boundary of a tracklet's admissible region."""

    s.epilog = """
        Example:\n
        \n
        # Evaluate 100 points in the boundary of 2014AA.txt with 10 workers\n
        # and 5 threads each\n
        julia -p 10 -t 5 --project scripts/mmov.jl -i 2014AA.txt -N 100\n
        \n
    """

    @add_arg_table! s begin
        "--input", "-i"
            help = "input optical astrometry file"
            arg_type = String
        "--varorder", "-v"
            help = "jet transport order"
            arg_type = Int
            default = 2
        "--N", "-N"
            help = "number of points"
            arg_type = Int
            default = 100
        "--maxiter", "-m"
            help = "maximum iterations per point"
            arg_type = Int
            default = 200
    end

    return parse_args(s)
end

function main()
    # Initial time
    init_time = now()
    # Parse arguments from commandline
    parsed_args = parse_commandline()
    # Input file
    input::String = parsed_args["input"]
    # Jet transport order
    varorder::Int = parsed_args["varorder"]
    # Number of points
    N::Int = parsed_args["N"]
    # Maximum iterations per points
    maxiter::Int = parsed_args["maxiter"]
    # Print header
    println("`NEOs.mmov` evaluation over the admissible region boundary")
    println("• Detected ", nworkers(), " worker(s) with ", Threads.nthreads(),
        " thread(s) each")
    println("• Input file: ", input)
    println("• Jet transport order: ", varorder)
    println("• Number of points: ", N)
    println("• Maximum iterations per point: ", maxiter)

    # Read optical astrometry
    radec = read_radec_mpc(input)
    # Orbit determination parameters
    params = Parameters(coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
        adamiter = maxiter)
    # Reduce tracklet
    tracklet = reduce_tracklets(radec)[1]
    # Admissible region
    A = AdmissibleRegion(tracklet, params)
    # Set jet transport variables
    set_variables(Float64, "dx"; order = varorder, numvars = 6)

    # Generate N points over the external boundary of A
    points = map(t -> arboundary(A, t, :outer, :log), LinRange(0, 3, N))
    # Distribute points over workers
    Np = round(Int, length(points) / nworkers())
    idxs = Iterators.partition(eachindex(points), Np)
    points_per_worker = [points[i] for i in idxs]
    # Evaluate `NEOs.mmov` in each point
    result = pmap(points_per_worker) do points
        # Orbit determination problem
        od = ODProblem(newtonian!, radec)
        # Pre-allocate output
        aes = Vector{Matrix{Float64}}(undef, length(points))
        Qs = Vector{Vector{Float64}}(undef, length(points))
        # Main loop
        for i in eachindex(points)
            x, v_ρ = points[i]
            porbit = mmov(od, 1, A, 10^x, v_ρ, params; scale = :log, adamorder = varorder)
            if iszero(porbit)
                aes[i] = Matrix{Float64}(undef, 0, 0)
                Qs[i] = Vector{Float64}(undef, 0)
            else
                aes[i] = porbit.aes
                Qs[i] = porbit.Qs
            end
        end
        filter!(!isempty, aes)
        filter!(!isempty, Qs)
        return aes, Qs
    end
    # Separate attributable elements from Qs
    aes = reduce(vcat, first.(result))
    Qs = reduce(vcat, last.(result))

    # Save results
    output = string(radec[1].tmpdesig, "MMOV.jld2")
    jldsave(output; A, aes, Qs)
    println("• Output saved to: ", output)
    # Final time
    final_time = now()
    println("• Run started ", init_time, " and finished ", final_time)

    return nothing
end

main()
