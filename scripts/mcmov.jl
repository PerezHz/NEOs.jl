using Distributed, ArgParse

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "mcmov.jl"
    # Desciption (for help screen)
    s.description = """
    Compute the manifold of variations of an optical astrometry
    file via Jet Transport Monte Carlo."""

    s.epilog = """
        Example:\n
        \n
        # Compute the MOV of 2014AA.txt with 10 workers and 5 threads each\n
        julia -p 10 -t 5 --project scripts/mcmov.jl -i 2014AA.txt\n
        \n
    """

    @add_arg_table! s begin
        "--input", "-i"
            help = "input optical astrometry file"
            arg_type = String
        "--varorder", "-v"
            help = "jet transport order"
            arg_type = Int
            default = 5
        "--Nx", "-x"
            help = "number of points in x"
            arg_type = Int
            default = 100
        "--Ny", "-y"
            help = "number of points in y"
            arg_type = Int
            default = 100
    end

    return parse_args(s)
end

@everywhere begin
    using NEOs, Dates, TaylorSeries, PlanetaryEphemeris, JLD2
    using NEOs: AdmissibleRegion, OpticalResidual, RadecMPC, PropagationBuffer,
                reduce_tracklets, argoldensearch, scaled_variables, attr2bary,
                propres!, nobs

    function mcmov(points::Vector{Tuple{T, T}}, A::AdmissibleRegion{T},
        radec::Vector{RadecMPC{T}}, params::NEOParameters{T},
        varorder::Int = 2) where {T <: Real}
        # Orbit determination problem
        od = ODProblem(newtonian!, radec)
        # JT variables
        dae = scaled_variables("dx", ones(T, 6), order = varorder)
        # Attributable elements (plain)
        ae = Vector{T}(undef, 6)
        ae[1:4] .= A.α, A.δ, A.v_α, A.v_δ
        ae[5:6] .= points[1]
        # Scaling factors
        scalings = Vector{T}(undef, 6)
        scalings[1:4] .= abs.(ae[1:4]) ./ 1e6
        x_min, x_max = log10.(A.ρ_domain)
        _, y_min = argoldensearch(A, 10^x_min, 10^x_max, :min, :outer, 1e-20)
        _, y_max = argoldensearch(A, 10^x_min, 10^x_max, :max, :outer, 1e-20)
        scalings[5:6] .= (x_max - x_min) / 100, (y_max - y_min) / 100
        # Attributable elements (jet transport)
        AE = ae .+ scalings .* dae
        # TDB epoch of admissible region
        _jd0_ = dtutc2jdtdb(A.date)
        # Propagation buffer
        buffer = PropagationBuffer(od, _jd0_, 1, nobs(od), AE, params)
        # Vector of residuals
        res = [zero(OpticalResidual{T, TaylorN{T}}) for _ in eachindex(radec)]
        # Origin
        x0 = zeros(T, 6)
        # Manifold of variations
        mov = Matrix{T}(undef, 13, length(points))
        # Iterate mov points
        for (i, point) in enumerate(points)
            # Attributable elements (plain)
            ae[5:6] .= point
            # Attributable elements (JT)
            AE .= ae .+ scalings .* dae
            # Topocentric range
            AE[5] = 10^AE[5]
            # Barycentric initial conditions (JT)
            q = attr2bary(A, AE, params)
            # Reference epoch in julian days (corrrected for light-time)
            jd0 = _jd0_ - constant_term(AE[5]) / c_au_per_day
            # Propagation and residuals
            propres!(res, od, jd0, q, params; buffer)
            # Least squares fit
            fit = newtonls(res, x0, 10, 1:4)
            # Objective function
            Q = nms(res)
            # Save results
            mov[1:6, i] .= fit.x
            mov[7:10, i] .= ae[1:4] .+ scalings[1:4] .* fit.x[1:4]
            mov[11:12, i] .= point
            mov[13, i] = Q(fit.x)
        end

        return mov
    end
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
    # Number of points in x (y)
    Nx::Int = parsed_args["Nx"]
    Ny::Int = parsed_args["Ny"]
    # Number of workers
    Nw = nworkers()
    # Print header
    println("Manifold of variations computation via Jet Transport Monte Carlo")
    println("• Detected ", Nw, " worker(s) with ", Threads.nthreads(), " thread(s) each")
    println("• Input file: ", input)
    println("• Jet transport order: ", varorder)
    println("• Number of points in x (y): ", Nx, " (", Ny, ")")

    # Read optical astrometry
    radec = read_radec_mpc(input)
    # Orbit determination parameters
    params = NEOParameters(coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007)
    # Reduce tracklet
    tracklet = reduce_tracklets(radec)[1]
    # Admissible region
    A = AdmissibleRegion(tracklet, params)
    # Set jet transport variables
    set_variables(Float64, "dx"; order = varorder, numvars = 6)

    # Global box
    x_min, x_max = log10.(A.ρ_domain)
    _, y_min = argoldensearch(A, A.ρ_domain[1], A.ρ_domain[2], :min, :outer, 1e-20)
    _, y_max = argoldensearch(A, A.ρ_domain[1], A.ρ_domain[2], :max, :outer, 1e-20)
    # Grid of domain points
    side_x = LinRange(x_min, x_max, Nx)
    side_y = LinRange(y_min, y_max, Ny)
    points = Iterators.product(side_x, side_y)
    # Distribute domain points over workers
    points_per_worker = [Tuple{Float64, Float64}[] for _ in 1:Nw]
    k = 1
    for point in points
        # Check point is inside AR
        if !((10^point[1], point[2]) in A)
             continue
        end
        # Save point in corresponding worker
        push!(points_per_worker[k], point)
        k = mod1(k+1, Nw)
    end
    println("• ", sum(length, points_per_worker), " domain points in the manifold of variations")
    # Manifold of variations
    movs = pmap(p -> mcmov(p, A, radec, params, varorder), points_per_worker)
    mov = reduce(hcat, movs)

    # Save results
    tmpdesig = radec[1].tmpdesig
    output = string(tmpdesig, "MCMOV.jld2")
    jldsave(output; mov = mov)
    println("• Output saved to: ", output)

    # Final time
    final_time = now()
    println("• Run started ", init_time, " and finished ", final_time)

    return nothing
end

main()