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
        "--output", "-o"
            help = "output .jld2 file"
            arg_type = String
            default = "MCMOV.jld2"
        "--tracklet", "-t"
            help = "initial tracklet"
            arg_type = Int
            default = 1
        "--varorder", "-v"
            help = "jet transport order"
            arg_type = Int
            default = 5
        "--scale", "-s"
            help = "horizontal scale (log / linear)"
            arg_type = Symbol
            default = :log
        "--Qmax", "-q"
            help = "maximum allowed nms"
            arg_type = Float64
            default = 1e10
        "--xmin"
            help = "lower horizontal bound"
            arg_type = Float64
            default = -Inf
        "--xmax"
            help = "upper horizontal bound"
            arg_type = Float64
            default = Inf
        "--ymin"
            help = "lower vertical bound"
            arg_type = Float64
            default = -Inf
        "--ymax"
            help = "upper vertical bound"
            arg_type = Float64
            default = Inf
        "--Nx"
            help = "number of points in x"
            arg_type = Int
            default = 100
        "--Ny"
            help = "number of points in y"
            arg_type = Int
            default = 100
    end

    s.epilog = """
        If the bounds `[xmin, xmax]×[ymin, ymax]` are oimitted, the script will fit
        a box to the AR. The horizontal bounds `[xmin, xmax]` must be consistent with
        `scale`.
    """

    return parse_args(s)
end

@everywhere begin
    using NEOs, Dates, TaylorSeries, PlanetaryEphemeris, JLD2
    using NEOs: AdmissibleRegion, OpticalResidual, RadecMPC, PropagationBuffer,
                reduce_tracklets, argoldensearch, scaled_variables, attr2bary,
                propres!, nobs, _lsmethods

    function mcmov(points::Vector{Tuple{T, T}}, A::AdmissibleRegion{T},
        radec::Vector{RadecMPC{T}}, bounds::Vector{Float64},
        params::NEOParameters{T}, varorder::Int = 2,
        scale::Symbol = :log, Qmax::Float64 = 1e10) where {T <: Real}
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
        xmin, xmax, ymin, ymax = bounds
        scalings[5:6] .= (xmax - xmin) / 100, (ymax - ymin) / 100
        # Attributable elements (jet transport)
        AE = ae .+ scalings .* dae
        # TDB epoch of admissible region
        _jd0_ = dtutc2jdtdb(A.date)
        # Propagation buffer
        _buffer_ = PropagationBuffer(od, _jd0_, 1, nobs(od), AE(), params)
        buffer = PropagationBuffer(od, _jd0_, 1, nobs(od), AE, params)
        # Vector of residuals
        _res_ = [zero(OpticalResidual{T, T}) for _ in eachindex(radec)]
        res = [zero(OpticalResidual{T, TaylorN{T}}) for _ in eachindex(radec)]
        # Origin
        x0 = zeros(T, 6)
        # Least squares cache and methods
        lscache, lsmethods = LeastSquaresCache(x0, 1:4, 20), _lsmethods(res, x0, 1:4)
        # Manifold of variations
        mov = fill(T(Inf), 13, length(points))
        # Iterate mov points
        for (i, point) in enumerate(points)
            # Attributable elements (plain)
            ae[5:6] .= point
            # Attributable elements (JT)
            AE .= ae .+ scalings .* dae
            # Topocentric range
            if scale == :log
                AE[5] = 10^AE[5]
            end
            # Barycentric initial conditions (JT)
            q = attr2bary(A, AE, params)
            # Reference epoch in julian days (corrrected for light-time)
            jd0 = _jd0_ - constant_term(AE[5]) / c_au_per_day
            # Propagation and residuals
            propres!(_res_, od, jd0, q(), params; buffer = _buffer_)
            if isempty(_res_)
                _res_ = [zero(OpticalResidual{T, T}) for _ in eachindex(radec)]
                continue
            end
            nms(_res_) > Qmax && continue
            propres!(res, od, jd0, q, params; buffer)
            # Least squares fit
            fit = tryls(res, x0, lscache, lsmethods)
            !fit.success && continue
            # Objective function
            Q = nms(res)
            # Save results
            mov[1:6, i] .= fit.x
            mov[7:10, i] .= ae[1:4] .+ scalings[1:4] .* fit.x[1:4]
            mov[11:12, i] .= point
            mov[13, i] = Q(fit.x)
        end
        # Select successful points
        mask = @. !isinf(mov[13, :])

        return mov[:, mask]
    end
end

function main()
    # Initial time
    init_time = now()
    # Parse arguments from commandline
    parsed_args = parse_commandline()
    # Input and output files
    input::String, output::String = parsed_args["input"], parsed_args["output"]
    # Initial tracklet
    idxtrk::Int = parsed_args["tracklet"]
    # Jet transport order
    varorder::Int = parsed_args["varorder"]
    # Horizontal scale
    scale::Symbol = Symbol(parsed_args["scale"])
    # Maximum allowed nms
    Qmax::Float64 = parsed_args["Qmax"]
    # Horizontal and vertical bounds
    xmin::Float64, xmax::Float64 = parsed_args["xmin"], parsed_args["xmax"]
    ymin::Float64, ymax::Float64 = parsed_args["ymin"], parsed_args["ymax"]
    # Number of points in x (y)
    Nx::Int, Ny::Int = parsed_args["Nx"], parsed_args["Ny"]
    # Number of workers
    Nw = nworkers()
    # Print header
    println("Manifold of variations computation via Jet Transport Monte Carlo")
    println("• Detected ", Nw, " worker(s) with ", Threads.nthreads(), " thread(s) each")
    println("• Input file: ", input)

    # Read optical astrometry
    radec = read_radec_mpc(input)
    # Orbit determination parameters
    params = NEOParameters(
        coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007, # Propagation
        jtlsiter = 20, lsiter = 20, significance = 0.99,       # Least squares
    )
    # Reduce tracklets
    tracklets = reduce_tracklets(radec)
    idxtrk = clamp(idxtrk, 1, length(tracklets))
    println("• Initial tracklet: ", idxtrk)
    # Admissible region
    A = AdmissibleRegion(tracklets[idxtrk], params)
    # Set jet transport variables
    set_variables(Float64, "dx"; order = varorder, numvars = 6)
    println("• Jet transport order: ", varorder)

    # Global box
    println("• Horizontal scale: ", scale)
    println("• Maximum allowed nms: ", Qmax)
    ρmin, ρmax = A.ρ_domain
    if scale == :log
        xmin, xmax = max(xmin, log10(ρmin)), min(xmax, log10(ρmax))
    elseif scale == :linear
        xmin, xmax = max(xmin, ρmin), min(xmax, ρmax)
    end
    ymin = max(ymin, argoldensearch(A, ρmin, ρmax, :min, :outer, 1e-20)[2])
    ymax = min(ymax, argoldensearch(A, ρmin, ρmax, :max, :outer, 1e-20)[2])
    bounds = [xmin, xmax, ymin, ymax]
    # Grid of domain points
    side_x = LinRange(xmin, xmax, Nx)
    side_y = LinRange(ymin, ymax, Ny)
    points = Iterators.product(side_x, side_y)
    println("• Number of points in x (y): ", Nx, " (", Ny, ")")
    # Distribute domain points over workers
    points_per_worker = [Tuple{Float64, Float64}[] for _ in 1:Nw]
    k = 1
    for point in points
        # Check point is inside AR
        mask = scale == :log ? ((10^point[1], point[2]) in A) :
            ((point[1], point[2]) in A)
        !mask && continue
        # Save point in corresponding worker
        push!(points_per_worker[k], point)
        k = mod1(k+1, Nw)
    end
    Np = sum(length, points_per_worker)
    println("• ", Np, " points in the manifold of variations")
    # Manifold of variations
    movs = pmap(points_per_worker) do p
        return mcmov(p, A, radec, bounds, params, varorder, scale, Qmax)
    end
    mov = reduce(hcat, movs)
    println("• ", count(!isnan, view(mov, 13, :)), "/", Np, " points converged")

    # Save results
    jldsave(output; mov = mov)
    println("• Output saved to: ", output)

    # Final time
    final_time = now()
    println("• Run started ", init_time, " and finished ", final_time)

    return nothing
end

main()