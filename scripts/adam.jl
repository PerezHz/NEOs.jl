using Distributed, ArgParse

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "adam.jl"
    # Desciption (for help screen)
    s.description = """
    Evaluate `NEOs.adam` over `N` equally spaced points
    in the boundary of a tracklet's admissible region."""

    s.epilog = """
        Example:\n
        \n
        # Evaluate 100 points in the boundary of 2014AA.txt with 10 workers\n
        # and 5 threads each\n
        julia -p 10 -t 5 --project scripts/adam.jl -i 2014AA.txt -N 100\n
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

@everywhere begin
    using NEOs, Dates, TaylorSeries, PlanetaryEphemeris, JLD2
    using NEOs: RadecMPC, AdmissibleRegion, PropagationBuffer, OpticalResidual,
        attr2bary, propres!, boundary_projection, reduce_tracklets, arboundary

    function adam(radec::Vector{RadecMPC{T}}, A::AdmissibleRegion{T}, ρ::T, v_ρ::T,
        params::NEOParameters{T}; scale::Symbol = :linear, η::T = 25.0,
        μ::T = 0.75, ν::T = 0.9, ϵ::T = 1e-8, Qtol::T = 0.001, adamorder::Int = 2,
        dynamics::D = newtonian!) where {T <: Real, D}
        # Initial time of integration [julian days]
        jd0 = datetime2julian(A.date)
        # Maximum number of iterations
        maxiter = params.adamiter
        # Allocate memory
        aes = Matrix{T}(undef, 6, maxiter+1)
        Qs = fill(T(Inf), maxiter+1)
        # Initial attributable elements
        aes[:, 1] .= [A.α, A.δ, A.v_α, A.v_δ, ρ, v_ρ]
        # Scaling factors
        scalings = Vector{T}(undef, 6)
        scalings[1:4] .= abs.(aes[1:4, 1]) ./ 1e6
        if scale == :linear
            scalings[5] = (A.ρ_domain[2] - A.ρ_domain[1]) / 1_000
        elseif scale == :log
            scalings[5] = (log10(A.ρ_domain[2]) - log10(A.ρ_domain[1])) / 1_000
        end
        scalings[6] = (A.v_ρ_domain[2] - A.v_ρ_domain[1]) / 1_000
        # Jet transport variables
        set_variables(Float64, "dx"; order = adamorder, numvars = 6)
        dae = [scalings[i] * TaylorN(i, order = adamorder) for i in 1:6]
        # Propagation buffer
        t0, tf = datetime2days(date(radec[1])), datetime2days(date(radec[end]))
        tlim = (t0 - params.bwdoffset, tf + params.fwdoffset)
        buffer = PropagationBuffer(dynamics, jd0, tlim, aes[:, 1] .+ dae, params)
        # Vector of O-C residuals
        res = Vector{OpticalResidual{T, TaylorN{T}}}(undef, length(radec))
        # Origin
        x0 = zeros(T, 6)
        x1 = zeros(T, 6)
        # Gradient of objective function wrt (ρ, v_ρ)
        g_t = Vector{T}(undef, 2)
        # First momentum
        m = zeros(T, 2)
        _m_ = zeros(T, 2)
        # Second momentum
        n = zeros(T, 2)
        _n_ = zeros(T, 2)
        # Gradient descent
        for t in 1:maxiter
            # Current attributable elements (plain)
            ae = aes[:, t]
            # Attributable elements (JT)
            if scale == :linear
                AE = ae + dae
            elseif scale == :log
                AE = [ae[1] + dae[1], ae[2] + dae[2], ae[3] + dae[3],
                    ae[4] + dae[4], 10^(log10(ae[5]) + dae[5]), ae[6] + dae[6]]
            end
            # Barycentric state vector
            q = attr2bary(A, AE, params)
            # Propagation and residuals
            # TO DO: `ρ::TaylorN` is too slow for `adam` due to evaluations
            # within the dynamical model
            propres!(res, radec, jd0 - ae[5]/c_au_per_day, q, params; buffer, dynamics)
            iszero(length(res)) && break
            # Least squares fit
            fit = tryls(res, x0, 5, 1:4)
            x1 .= fit.x
            # Current Q
            Q = nms(res)
            Q(x1) < 0 && break
            Qs[t] = Q(x1)
            # Convergence condition
            t > 1 && abs(Qs[t] - Qs[t-1]) / Qs[t] < Qtol && break
            # Gradient of objective function wrt (ρ, v_ρ)
            g_t[1] = differentiate(Q, 5)(x1)
            g_t[2] = differentiate(Q, 6)(x1)
            # First momentum
            m .= μ * m + (1 - μ) * g_t
            _m_ .= m / (1 - μ^t)
            # Second momentum
            n .= ν * n + (1 - ν) * g_t .^ 2
            _n_ .= n / (1 - ν^t)
            # Step
            x1[5:6] = x1[5:6] - η * _m_ ./ (sqrt.(_n_) .+ ϵ)
            # Update attributable elements
            aes[:, t+1] .= AE(x1)
            # Projection
            aes[5:6, t+1] .= boundary_projection(A, aes[5, t+1], aes[6, t+1])
        end

        return aes, Qs
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
    # Number of points
    N::Int = parsed_args["N"]
    # Maximum iterations per points
    maxiter::Int = parsed_args["maxiter"]
    # Print header
    println("`NEOS.adam` evaluation over the admissible region boundary")
    println("• Detected ", nworkers(), " worker(s) with ", Threads.nthreads(),
        " thread(s) each")
    println("• Input file: ", input)
    println("• Jet transport order: ", varorder)
    println("• Number of points: ", N)
    println("• Maximum iterations per point: ", maxiter)

    # Read optical astrometry
    radec = read_radec_mpc(input)
    # Orbit determination parameters
    params = NEOParameters(coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
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
    # Evaluate `NEOs.adam` in each point
    result = pmap(points_per_worker) do points
        aes = Vector{Matrix{Float64}}(undef, length(points))
        Qs = Vector{Vector{Float64}}(undef, length(points))
        for i in eachindex(points)
            x, v_ρ = points[i]
            aes[i], Qs[i] = adam(radec, A, 10^x, v_ρ, params; scale = :log,
                adamorder = varorder)
            j = findlast(!isinf, Qs[i])
            if isnothing(j)
                aes[i] = Matrix{Float64}(undef, 0, 0)
                Qs[i] = Vector{Float64}(undef, 0)
            else
                aes[i] = aes[i][:, 1:j]
                Qs[i] = Qs[i][1:j]
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
    output = string(radec[1].tmpdesig, "ADAM.jld2")
    jldsave(output; aes, Qs)
    println("• Output saved to: ", output)
    # Final time
    final_time = now()
    println("• Run started ", init_time, " and finished ", final_time)

    return nothing
end

main()
