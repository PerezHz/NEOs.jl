using Distributed, ArgParse

function parse_commandline()
    s = ArgParseSettings(add_version = true, version = "0.2")

    # Program name (for usage & help screen)
    s.prog = "mcmov.jl"
    # Desciption (for help screen)
    s.description = "Sample the manifold of variations (MOV) of a NEOCP \
    object using jet transport-assisted monte carlo."

    s.epilog = """
        Example:\n
        \n
        # Sample the MOV of P22hRXJ with 10 workers and 5 threads each\n
        julia -p 10 -t 5 --project mcmov.jl -i P22hRXJ -t 0000000Hz3SX -s linear --scout\n
        \n
    """

    @add_arg_table! s begin
        "--input", "-i"
            help = "NEOCP designation (trksub)"
            arg_type = String
        "--output", "-o"
            help = "Output file"
            arg_type = String
        "--trkids", "-t"
            help = "Trkids to include in OD"
            nargs = '+'
        "--scale", "-s"
            help = "Horizontal scale (log / linear)"
            arg_type = String
            default = "log"
        "--Nx"
            help = "Number of points in x"
            arg_type = Int
            default = 100
        "--Ny"
            help = "Number of points in y"
            arg_type = Int
            default = 100
        "--maxchi"
            help = "χ value threshold"
            arg_type = Float64
            default = 5.0
        "--nominal"
            help = "Compute a nominal orbit"
            arg_type = Bool
            default = true
        "--refine"
            help = "Refine the first grid"
            action = :store_true
        "--scout"
            help = "Fetch JPL Scout data and save it into a .csv file"
            action = :store_true
        "--neoscan"
            help = "Fetch NEODyS NEOScan data and save it into a .mov_sample file"
            action = :store_true
        "--neocp"
            help = "Fetch MPC NEOCP data and save it into a .orb file"
            action = :store_true
    end

    return parse_args(s)
end

@everywhere using Printf

@everywhere begin
    using NEOs, Dates, TaylorSeries, PlanetaryEphemeris, JLD2, HTTP, Statistics,
          ChunkSplitters, StaticArraysCore, DataFrames, CSV, JSON3
    using NEOs: PropresBuffer, PropagationBuffer, OpticalBuffer, OpticalADES,
                AbstractOrbit, KeplerianElements, parse_optical_rwo, argoldensearch,
                evaldeltas, init_optical_residuals, indices, _lsmethods, μ_S,
                equatorial2ecliptic, _propagate
    import NEOs: keplerian, initialcondition

    const VariantOrbit{T} = MMOVOrbit{typeof(newtonian!), T, T, Vector{OpticalADES{T}}}

    const NEOCP_ORBITS_HEADER = "Object   H     G    Epoch    M         Peri.      \
          Node       Incl.        e           n         a                     NObs \
          NOpp   Arc    r.m.s.       Orbit ID"

    computationtime(x::DateTime, y::DateTime) = (y - x).value / 60_000

    initialcondition(x::AbstractOrbit) = x(), epoch(x) + PE.J2000

    generate_grid(B::AbstractVector, Nx::Int, Ny::Int) =
        Iterators.product(LinRange(B[1], B[2], Nx), LinRange(B[3], B[4], Ny))

    function printitle(s::AbstractString, d::AbstractString)
        l = repeat(d, length(s))
        println(l, '\n', s, '\n', l)
    end

    function chi(x::AbstractVector{VariantOrbit{T}}) where {T <: Real}
        Qmin, i = findmin(nms, x)
        nobs = 2 * noptical(x[i])
        χs = @. sqrt(nobs * ( nms(x) - Qmin ))
        return χs
    end

    function fetch_scout_orbits(input::AbstractString)
        uri = HTTP.URI(
            scheme = "https",
            host   = "ssd-api.jpl.nasa.gov",
            path   = "/scout.api",
            query  = "tdes=$(input)&orbits=1"
        )
        response_scout = HTTP.get(string(uri) #=, require_ssl_verification = false=#)
        orbits_data = JSON3.read(response_scout.body)["orbits"]["data"]
        rows = [collect(row) for row in orbits_data]
        orbits_fields = JSON3.read(response_scout.body)["orbits"]["fields"]
        mat = vcat([permutedims(r) for r in rows]...)  # or hcat(rows...)' as another option
        df = DataFrame(mat, orbits_fields)
        CSV.write("$(input).csv", df)
        println("• Fetched Scout data; saved to: $(input).csv")
        return nothing
    end

    function fetch_neoscan_orbits(input::AbstractString)
        uri = HTTP.URI(
            scheme = "https",
            host   = "newton.spacedys.com",
            path   = "/neodys/NEOScan/scan_neocp/$input/$input.mov_sample"
        )
        response_neoscan = HTTP.get(string(uri) #=, require_ssl_verification = false=#)
        write("$input.mov_sample", response_neoscan.body)
        println("• Fetched NEOScan data; saved to: $(input).mov_sample")
        return nothing
    end

    function fetch_neocp_orbits(input::AbstractString)
        url_neocp = "https://cgi.minorplanetcenter.net/cgi-bin/showobsorbs.cgi"
        data = Dict("Obj" => input, "orb" => "y")
        response_neocp = HTTP.post(url_neocp, [], data #=, require_ssl_verification = false=#)
        text = String(response_neocp.body)
        lines = split(text, '\n')[2:end-2]
        write("$input.orb", join(lines, '\n'))
        println("• Fetched NEOCP data; saved to: $(input).orb")
        return nothing
    end

    function fetch_neodys_weights(desig::AbstractString)
        url = "https://newton.spacedys.com/neodys/NEOScan/scan_neocp/$desig/$desig.rwo"
        resp = HTTP.get(url)
        text = String(resp.body)
        optical = parse_optical_rwo(text)
        σs = rms.(optical)
        return @. tuple(1 / first(σs), 1 / last(σs))
    end

    function global_box(A::AdmissibleRegion, scale::Symbol)
        ρmin, ρmax = A.ρ_domain
        if scale == :log
            xmin, xmax = log10(ρmin), log10(ρmax)
        elseif scale == :linear
            xmin, xmax = ρmin, ρmax
        end
        ymin = argoldensearch(A, ρmin, ρmax, :min, :outer, 1e-20)[2]
        ymin = min(ymin, A.v_ρ_domain[1])
        ymax = argoldensearch(A, ρmin, ρmax, :max, :outer, 1e-20)[2]
        ymax = max(ymax, A.v_ρ_domain[2])
        bounds = [xmin, xmax, ymin, ymax]
        return bounds
    end

    function refined_box(mask::AbstractVector{Bool}, scale::Symbol,
                         points::Vector{Vector{NTuple{2, T}}}) where {T <: Real}
        xmin, xmax = typemax(T), typemin(T)
        ymin, ymax = typemax(T), typemin(T)
        for (i, point) in enumerate(Iterators.flatten(points))
            if mask[i]
                ρ, v_ρ = point
                xmin, xmax = min(xmin, ρ), max(xmax, ρ)
                ymin, ymax = min(ymin, v_ρ), max(ymax, v_ρ)
            end
        end
        if scale == :log
            xmin, xmax = log10(xmin), log10(xmax)
        end
        return [xmin, xmax, ymin, ymax]
    end

    function distribute_points(A::AdmissibleRegion, Nworkers::Int, scale::Symbol, points)
        points_per_worker = [NTuple{2, Float64}[] for _ in 1:Nworkers]
        k = 1
        for point in points
            # Check if point is inside the admissible region
            if scale == :log && (10^point[1], point[2]) in A
                push!(points_per_worker[k], (10^point[1], point[2]))
                k = mod1(k+1, Nworkers)
            elseif scale == :linear && (point[1], point[2]) in A
                push!(points_per_worker[k], (point[1], point[2]))
                k = mod1(k+1, Nworkers)
            end
        end
        return points_per_worker
    end

    function mcmov(points::AbstractVector{NTuple{2, T}}) where {T <: Real}
        # Attributable elements (plain)
        ae = Vector{T}(undef, 6)
        ae[1:4] .= A.ra, A.dec, A.vra, A.vdec
        ae[5:6] .= points[1]
        # Scaling factors
        scalings = Vector{T}(undef, 6)
        @. scalings[1:4] = abs(ae[1:4]) / 1E6
        xmin, xmax, ymin, ymax = bounds
        scalings[5:6] .= (xmax - xmin) / 1E3, (ymax - ymin) / 1E3
        # Jet transport variables
        dae = scaled_variables("dx", scalings, order = 2)
        variables = collect(1:6)
        # Attributable elements (jet transport)
        AE = Vector{TaylorN{T}}(undef, 6)
        @. AE[1:4] = ae[1:4] + dae[1:4]
        if scale == :linear
            AE[5] = ae[5] + dae[5]
        elseif scale == :log
            AE[5] = 10^(log10(ae[5]) + dae[5])
        end
        AE[6] = ae[6] + dae[6]
        # TDB epoch of admissible region
        jd0 = dtutc2jdtdb(A.date)
        # Initialize buffer and set of residuals
        idxs = indices(od.tracklets)
        buffer = PropresBuffer(od, AE, jd0, idxs, params)
        res = init_optical_residuals(TaylorN{T}, od, idxs)
        # Number of observations
        nobs = notoutobs(res)
        # Origin
        x0 = zeros(T, 6)
        # Least squares cache and methods
        lscache = LeastSquaresCache(x0, 1:4, 20)
        lsmethods = _lsmethods(res, x0, 1:4)
        # Manifold of variations
        orbits = [zero(VariantOrbit{T}) for _ in eachindex(points)]
        # Iterate mov points
        for (i, point) in enumerate(points)
            # Attributable elements (plain)
            ae[5:6] .= point
            # Attributable elements (JT)
            @. AE[1:4] = ae[1:4] + dae[1:4]
            if scale == :linear
                AE[5] = ae[5] + dae[5]
            elseif scale == :log
                AE[5] = 10^(log10(ae[5]) + dae[5])
            end
            AE[6] = ae[6] + dae[6]
            # Barycentric initial conditions (JT)
            q = attr2bary(A, AE, params)
            # Propagation and residuals
            bwd, fwd = propres!(res, od, q, jd0 - ae[5] / c_au_per_day, params; buffer, idxs)
            if isempty(res)
                res = init_optical_residuals(TaylorN{T}, od, idxs)
                continue
            end
            # Least squares fit
            fit = tryls(res, x0, lscache, lsmethods)
            !fit.success && continue
            # Current Q
            Q = nms(res)
            Q(fit.x) < 0 && continue
            # Covariance matrix
            C = (nobs/2) * TS.hessian(Q, fit.x)
            covariance = inv(C)
            # Residuals space to barycentric coordinates jacobian
            jacobian = Matrix(TS.jacobian(q - constant_term(q), fit.x))
            # Update orbit
            orbits[i] = evaldeltas(MMOVOrbit(
                newtonian!, variables, od.optical, od.tracklets, bwd, fwd,
                res, covariance, jacobian, [AE(fit.x);;], [Q(fit.x)]
            ), fit.x)
        end

        return orbits
    end

    function radec_next_day(orbits::AbstractVector{VariantOrbit{T}}) where {T <: Real}
        radec = Vector{NTuple{2, T}}(undef, length(orbits))
        q0, jd0 = initialcondition(orbits[1])
        t0 = minimum(epoch, orbits) - params.bwdoffset
        tf = dtutc2days(DAY_AFTER_EPOCH) + params.fwdoffset
        pbuffer = PropagationBuffer(newtonian!, q0, jd0, (t0, tf), params);
        for (i, orbit) in enumerate(orbits)
            q0, jd0 = initialcondition(orbit)
            tmax = ( tf + PE.J2000 - jd0 ) / yr
            fwd = _propagate(newtonian!, q0, jd0, tmax, pbuffer, params)
            obuffer = OpticalBuffer(zero(T))
            radec[i] = compute_radec(OBSERVER, DAY_AFTER_EPOCH, obuffer;
                xvs = params.eph_su, xve = params.eph_ea,
                xva = (orbit.bwd, fwd)
            )
        end
        return radec
    end

    function keplerian(orbit::VariantOrbit{T}, t::T,
                       params::Parameters{T}) where {T <: Real}
        # Reference epoch [MJD TDB]
        mjd0 = t + MJD2000
        # Scalar initial condition
        q0 = equatorial2ecliptic(orbit(t) - params.eph_su(t))
        # Osculating orbital elements
        elements = cartesian2keplerian(q0, mjd0; μ = μ_S)
        Γ_kep = SMatrix{6, 6}(fill(NaN, 6, 6))
        kep = KeplerianElements{T, T}(μ_S, mjd0, :ecliptic, elements, Γ_kep)

        return kep
    end

    function write_neocp_orbits(norbits::AbstractVector{<:AbstractOrbit},
                                vorbits::AbstractVector{<:AbstractOrbit},
                                filename::AbstractString)
        open(filename, "w") do file
            # Header
            write(file, NEOCP_ORBITS_HEADER, '\n')
            # Orbits
            for (j, orbit) in enumerate(Iterators.flatten((norbits, vorbits)))
                # Absolute magnitude
                H, _ = absolutemagnitude(orbit, params)
                # Slope parameter
                G = params.slope
                # Orbital elements
                kep = keplerian(orbit, REFERENCE_EPOCH, params)
                M = mod(meananomaly(kep), 360)
                ω = mod(argperi(kep), 360)
                Ω = mod(longascnode(kep), 360)
                i = mod(inclination(kep), 180)
                e = eccentricity(kep)
                n = meanmotion(kep)
                a = semimajoraxis(kep)
                # Number of observations
                nobs = noptical(orbit)
                # Arc length [days]
                arc = floor(Int, numberofdays(orbit.optical))
                # RMS
                Q = nrms(orbit)
                # Orbit ID
                id = isone(j) ? "NEOCPNomin" : "NEOCPV" * lpad(j-1, 4, '0')
                # Assemble line
                s = string(
                    rpad(input, 8),
                    rpad(@sprintf("%.1f", H), 6),
                    # ' ' ^ 6,
                    rpad(@sprintf("%.2f", G), 6),
                    ' ' ^ 6, # Epoch
                    rpad(@sprintf("%9.5f", M), 11),
                    rpad(@sprintf("%9.5f", ω), 11),
                    rpad(@sprintf("%9.5f", Ω), 11),
                    rpad(@sprintf("%9.5f", i), 11),
                    rpad(@sprintf("%9.7f", e), 11),
                    rpad(@sprintf("%10.8f", n), 11),
                    rpad(@sprintf("%11.7f", a), 11),
                    lpad(nobs, 19),
                    "   1 ", # Number of oppsitions
                    lpad(string(arc, " days "), 10),
                    rpad(@sprintf("%.2f", Q), 13),
                    id,
                    '\n'
                )
                write(file, s)
            end
        end
        return nothing
    end
end

function main()
    # Parse arguments from commandline
    parsed_args = parse_commandline()

    # Print header
    printitle("Manifold of variations sampling for NEOCP objects", "=")
    printitle("Parameters", "-")

    # Number of workers and threads
    Nworkers, Nthreads = nworkers(), Threads.nthreads()
    println("• Detected ", Nworkers, " worker(s) with ", Nthreads, " thread(s) each")

    # Input asteroid desgination
    input::String = parsed_args["input"]
    println("• Input NEOCP designation (trksub): ", input)

    # Output file
    output::String = isnothing(parsed_args["output"]) ? input*".neosjl" : parsed_args["output"]
    println("• Output file: ", output)

    # Horizontal scale
    scale_str::String = parsed_args["scale"]
    scale::Symbol = Symbol(scale_str)
    println("• Horizontal scale: ", scale)

    # Number of points in x (y)
    Nx::Int, Ny::Int = parsed_args["Nx"], parsed_args["Ny"]
    println("• Number of points in x (y): ", Nx, " (", Ny, ")")

    # χ value threshold
    χ_max::Float64 = parsed_args["maxchi"]
    println("• χ threshold: ", χ_max)

    # Compute a nominal orbit?
    compute_nominal = parsed_args["nominal"]
    println("• Compute nominal orbit?: ", compute_nominal)

    # Refine the first grid?
    refine_grid = parsed_args["refine"]
    println("• Refine the first grid?: ", refine_grid)

    # Fetching JPl Scout data?
    fetch_scout = parsed_args["scout"]
    println("• Fetch Scout data?: ", fetch_scout)

    # Fetch NEODyS NEOScan data?
    fetch_neoscan = parsed_args["neoscan"]
    println("• Fetch NEOScan data?: ", fetch_neoscan)

    # Fetch MPC NEOCP data?
    fetch_neocp = parsed_args["neocp"]
    println("• Fetch NEOCP data?: ", fetch_neocp)

    # Initial time
    initial_time = now()
    printitle("Computation", "-")
    println("• Run started at ", initial_time)

    # Get JPL Scout data for same object, if requested by user
    fetch_scout && fetch_scout_orbits(input)
    # Get NEODyS NEOScan data for same object, if requested by user
    fetch_neoscan && fetch_neoscan_orbits(input)
    # Get MPC NEOCP data for same object, if requested by user
    fetch_neocp && fetch_neocp_orbits(input)

    # Fetch optical astrometry
    optical_all = fetch_optical_ades(input, NEOCP)
    trkids::Vector{String} = parsed_args["trkids"]
    # If `trkids` is empty use all the astrometry; else, use only trkids contained in `trkids`
    optical = isempty(trkids) ? optical_all : filter(x->x.trkid in trkids, optical_all)
    println("• `trkids` included in run: ", unique(map(x->x.trkid, optical_all)))
    # Orbit determination problem
    od = ODProblem(newtonian!, optical)
    od.weights.weights .= fetch_neodys_weights(input)
    # Parameters
    params = Parameters(
        maxsteps = 1_000, order = 15, abstol = 1E-12, parse_eqs = true,
        coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
        jtlsorder = 2, jtlsmask = false, jtlsiter = 20, lsiter = 20,
        significance = 0.99, jtlsproject = true, outrej = false,
    )
    # Admissible region
    A = AdmissibleRegion(od.tracklets[1], params)
    # Backward offset must take into consideration the -ρ/c relativistic
    # correction to the epoch
    bwdoffset = params.bwdoffset + A.ρ_domain[2] / c_au_per_day
    params = Parameters(params; bwdoffset)

    # Global box
    bounds = global_box(A, scale)
    println("• Global box: ", bounds)

    # Declare constants in all workers
    @everywhere begin
        const input = $input
        const scale = Symbol($scale_str)
        const od = $od
        const params = $params
        const A = $A
        const bounds = $bounds
        const TRACKLET = first(od.tracklets)
        const REFERENCE_EPOCH = dtutc2days(TRACKLET)
        const DAY_AFTER_EPOCH = date(TRACKLET) + Day(1)
        const OBSERVER = observatory(TRACKLET)
    end

    # Grid of domain points
    points = generate_grid(bounds, Nx, Ny)
    # Distribute domain points over workers
    points_per_worker = distribute_points(A, Nworkers, scale, points)
    Npoints = sum(length, points_per_worker)
    println("• ", Npoints, " points in the manifold of variations")

    # Manifold of variations
    orbits = reduce(vcat, pmap(mcmov, points_per_worker))
    # Eliminate orbits with χ > 5
    χs = chi(orbits)
    mask = χs .≤ χ_max
    keepat!(orbits, mask)
    println("• ", length(orbits), " / ", Npoints, " points with χ ≤ χ_max = $(χ_max)")

    if refine_grid
        # Global box
        bounds .= refined_box(mask, scale, points_per_worker)
        @everywhere bounds .= $bounds
        println("• Refined box: ", bounds)
        # Grid of domain points
        points = generate_grid(bounds, Nx, Ny)
        # Distribute domain points over workers
        points_per_worker = distribute_points(A, Nworkers, scale, points)
        Npoints = sum(length, points_per_worker)
        println("• ", Npoints, " points in the manifold of variations")
        # Manifold of variations
        orbits = reduce(vcat, pmap(mcmov, points_per_worker))
        # Eliminate orbits with χ > 5
        χs = chi(orbits)
        keepat!(orbits, χs .≤ χ_max)
        println("• ", length(orbits), " / ", Npoints, " points with χ ≤ χ_max = $(χ_max)")
    end

    if compute_nominal
        # Right ascension and declination one day after REFERENCE_EPOCH
        radec = reduce(vcat, pmap(radec_next_day, chunks(orbits, n = Nworkers)))
        αs, δs = @. first(radec), last(radec)
        # Find closest orbit to the median
        αmedian, δmedian = median(αs), median(δs)
        i = argmin(@. hypot(αs - αmedian, δs - δmedian))
        # Sort orbits by nms
        orbits[1], orbits[i] = orbits[i], orbits[1]
        sort!(view(orbits, 2:length(orbits)), by = nms)
        # Nominal orbit
        norbit = jtls(od, orbits[1], params)
        norbits = iszero(norbit) ? view(orbits, 1:1) : [norbit]
    else
        # Sort orbits by nms
        sort!(orbits, by = nms)
        # Nominal orbit
        norbits = view(orbits, 1:1)
    end

    # Save results
    write_neocp_orbits(norbits, view(orbits, 2:length(orbits)), output)
    println("• Output saved to: ", output)

    # Final time
    final_time = now()
    println("• Run started ", initial_time, " and finished ", final_time)
    computation_time = computationtime(initial_time, final_time)
    println("• Total computation time was: ", @sprintf("%.2f", computation_time), " min")

    return nothing
end

main()