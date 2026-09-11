using ArgParse, Distributed, ChunkSplitters, StaticArraysCore
using HTTP, JSON, DataFrames, CSV, Printf, Statistics
@everywhere using NEOs, Dates, PlanetaryEphemeris, TaylorSeries
@everywhere using NEOs: PropresBuffer, PropagationBuffer, OpticalBuffer, OpticalADES,
                  AbstractOrbit, KeplerianElements, ObservatoryMPC, parse_optical_rwo,
                  argoldensearch, evaldeltas, init_optical_residuals, indices,
                  _lsmethods, μ_S, equatorial2ecliptic, _propagate, designation
@everywhere import NEOs: initialcondition, keplerian

const NEOCP_ORBITS_HEADER = "Object   H     G    Epoch    M         Peri.      \
      Node       Incl.        e           n         a                     NObs \
      NOpp   Arc    r.m.s.       Orbit ID"

@everywhere const OD{D, T} = ODProblem{D, T, Vector{OpticalADES{T}},
    Nothing, Veres17{T}, Eggl20{T}}
@everywhere const VariantOrbit{T} = MMOVOrbit{typeof(newtonian!),
    T, T, Vector{OpticalADES{T}}}

function parse_commandline(dict::AbstractDict = Dict())
    s = ArgParseSettings(add_version = true, version = "0.3")

    # Program name (for usage & help screen)
    s.prog = "mcmov.jl"
    # Desciption (for help screen)
    s.description = "Sample the manifold of variations (MOV) of a NEOCP \
    object using jet transport-assisted Monte Carlo."

    s.epilog = """
        Example:\n
        \n
        # Sample the MOV of P22hRXJ with 10 workers and 5 threads each\n
        julia -p 10 -t 5 --project mcmov.jl -i P22hRXJ -t 0000000Hz3SX -s linear --scout\n
        \n
    """

    @add_arg_table! s begin
        # General parameters
        "--input", "-i"
            help = "Input (ADES) file / (NEOCP) designation"
            arg_type = String
        "--output", "-o"
            help = "Output file"
            arg_type = String
        "--trkids", "-t"
            help = "Trkids to include in OD"
            nargs = '+'
        # Manifold of variations parameters
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
        "--refine"
            help = "Refine the first grid"
            action = :store_true
        # Nominal orbits parameters
        "--minimum"
            help = "Find the orbit that minimizes the target function"
            action = :store_true
        "--penalty"
            help = "Find the orbit that minimizes the penalized target function"
            action = :store_true
        "--median"
            help = "Find the orbit closest to the median in RA and DEC"
            action = :store_true
        # Third-party results parameters
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

    args = parse_args(s)

    for pair in dict
        args[first(pair)] = last(pair)
    end

    return args
end

computationtime(x::DateTime, y::DateTime) = @sprintf("%.2f", (y - x).value / 60_000)

printitle(s::AbstractString, d::AbstractString) = println(d ^ length(s), '\n', s,
    '\n', d ^ length(s))

@everywhere initialcondition(x::AbstractOrbit) = x(), epoch(x) + PE.J2000

function chi(x::AbstractVector{VariantOrbit{T}}) where {T <: Real}
    Qmin, i = findmin(nms, x)
    nobs = 2 * noptical(x[i])
    χs = @. sqrt(nobs * ( nms(x) - Qmin ))
    return χs
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
                     points::Vector{NTuple{2, T}}) where {T <: Real}
    xmin, xmax = typemax(T), typemin(T)
    ymin, ymax = typemax(T), typemin(T)
    for (i, point) in enumerate(points)
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

function generate_grid(A::AdmissibleRegion{T}, B::AbstractVector{T},
                       scale::Symbol, Nx::Int, Ny::Int) where {T <: Real}
    points = Vector{NTuple{2, T}}(undef, 0)
    for point in Iterators.product(LinRange(B[1], B[2], Nx), LinRange(B[3], B[4], Ny))
        # Check if point is inside the admissible region
        if scale == :log && (10^point[1], point[2]) in A
            push!(points, (10^point[1], point[2]))
        elseif scale == :linear && (point[1], point[2]) in A
            push!(points, (point[1], point[2]))
        end
    end
    return points
end

function fetch_scout_orbits(desig::AbstractString, write_output::Bool)
    uri = HTTP.URI(
        scheme = "https",
        host   = "ssd-api.jpl.nasa.gov",
        path   = "/scout.api",
        query  = "tdes=$desig&orbits=1"
    )
    response_scout = HTTP.get(string(uri), require_ssl_verification = false)
    text_scout = String(response_scout.body)
    dict_scout = JSON.parse(text_scout)
    orbits_data = dict_scout["orbits"]["data"]
    rows = [collect(row) for row in orbits_data]
    orbits_fields = dict_scout["orbits"]["fields"]
    mat = vcat([permutedims(r) for r in rows]...)  # or hcat(rows...)' as another option
    df = DataFrame(mat, orbits_fields)
    println("• Fetched Scout data")
    if write_output
        CSV.write("$desig.csv", df)
        println("• Scout data saved to: $desig.csv")
        return ""
    else
        io = IOBuffer()
        CSV.write(io, df)
        return String(take!(io))
    end
end

function fetch_neoscan_orbits(desig::AbstractString, write_output::Bool)
    uri = HTTP.URI(
        scheme = "https",
        host   = "newton.spacedys.com",
        path   = "/neodys/NEOScan/scan_neocp/$desig/$desig.mov_sample"
    )
    response_neoscan = HTTP.get(string(uri) #=, require_ssl_verification = false=#)
    println("• Fetched NEOScan data")
    if write_output
        write("$desig.mov_sample", response_neoscan.body)
        println("• NEOScan data saved to: $(desig).mov_sample")
        return ""
    else
        return String(response_neoscan.body)
    end
end

function fetch_neocp_orbits(desig::AbstractString, write_output::Bool)
    url_neocp = "https://cgi.minorplanetcenter.net/cgi-bin/showobsorbs.cgi"
    data = Dict("Obj" => desig, "orb" => "y")
    response_neocp = HTTP.post(url_neocp, [], data #=, require_ssl_verification = false=#)
    text = String(response_neocp.body)
    lines = split(text, '\n')[2:end-2]
    text = join(lines, '\n')
    println("• Fetched NEOCP data")
    if write_output
        write("$desig.orb", text)
        println("• NEOCP data saved to: $(desig).orb")
        return ""
    else
        return text
    end
end

function fetch_neodys_weights(desig::AbstractString)
    url = "https://newton.spacedys.com/neodys/NEOScan/scan_neocp/$desig/$desig.rwo"
    resp = HTTP.get(url)
    text = String(resp.body)
    optical = parse_optical_rwo(text)
    σs = rms.(optical)
    return @. tuple(1 / first(σs), 1 / last(σs))
end

@everywhere function radec_next_day(day_after_epoch::DateTime,
                                    observer::ObservatoryMPC{T},
                                    orbits::AbstractVector{VariantOrbit{T}},
                                    params::Parameters{T}) where {T <: Real}
    radec = Vector{NTuple{2, T}}(undef, length(orbits))
    q0, jd0 = initialcondition(orbits[1])
    t0 = minimum(epoch, orbits) - params.bwdoffset
    tf = dtutc2days(day_after_epoch) + params.fwdoffset
    pbuffer = PropagationBuffer(newtonian!, q0, jd0, (t0, tf), params);
    for (i, orbit) in enumerate(orbits)
        q0, jd0 = initialcondition(orbit)
        tmax = ( tf + PE.J2000 - jd0 ) / yr
        fwd = _propagate(newtonian!, q0, jd0, tmax, pbuffer, params)
        obuffer = OpticalBuffer(zero(T))
        radec[i] = compute_radec(observer, day_after_epoch, obuffer;
            xvs = params.eph_su, xve = params.eph_ea,
            xva = (orbit.bwd, fwd)
        )
    end
    return radec
end

function keplerian(orbit::AbstractOrbit{D, T, T}, t::T,
                   params::Parameters{T}) where {D, T <: Real}
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

function neocp_orbits_format(desig::AbstractString,
                             reference_epoch::Real,
                             ids::AbstractVector{String},
                             orbits::AbstractVector{<:VariantOrbit},
                             params::Parameters)
    orbits_lines = Vector{String}(undef, length(orbits) + 1)
    orbits_lines[1] = NEOCP_ORBITS_HEADER
    for (j, orbit) in enumerate(orbits)
        # Absolute magnitude
        H, _ = absolutemagnitude(orbit, params)
        # Slope parameter
        G = params.slope
        # Orbital elements
        kep = keplerian(orbit, reference_epoch, params)
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
        # Assemble line
        orbits_lines[j + 1] = string(
            rpad(desig, 8),
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
            ids[j]
        )
    end
    return join(orbits_lines, '\n')
end

@everywhere function mcmov(
        od::OD{typeof(newtonian!), T}, A::AdmissibleRegion{T},
        points::AbstractVector{NTuple{2, T}}, bounds::AbstractVector{T},
        scale::Symbol, params::Parameters{T}
    ) where {T <: Real}
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
    # Admissible region epoch [julian days TDB]
    _jd0_ = dtutc2jdtdb(A.date)
    # Initialize buffer and set of residuals
    nobs = 2 * noptical(od)
    idxs = indices(od.tracklets)
    buffer = PropresBuffer(od, AE, _jd0_, idxs, params)
    res = init_optical_residuals(TaylorN{T}, od, idxs)
    # Least squares cache and methods
    x0 = zeros(T, 6)
    lscache = LeastSquaresCache(x0, 1:4, 25)
    lsmethods = _lsmethods(res, x0, 1:4)
    Qtol, Mtol, penalty = params.lsQtol, params.lsMtol, nothing
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
        jd0 = _jd0_ - ae[5] / c_au_per_day
        bwd, fwd = propres!(res, od, q, jd0, params; buffer, idxs)
        if isempty(res)
            res = init_optical_residuals(TaylorN{T}, od, idxs)
            continue
        end
        # Least squares fit
        fit = tryls(res, x0, lscache, lsmethods; penalty, Qtol, Mtol)
        !fit.success && continue
        # Current Q
        Q = nms(res)
        Q(fit.x) < 0 && continue
        # Covariance matrix
        C = (nobs/2) * TS.hessian(Q, fit.x)
        Γ = project(q, fit.x, inv(C))
        # Update orbit
        orbits[i] = evaldeltas(MMOVOrbit(
            newtonian!, variables, od.optical, od.tracklets, bwd, fwd,
            res, Γ, [AE(fit.x);;], [Q(fit.x)]
        ), fit.x)
    end

    return orbits
end

function main(dict::AbstractDict = Dict(); write_output::Bool = true)

    #=================
    General parameters
    =================#

    parsed_args = parse_commandline(dict)

    printitle("Manifold of variations sampling for NEOCP objects", "=")
    printitle("Parameters", "-")

    Nworkers, Nthreads = nworkers(), Threads.nthreads()
    println("• Detected $Nworkers worker(s) with $Nthreads thread(s) each")

    # Load optical astrometry
    input::String = parsed_args["input"]
    optical_all = if isfile(input)
        read_optical_ades(input)
    else
        fetch_optical_ades(input, NEOCP)
    end
    desig = designation(last(optical_all))
    println("• Input (NEOCP) designation: ", desig)

    # If `trkids` is empty use all the astrometry; else, use only trkids contained in `trkids`
    trkids::Vector{String} = parsed_args["trkids"]
    optical = isempty(trkids) ? optical_all : filter(x -> x.trkid in trkids, optical_all)
    println("• `trkids` included in run: ", unique(map(x -> x.trkid, optical_all)))

    if write_output
        orbits_output = something(parsed_args["output"], desig) * ".neos"
        println("• Orbits output file: ", orbits_output)
    else
        orbits_output = ""
    end

    #================================
    Manifold of variations parameters
    ================================#

    scale_str::String = parsed_args["scale"]
    scale::Symbol = Symbol(scale_str)
    @assert scale in (:linear, :log) "Possible values for argument `scale` are: \
        `linear` or `log`"
    println("• Horizontal scale: ", scale)

    Nx::Int, Ny::Int = parsed_args["Nx"], parsed_args["Ny"]
    println("• Number of points in x (y): $Nx ($Ny)")

    χ_max::Float64 = parsed_args["maxchi"]
    println("• χ threshold: ", χ_max)

    refine_grid = parsed_args["refine"]
    println("• Refine the first grid?: ", refine_grid)

    #=======================
    Nominal orbits parameters
    =======================#

    compute_minimum = parsed_args["minimum"]
    println("• Find the orbit that minimizes the target function?: ", compute_minimum)

    compute_penalty = parsed_args["penalty"]
    println("• Find the orbit that minimizes the penalized target function?: ", compute_penalty)

    compute_median = parsed_args["median"]
    println("• Find the orbit closest to the median in RA and DEC?: ", compute_median)

    #=============================
    Third-party results parameters
    =============================#

    fetch_scout = parsed_args["scout"]
    println("• Fetch Scout data?: ", fetch_scout)

    fetch_neoscan = parsed_args["neoscan"]
    println("• Fetch NEOScan data?: ", fetch_neoscan)

    fetch_neocp = parsed_args["neocp"]
    println("• Fetch NEOCP data?: ", fetch_neocp)

    #=====================
    Manifold of variations
    =====================#

    # Initial time
    initial_time = now()
    printitle("Computation", "-")
    println("• Run started at ", initial_time)

    # Orbit determination problem
    od = ODProblem(newtonian!, optical)
    od.weights.weights .= fetch_neodys_weights(desig)
    # Parameters
    params = Parameters(
        maxsteps = 1_000, order = 15, abstol = 1E-12, parse_eqs = true,
        coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
        jtlsorder = 2, jtlsmask = false, jtlsiter = 20, lsiter = 20,
        significance = 0.99, jtlsproject = true, outrej = false,
    )
    # Admissible region
    tracklet = first(od.tracklets)
    A = AdmissibleRegion(tracklet, params)
    # Backward offset must take into consideration the -ρ/c relativistic
    # correction to the epoch
    bwdoffset = params.bwdoffset + A.ρ_domain[2] / c_au_per_day
    params = Parameters(params; bwdoffset)

    # Global box
    bounds = global_box(A, scale)
    println("• Global box: ", bounds)

    # Grid of domain points
    points = generate_grid(A, bounds, scale, Nx, Ny)
    Npoints = length(points)
    println("• $Npoints points in the manifold of variations")

    # Manifold of variations
    orbits = reduce(vcat, pmap(x -> mcmov(od, A, x, bounds, scale, params),
        chunks(points, n = Nworkers)))
    # Eliminate orbits with χ > 5
    χs = chi(orbits)
    mask = χs .≤ χ_max
    keepat!(orbits, mask)
    Norbits = length(orbits)
    println("• $Norbits / $Npoints points with χ ≤ χ_max = $χ_max")

    if refine_grid
        # Global box
        bounds .= refined_box(mask, scale, points)
        println("• Refined box: ", bounds)
        # Grid of domain points
        points = generate_grid(A, bounds, scale, Nx, Ny)
        Npoints = length(points)
        println("• $Npoints points in the manifold of variations")
        # Manifold of variations
        orbits = reduce(vcat, pmap(x -> mcmov(od, A, x, bounds, scale, params),
            chunks(points, n = Nworkers)))
        # Eliminate orbits with χ > 5
        χs = chi(orbits)
        mask = χs .≤ χ_max
        keepat!(orbits, mask)
        Norbits = length(orbits)
        println("• $Norbits / $Npoints points with χ ≤ χ_max = $χ_max")
    end

    #=============
    Nominal orbits
    =============#

    nominal_ids = Vector{String}(undef, 0)
    nominal_orbits = Vector{VariantOrbit{Float64}}(undef, 0)

    if compute_minimum
        minimum_orbit = argmin(nms, orbits)
        push!(nominal_ids, "NEOCPMinim")
        push!(nominal_orbits, minimum_orbit)
        printitle("Orbit that minimizes the target function", "-")
        println(summary(minimum_orbit))
    end

    if compute_penalty
        params = Parameters(params; lspenalty = 0.05)
        ρ = if scale === :linear
            (A.ρ_domain[1] + A.ρ_domain[2]) / 2
        else
            10^(log10(A.ρ_domain[1] * A.ρ_domain[2]) / 2)
        end
        v_ρ = (A.v_ρ_domain[1] + A.v_ρ_domain[2]) / 2
        penalty_orbit = mmov(od, A, ρ, v_ρ, params; scale)
        params = Parameters(params; lspenalty = 0.00)
        push!(nominal_ids, "NEOCPPenal")
        push!(nominal_orbits, penalty_orbit)
        printitle("Orbit that minimizes the penalized target function", "-")
        println(summary(penalty_orbit))
    end

    if compute_median
        day_after_epoch = date(tracklet) + Day(1)
        observer = observatory(tracklet)
        radec = reduce(vcat, pmap(x -> radec_next_day(day_after_epoch, observer, x, params),
            chunks(orbits, n = Nworkers)))
        αs, δs = @. first(radec), last(radec)
        αmedian, δmedian = median(αs), median(δs)
        i = argmin(@. hypot(αs - αmedian, δs - δmedian))
        median_orbit = orbits[i]
        push!(nominal_ids, "NEOCPMedia")
        push!(nominal_orbits, penalty_orbit)
        printitle("Orbit closest to the median in RA and DEC", "-")
        println(summary(median_orbit))
    end

    #==================
    Third-party results
    ==================#

    scout_string = fetch_scout ? fetch_scout_orbits(desig, write_output) : ""
    neoscan_string = fetch_neoscan ? fetch_neoscan_orbits(desig, write_output) : ""
    neocp_string = fetch_neocp ? fetch_neocp_orbits(desig, write_output) : ""

    #===========
    Save results
    ===========#

    # Sort by nms, then put the nominal orbits first
    sort!(orbits, by = nms)
    prepend!(orbits, nominal_orbits)
    unique!(orbits)
    N_nominal, N_orbits = length(nominal_orbits), length(orbits)
    ids = ["NEOCPV" * lpad(i-N_nominal, 4, '0') for i in N_nominal+1:N_orbits]
    prepend!(ids, nominal_ids)
    # Print orbits in NEOCP format
    reference_epoch = dtutc2days(tracklet)
    orbits_string = neocp_orbits_format(desig, reference_epoch, ids, orbits, params)
    if write_output
        write(orbits_output, orbits_string)
        println("• Orbits saved to: ", orbits_output)
    end

    # Final time
    final_time = now()
    println("• Run started $initial_time and finished $final_time")
    computation_time = computationtime(initial_time, final_time)
    println("• Total computation time was: $computation_time min")

    return orbits_string, scout_string, neoscan_string, neocp_string
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end