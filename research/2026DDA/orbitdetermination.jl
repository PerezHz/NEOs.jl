using Distributed, ArgParse
using Random: shuffle!

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "orbitdetermination.jl"
    # Desciption (for help screen)
    s.description = "Large scale orbit determination for NEOs via jet transport"

    s.epilog = """
        Example:\n
        \n
        julia -p 10 -t 5 --project orbitdetermination.jl -i names.txt -o orbits/\n
        \n
    """

    @add_arg_table! s begin
        "--input", "-i"
            help = "input asteroids file"
            arg_type = String
        "--output", "-o"
            help = "output directory"
            arg_type = String
        "--neocc"
            help = "download NEOCC orbit"
            action = :store_true
        "--neodys"
            help = "download NEODyS orbit"
            action = :store_true
        "--jpl"
            help = "download JPL orbit"
            action = :store_true
    end

    return parse_args(s)
end

@everywhere using Printf
@everywhere using Accessors: @reset
@everywhere using Parameters: @unpack

@everywhere begin
    using NEOs, PlanetaryEphemeris, TaylorSeries, Dates, Statistics
    using JLD2, HTTP, JSON, ThreadPinning

    using NEOs: NEOCC_URL, AbstractWeightingScheme, OpticalADES, OpticalTracklet,
          AdmissibleRegion, OpticalResidual, indices, set_od_order, updateorbit,
          rexveres17, σsveres17

    import NEOs: weights, corr, getid, update!

    const MINUTC = days2dtutc(0.0)
    const NEODyS2_URL = "https://newton.spacedys.com/~neodys2/"
    const SBDB_JPL_API = "https://ssd-api.jpl.nasa.gov/sbdb.api"
    const Orbit = LeastSquaresOrbit{typeof(newtonian!), Float64, Float64,
        Vector{OpticalADES{Float64}}, Nothing, Nothing}

    computationtime(x::DateTime, y::DateTime) = @sprintf("%.2f", (y - x).value / 60_000)

    timestring(x::DateTime) = string("(", computationtime(x, now()), " min)")

    printitle(s::AbstractString, d::AbstractString) = println(d ^ length(s),
        '\n', s, '\n', d ^ length(s))

    chi(x::OpticalResidual) = sqrt(chi2(x))
    logchi(x::OpticalResidual) = log(chi(x))

    daysbetween(x::DateTime, y::DateTime) = abs(x - y)
    daysbetween(x::DateTime, y::OpticalTracklet) = daysbetween(x, date(y))
    daysbetween(x::OpticalTracklet, y::DateTime) = daysbetween(date(x), y)
    daysbetween(x::OpticalTracklet, y::OpticalTracklet) = daysbetween(date(x), date(y))

    isodvalid(od::ODProblem, orbit::LeastSquaresOrbit, params::Parameters) =
        noptical(od) == noptical(orbit) && critical_value(orbit) < params.significance &&
        all(!isnan, sigmas(orbit))

    function meanepoch(x::ODProblem)
        t = Vector{Float64}(undef, noptical(x))
        w = Vector{Float64}(undef, noptical(x))
        for i in eachindex(x.optical)
            t[i] = dtutc2days(x.optical[i])
            δ = dec(x.optical[i])
            σα, σδ = inv.(x.weights.weights[i])
            w[i] = 1 / (σα^2 * cos(δ)^2 + σδ^2)
        end
        return mean(t, weights(w))
    end

    meandate(x::ODProblem) = days2dtutc(meanepoch(x))

    # Error model
    mutable struct ModifiedVeres17{T} <: AbstractWeightingScheme{T}
        weights::Vector{NTuple{2, T}}
        corr::Vector{T}
    end

    function ModifiedVeres17(optical::AbstractVector{OpticalADES{T}}) where {T <: Real}
        weights = w8smveres17(optical)
        corrs = corr.(optical)
        return ModifiedVeres17{T}(weights, corrs)
    end

    weights(x::ModifiedVeres17) = x.weights
    corr(x::ModifiedVeres17) = x.corr
    getid(::ModifiedVeres17) = "Modified Veres et al. (2017)"

    function update!(x::ModifiedVeres17{T}, optical::AbstractVector{OpticalADES{T}}) where {T <: Real}
        x.weights = w8smveres17(optical)
        x.corr = corr.(optical)
        return nothing
    end

    function σsmveres17(obs::OpticalADES{T}) where {T <: Real}
        σ = σsveres17(obs)
        σα = isnan(obs.rmsra) ? σ : max(obs.rmsra, σ)
        σδ = isnan(obs.rmsdec) ? σ : max(obs.rmsdec, σ)
        return σα, σδ
    end

    function w8smveres17(optical::AbstractVector{OpticalADES{T}}) where {T <: Real}
        σs = σsmveres17.(optical)
        rex = rexveres17(optical)
        return @. tuple(1 / (rex * first(σs)), 1 / (rex * last(σs)))
    end

    function on_error(x)
        println("• Orbit determination failed for an unknown asteroid")
        return false
    end

    function fetch_neocc_orbit(neo::AbstractString)
        # Assemble URL
        url = NEOCC_URL * "PSDB-portlet/download?file=$neo.eq0"
        # HTTP response
        resp = HTTP.get(url)
        # Convert to String
        text = String(resp.body)
        return text
    end

    function fetch_neodys_orbit(neo::AbstractString)
        # Assemble URL
        url = NEODyS2_URL * "epoch/$neo.eq0"
        # HTTP response
        resp = HTTP.get(url)
        # Convert to String
        text = String(resp.body)
        return text
    end

    function fetch_jpl_orbit(neo::AbstractString)
        # Assemble URL
        url = SBDB_JPL_API * "?sstr=$neo&full-prec=1&cov=mat&phys-par=1"
        # HTTP response
        resp = HTTP.get(url; require_ssl_verification = false)
        # Convert to String
        text = String(resp.body)
        return text
    end

    # Naive initial conditions for iod
    function initcond(A::AdmissibleRegion)
        v_ρ = sum(A.v_ρ_domain) / 2
        return [
            (A.ρ_domain[1], v_ρ, :log),
            (10^(sum(log10, A.ρ_domain) / 2), v_ρ, :log),
            (sum(A.ρ_domain) / 2, v_ρ, :linear),
            (A.ρ_domain[2], v_ρ, :linear),
        ]
    end

    function preprocessing!(optical::AbstractVector{OpticalADES{T}}) where {T <: Real}
        filter!(!isdeprecated, optical)
        if numberofdays(optical) > 30
            L0, Lf = 1, length(optical)
            j0 = something(findfirst(isdiscovery, optical), L0)
            d0 = date(optical[j0])
            jf = something(findlast(x -> daysbetween(d0, date(x)) ≤ Day(30), optical), Lf)
            keepat!(optical, j0:jf)
        end
        @unpack provid = optical[1]
        provid = replace(provid, " " => "")
        if provid == "2022PW40"
            for (i, x) in enumerate(optical)
                @reset x.trkid = x.trkid * string(ceil(Int, i / 3))
                optical[i] = x
            end
        elseif provid == "2015EQ"
            filter!(x -> observatorycode(x) != "703", optical)
        elseif provid == "2017UR52"
            filter!(x -> observatorycode(x) != "568", optical)
        elseif provid == "2025OB22"
            filter!(x -> observatorycode(x) != "U68", optical)
        elseif provid == "2020XL1"
            deleteat!(optical, 1:2)
        end
        return nothing
    end

    function parameters(OD::ODProblem)
        Nt = length(OD.tracklets)
        if Nt == 1
            χ2_rec, χ2_rej, fudge, max_per = 7.378, 7.824, 400.0, 10.0
        elseif Nt == 2
            χ2_rec, χ2_rej, fudge, max_per = 5.991, 7.378, 200.0, 20.0
        else
            χ2_rec, χ2_rej, fudge, max_per = 4.605, 5.991, 100.0, 30.0
        end
        return Parameters(;
            maxsteps = 100_000, order = 15, abstol = 1E-12, parse_eqs = true,
            coeffstol = Inf, bwdoffset = 0.05, fwdoffset = 0.05,
            gaussorder = 2, safegauss = false, refscale = :log,
            tsaorder = 2, adamiter = 500, adamQtol = 1E-5,
            jtlsorder = 2, jtlsmask = false, jtlsiter = 20, lsiter = 10,
            jtlsproject = true, significance = 0.99, verbose = true,
            outrej = true, χ2_rec, χ2_rej, fudge, max_per
        )
    end

    function initialtracklets(OD::ODProblem)
        # Initial tracklet
        i = findfirst(x -> !isempty(x.disc), OD.optical)
        j = isnothing(i) ? 1 : findfirst(x -> i in x.indices, OD.tracklets)
        # Sort OD.tracklets by the time to the initial tracklet
        tracklets = sort(OD.tracklets, by = Base.Fix1(daysbetween, OD.tracklets[j]))
        return tracklets
    end

    function newobservations(OD::ODProblem, od::ODProblem, res::AbstractVector)
        # Maximum chi
        trks = setdiff(OD.tracklets, od.tracklets)
        mags = Vector{Int}(undef, length(trks))
        for i in eachindex(mags)
            idxs = indices(trks[i])
            x = maximum(logchi, view(res, idxs))
            mags[i] = ceil(Int, x)
        end
        n = max(2, minimum(mags, init = typemax(Int)))
        cmax = exp(n)
        # New observations
        idxs = Int.(indexin(od.optical, OD.optical))
        if n ≥ 10
            trk = argmin(Base.Fix1(daysbetween, meandate(od)), trks)
            union!(idxs, indices(trk))
        else
            mask = @. mags ≤ n
            union!(idxs, indices(view(trks, mask)))
        end
        sort!(idxs)
        return idxs, cmax
    end

    function initial_orbit_determination(od::ODProblem, params::Parameters)
        # Unpack
        significance = params.significance
        # Pre-allocate orbit
        orbit = zero(Orbit)
        # Iterate tracklets
        for i in eachindex(od.tracklets)
            # Minimization over the MOV requires a minimum of 2 observations
            nobs(od.tracklets[i]) < 2 && continue
            # Admissible region
            A = AdmissibleRegion(od.tracklets[i], params)
            # List of naive initial conditions
            I0 = initcond(A)
            # Iterate naive initial conditions
            for (ρ, v_ρ, scale) in I0
                # Minimization over the MOV
                porbit = mmov(od, A, ρ, v_ρ, params; i, scale)
                # Failed to converge
                iszero(porbit) && continue
                # Jet Transport Least Squares
                _orbit_ = jtls(od, porbit, params, true)
                # Update orbit
                C1, C2 = critical_value(porbit), critical_value(_orbit_)
                if (significance < C1) && (significance < C2)
                    orbit = updateorbit(orbit, _orbit_, od.optical)
                    continue
                else
                    N1, N2 = length(porbit.Qs), length(_orbit_.Qs)
                    M1 = string("* Minimization over the MOV converged in $N1 \
                        iterations to:\n\n", summary(porbit))
                    if C1 < C2
                        orbit = LeastSquaresOrbit(od, porbit(), epoch(porbit) + J2000, params)
                        M2 = ""
                    else
                        orbit = _orbit_
                        M2 = string("\n", "* Jet Transport Least Squares converged in $N2 \
                            iterations to: \n\n", summary(orbit))
                    end
                    println(M1, M2)
                    return orbit
                end
            end
        end
        # Last resort: Gauss method
        orbit = updateorbit(orbit, gaussiod(od, params), od.optical)
        return orbit
    end

    function orbit_determination(od::ODProblem, orbit::LeastSquaresOrbit, params::Parameters)
        # Try #1: linkage
        params = Parameters(params; outrej = false, verbose = false)
        orbit1 = linkage(od, orbit, params)
        params = Parameters(params; outrej = true, verbose = true)
        if isodvalid(od, orbit1, params)
            orbit1 = jtls(od, orbit1, params)
            if isodvalid(od, orbit1, params)
                println(
                    "* Linkage converged to: \n\n",
                    summary(orbit1)
                )
                return orbit1
            else
                @warn("Linkage did not converge within the given parameters or \
                    could not fit all the astrometry")
            end
        end
        # Try #2: repeat initial orbit determination
        params = Parameters(params; outrej = false)
        orbit2 = initial_orbit_determination(od, params)
        params = Parameters(params; outrej = true)
        isodvalid(od, orbit2, params) && return orbit2
        orbit2 = initial_orbit_determination(od, params)
        return orbit2
    end

    function orbit_determination(neo::AbstractString; niter::Int = 20,
                                 dtmin::Second = Second(60))
        # Random start delay
        sleep(Second(rand(1:Nworkers)))
        # Fetch optical astrometry
        optical = fetch_optical_ades(neo, MPC)
        if any(x -> date(x) < MINUTC, optical)
            println("• Asteroid $neo skipped: observations before J2000")
            sleep(dtmin)
            return false
        end
        preprocessing!(optical)
        # Create output directory
        dirname = joinpath(output, neo)
        mkdir(dirname)
        # Fetch third party orbits
        if neocc
            orbitNEOCC = fetch_neocc_orbit(neo)
            filename = joinpath(dirname, neo * "NEOCC.txt")
            write(filename, orbitNEOCC)
        end
        if neodys
            orbitNEODyS = fetch_neodys_orbit(neo)
            filename = joinpath(dirname, neo * "NEODyS.txt")
            write(filename, orbitNEODyS)
        end
        if jpl
            orbitJPL = fetch_jpl_orbit(neo)
            filename = joinpath(dirname, neo * "JPL.json")
            JSON.json(filename, JSON.parse(orbitJPL))
        end
        # Outer initial time
        outer_initial_time = now()
        # Redirect output to a file
        filename = joinpath(dirname, "nohupOrbit.out")
        flag = redirect_stdio(; stdout = filename, stderr = filename) do
            # Inner initial time
            inner_initial_time = now()
            printitle("Orbit determination for asteroid $neo", "-")
            println("• Detected 1 worker with ",  Threads.nthreads(), " thread(s)")
            println("• Run started at ", inner_initial_time)
            # Orbit determination
            flag = orbit_determination(neo, optical; niter)
            # Inner final time
            inner_final_time = now()
            println("• Run finished at ", inner_final_time)
            inner_computation_time = computationtime(inner_initial_time, inner_final_time)
            println("• Computation time was ", inner_computation_time, " min")
            return flag
        end
        !flag && rm(dirname, recursive = true)
        # Outer final time
        outer_time_string = timestring(outer_initial_time)
        verb = flag ? "succeeded" : "failed"
        println("• Orbit determination $verb for asteroid $neo $outer_time_string")
        # Sleep if computation time was less than dtmin
        dt = floor(now() - outer_initial_time, Second)
        dt < dtmin && sleep(dtmin - dt)
        return flag
    end

    function orbit_determination(neo::AbstractString, optical::AbstractVector;
                                 niter::Int = 20)
        # Global orbit determination problem and parameters
        OD = ODProblem(newtonian!, optical, weights = ModifiedVeres17, debias = Eggl20)
        params = parameters(OD)
        # Initial orbit determination
        od = deepcopy(OD)
        orbit = zero(Orbit)
        tracklets = initialtracklets(OD)
        printitle("Iteration 1/$niter (cmax = 0.0)", "*")
        for i in eachindex(tracklets)
            trks = view(tracklets, 1:i)
            nobs(trks) < 3 && continue
            idxs = indices(trks)
            NEOs.update!(od, optical[idxs])
            orbit = initial_orbit_determination(od, params)
            isodvalid(od, orbit, params) && break
        end
        # Observations included in each iteration
        obs = Vector{Set{OpticalADES{Float64}}}(undef, niter+1)
        obs[1] = Set(od.optical)
        # Main loop
        cmax = 0.0
        for i in 1:niter
            # Orbit determination
            if i > 1
                printitle("Iteration $i/$niter (cmax = $cmax)", "*")
                orbit = orbit_determination(od, orbit, params)
            end
            iszero(orbit) && break
            # Break condition
            if noptical(orbit) == length(optical)
                # Last resort: direct Gauss method to a all observations
                if !isodvalid(od, orbit, params)
                    printitle("Last resort", "*")
                    orbit = gaussiod(od, params)
                end
                # Shift epoch to the middle of the observational arc
                tmean = meanepoch(od)
                orbit = shiftepoch(orbit, tmean + PE.J2000, params)
                iszero(orbit) && break
                _orbit_ = jtls(od, orbit, params)
                orbit = isodvalid(od, _orbit_, params) ? _orbit_ : orbit
                iszero(orbit) && break
                printitle("Final orbit", "*")
                println(summary(orbit))
                # Save output
                dirname = joinpath(output, neo)
                filename = joinpath(dirname, neo * "Orbit.jld2")
                jldsave(filename; orbit)
                return true
            end
            # Astrometric residuals wrt OD
            _, _, res = propres(OD, orbit(), epoch(orbit) + PE.J2000, params)
            isempty(res) && break
            # Update orbit determination problem
            idxs, cmax = newobservations(OD, od, res)
            NEOs.update!(od, optical[idxs])
            obs[i+1] = Set(setdiff(od.optical, view(obs, 1:i)...))
        end
        return false
    end
end

function main()
    # Parse arguments from commandline
    parsed_args = parse_commandline()

    # Print header
    printitle("Large scale orbit determination for NEOs via jet transport", "=")

    # Number of workers and threads
    Nworkers, Nthreads = nworkers(), Threads.nthreads()
    println("• Detected ", Nworkers, " worker(s) with ", Nthreads, " thread(s) each")

    # Pin threads to sockets
    distributed_pinthreads(:sockets)

    # Input asteroids file
    input::String = parsed_args["input"]
    println("• Input asteroids file: ", input)

    # Output directory
    output::String = parsed_args["output"]
    println("• Output directory: ", output)

    # Download NEOCC orbit?
    neocc::Bool = parsed_args["neocc"]
    println("• Download NEOCC orbit?: ", neocc)

    # Download NEODyS orbit?
    neodys::Bool = parsed_args["neodys"]
    println("• Download NEODyS orbit?: ", neodys)

    # Download JPL orbit?
    jpl::Bool = parsed_args["jpl"]
    println("• Download JPL orbit?: ", jpl)

    @everywhere begin
        # Declare constants in all workers
        const output = $output
        const neocc = $neocc
        const neodys = $neodys
        const jpl = $jpl
        const Nworkers = $Nworkers
        const Nthreads = $Nthreads
        # Update observatories and catalogues
        update_catalogues_mpc()
        update_observatories_mpc()
        # Set jet transport variables
        set_od_order(Float64, 2, 6)
    end

    # Global initial time
    global_initial_time = now()
    println("• Run started at ", global_initial_time)

    # Parse asteroid designations
    neos = readlines(input)
    println("• Detected ", length(neos), " NEOs")

    # Distributed orbit determination
    shuffle!(neos)
    mask = pmap(orbit_determination, neos; on_error)
    println("• Orbit determination succeeded in ", count(mask), " out of ",
            length(neos), " asteroids")

    # Final time
    global_final_time = now()
    println("• Run started ", global_initial_time, " and finished ", global_final_time)
    global_computation_time = computationtime(global_initial_time, global_final_time)
    println("• Total computation time was: ", global_computation_time, " min")

    return nothing
end

main()