using Distributed, ArgParse, Dates, JSON, JLD2
using Plots, Colors, Printf

ENV["GKSwstype"] = "100"

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "largescale.jl"
    # Desciption (for help screen)
    s.description = "Orbid determination large scale comparison with ESA and JPL"

    @add_arg_table! s begin
        "--directory", "-d"
            help = "I/O directory"
            arg_type = String
    end

    s.epilog = """
        Example:\n
        \n
        julia -p 10 -t 5 --project largescale.jl -d orbits/\n
        \n
    """

    return parse_args(s)
end

const UNAM_CYAN = RGB(6/255  , 151/255, 213/255)
const UNAM_BLUE = RGB(22/255 , 40/255 , 102/255)
const UNAM_GOLD = RGB(220/255, 167/255, 29/255 )

computationtime(x::DateTime, y::DateTime) = @sprintf("%.2f", (y - x).value / 60_000)

printitle(s::AbstractString, d::AbstractString) = println(d ^ length(s),
        '\n', s, '\n', d ^ length(s))

function parseESA(x::AbstractString)
    y = split(x)
    return string(y[1]), parse.(Float64, view(y, 2:length(y)))
end

function parseJPL(x::AbstractVector)
    a = string(replace(x[1], " " => ""))
    b = parse.(Float64, view(x, 2:length(x)))
    return a, b
end

@everywhere begin
    using NEOs, PlanetaryEphemeris, TaylorSeries, LinearAlgebra
    using NEOs: PropagationBuffer, OpticalADES, KeplerianElements, μ_S, set_od_order,
          elements, covariance, equatorial2ecliptic, evaldeltas

    const Orbit = LeastSquaresOrbit{typeof(newtonian!), Float64, Float64,
        Vector{OpticalADES{Float64}}, Nothing, Nothing}

    mahalanobis(x::AbstractVector, μ::AbstractVector, Σ::AbstractMatrix) =
        sqrt((x - μ)' * inv(Diagonal(Σ)) * (x - μ))

    function mahalanobis(orbitNEOs::LeastSquaresOrbit, orbitESA::AbstractVector,
                         orbitJPL::AbstractVector)
        # Reference epochs [JDTDB]
        jd0JPL = orbitJPL[1]
        jd0NEOs = epoch(orbitNEOs) + J2000
        jd0ESA = orbitESA[1] + (J2000 - MJD2000)
        # Reference epochs [MJDTDB]
        mjd0ESA = orbitESA[1]
        mjd0JPL = orbitJPL[1] + (MJD2000 - J2000)
        # Keplerian elements for ESA
        elementsESA = orbitESA[[2, 3, 4, 6, 5, 7]]
        # Keplerian elements for NEOs at ESA's epoch
        ndays = jd0ESA - jd0NEOs
        ndays += sign(ndays) * max(params.fwdoffset, params.bwdoffset)
        q0 = orbitNEOs() + 1E-8 * get_variables(Float64, 2)
        prop = NEOs._propagate(gravityonly!, q0, jd0NEOs, ndays / yr, bufferTN, params)
        q0 = equatorial2ecliptic(prop(jd0ESA - J2000) - params.eph_su(jd0ESA - J2000))
        ele = cartesian2keplerian(q0, mjd0ESA)
        Γ_ele = project(ele, covariance(orbitNEOs))
        kep = evaldeltas(KeplerianElements{Float64, TaylorN{Float64}}(
            μ_S, mjd0ESA, :ecliptic, ele, Γ_ele), zeros(Float64, 6))
        elementsNEOs = elements(kep)
        covarianceNEOs = covariance(kep)
        # Mahalanobis distance to ESA
        MESA = mahalanobis(elementsESA, elementsNEOs, covarianceNEOs)
        # Keplerian elements for JPL
        elementsJPL = orbitJPL[[3, 2, 5, 7, 6, 8]]
        # Keplerian elements for NEOs at JPL's epoch
        ndays = jd0JPL - jd0NEOs
        ndays += sign(ndays) * max(params.fwdoffset, params.bwdoffset)
        q0 = orbitNEOs() + 1E-8 * get_variables(Float64, 2)
        prop = NEOs._propagate(gravityonly!, q0, jd0NEOs, ndays / yr, bufferTN, params)
        q0 = equatorial2ecliptic(prop(jd0JPL - J2000) - params.eph_su(jd0JPL - J2000))
        ele = cartesian2keplerian(q0, mjd0JPL)
        Γ_ele = project(ele, covariance(orbitNEOs))
        kep = evaldeltas(KeplerianElements{Float64, TaylorN{Float64}}(
            μ_S, mjd0JPL, :ecliptic, ele, Γ_ele), zeros(Float64, 6))
        elementsNEOs = elements(kep)
        covarianceNEOs = covariance(kep)
        # Mahalanobis distance to JPL
        MJPL = mahalanobis(elementsJPL, elementsNEOs, covarianceNEOs)

        return MESA, MJPL
    end

    function orbitstatistics(orbitNEOs::LeastSquaresOrbit, orbitESA::AbstractVector,
                             orbitJPL::AbstractVector)
        try
            Q = nrms(orbitNEOs)
            S = minimum(snr(orbitNEOs); init = Inf)
            U = uncertaintyparameter(orbitNEOs, params)
            MESA, MJPL = mahalanobis(orbitNEOs, orbitESA, orbitJPL)
            return Q, S, U, MESA, MJPL
        catch
            return NaN, NaN, -1, NaN, NaN, NaN
        end
    end
end

function main()
    # Parse arguments from commandline
    parsed_args = parse_commandline()

    # Print header
    printitle("Orbid determination large scale comparison with ESA and JPL", "=")

    # Number of workers and threads
    println("• Detected ", nworkers(), " worker(s) with ", Threads.nthreads(),
            " thread(s) each")

    # I/O directory
    directory::String = parsed_args["directory"]
    println("• I/O directory: ", directory)

    # Global initial time
    global_initial_time = now()
    println("• Run started at ", global_initial_time)

    @everywhere begin
        # Set jet transport variables
        set_od_order(Float64, 2, 6)
        # Parameters
        const params = Parameters(
            maxsteps = 10_000, order = 25, abstol = 1E-20, parse_eqs = true,
            coeffstol = Inf, bwdoffset = 0.05, fwdoffset = 0.05,
            gaussorder = 2, safegauss = false, refscale = :log,
            tsaorder = 2, adamiter = 500, adamQtol = 1E-5,
            jtlsorder = 2, jtlsmask = false, jtlsiter = 20, lsiter = 10,
            significance = 0.99, verbose = true, outrej = true,
            χ2_rec = 7.0, χ2_rej = 8.0, fudge = 100.0, max_per = 20.0
        )
        # Propagation buffers
        jd0 = 2459200.5
        tlim = (0.0, 9861.5)
        q00 = [−0.18034828526, 0.94069105951, 0.34573599029,
               −0.0162659397882, 4.39154800E−5, −0.000395204013]
        q0 = q00 + 1E-8 * get_variables(Float64, 2)
        const bufferTN = PropagationBuffer(gravityonly!, q0, jd0, tlim, params)
    end

    # Parse successful NEOs
    filename = joinpath(directory, "nohupOD.out")
    text = read(filename, String)
    re = r"Orbit determination succeeded for asteroid (?<neo>[\d\w]+) \((?<t>\d*\.\d*) min\)"
    ms = collect(eachmatch(re, text))
    neos = first.(ms)
    sort!(neos)

    # Parse computation times and NEOs orbits
    re = r"Computation time was (?<t>\d*\.\d*) min"
    ts = Vector{Float64}(undef, length(neos))
    orbitsNEOs = Vector{Orbit}(undef, length(neos))
    for (i, neo) in enumerate(neos)
        filename = joinpath(directory, neo, "nohupOrbit.out")
        text = read(filename, String)
        m = match(re, text)
        ts[i] = parse(Float64, first(m))
        filename = joinpath(directory, neo, neo * "Orbit.jld2")
        orbitsNEOs[i] = JLD2.load(filename, "orbit")
    end

    # Parse ESA orbits file
    filename = joinpath(directory, "neo_km.cat")
    lines = readlines(filename)[7:end]
    neosESA = Vector{String}(undef, length(lines))
    _orbitsESA_ = Vector{Vector{Float64}}(undef, length(lines))
    for (i, line) in enumerate(lines)
        neosESA[i], _orbitsESA_[i] = parseESA(line)
    end
    perm = sortperm(neosESA)
    permute!(neosESA, perm)
    permute!(_orbitsESA_, perm)
    orbitsESA = Vector{Vector{Float64}}(undef, length(neos))
    for (i, neo) in enumerate(neos)
        idxs = searchsorted(neosESA, neo)
        orbitsESA[i] = _orbitsESA_[first(idxs)]
    end

    # Parse JPL orbits file
    filename = joinpath(directory, "neos.json")
    dict = JSON.parse(read(filename))
    neosJPL = Vector{String}(undef, length(dict["data"]))
    _orbitsJPL_ = Vector{Vector{Float64}}(undef, length(dict["data"]))
    for (i, line) in enumerate(dict["data"])
        neosJPL[i], _orbitsJPL_[i] = parseJPL(line)
    end
    perm = sortperm(neosJPL)
    permute!(neosJPL, perm)
    permute!(_orbitsJPL_, perm)
    orbitsJPL = Vector{Vector{Float64}}(undef, length(neos))
    for (i, neo) in enumerate(neos)
        idxs = searchsorted(neosJPL, neo)
        orbitsJPL[i] = _orbitsJPL_[first(idxs)]
    end

    # Orbit statistics
    stats = pmap(orbitstatistics, orbitsNEOs, orbitsESA, orbitsJPL)
    Qs = first.(stats)
    Ss = getindex.(stats, 2)
    Us = getindex.(stats, 3)
    MsESA = getindex.(stats, 4)
    MsJPL = getindex.(stats, 5)
    # Save orbit statstics
    filename = joinpath(directory, "stats.jld2")
    jldsave(filename; ts, Qs, Ss, Us, MsESA, MsJPL)
    println("• Saved orbit statistics: ", filename)

    # Plots.jl settings
    default(dpi = 300, titlefont = (12), guidefont = (12, :black),
            tickfont = (8, :black), framestyle = :box, yminorgrid = false,
            xminorgrid = false, legend = true #=margin = 0.25Measures.mm,=#)

    # Computation times histogram
    bins = 0:5:60
    histogram(ts .* 60, bins = bins, color = UNAM_CYAN, label = "", linewidth = 0.5,
        xlabel = "Computation time [s]", ylabel = "Count", xlim = extrema(bins),
        xticks = bins, ylim = (0, 8_000), yticks = 0:1_000:8_000)
    filename = joinpath(directory, "THIST.png")
    savefig(filename)
    println("• Saved computation times histogram to: ", filename)

    # NRMS histogram
    bins = 0:0.1:1
    histogram(Qs, bins = bins, color = UNAM_CYAN, label = "", linewidth = 0.5,
        xlabel = "NRMS", ylabel = "Count", xlim = extrema(bins), xticks = bins,
        ylim = (0, 6_000), yticks = 0:1_000:6_000)
    filename = joinpath(directory, "QHIST.png")
    savefig(filename)
    println("• Saved NRMS histogram to: ", filename)

    # SNR histogram
    bins = -3:0.5:7
    histogram(log10.(Ss), bins = bins, color = UNAM_CYAN, label = "", linewidth = 0.5,
        xlabel = "log₁₀(SNR)", ylabel = "Count", xlim = extrema(bins), xticks = -3:7,
        ylim = (0, 7_000), yticks = 0:1_000:7_000)
    filename = joinpath(directory, "SHIST.png")
    savefig(filename)
    println("• Saved SNR histogram to: ", filename)

    # Uncertainty parameters distance
    bins = -0.5:9.5
    histogram(Us, bins = bins, color = UNAM_CYAN, label = "", linewidth = 0.5,
        xlabel = "Uncertainty parameter", ylabel = "Count", xlim = extrema(bins), xticks = 0:9,
        ylim = (0, 8_000), yticks = 0:1_000:8_000)
    filename = joinpath(directory, "UHIST.png")
    savefig(filename)
    println("• Saved uncertainty parameters histogram to: ", filename)

    # ESA Mahalanobis distance
    bins = 0:0.25:3
    histogram(MsESA, bins = bins, color = UNAM_CYAN, label = "", linewidth = 0.5,
        xlabel = "Mahalanobis distance", ylabel = "Count", xlim = extrema(bins), xticks = bins,
        ylim = (0, 9_000), yticks = 0:1_000:9_000)
    filename = joinpath(directory, "MESAHIST.png")
    savefig(filename)
    println("• Saved ESA Mahalanobis distance histogram to: ", filename)

    # JPL Mahalanobis distance
    bins = 0:0.25:3
    histogram(MsJPL, bins = bins, color = UNAM_CYAN, label = "", linewidth = 0.5,
        xlabel = "Mahalanobis distance", ylabel = "Count", xlim = extrema(bins), xticks = bins,
        ylim = (0, 9_000), yticks = 0:1_000:9_000)
    filename = joinpath(directory, "MJPLHIST.png")
    savefig(filename)
    println("• Saved JPL Mahalanobis distance histogram to: ", filename)

    # Final time
    global_final_time = now()
    println("• Run started ", global_initial_time, " and finished ", global_final_time)
    global_computation_time = computationtime(global_initial_time, global_final_time)
    println("• Total computation time was: ", global_computation_time, " min")

    return nothing
end

main()