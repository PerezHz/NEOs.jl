using ArgParse, NEOs, PlanetaryEphemeris, DataFrames, Statistics, Dates, JLD2
using Plots, Colors, Printf

using NEOs: AbstractWeightingScheme, AbstractOpticalAstrometry, AbstractOpticalVector,
      OpticalADES, OpticalResidual, σsveres17, rexveres17, skipnanmean, indices

import NEOs: weights, corr, getid, update!

ENV["GKSwstype"] = "100"

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "3I.jl"
    # Desciption (for help screen)
    s.description = "Orbit determination for 3I/ATLAS"

    @add_arg_table! s begin
        "--directory", "-d"
            help = "I/O directory"
            arg_type = String
    end

    s.epilog = """
        Example:\n
        \n
        julia -t 5 --project 3I.jl -d data/\n
        \n
    """

    return parse_args(s)
end

const HIGH_FIDELITY_OBSERVATORIES = search_observatory_code.([
    "W68", "I41", "250", "X11", "T15", "T14", "I33", "705"
])

const UNAM_CYAN = RGB(6/255  , 151/255, 213/255)
const UNAM_BLUE = RGB(22/255 , 40/255 , 102/255)
const UNAM_GOLD = RGB(220/255, 167/255, 29/255 )

computationtime(x::DateTime, y::DateTime) = @sprintf("%.2f", (y - x).value / 60_000)

printitle(s::AbstractString, d::AbstractString) = println(d ^ length(s),
        '\n', s, '\n', d ^ length(s))

chi(x::OpticalResidual) = sqrt(chi2(x))
logchi(x::OpticalResidual) = log(chi(x))

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

"""
    Thoss26{T} <: AbstractWeightingScheme{T}

Thoss et al. (2026) optical astrometry weighting scheme
for 3I/ATLAS.

!!! reference
    See:
    - https://doi.org/10.48550/arXiv.2603.15735.
"""
mutable struct Thoss26{T} <: AbstractWeightingScheme{T}
    weights::Vector{NTuple{2, T}}
    corr::Vector{T}
end

# Constructor
function Thoss26(optical::AbstractOpticalVector{T}) where {T <: Real}
    weights = w8sthoss26(optical)
    corrs = corr.(optical)
    return Thoss26{T}(weights, corrs)
end

# Override weights
weights(x::Thoss26) = x.weights

# Override corr
corr(x::Thoss26) = x.corr

# Override getid
getid(::Thoss26) = "Thoss et al. (2026)"

# Override update!
function update!(x::Thoss26{T}, optical::AbstractOpticalVector{T}) where {T <: Real}
    x.weights = w8sthoss26(optical)
    x.corr = corr.(optical)
    return nothing
end

function σsthoss26(obs::AbstractOpticalAstrometry{T}) where {T <: Real}
    σα, σδ = rms(obs)
    code = observatory(obs).code
    σ0 = code in ("250", "X11", "T15", "T14", "I33", "705") ? 0.1 : 1.0
    return (max(σα, σ0), max(σδ, σ0))
end

function w8sthoss26(optical::AbstractOpticalVector{T}) where {T <: Real}
    σs = σsthoss26.(optical)
    rex = rexveres17(optical)
    return @. tuple(1 / (rex * first(σs)), 1 / (rex * last(σs)))
end

function meantime(x::ODProblem)
    t = Vector{Float64}(undef, noptical(x))
    w = Vector{Float64}(undef, noptical(x))
    for i in eachindex(x.optical)
        t[i] = dtutc2days(x.optical[i])
        δ = dec(x.optical[i])
        σα, σδ = x.weights.weights[i]
        w[i] = 1 / (σα^2 * cos(δ)^2 + σδ^2)
    end
    return mean(t, weights(w))
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
    j0, jf = indexin(@view(od.optical[[begin, end]]), OD.optical)
    mask = @. chi(res) < cmax
    ja, jb = findfirst(mask), findlast(mask)
    return min(j0, ja):max(jf, jb), cmax
end

function main(; niter::Int = 20)

    # Parse arguments from commandline
    parsed_args = parse_commandline()

    # Print header
    printitle("Orbit determination for 3I/ATLAS", "=")

    # Number of workers and threads
    println("• Detected 1 worker with ", Threads.nthreads(), " thread(s)")

    # I/O directory
    directory::String = parsed_args["directory"]
    println("• I/O directory: ", directory)

    # Global initial time
    global_initial_time = now()
    println("• Run started at ", global_initial_time)

    # Fetch optical astrometry
    optical = fetch_optical_ades("3I", MPC)
    # Get high-fidelity observatories
    df = DataFrame(optical)
    df.veres = σsveres17.(optical)
    gdf = groupby(df, :stn)
    cdf = combine(gdf, nrow, [:rmsra, :rmsdec, :veres] .=> skipnanmean, renamecols = false)
    sort!(cdf, [:veres, :rmsra, :rmsdec])
    i = findlast(<(1), cdf.veres)
    stations = cdf[1:i, :].stn
    union!(stations, HIGH_FIDELITY_OBSERVATORIES)
    obscodes = Set(getfield.(stations, :code))
    # Note: observations from Catalina Sky Survey appear to be biased
    pop!(obscodes, "703")
    # Filter observations
    filter!(x -> observatory(x).code in obscodes, optical)
    sort!(optical)

    # Parameters
    params = Parameters(
        maxsteps = 10_000, order = 25, abstol = 1E-20, parse_eqs = true,
        coeffstol = Inf, bwdoffset = 0.05, fwdoffset = 0.05,
        gaussorder = 2, safegauss = false, refscale = :log,
        tsaorder = 2, adamiter = 500, adamQtol = 1E-5, mmovproject = false,
        jtlsorder = 2, jtlsmask = false, jtlsiter = 20, lsiter = 10,
        jtlsproject = false, significance = 0.99, verbose = true,
        outrej = true, χ2_rec = 7.0, χ2_rej = 8.0, fudge = 100.0, max_per = 20.0,
        marsden_coeffs = (0.0, 0.0, 0.0), marsden_scalings = (1E-9, 1E-9, 1E-9)
    )
    # Orbit determination problem
    OD = ODProblem(gravityonly!, optical, weights = Thoss26, debias = Eggl20)

    # Initial orbit determination
    printitle("Iteration 1/$niter (cmax = 0.0)", "*")
    i = findfirst(x -> !isempty(x.disc), optical)
    i = findfirst(x -> i in x.indices, OD.tracklets)
    idxs = indices(OD.tracklets[i])
    od = ODProblem(gravityonly!, optical[idxs], weights = Thoss26, debias = Eggl20)
    orbit = tsaiod(od, params; initcond)

    # Observations included in each iteration
    obs = Vector{Set{OpticalADES{Float64}}}(undef, niter+1)
    obs[1] = Set(od.optical)
    # Main loop
    cmax = 0.0
    for i in 1:niter
        # Orbit determination
        if i > 1
            printitle("Iteration $i/$niter (cmax = $cmax)", "*")
            orbit = orbitdetermination(od, orbit, params)
        end
        iszero(orbit) && break
        # Break condition
        (noptical(orbit) == length(optical)) && break
        # Astrometric residuals wrt OD
        _, _, res = propres(OD, orbit(), epoch(orbit) + PE.J2000, params)
        isempty(res) && break
        # Update orbit determination problem
        idxs, cmax = newobservations(OD, od, res)
        NEOs.update!(od, optical[idxs])
        obs[i+1] = Set(setdiff(od.optical, view(obs, 1:i)...))
    end

    # Fit non-gravitational parameters
    printitle("Final orbit", "*")
    ODNG = ODProblem(nongravs!, optical, weights = Thoss26, debias = Eggl20)
    orbitNG = orbitdetermination(ODNG, orbit, params)

    # Save orbit
    filename = joinpath(directory, "3I.jld2")
    jldsave(filename; orbitNG)
    println("• Saved orbit to: ", filename)

    # Residuals plot
    default(dpi = 300, titlefont = (12), guidefont = (12, :black),
            tickfont = (8, :black), framestyle = :box, yminorgrid = false,
            xminorgrid = false, legend = true #=margin = 0.25Measures.mm,=#)

    α = @. ra(orbitNG.ores) / wra(orbitNG.ores)
    δ = @. dec(orbitNG.ores) / wdec(orbitNG.ores)
    omask = @. isoutlier(orbitNG.ores)
    layout = @layout[
        a{} _;
        b{0.8w,0.8h} c{};
    ]
    plot(layout = layout, link = :both, size = (500, 500), margin = -8Plots.px)
    histogram!(α[.!omask], bins = -2:0.1:2, subplot = 1, color = UNAM_CYAN,
        label = "", linewidth = 0.5, xlim = (-2, 2), xticks = -2:0.5:2, ylim = (0, 200),
        framestyle = :box, grid = true, xformatter = Returns(""), yformatter = Returns(""))
    scatter!(
        α[.!omask], δ[.!omask], subplot = 2, color = UNAM_CYAN, label = "Included",
        markersize = 3, markerstrokewidth = 0,
        xlabel = "R.A. Residual [arcsec]", ylabel = "Dec. Residual [arcsec]",
        xlim = (-2, 2), ylim = (-2, 2), aspect_ratio = 1,
        xticks = -2:0.5:2, yticks = -2:0.5:2
    )
    scatter!(α[omask], δ[omask], subplot = 2, color = :red, label = "Rejected",
        markersize = 3, markerstrokewidth = 0)
    histogram!(δ[.!omask], bins = -2:0.1:2, subplot = 3, color = UNAM_CYAN,
        orientation = :horizontal, label = "", linewidth = 0.5, ylim = (-2, 2),
        yticks = -2:0.5:2, xlim = (0, 300), framestyle = :box, grid = true,
        xformatter = Returns(""), yformatter = Returns(""))

    # Save plot
    filename = joinpath(directory, "3I.png")
    savefig(filename)
    println("• Saved plot to: ", filename)

    # Final time
    global_final_time = now()
    println("• Run started ", global_initial_time, " and finished ", global_final_time)
    global_computation_time = computationtime(global_initial_time, global_final_time)
    println("• Total computation time was: ", global_computation_time, " min")

    return nothing
end

main()