using ArgParse, NEOs, Dates, DelimitedFiles, PlanetaryEphemeris, JLD2
using Plots, Colors, Printf
using NEOs: OpticalADES, RadarJPL, AbstractOrbit, rexveres17, isoccultation

ENV["GKSwstype"] = "100"

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "apophis.jl"
    # Desciption (for help screen)
    s.description = "Orbit determination for Apophis"

    @add_arg_table! s begin
        "--directory", "-d"
            help = "I/O directory"
            arg_type = String
    end

    s.epilog = """
        Example:\n
        \n
        julia -t 5 --project apophis.jl -d data/\n
        \n
    """

    return parse_args(s)
end

const ODOrbit = LeastSquaresOrbit{typeof(gravityonly!), Float64, Float64,
        Vector{OpticalADES{Float64}}, Nothing, Nothing}

const UNAM_CYAN = RGB(6/255  , 151/255, 213/255)
const UNAM_BLUE = RGB(22/255 , 40/255 , 102/255)
const UNAM_GOLD = RGB(220/255, 167/255, 29/255 )

computationtime(x::DateTime, y::DateTime) = @sprintf("%.2f", (y - x).value / 60_000)

printitle(s::AbstractString, d::AbstractString) = println(d ^ length(s),
        '\n', s, '\n', d ^ length(s))

istholen13(x::OpticalADES) = catalogue(x).code == 'L' &&
    (x.prog == "0G" && observatory(x).code == "695") ||
    (x.prog == "02" && observatory(x).code == "568")

isapophis(x::OpticalADES) = isempty(x.deprecated) && abs(x.rmscorr) != 1 &&
    !(observatory(x).code == "217" && Date(date(x)) == Date(2021, 1, 28)) &&
    !(observatory(x).code == "L28" && Date(2021, 2, 19) ≤ date(x) ≤ Date(2021, 2, 22))

function initcond(A::AdmissibleRegion)
    v_ρ = sum(A.v_ρ_domain) / 2
    return [
        (A.ρ_domain[1], v_ρ, :log),
        (10^(sum(log10, A.ρ_domain) / 2), v_ρ, :log),
        (sum(A.ρ_domain) / 2, v_ρ, :linear),
        (A.ρ_domain[2], v_ρ, :linear),
    ]
end

function main()

    # Parse arguments from commandline
    parsed_args = parse_commandline()

    # Print header
    printitle("Orbit determination for Apophis", "=")

    # Number of workers and threads
    println("• Detected 1 worker with ", Threads.nthreads(), " thread(s)")

    # I/O directory
    directory::String = parsed_args["directory"]
    println("• I/O directory: ", directory)

    # Global initial time
    global_initial_time = now()
    println("• Run started at ", global_initial_time)

    # Read optical and radar astrometry
    optical = read_optical_ades(joinpath(directory, "99942.xml"))
    radar = read_radar_jpl(joinpath(directory, "99942.json"))
    # Filter out:
    # - deprecated observations
    # - observations with an absolute correlation of 1
    # - biased observations from observatory 217 on January 28, 2021
    # - biased observations from observatory L28 on February 19-21, 2021
    filter!(isapophis, optical)

    # Orbit determination problem
    OD = ODProblem(gravityonly!, optical, radar; weights = Veres17, debias = Eggl20)
    # Custom weights for Tholen et al. (2013) 432 observations
    idxsTholen13 = findall(istholen13, optical)
    rexTholen13 = rexveres17(view(optical, idxsTholen13))
    σsTholen13 = readdlm(joinpath(directory, "tholenetal2013_opterror.dat"), ',')
    OD.weights.weights[idxsTholen13] = @. tuple(
        1 / sqrt((σsTholen13[:,1]^2 + σsTholen13[:,3]^2 + σsTholen13[:,5]^2) * rexTholen13^2),
        1 / sqrt((σsTholen13[:,2]^2 + σsTholen13[:,4]^2 + σsTholen13[:,6]^2) * rexTholen13^2)
    )
    # Custom weights for occultation observations
    idxsOCC = findall(x -> isoccultation(observatory(x)), optical)
    σsOCC = [
        (0.01, 0.01) # (0.0008, 0.0007),
        (0.01, 0.01) # (0.0011, 0.0009),
        (0.01, 0.01) # (0.0004, 0.0005),
        (0.01, 0.01) # (0.0011, 0.0005),
        (0.01, 0.01) # (0.0003, 0.0003),
        # (0.01, 0.01) # (0.0002, 0.0004),
        (0.01, 0.01) # (0.0005, 0.0003)
    ]
    OD.weights.weights[idxsOCC] .= @. tuple(1 / (first(σsOCC)), 1 / (last(σsOCC)))
    # Set correlations in occultation observations
    @. OD.weights.corr[idxsOCC] = corr(optical[idxsOCC])

    # Initial orbit determination with the observations from June 2004
    idxs0 = findall(x -> Date(2004, 5) < date(x) < Date(2004, 7), OD.optical)
    od0 = ODProblem(newtonian!, OD.optical[idxs0]; weights = UniformWeights, debias = Eggl20)
    params = Parameters(
        maxsteps = 2_000, order = 25, abstol = 1E-20, parse_eqs = true,
        coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
        gaussorder = 2, safegauss = false, refscale = :log,
        tsaorder = 2, adamiter = 500, adamQtol = 1E-5,
        jtlsorder = 2, jtlsmask = false, jtlsiter = 20, lsiter = 10,
        significance = 0.99, outrej = false, verbose = true
    )
    orbit0 = initialorbitdetermination(od0, params; initcond)

    # Add the observations in May 2004
    idxs1 = findall(x -> date(x) < Date(2004, 7), OD.optical)
    od1 = ODProblem(newtonian!, OD.optical[idxs1]; weights = UniformWeights, debias = Eggl20)
    for i in 1:6
        od1.weights.weights[i] = od1.weights.weights[i] ./ 50
    end
    orbit1 = orbitdetermination(od1, orbit0, params)
    @. od1.weights.weights = OD.weights.weights[idxs1]
    @. od1.weights.corr = OD.weights.corr[idxs1]
    orbit1 = orbitdetermination(od1, orbit1, params)

    # Add the observations in December 18, 2004
    idxs2 = findall(x -> date(x) < Date(2004, 12, 19), OD.optical)
    od2 = ODProblem(gravityonly!, OD.optical[idxs2]; weights = Veres17, debias = Eggl20)
    @. od2.weights.weights = OD.weights.weights[idxs2]
    @. od2.weights.corr = OD.weights.corr[idxs2]
    orbit2 = orbitdetermination(od2, orbit1, params)

    # Add the remaining observations in December 2004
    idxs3 = findall(x -> date(x) < Date(2005), OD.optical)
    od3 = ODProblem(gravityonly!, OD.optical[idxs3]; weights = Veres17, debias = Eggl20)
    @. od3.weights.weights = OD.weights.weights[idxs3]
    @. od3.weights.corr = OD.weights.corr[idxs3]
    orbit3 = orbitdetermination(od3, orbit2, params)

    # Add all optical and radar observations before February 2005
    idxs4 = findall(x -> date(x) < Date(2005, 2), OD.optical)
    idxs4b = findall(x -> date(x) < Date(2005, 2), OD.radar)
    od4 = ODProblem(gravityonly!, OD.optical[idxs4], OD.radar[idxs4b];
        weights = Veres17, debias = Eggl20)
    @. od4.weights.weights = OD.weights.weights[idxs4]
    @. od4.weights.corr = OD.weights.corr[idxs4]
    orbit4 = orbitdetermination(od4, orbit3, params)

    # Add all optical and radar observations before 2006
    idxs5 = findall(x -> date(x) < Date(2006), OD.optical)
    idxs5b = findall(x -> date(x) < Date(2006), OD.radar)
    od5 = ODProblem(gravityonly!, OD.optical[idxs5], OD.radar[idxs5b];
        weights = Veres17, debias = Eggl20)
    @. od5.weights.weights = OD.weights.weights[idxs5]
    @. od5.weights.corr = OD.weights.corr[idxs5]
    orbit5 = orbitdetermination(od5, orbit4, params)

    # Add all optical and radar observations before 2009
    idxs6 = findall(x -> date(x) < Date(2009), OD.optical)
    idxs6b = findall(x -> date(x) < Date(2009), OD.radar)
    od6 = ODProblem(gravityonly!, OD.optical[idxs6], OD.radar[idxs6b];
        weights = Veres17, debias = Eggl20)
    @. od6.weights.weights = OD.weights.weights[idxs6]
    @. od6.weights.corr = OD.weights.corr[idxs6]
    orbit6 = orbitdetermination(od6, orbit5, params)

    # Add all optical and radar observations before 2014
    idxs7 = findall(x -> date(x) < Date(2014), OD.optical)
    idxs7b = findall(x -> date(x) < Date(2014), OD.radar)
    od7 = ODProblem(gravityonly!, OD.optical[idxs7], OD.radar[idxs7b];
        weights = Veres17, debias = Eggl20)
    @. od7.weights.weights = OD.weights.weights[idxs7]
    @. od7.weights.corr = OD.weights.corr[idxs7]
    orbit7 = orbitdetermination(od7, orbit6, params)

    # Add all optical and radar observations before 2021
    params = Parameters(params; marsden_coeffs = (0.0, 0.0, 0.0),
        marsden_scalings = (1E-14, 0.0, 0.0))
    idxs8 = findall(x -> date(x) < Date(2021), OD.optical)
    idxs8b = findall(x -> date(x) < Date(2021), OD.radar)
    od8 = ODProblem(nongravs!, OD.optical[idxs8], OD.radar[idxs8b];
        weights = Veres17, debias = Eggl20)
    @. od8.weights.weights = OD.weights.weights[idxs8]
    @. od8.weights.corr = OD.weights.corr[idxs8]
    orbit8 = orbitdetermination(od8, orbit7, params)

    # Add all optical and radar observations before 2022
    idxs9 = findall(x -> date(x) < Date(2022), OD.optical)
    idxs9b = findall(x -> date(x) < Date(2022), OD.radar)
    od9 = ODProblem(nongravs!, OD.optical[idxs9], OD.radar[idxs9b];
        weights = Veres17, debias = Eggl20)
    @. od9.weights.weights = OD.weights.weights[idxs9]
    @. od9.weights.corr = OD.weights.corr[idxs9]
    orbit9 = orbitdetermination(od9, orbit8, params)

    # Add the occultation observation from 2022
    od10 = ODProblem(nongravs!, OD.optical, OD.radar;
        weights = Veres17, debias = Eggl20)
    @. od10.weights.weights = OD.weights.weights
    @. od10.weights.corr = OD.weights.corr
    orbit10 = orbitdetermination(od10, orbit9, params)

    # Outlier rejection
    params = Parameters(params; outrej = true, χ2_rec = sqrt(9.21),
        χ2_rej = sqrt(10), fudge = 100.0, max_per = 34.0)
    orbit = orbitdetermination(od10, orbit10, params)

    # Save orbit
    filename = joinpath(directory, "apophis.jld2")
    jldsave(filename; orbit)
    println("• Saved orbit to: ", filename)

    # Residuals plot
    default(dpi = 300, titlefont = (12), guidefont = (12, :black),
            tickfont = (8, :black), framestyle = :box, yminorgrid = false,
            xminorgrid = false, legend = true #=margin = 0.25Measures.mm,=#)

    α = @. ra(orbit.ores) / wra(orbit.ores)
    δ = @. dec(orbit.ores) / wdec(orbit.ores)
    omask = @. isoutlier(orbit.ores)
    layout = @layout[
        a{} _;
        b{0.8w,0.8h} c{};
    ]
    plot(layout = layout, link = :both, size = (500, 500), margin = -8Plots.px)
    histogram!(α[.!omask], bins = -4:0.1:4, subplot = 1, color = UNAM_CYAN,
        label = "", linewidth = 0.5, xlim = (-4, 4), xticks = -4:4, ylim = (0, 1800),
        framestyle = :box, grid = true, xformatter = Returns(""), yformatter = Returns(""))
    scatter!(
        α[.!omask], δ[.!omask], subplot = 2, color = UNAM_CYAN, label = "Included",
        markersize = 3, markerstrokewidth = 0,
        xlabel = "R.A. Residual [arcsec]", ylabel = "Dec. Residual [arcsec]",
        xlim = (-4, 4), ylim = (-4, 4), aspect_ratio = 1,
        xticks = -4:4, yticks = -4:4
    )
    scatter!(α[omask], δ[omask], subplot = 2, color = :red, label = "Rejected",
        markersize = 3, markerstrokewidth = 0)
    histogram!(δ[.!omask], bins = -4:0.1:4, subplot = 3, color = UNAM_CYAN,
        orientation = :horizontal, label = "", linewidth = 0.5, ylim = (-4, 4),
        yticks = -4:4, xlim = (0, 2250), framestyle = :box, grid = true,
        xformatter = Returns(""), yformatter = Returns(""))

    # Save plot
    filename = joinpath(directory, "apophis.png")
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