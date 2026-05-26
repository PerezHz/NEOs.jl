using ArgParse, Printf, Plots, Colors, Measures, Statistics
using NEOs, PlanetaryEphemeris, TaylorSeries, Dates, JLD2

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "2024YR4Plots.jl"
    # Desciption (for help screen)
    s.description = "Impact monitoring plots for asteroid 2024 YR4"

    s.epilog = """
        Example:\n
        \n
        julia -t 5 --project 2024YR4Plots.jl -d 2024YR4/\n
        \n
    """

    @add_arg_table! s begin
        "--directory", "-d"
            help = "I/O directory"
            arg_type = String
    end

    return parse_args(s)
end

const UNAM_CYAN = RGB(6/255  , 151/255, 213/255)
const UNAM_BLUE = RGB(22/255 , 40/255 , 102/255)
const UNAM_GOLD = RGB(220/255, 167/255, 29/255 )

computationtime(x::DateTime, y::DateTime) = @sprintf("%.2f", (y - x).value / 60_000)

printitle(s::AbstractString, d::AbstractString) = println(d ^ length(s),
    '\n', s, '\n', d ^ length(s))

cumip(VIs) = sum(impact_probability, VIs; init = 0.0)

function main()
    # Parse arguments from commandline
    parsed_args = parse_commandline()

    # Print header
    printitle("Impact monitoring plots for asteroid 2024 YR4", "=")

    # Number of workers and threads
    println("• Detected 1 worker(s) with ", Threads.nthreads(), " thread(s) each")

    # I/O directory
    directory::String = parsed_args["directory"]
    println("• I/O directory: ", directory)

    # Global initial time
    global_initial_time = now()
    println("• Run started at ", global_initial_time)

    # NEOs.jl

    # Computation dates
    filename = joinpath(directory, "2024YR4Orbits.jld2")
    dates_NEOs = JLD2.load(filename, "dates")

    # Virtual impactors
    filename = joinpath(directory, "2024YR4IMMoon.jld2")
    VIss_Moon = JLD2.load(filename, "VIss")

    filename = joinpath(directory, "2024YR4IMEarth.jld2")
    VIss_Earth = JLD2.load(filename, "VIss")

    # Cumulative impact probabilities
    CIPs_NEOs_Moon = cumip.(VIss_Moon)
    CIPs_NEOs_Earth = cumip.(VIss_Earth)

    # NEOCC

    # Parse file
    filename = joinpath(directory, "NEOCCIM.txt")
    lines = readlines(filename)[2:end]
    # Computation dates, cumulative and maximum impact probabilities
    dates_NEOCC = Vector{Date}(undef, length(lines))
    CIPs_NEOCC = Vector{Float64}(undef, length(lines))
    for (i, l) in enumerate(lines)
        dates_NEOCC[i] = Date(l[1:10])
        CIPs_NEOCC[i] = eval(Meta.parse(l[32:41]))
    end
    perm = sortperm(dates_NEOCC)
    permute!(dates_NEOCC, perm)
    permute!(CIPs_NEOCC, perm)

    # Plots.jl settings
    default(
        dpi = 300, titlefont = (12), guidefont = (12, :black),
        tickfont = (8, :black), framestyle = :box, yminorgrid = false, xminorgrid = false,
        legend = true #=margin = 0.25Measures.mm,=#
    )

    # Plot: Cumulative impact probability vs time
    xticksv = Date(2025, 1):Month(1):Date(2025, 9)
    xtickss = Dates.format.(xticksv, "u/yy")
    xtickss[end] = "Mar/26"
    yticksv = -6:-1
    ytickss = @. string("10⁻", TS.superscriptify(abs(yticksv)))
    x_NEOCC = vcat(dates_NEOCC[1:end-1], Date(2025, 9, 5))
    y_NEOCC = vcat(log10.(CIPs_NEOCC[1:end-1]), -10)
    x_NEOs = vcat(dates_NEOs[1:end-1], Date(2025, 9, 5))
    y_NEOs_Moon = vcat(log10.(CIPs_NEOs_Moon[1:end-1]), -10)
    y_NEOs_Earth = vcat(log10.(CIPs_NEOs_Earth[1:end-1]), -10)
    P = plot(
        xlabel = "Computation time", ylabel = "Impact probability",
        xticks = (xticksv, xtickss), yticks = (yticksv, ytickss),
        xlim = (Date(2024, 12, 15), Date(2025, 9, 15)), ylim = (-6, -1),
        margin = 1.0Plots.mm, legend = (0.07, 0.14), legendfontsize = 8
    )
    plot!(P, x_NEOs, y_NEOs_Moon, label = "", linetype = :steppost,
          color = UNAM_CYAN, linewidth = 1.75)
    scatter!(P, [], [], label = "NEOs.jl (Moon)", color = UNAM_CYAN)
    plot!(P, x_NEOs, y_NEOs_Earth, label = "", linetype = :steppost,
          color = UNAM_BLUE, linewidth = 1.75)
    scatter!(P, [], [], label = "NEOs.jl (Earth)", color = UNAM_BLUE)
    plot!(P, x_NEOCC, y_NEOCC, label = "", linetype = :steppost,
          color = UNAM_GOLD, linewidth = 1.75)
    scatter!(P, [], [], label = "NEOCC (Earth)", color = UNAM_GOLD)
    filename = joinpath(directory, "2024YR4CIP.png")
    savefig(filename)
    println("• Cumulative impact probability plot saved to: ", filename)

    # Final time
    global_final_time = now()
    println("• Run started ", global_final_time, " and finished ", global_final_time)
    global_computation_time = computationtime(global_initial_time, global_final_time)
    println("• Total computation time was: ", global_computation_time, " min")

    return nothing
end

main()