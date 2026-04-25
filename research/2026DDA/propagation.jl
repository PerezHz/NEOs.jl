using ArgParse, Dates, NEOs, PlanetaryEphemeris, TaylorSeries, LinearAlgebra, JLD2
using Plots, Colors, Printf

ENV["GKSwstype"] = "100"

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "propagation.jl"
    # Desciption (for help screen)
    s.description = "Propagate JPL's #220 orbit for Apophis"

    @add_arg_table! s begin
        "--directory", "-d"
            help = "I/O directory"
            arg_type = String
    end

    s.epilog = """
        Example:\n
        \n
        julia -t 5 --project propagation.jl -d data/\n
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

function superscriptify(i::Int)
    if i < 0
        return string("⁻", TS.superscriptify(abs(i)))
    else
        return TS.superscriptify(abs(i))
    end
end

function main()
    # Parse arguments from commandline
    parsed_args = parse_commandline()

    # Print header
    printitle("Propagate JPL's #220 orbit for Apophis", "=")

    # Number of workers and threads
    println("• Detected 1 worker with ", Threads.nthreads(), " thread(s)")

    # I/O directory
    directory::String = parsed_args["directory"]
    println("• I/O directory: ", directory)

    # Global initial time
    global_initial_time = now()
    println("• Run started at ", global_initial_time)

    # Parse Horizons file
    filename = joinpath(directory, "horizons_results.txt")
    lines = readlines(filename)
    j0, jf = findall(Base.Fix2(startswith, "\$\$"), lines)
    keepat!(lines, j0+1:jf-1)
    jds = Vector{Float64}(undef, length(lines)÷3)
    jplrv = Matrix{Float64}(undef, 6, length(lines)÷3)
    jplsig = Matrix{Float64}(undef, 6, length(lines)÷3)
    for (i, sublines) in enumerate(Iterators.partition(lines, 3))
       jds[i] = parse(Float64, sublines[1][1:18])
       jplrv[:, i] = parse.(Float64, split(sublines[2][10:end]))
       jplsig[:, i] = parse.(Float64, split(sublines[3][10:end]))
    end

    # Parameters
    params = Parameters(
        order = 25, abstol = 1E-20, parse_eqs = true
    )
    # Time span
    jd0, jdf = jds[1], jds[end]
    nyears = (jdf - jd0) / yr
    # Initial condition
    q00 = jplrv[:, 1] + params.eph_su(jd0 - J2000)
    q00 = vcat(q00, -2.901766720242E-14, 4.999999873689E-13, 0.0)

    # Warmup propagation
    params = Parameters(params; maxsteps = 1)
    fwd = NEOs.propagate(nongravs!, q00, jd0, nyears, params)
    # Propagation
    params = Parameters(params; maxsteps = 10_000)
    fwd = NEOs.propagate(nongravs!, q00, jd0, nyears, params)

    # Save propagation
    filename = joinpath(directory, "propagation.jld2")
    jldsave(filename; fwd)
    println("• Saved propagation to: ", filename)

    # Compute difference between JPL and NEOs
    ts = jds .- J2000
    dx = Matrix{Float64}(undef, 4, length(ts))
    for (i, t) in enumerate(ts)
        jpl_minus_neos = auday2kmsec(jplrv[:, i] - (fwd(t)[1:6] - params.eph_su(t)))
        jplsigma = auday2kmsec(jplsig[:, i])
        dx[1, i] = norm(jpl_minus_neos[1:3])
        dx[2, i] = norm(jplsigma[1:3])
        dx[3, i] = norm(jpl_minus_neos[4:6])
        dx[4, i] = norm(jplsigma[4:6])
    end

    # Plot
    default(dpi = 300, titlefont = (12), guidefont = (12, :black),
            tickfont = (8, :black), framestyle = :box, yminorgrid = false,
            xminorgrid = false, legend = true #=margin = 0.25Measures.mm,=#)

    dates = @. julian2datetime(ts + J2000)
    xticksv = DateTime(2000):Year(10):DateTime(2100)
    xtickss = Dates.format.(xticksv, "yyyy")
    yticksv = -10:2:6
    ytickss = @. string("10", superscriptify(yticksv))
    begin
        plot(
            xlabel = "Year [TDB]", ylabel = "Absolute difference",
            xticks = (xticksv, xtickss), yticks = (yticksv, ytickss),
            framestyle = :box, dpi = 300
        )
        plot!(dates, log10.(dx[1, :]), linewidth = 1.5, color = UNAM_CYAN,
            linestyle = :solid, label = "Position difference [km]")
        plot!(dates, log10.(dx[2, :]), linewidth = 1.5, color = UNAM_BLUE,
            linestyle = :solid, label = "Position uncertainty [km]")
        plot!(dates, log10.(dx[3, :]), linewidth = 1.5, color = :gold,
            linestyle = :solid, label = "Velocity difference [km/s]")
        plot!(dates, log10.(dx[4, :]), linewidth = 1.5, color = UNAM_GOLD,
            linestyle = :solid, label = "Velocity uncertainty [km/s]")

        # Save plot
        filename = joinpath(directory, "propagation.png")
        savefig(filename)
        println("• Saved plot to: ", filename)
    end

    # Final time
    global_final_time = now()
    println("• Run started ", global_initial_time, " and finished ", global_final_time)
    global_computation_time = computationtime(global_initial_time, global_final_time)
    println("• Total computation time was: ", global_computation_time, " min")

    return nothing
end

main()