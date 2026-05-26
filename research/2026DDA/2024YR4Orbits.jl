using ArgParse, Printf
using NEOs, PlanetaryEphemeris, Dates, Statistics, JLD2, JSON
using NEOs: OpticalRWO, parse_optical_mpc80

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "2024YR4Orbits.jl"
    # Desciption (for help screen)
    s.description = "Progressive orbit determination for asteroid 2024 YR4"

    s.epilog = """
        Example:\n
        \n
        julia -t 5 --project 2024YR4Orbits.jl -i 2024YR4/ -o 2024YR4Orbits.jld2\n
        \n
    """

    @add_arg_table! s begin
        "--input", "-i"
            help = "input directory"
            arg_type = String
        "--output", "-o"
            help = "output orbit determination file"
            arg_type = String
    end

    return parse_args(s)
end

const Orbit{D, T} = LeastSquaresOrbit{D, T, T, Vector{OpticalRWO{T}}, Nothing, Nothing}

computationtime(x::DateTime, y::DateTime) = @sprintf("%.2f", (y - x).value / 60_000)

printitle(s::AbstractString, d::AbstractString) = println(d ^ length(s),
    '\n', s, '\n', d ^ length(s))

isodvalid(od::ODProblem, orbit::LeastSquaresOrbit, params::Parameters) =
    noptical(od) == noptical(orbit) && critical_value(orbit) < params.significance &&
    all(!isnan, sigmas(orbit))

function meanepoch(x::Orbit)
    t = Vector{Float64}(undef, noptical(x))
    w = Vector{Float64}(undef, noptical(x))
    for i in eachindex(x.optical)
        t[i] = dtutc2days(x.optical[i])
        δ = dec(x.optical[i])
        σα, σδ = 1 / wra(x.ores[i]), 1 / wdec(x.ores[i])
        w[i] = 1 / (σα^2 * cos(δ)^2 + σδ^2)
    end
    return mean(t, weights(w))
end

function main()
    # Parse arguments from commandline
    parsed_args = parse_commandline()

    # Print header
    printitle("Progressive orbit determination for asteroid 2024 YR4", "=")

    # Number of workers and threads
    println("• Detected 1 worker(s) with ", Threads.nthreads(), " thread(s) each")

    # Input directory
    input::String = parsed_args["input"]
    println("• Input directory: ", input)

    # Output orbit determination file
    output::String = parsed_args["output"]
    println("• Output orbit determination file: ", output)

    # Global initial time
    global_initial_time = now()
    println("• Run started at ", global_initial_time)

    # Load optical astrometry
    filename = joinpath(input, "2024YR4.rwo")
    optical = read_optical_rwo(filename)
    sort!(optical)
    # Load dates of submission to the MPC
    filename = joinpath(input, "2024YR4.json")
    text = read(filename, String)
    dict = JSON.parse(text)
    _optical_ = map(dict) do x
        d = DateTime(x["submission_id"][1:23])
        s = x["obs80"]
        s = length(s) ≤ 80 ? s : string(s[1:80], '\n', s[81:end])
        o = first(parse_optical_mpc80(string(s, '\n')))
        return (d, o)
    end
    sort!(_optical_, by = last)
    submission_dates = first.(_optical_)
    # Eliminate observations not used by ESA
    idxs = findall(x -> x.sel_A == 1, optical)
    keepat!(optical, idxs)
    keepat!(submission_dates, idxs)
    # Compute orbit computation dates
    dates = sort!(unique(Date.(submission_dates)))

    # Parameters
    params = Parameters(
        maxsteps = 10_000, order = 15, abstol = 1E-12, parse_eqs = true,
        coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
        gaussorder = 2, safegauss = false, refscale = :log,
        tsaorder = 2, adamiter = 500, adamQtol = 1E-5,
        jtlsorder = 2, jtlsmask = false, jtlsiter = 20, lsiter = 20,
        jtlsproject = false, significance = 0.99, verbose = true,
        outrej = false, χ2_rec = 7.0, χ2_rej = 8.0, fudge = 100.0,
        max_per = 33.3
    )

    # Main loop
    od = ODProblem(gravityonly!, optical, weights = SourceWeights, debias = SourceDebiasing)
    orbits = Vector{Orbit{typeof(gravityonly!), Float64}}(undef, length(dates))
    for (i, d) in enumerate(dates)
        # Iteration initial time
        iteration_initial_time = now()
        printitle("Iteration #$i", "-")
        println("• Started ", iteration_initial_time)
        println("• Current date is ", d)
        # Current orbit determination problem
        idxs = findall(x -> Date(x) ≤ d, submission_dates)
        NEOs.update!(od, optical[idxs])
        # Current orbit
        if isone(i)
            orbits[i] = initialorbitdetermination(od, params)
        else
            orbits[i] = linkage(od, orbits[i-1], params)
        end
        # Shift epoch to the middle of the observational arc
        tmean = meanepoch(orbits[i])
        orbits[i] = shiftepoch(orbits[i], tmean + PE.J2000, params)
        orbits[i] = jtls(od, orbits[i], params)
        N = length(orbits[i].Qs)
        params.verbose && println(
            "* Shiftepoch converged in $N iterations to: \n\n",
            summary(orbits[i])
        )
        # Iteration final time
        iteration_final_time = now()
        println("• Finished ", iteration_final_time)
        iteration_computation_time = computationtime(iteration_initial_time, iteration_final_time)
        println("• Iteration computation time was ", iteration_computation_time, " min\n")
    end

    # Save output
    jldsave(output; dates, orbits)
    println("• Output saved to: ", output)

    # Final time
    global_final_time = now()
    println("• Run started ", global_final_time, " and finished ", global_final_time)
    global_computation_time = computationtime(global_initial_time, global_final_time)
    println("• Total computation time was: ", global_computation_time, " min")

    return nothing
end

main()