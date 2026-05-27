using ArgParse
using NEOs, PlanetaryEphemeris, JLD2, Dates, Statistics, Printf
using NEOs: AbstractOpticalAstrometry, AbstractOpticalVector, OpticalADES,
            OpticalMPC80, AbstractOrbit, log10chi
import NEOs: indices, numberofdays, noptical

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "findorb.jl"
    # Desciption (for help screen)
    s.description = "Find an orbit from a set of optical astrometry"

    s.epilog = """
        Example:\n
        \n
        julia -t 5 --project findorb.jl -i 2024YR4\n
        \n
    """

    @add_arg_table! s begin
        "--input", "-i"
            help = "input designation/astrometry file"
            arg_type = String
        "--output", "-o"
            help = "output .jld2 file"
            arg_type = String
        "--format", "-f"
            help = "input format: auto, ades, mpc80, or obs80"
            arg_type = String
            default = "auto"
    end

    return parse_args(s)
end

const SingleApparitionOrbit{O <: AbstractOpticalVector{Float64}} =
    LeastSquaresOrbit{typeof(newtonian!), Float64, Float64, O, Nothing, Nothing}

const MultipleApparitionOrbit{O <: AbstractOpticalVector{Float64}} =
    LeastSquaresOrbit{typeof(gravityonly!), Float64, Float64, O, Nothing, Nothing}

function normalize_astrometry_format(format::AbstractString)
    fmt = lowercase(strip(format))
    if fmt in ("auto", "ades", "xml")
        return fmt == "xml" ? "ades" : fmt
    elseif fmt in ("mpc80", "obs80")
        return "mpc80"
    else
        throw(ArgumentError("Unknown input format: $format. Use auto, ades, mpc80, or obs80."))
    end
end

function detect_astrometry_format(filename::AbstractString)
    for line in eachline(filename)
        stripped = strip(line)
        isempty(stripped) && continue
        return startswith(stripped, '<') ? "ades" : "mpc80"
    end
    throw(ArgumentError("Cannot detect astrometry format from empty file: $filename"))
end

function load_optical_astrometry(input::AbstractString, format::AbstractString)
    fmt = normalize_astrometry_format(format)
    if isfile(input)
        fmt = fmt == "auto" ? detect_astrometry_format(input) : fmt
        optical = fmt == "ades" ? read_optical_ades(input) : read_optical_mpc80(input)
    else
        fmt = fmt == "auto" ? "ades" : fmt
        optical = fmt == "ades" ? fetch_optical_ades(input, MPC) : fetch_optical_mpc80(input, MPC)
    end
    return optical, fmt
end

struct Apparition{T <: Real, O <: AbstractOpticalAstrometry{T}, V <: AbstractVector{O},
                  I <: AbstractVector{Int}, B}
    optical::SubArray{O, 1, V, Tuple{I}, B}
end

const AbstractApparitionVector{T} = AbstractVector{Apparition{T, O, V, I, B}} where {O, V, I, B}

indices(x::Apparition) = first(x.optical.indices)
NEOs.optical(x::Apparition) = collect(x.optical)
NEOs.optical(x::AbstractApparitionVector) = sort!(mapreduce(NEOs.optical, vcat, x))
numberofdays(x::Apparition) = numberofdays(x.optical)
noptical(x::Apparition) = length(x.optical)
noptical(x::AbstractApparitionVector) = sum(noptical, x)

function apparitions(optical::AbstractOpticalVector{T},
                     gap::Period = Day(30)) where {T <: Real}
    sort!(optical)
    apps = [[1]]
    for i in 2:length(optical)
        if date(optical[i]) - date(optical[i-1]) > gap
            push!(apps, [i])
        else
            push!(apps[end], i)
        end
    end
    return [Apparition(view(optical, i)) for i in apps]
end

computationtime(x::DateTime, y::DateTime) = @sprintf("%.2f", (y - x).value / 60_000)

printitle(s::AbstractString, d::AbstractString) = println(d ^ length(s),
          '\n', s, '\n', d ^ length(s))

isodvalid(od::ODProblem, orbit::LeastSquaresOrbit, params::Parameters) =
        noptical(od) == noptical(orbit) && critical_value(orbit) < params.significance &&
        all(!isnan, sigmas(orbit))

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

function meanepoch(x::AbstractOrbit)
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

function singleapparition(apps::AbstractApparitionVector, params::Parameters)
    # Single apparition orbit determination
    optical = NEOs.optical(first(apps))
    orbitSA = zero(SingleApparitionOrbit{typeof(optical)})
    sort!(apps, by = numberofdays, rev = true)
    od = ODProblem(newtonian!, NEOs.optical(apps[1]), weights = Veres17,
                   debias = Eggl20)
    for i in 1:2
        for app in apps
            NEOs.update!(od, NEOs.optical(app))
            if i == 1
                orbitSA = gaussiod(od, params)
            else
                orbitSA = tsaiod(od, params; initcond)
            end
            isodvalid(od, orbitSA, params) && break
        end
        isodvalid(od, orbitSA, params) && break
    end
    return orbitSA
end

function bridge(apps::AbstractApparitionVector, orbitSA::SingleApparitionOrbit,
                params::Parameters)
    # Bridge between single and multiple apparitions
    OD = ODProblem(gravityonly!, NEOs.optical(apps), weights = Veres17, debias = Eggl20)
    _, _, res = propres(OD, orbitSA(), epoch(orbitSA) + PE.J2000, params)
    mags = Vector{Float64}(undef, length(apps))
    for (i, app) in enumerate(apps)
        if issubset(NEOs.optical(app), orbitSA.optical)
            mags[i] = zero(Float64)
        else
            mags[i] = maximum(log10chi, view(res, indices(app)))
        end
    end
    perm = sortperm(mags)
    permute!(mags, perm)
    permute!(apps, perm)
    # Step #1: Linkage with newtonian!
    i = findfirst(>(0), mags)
    params = Parameters(params; outrej = false)
    od = ODProblem(newtonian!, NEOs.optical(view(apps, 1:i)), weights = Veres17,
                   debias = Eggl20)
    orbitMID = linkage(od, orbitSA, params)
    # Step #2: JTLS with gravityonly!
    NEOs.update!(OD, od.optical)
    orbitMA = jtls(OD, orbitMID, params)
    # Step #3: Outlier rejection
    params = Parameters(params; outrej = true, χ2_rec = sqrt(9.21), χ2_rej = sqrt(10),
                        fudge = 100.0, max_per = 33.3)
    orbitMA = jtls(OD, orbitMA, params)
    return orbitMA
end

function multipleapparition(apps::AbstractApparitionVector, orbitMA::MultipleApparitionOrbit,
                            params::Parameters)
    # Multiple apparition orbit determination
    OD = ODProblem(gravityonly!, NEOs.optical(apps), weights = Veres17, debias = Eggl20)
    _, _, res = propres(OD, orbitMA(), epoch(orbitMA) + PE.J2000, params)
    mags = Vector{Float64}(undef, length(apps))
    for (i, app) in enumerate(apps)
        if issubset(NEOs.optical(app), orbitMA.optical)
            mags[i] = zero(Float64)
        else
            mags[i] = maximum(log10chi, view(res, indices(app)))
        end
    end
    perm = sortperm(mags)
    permute!(mags, perm)
    permute!(apps, perm)
    for i in eachindex(mags)
        iszero(mags[i]) && continue
        NEOs.update!(OD, NEOs.optical(view(apps, 1:i)))
        orbitMA = linkage(OD, orbitMA, params)
        # Break condition
        # isodvalid(OD, orbitMA, params) && break
        noptical(orbitMA) == noptical(apps) && break
    end
    return orbitMA
end

function main()
    # Parse arguments from commandline
    parsed_args = parse_commandline()

    # Print header
    printitle("Find an orbit from a set of optical astrometry", "=")

    # Number of workers and threads
    println("• Detected 1 worker with ", Threads.nthreads(), " thread(s)")

    # Input designation/astrometry file
    input::String = parsed_args["input"]
    println("• Input designation/astrometry file: ", input)

    # Output .jld2 file
    output::String = parsed_args["output"]
    println("• Output .jld2 file: ", output)

    # Input astrometry format
    format::String = parsed_args["format"]
    println("• Requested input astrometry format: ", format)

    # Global initial time
    global_initial_time = now()
    println("• Run started at ", global_initial_time)

    # Load optical astrometry
    optical, format = load_optical_astrometry(input, format)
    println("• Loaded ", length(optical), " ", uppercase(format), " optical observations")
    filter!(!isdeprecated, optical)

    # Parameters
    params = Parameters(
        maxsteps = 20_000, order = 15, abstol = 1E-12, parse_eqs = true,
        coeffstol = Inf, bwdoffset = 0.05, fwdoffset = 0.05,
        gaussorder = 2, safegauss = false, refscale = :log,
        tsaorder = 2, adamiter = 500, adamQtol = 1E-5,
        jtlsorder = 2, jtlsmask = false, jtlsiter = 20, lsiter = 10,
        jtlsproject = true, significance = 0.99, verbose = true,
        outrej = true, χ2_rec = 7.0, χ2_rej = 8.0, fudge = 100.0,
        max_per = 33.3
    )

    # Split observational arc into apparitions
    apps = apparitions(optical, Day(238))
    # Single apparition orbit determination
    orbit = singleapparition(apps, params)

    if length(apps) > 1
        # Bridge between single and multiple apparitions
        orbit = bridge(apps, orbit, params)
    end
    if length(apps) > 2
        # Multiple apparition orbit determination
        orbit = multipleapparition(apps, orbit, params)
    end

    # Shift epoch to the middle of the observational arc
    tmean = meanepoch(orbit)
    orbit = shiftepoch(orbit, tmean + PE.J2000, params)
    printitle("Final orbit", "*")
    println(summary(orbit))

    # Print heliocentric ecliptic Keplerian elements (plus q, tp)
    kep = keplerian(orbit, params)
    H, dH = absolutemagnitude(orbit, params)
    printitle("Keplerian elements", "*")
    println(kep)
    println("q  = ", @sprintf("%+.12E", pericenter(kep)), " au")
    println("tp = ", @sprintf("%+.12E", timeperipass(kep)), " MJD TDB")
    println("H  = ", @sprintf("%.3f", H), " +/- ", @sprintf("%.3f", dH), " mag")
    println("")

    # Save orbit
    jldsave(output; orbit)
    println("Final orbit saved to: ", output)

    # Final time
    global_final_time = now()
    println("• Run started ", global_initial_time, " and finished ", global_final_time)
    global_computation_time = computationtime(global_initial_time, global_final_time)
    println("• Total computation time was: ", global_computation_time, " min")

    return nothing
end

main()
