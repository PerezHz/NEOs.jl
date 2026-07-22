using ArgParse
using NEOs, PlanetaryEphemeris, JLD2, Dates, Statistics, Printf
using NEOs: AbstractOpticalAstrometry, AbstractOpticalVector, OpticalADES,
            OpticalMPC80, AbstractOrbit, log10chi
import NEOs: indices, numberofdays, noptical

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "findorbit.jl"
    # Desciption (for help screen)
    s.description = "Find an orbit from a set of optical astrometry"

    s.epilog = """
        Example:\n
        \n
        julia -t 5 --project findorbit.jl -i 2024YR4\n
        julia -t 5 --project findorbit.jl -i obs.obs -o obs.jld2\n
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
        "--epoch", "-e"
            help = "solution epoch as a Julian Date in the TDB time scale"
            arg_type = Float64
        "--gap-days", "-g"
            help = "maximum gap in days used to split observations into apparitions"
            arg_type = Int
            default = 15
        "--nongravs"
            help = "use nongravs! with A2 for the final refinement"
            action = :store_true
        "--cometary-nongravs"
            help = "use nongravs! with A1, A2 and A3 plus the standard Marsden water-ice radial law for the final refinement"
            action = :store_true
        "--exclude-fit", "-x"
            help = "comma-separated 1-based input file line numbers/ranges to exclude from the fit but keep in the saved orbit residuals, e.g. 1,4,10-12"
            arg_type = String
            default = ""
    end

    return parse_args(s)
end

const SingleApparitionOrbit{O <: AbstractOpticalVector{Float64}} =
    LeastSquaresOrbit{typeof(newtonian!), Float64, Float64, O, Nothing, Nothing}

const MultipleApparitionOrbit{O <: AbstractOpticalVector{Float64}} =
    LeastSquaresOrbit{typeof(gravityonly!), Float64, Float64, O, Nothing, Nothing}

const A2_NONGRAV_SCALINGS = (1E-14, 0.0, 0.0)
const COMETARY_NONGRAV_SCALINGS = (1E-8, 1E-8, 1E-8)
const COMETARY_MARSDEN_RADIAL = (0.111262, 2.808, 2.15, 5.093, 4.6142)
const MAX_RMS_REGRESSION_FACTOR = 1.5
const MAX_RMS_REGRESSION_ARCSEC = 0.5

function fetch_astrometry_format(format::AbstractString)
    fmt = lowercase(strip(format))
    if fmt in ("auto", "ades", "xml")
        return :ades
    elseif fmt in ("mpc80", "obs80")
        return :mpc80
    else
        throw(ArgumentError("Unknown input format: $format. Use auto, ades, mpc80, or obs80."))
    end
end

fetch_optical_astrometry(input::AbstractString, ::Val{:ades}) =
    fetch_optical_ades(input, MPC)

fetch_optical_astrometry(input::AbstractString, ::Val{:mpc80}) =
    fetch_optical_mpc80(input, MPC)

astrometry_format(::AbstractVector{<:OpticalADES}) = "ades"
astrometry_format(::AbstractVector{<:OpticalMPC80}) = "mpc80"

function looks_like_astrometry_file(input::AbstractString)
    return !isempty(splitext(input)[2]) || occursin("/", input) || occursin("\\", input)
end

function load_optical_astrometry(input::AbstractString, format::AbstractString)
    if isfile(input)
        optical = read_optical_astrometry(input; format)
    elseif looks_like_astrometry_file(input)
        throw(ArgumentError("Input astrometry file not found: $input"))
    else
        fmt = fetch_astrometry_format(format)
        optical = fetch_optical_astrometry(input, Val(fmt))
    end
    return optical, astrometry_format(optical)
end

function parseexcludedindices(spec::AbstractString, n::Int;
                              unit::AbstractString = "indices")
    idxs = Int[]
    for rawtoken in split(strip(spec), ',')
        token = replace(strip(rawtoken), ':' => '-')
        isempty(token) && continue
        bounds = split(token, '-')
        if length(bounds) == 1
            push!(idxs, parse(Int, strip(only(bounds))))
        elseif length(bounds) == 2
            firstidx, lastidx = parse.(Int, strip.(bounds))
            firstidx <= lastidx || throw(ArgumentError(
                "Invalid --exclude-fit range: $token"
            ))
            append!(idxs, firstidx:lastidx)
        else
            throw(ArgumentError("Invalid --exclude-fit token: $token"))
        end
    end
    sort!(unique!(idxs))
    if any(i -> i < 1 || i > n, idxs)
        throw(ArgumentError("--exclude-fit $unit must be between 1 and $n"))
    elseif length(idxs) == n && n > 0
        throw(ArgumentError("--exclude-fit cannot exclude every observation"))
    end
    return idxs
end

function ismpc80twolinerline(line::AbstractString)
    length(line) >= 15 && line[15] == 'S' && return true
    length(line) >= 80 || return false
    obscode = view(line, 78:80)
    return obscode == "247" || obscode == "270"
end

function mpc80filerecords(filename::AbstractString)
    lines = readlines(filename)
    records = Dict{Int, String}()
    i = 1
    while i <= length(lines)
        line = lines[i]
        if isempty(strip(line))
            i += 1
        elseif ismpc80twolinerline(line)
            i < length(lines) || throw(ArgumentError(string(
                "MPC80 two-line observation starts at input file line $i, ",
                "but the second line is missing"
            )))
            record = string(line, '\n', lines[i+1], '\n')
            records[i] = record
            records[i+1] = record
            i += 2
        else
            records[i] = string(line, '\n')
            i += 1
        end
    end
    return records, length(lines)
end

function parsempc80record(record::AbstractString, linenumber::Int)
    optical = NEOs.parse_optical_mpc80(record)
    length(optical) == 1 || throw(ArgumentError(
        "Input file line $linenumber does not map to exactly one MPC80 observation"
    ))
    return only(optical)
end

function lineexcludedindices(spec::AbstractString, input::AbstractString,
                             optical::AbstractOpticalVector, format::AbstractString)
    isempty(strip(spec)) && return Int[]
    isfile(input) || throw(ArgumentError(
        "--exclude-fit uses input file line numbers and requires a local astrometry file"
    ))
    format == "mpc80" || throw(ArgumentError(
        "--exclude-fit file line numbers are currently supported for MPC80/OBS80 files only"
    ))
    records, nlines = mpc80filerecords(input)
    linenumbers = parseexcludedindices(spec, nlines; unit = "file line numbers")
    badlines = [line for line in linenumbers if !haskey(records, line)]
    isempty(badlines) || throw(ArgumentError(string(
        "--exclude-fit line(s) do not contain MPC80 optical astrometry: ",
        join(badlines, ",")
    )))
    requested = [line => parsempc80record(records[line], line) for line in linenumbers]
    excluded = Int[]
    ignored = Int[]
    for (line, obs) in requested
        idx = findfirst(==(obs), optical)
        if isnothing(idx)
            push!(ignored, line)
        else
            push!(excluded, idx)
        end
    end
    sort!(unique!(excluded))
    if length(excluded) == length(optical) && !isempty(optical)
        throw(ArgumentError("--exclude-fit cannot exclude every observation"))
    elseif !isempty(ignored)
        println("• Ignored ", length(ignored), " --exclude-fit file line(s) that ",
                "are not present after parsing/filtering: ", join(ignored, ","))
    end
    return excluded
end

function includedindices(n::Int, excluded::AbstractVector{Int})
    excludedset = Set(excluded)
    return [i for i in 1:n if !(i in excludedset)]
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

apparitionrank(x::Apparition) = (noptical(x) >= 5, numberofdays(x), noptical(x))

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

function isrmsregression(previous::LeastSquaresOrbit, candidate::LeastSquaresOrbit)
    previous_rms, candidate_rms = NEOs.opticalrms(previous), NEOs.opticalrms(candidate)
    return isfinite(previous_rms) && isfinite(candidate_rms) &&
           candidate_rms > MAX_RMS_REGRESSION_ARCSEC &&
           candidate_rms > MAX_RMS_REGRESSION_FACTOR * previous_rms
end

function printrmsregression(label::AbstractString, previous::LeastSquaresOrbit,
                            candidate::LeastSquaresOrbit)
    @printf("• %s rejected: RMS would increase from %.6f to %.6f arcsec\n",
            label, NEOs.opticalrms(previous), NEOs.opticalrms(candidate))
    return nothing
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
    sort!(apps, by = apparitionrank, rev = true)
    od = ODProblem(newtonian!, NEOs.optical(apps[1]), weights = Veres17,
                   debias = Eggl20)
    for app in apps
        NEOs.update!(od, NEOs.optical(app))
        for i in 1:2
            orbitSA = i == 1 ? gaussiod(od, params) : tsaiod(od, params; initcond)
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
    i = something(findfirst(>(0), mags), length(apps))
    params = Parameters(params; outrej = false)
    od = ODProblem(newtonian!, NEOs.optical(view(apps, 1:i)), weights = Veres17,
                   debias = Eggl20)
    orbitMID = linkage(od, orbitSA, params)
    iszero(orbitMID) && return orbitSA
    # Step #2: JTLS with gravityonly!
    NEOs.update!(OD, od.optical)
    orbitMA = jtls(OD, orbitMID, params)
    iszero(orbitMA) && return orbitMID
    # Step #3: Outlier rejection
    params = Parameters(params; outrej = true, χ2_rec = sqrt(9.21), χ2_rej = sqrt(10),
                        fudge = 100.0, max_per = 33.3)
    orbitRJ = jtls(OD, orbitMA, params)
    return iszero(orbitRJ) ? orbitMA : orbitRJ
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
        orbit = linkage(OD, orbitMA, params)
        iszero(orbit) && break
        if isrmsregression(orbitMA, orbit)
            printrmsregression("Apparition linkage", orbitMA, orbit)
            break
        end
        orbitMA = orbit
        # Break condition
        # isodvalid(OD, orbitMA, params) && break
        noptical(orbitMA) == noptical(apps) && break
    end
    return orbitMA
end

function finalrefinement(dynamics, orbit::LeastSquaresOrbit,
                         optical::AbstractOpticalVector, params::Parameters;
                         marsden_scalings = params.marsden_scalings)
    label = dynamics === nongravs! ? "Non-gravitational" : "Gravity-only"
    scalings = dynamics === nongravs! ? marsden_scalings : params.marsden_scalings
    paramsFR = Parameters(params; marsden_scalings = scalings,
                          jtlsproject = false, outrej = false)
    od = ODProblem(dynamics, optical, weights = Veres17, debias = Eggl20)
    seed = withoutoutliers(orbit)
    orbitFR = noptical(seed) < length(optical) ?
              linkage(od, seed, paramsFR; maxiter = 20) :
              jtls(od, seed, paramsFR)
    if iszero(orbitFR)
        println("• $label final refinement failed; keeping previous orbit")
        return orbit
    elseif noptical(orbitFR) < noptical(orbit)
        println("• $label final refinement fit fewer observations than the seed; keeping previous orbit")
        return orbit
    end
    candidate = orbitFR
    paramsFR = Parameters(paramsFR; outrej = params.outrej)
    orbitFR = jtls(od, orbitFR, paramsFR)
    if !iszero(orbitFR) && isfinite(nrms(orbitFR)) &&
       (noptical(orbitFR) > noptical(candidate) ||
        (noptical(orbitFR) == noptical(candidate) && nrms(orbitFR) < nrms(candidate)))
        candidate = orbitFR
    end
    if isrmsregression(orbit, candidate)
        printrmsregression("$label final refinement", orbit, candidate)
        return orbit
    end
    if noptical(candidate) == length(optical)
        println("• $label final refinement accepted")
    else
        println("• $label final refinement accepted with ", noptical(candidate),
                " of ", length(optical), " observations")
    end
    return candidate
end

function withoutlier(res::NEOs.OpticalResidual{T, U}, outlier::Bool) where {T, U}
    return NEOs.OpticalResidual{T, U}(
        ra(res), dec(res), wra(res), wdec(res), dra(res), ddec(res), corr(res), outlier
    )
end

function withoutoutliers(orbit::LeastSquaresOrbit)
    return LeastSquaresOrbit(
        orbit.dynamics, orbit.variables, orbit.optical, orbit.tracklets, orbit.radar,
        orbit.bwd, orbit.fwd, [withoutlier(res, false) for res in orbit.ores],
        orbit.rres, orbit.fit, orbit.qs, orbit.Qs
    )
end

function fitoutlierflags(orbit::LeastSquaresOrbit, optical::AbstractOpticalVector)
    flags = trues(length(optical))
    used = falses(length(orbit.optical))
    for i in eachindex(optical)
        j = findfirst(eachindex(orbit.optical)) do k
            return !used[k] && optical[i] == orbit.optical[k]
        end
        if !isnothing(j)
            used[j] = true
            flags[i] = isoutlier(orbit.ores[j])
        end
    end
    return flags
end

function attachfullresiduals(orbit::LeastSquaresOrbit, optical::AbstractOpticalVector,
                             excluded::AbstractVector{Int}, params::Parameters)
    od = ODProblem(orbit.dynamics, optical, weights = Veres17, debias = Eggl20)
    q0 = NEOs.initialcondition(orbit, NEOs.dof(od), params)
    bwd, fwd, ores = propres(od, q0, epoch(orbit) + PE.J2000, params)
    length(ores) == length(optical) || error("Could not compute residuals for full input arc")
    flags = fitoutlierflags(orbit, optical)
    flags[excluded] .= true
    for i in eachindex(ores)
        ores[i] = withoutlier(ores[i], flags[i])
    end
    return LeastSquaresOrbit(
        orbit.dynamics, orbit.variables, optical, reduce_tracklets(optical), orbit.radar,
        bwd, fwd, ores, orbit.rres, orbit.fit, orbit.qs, orbit.Qs
    )
end

function residualchi(res::NEOs.OpticalResidual)
    ρ = corr(res)
    abs(ρ) >= 1 && return NaN
    ξα, ξδ = ra(res), dec(res)
    return sqrt((ξα^2 + ξδ^2 - 2ρ * ξα * ξδ) / (1 - ρ^2))
end

function printexcludedresiduals(optical::AbstractOpticalVector,
                                ores::AbstractVector{<:NEOs.OpticalResidual},
                                excluded::AbstractVector{Int})
    isempty(excluded) && return nothing
    printitle("Held-out residuals", "*")
    println(rpad("idx", 8), rpad("date", 28), rpad("RA", 16),
            rpad("Dec", 16), "chi")
    for i in excluded
        α, δ = residualarcsec(ores[i])
        @printf("%-8d%-28s%+15.5f %+15.5f %12.5f\n",
                i, string(date(optical[i])), α, δ, residualchi(ores[i]))
    end
    println("")
    return nothing
end

function residualarcsec(res::NEOs.OpticalResidual)
    return ra(res) / wra(res), dec(res) / wdec(res)
end

function mpeclikeresidualvalue(x::Real)
    if !isfinite(x)
        return "  NaN"
    end
    signchar = x > 0 ? '+' : x < 0 ? '-' : ' '
    ax = abs(x)
    if ax < 60
        return @sprintf("%4.1f%c", ax, signchar)
    else
        dx = ax / 3600
        return dx < 10 ? @sprintf("%4.2f%c", dx, signchar) :
                         @sprintf("%5.2f%c", dx, signchar)
    end
end

function mpeclikedatecode(x::AbstractOpticalAstrometry)
    return Dates.format(date(x), dateformat"yymmdd")
end

function mpeclikeresidualentry(x::AbstractOpticalAstrometry, res::NEOs.OpticalResidual)
    α, δ = residualarcsec(res)
    αs, δs = mpeclikeresidualvalue(α), mpeclikeresidualvalue(δ)
    values = isoutlier(res) ? "($αs  $δs)" : " $αs  $δs "
    return string(rpad(mpeclikedatecode(x), 6), "  ",
                  rpad(string(observatorycode(x)), 3), " ", values)
end

function printmpeclikeresiduals(orbit::LeastSquaresOrbit)
    printitle("Residuals in seconds of arc", "*")
    perm = sortperm(collect(eachindex(orbit.optical)), by = i -> date(orbit.optical[i]))
    entries = [mpeclikeresidualentry(orbit.optical[i], orbit.ores[i]) for i in perm]
    nrows = ceil(Int, length(entries) / 3)
    for i in 1:nrows
        row = String[]
        for j in i:nrows:length(entries)
            push!(row, entries[j])
        end
        println(join(rpad.(row, 31)))
    end
    println("")
    return nothing
end

function orbitapptype(apps::AbstractApparitionVector, params::Parameters)
    napps = length(apps)
    napps > 0 || throw(ArgumentError("At least one apparition is required"))
    return orbitapptype(Val(min(napps, 3)), apps, params)
end

function orbitapptype(::Val{1}, apps::AbstractApparitionVector, params::Parameters)
    return singleapparition(apps, params)
end

function orbitapptype(::Val{2}, apps::AbstractApparitionVector, params::Parameters)
    orbitSA = singleapparition(apps, params)
    return bridge(apps, orbitSA, params)
end

function orbitapptype(::Val{3}, apps::AbstractApparitionVector, params::Parameters)
    orbitSA = singleapparition(apps, params)
    orbitMA = bridge(apps, orbitSA, params)
    orbitMA isa MultipleApparitionOrbit || return orbitMA
    return multipleapparition(apps, orbitMA, params)
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

    # Solution epoch
    solution_epoch = parsed_args["epoch"]
    if isnothing(solution_epoch)
        println("• Requested solution epoch: weighted mean epoch of observations")
    else
        println("• Requested solution epoch: ",
                @sprintf("%.12f", solution_epoch), " JDTDB")
    end

    # Apparition split gap
    gap_days::Int = parsed_args["gap-days"]
    gap_days > 0 || throw(ArgumentError("Apparition split gap must be positive"))
    println("• Apparition split gap: ", gap_days, " days")

    # Final refinement dynamical model
    use_nongravs::Bool = parsed_args["nongravs"]
    use_cometary_nongravs::Bool = parsed_args["cometary-nongravs"]
    if use_nongravs && use_cometary_nongravs
        throw(ArgumentError("--nongravs and --cometary-nongravs are mutually exclusive"))
    end
    final_dynamics = (use_nongravs || use_cometary_nongravs) ? nongravs! : gravityonly!
    nongrav_scalings = use_cometary_nongravs ? COMETARY_NONGRAV_SCALINGS :
                       use_nongravs ? A2_NONGRAV_SCALINGS : (0.0, 0.0, 0.0)
    nongrav_label = use_cometary_nongravs ? "cometary A1/A2/A3" :
                    use_nongravs ? "A2 only" : "none"
    println("• Final refinement dynamical model: ", final_dynamics,
            " (nongravs: ", nongrav_label, ")")

    # Observations excluded from fit
    exclude_fit_spec::String = parsed_args["exclude-fit"]
    if isempty(strip(exclude_fit_spec))
        println("• Observations excluded from fit: none")
    else
        println("• Requested observations excluded from fit: ", exclude_fit_spec)
    end

    # Global initial time
    global_initial_time = now()
    println("• Run started at ", global_initial_time)

    # Load optical astrometry
    optical, format = load_optical_astrometry(input, format)
    println("• Loaded ", length(optical), " ", uppercase(format), " optical observations")
    filter!(!isdeprecated, optical)
    sort!(optical)
    fit_excluded = lineexcludedindices(exclude_fit_spec, input, optical, format)
    fit_optical = isempty(fit_excluded) ? optical :
                  optical[includedindices(length(optical), fit_excluded)]
    if !isempty(fit_excluded)
        println("• Excluding ", length(fit_excluded),
                " observation(s) from the fit; residuals will be kept in output")
    end
    arc_days = numberofdays(optical)
    split_gap_days = arc_days <= 2 * gap_days ? max(gap_days, ceil(Int, arc_days) + 1) :
                     gap_days
    if split_gap_days != gap_days
        println("• Effective apparition split gap: ", split_gap_days,
                " days (single short arc)")
    end

    # Parameters
    params = Parameters(
        maxsteps = 20_000, order = 15, abstol = 1E-12, parse_eqs = true,
        coeffstol = Inf, bwdoffset = 0.2, fwdoffset = 0.2,
        marsden_radial = use_cometary_nongravs ? COMETARY_MARSDEN_RADIAL :
                         (1.0, 1.0, 2.0, 0.0, 0.0),
        gaussorder = 2, safegauss = false, refscale = :log,
        tsaorder = 2, adamiter = 500, adamQtol = 1E-5,
        jtlsorder = 2, jtlsmask = false, jtlsiter = 20, lsiter = 10,
        jtlsproject = true, significance = 0.99, verbose = true,
        outrej = true, χ2_rec = 4.0, χ2_rej = 5.0, fudge = 100.0,
        max_per = 33.3
    )

    # Split observational arc into apparitions
    apps = apparitions(fit_optical, Day(split_gap_days))
    # Compute orbit by apparition type (single/multiple apparition)
    orbit = orbitapptype(apps, params)
    orbit = finalrefinement(final_dynamics, orbit, fit_optical, params;
                            marsden_scalings = nongrav_scalings)

    # Shift epoch to requested epoch, or to the middle of the observational arc
    jdsolution = isnothing(solution_epoch) ? meanepoch(orbit) + PE.J2000 : solution_epoch
    try
        orbit = shiftepoch(orbit, jdsolution, params; beyondarc = true)
    catch err
        println("• Could not shift final orbit epoch; keeping orbit at fitted epoch (",
                typeof(err), ")")
    end
    # Compute elements and H before attaching held-out observations because H currently
    # uses every optical observation and does not inspect residual outlier flags.
    kep = keplerian(orbit, params)
    H, dH = absolutemagnitude(orbit, params)
    full_residuals_attached = false
    if !isempty(fit_excluded)
        try
            orbit = attachfullresiduals(orbit, optical, fit_excluded, params)
            full_residuals_attached = true
        catch err
            println("• Could not compute residuals for held-out observations; ",
                    "saving fit-only orbit (", typeof(err), ")")
        end
    end
    printitle("Final orbit", "*")
    println(summary(orbit))

    # Print heliocentric ecliptic Keplerian elements (plus q, tp)
    printitle("Keplerian elements", "*")
    println(kep)
    println("q  = ", @sprintf("%+.12E", pericenter(kep)), " au")
    println("tp = ", @sprintf("%+.12E", timeperipass(kep)), " MJD TDB")
    println("H  = ", @sprintf("%.3f", H), " +/- ", @sprintf("%.3f", dH), " mag")
    println("")
    printmpeclikeresiduals(orbit)

    # Save orbit
    if isempty(fit_excluded)
        jldsave(output; orbit)
    elseif full_residuals_attached
        heldout_optical = orbit.optical[fit_excluded]
        heldout_ores = orbit.ores[fit_excluded]
        jldsave(output; orbit, fit_excluded, heldout_optical, heldout_ores)
    else
        heldout_optical = optical[fit_excluded]
        heldout_ores = typeof(orbit.ores)()
        jldsave(output; orbit, fit_excluded, heldout_optical, heldout_ores)
    end
    println("Final orbit saved to: ", output)

    # Final time
    global_final_time = now()
    println("• Run started ", global_initial_time, " and finished ", global_final_time)
    global_computation_time = computationtime(global_initial_time, global_final_time)
    println("• Total computation time was: ", global_computation_time, " min")

    return nothing
end

main()
