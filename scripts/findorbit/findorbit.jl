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
            help = "comma-separated 1-based observation indices/ranges to exclude from the fit but keep in the saved orbit residuals, e.g. 1,4,10-12"
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

function parseexcludedindices(spec::AbstractString, n::Int)
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
        throw(ArgumentError("--exclude-fit indices must be between 1 and $n"))
    elseif length(idxs) == n && n > 0
        throw(ArgumentError("--exclude-fit cannot exclude every observation"))
    end
    return idxs
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
    orbitFR = noptical(orbit) < length(optical) ?
              linkage(od, orbit, paramsFR; maxiter = 20) :
              jtls(od, orbit, paramsFR)
    if iszero(orbitFR)
        println("• $label final refinement failed; keeping previous orbit")
        return orbit
    elseif noptical(orbitFR) < noptical(orbit)
        println("• $label final refinement fit fewer observations than the seed; keeping previous orbit")
        return orbit
    elseif isrmsregression(orbit, orbitFR)
        printrmsregression("$label final refinement", orbit, orbitFR)
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
    bwd, fwd, ores = propres(od, orbit(), epoch(orbit) + PE.J2000, params)
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

function elefield!(line::Vector{Char}, range::UnitRange{Int}, value::AbstractString;
                   align::Symbol = :left)
    width = length(range)
    s = length(value) > width ? value[1:width] :
        align === :right ? lpad(value, width) : rpad(value, width)
    for (i, c) in enumerate(s)
        line[first(range) + i - 1] = c
    end
    return line
end

function elefield!(line::Vector{Char}, start0::Int, stop0::Int, value::AbstractString;
                   align::Symbol = :left)
    return elefield!(line, start0 + 1:stop0, value; align)
end

packeddatecode(x::Integer) = x <= 9 ? string(x) : string(Char('A' + x - 10))

function packedcentury(year::Integer)
    century = year ÷ 100
    return string(Char('I' + century - 18))
end

function packeddate5(dt::DateTime)
    return string(packedcentury(Dates.year(dt)),
                  lpad(mod(Dates.year(dt), 100), 2, '0'),
                  packeddatecode(Dates.month(dt)),
                  packeddatecode(Dates.day(dt)))
end

function dayfraction(dt::DateTime)
    midnight = DateTime(Dates.Date(dt))
    return (dt - midnight).value / 86_400_000
end

function packeddate12(dt::DateTime)
    frac = round(Int, dayfraction(dt) * 10_000_000)
    if frac == 10_000_000
        dt += Day(1)
        frac = 0
    end
    return string(packedcentury(Dates.year(dt)),
                  lpad(mod(Dates.year(dt), 100), 2, '0'),
                  packeddatecode(Dates.month(dt)),
                  packeddatecode(Dates.day(dt)),
                  lpad(frac, 7, '0'))
end

packeddate5(mjd::Real) = packeddate5(julian2datetime(mjd + 2_400_000.5))
packeddate12(mjd::Real) = packeddate12(julian2datetime(mjd + 2_400_000.5))

function monthabbr(dt::DateTime)
    months = ("Jan.", "Feb.", "Mar.", "Apr.", "May ", "June",
              "July", "Aug.", "Sept.", "Oct.", "Nov.", "Dec.")
    return months[Dates.month(dt)]
end

function daydecimal(dt::DateTime)
    return Dates.day(dt) + dayfraction(dt)
end

function impliedfield(x::Real, scale::Real, width::Int)
    n = round(Int, x * scale)
    return lpad(n, width, '0')
end

anglefield(x::Real) = impliedfield(mod(x, 360.0), 1e7, 10)
eccfield(x::Real) = impliedfield(x, 1e9, 10)
qfield(x::Real) = x < 10 ? impliedfield(x, 1e9, 10) : @sprintf("%10.7f", x)

function nonemptyproperty(x, names::Symbol...)
    for name in names
        hasproperty(x, name) || continue
        value = strip(string(getproperty(x, name)))
        isempty(value) || return value
    end
    return ""
end

function opticalids(x::AbstractOpticalAstrometry)
    number = nonemptyproperty(x, :number)
    desig = nonemptyproperty(x, :desig)
    if length(number) == 1 && isletter(only(number)) && !isempty(desig)
        return "", string(number, desig)
    end
    permid = nonemptyproperty(x, :permid)
    provid = nonemptyproperty(x, :provid)
    trksub = nonemptyproperty(x, :trksub)
    permanent = !isempty(number) ? number : !isempty(permid) ? maybe_packnum(permid) : ""
    provisional = !isempty(desig) ? desig : !isempty(provid) ? maybe_packdesig(provid) :
                  trksub
    return permanent, provisional
end

function maybe_packnum(s::AbstractString)
    try
        return packnum(s)
    catch
        return s
    end
end

function maybe_packdesig(s::AbstractString)
    try
        return packdesig(s)
    catch
        return s
    end
end

function packeddesignation(x::AbstractOpticalAstrometry)
    permanent, provisional = opticalids(x)
    isempty(permanent) || return permanent
    isempty(provisional) || return provisional
    return "UNKNOWN"
end

packeddesignation(orbit::LeastSquaresOrbit) = packeddesignation(first(orbit.optical))

function setele255designation!(line::Vector{Char}, orbit::LeastSquaresOrbit)
    permanent, provisional = opticalids(first(orbit.optical))
    if !isempty(permanent)
        elefield!(line, 0, 5, permanent; align = :right)
    elseif !isempty(provisional)
        elefield!(line, 4, 12, provisional)
    end
    return line
end

function ele255name(orbit::LeastSquaresOrbit)
    id = packeddesignation(orbit)
    if length(id) == 8 && isletter(id[1])
        try
            return unpackdesig(id[2:end])
        catch
            return id
        end
    elseif length(id) == 7
        try
            return unpackdesig(id)
        catch
            return id
        end
    end
    return id
end

function orbitcoefficient(orbit::LeastSquaresOrbit, variable::Int)
    j = findfirst(==(variable), variables(orbit))
    return isnothing(j) ? 0.0 : orbit()[j]
end

hasorbitvariable(orbit::LeastSquaresOrbit, variable::Int) =
    !isnothing(findfirst(==(variable), variables(orbit)))

function usedopticaldates(orbit::LeastSquaresOrbit)
    idxs = [i for i in eachindex(orbit.optical) if !isoutlier(orbit.ores[i])]
    isempty(idxs) && return minmaxdates(orbit.optical)
    return extrema(date(orbit.optical[i]) for i in idxs)
end

function cometngfield(x::Real, width::Int, precision::Int)
    return @sprintf("%+*.*f", width, precision, x * 1e8)
end

function mpecele255(orbit::LeastSquaresOrbit, kep, H::Real, params::Parameters)
    line = fill(' ', 255)
    setele255designation!(line, orbit)
    tp_mjd = timeperipass(kep)
    tp_dt = julian2datetime(tp_mjd + 2_400_000.5)
    firstobs, lastobs = usedopticaldates(orbit)
    A2 = orbitcoefficient(orbit, 7)
    A1 = orbitcoefficient(orbit, 8)
    A3 = orbitcoefficient(orbit, 9)
    isng = orbit.dynamics === nongravs! || hasorbitvariable(orbit, 7) ||
           hasorbitvariable(orbit, 8) || hasorbitvariable(orbit, 9)

    elefield!(line, 12, 24, packeddate12(tp_mjd))
    elefield!(line, 24, 34, anglefield(argperi(kep)))
    elefield!(line, 34, 44, anglefield(longascnode(kep)))
    elefield!(line, 44, 54, anglefield(inclination(kep)))
    elefield!(line, 54, 64, qfield(pericenter(kep)))
    elefield!(line, 64, 74, eccfield(eccentricity(kep)))
    elefield!(line, 75, 80, packeddate5(epoch(kep)))
    elefield!(line, 80, 109, ele255name(orbit))
    elefield!(line, 109, 120, "NEOs.jl")
    elefield!(line, 121, 125, string(Dates.year(now())))
    isfinite(H) && elefield!(line, 130, 134, @sprintf("%4.1f", H); align = :right)
    elefield!(line, 135, 139, "10")
    elefield!(line, 148, 150, lpad(mod(Dates.year(tp_dt), 100), 2, '0'))
    elefield!(line, 151, 156, monthabbr(tp_dt))
    elefield!(line, 156, 163, @sprintf("%7.4f", daydecimal(tp_dt)); align = :right)
    elefield!(line, 165, 174, "NEOs.jl")
    elefield!(line, 177, 181, string(min(notout(orbit.ores), 9999)); align = :right)
    isng && elefield!(line, 181, 182, "*")
    elefield!(line, 183, 188, packeddate5(firstobs))
    elefield!(line, 189, 194, packeddate5(lastobs))
    elefield!(line, 195, 196, "h")
    elefield!(line, 197, 200, "M-v")
    elefield!(line, 201, 205, "003E")
    elefield!(line, 206, 207, "9")

    if isng
        elefield!(line, 221, 222, "2")
        hasorbitvariable(orbit, 8) && elefield!(line, 224, 231, cometngfield(A1, 7, 4))
        hasorbitvariable(orbit, 7) && elefield!(line, 233, 242, cometngfield(A2, 9, 6))
        hasorbitvariable(orbit, 9) && elefield!(line, 244, 251, cometngfield(A3, 7, 4))
    elseif iselliptic(kep)
        recip = @sprintf("%9.6f", 1 / semimajoraxis(kep))
        elefield!(line, 230, 239, recip; align = :right)
        elefield!(line, 250, 251, "9")
    end
    elefield!(line, 251, 255, @sprintf("%4.1f", min(NEOs.opticalrms(orbit), 99.9));
              align = :right)

    return String(line)
end

function printele255(ele255::AbstractString)
    printitle("MPC ele255 orbit line", "*")
    println(ele255)
    println("length = ", length(ele255), " characters")
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
    fit_excluded = parseexcludedindices(exclude_fit_spec, length(optical))
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
    ele255 = mpecele255(orbit, kep, H, params)
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
    use_cometary_nongravs && printele255(ele255)

    # Save orbit
    if isempty(fit_excluded)
        jldsave(output; orbit, ele255)
    elseif full_residuals_attached
        heldout_optical = orbit.optical[fit_excluded]
        heldout_ores = orbit.ores[fit_excluded]
        jldsave(output; orbit, fit_excluded, heldout_optical, heldout_ores, ele255)
    else
        heldout_optical = optical[fit_excluded]
        heldout_ores = typeof(orbit.ores)()
        jldsave(output; orbit, fit_excluded, heldout_optical, heldout_ores, ele255)
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
