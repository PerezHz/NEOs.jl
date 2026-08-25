using ArgParse
using NEOs, PlanetaryEphemeris, JLD2, Dates, LinearAlgebra, Statistics, Printf
using NEOs: AbstractOpticalAstrometry, AbstractOpticalVector, AbstractApparitionVector,
            OpticalADES, OpticalMPC80, AbstractOrbit, log10chi, indices

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
        "--no-mags-res"
            help = "do not print magnitude residuals in the final residual table"
            action = :store_true
        "--no-track-res"
            help = "do not print along-track and cross-track residuals in the final residual table"
            action = :store_true
        "--force-all-fit"
            help = "keep the final refinement that fits every non-manually-excluded observation, even if its RMS is worse than the linked seed"
            action = :store_true
        "--max-rms-jump-factor"
            help = "maximum allowed RMS increase factor when accepting a linked/refined orbit"
            arg_type = Float64
            default = 1.5
        "--max-rms-jump-arcsec"
            help = "minimum candidate RMS in arcsec before the RMS jump guard rejects a linked/refined orbit"
            arg_type = Float64
            default = 0.5
        "--ephem-days"
            help = "post-fit ephemeris horizon in days after the latest observation; 0 disables it"
            arg_type = Float64
            default = 0.0
        "--ephem-step"
            help = "post-fit ephemeris step in days"
            arg_type = Float64
            default = 1.0
        "--ephem-obscode"
            help = "MPC observatory code for the post-fit ephemeris"
            arg_type = String
            default = "500"
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
const JUPITER_SEMIMAJOR_AXIS_AU = 5.20336301
const C3_AUDAY2_TO_KMSEC2 = (au / daysec)^2
const MAX_RMS_REGRESSION_FACTOR = 1.5
const MAX_RMS_REGRESSION_ARCSEC = 0.5
const MAX_ACCEPTABLE_NRMS = 10.0
const MOON_EPH = NEOs.selecteph(NEOs.sseph, NEOs.mo)

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

apparitionrank(x) = (noptical(x) >= 5, noptical(x), numberofdays(x))

computationtime(x::DateTime, y::DateTime) = @sprintf("%.2f", (y - x).value / 60_000)

printitle(s::AbstractString, d::AbstractString) = println(d ^ length(s),
          '\n', s, '\n', d ^ length(s))

isodvalid(od::ODProblem, orbit::LeastSquaresOrbit, params::Parameters) =
        noptical(od) == noptical(orbit) && critical_value(orbit) < params.significance &&
        all(!isnan, sigmas(orbit))

function isrmsregression(previous::LeastSquaresOrbit, candidate::LeastSquaresOrbit;
                         max_factor::Real = MAX_RMS_REGRESSION_FACTOR,
                         max_arcsec::Real = MAX_RMS_REGRESSION_ARCSEC)
    previous_rms, candidate_rms = NEOs.opticalrms(previous), NEOs.opticalrms(candidate)
    return isfinite(previous_rms) && isfinite(candidate_rms) &&
           candidate_rms > max_arcsec &&
           candidate_rms > max_factor * previous_rms
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

function singleapparitionseeds(apps::AbstractApparitionVector, params::Parameters)
    # Single apparition orbit determination
    optical = NEOs.optical(first(apps))
    orbitSA = zero(SingleApparitionOrbit{typeof(optical)})
    seeds = SingleApparitionOrbit[]
    sort!(apps, by = apparitionrank, rev = true)
    od = ODProblem(newtonian!, NEOs.optical(apps[1]), weights = Veres17,
                   debias = Eggl20)
    for app in apps
        NEOs.update!(od, NEOs.optical(app))
        for i in 1:2
            orbitSA = i == 1 ? gaussiod(od, params) : tsaiod(od, params; initcond)
            isodvalid(od, orbitSA, params) && break
        end
        isodvalid(od, orbitSA, params) && push!(seeds, orbitSA)
    end
    return seeds
end

function singleapparition(apps::AbstractApparitionVector, params::Parameters)
    optical = NEOs.optical(first(apps))
    seeds = singleapparitionseeds(apps, params)
    isempty(seeds) && return zero(SingleApparitionOrbit{typeof(optical)})
    orbitSA = first(seeds)
    return orbitSA
end

function bridge(apps::AbstractApparitionVector, orbitSA::SingleApparitionOrbit,
                params::Parameters;
                max_rms_jump_factor::Real = MAX_RMS_REGRESSION_FACTOR,
                max_rms_jump_arcsec::Real = MAX_RMS_REGRESSION_ARCSEC)
    # Bridge between single and multiple apparitions
    OD = ODProblem(gravityonly!, NEOs.optical(apps), weights = Veres17, debias = Eggl20)
    _, _, res = propres(OD, orbitSA(), epoch(orbitSA) + PE.J2000, params)
    length(res) == length(OD.optical) || return orbitSA
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
    idxs = isnothing(i) ? eachindex(apps) : 1:i
    params = Parameters(params; outrej = false)
    od = ODProblem(newtonian!, NEOs.optical(view(apps, idxs)), weights = Veres17,
                   debias = Eggl20)
    orbitMID = linkage(od, orbitSA, params)
    iszero(orbitMID) && return orbitSA
    if isrmsregression(orbitSA, orbitMID; max_factor = max_rms_jump_factor,
                       max_arcsec = max_rms_jump_arcsec)
        printrmsregression("Initial apparition linkage", orbitSA, orbitMID)
        return orbitSA
    elseif isfitloss(orbitSA, orbitMID)
        printfitloss("Initial apparition linkage", orbitSA, orbitMID)
        return orbitSA
    end
    # Step #2: JTLS with gravityonly!
    NEOs.update!(OD, od.optical)
    orbitMA = jtls(OD, orbitMID, params)
    iszero(orbitMA) && return orbitMID
    if isrmsregression(orbitMID, orbitMA; max_factor = max_rms_jump_factor,
                       max_arcsec = max_rms_jump_arcsec)
        printrmsregression("Initial gravity-only refinement", orbitMID, orbitMA)
        return orbitMID
    elseif isfitloss(orbitMID, orbitMA)
        printfitloss("Initial gravity-only refinement", orbitMID, orbitMA)
        return orbitMID
    end
    # Step #3: Outlier rejection
    params = Parameters(params; outrej = true, χ2_rec = 7.0, χ2_rej = 8.0,
                        fudge = 100.0, max_per = 33.3)
    orbitRJ = jtls(OD, orbitMA, params)
    return iszero(orbitRJ) ? orbitMA : orbitRJ
end

function multipleapparition(apps::AbstractApparitionVector, orbitMA::MultipleApparitionOrbit,
                            params::Parameters;
                            max_rms_jump_factor::Real = MAX_RMS_REGRESSION_FACTOR,
                            max_rms_jump_arcsec::Real = MAX_RMS_REGRESSION_ARCSEC)
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
    accepted = [apps[i] for i in eachindex(apps) if iszero(mags[i])]
    # Keep accepted inliers stable while growing the arc; final refinement
    # performs global outlier rejection after all linkable apparitions are present.
    linkage_params = Parameters(params; outrej = false)
    for i in eachindex(mags)
        iszero(mags[i]) && continue
        trial_apps = push!(copy(accepted), apps[i])
        NEOs.update!(OD, NEOs.optical(trial_apps))
        orbit = linkage(OD, orbitMA, linkage_params)
        iszero(orbit) && continue
        if isrmsregression(orbitMA, orbit; max_factor = max_rms_jump_factor,
                           max_arcsec = max_rms_jump_arcsec)
            printrmsregression("Apparition linkage", orbitMA, orbit)
            continue
        elseif isfitloss(orbitMA, orbit)
            printfitloss("Apparition linkage", orbitMA, orbit)
            continue
        end
        push!(accepted, apps[i])
        orbitMA = orbit
        # Break condition
        # isodvalid(OD, orbitMA, params) && break
        noptical(orbitMA) == noptical(apps) && break
    end
    return orbitMA
end

function finalrefinement(dynamics, orbit::LeastSquaresOrbit,
                         optical::AbstractOpticalVector, params::Parameters;
                         marsden_scalings = params.marsden_scalings,
                         force_all_fit::Bool = false,
                         max_rms_jump_factor::Real = MAX_RMS_REGRESSION_FACTOR,
                         max_rms_jump_arcsec::Real = MAX_RMS_REGRESSION_ARCSEC)
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
    if !isfinite(nrms(candidate)) || nrms(candidate) >= MAX_ACCEPTABLE_NRMS
        println("• $label final refinement failed fit-quality check; keeping previous orbit")
        return orbit
    end
    if isrmsregression(orbit, candidate; max_factor = max_rms_jump_factor,
                       max_arcsec = max_rms_jump_arcsec) &&
       !(force_all_fit && noptical(candidate) == length(optical))
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

function nfitoptical(orbit::LeastSquaresOrbit)
    return count(!isoutlier, orbit.ores)
end

function retainedfitoptical(previous::LeastSquaresOrbit, candidate::LeastSquaresOrbit)
    flags = fitoutlierflags(candidate, previous.optical)
    return count(eachindex(previous.optical)) do i
        !isoutlier(previous.ores[i]) && !flags[i]
    end
end

function isfitloss(previous::LeastSquaresOrbit, candidate::LeastSquaresOrbit)
    return retainedfitoptical(previous, candidate) < nfitoptical(previous)
end

function printfitloss(label::AbstractString, previous::LeastSquaresOrbit,
                      candidate::LeastSquaresOrbit)
    lost = nfitoptical(previous) - retainedfitoptical(previous, candidate)
    println("• $label rejected: would discard $lost previously fitted observation(s)")
    return nothing
end

function missingopticalindices(orbit::LeastSquaresOrbit, optical::AbstractOpticalVector)
    used = falses(length(orbit.optical))
    missing = Int[]
    for i in eachindex(optical)
        j = findfirst(eachindex(orbit.optical)) do k
            return !used[k] && optical[i] == orbit.optical[k]
        end
        if isnothing(j)
            push!(missing, i)
        else
            used[j] = true
        end
    end
    return missing
end

function attachfullresiduals(orbit::LeastSquaresOrbit, optical::AbstractOpticalVector,
                             excluded::AbstractVector{Int}, params::Parameters)
    od = ODProblem(orbit.dynamics, optical, weights = Veres17, debias = Eggl20)
    q0 = NEOs.initialcondition(orbit, epoch(orbit) + PE.J2000, NEOs.dof(od), params)
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

dot3(x::AbstractVector, y::AbstractVector) = x[1]*y[1] + x[2]*y[2] + x[3]*y[3]
norm3(x::AbstractVector) = sqrt(dot3(x, x))

function asteroidstatekmsec(orbit::LeastSquaresOrbit, et)
    t = et / PlanetaryEphemeris.daysec
    sol = NEOs.cte(t) <= epoch(orbit) ? orbit.bwd : orbit.fwd
    return PlanetaryEphemeris.auday2kmsec(sol(t))
end

denseephstatekmsec(eph, et) =
    PlanetaryEphemeris.auday2kmsec(eph(et / PlanetaryEphemeris.daysec))

# Differentiate the topocentric light-time RA/Dec model with respect to receive ET.
# Returned rates are arcsec/sec; only their direction is used for along/cross residuals.
function computedradecrate(x::AbstractOpticalAstrometry, orbit::LeastSquaresOrbit,
                           params::Parameters; niter::Int = 5)
    xva = et -> asteroidstatekmsec(orbit, et)
    return computedradecrate(observatory(x), date(x), xva, params; niter)
end

function computedradecrate(obs, date::DateTime, xva, params::Parameters;
                           niter::Int = 5)
    et0 = dtutc2et(date)
    et_r_secs = NEOs.Taylor1([et0, one(et0)], 1)

    rv_s_t_r = denseephstatekmsec(params.eph_su, et_r_secs)
    r_s_t_r = rv_s_t_r[1:3]
    rv_e_t_r = denseephstatekmsec(params.eph_ea, et_r_secs)
    r_e_t_r = rv_e_t_r[1:3]
    rv_a_t_r = xva(et_r_secs)
    r_a_t_r = rv_a_t_r[1:3]
    RV_r = obsposvelECI(obs, et_r_secs)
    R_r = RV_r[1:3]
    r_r_t_r = r_e_t_r + R_r

    ρ_vec_r = r_a_t_r - r_r_t_r
    ρ_r = norm3(ρ_vec_r)
    τ_D = ρ_r / NEOs.clightkms
    et_b_secs = et_r_secs - τ_D
    e_D_vec = r_r_t_r - r_s_t_r
    e_D = norm3(e_D_vec)

    local r_a_t_b, r_s_t_b
    for _ in 1:niter
        rv_a_t_b = xva(et_b_secs)
        r_a_t_b = rv_a_t_b[1:3]
        v_a_t_b = rv_a_t_b[4:6]
        ρ_vec_r = r_a_t_b - r_r_t_r
        ρ_r = norm3(ρ_vec_r)
        rv_s_t_b = denseephstatekmsec(params.eph_su, et_b_secs)
        r_s_t_b = rv_s_t_b[1:3]
        p_D_vec = r_a_t_b - r_s_t_b
        p_D = norm3(p_D_vec)
        Δτ_rel_D = NEOs.shapiro_delay(e_D, p_D, ρ_r)
        p_dot_23 = dot3(ρ_vec_r, v_a_t_b) / ρ_r
        Δt_2 = (τ_D - ρ_r / NEOs.clightkms - Δτ_rel_D) /
               (1.0 - p_dot_23 / NEOs.clightkms)
        τ_D -= Δt_2
        et_b_secs = et_r_secs - τ_D
    end

    rv_a_t_b = xva(et_b_secs)
    r_a_t_b = rv_a_t_b[1:3]
    ρ_vec_r = r_a_t_b - r_r_t_r
    ρ_r = norm3(ρ_vec_r)
    rv_s_t_b = denseephstatekmsec(params.eph_su, et_b_secs)
    r_s_t_b = rv_s_t_b[1:3]

    E_H_vec = r_r_t_r - r_s_t_r
    U_vec = ρ_vec_r
    U_norm = ρ_r
    u_vec = U_vec / U_norm
    Q_vec = r_a_t_b - r_s_t_b
    q_vec = Q_vec / norm3(Q_vec)
    E_H = norm3(E_H_vec)
    e_vec = E_H_vec / E_H
    g1 = NEOs.g1coeff / (E_H / PlanetaryEphemeris.au)
    g2 = 1 + dot3(q_vec, e_vec)
    u1_vec = U_norm * (u_vec + (g1 / g2) *
        (dot3(u_vec, q_vec) * e_vec - dot3(e_vec, u_vec) * q_vec))
    u1_norm = norm3(u1_vec)

    αdot_rad = (u1_vec[1][0] * u1_vec[2][1] - u1_vec[2][0] * u1_vec[1][1]) /
               (u1_vec[1][0]^2 + u1_vec[2][0]^2)
    δ_rad = asin(u1_vec[3] / u1_norm)
    return rad2arcsec(αdot_rad), rad2arcsec(δ_rad[1])
end

function alongcrossresidualarcsec(x::AbstractOpticalAstrometry, res::NEOs.OpticalResidual,
                                  orbit::LeastSquaresOrbit, params::Parameters)
    α, δ = residualarcsec(res)
    αdot, δdot = computedradecrate(x, orbit, params)
    vα = αdot * cos(dec(x))
    vδ = δdot
    speed = hypot(vα, vδ)
    isfinite(speed) && !iszero(speed) || return NaN, NaN
    uα, uδ = vα / speed, vδ / speed
    return α * uα + δ * uδ, -α * uδ + δ * uα
end

function computedradecarcsec(x::AbstractOpticalAstrometry, orbit::LeastSquaresOrbit,
                             params::Parameters)
    return computedradecarcsec(observatory(x), date(x), orbit, params)
end

function computedradecarcsec(obs, date::DateTime, orbit::LeastSquaresOrbit,
                             params::Parameters)
    return compute_radec(obs, date;
        xvs = et -> denseephstatekmsec(params.eph_su, et),
        xve = et -> denseephstatekmsec(params.eph_ea, et),
        xva = et -> asteroidstatekmsec(orbit, et)
    )
end

function tangentoffsetsarcsec(α::AbstractVector{<:Real}, δ::AbstractVector{<:Real},
                              δref::Real)
    ξ = [NEOs.anglediff(a, first(α)) * cos(δref) for a in α]
    η = [d - first(δ) for d in δ]
    return ξ, η
end

function linearslope(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    length(x) >= 2 || return NaN
    x̄, ȳ = mean(x), mean(y)
    dx = x .- x̄
    denom = sum(abs2, dx)
    iszero(denom) && return NaN
    return sum(dx .* (y .- ȳ)) / denom
end

function positionangledeg(vα::Real, vδ::Real)
    speed = hypot(vα, vδ)
    isfinite(speed) && !iszero(speed) || return NaN
    return mod(rad2deg(atan(vα, vδ)), 360)
end

function elevationdeg(obs, date::DateTime, orbit::LeastSquaresOrbit, params::Parameters)
    αas, δas = computedradecarcsec(obs, date, orbit, params)
    return elevationdeg(obs, date, αas, δas)
end

function elevationdeg(obs, date::DateTime, αas::Real, δas::Real)
    α, δ = arcsec2rad(αas), arcsec2rad(δas)
    lineofsight = [cos(δ) * cos(α), cos(δ) * sin(α), sin(δ)]
    zenith = obsposvelECI(obs, dtutc2et(date))[1:3]
    z_norm = norm3(zenith)
    isfinite(z_norm) && !iszero(z_norm) || return NaN
    sin_elevation = clamp(dot3(lineofsight, zenith) / z_norm, -1, 1)
    return rad2deg(asin(sin_elevation))
end

function vectorangledeg(x::AbstractVector, y::AbstractVector)
    nx, ny = norm3(x), norm3(y)
    all(isfinite, (nx, ny)) && !iszero(nx) && !iszero(ny) || return NaN
    return rad2deg(acos(clamp(dot3(x, y) / (nx * ny), -1, 1)))
end

function ephemerisgeometry(t::Real, orbit::LeastSquaresOrbit, params::Parameters)
    return ephemerisgeometry(t, orbit(t)[1:3], params)
end

function ephemerisgeometry(t::Real, xa::AbstractVector, params::Parameters)
    xe = params.eph_ea(t)[1:3]
    xs = params.eph_su(t)[1:3]
    xm = MOON_EPH(t)[1:3]
    geocentric = xa - xe
    heliocentric = xa - xs
    return (
        solar_elongation = vectorangledeg(geocentric, xs - xe),
        lunar_elongation = vectorangledeg(geocentric, xm - xe),
        geocentric_distance = norm3(geocentric),
        heliocentric_distance = norm3(heliocentric),
        phase_angle = vectorangledeg(-geocentric, -heliocentric),
    )
end

function postfitpropagation(orbit::LeastSquaresOrbit, tmin::Real, tmax::Real,
                            params::Parameters; jet::Bool = false)
    t0 = epoch(orbit)
    variables = NEOs.variables(orbit)
    Ndof = NEOs.dof(orbit)
    q0 = NEOs.initialcondition(orbit, epoch(orbit) + PE.J2000, Ndof, params)
    if jet
        T = eltype(q0)
        NEOs.set_od_order(T, 1, length(variables))
        q0 += NEOs.jtperturbation(ones(T, Ndof), variables, Ndof, 1, params)
    end
    margin = 0.1
    tlo = min(tmin, t0) - margin
    thi = max(tmax, t0) + margin
    bwd = NEOs.propagate(NEOs.dynamicalmodel(orbit), q0, t0 + PE.J2000,
                         (tlo - t0) / PlanetaryEphemeris.yr, params)
    fwd = NEOs.propagate(NEOs.dynamicalmodel(orbit), q0, t0 + PE.J2000,
                         (thi - t0) / PlanetaryEphemeris.yr, params)
    return (; epoch = t0, bwd, fwd)
end

function postfitstate(eph, t)
    sol = NEOs.cte(t) <= eph.epoch ? eph.bwd : eph.fwd
    return sol(t)
end

postfitstatekmsec(eph, et) =
    PlanetaryEphemeris.auday2kmsec(postfitstate(eph, et / PlanetaryEphemeris.daysec))

function skyplaneuncertainty(obs, date::DateTime, eph, orbit::LeastSquaresOrbit,
                             params::Parameters)
    αas, δas = compute_radec(obs, date;
        xvs = et -> denseephstatekmsec(params.eph_su, et),
        xve = et -> denseephstatekmsec(params.eph_ea, et),
        xva = et -> postfitstatekmsec(eph, et),
    )
    δ0 = NEOs.cte(δas)
    tangent = [αas * cos(arcsec2rad(δ0)), δas]
    x0 = zeros(length(NEOs.variables(orbit)))
    Γ = Matrix(NEOs.project(tangent, x0, NEOs.covariance(orbit)))
    eig = eigen(Symmetric(Γ))
    order = sortperm(eig.values; rev = true)
    major_index, minor_index = order[1], order[2]
    major = sqrt(max(eig.values[major_index], 0.0))
    minor = sqrt(max(eig.values[minor_index], 0.0))
    major_axis = eig.vectors[:, major_index]
    pa = mod(rad2deg(atan(major_axis[1], major_axis[2])), 180)
    return major, minor, pa
end

function formathms(αas::Real)
    total_ms = mod(round(Int, αas * 1_000 / 15), 86_400_000)
    h, rem_ms = divrem(total_ms, 3_600_000)
    m, rem_ms = divrem(rem_ms, 60_000)
    return @sprintf("%02d:%02d:%06.3f", h, m, rem_ms / 1_000)
end

function formatdms(δas::Real)
    sign = signbit(δas) ? '-' : '+'
    total_ms = round(Int, abs(δas) * 1_000)
    d, rem_ms = divrem(total_ms, 3_600_000)
    m, rem_ms = divrem(rem_ms, 60_000)
    return @sprintf("%c%02d:%02d:%05.2f", sign, d, m, rem_ms / 1_000)
end

function ephemerisnumber(x::Real, width::Int, digits::Int)
    isfinite(x) || return " " ^ width
    value = @sprintf("%.*f", digits, x)
    length(value) <= width && return lpad(value, width)
    return lpad(@sprintf("%.2e", x), width)
end

function ephemerisoffsets(days::Real, step::Real)
    offsets = collect(0.0:step:days)
    isempty(offsets) && push!(offsets, 0.0)
    if days - last(offsets) > 10eps(Float64) * max(1.0, abs(days))
        push!(offsets, Float64(days))
    end
    return offsets
end

function printpostfitephemeris(orbit::LeastSquaresOrbit, params::Parameters, H::Real,
                               days::Real, step::Real, obs)
    days > 0 || return nothing
    _, start = minmaxdates(orbit)
    dates = [start + Millisecond(round(Int, offset * PlanetaryEphemeris.daysec * 1_000))
             for offset in ephemerisoffsets(days, step)]
    times = dtutc2days.(dates)
    scalar_eph = postfitpropagation(orbit, minimum(times), maximum(times), params)
    jet_eph = postfitpropagation(orbit, minimum(times), maximum(times), params; jet = true)
    xva = et -> postfitstatekmsec(scalar_eph, et)

    printitle("Post-fit ephemeris", "*")
    println("Observatory: ", obs.code, " (", obs.name, ")",
            "; times are UTC; RA/Dec are astrometric J2000.")
    println("Rate is arcsec/min; elevation, elongations and PA are degrees; distances are au.")
    println("SigMaj/SigMin are formal linearized 1-sigma sky-plane axes in arcsec; SigPA is east of north.")
    println(
        rpad("Date UTC", 20), "  ", rpad("RA", 12), "  ", rpad("Dec", 12), "  ",
        lpad("Elev", 6), "  ", lpad("Rate", 8), "  ", lpad("PA", 6), "  ",
        lpad("SolEl", 6), "  ", lpad("LunEl", 6), "  ", lpad("GeoDist", 9), "  ",
        lpad("r", 8), "  ", lpad("V", 6), "  ", lpad("SigMaj", 10), "  ",
        lpad("SigMin", 10), "  ", lpad("SigPA", 6),
    )
    for (date, t) in zip(dates, times)
        αas, δas = compute_radec(obs, date;
            xvs = et -> denseephstatekmsec(params.eph_su, et),
            xve = et -> denseephstatekmsec(params.eph_ea, et),
            xva,
        )
        αdot, δdot = computedradecrate(obs, date, xva, params)
        vα = 60 * αdot * cos(arcsec2rad(δas))
        vδ = 60 * δdot
        rate = motionspeed(vα, vδ)
        pa = positionangledeg(vα, vδ)
        geometry = ephemerisgeometry(t, postfitstate(scalar_eph, t)[1:3], params)
        elevation = elevationdeg(obs, date, αas, δas)
        phase = deg2rad(geometry.phase_angle)
        Φ = NEOs.phase_integral(phase, params.slope)
        V = H + 5 * log10(geometry.geocentric_distance * geometry.heliocentric_distance) -
            2.5 * log10(Φ)
        σmajor, σminor, σpa = skyplaneuncertainty(obs, date, jet_eph, orbit, params)
        println(
            rpad(Dates.format(date, "yyyy-mm-ddTHH:MM:SS"), 20), "  ",
            rpad(formathms(αas), 12), "  ", rpad(formatdms(δas), 12), "  ",
            ephemerisnumber(elevation, 6, 1), "  ", ephemerisnumber(rate, 8, 2), "  ",
            ephemerisnumber(pa, 6, 1), "  ",
            ephemerisnumber(geometry.solar_elongation, 6, 1), "  ",
            ephemerisnumber(geometry.lunar_elongation, 6, 1), "  ",
            ephemerisnumber(geometry.geocentric_distance, 9, 6), "  ",
            ephemerisnumber(geometry.heliocentric_distance, 8, 4), "  ",
            ephemerisnumber(V, 6, 2), "  ",
            ephemerisnumber(σmajor, 10, 3), "  ",
            ephemerisnumber(σminor, 10, 3), "  ",
            ephemerisnumber(σpa, 6, 1),
        )
    end
    println("")
    return nothing
end

function alongcrossrate(vαobs::Real, vδobs::Real, vαeph::Real, vδeph::Real)
    all(isfinite, (vαobs, vδobs, vαeph, vδeph)) || return NaN, NaN
    speed = hypot(vαeph, vδeph)
    iszero(speed) && return NaN, NaN
    uα, uδ = vαeph / speed, vδeph / speed
    dvα, dvδ = vαobs - vαeph, vδobs - vδeph
    return dvα * uα + dvδ * uδ, -dvα * uδ + dvδ * uα
end

motionspeed(vα::Real, vδ::Real) = all(isfinite, (vα, vδ)) ? hypot(vα, vδ) : NaN

function motionratevalue(x::Real, width::Int = 11)
    isfinite(x) || return " " ^ width
    return lpad(@sprintf("%+.2f", x), width)
end

function motionanglevalue(x::Real, width::Int = 7)
    isfinite(x) || return " " ^ width
    return lpad(@sprintf("%.1f", x), width)
end

function motionelevationvalue(x::Real, width::Int = 7)
    isfinite(x) || return " " ^ width
    return lpad(@sprintf("%+.1f", x), width)
end

function motiondistancevalue(x::Real, width::Int = 9)
    isfinite(x) || return " " ^ width
    return lpad(@sprintf("%.5f", x), width)
end

function motionarcvalue(x::Real, width::Int = 9)
    isfinite(x) || return " " ^ width
    return lpad(@sprintf("%.2f", x), width)
end

motiontrackletindices(orbit::LeastSquaresOrbit) = indices.(orbit.tracklets)


function trackletmotion(idxs::AbstractVector{Int}, orbit::LeastSquaresOrbit,
                        params::Parameters)
    idxs = sort(collect(idxs), by = i -> date(orbit.optical[i]))
    optical = orbit.optical[idxs]
    dates = date.(optical)
    jds = dtutc2jdtdb.(dates)
    t0 = mean(jds)
    tmin = @. (jds - t0) * 1_440
    span = (maximum(dates) - minimum(dates)).value / 60_000
    δref = mean(dec.(optical))

    αobs = rad2arcsec.(ra.(optical))
    δobs = rad2arcsec.(dec.(optical))
    ξobs, ηobs = tangentoffsetsarcsec(αobs, δobs, δref)
    vαobs = linearslope(tmin, ξobs)
    vδobs = linearslope(tmin, ηobs)

    if length(optical) == 1
        αdot, δdot = computedradecrate(only(optical), orbit, params)
        vαeph, vδeph = 60 * αdot * cos(dec(only(optical))), 60 * δdot
    else
        eph = [computedradecarcsec(x, orbit, params) for x in optical]
        αeph = first.(eph)
        δeph = last.(eph)
        ξeph, ηeph = tangentoffsetsarcsec(αeph, δeph, δref)
        vαeph = linearslope(tmin, ξeph)
        vδeph = linearslope(tmin, ηeph)
    end

    dAT, dCT = alongcrossrate(vαobs, vδobs, vαeph, vδeph)
    geometry = ephemerisgeometry(t0 - PE.J2000, orbit, params)
    return (
        date = jdtdb2dtutc(t0),
        observatory = observatorycode(first(optical)),
        nobs = length(idxs),
        nout = count(isoutlier, orbit.ores[idxs]),
        span = span,
        elevation = elevationdeg(observatory(first(optical)), jdtdb2dtutc(t0),
                                 orbit, params),
        solar_elongation = geometry.solar_elongation,
        lunar_elongation = geometry.lunar_elongation,
        geocentric_distance = geometry.geocentric_distance,
        arcobs = motionspeed(vαobs, vδobs) * span,
        arceph = motionspeed(vαeph, vδeph) * span,
        vαobs = vαobs,
        vδobs = vδobs,
        speedobs = motionspeed(vαobs, vδobs),
        vαeph = vαeph,
        vδeph = vδeph,
        speedeph = motionspeed(vαeph, vδeph),
        dAT = dAT,
        dCT = dCT,
        PAobs = positionangledeg(vαobs, vδobs),
        PAeph = positionangledeg(vαeph, vδeph),
    )
end

function printmotionbytracklet(orbit::LeastSquaresOrbit, params::Parameters)
    printitle("Motion by tracklet", "*")
    println("Rates and dAT/dCT are arcsec/min; POS arclengths are arcsec; elevation, elongations and PA are degrees; GeoDist is au.")
    println("RAcosD is dRA*cosDec; dAT/dCT are observed-minus-ephemeris rates.")
    header = string(
        rpad("Date", 6), "  ", rpad("Obs", 3), "  ", lpad("N", 3), "  ",
        lpad("Out", 3), "  ", lpad("Span", 8), "  ", lpad("Elev", 7), "  ",
        lpad("SolEl", 7), "  ", lpad("LunEl", 7), "  ", lpad("GeoDist", 9), "  ",
        lpad("Arc obs", 9), "  ", lpad("Arc eph", 9), "  ",
        lpad("RAcosD obs", 11), "  ", lpad("Dec obs", 9), "  ", lpad("Speed obs", 9), "  ",
        lpad("RAcosD eph", 11), "  ", lpad("Dec eph", 9), "  ",
        lpad("Speed eph", 9), "  ", lpad("dAT", 9), "  ", lpad("dCT", 9), "  ",
        lpad("PAobs", 7), "  ", lpad("PAeph", 7)
    )
    println(header)
    for idxs in motiontrackletindices(orbit)
        m = trackletmotion(idxs, orbit, params)
        @printf("%-6s  %-3s  %3d  %3d  %8.2f  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s\n",
                Dates.format(m.date, dateformat"yymmdd"), m.observatory, m.nobs, m.nout,
                m.span, motionelevationvalue(m.elevation),
                motionanglevalue(m.solar_elongation), motionanglevalue(m.lunar_elongation),
                motiondistancevalue(m.geocentric_distance), motionarcvalue(m.arcobs),
                motionarcvalue(m.arceph), motionratevalue(m.vαobs),
                motionratevalue(m.vδobs, 9), motionratevalue(m.speedobs, 9), motionratevalue(m.vαeph),
                motionratevalue(m.vδeph, 9), motionratevalue(m.speedeph, 9),
                motionratevalue(m.dAT, 9), motionratevalue(m.dCT, 9),
                motionanglevalue(m.PAobs), motionanglevalue(m.PAeph))
    end
    println("")
    return nothing
end

residualsignchar(x::Real, ax::Real) = iszero(ax) ? ' ' : x > 0 ? '+' : x < 0 ? '-' : ' '

function mpeclikeresidualvalue(x::Real)
    if !isfinite(x)
        return "  NaN"
    end
    ax = abs(x)
    if ax < 60
        rounded = round(ax; digits = 1)
        return @sprintf("%4.1f%c", rounded, residualsignchar(x, rounded))
    else
        rounded = round(ax / 3600; digits = 2)
        return rounded < 10 ? @sprintf("%4.2f%c", rounded, residualsignchar(x, rounded)) :
                              @sprintf("%5.2f%c", rounded, residualsignchar(x, rounded))
    end
end

function magnituderesiduals(orbit::LeastSquaresOrbit, H::Real, params::Parameters)
    res = fill(NaN, length(orbit.optical))
    isfinite(H) || return res
    for i in eachindex(orbit.optical)
        optical = orbit.optical[i]
        observed_v = mag(optical) + vconversion(optical)
        isfinite(observed_v) || continue
        t = dtutc2days(optical)
        xa = orbit(t)[1:3]
        xe = params.eph_ea(t)[1:3]
        xs = params.eph_su(t)[1:3]
        d_BS = norm(xa - xs)
        d_OS = norm(xe - xs)
        d_BO = norm(xa - xe)
        phase = acos(clamp((d_BO^2 + d_BS^2 - d_OS^2) / (2 * d_BO * d_BS), -1, 1))
        phase_integral = NEOs.phase_integral(phase, params.slope)
        isfinite(phase_integral) && phase_integral > 0 || continue
        computed_v = H + 5 * log10(d_BS * d_BO) - 2.5 * log10(phase_integral)
        res[i] = observed_v - computed_v
    end
    return res
end

function mpeclikemagresidualvalue(x::Real)
    isfinite(x) || return "      "
    signchar = x < 0 ? '-' : '+'
    return lpad(@sprintf("%c%.2f", signchar, abs(x)), 6)
end

function mpeclikedatecode(x::AbstractOpticalAstrometry)
    return Dates.format(date(x), dateformat"yymmdd")
end

function mpeclikeresidualheader(include_track::Bool, include_mags::Bool)
    parts = String[lpad("RA", 5), lpad("Dec", 5)]
    if include_track
        append!(parts, (lpad("AT", 5), lpad("CT", 5)))
    end
    if include_mags
        push!(parts, lpad("Mag", 6))
    end
    return string(rpad("Date", 6), "  ", rpad("Obs", 3), "  ", join(parts, "  "))
end

function mpeclikeresidualentry(x::AbstractOpticalAstrometry, res::NEOs.OpticalResidual,
                               radec::Tuple{<:Real, <:Real},
                               atct::Union{Nothing, Tuple{<:Real, <:Real}} = nothing,
                               magres::Union{Nothing, Real} = nothing)
    αs, δs = mpeclikeresidualvalue.(radec)
    parts = String[αs, δs]
    if !isnothing(atct)
        ats, cts = mpeclikeresidualvalue.(atct)
        append!(parts, (ats, cts))
    end
    if !isnothing(magres)
        push!(parts, mpeclikemagresidualvalue(magres))
    end
    values = join(parts, "  ")
    values = isoutlier(res) ? "($values)" : " $values "
    return string(rpad(mpeclikedatecode(x), 6), "  ",
                  rpad(string(observatorycode(x)), 3), " ", values)
end

function printmpeclikeresiduals(orbit::LeastSquaresOrbit, H::Real, params::Parameters;
                                include_mags::Bool = true,
                                include_track::Bool = true)
    title = if include_track
        include_mags ? "Residuals in seconds of arc (RA, Dec, AT, CT) and magnitudes" :
                       "Residuals in seconds of arc (RA, Dec, AT, CT)"
    else
        include_mags ? "Residuals in seconds of arc and magnitudes" :
                       "Residuals in seconds of arc"
    end
    printitle(title, "*")
    perm = sortperm(collect(eachindex(orbit.optical)), by = i -> date(orbit.optical[i]))
    magres = include_mags ? magnituderesiduals(orbit, H, params) : nothing
    radec = [residualarcsec(orbit.ores[i]) for i in perm]
    atct = include_track ?
        [alongcrossresidualarcsec(orbit.optical[i], orbit.ores[i], orbit, params) for i in perm] :
        nothing
    entries = [mpeclikeresidualentry(orbit.optical[i], orbit.ores[i], radec[j],
                                     include_track ? atct[j] : nothing,
                                     include_mags ? magres[i] : nothing)
               for (j, i) in enumerate(perm)]
    if isempty(entries)
        println("")
        return nothing
    end
    entry_width = include_mags ? (include_track ? 54 : 40) : (include_track ? 45 : 31)
    nrows = ceil(Int, length(entries) / 3)
    ncols = ceil(Int, length(entries) / nrows)
    header = mpeclikeresidualheader(include_track, include_mags)
    println(join(rpad.(fill(header, ncols), entry_width)))
    for i in 1:nrows
        row = String[]
        for j in i:nrows:length(entries)
            push!(row, entries[j])
        end
        println(join(rpad.(row, entry_width)))
    end
    println("")
    return nothing
end

function pystringliteral(s::AbstractString)
    return string('"', replace(s, "\\" => "\\\\", "\"" => "\\\""), '"')
end

function pyele220arg(x::AbstractString)
    return pystringliteral(x)
end

function pyele220arg(x::Integer)
    return string(x)
end

function pyele220arg(x::Real)
    isfinite(x) && return @sprintf("%.15g", x)
    isnan(x) && return "float(\"nan\")"
    return x > 0 ? "float(\"inf\")" : "float(\"-inf\")"
end

function ele220inputname(input::AbstractString, format::AbstractString,
                         fallback::AbstractString)
    if isfile(input) && format == "mpc80"
        line = rpad(readline(input), 12)
        return strip(line[6:12])
    end
    return fallback
end

function notoutoptical(orbit::LeastSquaresOrbit)
    return [orbit.optical[i] for i in eachindex(orbit.optical) if !isoutlier(orbit.ores[i])]
end

function ele220argumentlist(name::AbstractString, orbit::LeastSquaresOrbit, kep, H::Real,
                            params::Parameters, gap_days::Int)
    optical = notoutoptical(orbit)
    isempty(optical) && throw(ArgumentError("Cannot build ele220 arguments with zero fitted observations"))
    tbegin, tend = extrema(date.(optical))
    args = Any[
        name,
        timeperipass(kep) + 2_400_000.5,
        mod(meananomaly(kep), 360),
        mod(argperi(kep), 360),
        mod(longascnode(kep), 360),
        inclination(kep),
        pericenter(kep),
        eccentricity(kep),
        H,
        NEOs.opticalrms(orbit),
        length(optical),
        length(apparitions(collect(optical), Day(gap_days))),
        dtutc2jdtdb(tbegin),
        dtutc2jdtdb(tend),
        epoch(kep) + 2_400_000.5,
        ishyperbolic(kep) ? 9 : uncertaintyparameter(orbit, params),
    ]
    return join(pyele220arg.(args), ", ")
end

function printele220arguments(name::AbstractString, orbit::LeastSquaresOrbit, kep, H::Real,
                              params::Parameters, gap_days::Int)
    printitle("MPC ele220 arguments", "*")
    println(ele220argumentlist(name, orbit, kep, H, params, gap_days))
    println("")
    return nothing
end

function propagatedsigma(kep, gradient::AbstractVector{<:Real})
    variance = dot(gradient, NEOs.covariance(kep) * gradient)
    if variance < 0
        return abs(variance) <= 100 * eps(typeof(float(variance))) ? 0.0 : NaN
    end
    return sqrt(variance)
end

function pericentersigma(kep)
    gradient = zeros(length(elements(kep)))
    if iselliptic(kep)
        a, e = semimajoraxis(kep), eccentricity(kep)
        gradient[1] = 1 - e
        gradient[2] = -a
    else
        gradient[1] = 1
    end
    return propagatedsigma(kep, gradient)
end

function timeperipasssigma(kep)
    gradient = zeros(length(elements(kep)))
    if iselliptic(kep)
        a, M, n = semimajoraxis(kep), meananomaly(kep), meanmotion(kep)
        gradient[1] = -3 * M / (2 * a * n)
        gradient[6] = -1 / n
    else
        gradient[6] = 1
    end
    return propagatedsigma(kep, gradient)
end

function tisserandjupiter(kep; a_jupiter::Real = JUPITER_SEMIMAJOR_AXIS_AU)
    a, e, inc = semimajoraxis(kep), eccentricity(kep), inclination(kep)
    radicand = (a / a_jupiter) * (1 - e^2)
    radicand < 0 && return NaN
    return a_jupiter / a + 2 * cosd(inc) * sqrt(radicand)
end

function tisserandjupitersigma(kep; a_jupiter::Real = JUPITER_SEMIMAJOR_AXIS_AU)
    a, e, inc = semimajoraxis(kep), eccentricity(kep), inclination(kep)
    radicand = (a / a_jupiter) * (1 - e^2)
    radicand <= 0 && return NaN
    root = sqrt(radicand)
    dT_da = -a_jupiter / a^2 + cosd(inc) * (1 - e^2) / (a_jupiter * root)
    dT_de = -2 * a * e * cosd(inc) / (a_jupiter * root)
    dT_di = -2 * root * sind(inc) * (pi / 180)
    gradient = zeros(length(elements(kep)))
    if iselliptic(kep)
        gradient[1] = dT_da
        gradient[2] = dT_de
    else
        q = pericenter(kep)
        da_dq = 1 / (1 - e)
        da_de = q / (1 - e)^2
        gradient[1] = dT_da * da_dq
        gradient[2] = dT_de + dT_da * da_de
    end
    gradient[3] = dT_di
    return propagatedsigma(kep, gradient)
end

function heliocentricc3(kep)
    return -gm(kep) / semimajoraxis(kep)
end

function heliocentricc3sigma(kep)
    μ, a = gm(kep), semimajoraxis(kep)
    gradient = zeros(length(elements(kep)))
    if iselliptic(kep)
        gradient[1] = μ / a^2
    else
        q, e = pericenter(kep), eccentricity(kep)
        gradient[1] = -μ * (e - 1) / q^2
        gradient[2] = μ / q
    end
    return propagatedsigma(kep, gradient)
end

function geocentricc3(attr)
    _, δdeg, vαdeg, vδdeg, ρ, vρ = elements(attr)
    δ = deg2rad(δdeg)
    vα, vδ = deg2rad(vαdeg), deg2rad(vδdeg)
    η2 = vα^2 * cos(δ)^2 + vδ^2
    return vρ^2 + ρ^2 * η2 - 2 * gm(attr) / ρ
end

function geocentricc3sigma(attr)
    _, δdeg, vαdeg, vδdeg, ρ, vρ = elements(attr)
    δ = deg2rad(δdeg)
    vα, vδ = deg2rad(vαdeg), deg2rad(vδdeg)
    η2 = vα^2 * cos(δ)^2 + vδ^2
    gradient = zeros(length(elements(attr)))
    gradient[2] = -2 * ρ^2 * vα^2 * cos(δ) * sin(δ) * (pi / 180)
    gradient[3] = 2 * ρ^2 * vα * cos(δ)^2 * (pi / 180)
    gradient[4] = 2 * ρ^2 * vδ * (pi / 180)
    gradient[5] = 2 * ρ * η2 + 2 * gm(attr) / ρ^2
    gradient[6] = 2 * vρ
    return propagatedsigma(attr, gradient)
end

function c3classification(c3::Real, σ::Real)
    isfinite(c3) && isfinite(σ) && σ >= 0 || return "classification unavailable"
    c3 + σ < 0 && return "two-body bound at >1σ"
    c3 - σ > 0 && return "two-body unbound at >1σ"
    return "bound/unbound indeterminate at 1σ"
end

function printc3(label::AbstractString, c3::Real, σ::Real)
    c3_km, σ_km = C3_AUDAY2_TO_KMSEC2 .* (c3, σ)
    println(rpad(label, 10), "= ", @sprintf("%+.8E", c3_km), " +/- ",
            @sprintf("%.8E", σ_km), " km²/s² (", c3classification(c3, σ), ")")
    return nothing
end

keplerianelementnames(kep) =
    iselliptic(kep) ? ["a", "e", "i", "ω", "Ω", "M"] :
                      ["q", "e", "i", "ω", "Ω", "tp"]

function snrratio(nominal::Real, uncertainty::Real)
    isfinite(nominal) && isfinite(uncertainty) && uncertainty >= 0 || return NaN
    iszero(uncertainty) && return iszero(nominal) ? NaN : Inf
    return abs(nominal) / uncertainty
end

function snrvalue(x::Real)
    isnan(x) && return lpad("NaN", 16)
    isinf(x) && return lpad("Inf", 16)
    return lpad(@sprintf("%.6E", x), 16)
end

function printkepleriansnrs(kep)
    rows = collect(zip(keplerianelementnames(kep), elements(kep), sigmas(kep)))
    names = first.(rows)
    if !("q" in names)
        push!(rows, ("q", pericenter(kep), pericentersigma(kep)))
    end
    if !("tp" in names)
        push!(rows, ("tp", timeperipass(kep), timeperipasssigma(kep)))
    end
    push!(rows, ("Tj", tisserandjupiter(kep), tisserandjupitersigma(kep)))

    printitle("Keplerian SNRs", "*")
    println("SNR = |nominal value| / uncertainty")
    println(rpad("Variable", 12), lpad("SNR", 16))
    for (name, nominal, uncertainty) in rows
        println(rpad(name, 12), snrvalue(snrratio(nominal, uncertainty)))
    end
    println("")
    return nothing
end

function orbitapptype(apps::AbstractApparitionVector, params::Parameters;
                      max_rms_jump_factor::Real = MAX_RMS_REGRESSION_FACTOR,
                      max_rms_jump_arcsec::Real = MAX_RMS_REGRESSION_ARCSEC)
    napps = length(apps)
    napps > 0 || throw(ArgumentError("At least one apparition is required"))
    return orbitapptype(Val(min(napps, 3)), apps, params;
                        max_rms_jump_factor, max_rms_jump_arcsec)
end

function orbitapptype(::Val{1}, apps::AbstractApparitionVector, params::Parameters;
                      max_rms_jump_factor::Real = MAX_RMS_REGRESSION_FACTOR,
                      max_rms_jump_arcsec::Real = MAX_RMS_REGRESSION_ARCSEC)
    return singleapparition(apps, params)
end

function preferorbit(candidate::LeastSquaresOrbit, best::Union{Nothing, LeastSquaresOrbit})
    candidate_nrms = NEOs.nrms(candidate)
    isfinite(candidate_nrms) && candidate_nrms < MAX_ACCEPTABLE_NRMS || return best
    if isnothing(best) || !isfinite(NEOs.nrms(best)) || NEOs.nrms(best) >= MAX_ACCEPTABLE_NRMS
        return candidate
    end
    candidate_nfit, best_nfit = nfitoptical(candidate), nfitoptical(best)
    if candidate_nfit != best_nfit
        return candidate_nfit > best_nfit ? candidate : best
    end
    return NEOs.opticalrms(candidate) < NEOs.opticalrms(best) ? candidate : best
end

function orbitapptype(::Val{2}, apps::AbstractApparitionVector, params::Parameters;
                      max_rms_jump_factor::Real = MAX_RMS_REGRESSION_FACTOR,
                      max_rms_jump_arcsec::Real = MAX_RMS_REGRESSION_ARCSEC)
    optical = NEOs.optical(first(apps))
    seeds = singleapparitionseeds(copy(apps), params)
    isempty(seeds) && return zero(SingleApparitionOrbit{typeof(optical)})
    best = nothing
    for orbitSA in seeds
        candidate = bridge(copy(apps), orbitSA, params;
                           max_rms_jump_factor, max_rms_jump_arcsec)
        iszero(candidate) && continue
        best = preferorbit(candidate, best)
    end
    return isnothing(best) ? first(seeds) : best
end

function orbitapptype(::Val{3}, apps::AbstractApparitionVector, params::Parameters;
                      max_rms_jump_factor::Real = MAX_RMS_REGRESSION_FACTOR,
                      max_rms_jump_arcsec::Real = MAX_RMS_REGRESSION_ARCSEC)
    optical = NEOs.optical(first(apps))
    seeds = singleapparitionseeds(copy(apps), params)
    isempty(seeds) && return zero(SingleApparitionOrbit{typeof(optical)})
    best = nothing
    for orbitSA in seeds
        trialapps = copy(apps)
        orbitMA = bridge(trialapps, orbitSA, params;
                         max_rms_jump_factor, max_rms_jump_arcsec)
        iszero(orbitMA) && continue
        candidate = orbitMA isa MultipleApparitionOrbit ?
                    multipleapparition(trialapps, orbitMA, params;
                                       max_rms_jump_factor, max_rms_jump_arcsec) :
                    orbitMA
        iszero(candidate) && continue
        best = preferorbit(candidate, best)
        if noptical(candidate) == noptical(apps) &&
           critical_value(candidate) < params.significance &&
           all(!isnan, sigmas(candidate))
            return candidate
        end
    end
    return isnothing(best) ? first(seeds) : best
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
    output::Union{Nothing, String} = parsed_args["output"]
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

    # Residual table options
    print_mag_residuals::Bool = !parsed_args["no-mags-res"]
    print_track_residuals::Bool = !parsed_args["no-track-res"]
    force_all_fit::Bool = parsed_args["force-all-fit"]
    force_all_fit && println("• Final refinement will keep all fitted observations if it converges")
    max_rms_jump_factor::Float64 = parsed_args["max-rms-jump-factor"]
    max_rms_jump_arcsec::Float64 = parsed_args["max-rms-jump-arcsec"]
    max_rms_jump_factor >= 1 ||
        throw(ArgumentError("--max-rms-jump-factor must be at least 1"))
    max_rms_jump_arcsec >= 0 ||
        throw(ArgumentError("--max-rms-jump-arcsec must be non-negative"))
    println("• RMS jump guard: factor ", max_rms_jump_factor,
            ", active above ", max_rms_jump_arcsec, " arcsec")

    # Post-fit ephemeris options
    ephem_days::Float64 = parsed_args["ephem-days"]
    ephem_step::Float64 = parsed_args["ephem-step"]
    ephem_obscode::String = parsed_args["ephem-obscode"]
    ephem_days >= 0 || throw(ArgumentError("--ephem-days must be non-negative"))
    ephem_step > 0 || throw(ArgumentError("--ephem-step must be positive"))
    ephem_observatory = ephem_days > 0 ? search_observatory_code(ephem_obscode) : nothing
    if ephem_days > 0
        NEOs.isunknown(ephem_observatory) && throw(ArgumentError(
            "Unknown MPC observatory code for --ephem-obscode: $ephem_obscode"
        ))
        println("• Post-fit ephemeris: ", ephem_days, " days at ", ephem_step,
                "-day steps from observatory ", ephem_obscode)
    else
        println("• Post-fit ephemeris: disabled")
    end

    # Global initial time
    global_initial_time = now()
    println("• Run started at ", global_initial_time)

    # Load optical astrometry
    optical, format = load_optical_astrometry(input; format)
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
        coeffstol = Inf, bwdoffset = 0.05, fwdoffset = 0.05,
        marsden_radial = use_cometary_nongravs ? COMETARY_MARSDEN_RADIAL :
                         (1.0, 1.0, 2.0, 0.0, 0.0),
        gaussorder = 2, safegauss = true, refscale = :log,
        tsaorder = 2, adamiter = 500, adamQtol = 1E-5,
        jtlsorder = 2, jtlsmask = false, jtlsiter = 20, lsiter = 30,
        jtlsproject = true, significance = 0.99, verbose = true,
        outrej = true, χ2_rec = 4.0, χ2_rej = 5.0, fudge = 100.0,
        max_per = 33.3
    )

    # Split observational arc into apparitions
    apps = apparitions(fit_optical, Day(split_gap_days))
    # Compute orbit by apparition type (single/multiple apparition)
    orbit = orbitapptype(apps, params; max_rms_jump_factor, max_rms_jump_arcsec)
    iszero(orbit) && error(
        "Unable to determine a preliminary orbit from the fitted observations"
    )
    orbit = finalrefinement(final_dynamics, orbit, fit_optical, params;
                            marsden_scalings = nongrav_scalings,
                            force_all_fit, max_rms_jump_factor, max_rms_jump_arcsec)

    # JTLS fits at the weighted mean observational epoch; only shift when requested.
    if !isnothing(solution_epoch)
        try
            orbit = shiftepoch(orbit, solution_epoch, params; beyondarc = true)
        catch err
            println("• Could not shift final orbit epoch; keeping orbit at fitted epoch (",
                    typeof(err), ")")
        end
    end
    # Compute elements and H before attaching held-out observations because H currently
    # uses every optical observation and does not inspect residual outlier flags.
    kep = keplerian(orbit, params)
    attr = attributable(orbit, params)
    H, dH = absolutemagnitude(orbit, params)
    full_residuals_attached = false
    auto_excluded = missingopticalindices(orbit, optical)
    residual_excluded = sort!(unique!(vcat(fit_excluded, auto_excluded)))
    if !isempty(auto_excluded)
        println("• Keeping ", length(auto_excluded),
                " automatically skipped observation(s) as outliers in residual output")
    end
    if !isempty(residual_excluded)
        try
            orbit = attachfullresiduals(orbit, optical, residual_excluded, params)
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
    println("q  = ", @sprintf("%+.12E", pericenter(kep)), " +/- ",
            @sprintf("%.12E", pericentersigma(kep)), " au")
    println("tp = ", @sprintf("%+.12E", timeperipass(kep)), " +/- ",
            @sprintf("%.12E", timeperipasssigma(kep)), " MJD TDB")
    println("Tj = ", @sprintf("%.8f", tisserandjupiter(kep)), " +/- ",
            @sprintf("%.8f", tisserandjupitersigma(kep)))
    printc3("C3(Sun)", heliocentricc3(kep), heliocentricc3sigma(kep))
    printc3("C3(Earth)", geocentricc3(attr), geocentricc3sigma(attr))
    println("H  = ", @sprintf("%.3f", H), " +/- ", @sprintf("%.3f", dH), " mag")
    println("")
    printkepleriansnrs(kep)
    ele220_name = ele220inputname(input, format, designation(orbit))
    printele220arguments(ele220_name, orbit, kep, H, params, gap_days)
    try
        printpostfitephemeris(orbit, params, H, ephem_days, ephem_step,
                              ephem_observatory)
    catch err
        print("• Could not compute post-fit ephemeris: ")
        showerror(stdout, err)
        println("")
    end
    printmotionbytracklet(orbit, params)
    printmpeclikeresiduals(orbit, H, params; include_mags = print_mag_residuals,
                           include_track = print_track_residuals)

    # Save orbit
    if !isnothing(output)
        if isempty(residual_excluded)
            jldsave(output; orbit)
        elseif full_residuals_attached
            heldout_optical = orbit.optical[residual_excluded]
            heldout_ores = orbit.ores[residual_excluded]
            jldsave(output; orbit, fit_excluded = residual_excluded,
                    heldout_optical, heldout_ores)
        else
            heldout_optical = optical[residual_excluded]
            heldout_ores = typeof(orbit.ores)()
            jldsave(output; orbit, fit_excluded = residual_excluded,
                    heldout_optical, heldout_ores)
        end
        println("Final orbit saved to: ", output)
    end

    # Final time
    global_final_time = now()
    println("• Run started ", global_initial_time, " and finished ", global_final_time)
    global_computation_time = computationtime(global_initial_time, global_final_time)
    println("• Total computation time was: ", global_computation_time, " min")

    return nothing
end

main()
