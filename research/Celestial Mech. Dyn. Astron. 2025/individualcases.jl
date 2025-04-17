using ArgParse, NEOs, XML, Dates, LinearAlgebra, PlanetaryEphemeris
using TaylorSeries, StatsBase, JLD2
using Plots, LaTeXStrings, Printf, Measures, PrettyTables
using NEOs: RadecMPC, ObservatoryMPC, AdmissibleRegion, scaled_variables,
    _gausstriplet!, cte, heliocentric_energy, euclid3D, jtls, wra, wdec,
    iodinitcond, adam, attr2bary, propres, bary2topo, boundary_projection

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "individualcases.jl"
    # Desciption (for help screen)
    s.description = """
    Reproduce the plots and tables in Ramírez-Montoya et al, 2024, section 5.1."""

    @add_arg_table! s begin
        "--input1"
            help = "2024 PDC25 .xml astrometry file"
            arg_type = String
            default = ""
        "--input2"
            help = "2014 AA MOV .jld2 file"
            arg_type = String
            default = ""
        "--input3"
            help = "2020 CV1 MOV .jld2 file"
            arg_type = String
            default = ""
        "--input4"
            help = "2024 MK MOV .jld2 file"
            arg_type = String
            default = ""
    end

    s.epilog = """
        notes:\n
        • If a --input is empty, the corresponding NEO will be omitted\n
        • The MOV .jld2 files can be computed with the mcmov.jl script\n
        \n
        example:\n
        \n
        julia -t 5 --project scripts/individualcases.jl --input1 2024PDC25.xml
        --input2 2014AAMOV.jld2 --input3 2020CV1MOV.jld2 --input4 2024MKMOV.jld2\n
        \n
    """

    return parse_args(s)
end

function parsexml(filename::String)
    # Read XML file as Node
    n = read(filename, Node)
    # Get optical astrometry in ADES format
    obs = n.children[2].children
    # Converte observations to NEOs format
    radec = Vector{RadecMPC{Float64}}(undef, length(obs))
    for i in eachindex(obs)
        o = obs[i].children
        observatory = search_obs_code(o[3].children[1].value)
        date = DateTime(o[4].children[1].value[1:end-1])
        α = deg2rad(parse(Float64, o[5].children[1].value))
        δ = deg2rad(parse(Float64, o[6].children[1].value))
        catalogue = unknowncat()
        mag = parse(Float64, o[10].children[1].value)
        band = o[12].children[1].value
        radec[i] = RadecMPC{Float64}("", "K24P25DC", "", "", "C", date,
            α, δ, "", mag, band, catalogue, "", observatory)
    end
    sort!(radec)

    return radec
end

function neo2024PDC25(filename::String)
    println("• 2024 PDC25")
    # Load astrometry
    radec = parsexml(filename)
    # Orbit determination parameters and problem
    params = Parameters(
        coeffstol = Inf, bwdoffset = 0.042, fwdoffset = 0.042, # Propagation
        adamiter = 500, adamQtol = 1e-5,                       # ADAM
        jtlsiter = 20, lsiter = 10, significance = 0.99,       # Least squares
        outrej = true, χ2_rec = 7.0, χ2_rej = 8.0,             # Outlier rejection
        fudge = 100.0, max_per = 20.0
    )
    od = ODProblem(newtonian!, radec[1:6], UniformWeights, ZeroDebiasing)
    od.w8s.w8s .= [(1 / 0.2^2, 1 / 0.2^2) for _ in 1:6]
    # Declare jet transport variables
    varorder = max(params.tsaorder, params.gaussorder, params.jtlsorder)
    dq = scaled_variables("dx", fill(1e-6, 6); order = varorder)
    # Choose a three observations for Gauss Method
    observatories = Vector{ObservatoryMPC{Float64}}(undef, 3)
    dates = Vector{DateTime}(undef, 3)
    α, δ = Vector{Float64}(undef, 3), Vector{Float64}(undef, 3)
    _gausstriplet!(observatories, dates, α, δ, od.tracklets)
    _jd0_ = dtutc2jdtdb(dates[2])
    # Gauss method
    solG = gauss_method(observatories, dates, α .+ dq[1:3], δ .+ dq[4:6], params)
    println("Roots of Gauss polynomial")
    for i in eachindex(solG)
        r_2 = euclid3D(cte.(solG[i].statevect))
        println("r_2 = ", @sprintf("%.3f", r_2), " au")
    end
    # Filter Gauss solutions by heliocentric energy
    filter!(x -> cte(cte(heliocentric_energy(x.statevect))) <= 0, solG)
    println("Roots with negative heliocentric energy")
    for i in eachindex(solG)
        r_2 = euclid3D(cte.(solG[i].statevect))
        println("r_2 = ", @sprintf("%.3f", r_2), " au")
    end
    # First root
    jd0 = _jd0_ - cte(solG[1].ρ[2]) / c_au_per_day
    q = solG[1].statevect .+ params.eph_su(jd0 - PE.J2000);
    sol1 = jtls(od, jd0, q, od.tracklets, params, true)
    println("First root converges via JT least squares to an orbit with ",
        "nrms = ", @sprintf("%.3f", nrms(sol1)))
    # Second root
    jd0 = _jd0_ - cte(solG[2].ρ[2]) / c_au_per_day
    q = solG[2].statevect .+ params.eph_su(jd0 - PE.J2000);
    sol2 = jtls(od, jd0, q, od.tracklets, params, true)
    println("Second root converges via JT least squares to an orbit with ",
        "nrms = ", @sprintf("%.3f", nrms(sol2)))
    # Include remaining observations
    NEOs.update!(od, radec)
    od.w8s.w8s .= [(1 / 0.2^2, 1 / 0.2^2) for _ in eachindex(radec)]
    sol12 = orbitdetermination(od, sol1, params)
    println("First root converges via JT least squares to an orbit with ",
        "nrms = ", @sprintf("%.3f", nrms(sol12)), " when considering the full ",
        "set of optical astrometry")
    sol22 = orbitdetermination(od, sol2, params)
    println("Second root converges via JT least squares to an orbit with ",
        "nrms = ", @sprintf("%.3f", nrms(sol22)), " when considering the full ",
        "set of optical astrometry")
    # Residuals statistics
    ma = @sprintf("%.1e", mean(ra.(sol22.res), weights(wra.(sol22.res))))
    md = @sprintf("%.1e", mean(dec.(sol22.res), weights(wdec.(sol22.res))))
    println("Residuals statistics of latter orbital fit:")
    println("Mean weighted right ascension residual: ", ma)
    println("Mean weighted declination residual: ", md)

    return nothing
end

function adamplot(A::AdmissibleRegion{T}, aes::Matrix{T}, mov::Matrix{T},
    filename::String; scale::Symbol = :log, Qmax::T = 1e10,
    bounds = nothing) where {T <: Real}
    mask = mov[13, :] .< Qmax
    xlabel = scale == :log ? L"\log_{10}(\rho \ \textrm{[au]})" :
        L"\rho \ \textrm{[au]}"
    scatter(mov[11, mask], mov[12, mask], zcolor = log10.(mov[13, mask]),
        label = "", color = :deep, markershape = :square, markersize = 2.3,
        markerstrokewidth = 0, xlabel = xlabel,
        ylabel = L"\dot{\rho} \ \textrm{[au/day]}", colorbar_title = L"\log_{10}(Q)",
        dpi = 300, fontfamily = "Computer Modern", titlefont = (14),
        guidefont = (12, :black), tickfont = (10, :black), framestyle = :box,
        yminorgrid = true, xminorgrid = true, margin = 0.5Measures.mm,)
    plot!(A, N = 10_000, ρscale = scale, label = "", c = :black, linewidth = 3)
    plot!(A, boundary = :inner, N = 10_000, ρscale = scale, color = :black,
        label = "", linewidth = 2)
    xs = scale == :log ? log10.(aes[5, :]) : aes[5, :]
    shapes, sizes = fill(:circle, length(xs)), fill(4.0, length(xs))
    shapes[1], shapes[end] = :star5, :square
    sizes[1], sizes[end] = 5.5, 4.5
    scatter!(xs, aes[6, :], markersize = sizes, color = :white,
        markershape = shapes, markerstrokecolor = :black, label = "")
    if !isnothing(bounds)
        xlims!(bounds[1], bounds[2])
        ylims!(bounds[3], bounds[4])
    end
    savefig(filename)
end

function neo2014AA(filename::String)
    println("• 2014 AA")
    # Load astrometry
    radec = fetch_radec_mpc("2014 AA")
    # Orbit determination parameters and problem
    params = Parameters(
        coeffstol = Inf, bwdoffset = 0.042, fwdoffset = 0.042, # Propagation
        adamiter = 500, adamQtol = 1e-5,                       # ADAM
        jtlsiter = 20, lsiter = 10, significance = 0.99,       # Least squares
        outrej = true, χ2_rec = 7.0, χ2_rej = 8.0,             # Outlier rejection
        fudge = 100.0, max_per = 20.0
    )
    od = ODProblem(newtonian!, radec)
    # Declare jet transport variables
    varorder = max(params.tsaorder, params.gaussorder, params.jtlsorder)
    scaled_variables("dx", ones(Float64, 6); order = varorder)
    # List of naive initial conditions
    A = AdmissibleRegion(od.tracklets[1], params)
    I0 = iodinitcond(A)
    # ADAM minimization over manifold of variations
    ρ, v_ρ, scale = I0[1]
    aes, Qs = adam(od, 1, A, ρ, v_ρ, params; scale)
    ae, Q = aes[:, end], Qs[end]
    println("ADAM converges in ", length(Qs), " iterations to an orbit with ",
        "nrms = ", @sprintf("%.3f", sqrt(Q)))
    # Jet transport least squares
    jd0 = dtutc2jdtdb(A.date) - ae[5] / c_au_per_day
    q0 = attr2bary(A, ae, params)
    q = [q0[k] + (abs(q0[k]) / 10^5) * TaylorN(k, order = varorder) for k in 1:6]
    sol = jtls(od, jd0, q, od.tracklets, params, true)
    println("Two further JT least squares iterations suffice to produce an orbit ",
        "with nrms = ", @sprintf("%.3f", nrms(sol)))
    dist2jpl = @sprintf("%.1f", maximum(jplcompare("2014 AA", sol)))
    println("Distance to JPL's orbit in units of the maximum uncertainty: ",
        dist2jpl, "σ")
    # Table of osculating elements
    jd0 = epoch(sol) + J2000
    q0 .= sol()
    q .= [q0[k] + (abs(q0[k]) / 10^5) * TaylorN(k, order = varorder) for k in 1:6]
    params = Parameters(params; bwdoffset = 0.5)
    bwd, _, _ = propres(od, jd0, q, params)
    jd0, t = 2456658.5, 2456658.5 - J2000
    osc = pv2kep(bwd(t) - params.eph_su(t); jd = jd0, frame = :ecliptic)
    ele = [osc.a, osc.e, osc.i, osc.Ω, osc.ω, osc.M + 360]
    Γ_E = project(ele, sol.fit)
    σ_E = sqrt.(diag(Γ_E))
    ss = [@sprintf("%.3f", cte(e[1])) * " ± " * @sprintf("%.3f", e[2]) for e in zip(ele, σ_E)]
    jplele = [1.161570914329746, 0.210903385513902, 1.4109467942359,
        101.613014674528, 52.3157504212017, 324.0051879996831]
    jplsig = [0.0042653, 0.0042378, 0.029249, 0.0096696, 0.19112, 0.18696]
    jplss = [@sprintf("%.3f", e[1]) * " ± " * @sprintf("%.3f", e[2]) for e in zip(jplele, jplsig)]
    println("Comparison of the orbital elements obtained for 2014 AA at epoch ",
        @sprintf("%.2f", jd0), " JD using ADAM-JT with the JPL ",
        " result (Farnocchia et al., 2016):")
    pretty_table(
    [
        "JPL"       jplss...
        "ADAM-JT"   ss...
    ],
    header = ["", "a [au]", "e", "i [deg]", "Ω [deg]", "ω [deg]", "M [deg]"],
    alignment = :c
    )
    # ADAM Plot
    mov = JLD2.load(filename, "mov")
    filename = joinpath(dirname(filename), "2014AA.pdf")
    adamplot(A, aes, mov, filename; scale = :log, Qmax = 1e5)
    println("ADAM plot saved to: ", filename)
end

function neo2020CV1(filename::String)
    println("• 2020 CV1")
    # Load astrometry
    radec = fetch_radec_mpc("2020 CV1")
    # Orbit determination parameters and problem
    params = Parameters(
        coeffstol = Inf, bwdoffset = 0.042, fwdoffset = 0.042, # Propagation
        adamiter = 500, adamQtol = 1e-5,                       # ADAM
        jtlsiter = 20, lsiter = 10, significance = 0.99,       # Least squares
        outrej = true, χ2_rec = 7.0, χ2_rej = 8.0,             # Outlier rejection
        fudge = 100.0, max_per = 20.0
    )
    od = ODProblem(newtonian!, radec[7:end])
    # Declare jet transport variables
    varorder = max(params.tsaorder, params.gaussorder, params.jtlsorder)
    scaled_variables("dx", ones(Float64, 6); order = varorder)
    # List of naive initial conditions
    A = AdmissibleRegion(od.tracklets[1], params)
    I0 = [(mean(A.ρ_domain), mean(A.v_ρ_domain), :linear)]
    # ADAM minimization over manifold of variations
    ρ, v_ρ, scale = I0[1]
    aes, Qs = adam(od, 1, A, ρ, v_ρ, params; scale)
    ae, Q = aes[:, end], Qs[end]
    println("ADAM converges in ", length(Qs), " iterations to an orbit with ",
        "nrms = ", @sprintf("%.3f", sqrt(Q)))
    # Jet transport least squares
    jd0 = dtutc2jdtdb(A.date) - ae[5] / c_au_per_day
    q0 = attr2bary(A, ae, params)
    q = [q0[k] + (abs(q0[k]) / 10^5) * TaylorN(k, order = varorder) for k in 1:6]
    sol = jtls(od, jd0, q, od.tracklets, params, true)
    println("Four further JT least squares iterations suffice to produce an orbit ",
        "with nrms = ", @sprintf("%.3f", nrms(sol)))
    dist2jpl = @sprintf("%.1f", maximum(jplcompare("2020 CV1", sol)))
    println("Distance to JPL's orbit in units of the maximum uncertainty: ",
        dist2jpl, "σ")
    # ADAM Plot
    mov = JLD2.load(filename, "mov")
    filename = joinpath(dirname(filename), "2020CV1.pdf")
    adamplot(A, aes, mov, filename; scale = :linear, Qmax = 1e10)
    println("ADAM plot saved to: ", filename)
end

function neo2024MK(filename::String)
    println("• 2024 MK")
    # Load astrometry
    radec = fetch_radec_mpc("2024 MK")
    # Orbit determination parameters and problem
    params = Parameters(
        coeffstol = Inf, bwdoffset = 0.042, fwdoffset = 0.042, # Propagation
        adamiter = 500, adamQtol = 1e-5,                       # ADAM
        jtlsiter = 20, lsiter = 10, significance = 0.99,       # Least squares
        outrej = true, χ2_rec = 7.0, χ2_rej = 8.0,             # Outlier rejection
        fudge = 100.0, max_per = 20.0
    )
    od = ODProblem(newtonian!, radec[10:21])
    # Declare jet transport variables
    varorder = max(params.tsaorder, params.gaussorder, params.jtlsorder)
    dq = scaled_variables("dx", ones(Float64, 6); order = varorder)
    # Choose a three observations for Gauss Method
    observatories = Vector{ObservatoryMPC{Float64}}(undef, 3)
    dates = Vector{DateTime}(undef, 3)
    α, δ = Vector{Float64}(undef, 3), Vector{Float64}(undef, 3)
    _gausstriplet!(observatories, dates, α, δ, od.tracklets)
    # Gauss method
    solG = gauss_method(observatories, dates, α .+ dq[1:3], δ .+ dq[4:6], params)
    println("The Gauss polynomial has a single root:")
    for i in eachindex(solG)
        r_2 = euclid3D(cte.(solG[i].statevect))
        println("r_2 = ", @sprintf("%.3f", r_2), " au")
    end
    print("which fails to converge in the LS step, because ")
    # ADAM refinement
    jd0 = dtutc2jdtdb(dates[2]) - cte(solG[1].ρ[2]) / c_au_per_day
    q = solG[1].statevect .+ params.eph_su(jd0 - J2000)
    bwd, fwd, res = propres(od, jd0, q(), params)
    println(" nrms = ", @sprintf("%.1e", nrms(res)))
    A = AdmissibleRegion(od.tracklets[2], params)
    At0 = dtutc2days(A.date)
    q0 = bwdfwdeph(At0, bwd, fwd, false, false)
    ρ, v_ρ = bary2topo(A, q0)
    ρ, v_ρ = boundary_projection(A, ρ, v_ρ)
    aes, Qs = adam(od, 2, A, ρ, v_ρ, params; scale = :log)
    ae, Q = aes[:, end], Qs[end]
    r_2 = euclid3D(solG[1].R_vec[2, :] + ae[5] * cte.(solG[1].ρ_vec[2, :]))
    println("Gauss preliminary orbit is refined with ", length(Qs),
        " iterations of ADAM minimization over the MOV to an orbit with ",
        "r_2 = ", @sprintf("%.3f", r_2), " au and nrms = ", @sprintf("%.3f", sqrt(Q)))
    # Jet transport least squares
    jd0 = dtutc2jdtdb(A.date) - ae[5] / c_au_per_day
    q0 = attr2bary(A, ae, params)
    q .= [q0[k] + (abs(q0[k]) / 10^5) * TaylorN(k, order = varorder) for k in 1:6]
    sol = jtls(od, jd0, q, od.tracklets, params, true)
    println("Two further JT least squares iterations suffice to produce an orbit ",
        "with nrms = ", @sprintf("%.3f", nrms(sol)))
    dist2jpl = @sprintf("%.1f", maximum(jplcompare("2024 MK", sol)))
    println("Distance to JPL's orbit in units of the maximum uncertainty: ",
        dist2jpl, "σ")
    # ADAM Plot
    mov = JLD2.load(filename, "mov")
    filename = joinpath(dirname(filename), "2024MK.pdf")
    adamplot(A, aes, mov, filename; scale = :log, Qmax = 1e10,
        bounds = [-2.2, 0.0, -0.01, 0.002])
    println("ADAM plot saved to: ", filename)
end

function main()
    # Initial time
    init_time = now()
    # Parse arguments from commandline
    parsed_args = parse_commandline()
    # Input files
    input1::String = parsed_args["input1"]
    input2::String = parsed_args["input2"]
    input3::String = parsed_args["input3"]
    input4::String = parsed_args["input4"]
    # Print header
    println("Reproduce the plots and tables in Ramírez-Montoya et al, 2024, section 5.1.")
    println("• Input file for 2024 PDC25: ", isempty(input1) ? "omitted" : input1)
    println("• Input file for 2014 AA: ", isempty(input2) ? "omitted" : input2)
    println("• Input file for 2020 CV1: ", isempty(input3) ? "omitted" : input3)
    println("• Input file for 2024 MK: ", isempty(input4) ? "omitted" : input4)

    # 2024 PDC25
    !isempty(input1) && neo2024PDC25(input1)
    # 2014 AA
    !isempty(input2) && neo2014AA(input2)
    # 2020 CV1
    !isempty(input3) && neo2020CV1(input3)
    # 2024 MK
    !isempty(input4) && neo2024MK(input4)

    # Final time
    final_time = now()
    println("• Run started ", init_time, " and finished ", final_time)

    return nothing
end

main()