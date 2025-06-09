using ArgParse, NEOs, XML, Dates, LinearAlgebra, PlanetaryEphemeris
using TaylorSeries, StatsBase, JLD2
using Plots, LaTeXStrings, Printf, Measures, PrettyTables
using NEOs: RadecMPC, AdmissibleRegion, cte, isclosed, euclid3D, jtls,
      wra, wdec, iodinitcond, propres, set_od_order, covariance

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
        adamiter = 500, adamQtol = 1e-5,                       # MMOV
        jtlsiter = 20, lsiter = 10, significance = 0.99,       # Least squares
        outrej = true, χ2_rec = 7.0, χ2_rej = 8.0,             # Outlier rejection
        fudge = 100.0, max_per = 20.0
    )
    od = ODProblem(newtonian!, radec[1:6], UniformWeights, ZeroDebiasing)
    od.w8s.w8s .= Ref((1 / 0.2^2, 1 / 0.2^2))
    # Gauss method
    porbits = gaussmethod(od, params)
    println("Roots of Gauss polynomial")
    for porbit in porbits
        println("r_2 = ", @sprintf("%.3f", porbit.r_2), " au")
    end
    # Filter Gauss solutions by heliocentric energy
    filter!(isclosed, porbits)
    println("Roots with negative heliocentric energy")
    for porbit in porbits
        println("r_2 = ", @sprintf("%.3f", porbit.r_2), " au")
    end
    # First root
    orbit1 = jtls(od, porbits[1], params, true)
    println("First root converges via JT least squares to an orbit with ",
        "nrms = ", @sprintf("%.3f", nrms(orbit1)))
    # Second root
    orbit2 = jtls(od, porbits[2], params, true)
    println("Second root converges via JT least squares to an orbit with ",
        "nrms = ", @sprintf("%.3f", nrms(orbit2)))
    # Include remaining observations
    NEOs.update!(od, radec)
    od.w8s.w8s .= Ref((1 / 0.2^2, 1 / 0.2^2))
    orbit12 = orbitdetermination(od, orbit1, params)
    println("First root converges via JT least squares to an orbit with ",
        "nrms = ", @sprintf("%.3f", nrms(orbit12)), " when considering the full ",
        "set of optical astrometry")
    orbit22 = orbitdetermination(od, orbit2, params)
    println("Second root converges via JT least squares to an orbit with ",
        "nrms = ", @sprintf("%.3f", nrms(orbit22)), " when considering the full ",
        "set of optical astrometry")
    # Residuals statistics
    ma = @sprintf("%.1e", mean(ra.(orbit22.res), weights(wra.(orbit22.res))))
    md = @sprintf("%.1e", mean(dec.(orbit22.res), weights(wdec.(orbit22.res))))
    println("Residuals statistics of latter orbital fit:")
    println("Mean weighted right ascension residual: ", ma)
    println("Mean weighted declination residual: ", md)

    return nothing
end

function mmovplot(A::AdmissibleRegion{T}, aes::Matrix{T}, mov::Matrix{T},
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
        adamiter = 500, adamQtol = 1e-5,                       # MMOV
        jtlsiter = 20, lsiter = 10, significance = 0.99,       # Least squares
        outrej = true, χ2_rec = 7.0, χ2_rej = 8.0,             # Outlier rejection
        fudge = 100.0, max_per = 20.0
    )
    od = ODProblem(newtonian!, radec)
    # List of naive initial conditions
    A = AdmissibleRegion(od.tracklets[1], params)
    I0 = iodinitcond(A)
    # Minimization over manifold of variations
    ρ, v_ρ, scale = I0[1]
    porbit = mmov(od, 1, A, ρ, v_ρ, params; scale)
    aes, Qs = porbit.aes, porbit.Qs
    println("MMOV converges in ", length(Qs), " iterations to an orbit with ",
        "nrms = ", @sprintf("%.3f", nrms(porbit)))
    # Jet transport least squares
    orbit = jtls(od, porbit, params, true)
    println("Two further JT least squares iterations suffice to produce an orbit ",
        "with nrms = ", @sprintf("%.3f", nrms(orbit)))
    dist2jpl = @sprintf("%.1f", maximum(jplcompare("2014 AA", orbit)))
    println("Distance to JPL's orbit in units of the maximum uncertainty: ",
        dist2jpl, "σ")
    # Table of osculating elements
    set_od_order(Float64, 2)
    t = epoch(orbit)
    jd0 = t + PE.J2000
    q0 = orbit(t) + diag(orbit.J) .* get_variables(Float64, 2)
    params = Parameters(params; bwdoffset = 0.5)
    bwd, _, _ = propres(od, jd0, q0, params)
    jd0, t = 2456658.5, 2456658.5 - J2000
    x0 = zeros(Float64, 6)
    osc = pv2kep(bwd(t) - params.eph_su(t); jd = jd0, frame = :ecliptic)
    ele = [osc.a, osc.e, osc.i, osc.Ω, osc.ω, mod(osc.M, 360)]
    t_car2kep = TS.jacobian(ele, x0)
    Γ_E = t_car2kep * covariance(orbit) * t_car2kep'
    σ_E = sqrt.(diag(Γ_E))
    ss = [@sprintf("%.3f", cte(e[1])) * " ± " * @sprintf("%.3f", e[2]) for e in zip(ele, σ_E)]
    jplele = [1.161570914329746, 0.210903385513902, 1.4109467942359,
        101.613014674528, 52.3157504212017, 324.0051879996831]
    jplsig = [0.0042653, 0.0042378, 0.029249, 0.0096696, 0.19112, 0.18696]
    jplss = [@sprintf("%.3f", e[1]) * " ± " * @sprintf("%.3f", e[2]) for e in zip(jplele, jplsig)]
    println("Comparison of the orbital elements obtained for 2014 AA at epoch ",
        @sprintf("%.2f", jd0), " JD using MMOV with the JPL ",
        " result (Farnocchia et al., 2016):")
    pretty_table(
    [
        "JPL"       jplss...
        "MMOV"   ss...
    ],
    header = ["", "a [au]", "e", "i [deg]", "Ω [deg]", "ω [deg]", "M [deg]"],
    alignment = :c
    )
    # MMOV Plot
    mov = JLD2.load(filename, "mov")
    filename = joinpath(dirname(filename), "2014AA.pdf")
    mmovplot(A, aes, mov, filename; scale = :log, Qmax = 1e5)
    println("MMOV plot saved to: ", filename)
end

function neo2020CV1(filename::String)
    println("• 2020 CV1")
    # Load astrometry
    radec = fetch_radec_mpc("2020 CV1")
    # Orbit determination parameters and problem
    params = Parameters(
        coeffstol = Inf, bwdoffset = 0.042, fwdoffset = 0.042, # Propagation
        adamiter = 500, adamQtol = 1e-5,                       # MMOV
        jtlsiter = 20, lsiter = 10, significance = 0.99,       # Least squares
        outrej = true, χ2_rec = 7.0, χ2_rej = 8.0,             # Outlier rejection
        fudge = 100.0, max_per = 20.0
    )
    od = ODProblem(newtonian!, radec[7:end])
    # List of naive initial conditions
    A = AdmissibleRegion(od.tracklets[1], params)
    I0 = [(mean(A.ρ_domain), mean(A.v_ρ_domain), :linear)]
    # Minimization over manifold of variations
    ρ, v_ρ, scale = I0[1]
    porbit = mmov(od, 1, A, ρ, v_ρ, params; scale)
    aes, Qs = porbit.aes, porbit.Qs
    println("MMOV converges in ", length(Qs), " iterations to an orbit with ",
        "nrms = ", @sprintf("%.3f", nrms(porbit)))
    # Jet transport least squares
    orbit = jtls(od, porbit, params, true)
    println("Four further JT least squares iterations suffice to produce an orbit ",
        "with nrms = ", @sprintf("%.3f", nrms(orbit)))
    dist2jpl = @sprintf("%.1f", maximum(jplcompare("2020 CV1", orbit)))
    println("Distance to JPL's orbit in units of the maximum uncertainty: ",
        dist2jpl, "σ")
    # MMOV Plot
    mov = JLD2.load(filename, "mov")
    filename = joinpath(dirname(filename), "2020CV1.pdf")
    mmovplot(A, aes, mov, filename; scale = :linear, Qmax = 1e10)
    println("MMOV plot saved to: ", filename)
end

function neo2024MK(filename::String)
    println("• 2024 MK")
    # Load astrometry
    radec = fetch_radec_mpc("2024 MK")
    # Orbit determination parameters and problem
    params = Parameters(
        coeffstol = Inf, bwdoffset = 0.042, fwdoffset = 0.042, # Propagation
        adamiter = 500, adamQtol = 1e-5,                       # MMOV
        jtlsiter = 20, lsiter = 10, significance = 0.99,       # Least squares
        outrej = true, χ2_rec = 7.0, χ2_rej = 8.0,             # Outlier rejection
        fudge = 100.0, max_per = 20.0
    )
    od = ODProblem(newtonian!, radec[10:21])
    # Gauss method
    porbits = gaussmethod(od, params)
    println("The Gauss polynomial has a single root:")
    for porbit in porbits
        println("r_2 = ", @sprintf("%.3f", porbit.r_2), " au")
    end
    print("which fails to converge in the LS step, because nrms = ",
        @sprintf("%.1e", nrms(porbits[1])))
    # MMOV refinement
    porbit = mmov(od, porbits[1], 2, :log, params)
    aes, Qs = porbit.aes, porbit.Qs
    r_2 = euclid3D(porbit() - params.eph_su(epoch(porbit)))
    println("Gauss preliminary orbit is refined with ", length(Qs),
        " iterations of minimization over the MOV to an orbit with ",
        "r_2 = ", @sprintf("%.3f", r_2), " au and nrms = ", @sprintf("%.3f", nrms(porbit)))
    # Jet transport least squares
    orbit = jtls(od, porbit, params, true)
    println("Two further JT least squares iterations suffice to produce an orbit ",
        "with nrms = ", @sprintf("%.3f", nrms(orbit)))
    dist2jpl = @sprintf("%.1f", maximum(jplcompare("2024 MK", orbit)))
    println("Distance to JPL's orbit in units of the maximum uncertainty: ",
        dist2jpl, "σ")
    # MMOV Plot
    mov = JLD2.load(filename, "mov")
    filename = joinpath(dirname(filename), "2024MK.pdf")
    A = AdmissibleRegion(od.tracklets[2], params)
    mmovplot(A, aes, mov, filename; scale = :log, Qmax = 1e10,
        bounds = [-2.2, 0.0, -0.01, 0.002])
    println("MMOV plot saved to: ", filename)
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

    # Set orbit determination order
    set_od_order(Float64, 6)
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