using ArgParse, NEOs, JLD2, StatsBase, Dates
using Plots, LaTeXStrings, Measures, Printf, PrettyTables
using NEOs: NEOSolution, numberofdays

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "largescale.jl"
    # Desciption (for help screen)
    s.description = """
    Reproduce the plots and tables in Ramírez-Montoya et al, 2024, section 5.2."""

    @add_arg_table! s begin
        "--input1"
            help = "input directory for sample 1"
            arg_type = String
            default = "orbits5000"
        "--input2"
            help = "input directory for sample 2"
            arg_type = String
            default = "orbits2024"
    end

    s.epilog = """
        example:\n
        \n
        julia --project scripts/largescale.jl --input1 orbits5000 --input2 orbits2024\n
        \n
    """

    return parse_args(s)
end

function qhist(Qs::Vector{T}, mask::BitVector, filename::String;
    ymax::T = 600.0, ysub::T = 550.0, sub::String = "(a)") where {T <: Real}
    p = @sprintf("%.2f", count(mask) * 100 / length(Qs))
    l = "$(count(mask))/$(length(Qs)) NEOs ($p%)"
    histogram(Qs[mask], label = l, xlabel = "NRMS", ylabel = "Number of NEOs",
        xlim = (0, 1.5), xticks = 0:0.1:1.5, ylim = (0, ymax))
    annotate!(0.1, ysub, sub)
    savefig(filename)
end

function shist(snrs::Vector{T}, mask::BitVector, filename::String;
    ymax::T = 600.0, ysub::T = 550.0, sub::String = "(a)") where {T <: Real}
    p = @sprintf("%.2f", count(mask) * 100 / length(snrs))
    l = "$(count(mask))/$(length(snrs)) NEOs ($p%)"
    histogram(log10.(snrs[mask]), label = l,
        xlabel = L"\log_{10}(\textrm{minimum \ SNR})",
        ylabel = "Number of NEOs", xlim = (-4, 5), ylim = (0, ymax),
        xticks = -4:5, yticks = 0:50:ymax, legend = :topleft, bins = -4:0.2:5)
    annotate!(4.5, ysub, sub)
    savefig(filename)
end

function neosmarginalhist(αs::Vector{T}, δs::Vector{T}, filename::String;
    sub::String = "(a)", rmax::T = 6.0, hmax::T = 30_000.0,
    bins::AbstractRange{T} = -6:0.4:6, sticks::AbstractRange{S} = -6:2:6,
    hticks::AbstractRange{U} = 6_000:6_000:30_000,
    slabel::AbstractVector{T} = [-5, 5, -5.5, 32_500, 35_000, -5.75]) where {T <: Real,
    S <: Real, U <: Real}
    layout = @layout [a              _
        b{0.75w, 0.75h} c]
    plot(layout = layout, link = :both, size = (600, 600), legend = false,
        guidefont = (14, :black), tickfont = (12, :black))
    scatter!(αs, δs, subplot = 2, framestyle = :box,
        xlims = (-rmax, rmax), ylims = (-rmax, rmax), alpha = 0.25,
        xticks = (sticks, string.(sticks)), yticks = (sticks, string.(sticks)),
        markerstrokewidth = 0,
        xlabel = L"\alpha\cos(\delta) \ \textrm{[arcsec]}",
        ylabel = L"\delta \ \textrm{[arcsec]}")
    annotate!(slabel[1], slabel[2], sub, subplot = 2)
    histogram!(αs, subplot = 1, orientation = :v,
        framestyle = :box, legend = false, margin = -12mm,
        bins = bins, alpha = 0.5,
        xlim = (-rmax, rmax), ylim = (0, hmax),
        xticks = (sticks, fill("", length(sticks))),
        yticks = (hticks, string.(hticks .÷ 1_000)),
        ylabel = "Counts",
        top_margin = 3mm)
    annotate!([(slabel[3], slabel[4], Plots.text(L"\times10^{3}", 11, :black, :center))],
        subplot = 1)
    histogram!(δs, subplot = 3, orientation = :h,
        framestyle = :box, legend = false, margin = -12mm,
        bins = bins, alpha = 0.5,
        ylim = (-rmax, rmax), xlim = (0, hmax),
        yticks = (sticks, fill("", length(sticks))),
        xticks = (hticks, string.(hticks .÷ 1_000)),
        xlabel = "Counts",
        right_margin = 7mm)
    annotate!([(slabel[5], slabel[6], Plots.text(L"\times10^{3}", 11, :black, :center))],
        subplot = 3)
    savefig(filename)
end

function main()
    # Initial time
    init_time = now()
    # Parse arguments from commandline
    parsed_args = parse_commandline()
    # Input directories
    input1::String = parsed_args["input1"]
    input2::String = parsed_args["input2"]
    # Print header
    println("Reproduce the plots and tables in Ramírez-Montoya et al, 2024, section 5.2.")
    println("• Input directory for sample 1: ", input1)
    println("• Input directory for sample 2: ", input2)

    # Parse NEOs'
    neos1 = readlines(joinpath(input1, "names.txt"))
    neos2 = readlines(joinpath(input2, "names.txt"))
    # Load orbits
    sols1 = map(neos1) do neo
        filename = joinpath(input1, replace(neo, " " => "") * ".jld2")
        !isfile(filename) && return zero(NEOSolution{Float64, Float64})
        return JLD2.load(filename, "sol")
    end
    sols2 = map(neos2) do neo
        filename = joinpath(input2, replace(neo, " " => "") * ".jld2")
        !isfile(filename) && return zero(NEOSolution{Float64, Float64})
        return JLD2.load(filename, "sol")
    end

    # NRMS
    Q1, Q2 = @. nrms(sols1), nrms(sols2)
    # Chi-squared critical values
    C1, C2 = @. critical_value(sols1), critical_value(sols2)
    # SNR
    S1, S2 = @. minimum(snr(sols1)), minimum(snr(sols2))
    # NEOs that converged within 0.99 confidence
    mask1 = @. (C1 ≤ 0.99) && !isinf(S1)
    mask2 = @. (C2 ≤ 0.99) && !isinf(S2)
    # Right ascension and declination residuals
    αs1 = reduce(vcat, map(sols1[mask1]) do s
        m = @. !isoutlier(s.res)
        return @. ra(s.res[m])
    end)
    δs1 = reduce(vcat, map(sols1[mask1]) do s
        m = @. !isoutlier(s.res)
        return @. dec(s.res[m])
    end)
    αs2 = reduce(vcat, map(sols2[mask2]) do s
        m = @. !isoutlier(s.res)
        return @. ra(s.res[m])
    end)
    δs2 = reduce(vcat, map(sols2[mask2]) do s
        m = @. !isoutlier(s.res)
        return @. dec(s.res[m])
    end)

    # Convergence statistics
    println("• Residuals statistics: ")
    p1 = @sprintf("%.2f", count(mask1) * 100 / length(neos1))
    p2 = @sprintf("%.2f", count(mask2) * 100 / length(neos2))
    l1 = "$(count(mask1))/$(length(neos1)) NEOs ($p1%)"
    l2 = "$(count(mask2))/$(length(neos2)) NEOs ($p2%)"
    println("NEOs within 99% confidence in sample 1: ", l1)
    println("NEOs within 99% confidence in sample 2: ", l2)
    mask1S, mask2S = @. mask1 && S1 < 1, mask2 && S2 < 1
    p1 = @sprintf("%.2f", count(mask1S) * 100 / count(mask1)) * "%"
    p2 = @sprintf("%.2f", count(mask2S) * 100 / count(mask2)) * "%"
    println("In sample 1: ", p1, " of the NEOs within 99% confidence",
        " have a minimal SNR lower that 1")
    println("In sample 2: ", p2, " of the NEOs within 99% confidence",
        " have a minimal SNR lower that 1")
    m1 = @sprintf("%.2f", mean(numberofdays(s.tracklets) for s in sols1[mask1]))
    m1S = @sprintf("%.2f", mean(numberofdays(s.tracklets) for s in sols1[mask1S]))
    m2 = @sprintf("%.2f", mean(numberofdays(s.tracklets) for s in sols2[mask2]))
    m2S = @sprintf("%.2f", mean(numberofdays(s.tracklets) for s in sols2[mask2S]))
    println("On average, the observations used in sample 1 span a time of ",
        m1, " days for the NEOs under 99% confidence, but only ", m1S,
        " days for the NEOs with a minimal SNR smaller than 1")
    println("On average, the observations used in sample 2 span a time of ",
        m2, " days for the NEOs under 99% confidence, but only ", m2S,
        " days for the NEOs with a minimal SNR smaller than 1")
    # Residuals statistics
    println("• Residuals statistics: ")
    pretty_table(
    [
        "R.A."  "Median"    median(αs1)     median(αs2)
        ""      "Mean"      mean(αs1)       mean(αs2)
        ""      "Variance"  var(αs1)        var(αs2)
        ""      "Skewness"  skewness(αs1)   skewness(αs2)
        ""      "Kurtosis"  kurtosis(αs1)   kurtosis(αs2)
        "Dec."  "Median"    median(δs1)     median(δs2)
        ""      "Mean"      mean(δs1)       mean(δs2)
        ""      "Variance"  var(δs1)        var(δs2)
        ""      "Skewness"  skewness(δs1)   skewness(δs2)
        ""      "Kurtosis"  kurtosis(δs1)   kurtosis(δs2)
    ],
    header = ["", "", "Sample 1", "Sample 2"], formatters = ft_printf("%.3e", 3:4)
    )

    # Plots.jl setting
    default(fontfamily = "Computer Modern", titlefont = (14),
        guidefont = (12, :black), tickfont = (10, :black),
        framestyle = :box, yminorgrid = true, xminorgrid = true,
        margin = 0.5Measures.mm, dpi = 300, legend = true)
    # NRMS histogram
    filename1, filename2 = joinpath(input1, "QHIST.pdf"), joinpath(input2, "QHIST.pdf")
    qhist(Q1, mask1, filename1; ymax = 600.0, ysub = 550.0, sub = "(a)")
    qhist(Q2, mask2, filename2; ymax = 200.0, ysub = 182.5, sub = "(b)")
    println("• Saved NRMS histograms to: \n", filename1, "\n", filename2)
    # SNR histogram
    filename1, filename2 = joinpath(input1, "SHIST.pdf"), joinpath(input2, "SHIST.pdf")
    shist(S1, mask1, filename1; ymax = 450.0, ysub = 425.0, sub = "(a)")
    shist(S2, mask2, filename2; ymax = 150.0, ysub = 142.5, sub = "(b)")
    println("• Saved SNR histograms to: \n", filename1, "\n", filename2)
    # Residuals marginal histograms
    filename1, filename2 = joinpath(input1, "RHIST.pdf"), joinpath(input2, "RHIST.pdf")
    neosmarginalhist(αs1, δs1, filename1;
        sub = "(a)", rmax = 6.0, hmax = 30_000.0, bins = -6:0.4:6,
        sticks = -6.0:2.0:6.0, hticks = 6_000:6_000:30_000,
        slabel = [-5, 5, -5.5, 32_500, 35_000, -5.75])
    neosmarginalhist(αs2, δs2, filename2;
        sub = "(b)", rmax = 5.0, hmax = 10_000.0, bins = -5:0.5:5,
        sticks = -5.0:2.5:5.0, hticks = 2_000:2_000:10_000,
        slabel = [-4, 4.15, -4.5, 10_850, 11_600, -4.75])
    println("• Saved residuals marginal histograms to: \n", filename1, "\n", filename2)

    # Final time
    final_time = now()
    println("• Run started ", init_time, " and finished ", final_time)

    return nothing
end

main()