# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using Dates
using PlanetaryEphemeris
using TaylorIntegration
using Plots
using Test

const NEOs_DATA = joinpath(pkgdir(NEOs), "data")
const TEST_DATA = joinpath(pkgdir(NEOs), "test", "data")

# Load optical astrometry
obs99942 = read_optical_mpc80(joinpath(NEOs_DATA, "99942_2004_2020.dat"))
filter!(x -> Date(2005, 1, 27) < date(x) < Date(2005, 1, 31), obs99942)
obs895907 = read_optical_mpc80(joinpath(TEST_DATA, "895907.txt"))
filter!(x -> Date(2016, 1, 9) < date(x) < Date(2016, 1, 20), obs895907)
# Reduce optical tracklets
trks99942 = reduce_tracklets(obs99942)
trks895907 = reduce_tracklets(obs895907)
# Parameters
params = Parameters(
    coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
    gaussorder = 2, safegauss = false,
    mmovorder = 2, mmoviter = 500, mmovQtol = 1e-5, jtlsorder = 2,
    jtlsmask = false, jtlsiter = 20, lsiter = 10, significance = 0.99,
    outrej = true, χ2_rec = 7.0, χ2_rej = 8.0,
    fudge = 100.0, max_per = 34.0,
)
# Preliminary orbit
loadjpleph()
od = ODProblem(newtonian!, obs99942)
jd0 = datetime2julian(DateTime(2005, 1, 29))
q00 = kmsec2auday(apophisposvel199(julian2etsecs(jd0)))
orbit = LeastSquaresOrbit(od, q00, jd0, params)

@testset "RecipesBase.jl extension" begin

    @testset "OpticalResidual" begin
        @test plot(
            orbit.ores, xlabel = "αcos(δ)", ylabel = "δ",
            xlim = (-2.5, 2.5), ylim = (-2.5, 2.5),
            xticks = -2.5:0.5:2.5, yticks = -2.5:0.5:2.5,
            aspect_ratio = 1, label = "", framestyle = :box
        ) isa Plots.Plot
    end

    @testset "AdmissibleRegion" begin

        # Common keyword arguments
        N = 1_000
        framestyle = :box
        Hs = vcat(34.5, 32:-2:14)

        @testset "One component" begin
            A = AdmissibleRegion(trks99942[1], params)
            @test_throws AssertionError plot(A, N = 0)
            @test_throws AssertionError plot(A, ρscale = :invalid)
            @test plot(
                A; ρscale = :log, N, Hs, framestyle,
                xlabel = "log₁₀(ρ)", ylabel = "v_ρ",
                xlim = (-4, 0), ylim = (-0.01, 0.04),
                xticks = -4:0, yticks = -0.01:0.01:0.04
            ) isa Plots.Plot
            @test plot(
                A; ρscale = :linear, N, Hs, framestyle,
                xlabel = "ρ", ylabel = "v_ρ",
                xlim = (-0.01, 1.0), ylim = (-0.01, 0.04),
                xticks = 0:0.1:1.0, yticks = -0.01:0.01:0.04
            ) isa Plots.Plot
        end

        @testset "Two components" begin
            A = AdmissibleRegion(trks895907[1], params)
            @test_throws AssertionError plot(A, N = 0)
            @test_throws AssertionError plot(A, ρscale = :invalid)
            @test plot(
                A; ρscale = :log, N, Hs, framestyle,
                xlabel = "log₁₀(ρ)", ylabel = "v_ρ",
                xlim = (-3, 2), ylim = (-0.03, 0.02),
                xticks = -3:2, yticks = -0.03:0.01:0.02
            ) isa Plots.Plot
            @test plot(
                A; ρscale = :linear, N, Hs, framestyle,
                xlabel = "ρ", ylabel = "v_ρ",
                xlim = (-1, 50), ylim = (-0.03, 0.02),
                xticks = 0:10:50, yticks = -0.03:0.01:0.02
            ) isa Plots.Plot
        end

    end

    @testset "TaylorSolution / AbstractOrbit" begin
        N = 1_000
        t0, tf = lasttime(orbit.bwd), lasttime(orbit.fwd)
        projections = (:x, :y, :z, :xy, :xz, :yz, :xyz)
        @test_throws AssertionError plot(orbit, t0, tf, N = 0)
        @test_throws AssertionError plot(orbit, t0, tf, projection = :invalid)
        for projection in projections
            @test plot(orbit, t0, tf; N, projection) isa Plots.Plot
        end
    end
end