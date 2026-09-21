# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using Dates
using PlanetaryEphemeris
using TaylorIntegration
using Plots
using Test

# Load optical astrometry
optical = read_optical_mpc80(joinpath(pkgdir(NEOs), "data",
    "99942_2004_2020.dat"))
filter!(x -> Date(2005, 1, 27) < date(x) < Date(2005, 1, 31), optical)
# Parameters
params = Parameters(
    coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
    gaussorder = 2, safegauss = false,
    tsaorder = 2, adamiter = 500, adamQtol = 1e-5, jtlsorder = 2,
    jtlsmask = false, jtlsiter = 20, lsiter = 10, significance = 0.99,
    outrej = true, χ2_rec = 7.0, χ2_rej = 8.0,
    fudge = 100.0, max_per = 34.0,
)
# Orbit determination problem (only optical astrometry)
od = ODProblem(newtonian!, optical)
# Admissible region
A = AdmissibleRegion(od.tracklets[1], params)
# Preliminary orbit
loadjpleph()
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
        N = 1_000
        framestyle = :box
        Hs = vcat(34.5, 32:-2:14)
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