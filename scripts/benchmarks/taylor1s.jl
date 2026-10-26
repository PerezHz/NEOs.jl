# Taylor1S vs nested Taylor1 benchmarks: (1) in-place kernels, (2) closeapproaches
#
# Example:
#   julia --project taylor1s.jl -i ../../risklist/2022UE61/2022UE61Orbit.jld2
using ArgParse, NEOs, PlanetaryEphemeris, TaylorSeries, Dates, JLD2, Printf, BenchmarkTools
using NEOs: CloseApproachesBuffer, Taylor1S, initialcondition, radius

function parse_commandline()
    s = ArgParseSettings()
    s.prog = "taylor1s.jl"
    s.description = "Taylor1S vs nested Taylor1 benchmarks"
    @add_arg_table! s begin
        "--input", "-i"
            help = "input orbit file (JLD2 with key \"orbit\"); if absent, only \
                the kernels are benchmarked"
            arg_type = String
            default = ""
        "--nva"
            help = "number of virtual asteroids to propagate"
            arg_type = Int
            default = 3
        "--ctol"
            help = "convergence tolerance of closeapproaches"
            arg_type = Float64
            default = 1.0
        "--samples", "-s"
            help = "number of samples per closeapproaches benchmark"
            arg_type = Int
            default = 5
    end
    return parse_args(s)
end

# Copy of a Taylor1 as a Taylor1S of the same order
t1s(x::Taylor1) = Taylor1S(Tuple(x.coeffs))

# Deterministic test series: nested Taylor1 (time order n, σ order m) and its Taylor1S copy
function testseries(n::Int, m::Int, s::Int, a0::Float64)
    inner(k) = Taylor1([cos(s * k + j) / (1 + k + j) for j in 0:m], m)
    x = Taylor1([inner(k) + a0 * (k == 0) for k in 0:n], n)
    X = Taylor1(t1s.(x.coeffs), n)
    return x, X
end

# Run kernel `f!(c, x, y, aux, k)` for all orders
function runkernel!(f!::F, c, x, y, aux) where {F}
    for k in eachindex(c)
        f!(c, x, y, aux, k)
    end
    return nothing
end

const KERNELS = (
    ("mul!",          (c, x, y, aux, k) -> TaylorSeries.mul!(c, x, y, k)),
    ("pow! (r = 2)",  (c, x, y, aux, k) -> TaylorSeries.pow!(c, x, aux, 2, k)),
    ("div!",          (c, x, y, aux, k) -> TaylorSeries.div!(c, x, y, k)),
    ("div! (scalar)", (c, x, y, aux, k) -> TaylorSeries.div!(c, 0.7, y, k)),
    ("sqrt!",         (c, x, y, aux, k) -> TaylorSeries.sqrt!(c, x, aux, k)),
    ("pow! (r = -2)", (c, x, y, aux, k) -> TaylorSeries.pow!(c, x, aux, -2.0, k)),
    ("pow! (r = 1.5)",(c, x, y, aux, k) -> TaylorSeries.pow!(c, x, aux, 1.5, k)),
)

function benchmark_kernels(; n::Int = 20, ms = (2, 4, 6, 8))
    println("\n# Kernels (time order $n)\n")
    @printf("| %-15s | %2s | %12s | %13s | %8s | %5s | %9s |\n", "kernel", "m",
            "nested (ns)", "Taylor1S (ns)", "speed-up", "==", "isequal")
    println("|", "-"^17, "|", "-"^4, "|", "-"^14, "|", "-"^15, "|", "-"^10, "|",
            "-"^7, "|", "-"^11, "|")
    for m in ms
        T1S = Taylor1S{m + 1, Float64}
        x, X = testseries(n, m, 3, 2.0)
        y, Y = testseries(n, m, 5, 1.5)
        for (name, f!) in KERNELS
            cN, auxN, cS, auxS = zero(x), zero(x), zero(X), zero(X)
            runkernel!(f!, cN, x, y, auxN)
            runkernel!(f!, cS, X, Y, auxS)
            eq = all(T1S(cN[k]) == cS[k] for k in eachindex(cS))
            ieq = all(isequal(T1S(cN[k]), cS[k]) for k in eachindex(cS))
            tN = @belapsed runkernel!($f!, $cN, $x, $y, $auxN)
            tS = @belapsed runkernel!($f!, $cS, $X, $Y, $auxS)
            @printf("| %-15s | %2d | %12.1f | %13.1f | %8.2f | %5s | %9s |\n",
                    name, m, 1E9 * tN, 1E9 * tS, tN / tS, eq, ieq)
        end
        # Horner evaluation of a state vector at a time step
        qN, qS = fill(x, 6), fill(X, 6)
        dN, dS = [zero(x[0]) for _ in 1:6], [zero(X[0]) for _ in 1:6]
        TaylorSeries.evaluate!(qN, 0.25, dN)
        TaylorSeries.evaluate!(qS, 0.25, dS)
        eq, ieq = all(T1S.(dN) .== dS), all(isequal.(T1S.(dN), dS))
        tN = @belapsed TaylorSeries.evaluate!($qN, 0.25, $dN)
        tS = @belapsed TaylorSeries.evaluate!($qS, 0.25, $dS)
        @printf("| %-15s | %2d | %12.1f | %13.1f | %8.2f | %5s | %9s |\n",
                "evaluate! (x6)", m, 1E9 * tN, 1E9 * tS, tN / tS, eq, ieq)
    end
    return nothing
end

# Maximum absolute difference between the coefficients of a Taylor1S close approach `x`
# and a nested Taylor1 close approach `y`
cadiff(x, y) = maximum(f -> maximum(abs, (getfield(x, f) - t1s(getfield(y, f))).c),
                       (:t, :a, :x, :y, :z))

function benchmark_closeapproaches(input::String, nva::Int, ctol::Float64, samples::Int)
    # Same setup as the impact monitoring of 2022 UE61
    σmax, lovorder, lovtol, lovsteps = 5.0, 4, 1E-8, 1_000
    IP, Δσmax, vaorder = 1E-6, 0.01, 4
    params = Parameters(
        maxsteps = 20_000, order = 20, abstol = eps(), parse_eqs = true,
        coeffstol = Inf, bwdoffset = 0.05, fwdoffset = 0.05,
        gaussorder = 2, safegauss = true, refscale = :log,
        mmovorder = 2, mmoviter = 500, mmovQtol = 1E-5,
        jtlsorder = 2, jtlsmask = false, jtlsiter = 20, lsiter = 10,
        jtlsproject = true, significance = 0.99, verbose = true,
        outrej = false, χ2_rec = 9.21, χ2_rej = 10.0, fudge = 100.0,
        max_per = 33.3, threads = false
    )
    orbit = JLD2.load(input, "orbit")
    IM = IMProblem(orbit, ImpactTarget(:earth))
    lov = lineofvariations(IM, params; coord = :cartesian, σmax, lovorder, lovtol, lovsteps)
    R_TP, R_P = 0.2, radius(IM.target)
    VAs = virtualasteroids(lov, :DelVigna19, vaorder; R_TP, R_P, IP, Δσmax)
    jd0 = epoch(lov) + PE.J2000
    nyears = (datetime2julian(DateTime(2100, 1, 1, 12)) - jd0) / yr
    q0 = orbit() .+ 1E-8 * Taylor1(vaorder)
    nested = CloseApproachesBuffer(IM, q0, nyears, params)
    t1sbuf = CloseApproachesBuffer(IM, t1s.(q0), nyears, params)
    @assert nested.prop.cache.parse_eqs === t1sbuf.prop.cache.parse_eqs === Val(true)
    idxs = nva == 1 ? [1] : unique(round.(Int, range(1, length(VAs), length = nva)))

    println("\n# closeapproaches ($(length(VAs)) VAs, ctol = $ctol, ",
            "order = $(params.order), vaorder = $vaorder)\n")
    @printf("| %5s | %4s | %5s | %9s | %11s | %12s | %8s | %12s | %14s |\n",
            "VA", "CAs", "==", "max diff", "nested (s)", "Taylor1S (s)", "speed-up",
            "nested (MiB)", "Taylor1S (MiB)")
    println("|", "-"^7, "|", "-"^6, "|", "-"^7, "|", "-"^11, "|", "-"^13, "|", "-"^14,
            "|", "-"^10, "|", "-"^14, "|", "-"^16, "|")
    for i in idxs
        VA = VAs[i]
        VAS = VirtualAsteroid(epoch(VA), sigma(VA), VA.domain, t1s.(initialcondition(VA)))
        CAsN = closeapproaches(IM, VA, nyears, params; R_TP, ctol, buffer = nested)
        CAsS = closeapproaches(IM, VAS, nyears, params; R_TP, ctol, buffer = t1sbuf)
        d = length(CAsN) == length(CAsS) ?
            maximum(cadiff.(CAsS, CAsN); init = 0.0) : NaN
        eq = iszero(d)
        bN = @benchmark closeapproaches($IM, $VA, $nyears, $params; R_TP = $R_TP,
            ctol = $ctol, buffer = $nested) samples = samples evals = 1 seconds = 1E4
        bS = @benchmark closeapproaches($IM, $VAS, $nyears, $params; R_TP = $R_TP,
            ctol = $ctol, buffer = $t1sbuf) samples = samples evals = 1 seconds = 1E4
        tN, tS = 1E-9 * minimum(bN).time, 1E-9 * minimum(bS).time
        @printf("| %5d | %4d | %5s | %9.2e | %11.3f | %12.3f | %8.2f | %12.1f | %14.1f |\n",
                i, length(CAsN), eq, d, tN, tS, tN / tS, bN.memory / 2^20, bS.memory / 2^20)
    end
    return nothing
end

function main()
    args = parse_commandline()
    println("• Run started at ", now(), " with ", Threads.nthreads(), " thread(s)")
    benchmark_kernels()
    isempty(args["input"]) || benchmark_closeapproaches(args["input"], args["nva"],
        args["ctol"], args["samples"])
    println("\n• Run finished at ", now())
    return nothing
end

main()
