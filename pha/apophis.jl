# This file is part of the NEOs.jl package; MIT licensed

### This script can be run either as a standalone script or via ArgParse. In any case,
### this folder's Project.toml environment has to be already active and instantiated
### before running this script. Uncomment the following three lines to activate and
### instantiate this folder's environment:

# import Pkg
# Pkg.activate(".")
# Pkg.instantiate()

### If ran as a standalone and the environment defined in this folder's Project.toml is
### already active and instantiated, then this script can be run from this folder with the
### default settings simply as:
### $ julia -t <number-of-threads> --project=. apophis.jl
### Finally, this script can be run via the ArgParse mechanism. Help can be displayed doing:
### $ julia --project=. apophis.jl --help

using ArgParse
using NEOs
using Dates
using TaylorIntegration
using JLD2
using PlanetaryEphemeris
using DataFrames
using DelimitedFiles
using LinearAlgebra: diag
using Statistics
using StatsBase

# Load JPL ephemeris
loadjpleph()

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "apophis.jl"
    # Desciption (for help screen)
    s.description = "Propagates Apophis orbit via jet transport"

    @add_arg_table! s begin
        "--jd0"
            help = "Initial date"
            arg_type = DateTime
            default = DateTime(2020, 12, 17)
        "--varorder"
            help = "Order of the jet transport perturbation"
            arg_type = Int
            default = 5
        "--maxsteps"
            help = "Maximum number of steps during integration"
            arg_type = Int
            default = 10_000
        "--nyears_bwd"
            help = "Years in backward integration"
            arg_type = Float64
            default = -17.0
        "--nyears_fwd"
            help = "Years in forward integration"
            arg_type = Float64
            default = 9.0
        "--order"
            help = "Order of Taylor polynomials expansions during integration"
            arg_type = Int
            default = 25
        "--abstol"
            help = "Absolute tolerance"
            arg_type = Float64
            default = 1.0E-20
        "--parse_eqs"
            help = "Whether to use the taylorized method of jetcoeffs or not"
            arg_type = Bool
            default = true
    end

    s.epilog = """
        Examples (run from the `pha` folder):\n
        \n
        # Multi-threaded\n
        julia -t 4 --project=. apophis.jl --maxsteps 100 --nyears_bwd -0.02 --nyears_fwd 0.02 --parse_eqs true\n
        \n
        # Single-threaded\n
        julia --project=. apophis.jl --maxsteps 100 --nyears_bwd -0.02 --nyears_fwd 0.02 --parse_eqs true\n
        \n
    """

    return parse_args(s)
end

function print_header(header::String, level::Int = 1)
    L = length(header)
    if level == 1
        c = "="
    else
        c = "-"
    end
    println(repeat(c, L))
    println(header)
    println(repeat(c, L))
end

function main(dynamics::D, maxsteps::Int, jd0_datetime::DateTime, nyears_bwd::T, nyears_fwd::T, order::Int, varorder::Int,
              abstol::T, parse_eqs::Bool) where {T <: Real, D}

    # Initial conditions from Apophis JPL solution #197
    q00 = kmsec2auday(apophisposvel197(datetime2et(jd0_datetime)))

    # Perturbation to nominal initial condition (Taylor1 jet transport)
    # vcat(fill(1e-8, 6), 1e-14, 1e-15) are the scaling factors for jet transport perturbation,
    # these are needed to ensure expansion coefficients remain small.
    # The magnitudes correspond to the typical order of magnitude of errors in
    # position/velocity (1e-8), Yarkovsky (1e-14) and radiation pressure (1e-15)
    scalings = vcat(fill(1e-8, 6), 1e-14, 1e-15)
    if varorder == 0
        dq = zeros(8)
    else
        dq = NEOs.scaled_variables("δx", scalings, order = varorder)
    end

    q0 = vcat(q00, 0.0, 0.0) .+ dq

    # Initial date (in Julian days)
    jd0 = datetime2julian(jd0_datetime)

    print_header("Integrator warmup", 2)
    _ = NEOs.propagate(dynamics, 1, jd0, nyears_fwd, q0, Val(true);
                       order, abstol, parse_eqs)
    _ = NEOs.propagate_root(dynamics, 1, jd0, nyears_fwd, q0, Val(true);
                       order, abstol, parse_eqs)
    println()

    print_header("Main integration", 2)
    tmax = nyears_bwd*yr
    println("• Initial time of integration: ", string(jd0_datetime))
    println("• Final time of integration: ", julian2datetime(jd0 + tmax))

    sol_bwd = NEOs.propagate(dynamics, maxsteps, jd0, nyears_bwd, q0, Val(true);
                             order, abstol, parse_eqs)
    jldsave("Apophis_bwd.jld2"; sol_bwd)
    # sol_bwd = JLD2.load("Apophis_bwd.jld2", "sol_bwd")

    tmax = nyears_fwd*yr
    println("• Initial time of integration: ", string(jd0_datetime))
    println("• Final time of integration: ", julian2datetime(jd0 + tmax))

    sol_fwd, tvS, xvS, gvS = NEOs.propagate_root(dynamics, maxsteps, jd0, nyears_fwd, q0, Val(true);
                             order, abstol, parse_eqs)
    jldsave("Apophis_fwd.jld2"; sol_fwd, tvS, xvS, gvS)
    # sol_fwd = JLD2.load("Apophis_fwd.jld2", "sol_bwd")
    println()

    # load Solar System ephemeris
    sseph::TaylorInterpolant{Float64,Float64,2} = loadpeeph(NEOs.sseph, sol_bwd.t0+sol_bwd.t[end], sol_fwd.t0+sol_fwd.t[end])
    eph_su::TaylorInterpolant{Float64,Float64,2} = selecteph(sseph, su)
    eph_ea::TaylorInterpolant{Float64,Float64,2} = selecteph(sseph, ea)

    # NEO
    # Change t, x, v units, resp., from days, au, au/day to sec, km, km/sec
    xva_bwd(et) = auday2kmsec(sol_bwd(et/daysec)[1:6])
    xva_fwd(et) = auday2kmsec(sol_fwd(et/daysec)[1:6])
    xva(et) = bwdfwdeph(et, sol_bwd, sol_fwd)
    # Earth
    # Change x, v units, resp., from au, au/day to km, km/sec
    xve(et) = auday2kmsec(eph_ea(et/daysec))
    # Sun
    # Change x, v units, resp., from au, au/day to km, km/sec
    xvs(et) = auday2kmsec(eph_su(et/daysec))


    radec_2004_2020 = read_radec_mpc(joinpath(pkgdir(NEOs), "data", "99942_2004_2020.dat"))
    radec_2020_2021 = read_radec_mpc(joinpath(pkgdir(NEOs), "data", "99942_2020_2021.dat"))
    radec = vcat(radec_2004_2020,radec_2020_2021)

    deldop_2005_2013 = read_radar_jpl(joinpath(pkgdir(NEOs), "data", "99942_RADAR_2005_2013.dat"))
    deldop_2021 = read_radar_jpl(joinpath(pkgdir(NEOs), "data", "99942_RADAR_2021.dat"))
    deldop = vcat(deldop_2005_2013,deldop_2021)

    # Compute optical residuals
    res_radec_all, w_radec_all = NEOs.residuals(radec; xvs, xve, xva)
    jldsave("Apophis_res_w_radec.jld2"; res_radec_all, w_radec_all)
    # JLD2.@load "Apophis_res_w_radec.jld2"

    # Compute radar residuals
    res_del, w_del, res_dop, w_dop = NEOs.residuals(deldop; xvs, xve, xva, niter=10, tord=10)
    jldsave("Apophis_res_w_deldop.jld2"; res_del, w_del, res_dop, w_dop)
    # JLD2.@load "Apophis_res_w_deldop.jld2"

    ### Process optical astrometry (filter, weight, debias)

    # filter out biased observations from observatory 217 on 28-Jan-2021
    df_radec = DataFrame(radec)
    # add residuals and weights to optical astrometry DataFrame
    df_radec[!, :res_α] .= res_radec_all[1:round(Int,length(res_radec_all)/2)]
    df_radec[!, :res_δ] .= res_radec_all[1+round(Int,length(res_radec_all)/2):end]
    df_radec[!, :w_α] .= w_radec_all[1:round(Int,length(res_radec_all)/2)]
    df_radec[!, :w_δ] .= w_radec_all[1+round(Int,length(res_radec_all)/2):end]
    filter!(
        x->(Date(x.date) != Date(2021, 1, 28)),
        df_radec
    )

    # read astrometric errors from Tholen et al. (2013)
    tho13_errors = readdlm(joinpath(pkgdir(NEOs), "data", "tholenetal2013_opterror.dat"), ',')
    # compute weights
    w_α_tho13 = 1 ./ (tho13_errors[:,1].^2 .+ tho13_errors[:,3].^2 .+ tho13_errors[:,5].^2)
    w_δ_tho13 = 1 ./ (tho13_errors[:,2].^2 .+ tho13_errors[:,4].^2 .+ tho13_errors[:,6].^2)
    # Tholen et al. (2013) obs table
    radec_tho13 = DataFrame(read_radec_mpc(joinpath(pkgdir(NEOs), "test", "data", "99942_Tholen_etal_2013.dat")))
    # vector of RA values from Tholen et al. (2013) observations (used for filtering)
    tho13_α = radec_tho13[!,:α]
    # set weights in Tholen et al. (2013) astrometry corresponding to associated uncertainties
    df_radec[in.(df_radec.α, Ref(tho13_α)),:w_α] = w_α_tho13
    df_radec[in.(df_radec.α, Ref(tho13_α)),:w_δ] = w_δ_tho13

    # Relaxation factor (account for correlations in optical astrometry data)
    # for each observation batch, count the number of observations made in
    # the same night by the same observatory
    # Ref: Veres et al. (2017)
    date_site_v = select(df_radec, :date => ByRow(Date), :observatory)
    Nv = [count(x->x.date_Date==i.date_Date && x.observatory==i.observatory, eachrow(date_site_v)) for i in eachrow(date_site_v)]
    relax_factor = map(x->x>4.0 ? x/4.0 : 1.0, Nv)
    # inflate uncertainties (i.e., relax weights) by relaxation factor
    df_radec[!, :w_α] .= (df_radec.w_α)./relax_factor
    df_radec[!, :w_δ] .= (df_radec.w_δ)./relax_factor

    # update optical residuals and weights
    res_radec = vcat(df_radec.res_α, df_radec.res_δ)
    w_radec = vcat(df_radec.w_α, df_radec.w_δ)

    ### Perform orbital fit to optical and radar astrometry data

    res = vcat(res_radec, res_del, res_dop)
    w = vcat(w_radec, w_del, w_dop)

    success, δx_OR8, Γ_OR8 = newtonls(res, w, zeros(get_numvars()), 10)
    x_OR8 = sol_fwd(sol_fwd.t0)(δx_OR8)
    σ_OR8 = sqrt.(diag(Γ_OR8)).*scalings

    nradec = length(res_radec)
    res_ra = view(res_radec, 1:nradec÷2)
    res_dec = view(res_radec, 1+nradec÷2:nradec)
    w_ra = view(w_radec, 1:nradec÷2)
    w_dec = view(w_radec, 1+nradec÷2:nradec)

    ### Print results

    print_header("Orbital fit (8-DOF) and post-fit statistics", 2)

    # orbital fit
    println("Success flag                                                      : ", success, "\n")
    println("Nominal solution             [au,au,au,au/d,au/d,au/d,au/d²,au/d²]: ", x_OR8, "\n")
    println("1-sigma formal uncertainties [au,au,au,au/d,au/d,au/d,au/d²,au/d²]: ", σ_OR8, "\n")

    # post-fit statistics
    println("Normalized RMS (optical-only)               [adimensional] : ", nrms(res_radec(δx_OR8),w_radec))
    println("Normalized RMS (radar-only)                 [adimensional] : ", nrms(vcat(res_del,res_dop)(δx_OR8),vcat(w_del,w_dop)))
    println("Normalized RMS (combined optical and radar) [adimensional] : ", nrms(res(δx_OR8),w), "\n")
    println("Mean weighted right-ascension residual      [arcseconds]   : ", mean(res_ra(δx_OR8), weights(w_ra)))
    println("Mean weighted declination residual          [arcseconds]   : ", mean(res_dec(δx_OR8), weights(w_dec)))
    println("Mean weighted time-delay residual           [micro-seconds]: ", mean(res_del(δx_OR8), weights(w_del)))
    println("Mean weighted Doppler-shift residual        [Hz]           : ", mean(res_dop(δx_OR8), weights(w_dop)), "\n")
    println("Chi-squared statistic (χ²):                 [adimensional] : ", chi2(res(δx_OR8),w))

    return sol_bwd, sol_fwd, res_radec, res_del, res_dop, w_radec, w_del, w_dop
end

function main()

    # Parse arguments from commandline
    parsed_args = parse_commandline()

    print_header("Asteroid Apophis")
    println()
    print_header("General parameters", 2)

    # Number of threads
    N_threads = Threads.nthreads()
    println("• Number of threads: ", N_threads)

    # Dynamical function
    dynamics = RNp1BP_pN_A_J23E_J2S_ng_eph_threads!

    println("• Dynamical function: ", dynamics)

    # Maximum number of steps
    maxsteps = parsed_args["maxsteps"]
    println("• Maximum number of steps: ", maxsteps)

    # Initial date
    jd0_datetime = parsed_args["jd0"]

    # Number of years in backward integration
    nyears_bwd = parsed_args["nyears_bwd"]

    # Number of years in forward integration
    nyears_fwd = parsed_args["nyears_fwd"]

    # Order of Taylor polynomials
    order = parsed_args["order"]
    println("• Order of Taylor polynomials: ", order)

    # Order of jet transport perturbation
    varorder = parsed_args["varorder"]
    println("• Order of jet transport perturbation: ", varorder)

    # Absolute tolerance
    abstol = parsed_args["abstol"]
    println("• Absolute tolerance: ", abstol)

    # Whether to use @taylorize or not
    parse_eqs = parsed_args["parse_eqs"]
    println("• Use @taylorize: ", parse_eqs, "\n")

    main(dynamics, maxsteps, jd0_datetime,  nyears_bwd, nyears_fwd, order, varorder, abstol, parse_eqs)
end

main()
