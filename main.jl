#Multi-threaded:
#JULIA_NUM_THREADS=<number-of-threads> julia --project=@. main.jl
#Single-threaded:
#julia --project=@. main.jl
using Apophis
using Dates
using TaylorIntegration
using JLD
using SPICE: furnsh
@show Threads.nthreads()

#script parameters (TODO: use ArgParse.jl instead)
const varorder = 10 # varorder is the order corresponding to the jet transport perturbation
const nv = 7 #number of TaylorN variables
const objname = "Apophis"
const maxsteps = 10000
const nyears = 5.0 #24.0
const dense = false #true
const quadmath = true # use quadruple precision
const apophisjlpath = pkgdir(Apophis)
const radarobsfile = joinpath(apophisjlpath, "Apophis_JPL_data_2012_2013.dat")
# const radarobsfile = joinpath(apophisjlpath, "Apophis_JPL_data_2005_2006.dat")
# const dynamics = RNp1BP_pN_A_J23E_J2S_ng_eph!
const dynamics = RNp1BP_pN_A_J23E_J2S_ng_eph_threads!
const jd0 = datetime2julian(DateTime(2008,9,24,0,0,0)) #Julian date of integration initial time
@show jd0 == 2454733.5
const t0 = 0.0 # integration initial time

# path to local Solar System ephemeris file
ss_eph_file = joinpath(apophisjlpath, "jldeph", "ss16ast343_eph_p5y_et.jld")
#ss_eph_file = joinpath(apophisjlpath, "jldeph", "ss16ast343_eph_p24y_et.jld")
#ss_eph_file = joinpath(apophisjlpath, "jldeph", "ss16ast343_eph_m4y_et.jld")

#### dq: perturbation to nominal initial condition (Taylor1 jet transport)
dq = Taylor1.(zeros(7), varorder)
dq[end][1] = 1e-14

#### dq: perturbation to nominal initial condition (TaylorN jet transport)
# dq = set_variables("ξ", order=varorder, numvars=nv)
# for i in 1:6
#     dq[i][1][i] = 1e-8
# end
# if get_numvars() == 7
#     dq[7][1][7] = 1e-14
# end

####integrator warmup
propagate(objname, dynamics, 1, t0, nyears, ss_eph_file, output=false, dense=dense, dq=dq, quadmath=quadmath)
println("*** Finished warmup")

propagate(objname, dynamics, 5, t0, nyears, ss_eph_file, dense=dense, dq=dq, quadmath=quadmath)
println("*** Finished 2nd warmup with output")

####Full jet transport integration until ~2038: about 8,000 steps
##propagate(objname, dynamics, maxsteps, t0, nyears, ss_eph_file, dense=dense, dq=dq, quadmath=quadmath)
##println("*** Finished full jet transport integration")

# #### calculate computed values of time-delays and Doppler shifts
# ss16asteph = load(ss_eph_file, "ss16ast_eph")
# astfname = "Apophis_jt.0.jld"
# tx = load(astfname, "apophis")
# furnsh( joinpath(apophisjlpath, "jldeph", "naif0012.tls") ) # load leapseconds kernel
# furnsh( joinpath(apophisjlpath, "jldeph", "de430_1850-2150.bsp") ) # at least one SPK file must be loaded to read .tls file
# Apophis.compute_radar_obs("deldop.jld", radarobsfile, tx, ss16asteph)
