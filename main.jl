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
const varorder = 5 # varorder is the order corresponding to the jet transport perturbation
const nv = 7 #number of TaylorN variables
const objname = "Apophis"
const maxsteps = 10000
const nyears = 5.0 #-5.0 #24.0
const dense = true #false
const apophisjlpath = dirname(pathof(Apophis))
const radarobsfile = joinpath(apophisjlpath, "../Apophis_JPL_data_2012_2013.dat")
# const radarobsfile = joinpath(apophisjlpath, "../Apophis_JPL_data_2005_2006.dat")
# const dynamics = RNp1BP_pN_A_J23E_J2S_ng_eph!
const dynamics = RNp1BP_pN_A_J23E_J2S_ng_eph_threads!
const t0 = datetime2julian(DateTime(2008,9,24,0,0,0)) #starting time of integration
@show t0 == 2454733.5

# path to local Solar System ephemeris file
# ss_eph_file = joinpath(apophisjlpath, "../jpleph", "ss16ast343_eph_24yr_tx.jld")
ss_eph_file = joinpath(apophisjlpath, "../jpleph", "ss16ast343_eph_5yr_tx.jld")
# ss_eph_file = joinpath(apophisjlpath, "../jpleph", "ss16ast343_eph_minus5yr_tx_BACKWARDS.jld")

# dq: perturbation to nominal initial condition (Taylor1 jet transport)
dq = Taylor1.(zeros(7), varorder)
dq[end][1] = 1e-14

# # dq: perturbation to nominal initial condition (TaylorN jet transport)
# dq = set_variables("ξ", order=varorder, numvars=nv)
# for i in 1:6
#     dq[i][1][i] = 1e-8
# end
# if get_numvars() == 7
#     dq[7][1][7] = 1e-14
# end

#integrator warmup
propagate(objname, dynamics, 1, t0, nyears, ss_eph_file, output=false, dense=dense, dq=dq)
println("*** Finished warmup")

propagate(objname, dynamics, 5, t0, nyears, ss_eph_file, dense=dense, dq=dq)
println("*** Finished 2nd warmup")

#root-finding methods warmup (integrate until first root-finding event):
# propagate(objname, dynamics, 50, t0, nyears, ss_eph_file, dense=dense, dq=dq)
# println("*** Finished root-finding warmup")

#propagate(objname, dynamics, 100, t0, nyears, ss_eph_file, dense=dense, dq=dq)
#println("*** Finished root-finding test: several roots")

#Full jet transport integration until ~2038: about 8,000 steps
#propagate(objname, dynamics, maxsteps, t0, nyears, ss_eph_file, dense=dense, dq=dq)
#println("*** Finished full jet transport integration")

# # calculate computed values of time-delays and Doppler shifts
# ss16asteph, acc_eph, newtonianNb_Potential = Apophis.loadeph(ss_eph_file)
# astfname = "Apophis_jt.0.jld"
# t = load(astfname, "t")
# x = load(astfname, "x")
# tx = TaylorInterpolant(t, x)
# furnsh( joinpath(apophisjlpath, "../jpleph", "naif0012.tls") ) # load leapseconds kernel
# furnsh( joinpath(apophisjlpath, "../jpleph", "de431t.bsp") ) # at least one SPK file must be loaded to read .tls file
# Apophis.compute_radar_obs("deldop.jld", radarobsfile, tx, ss16asteph)
