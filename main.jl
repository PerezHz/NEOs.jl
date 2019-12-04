#Multi-threaded:
#JULIA_NUM_THREADS=8 julia --project=@. main.jl
#Single-threaded:
#julia --project=@. main.jl
using Apophis
using Dates
using TaylorSeries

@show Threads.nthreads()

#script parameters (TODO: use ArgParse.jl instead)
const varorder = 5 # varorder is the order corresponding to the jet transport perturbation
const objname = "Apophis"
const maxsteps = 10000
const nyears = 5.0
const dense = true#false
const apophisjlpath = dirname(pathof(Apophis))
const radarobsfile = joinpath(apophisjlpath, "../Apophis_JPL_data_2012_2013.dat")
# const dynamics = RNp1BP_pN_A_J23E_J2S_ng_eph!
const dynamics = RNp1BP_pN_A_J23E_J2S_ng_eph_threads!
const t0 = datetime2julian(DateTime(2008,9,24,0,0,0)) #starting time of integration
@show t0 == 2454733.5

# path to local Solar System ephemeris file
# ast_eph_file = joinpath(dirname(pathof(Apophis)), "../jpleph", "ss16ast343_eph_24yr_tx.jld")
ast_eph_file = joinpath(dirname(pathof(Apophis)), "../jpleph", "ss16ast343_eph_5yr_tx.jld")

# dq: perturbation to nominal initial condition (Taylor1 jet transport)
dq = Taylor1.(zeros(7), varorder)
dq[end][1] = 1e-14

# dq: perturbation to nominal initial condition (TaylorN jet transport)
# dq = set_variables("ξ", order=varorder, numvars=7)
# for i in 1:6
#     dq[i][1][i] = 1e-8
# end
# dq[end][1][end] = 1e-14

#integrator warmup
propagate(objname, dynamics, 1, t0, nyears, ast_eph_file, output=false, dense=dense, dq=dq)
println("*** Finished warmup")

#propagate(objname, dynamics, 2, t0, nyears, ast_eph_file, dense=dense, dq=dq)
#println("*** Finished 2nd warmup")

#root-finding methods warmup (integrate until first root-finding event):
# propagate(objname, dynamics, 50, t0, nyears, ast_eph_file, dense=dense, dq=dq)
# println("*** Finished root-finding warmup")

#propagate(objname, dynamics, 100, t0, nyears, ast_eph_file, dense=dense, dq=dq)
#println("*** Finished root-finding test: several roots")

#Full jet transport integration until ~2038: about 8,000 steps
propagate(objname, dynamics, maxsteps, t0, nyears, ast_eph_file, dense=dense, dq=dq, radarobsfile=radarobsfile)
println("*** Finished full jet transport integration")
