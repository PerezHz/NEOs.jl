#Multi-threaded:
#JULIA_NUM_THREADS=8 julia --project=@. main.jl
#Single-threaded:
#julia --project=@. main.jl
using Apophis
using Dates
using TaylorSeries

@show Threads.nthreads()

#script parameters (TODO: use ArgParse.jl instead)
const objname = "Apophis"
const maxsteps = 10000
const nyears = 24.0
const dense = true#false
# const dynamics = RNp1BP_pN_A_J23E_J2S_ng_eph!
const dynamics = RNp1BP_pN_A_J23E_J2S_ng_eph_threads!
const t0 = datetime2julian(DateTime(2008,9,24,0,0,0)) #starting time of integration
@show t0 == 2454733.5

# path to local Solar System ephemeris file
my_eph_file = joinpath(dirname(pathof(Apophis)), "../jpleph", "ss16ast343_eph_24yr_tx.jld")
# my_eph_file = joinpath(dirname(pathof(Apophis)), "../jpleph", "ss16ast343_eph_5yr_tx.jld")

varorder = 4 # varorder is the order corresponding to the jet transport perturbation
# dq: perturbation to nominal initial condition (Taylor1 jet transport)
# dq = Taylor1.(zeros(7), varorder)
# dq[end][1] = 1e-14

# dq: perturbation to nominal initial condition (TaylorN jet transport)
# ξv = set_variables("ξ", order=varorder, numvars=1)
# zeroxi = zero(ξv[1])
# dq = [zeroxi, zeroxi, zeroxi, zeroxi, zeroxi, zeroxi, 1e-14ξv[1]]
dq = set_variables("ξ", order=varorder, numvars=7)

#integrator warmup
propagate(objname, dynamics, 1, t0, nyears, my_eph_file, output=false, dense=dense, dq=dq)
println("*** Finished warmup")

propagate(objname, dynamics, 2, t0, nyears, my_eph_file, dense=dense, dq=dq)
println("*** Finished 2nd warmup")

#root-finding methods warmup (integrate until first root-finding event):
# propagate(objname, dynamics, 50, t0, nyears, my_eph_file, dense=dense, dq=dq)
# println("*** Finished root-finding warmup")

#propagate(objname, dynamics, 100, t0, nyears, my_eph_file, dense=dense, dq=dq)
#println("*** Finished root-finding test: several roots")

#Full jet transport integration until ~2038: about 8,000 steps
# propagate(objname, dynamics, 8000, t0, nyears, my_eph_file, dense=dense, dq=dq)
# println("*** Finished full jet transport integration")
