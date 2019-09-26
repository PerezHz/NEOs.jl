#Multi-threaded:
#JULIA_NUM_THREADS=8 julia --project=@. main.jl
#Single-threaded:
#julia --project=@. main.jl
using Apophis
using Dates

@show Threads.nthreads()

#script parameters (TODO: use ArgParse.jl instead)
const objname = "Apophis"
const objdata = "../Apophis_JPL_data.dat"
const newtoniter = 10
const maxsteps = 10000
const nyears = 24.0 # since t0 is 2008-9-24, this ends the integration on 2032-9-24; NOTE: this value is overriden when evaluating solution at JPL radar observation times
const jt = true
const dense = true#false
# const dynamics = RNp1BP_pN_A_J23E_J2S_ng_eph!
const dynamics = RNp1BP_pN_A_J23E_J2S_ng_eph_threads!
const t0 = datetime2julian(DateTime(2008,9,24,0,0,0)) #starting time of integration
@show t0 == 2454733.5

#integrator warmup
propagate(objname, objdata, dynamics, 1, newtoniter, t0, nyears, output=false, jt=jt, dense=dense)
println("*** Finished warmup")

propagate(objname, objdata, dynamics, 2, newtoniter, t0, nyears, jt=jt, dense=dense)
println("*** Finished 2nd warmup")

#root-finding methods warmup (integrate until first root-finding event):
# propagate(objname, objdata, dynamics, 50, newtoniter, t0, nyears, jt=jt, dense=dense)
# println("*** Finished root-finding warmup")

#propagate(objname, objdata, dynamics, 100, newtoniter, t0, nyears, jt=jt, dense=dense)
#println("*** Finished root-finding test: several roots")

#Full jet transport integration until ~2038: about 8,000 steps
# propagate(objname, objdata, dynamics, 8000, newtoniter, t0, nyears, jt=jt, dense=dense)
# println("*** Finished full jet transport integration")
