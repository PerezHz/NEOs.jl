#julia --project=@. main.jl
using Apophis
using Dates

#script parameters (TODO: use ArgParse.jl instead)
const objname = "Apophis"
const objdata = "Apophis_JPL_data.dat"
const newtoniter = 5
const maxsteps = 10000
const nyears = 24.0 # since t0 is 2008-9-24, this ends the integration on 2032-9-24; NOTE: this value is overriden when evaluating solution at JPL radar observation times
const radarobs = false#true
const jt = true#false#true
const dynamics = RNp1BP_pN_A_J234E_J2S_ng!
const t0 = datetime2julian(DateTime(2008,9,24,0,0,0)) #starting time of integration
@show t0 == 2454733.5

#integrator warmup
main(objname, objdata, dynamics, 2, newtoniter, t0, nyears, output=false, radarobs=radarobs, jt=jt)
println("*** Finished warmup")

#root-finding methods warmup (integrate until first root-finding event):
#main(objname, objdata, dynamics, 110, newtoniter, t0, nyears, radarobs=radarobs, jt=jt)
#println("*** Finished root-finding warmup")

#main(objname, objdata, dynamics, 300, newtoniter, t0, nyears, radarobs=radarobs, jt=jt)
#println("*** Finished root-finding test: several roots")

#Full jet transport integration until ~2038: about 8,000 steps
# main(objname, objdata, dynamics, 8000, newtoniter, t0, nyears, radarobs=radarobs, jt=jt)
# println("*** Finished full jet transport integration")
