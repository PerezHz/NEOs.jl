import Pkg
Pkg.activate("./")
using Apophis
using Dates

#script parameters (TODO: use ArgParse.jl instead)
const objname = "Apophis"
const objdata = "Apophis_JPL_data.dat"
const newtoniter = 5
const maxsteps = 10000
const nyears = 24.0 # since t0 is 2008-9-24, this ends the integration on 2032-9-24; NOTE: this value is overriden when evaluating solution at JPL radar observation times
const radarobs = false#true
const jt = false#true
const t0 = datetime2julian(DateTime(2008,9,24,0,0,0)) #starting time of integration
@show t0 == 2454733.5

#integrator warmup
main(objname, objdata, 1, newtoniter, t0, nyears, output=false, radarobs=radarobs, jt=jt)
println("*** Finished warmup")

#root-finding methods warmup (integrate until first root-finding event):
# main(objname, objdata, 100, newtoniter, t0, nyears, radarobs=radarobs, jt=jt)
# println("*** Finished root-finding warmup")

main(objname, objdata, 300, newtoniter, t0, nyears, radarobs=radarobs, jt=jt)
println("*** Finished root-finding test: several roots")

#Full jet transport integration until ~2038: about 8,000 steps
# main(objname, objdata, 8000, newtoniter, t0, nyears, radarobs=radarobs, jt=jt)
# println("*** Finished full jet transport integration")
