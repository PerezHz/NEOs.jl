include("src/Apophis.jl")
using .Apophis

#script parameters (TODO: use ArgParse.jl instead)
const newtoniter = 5
const maxsteps = 10000
const nyears = 24.0 # since t0 is 2008-9-24, this ends the integration on 2032-9-24; NOTE: this value is overriden when evaluating solution at NEODyS radar observation times
const radarobs = true

#integrator warmup
main(2, newtoniter, nyears, output=false, radarobs=radarobs)
# main(2, newtoniter)
println("*** Finished warmup")

#root-finding methods warmup (integrate until first root-finding event):
main(100, newtoniter, nyears, radarobs=radarobs)
println("*** Finished root-finding warmup")

#main(300, newtoniter, nyears, radarobs=radarobs)
#println("*** Finished root-finding test: several roots")

#Full jet transport integration until ~2038: about 8,000 steps
# main(8000, newtoniter, nyears, radarobs=radarobs)
# println("*** Finished full jet transport integration")
