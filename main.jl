include("src/Apophis.jl")
using .Apophis

#script parameters (TODO: use ArgParse.jl instead)
const newtoniter = 5
const maxsteps = 10000

#integrator warmup
main(2, newtoniter, false)
# main(2, newtoniter, true)
println("*** Finished warmup")

#root-finding methods warmup (integrate until first root-finding event):
main(100, newtoniter, true)
println("*** Finished root-finding warmup")

#main(300, newtoniter, true)
#println("*** Finished root-finding test: several roots")

#Full jet transport integration until ~2038: about 8,000 steps
# main(8000, newtoniter, true)
# println("*** Finished full jet transport integration")
