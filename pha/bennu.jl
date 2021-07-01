#Multi-threaded:
#julia -t <number-of-threads> --project=@. main.jl
#Single thread:
#julia --project=@. main.jl
using NEOs
using Dates
using TaylorIntegration
using JLD
using PlanetaryEphemeris
@show Threads.nthreads()

#script parameters (TODO: use ArgParse.jl instead)
const varorder = 1 #10 #5 #1 # varorder is the order corresponding to the jet transport perturbation
const nv = 7 #1 #number of TaylorN variables
const objname = "Bennu" ##"Apophis"
const maxsteps = 10000
const nyears = 10.0 #21.0 #-13.0 #39.5 #-18.0 #-18.0 #-5.0 #6.0 #21.0
const dense = false #true
const quadmath = false # use quadruple precision
const lyap = false #true # compute Lyapunov exponents
const neosjlpath = pkgdir(NEOs)
const dynamics = RNp1BP_pN_A_J23E_J2S_ng_eph_threads!

# integration parameters
const order = 25
const abstol = 1.0E-20

#### observation data files (ra/dec, del/dop)
const opticalobsfile = joinpath(neosjlpath, "data", "vokr15_mpc_formatted.dat") # "99942_2004_2020.dat")
const radarobsfile = joinpath(neosjlpath, "data", "99942_RADAR_2005_2013.dat")
# const opticalobsfile = joinpath(neosjlpath, "data", "99942_2020_2021.dat")
# const radarobsfile = joinpath(neosjlpath, "data", "99942_RADAR_2021.dat")

# path to local Solar System ephemeris file
#ss_eph_file = "../../PlanetaryEphemeris/notebooks/sseph343ast016_p70y_et_2000_2069.jld"
#ss_eph_file = joinpath(neosjlpath, "jldeph", "sseph343ast016_p30y_et_J2000.jld")
ss_eph_file = joinpath(neosjlpath, "jldeph", "sseph343ast016_m30y_et_2030.jld")



#### dq: perturbation to nominal initial condition (Taylor1 jet transport)
#dq = Taylor1.(zeros(7), varorder)
#dq[end][1] = 1e-14

# # dq: perturbation to nominal initial condition (TaylorN jet transport)
# dq = set_variables("δx", order=varorder, numvars=nv)
# for i in 1:6
#    dq[i][1][i] = 1e-8
# end
# if get_numvars() == 7
#    dq[7][1][7] = 1e-14
# end

#### WARNING: When calling with lyap=true, remember to setup TaylorN variables with order=1, numvars=7
#TNvars = set_variables("δx", order=1, numvars=7)
#dq=zeros(7)

### OR7 solution, post-2021
# jd0 = datetime2julian(DateTime(2020,12,17)) ###Julian date (TDB) of integration initial time
# q00 = [-0.18034747703273316, 0.9406910666200128, 0.3457360259054398, -0.016265942170279046, 4.392889725556651e-5, -0.00039519931615139716] ### s197 2020Dec17.0 (TDB)
# q0 = (vcat(q00, 0.0) .+ [-8.072551353250423e-7, -6.977674993203065e-9, -3.535944557003007e-8, 2.3820741189757113e-9, -1.340336958871951e-8, -4.6952597219089e-9, -2.895160537809007e-14]) .+ dq

### OR7 solution, pre-2021, equivalent to Vokrouhlicky et al (2015) solution
# jd0 = datetime2julian(DateTime(2008,9,24)) ###Julian date (TDB) of integration initial time
# q00 = x0_JPL_s197
# q0 = (vcat(q00, 0.0) .+ [7.16926378354607e-8, 1.4192184039812996e-7, 5.186569758679862

### BENNU, initial conditions from JPL solution \#97
### Note: JPL sol \#97 uses nongrav function g(r)=(1 au)/r^2.25, while we use g(r)=(1 au)/r^2
jd0 = datetime2julian(DateTime(2011,1,1)) # JDTDB = 2455562.5
q00 = [-1.1951358208617802, -0.20726185835689961, -0.11201678544935807, 8.881637772597003e-5, -0.013056288090844732, -0.007377624521045638]
q0 = vcat(q00, -4.614425051841E-14) #.+ dq

####integrator warmup
propagate(objname, dynamics, 1, jd0, nyears, ss_eph_file, output=false, dense=dense, q0=q0, quadmath=quadmath, lyap=lyap, order=order, abstol=abstol)
println("*** Finished warmup")

#propagate(objname, dynamics, 5 #=300=#, jd0, nyears, ss_eph_file, dense=dense, q0=q0, quadmath=quadmath, lyap=lyap, order=order, abstol=abstol)
#println("*** Finished 5 steps")

######Full jet transport integration until ~2038: about 8,000 steps
##### no obs ephemeris computation
propagate(objname, dynamics, maxsteps, jd0, nyears, ss_eph_file, dense=dense, q0=q0, quadmath=quadmath, lyap=lyap, order=order, abstol=abstol)

##### with obs ephemeris computation
#propagate(objname, dynamics, maxsteps, jd0, nyears, ss_eph_file, dense=dense, q0=q0, quadmath=quadmath, order=order, abstol=abstol, radarobsfile=radarobsfile, opticalobsfile=opticalobsfile, tord=10, niter=5)
println("*** Finished full jet transport integration")


########


# astfname = "Bennu_jt.jld"
# asteph = JLD.load(astfname, "asteph")
# ss16asteph_et = JLD.load(ss_eph_file, "ss16ast_eph")
# NEOs.loadjpleph()
# # NEOs.compute_radar_obs("deldop_2005_2013_o$(varorder)v$(nv).jdb", "../notebooks/99942_RADAR_2005_2013.dat", asteph, ss16asteph_et)
# # NEOs.compute_optical_obs("radec_2004_2020_o$(varorder)v$(nv).jdb", "99942_2004_2020.dat", asteph, ss16asteph_et, debias_table=debias_table)

# NEOs.compute_radar_obs("NEW_deldop_2005_2013_o$(varorder)v$(nv).jdb", "../data/99942_RADAR_2005_2013.dat", asteph, ss16asteph_et, tord=10)
# #NEOs.compute_optical_obs("NEW_radec_2004_2020_o$(varorder)v$(nv).jdb", "../data/vokr15_mpc_formatted.dat", asteph, ss16asteph_et, debias_table=debias_table)
