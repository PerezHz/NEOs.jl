#Multi-threaded:
#julia -t <number-of-threads> --project=@. main.jl
#Single thread:
#julia --project=@. main.jl
using Apophis
using Dates
using TaylorIntegration
using JLD
using SPICE: furnsh
@show Threads.nthreads()

#script parameters (TODO: use ArgParse.jl instead)
const varorder = 5 # 1 # varorder is the order corresponding to the jet transport perturbation
const nv = 7 #number of TaylorN variables
const objname = "Apophis"
const maxsteps = 10000
const nyears = 6.0 #-5.0 #21.0
const dense = false #true
const quadmath = false # use quadruple precision
const debias_table = "2018" # "2014", "hires2018"
const apophisjlpath = pkgdir(Apophis)
# const dynamics = RNp1BP_pN_A_J23E_J2S_ng_eph!
const dynamics = RNp1BP_pN_A_J23E_J2S_ng_eph_threads!
const jd0 = datetime2julian(DateTime(2008,9,24,0,0,0)) #Julian date of integration initial time
@show jd0 == 2454733.5
const t0 = 0.0 # integration initial time

#### observation data files (ra/dec, del/dop)
#const opticalobsfile = joinpath(apophisjlpath, "tholen13_mpc_formatted.dat")
#const radarobsfile = joinpath(apophisjlpath, "Apophis_JPL_data_2005_2006.dat")
const opticalobsfile = "" #joinpath(apophisjlpath, "vokr15_mpc_formatted.dat") # #""
const radarobsfile = "" #joinpath(apophisjlpath, "Apophis_JPL_data_2012_2013.dat") #""

# path to local Solar System ephemeris file
#ss_eph_file = joinpath(apophisjlpath, "jldeph", "ss16ast343_eph_m5y_et.jld")
ss_eph_file = joinpath(apophisjlpath, "jldeph", "ss16ast343_eph_p6y_et.jld")

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

#propagate(objname, dynamics, 300 #=5=#, t0, nyears, ss_eph_file, dense=dense, dq=dq, quadmath=quadmath)
#println("*** Finished 5 steps")

######Full jet transport integration until ~2038: about 8,000 steps
###propagate(objname, dynamics, maxsteps, t0, nyears, ss_eph_file, dense=dense, dq=dq) # no obs ephemeris computation
propagate(objname, dynamics, maxsteps, t0, nyears, ss_eph_file, dense=dense, dq=dq, quadmath=quadmath, radarobsfile=radarobsfile, opticalobsfile=opticalobsfile, debias_table=debias_table)
println("*** Finished full jet transport integration")





#const ss_eph_file_m5y = joinpath(apophisjlpath, "jldeph", "ss16ast343_eph_m5y_et_f64.jld")
#const ss_eph_file_p6y = joinpath(apophisjlpath, "jldeph", "ss16ast343_eph_p6y_et_f64.jld")
#const astfname_m5y = "Apophis_jt_m5y_o5v7.jld"
#const astfname_p6y = "Apophis_jt_p6y_o5v7.jld"
#ss16asteph_m5y = load(ss_eph_file_m5y, "ss16ast_eph")
#ss16asteph_p6y = load(ss_eph_file_p6y, "ss16ast_eph")
#apophis_m5y = load(astfname_m5y, "apophis")
#apophis_p6y = load(astfname_p6y, "apophis")
#furnsh( joinpath(apophisjlpath, "jldeph", "naif0012.tls") ) # load leapseconds kernel
#furnsh( joinpath(apophisjlpath, "jldeph", "de430_1850-2150.bsp") ) # at least one SPK file must be loaded to read .tls file
#Apophis.compute_radar_obs("deldop_m5y_o5v7.jld", radarobsfile0506, apophis_m5y, ss16asteph_m5y)
#Apophis.compute_radar_obs("deldop_p6y_o5v7.jld", radarobsfile1213, apophis_p6y, ss16asteph_p6y)
#Apophis.compute_optical_obs("radec_tho13_o5v7.jld", tho13obsfile, apophis_m5y, ss16asteph_m5y)
#Apophis.compute_optical_obs_v15("radec_vokr15_o5v7.jld", vokr15obsfile, apophis_p6y, ss16asteph_p6y)
