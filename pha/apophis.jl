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
varorder = 5 # varorder is the order corresponding to the jet transport perturbation
nv = 7 #number of TaylorN variables
objname = "Apophis"
maxsteps = 10000
nyears_bwd = -18.0 # years in backward integration
nyears_fwd = 9.0 # years in forward integration
dense = false #true
quadmath = false # use quadruple precision
lyap = false #true # compute Lyapunov exponents
neosjlpath = pkgdir(NEOs)
dynamics = RNp1BP_pN_A_J23E_J2S_ng_eph_threads!

# integration parameters
order = 25
abstol = 1.0E-20

### observation data files (ra/dec, del/dop)
optfile_bwd = joinpath(neosjlpath, "data", "99942_2004_2020.dat")
optfile_fwd = joinpath(neosjlpath, "data", "99942_2020_2021.dat")
radarfile_bwd = joinpath(neosjlpath, "data", "99942_RADAR_2005_2013.dat")
radarfile_fwd = joinpath(neosjlpath, "data", "99942_RADAR_2021.dat")

### path to local Solar System ephemeris file
ss_eph_file = "./sseph343ast016_p31y_et.jld"
ss16asteph_et = JLD.load(ss_eph_file, "ss16ast_eph")

### TaylorN variables setup
if lyap
    ### setup TaylorN variables with order=1, numvars=7
    TNvars = set_variables("δx", order=1, numvars=7)
    dq=zeros(7)
else
    #### dq: perturbation to nominal initial condition (Taylor1 jet transport)
    #dq = Taylor1.(zeros(7), varorder)
    #dq[end][1] = 1e-14
    # dq: perturbation to nominal initial condition (TaylorN jet transport)
    dq = set_variables("δx", order=varorder, numvars=nv)
    for i in 1:6
        dq[i][1][i] = 1e-8
    end
    if get_numvars() == 7
        dq[7][1][7] = 1e-14
    end
end

### initial conditions from Apophis JPL solution #197 at the Dec-17-2020.0 (TDB) epoch
jd0 = datetime2julian(DateTime(2020,12,17)) # JDTDB = 2459200.5
q00 = [-0.18034747703273316, 0.9406910666200128, 0.3457360259054398, -0.016265942170279046, 4.392889725556651e-5, -0.00039519931615139716] ### JPL solution #197 at 2020Dec17.0 (TDB)
q0 = vcat(q00, 0.0) .+ dq

####integrator warmup
propagate(objname, dynamics, 1, jd0, nyears_fwd, ss16asteph_et, output=false, dense=dense, q0=q0, quadmath=quadmath, lyap=lyap, order=order, abstol=abstol)
println("*** Finished warmup")

######Full jet transport integration
propagate(objname*"_bwd", dynamics, maxsteps, jd0, nyears_bwd, ss16asteph_et, dense=dense, q0=q0, quadmath=quadmath, lyap=lyap, order=order, abstol=abstol, radarobsfile=radarfile_bwd, opticalobsfile=optfile_bwd, tord=10, niter=5)
propagate(objname*"_fwd", dynamics, maxsteps, jd0, nyears_fwd, ss16asteph_et, dense=dense, q0=q0, quadmath=quadmath, lyap=lyap, order=order, abstol=abstol, radarobsfile=radarfile_fwd, opticalobsfile=optfile_fwd, tord=10, niter=5)
println("*** Finished asteroid ephemeris integration")
