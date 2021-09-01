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
objname = "Bennu"
maxsteps = 10000
nyears_bwd = -13.0 # years in backward integration
nyears_fwd = 10.0 # years in forward integration
dense = false #true
quadmath = false # use quadruple precision
lyap = false #true # compute Lyapunov exponents
neosjlpath = pkgdir(NEOs)
dynamics = RNp1BP_pN_A_J23E_J2S_ng_eph_threads!

# integration parameters
order = 25
abstol = 1.0E-20

### observation data files (ra/dec, del/dop)
optfile_bwd = joinpath(neosjlpath, "data", "101955_OPTICAL_1999_2006.dat")
optfile_fwd = joinpath(neosjlpath, "data", "101955_OPTICAL_2011_2018.dat")
radarfile_bwd = joinpath(neosjlpath, "data", "101955_RADAR_1999_2005.dat")
radarfile_fwd = joinpath(neosjlpath, "data", "101955_RADAR_2011.dat")

### path to local Solar System ephemeris file
ss_eph_file = "./sseph343ast016_p56y_et.jld"
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

### initial conditions from Bennu JPL solution \#97 ### Note: JPL sol \#97 uses nongrav function g(r)=(1 au)/r^2.25, while we use g(r)=(1 au)/r^2
jd0 = datetime2julian(DateTime(2011,1,1)) # JDTDB = 2455562.5
q00 = [-1.1951358208617802, -0.20726185835689961, -0.11201678544935807, 8.881637772597003e-5, -0.013056288090844732, -0.007377624521045638]
q0 = vcat(q00, 0.0) .+ dq ####vcat(q00, -4.614425051841E-14) #.+ dq

####integrator warmup
propagate(objname, dynamics, 1, jd0, nyears_fwd, ss16asteph_et, output=false, dense=dense, q0=q0, quadmath=quadmath, lyap=lyap, order=order, abstol=abstol)
println("*** Finished warmup")

######Full jet transport integration
propagate(objname*"_bwd", dynamics, maxsteps, jd0, nyears_bwd, ss16asteph_et, dense=dense, q0=q0, quadmath=quadmath, lyap=lyap, order=order, abstol=abstol, radarobsfile=radarfile_bwd, opticalobsfile=optfile_bwd, tord=10, niter=5)
propagate(objname*"_fwd", dynamics, maxsteps, jd0, nyears_fwd, ss16asteph_et, dense=dense, q0=q0, quadmath=quadmath, lyap=lyap, order=order, abstol=abstol, radarobsfile=radarfile_fwd, opticalobsfile=optfile_fwd, tord=10, niter=5)
println("*** Finished asteroid ephemeris integration")
