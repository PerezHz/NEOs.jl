#numerator of Apophis radial velocity wrt Earth
function rvelea(dx, x, params, t)
    ss16asteph_t = params[1](t) #ss16asteph(t)
    xe = ss16asteph_t[union(3ea-2:3ea,3(N-1+ea)-2:3(N-1+ea))]
    return (x[1]-xe[1])*(x[4]-xe[4]) + (x[2]-xe[2])*(x[5]-xe[5]) + (x[3]-xe[3])*(x[6]-xe[6])
end

function propagate(objname::String, dynamics::Function, maxsteps::Int, t0::T,
        tspan::T, ephfile::String; output::Bool=true, newtoniter::Int=10,
        dense::Bool=false, dq::Vector=zeros(7), radarobsfile::String="") where {T<:Real}

    # read Solar System ephemeris (Sun+8 planets+Moon+Pluto+16 main belt asteroids)
    # ephfile = "ss16ast343_eph_24yr_tx.jld"
    # ss16ast_eph_t = load(joinpath(jplephpath, ephfile), "ss16ast_eph_t")
    # ss16ast_eph_x = load(joinpath(jplephpath, ephfile), "ss16ast_eph_x")
    # ss16ast_eph_t = load(ephfile, "t")
    # ss16ast_eph_x = load(ephfile, "x")
    ss16ast_eph_t = load(ephfile, "ss16ast_eph_t")
    ss16ast_eph_x = load(ephfile, "ss16ast_eph_x")
    ss16asteph = TaylorInterpolant(ss16ast_eph_t, ss16ast_eph_x)
    #compute point-mass Newtonian accelerations from ephemeris: all bodies except Apophis
    # accelerations of "everybody else" are needed when evaluating Apophis post-Newtonian acceleration
    Nm1 = N-1
    acc_eph = TaylorInterpolant(ss16ast_eph_t, Matrix{eltype(ss16ast_eph_x)}(undef, length(ss16ast_eph_t)-1, 3Nm1))
    newtonianNb_Potential = TaylorInterpolant(ss16ast_eph_t, Matrix{eltype(ss16ast_eph_x)}(undef, length(ss16ast_eph_t)-1, Nm1))
    fill!(acc_eph.x, zero(ss16ast_eph_x[1]))
    fill!(newtonianNb_Potential.x, zero(ss16ast_eph_x[1]))
    _1_to_Nm1 = Base.OneTo(Nm1) # iterator over all bodies except Apophis
    for j in _1_to_Nm1
        for i in _1_to_Nm1
            # i == j && continue
            if i == j
            else
                X_ij = ss16ast_eph_x[:,3i-2] .- ss16ast_eph_x[:,3j-2]
                Y_ij = ss16ast_eph_x[:,3i-1] .- ss16ast_eph_x[:,3j-1]
                Z_ij = ss16ast_eph_x[:,3i  ] .- ss16ast_eph_x[:,3j  ]
                r_p2_ij = ( (X_ij.^2) .+ (Y_ij.^2) ) .+ (Z_ij.^2)
                r_p3d2_ij = r_p2_ij.^1.5
                r_ij = sqrt.(r_p2_ij)
                newtonianCoeff_ij =  μ[i]./r_p3d2_ij
                acc_eph.x[:,3j-2] .+= (X_ij.*newtonianCoeff_ij)
                acc_eph.x[:,3j-1] .+= (Y_ij.*newtonianCoeff_ij)
                acc_eph.x[:,3j  ] .+= (Z_ij.*newtonianCoeff_ij)
                newtonianNb_Potential.x[:,j] .+= (μ[i]./r_ij)
            end #if i != j
        end #for, i
    end #for, j
    params = (ss16asteph, acc_eph, newtonianNb_Potential)
    # get asteroid initial conditions
    q0 = initialcond(dq)

    @show tmax = t0+tspan*yr #final time of integration

    # do integration
    if dense
        if radarobsfile != ""
            asteroid_data = process_radar_data_jpl(radarobsfile)
            # TODO: check that first and last observation times are within interpolation interval
            @time interp = apophisinteg(dynamics, q0, t0, tmax, order, abstol, params; maxsteps=maxsteps, dense=dense)
            # Load TT-TDB DE430 ephemeris
            furnsh( joinpath(jplephpath, "TTmTDB.de430.19feb2015.bsp") )
            function apophis_et(et)
                return interp( etsecs2julian(et) )[1:6]
            end
            function earth_et(et)
                return ss16asteph( etsecs2julian(et) )[union(3*4-2:3*4,3*(27+4)-2:3*(27+4))]
            end
            function sun_et(et)
                return ss16asteph( etsecs2julian(et) )[union(3*1-2:3*1,3*(27+1)-2:3*(27+1))]
            end
            #compute time-delay and Doppler-shift "ephemeris" (i.e., predicted values according to ephemeris)
            vdel, vdop = delay_doppler(asteroid_data; xve=earth_et, xvs=sun_et, xva=apophis_et)
            sol = (t=interp.t[:], x=interp.x[:,:], vdel=vdel, vdop=vdop)
        else
            @time interp = apophisinteg(dynamics, q0, t0, tmax, order, abstol, params; maxsteps=maxsteps, dense=dense)
            sol = (t=interp.t[:], x=interp.x[:,:])
        end
    else
        @time sol_objs = apophisinteg(dynamics, rvelea, q0, t0, tmax, order, abstol, params; maxsteps=maxsteps, newtoniter=newtoniter)
        sol = (
            tv1 = sol_objs[1][:],
            xv1 = sol_objs[2][:,:],
            tvS1 = sol_objs[3][:],
            xvS1 = sol_objs[4][:,:],
            gvS1 = sol_objs[5][:]
        )
    end

    #write solution to .jld files
    if output
        filename = string(objname, "_jt.jld")
        println("Saving solution to file: $filename")
        #first, deal with `tv_jpl_integ`
        jldopen(filename, "w") do file
            #loop over solution variables
            for ind in eachindex(sol)
                varname = string(ind)
                println("Saving variable: ", varname)
                write(file, varname, sol[ind])
            end
        end
        println("Checking that all variables were saved correctly...")
        #loop over solution variables
        for ind in eachindex(sol)
            varname = string(ind)
            #read varname from files and assign recovered variable to recovered_sol_i
            recovered_sol_i = load(filename, varname)
            #check that varname was recovered succesfully
            @show recovered_sol_i == sol[ind]
        end
        println("Saved solution")
    end
    return nothing
end

function testjetcoeffs()
    # read Solar System ephemeris (Sun+8 planets+Moon+Pluto+16 main belt asteroids)
    ephfile = "ss16ast343_eph_24yr_tx.jld"
    ss16ast_eph_t = load(joinpath(jplephpath, ephfile), "ss16ast_eph_t")
    ss16ast_eph_x = load(joinpath(jplephpath, ephfile), "ss16ast_eph_x")
    ss16asteph = TaylorInterpolant(ss16ast_eph_t, ss16ast_eph_x)
    #compute point-mass Newtonian accelerations from ephemeris: all bodies except Apophis
    # accelerations of "everybody else" are needed when evaluating Apophis post-Newtonian acceleration
    @show N
    Nm1 = N-1
    acc_eph = TaylorInterpolant(ss16ast_eph_t, Matrix{eltype(ss16ast_eph_x)}(undef, length(ss16ast_eph_t)-1, 3Nm1))
    newtonianNb_Potential = TaylorInterpolant(ss16ast_eph_t, Matrix{eltype(ss16ast_eph_x)}(undef, length(ss16ast_eph_t)-1, Nm1))
    fill!(acc_eph.x, zero(ss16ast_eph_x[1]))
    fill!(newtonianNb_Potential.x, zero(ss16ast_eph_x[1]))
    _1_to_Nm1 = Base.OneTo(Nm1) # iterator over all bodies except Apophis
    for j in _1_to_Nm1
        for i in _1_to_Nm1
            # i == j && continue
            if i == j
            else
                X_ij = ss16ast_eph_x[:,3i-2] .- ss16ast_eph_x[:,3j-2]
                Y_ij = ss16ast_eph_x[:,3i-1] .- ss16ast_eph_x[:,3j-1]
                Z_ij = ss16ast_eph_x[:,3i  ] .- ss16ast_eph_x[:,3j  ]
                r_p2_ij = ( (X_ij.^2) .+ (Y_ij.^2) ) .+ (Z_ij.^2)
                r_p3d2_ij = r_p2_ij.^1.5
                r_ij = sqrt.(r_p2_ij)
                newtonianCoeff_ij =  μ[i]./r_p3d2_ij
                acc_eph.x[:,3j-2] .+= (X_ij.*newtonianCoeff_ij)
                acc_eph.x[:,3j-1] .+= (Y_ij.*newtonianCoeff_ij)
                acc_eph.x[:,3j  ] .+= (Z_ij.*newtonianCoeff_ij)
                newtonianNb_Potential.x[:,j] .+= (μ[i]./r_ij)
            end #if i != j
        end #for, i
    end #for, j
    params = (ss16asteph, acc_eph, newtonianNb_Potential)

    # test **WITHOUT** jet transport
    q0 = initialcond()
    t0 = datetime2julian(DateTime(2008,9,24))
    tT = t0 + Taylor1(order)
    q0T = Taylor1.(q0, order)
    dq0T = similar(q0T)
    xaux = similar(q0T)
    q0T1 = Taylor1.(q0, order)
    dq0T1 = similar(q0T)

    @show methods(TaylorIntegration.jetcoeffs!)

    TaylorIntegration.jetcoeffs!(RNp1BP_pN_A_J23E_J2S_ng_eph!, tT, q0T, dq0T, xaux, params)
    @time TaylorIntegration.jetcoeffs!(RNp1BP_pN_A_J23E_J2S_ng_eph!, tT, q0T, dq0T, xaux, params)
    TaylorIntegration.jetcoeffs!(Val(RNp1BP_pN_A_J23E_J2S_ng_eph!), tT, q0T1, dq0T1, params)
    @time TaylorIntegration.jetcoeffs!(Val(RNp1BP_pN_A_J23E_J2S_ng_eph!), tT, q0T1, dq0T1, params)

    # @show q0T
    # @show q0T1
    @show norm(q0T-q0T1, Inf)
    @show norm(dq0T-dq0T1, Inf)
    @show q0T==q0T1

    # test **WITH** jet transport
    __q0 = Taylor1.(q0,10)
    __q0[1:end-1] = Taylor1.(q0[1:end-1],10)
    __q0[end] = Taylor1([q0[end],1e-14],10) #note the 1e-14!!!
    q0TT = Taylor1.(__q0, order)
    dq0TT = similar(q0TT)
    xauxTT = similar(q0TT)
    q0TT1 = Taylor1.(__q0, order)
    dq0TT1 = similar(q0TT1)
    q0TT2 = Taylor1.(__q0, order)
    dq0TT2 = similar(q0TT2)

    TaylorIntegration.jetcoeffs!(RNp1BP_pN_A_J23E_J2S_ng_eph!, tT, q0TT, dq0TT, xauxTT, params)
    @time TaylorIntegration.jetcoeffs!(RNp1BP_pN_A_J23E_J2S_ng_eph!, tT, q0TT, dq0TT, xauxTT, params)
    TaylorIntegration.jetcoeffs!(Val(RNp1BP_pN_A_J23E_J2S_ng_eph!), tT, q0TT1, dq0TT1, params)
    @time TaylorIntegration.jetcoeffs!(Val(RNp1BP_pN_A_J23E_J2S_ng_eph!), tT, q0TT1, dq0TT1, params)
    TaylorIntegration.jetcoeffs!(Val(RNp1BP_pN_A_J23E_J2S_ng_eph_threads!), tT, q0TT2, dq0TT2, params)
    @time TaylorIntegration.jetcoeffs!(Val(RNp1BP_pN_A_J23E_J2S_ng_eph_threads!), tT, q0TT2, dq0TT2, params)

    # @show q0TT
    # @show q0TT1
    @show norm(q0TT-q0TT1, Inf)
    @show norm(dq0TT-dq0TT1, Inf)
    @show q0TT==q0TT1
    @show Threads.nthreads()
    @show norm(q0TT2-q0TT1, Inf)
    @show norm(dq0TT2-dq0TT1, Inf)
    @show q0TT2==q0TT1
end
