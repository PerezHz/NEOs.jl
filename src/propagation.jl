#numerator of Apophis radial velocity wrt Earth
function rvelea(dx, x, params, t)
    ss16asteph_t = params[1](t) #ss16asteph(t)
    xe = ss16asteph_t[union(3ea-2:3ea,3(N-1+ea)-2:3(N-1+ea))]
    return (x[1]-xe[1])*(x[4]-xe[4]) + (x[2]-xe[2])*(x[5]-xe[5]) + (x[3]-xe[3])*(x[6]-xe[6])
end

function loadeph(ephfile)
    # read Solar System ephemeris (Sun+8 planets+Moon+Pluto+16 main belt asteroids)
    ss16asteph_ = load(ephfile, "ss16ast_eph")
    ss16asteph_t0 = (ss16asteph_.t0 ./ daysec) - (jd0-J2000)
    ss16asteph_t = (ss16asteph_.t ./ daysec)
    ephord = ss16asteph_.x[1].order
    ss16asteph_x = map(x->x(Taylor1(ephord)*daysec), ss16asteph_.x)
    ss16asteph = TaylorInterpolant(ss16asteph_t0, ss16asteph_t, ss16asteph_x)
    #compute point-mass Newtonian accelerations from ephemeris: all bodies except Apophis
    # accelerations of "everybody else" are needed when evaluating Apophis post-Newtonian acceleration
    Nm1 = N-1
    acc_eph = TaylorInterpolant(ss16asteph.t0, ss16asteph.t, Matrix{eltype(ss16asteph.x)}(undef, length(ss16asteph.t)-1, 3Nm1))
    newtonianNb_Potential = TaylorInterpolant(ss16asteph.t0, ss16asteph.t, Matrix{eltype(ss16asteph.x)}(undef, length(ss16asteph.t)-1, Nm1))
    fill!(acc_eph.x, zero(ss16asteph.x[1]))
    fill!(newtonianNb_Potential.x, zero(ss16asteph.x[1]))
    _1_to_Nm1 = Base.OneTo(Nm1) # iterator over all bodies except Apophis
    for j in _1_to_Nm1
        for i in _1_to_Nm1
            # i == j && continue
            if i == j
            else
                X_ij = ss16asteph.x[:,3i-2] .- ss16asteph.x[:,3j-2]
                Y_ij = ss16asteph.x[:,3i-1] .- ss16asteph.x[:,3j-1]
                Z_ij = ss16asteph.x[:,3i  ] .- ss16asteph.x[:,3j  ]
                r_p2_ij = ( (X_ij.^2) .+ (Y_ij.^2) ) .+ (Z_ij.^2)
                r_ij = sqrt.(r_p2_ij)
                newtonianNb_Potential.x[:,j] .+= (μ[i]./r_ij)
            end #if i != j
        end #for, i
        acc_eph.x[:,3j-2] .= differentiate.(ss16asteph.x[:,3(Nm1+j)-2])
        acc_eph.x[:,3j-1] .= differentiate.(ss16asteph.x[:,3(Nm1+j)-1])
        acc_eph.x[:,3j  ] .= differentiate.(ss16asteph.x[:,3(Nm1+j)  ])
    end #for, j
    return ss16asteph, acc_eph, newtonianNb_Potential
end

function save2jldandcheck(objname, sol)
    outfilename = string(objname, "_jt.", myid()-1, ".jld")
    return __save2jldandcheck(outfilename, sol)
end

function __save2jldandcheck(outfilename, sol)
    println("Saving solution to file: $outfilename")
    jldopen(outfilename, "w") do file
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
        recovered_sol_i = load(outfilename, varname)
        #check that varname was recovered succesfully
        @show recovered_sol_i == sol[ind]
    end
    println("Saved solution")
    return outfilename
end

function taylor_minimum(pol::Taylor1{T}, x0::T; niters::Int=10) where {T<:Real}
    dpol = differentiate(pol)
    dpol2 = differentiate(dpol)
    xnewton::T = x0
    #@show xnewton
    for i in 1:niters
        xnewton -= dpol(xnewton)/dpol2(xnewton)
        #@show xnewton, dpol(xnewton)
    end
    return xnewton
end

function taylor_roots(pol::Taylor1{T}, x0::T; niters::Int=10) where {T<:Real}
    dpol = differentiate(pol)
    xnewton::T = x0
    #@show xnewton
    for i in 1:niters
        xnewton -= pol(xnewton)/dpol(xnewton)
        #@show xnewton, pol(xnewton)
    end
    return xnewton
end

function least_squares_A2(asteroid_data::Vector{RadarDataJPL{T}},
        vdel::Vector{Taylor1{U}}, vdop::Vector{Taylor1{U}}) where {T<:Number, U<:Number}
    delay_index = findall(x->x.delay_units=="us", asteroid_data)
    doppler_index = findall(x->x.doppler_units=="Hz", asteroid_data)
    tdelay_jpl_obs = [x.delay for x in asteroid_data][delay_index]
    dshift_jpl_obs = [x.doppler for x in asteroid_data][doppler_index]
    tdelay_jpl_obs_sigma = [x.delay_sigma for x in asteroid_data][delay_index]
    dshift_jpl_obs_sigma = [x.doppler_sigma for x in asteroid_data][doppler_index]
    res_del = tdelay_jpl_obs .- vdel # delay residuals a.a.f. of A2
    res_dop = dshift_jpl_obs .- vdop # Doppler residuals a.a.f. of A2
    W_del = 1 ./(tdelay_jpl_obs_sigma.^2) # delay weights
    W_dop = 1 ./(dshift_jpl_obs_sigma.^2) # Doppler weights
    res_deldop = vcat(res_del, res_dop) # delay + Doppler residuals as a function of A2
    W_deldop = vcat(W_del, W_dop) # delay + Doppler weights
    res_xWx_deldop = res_deldop .* W_deldop .* res_deldop
    Q_A2_deldop = sum(res_xWx_deldop)/length(res_xWx_deldop)
    A2_lsqfit_deldop = taylor_minimum(Q_A2_deldop, 0.0, niters=5)
    B_deldop = differentiate.(res_deldop) # design matrix
    C_deldop = ( transpose(B_deldop) )*( W_deldop.*B_deldop ) # normal matrix
    Γ_deldop = inv(C_deldop) # covariance matrix
    return A2_lsqfit_deldop, Γ_deldop(A2_lsqfit_deldop)
end

function propagate(objname::String, dynamics::Function, maxsteps::Int, t0::T,
        tspan::T, ephfile::String; output::Bool=true, newtoniter::Int=10,
        dense::Bool=false, dq::Vector=zeros(7), radarobsfile::String="") where {T<:Real}

    ss16asteph, acc_eph, newtonianNb_Potential = loadeph(ephfile)
    jd0 = datetime2julian(DateTime(2008, 9, 24))
    params = (ss16asteph, acc_eph, newtonianNb_Potential, jd0)
    # get asteroid initial conditions
    q0 = initialcond(dq)
    @show q0

    @show tmax = t0+tspan*yr #final time of integration

    # propagate orbit
    if dense
        @time interp = apophisinteg(dynamics, q0, t0, tmax, order, abstol, params; maxsteps=maxsteps, dense=dense)
        et0 = (jd0-J2000)*daysec
        etv = interp.t[:]*daysec
        interp_x_et = map(x->x(Taylor1(order)/daysec), interp.x[:,:])
        apophis = TaylorInterpolant(et0, etv, interp_x_et)
        sol = (apophis=apophis,
        )
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

    #write solution and predicted values of observations (if requested) to .jld files
    if output
        outfilename = save2jldandcheck(objname, sol)
        if dense
            # if requested by user, calculate computed (i.e., predicted) values of observations
            compute_radar_obs(outfilename, radarobsfile, interp, ss16asteph)
        end
    end

    return nothing
end

function compute_radar_obs(outfilename::String, radarobsfile::String, apophis_interp, ss16asteph; tc::Real=3.0)
    if radarobsfile != ""
        asteroid_data = process_radar_data_jpl(radarobsfile)
        # TODO: check that first and last observation times are within interpolation interval
        jd0 = datetime2julian(DateTime(2008,9,24))
        function apophis_et(et)
            return auday2kmsec(apophis_interp(et)[1:6])
        end
        function earth_et(et)
            return auday2kmsec(ss16asteph(et)[union(3*4-2:3*4,3*(N-1+4)-2:3*(N-1+4))])
        end
        function sun_et(et)
            return auday2kmsec(ss16asteph(et)[union(3*1-2:3*1,3*(N-1+1)-2:3*(N-1+1))])
        end
        #compute time-delay and Doppler-shift "ephemeris" (i.e., predicted values according to ephemeris)
        vdel, vdop = delay_doppler(asteroid_data, tc=tc, xve=earth_et, xvs=sun_et, xva=apophis_et)
        sol = (vdel=vdel, vdop=vdop)
        #save data to file
        __save2jldandcheck(outfilename, sol)
    end
    return nothing
end

# distributed computing (Monte-Carlo) version of `propagate`
function propagate_distributed(objname::String, dynamics::Function, maxsteps::Int,
        t0::T, tmax::T, aux; output::Bool=true, newtoniter::Int=10,
        dq::Vector=zeros(7), radarobsfile::String="") where {T<:Real}

    ss16asteph, acc_eph, newtonianNb_Potential, earth_et, sun_et = aux
    params = aux[1:3]

    # get asteroid initial conditions
    q0 = initialcond(dq)

    @show myid()

    # do integration
    if output && radarobsfile != ""
        asteroid_data = process_radar_data_jpl(radarobsfile)
        # TODO: check that first and last observation times are within interpolation interval
        @time interp = apophisinteg(dynamics, q0, t0, tmax, order, abstol, params; maxsteps=maxsteps, dense=true)
        function apophis_et(et)
            return interp( etsecs2julian(et) )[1:6]
        end
        #compute time-delay and Doppler-shift "ephemeris" (i.e., predicted values according to ephemeris)
        vdel, vdop = delay_doppler(asteroid_data; xve=earth_et, xvs=sun_et, xva=apophis_et)
        A2, Γ_A2 = least_squares_A2(asteroid_data, vdel, vdop)
        sol = (t=interp.t[:], x=interp.x[:,:], vdel=vdel, vdop=vdop, A2=A2, Γ_A2=Γ_A2)
    else
        @time interp = apophisinteg(dynamics, q0, t0, tmax, order, abstol, params; maxsteps=maxsteps, dense=true)
        sol = (t=interp.t[:], x=interp.x[:,:])
    end

    #write solution to .jld files
    if output
        save2jldandcheck(objname, sol)
    end
    return nothing
end

function parallel_run(objname::String, dynamics::Function, maxsteps::Int,
        t0::T, tmax::T, aux; output::Bool=true, newtoniter::Int=10,
        radarobsfile::String="") where {T<:Real}

    varorder = 5 # varorder is the order corresponding to the jet transport perturbation
    dxv = [1e-8randn(6) for w in workers()]

    dqv = Vector{typeof(Taylor1.(dxv[1], varorder))}(undef, length(dxv))
    for j in eachindex(dqv)
        # dqv[j]: perturbation to nominal initial condition (Taylor1 jet transport)
        dqv[j] = Taylor1.(zeros(7), varorder)
        for i in 1:6
            dqv[j][i][0] = dxv[j][i]
        end
        dqv[j][end][1] = 1e-14
        # dq: perturbation to nominal initial condition (TaylorN jet transport)
        # dq = set_variables("ξ", order=varorder, numvars=7)
        # for i in 1:6
        #     dq[i][1][i] = 1e-8
        # end
        # dq[end][1][end] = 1e-14
    end
    @show dqv

    f1 = x -> propagate_distributed(objname, dynamics, maxsteps, t0, tmax,
        aux, output=output, radarobsfile=radarobsfile, dq=x)
    pmap(f1, dqv[1:nworkers()])
end

function testjetcoeffs(ephfile)
    # read Solar System ephemeris (Sun+8 planets+Moon+Pluto+16 main belt asteroids)
    ss16asteph, acc_eph, newtonianNb_Potential = loadeph(ephfile)
    testjetcoeffs(ss16asteph, acc_eph, newtonianNb_Potential)
end

function testjetcoeffs(ss16asteph, acc_eph, newtonianNb_Potential)
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

    # Determine if specialized jetcoeffs! method runs without trouble (parsed but no threads)
    parse_eqs = true
    if parse_eqs
        try
            TaylorIntegration.jetcoeffs!(Val(RNp1BP_pN_A_J23E_J2S_ng_eph!), tT, q0T1, dq0T1, params)
        catch
            parse_eqs = false
        end
    end
    @show parse_eqs # should evaluate to true, if evaluation of parsed jetcoeffs! method worked
    q0T1 = Taylor1.(q0, order)
    dq0T1 = similar(q0T)
    # Determine if specialized jetcoeffs! method runs without trouble (parsed+threads)
    parse_eqs = true
    if parse_eqs
        try
            TaylorIntegration.jetcoeffs!(Val(RNp1BP_pN_A_J23E_J2S_ng_eph_threads!), tT, q0T1, dq0T1, params)
        catch
            parse_eqs = false
        end
    end
    @show parse_eqs # should evaluate to true, if evaluation of parsed jetcoeffs! method worked

    # q0T1 = Taylor1.(q0, order)
    # dq0T1 = similar(q0T)

    @show methods(TaylorIntegration.jetcoeffs!)

    TaylorIntegration.jetcoeffs!(RNp1BP_pN_A_J23E_J2S_ng_eph!, tT, q0T, dq0T, xaux, params)
    @time TaylorIntegration.jetcoeffs!(RNp1BP_pN_A_J23E_J2S_ng_eph!, tT, q0T, dq0T, xaux, params)
    TaylorIntegration.jetcoeffs!(Val(RNp1BP_pN_A_J23E_J2S_ng_eph!), tT, q0T1, dq0T1, params)
    @time TaylorIntegration.jetcoeffs!(Val(RNp1BP_pN_A_J23E_J2S_ng_eph!), tT, q0T1, dq0T1, params)

    # @show q0T
    # @show q0T1
    # Test equality of non-parsed vs parsed jetcoeffs! methods
    @show norm(q0T-q0T1, Inf)
    @show norm(dq0T-dq0T1, Inf)
    @show q0T==q0T1

    # test **WITH** jet transport
    # dq: perturbation to nominal initial condition (Taylor1 jet transport)
    dq = Taylor1.(zeros(7), 10)
    dq[end][1] = 1e-14 #note the 1e-14!!!
    __q0 = initialcond(dq)
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
