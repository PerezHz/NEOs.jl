function apophisstep!(f!, t::Taylor1{T}, x::Vector{Taylor1{U}},
        dx::Vector{Taylor1{U}}, xaux::Vector{Taylor1{U}}, abstol::T, params,
        parse_eqs::Bool=true) where {T<:Real, U<:Number}
    # Compute the Taylor coefficients
    TaylorIntegration.__jetcoeffs!(Val(parse_eqs), f!, t, x, dx, xaux, params)
    # Compute the step-size of the integration using `abstol`
    δt = TaylorIntegration.stepsize(x, abstol)
    # Force Apophis time-step to be no larger than planetary ephemeris time-step
    ind, _ = PlanetaryEphemeris.getinterpindex(params[1], t[0])
    Δt = abs((params[1].t[ind+1]+params[1].t0) - t[0])
    δt = min(δt, Δt)
    return δt
end

function lyap_apophisstep!(f!, t::Taylor1{T}, x::Vector{Taylor1{U}},
        dx::Vector{Taylor1{U}}, xaux::Vector{Taylor1{U}},
        δx::Array{TaylorN{Taylor1{U}},1}, dδx::Array{TaylorN{Taylor1{U}},1},
        jac::Array{Taylor1{U},2}, t0::T, t1::T, order::Int,
        abstol::T, _δv::Vector{TaylorN{Taylor1{U}}}, varsaux::Array{Taylor1{U},3},
        params, parse_eqs::Bool=true, jacobianfunc! =nothing) where {T<:Real, U<:Number}

    # Dimensions of phase-space: dof
    nx = length(x)
    dof = length(δx)

    # Compute the Taylor coefficients associated to trajectory
    TaylorIntegration.__jetcoeffs!(Val(parse_eqs), f!, t, view(x, 1:dof), view(dx, 1:dof), view(xaux, 1:dof), params)

    # Compute stability matrix
    TaylorIntegration.stabilitymatrix!(f!, t, x, δx, dδx, jac, _δv, params, jacobianfunc!)

    # Compute the Taylor coefficients associated to variational equations
    TaylorIntegration.lyap_jetcoeffs!(t, view(x, dof+1:nx), view(dx, dof+1:nx), jac, varsaux)

    # Compute the step-size of the integration using `abstol`
    δt = TaylorIntegration.stepsize(view(x, 1:dof), abstol)
    # Force Apophis time-step to be no larger than planetary ephemeris time-step
    ind, _ = PlanetaryEphemeris.getinterpindex(params[1], t[0])
    Δt = abs((params[1].t[ind+1]+params[1].t0) - t[0])
    δt = min(δt, Δt)

    δt = min(δt, t1-t0)

    return δt
end

function apophisinteg(f!, q0::Array{U,1}, t0::T, tmax::T, order::Int, abstol::T,
        params = nothing; maxsteps::Int=500, parse_eqs::Bool=true, dense::Bool=false) where {T<:Real, U<:Number}

    # Allocation
    tv = Array{T}(undef, maxsteps+1)
    dof = length(q0)
    xv = Array{U}(undef, dof, maxsteps+1)
    if dense
        xv_interp = Array{Taylor1{U}}(undef, dof, maxsteps+1)
    end

    # Initialize the vector of Taylor1 expansions
    t = Taylor1(T, order)
    x = Array{Taylor1{U}}(undef, dof)
    dx = Array{Taylor1{U}}(undef, dof)
    xaux = Array{Taylor1{U}}(undef, dof)
    dx .= Taylor1.(zeros(U), order)

    # Initial conditions
    @inbounds t[0] = t0
    x .= Taylor1.(q0, order)
    x0 = deepcopy(q0)
    @inbounds tv[1] = t0
    @inbounds xv[:,1] .= q0
    sign_tstep = copysign(1, tmax-t0)

    # Determine if specialized jetcoeffs! method exists
    parse_eqs = TaylorIntegration._determine_parsing!(parse_eqs, f!, t, x, dx, params)

    # Integration
    nsteps = 1
    while sign_tstep*t0 < sign_tstep*tmax
        δt = apophisstep!(f!, t, x, dx, xaux, abstol, params, parse_eqs) # δt is positive!
        # Below, δt has the proper sign according to the direction of the integration
        δt = sign_tstep * min(δt, sign_tstep*(tmax-t0))
        evaluate!(x, δt, x0) # new initial condition
        if dense
            xv_interp[:,nsteps] .= deepcopy(x)
        end
        for i in eachindex(x0)
            @inbounds x[i][0] = x0[i]
            @inbounds dx[i] = Taylor1( zero(x0[i]), order )
        end
        t0 += δt
        @inbounds t[0] = t0
        nsteps += 1
        @inbounds tv[nsteps] = t0
        @inbounds xv[:,nsteps] .= x0
        if nsteps > maxsteps
            @warn("""
            Maximum number of integration steps reached; exiting.
            """)
            break
        end
    end

    if dense
        return TaylorInterpolant(tv[1], view(tv.-tv[1],1:nsteps), view(transpose(view(xv_interp,:,1:nsteps-1)),1:nsteps-1,:))
    else
        return view(tv,1:nsteps), view(transpose(view(xv,:,1:nsteps)),1:nsteps,:)
    end
end

#TODO: allow backwards integration
# root-finding integration method
function apophisinteg(f!, g, q0::Array{U,1}, t0::T, tmax::T, order::Int,
        abstol::T, params = nothing; maxsteps::Int=500, parse_eqs::Bool=true,
        dense::Bool=false, eventorder::Int=0, newtoniter::Int=10,
        nrabstol::T=eps(T)) where {T <: Real,U <: Number}

    @assert order ≥ eventorder "`eventorder` must be less than or equal to `order`"

    # Allocation
    tv = Array{T}(undef, maxsteps+1)
    dof = length(q0)
    xv = Array{U}(undef, dof, maxsteps+1)
    if dense
        xv_interp = Array{Taylor1{U}}(undef, dof, maxsteps+1)
    end

    # Initialize the vector of Taylor1 expansions
    t = Taylor1(T, order)
    x = Array{Taylor1{U}}(undef, dof)
    dx = Array{Taylor1{U}}(undef, dof)
    xaux = Array{Taylor1{U}}(undef, dof)
    dx = Array{Taylor1{U}}(undef, dof)

    # Initial conditions
    @inbounds t[0] = t0
    x .= Taylor1.(q0, order)
    dx .= zero.(x)
    x0 = deepcopy(q0)
    @inbounds tv[1] = t0
    @inbounds xv[:,1] .= q0
    sign_tstep = copysign(1, tmax-t0)

    # Determine if specialized jetcoeffs! method exists
    parse_eqs = TaylorIntegration._determine_parsing!(parse_eqs, f!, t, x, dx, params)

    # Some auxiliary arrays for root-finding/event detection/Poincaré surface of section evaluation
    g_tupl = g(dx, x, params, t)
    g_tupl_old = g(dx, x, params, t)
    # g_val = zero(gg_tupl[2])
    # g_val_old = zero(g_val)
    slope = zero(x[1])
    dt_li = zero(x[1])
    dt_nr = zero(x[1])
    δt = zero(x[1])
    δt_old = zero(x[1])

    x_dx = vcat(x, dx)
    # g_dg = vcat(g_val, g_val_old)
    g_dg = vcat(g_tupl[2], g_tupl_old[2])
    x_dx_val = Array{U}(undef, length(x_dx) )
    # g_dg_val = vcat(evaluate(g_val), evaluate(g_val_old))
    g_dg_val = vcat(evaluate(g_tupl[2]), evaluate(g_tupl_old[2]))

    tvS = Array{U}(undef, maxsteps+1)
    xvS = similar(xv)
    gvS = similar(tvS)

    # Integration
    nsteps = 1
    nevents = 1 #number of detected events
    while sign_tstep*t0 < sign_tstep*tmax
        δt_old = δt
        δt = apophisstep!(f!, t, x, dx, xaux, abstol, params, parse_eqs) # δt is positive!
        # Below, δt has the proper sign according to the direction of the integration
        δt = sign_tstep * min(δt, sign_tstep*(tmax-t0))
        evaluate!(x, δt, x0) # new initial condition
        # g_val = g(dx, x, params, t)
        # nevents = findroot!(t, x, dx, g_val_old, g_val, eventorder,
        #     tvS, xvS, gvS, t0, δt_old, x_dx, x_dx_val, g_dg, g_dg_val,
        #     nrabstol, newtoniter, nevents)
        # g_val_old = deepcopy(g_val)
        g_tupl = g(dx, x, params, t)
        nevents = TaylorIntegration.findroot!(t, x, dx, g_tupl_old, g_tupl, eventorder,
        tvS, xvS, gvS, t0, δt_old, x_dx, x_dx_val, g_dg, g_dg_val,
        nrabstol, newtoniter, nevents)
        if dense
            xv_interp[:,nsteps] .= deepcopy(x)
        end
        g_tupl_old = deepcopy(g_tupl)
        for i in eachindex(x0)
            @inbounds x[i][0] = x0[i]
        end
        t0 += δt
        @inbounds t[0] = t0
        nsteps += 1
        @inbounds tv[nsteps] = t0
        @inbounds xv[:,nsteps] .= x0
        if nsteps > maxsteps
            @warn("""
            Maximum number of integration steps reached; exiting.
            """)
            break
        end
    end
    if dense
        return TaylorInterpolant(tv[1], view(tv.-tv[1],1:nsteps), view(transpose(view(xv_interp,:,1:nsteps-1)),1:nsteps-1,:)), view(tvS,1:nevents-1), view(transpose(view(xvS,:,1:nevents-1)),1:nevents-1,:), view(gvS,1:nevents-1)
    else
        return view(tv,1:nsteps), view(transpose(view(xv,:,1:nsteps)),1:nsteps,:), view(tvS,1:nevents-1), view(transpose(view(xvS,:,1:nevents-1)),1:nevents-1,:), view(gvS,1:nevents-1)
    end
end

function lyap_apophisinteg(f!, q0::Array{U,1}, t0::T, tmax::T,
        order::Int, abstol::T, params = nothing, jacobianfunc! =nothing;
        maxsteps::Int=500, parse_eqs::Bool=true) where {T<:Real, U<:Number}
    # Allocation
    tv = Array{T}(undef, maxsteps+1)
    dof = length(q0)
    xv = Array{U}(undef, dof, maxsteps+1)
    λ = similar(xv)
    λtsum = similar(q0)
    jt = Matrix{U}(I, dof, dof)
    _δv = Array{TaylorN{Taylor1{U}}}(undef, dof)

    # Initial conditions
    @inbounds tv[1] = t0
    @inbounds for ind in eachindex(q0)
        xv[ind,1] = q0[ind]
        λ[ind,1] = zero(U)
        λtsum[ind] = zero(U)
    end
    x0 = vcat(q0, reshape(jt, dof*dof))
    nx0 = length(x0)
    t00 = t0

    # Initialize the vector of Taylor1 expansions
    t = Taylor1(T, order)
    x = Array{Taylor1{U}}(undef, nx0)
    x .= Taylor1.( x0, order )
    @inbounds t[0] = t0

    # If user does not provide Jacobian, check number of TaylorN variables and initialize _δv
    if isa(jacobianfunc!, Nothing)
        @assert get_numvars() == dof "`length(q0)` must be equal to number of variables set by `TaylorN`"
        for ind in eachindex(q0)
            _δv[ind] = one(x[1])*TaylorN(Taylor1{U}, ind, order=1)
        end
    end

    #Allocate auxiliary arrays
    dx = Array{Taylor1{U}}(undef, nx0)
    xaux = Array{Taylor1{U}}(undef, nx0)
    δx = Array{TaylorN{Taylor1{U}}}(undef, dof)
    dδx = Array{TaylorN{Taylor1{U}}}(undef, dof)
    jac = Array{Taylor1{U}}(undef, dof, dof)
    varsaux = Array{Taylor1{U}}(undef, dof, dof, dof)
    fill!(jac, zero(x[1]))
    QH = Array{U}(undef, dof, dof)
    RH = Array{U}(undef, dof, dof)
    aⱼ = Array{U}(undef, dof )
    qᵢ = similar(aⱼ)
    vⱼ = similar(aⱼ)

    # Determine if specialized jetcoeffs! method exists
    parse_eqs = TaylorIntegration._determine_parsing!(parse_eqs, f!, t, view(x, 1:dof), view(dx, 1:dof), params)

    # Integration
    nsteps = 1
    while t0 < tmax
        δt = lyap_apophisstep!(f!, t, x, dx, xaux, δx, dδx, jac, t0, tmax,
            order, abstol, _δv, varsaux, params, parse_eqs, jacobianfunc!)
        evaluate!(x, δt, x0) # Update x0
        for ind in eachindex(jt)
            @inbounds jt[ind] = x0[dof+ind]
        end
        TaylorIntegration.modifiedGS!( jt, QH, RH, aⱼ, qᵢ, vⱼ )
        t0 += δt
        @inbounds t[0] = t0
        tspan = t0-t00
        nsteps += 1
        @inbounds tv[nsteps] = t0
        @inbounds for ind in eachindex(q0)
            xv[ind,nsteps] = x0[ind]
            λtsum[ind] += log(RH[ind,ind])
            λ[ind,nsteps] = λtsum[ind]/tspan
        end
        for ind in eachindex(QH)
            @inbounds x0[dof+ind] = QH[ind]
        end
        x .= Taylor1.( x0, order )
        if nsteps > maxsteps
            @warn("""
            Maximum number of integration steps reached; exiting.
            """)
            break
        end
    end

    return view(tv,1:nsteps),  view(transpose(xv),1:nsteps,:),  view(transpose(λ),1:nsteps,:)
end
