function apophisstep!(f!, t::Taylor1{T}, x::Vector{Taylor1{U}}, dx::Vector{Taylor1{U}},
        xaux::Vector{Taylor1{U}}, t0::T, t1::T, order::Int, abstol::T, params,
        parse_eqs::Bool=true) where {T<:Real, U<:Number}

    @assert t1 > t0

    # Compute the Taylor coefficients
    TaylorIntegration.__jetcoeffs!(Val(parse_eqs), f!, t, x, dx, xaux, params)

    # Compute the step-size of the integration using `abstol`
    δt = TaylorIntegration.stepsize(x, abstol)
    # Handle case: when Apophis time-step is larger than ephemeris time-step
    ind = findlast(x->x≤t0, params[1].t)
    Δt = params[1].t[ind+1] - t0
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

    # Determine if specialized jetcoeffs! method exists
    parse_eqs = parse_eqs && (length(methods(TaylorIntegration.jetcoeffs!)) > 2)
    if parse_eqs
        try
            TaylorIntegration.jetcoeffs!(Val(f!), t, x, dx, params)
        catch
            parse_eqs = false
        end
    end

    # Integration
    nsteps = 1
    while t0 < tmax
        δt = apophisstep!(f!, t, x, dx, xaux, t0, tmax, order, abstol, params, parse_eqs)
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
            @info("""
            Maximum number of integration steps reached; exiting.
            """)
            break
        end
    end

    if dense
        return TaylorInterpolant(view(tv,1:nsteps), view(transpose(view(xv_interp,:,1:nsteps-1)),1:nsteps-1,:))
    else
        return view(tv,1:nsteps), view(transpose(view(xv,:,1:nsteps)),1:nsteps,:)
    end
end

# root-finding integration method
function apophisinteg(f!, g, q0::Array{U,1}, trange::AbstractVector{T},
        order::Int, abstol::T, params = nothing; maxsteps::Int=500, parse_eqs::Bool=true,
        eventorder::Int=0, newtoniter::Int=10, nrabstol::T=eps(T)) where {T <: Real,U <: Number}

    @assert order ≥ eventorder "`eventorder` must be less than or equal to `order`"

    # Allocation
    nn = length(trange)
    dof = length(q0)
    x0 = similar(q0, eltype(q0), dof)
    fill!(x0, T(NaN))
    xv = Array{eltype(q0)}(undef, dof, nn)
    for ind in 1:nn
        @inbounds xv[:,ind] .= x0
    end

    # Initialize the vector of Taylor1 expansions
    t = Taylor1( T, order )
    x = Array{Taylor1{U}}(undef, dof)
    dx = Array{Taylor1{U}}(undef, dof)
    xaux = Array{Taylor1{U}}(undef, dof)
    for i in eachindex(q0)
        @inbounds x[i] = Taylor1( q0[i], order )
        @inbounds dx[i] = Taylor1( zero(q0[i]), order )
    end

    # Initial conditions
    @inbounds t[0] = trange[1]
    @inbounds t0, t1, tmax = trange[1], trange[2], trange[end]
    x0 = deepcopy(q0)
    x1 = similar(x0)
    x .= Taylor1.(q0, order)
    @inbounds xv[:,1] .= q0

    # Determine if specialized jetcoeffs! method exists
    parse_eqs = parse_eqs && (length(methods(TaylorIntegration.jetcoeffs!)) > 2)
    if parse_eqs
        try
            TaylorIntegration.jetcoeffs!(Val(f!), t, x, dx, params)
        catch
            parse_eqs = false
        end
    end

    #Some auxiliary arrays for root-finding/event detection/Poincaré surface of section evaluation
    g_val = zero(g(x, x, params, t))
    g_val_old = zero(g_val)
    slope = zero(U)
    dt_li = zero(U)
    dt_nr = zero(U)
    δt = zero(U)
    δt_old = zero(U)

    x_dx = vcat(x, dx)
    g_dg = vcat(g_val, g_val_old)
    x_dx_val = Array{U}(undef, length(x_dx) )
    g_dg_val = vcat(evaluate(g_val), evaluate(g_val_old))

    tvS = Array{U}(undef, maxsteps+1)
    xvS = similar(xv)
    gvS = similar(tvS)

    # Integration
    iter = 2
    nsteps = 1
    nevents = 1 #number of detected events
    while t0 < tmax
        δt_old = δt
        δt = apophisstep!(f!, t, x, dx, xaux, t0, tmax, order, abstol, params, parse_eqs)
        evaluate!(x, δt, x0) # new initial condition
        tnext = t0+δt
        # Evaluate solution at times within convergence radius
        while t1 < tnext
            evaluate!(x, t1-t0, x1)
            @inbounds xv[:,iter] .= x1
            iter += 1
            @inbounds t1 = trange[iter]
        end
        if δt == tmax-t0
            @inbounds xv[:,iter] .= x0
            break
        end
        g_val = g(dx, x, params, t)
        nevents = TaylorIntegration.findroot!(t, x, dx, g_val_old, g_val, eventorder,
            tvS, xvS, gvS, t0, δt_old, x_dx, x_dx_val, g_dg, g_dg_val,
            nrabstol, newtoniter, nevents)
        g_val_old = deepcopy(g_val)
        for i in eachindex(x0)
            @inbounds x[i][0] = x0[i]
        end
        t0 = tnext
        @inbounds t[0] = t0
        nsteps += 1
        if nsteps > maxsteps
            @info("""
            Maximum number of integration steps reached; exiting.
            """)
            break
        end
    end

    return transpose(xv), view(tvS,1:nevents-1), view(transpose(view(xvS,:,1:nevents-1)),1:nevents-1,:), view(gvS,1:nevents-1)
end
