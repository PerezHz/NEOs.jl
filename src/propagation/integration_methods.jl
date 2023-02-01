@doc raw"""
    neos_jetcoeffs!(eqsdiff!::Function, t::Taylor1{T}, x::AbstractArray{Taylor1{U}, N}, dx::AbstractArray{Taylor1{U}, N},
                    xaux::AbstractArray{Taylor1{U}, N}, params) where {T <: Real, U <: Number, N}

Specialized method of `TaylorIntegration.jetcoeffs!` for the integration of a NEO. 

See also [`TaylorIntegration.jetcoeffs!`](@ref). 
"""
function neos_jetcoeffs!(eqsdiff!::Function, t::Taylor1{T}, x::AbstractArray{Taylor1{U}, N}, dx::AbstractArray{Taylor1{U}, N},
                    xaux::AbstractArray{Taylor1{U}, N}, params) where {T <: Real, U <: Number, N}
    order = x[1].order
    for ord in 0:order-1
        ordnext = ord+1

        # Set `taux`, auxiliary Taylor1 variable to order `ord`
        @inbounds taux = Taylor1( t.coeffs )
        # Set `xaux`, auxiliary vector of Taylor1 to order `ord`
        for j in eachindex(x)
            @inbounds xaux[j] = Taylor1( x[j].coeffs[1:ordnext] )
        end
        
        # Equations of motion
        eqsdiff!(dx, xaux, params, taux)

        # Recursion relations
        for j in eachindex(x)
            @inbounds x[j].coeffs[ordnext+1] = dx[j].coeffs[ordnext]/ordnext
        end
    end
    nothing
end

@doc raw"""
    neosstep!(f!, t::Taylor1{T}, x::Vector{Taylor1{U}}, dx::Vector{Taylor1{U}}, 
              xaux::Vector{Taylor1{U}}, abstol::T, params) where {T <: Real, U <: Number}
    neosstep!(f!, t::Taylor1{T}, x::Vector{Taylor1{U}}, dx::Vector{Taylor1{U}}, 
              abstol::T, params, rv::RetAlloc{Taylor1{U}}) where {T <: Real, U <: Number}

Specialized methods of `TaylorIntegration.taylorstep!` for the integration of a NEO. 

See also [`TaylorIntegration.taylorstep!`](@ref).
"""
function neosstep!(f!, t::Taylor1{T}, x::Vector{Taylor1{U}}, dx::Vector{Taylor1{U}}, 
                   xaux::Vector{Taylor1{U}}, abstol::T, params) where {T <: Real, U <: Number}

    # Compute the Taylor coefficients
    neos_jetcoeffs!(f!, t, x, dx, xaux, params)

    # Compute the step-size of the integration using `abstol`
    δt = TaylorIntegration.stepsize(x, abstol)

    return δt
end

function neosstep!(f!, t::Taylor1{T}, x::Vector{Taylor1{U}}, dx::Vector{Taylor1{U}}, 
                   abstol::T, params, rv::TaylorIntegration.RetAlloc{Taylor1{U}}) where {T <: Real, U <: Number}

    # Compute the Taylor coefficients
    TaylorIntegration.__jetcoeffs!(Val(true), f!, t, x, dx, params, rv)

    # Compute the step-size of the integration using `abstol`
    δt = TaylorIntegration.stepsize(x, abstol)

    return δt
end

@doc raw"""
    lyap_neosstep!(f!, t::Taylor1{T}, x::Vector{Taylor1{U}}, dx::Vector{Taylor1{U}}, xaux::Vector{Taylor1{U}},
                   δx::Array{TaylorN{Taylor1{U}}, 1}, dδx::Array{TaylorN{Taylor1{U}}, 1}, jac::Array{Taylor1{U}, 2}, 
                   abstol::T, _δv::Vector{TaylorN{Taylor1{U}}}, varsaux::Array{Taylor1{U},3}, params, 
                   jacobianfunc! = nothing) where {T <: Real, U <: Number}
    lyap_neosstep!(f!, t::Taylor1{T}, x::Vector{Taylor1{U}}, dx::Vector{Taylor1{U}}, δx::Array{TaylorN{Taylor1{U}}, 1}, 
                   dδx::Array{TaylorN{Taylor1{U}}, 1}, jac::Array{Taylor1{U}, 2}, abstol::T, _δv::Vector{TaylorN{Taylor1{U}}},
                   varsaux::Array{Taylor1{U}, 3}, params, rv::TaylorIntegration.RetAlloc{Taylor1{U}}, 
                   jacobianfunc! = nothing) where {T <: Real, U <: Number}

Specialized method of `TaylorIntegration.lyap_taylorstep!` for the calculation of the Lyapunov spectrum of a NEO. 

See also [`TaylorIntegration.lyap_taylorstep`](@ref).
"""
function lyap_neosstep!(f!, t::Taylor1{T}, x::Vector{Taylor1{U}}, dx::Vector{Taylor1{U}}, xaux::Vector{Taylor1{U}},
                        δx::Array{TaylorN{Taylor1{U}}, 1}, dδx::Array{TaylorN{Taylor1{U}}, 1}, jac::Array{Taylor1{U}, 2}, 
                        abstol::T, _δv::Vector{TaylorN{Taylor1{U}}}, varsaux::Array{Taylor1{U},3}, params, 
                        jacobianfunc! = nothing) where {T <: Real, U <: Number}

    # Dimensions of phase-space: dof
    nx = length(x)
    dof = length(δx)

    # Compute the Taylor coefficients associated to trajectory
    neos_jetcoeffs!(f!, t, view(x, 1:dof), view(dx, 1:dof), view(xaux, 1:dof), params)

    # Compute stability matrix
    TaylorIntegration.stabilitymatrix!(f!, t, x, δx, dδx, jac, _δv, params, jacobianfunc!)

    # Compute the Taylor coefficients associated to variational equations
    TaylorIntegration.lyap_jetcoeffs!(t, view(x, dof+1:nx), view(dx, dof+1:nx), jac, varsaux)
    
    # Compute the step-size of the integration using `abstol`
    δt = TaylorIntegration.stepsize(view(x, 1:dof), abstol)

    return δt
end

function lyap_neosstep!(f!, t::Taylor1{T}, x::Vector{Taylor1{U}}, dx::Vector{Taylor1{U}}, δx::Array{TaylorN{Taylor1{U}}, 1}, 
                        dδx::Array{TaylorN{Taylor1{U}}, 1}, jac::Array{Taylor1{U}, 2}, abstol::T, _δv::Vector{TaylorN{Taylor1{U}}},
                        varsaux::Array{Taylor1{U}, 3}, params, rv::TaylorIntegration.RetAlloc{Taylor1{U}}, 
                        jacobianfunc! = nothing) where {T <: Real, U <: Number}

    # Dimensions of phase-space: dof
    nx = length(x)
    dof = length(δx)

    # Compute the Taylor coefficients associated to trajectory
    TaylorIntegration.__jetcoeffs!(Val(true), f!, t, view(x, 1:dof), view(dx, 1:dof), params, rv)

    # Compute stability matrix
    TaylorIntegration.stabilitymatrix!(f!, t, x, δx, dδx, jac, _δv, params, jacobianfunc!)

    # Compute the Taylor coefficients associated to variational equations
    TaylorIntegration.lyap_jetcoeffs!(t, view(x, dof+1:nx), view(dx, dof+1:nx), jac, varsaux)
    
    # Compute the step-size of the integration using `abstol`
    δt = TaylorIntegration.stepsize(view(x, 1:dof), abstol)

    return δt
end

const V_true = :(Val{true})
const V_false = :(Val{false})
const V_true_false = (V_true, V_false)

@doc raw"""
    neosinteg(f!, q0::Array{U, 1}, t0::T, tmax::T, order::Int, abstol::T, ::$V, params = nothing;
              maxsteps::Int = 500, parse_eqs::Bool = true) where {T <: Real, U <: Number}
    neosinteg(f!, g, q0::Array{U, 1}, t0::T, tmax::T, order::Int, abstol::T, ::$V, params = nothing; 
              maxsteps::Int = 500, parse_eqs::Bool = true, eventorder::Int = 0, newtoniter::Int = 10, 
              nrabstol::T = eps(T)) where {T <: Real, U <: Number}

Specialized methods of `TaylorIntegration.taylorinteg` for the integration of a NEO.

See also [`TaylorIntegration.taylorinteg`](@ref). 
""" neosinteg

for V in V_true_false
    @eval begin

        function neosinteg(f!, q0::Array{U, 1}, t0::T, tmax::T, order::Int, abstol::T, ::$V, params = nothing;
                           maxsteps::Int = 500, parse_eqs::Bool = true) where {T <: Real, U <: Number}

            # Initialize the vector of Taylor1 expansions
            dof = length(q0)
            t = t0 + Taylor1( T, order )
            x = Array{Taylor1{U}}(undef, dof)
            dx = Array{Taylor1{U}}(undef, dof)
            @inbounds for i in eachindex(q0)
            @inbounds x[i] = Taylor1( q0[i], order )
            @inbounds dx[i] = Taylor1( zero(q0[i]), order )
            end

            # Determine if specialized jetcoeffs! method exists
            parse_eqs, rv = TaylorIntegration._determine_parsing!(parse_eqs, f!, t, x, dx, params)

            if parse_eqs
                # Re-initialize the Taylor1 expansions
                t = t0 + Taylor1( T, order )
                x .= Taylor1.( q0, order )
                return _neosinteg!(f!, t, x, dx, q0, t0, tmax, abstol, rv, $V(), params, maxsteps = maxsteps)
            else
                return _neosinteg!(f!, t, x, dx, q0, t0, tmax, abstol, $V(), params, maxsteps=maxsteps)
            end

        end

        function _neosinteg!(f!, t::Taylor1{T}, x::Array{Taylor1{U}, 1}, dx::Array{Taylor1{U}, 1}, q0::Array{U, 1}, t0::T, 
                             tmax::T, abstol::T, ::$V, params; maxsteps::Int = 500) where {T<:Real, U<:Number}

            # Initialize the vector of Taylor1 expansions
            dof = length(q0)

            # Allocation
            tv = Array{T}(undef, maxsteps+1)
            xv = Array{U}(undef, dof, maxsteps+1)
            if $V == Val{true}
                polynV = Array{Taylor1{U}}(undef, dof, maxsteps+1)
            end
            xaux = Array{Taylor1{U}}(undef, dof)

            # Initial conditions
            @inbounds t[0] = t0
            # x .= Taylor1.(q0, order)
            x0 = deepcopy(q0)
            @inbounds tv[1] = t0
            @inbounds xv[:,1] .= q0
            if $V == Val{true}
                @inbounds polynV[:,1] .= deepcopy.(x)
            end
            sign_tstep = copysign(1, tmax-t0)

            # Integration
            nsteps = 1
            while sign_tstep*t0 < sign_tstep*tmax
                δt = neosstep!(f!, t, x, dx, xaux, abstol, params) # δt is positive!
                # Below, δt has the proper sign according to the direction of the integration
                δt = sign_tstep * min(δt, sign_tstep*(tmax-t0))
                evaluate!(x, δt, x0) # new initial condition
                if $V == Val{true}
                    # Store the Taylor polynomial solution
                    @inbounds polynV[:,nsteps+1] .= deepcopy.(x)
                end
                @inbounds Threads.@threads for i in eachindex(x0)
                    x[i][0] = x0[i]
                    dx[i][0] = zero(x0[i])
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

            if $V == Val{true}
                return TaylorInterpolant(tv[1], view(tv.-tv[1],1:nsteps), view(transpose(view(polynV,:,2:nsteps)),1:nsteps-1,:))
            elseif $V == Val{false}
                return view(tv,1:nsteps), view(transpose(view(xv,:,1:nsteps)),1:nsteps,:)
            end
        end

        function _neosinteg!(f!, t::Taylor1{T}, x::Array{Taylor1{U},1}, dx::Array{Taylor1{U},1}, q0::Array{U,1}, t0::T, 
            tmax::T, abstol::T, rv::TaylorIntegration.RetAlloc{Taylor1{U}}, ::$V, params; maxsteps::Int=500) where {T<:Real, U<:Number}

            # Initialize the vector of Taylor1 expansions
            dof = length(q0)

            # Allocation of output
            tv = Array{T}(undef, maxsteps+1)
            xv = Array{U}(undef, dof, maxsteps+1)
            if $V == Val{true}
                polynV = Array{Taylor1{U}}(undef, dof, maxsteps+1)
            end

            # Initial conditions
            @inbounds t[0] = t0
            x0 = deepcopy(q0)
            @inbounds tv[1] = t0
            @inbounds xv[:,1] .= q0
            if $V == Val{true}
                @inbounds polynV[:,1] .= deepcopy.(x)
            end
            sign_tstep = copysign(1, tmax-t0)

            # Integration
            nsteps = 1
            while sign_tstep*t0 < sign_tstep*tmax
                δt = neosstep!(f!, t, x, dx, abstol, params, rv) # δt is positive!
                # Below, δt has the proper sign according to the direction of the integration
                δt = sign_tstep * min(δt, sign_tstep*(tmax-t0))
                evaluate!(x, δt, x0) # new initial condition
                if $V == Val{true}
                    # Store the Taylor polynomial solution
                    @inbounds polynV[:,nsteps+1] .= deepcopy.(x)
                end

                @inbounds Threads.@threads for i in eachindex(x0)
                    x[i][0] = x0[i]
                    dx[i][0] = zero(x0[i])
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

            if $V == Val{true}
                return TaylorInterpolant(tv[1], view(tv.-tv[1],1:nsteps), view(transpose(view(polynV,:,2:nsteps)),1:nsteps-1,:))
            elseif $V == Val{false}
                return view(tv,1:nsteps), view(transpose(view(xv,:,1:nsteps)),1:nsteps,:)
            end
        end

        # TODO: allow backwards integration
        function neosinteg(f!, g, q0::Array{U, 1}, t0::T, tmax::T, order::Int, abstol::T, ::$V, params = nothing; 
                           maxsteps::Int = 500, parse_eqs::Bool = true, eventorder::Int = 0, newtoniter::Int = 10, 
                           nrabstol::T = eps(T)) where {T <: Real, U <: Number}
    
            @assert order ≥ eventorder "`eventorder` must be less than or equal to `order`"

            # Initialize the vector of Taylor1 expansions
            dof = length(q0)
            t = t0 + Taylor1( T, order )
            x = Array{Taylor1{U}}(undef, dof)
            dx = Array{Taylor1{U}}(undef, dof)
            @inbounds for i in eachindex(q0)
                x[i] = Taylor1( q0[i], order )
                dx[i] = Taylor1( zero(q0[i]), order )
            end

            # Determine if specialized jetcoeffs! method exists
            parse_eqs, rv = TaylorIntegration._determine_parsing!(parse_eqs, f!, t, x, dx, params)

            if parse_eqs
                # Re-initialize the Taylor1 expansions
                t = t0 + Taylor1( T, order )
                x .= Taylor1.( q0, order )
                return _neosinteg!(f!, g, t, x, dx, q0, t0, tmax, abstol, rv, $V(), params, maxsteps = maxsteps, 
                                   eventorder = eventorder, newtoniter = newtoniter, nrabstol = nrabstol)
            else
                return _neosinteg!(f!, g, t, x, dx, q0, t0, tmax, abstol, $V(), params, maxsteps = maxsteps, 
                                   eventorder = eventorder, newtoniter = newtoniter, nrabstol = nrabstol)
            end

        end 

        function _neosinteg!(f!, g, t::Taylor1{T}, x::Array{Taylor1{U}, 1}, dx::Array{Taylor1{U}, 1}, q0::Array{U, 1}, t0::T, 
                             tmax::T, abstol::T, ::$V, params; maxsteps::Int = 500, eventorder::Int = 0, newtoniter::Int = 10, 
                             nrabstol::T = eps(T)) where {T <: Real, U <: Number}
    
            # Initialize the vector of Taylor1 expansions
            dof = length(q0)

            # Allocation
            tv = Array{T}(undef, maxsteps+1)
            xv = Array{U}(undef, dof, maxsteps+1)
            if $V == Val{true}
                polynV = Array{Taylor1{U}}(undef, dof, maxsteps+1)
            end
            xaux = Array{Taylor1{U}}(undef, dof)

            # Initial conditions
            @inbounds t[0] = t0
            # x .= Taylor1.(q0, order)
            x0 = deepcopy(q0)
            @inbounds tv[1] = t0
            @inbounds xv[:,1] .= q0
            if $V == Val{true}
                @inbounds polynV[:,1] .= deepcopy.(x)
            end
            sign_tstep = copysign(1, tmax-t0)

            # Some auxiliary arrays for root-finding/event detection/Poincaré surface of section evaluation
            g_tupl = g(dx, x, params, t)
            g_tupl_old = g(dx, x, params, t)
            δt = zero(x[1])
            δt_old = zero(x[1])

            x_dx = vcat(x, dx)
            g_dg = vcat(g_tupl[2], g_tupl_old[2])
            x_dx_val = Array{U}(undef, length(x_dx) )
            g_dg_val = vcat(evaluate(g_tupl[2]), evaluate(g_tupl_old[2]))

            tvS = Array{U}(undef, maxsteps+1)
            xvS = similar(xv)
            gvS = similar(tvS)

            # Integration
            nsteps = 1
            nevents = 1 #number of detected events
            while sign_tstep*t0 < sign_tstep*tmax
                δt_old = δt
                δt = neosstep!(f!, t, x, dx, xaux, abstol, params) # δt is positive!
                # Below, δt has the proper sign according to the direction of the integration
                δt = sign_tstep * min(δt, sign_tstep*(tmax-t0))
                evaluate!(x, δt, x0) # new initial condition
                g_tupl = g(dx, x, params, t)
                nevents = TaylorIntegration.findroot!(t, x, dx, g_tupl_old, g_tupl, eventorder, tvS, xvS, gvS, t0, 
                                                      δt_old, x_dx, x_dx_val, g_dg, g_dg_val, nrabstol, newtoniter, nevents)
                if $V == Val{true}
                    # Store the Taylor polynomial solution
                    @inbounds polynV[:,nsteps+1] .= deepcopy.(x)
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

            if $V == Val{true}
                return TaylorInterpolant(tv[1], view(tv .- tv[1], 1:nsteps), view(transpose(view(polynV, :, 2:nsteps)), 1:nsteps-1, :)),
                       view(tvS, 1:nevents-1), view(transpose(view(xvS, :, 1:nevents-1)), 1:nevents-1, :), view(gvS, 1:nevents-1)
            elseif $V == Val{false}
                return view(tv, 1:nsteps), view(transpose(view(xv, :, 1:nsteps)), 1:nsteps, :), view(tvS, 1:nevents-1), 
                view(transpose(view(xvS, :, 1:nevents-1)), 1:nevents-1, :), view(gvS, 1:nevents-1)
            end
            
        end 

        function _neosinteg!(f!, g, t::Taylor1{T}, x::Array{Taylor1{U}, 1}, dx::Array{Taylor1{U}, 1}, q0::Array{U, 1}, t0::T, 
                             tmax::T, abstol::T, rv::TaylorIntegration.RetAlloc{Taylor1{U}}, ::$V, params; maxsteps::Int = 500, 
                             eventorder::Int = 0, newtoniter::Int = 10, nrabstol::T = eps(T)) where {T <: Real, U <: Number}

            # Initialize the vector of Taylor1 expansions
            dof = length(q0)

            # Allocation
            tv = Array{T}(undef, maxsteps+1)
            xv = Array{U}(undef, dof, maxsteps+1)
            if $V == Val{true}
                polynV = Array{Taylor1{U}}(undef, dof, maxsteps+1)
            end
            
            # Initial conditions
            @inbounds t[0] = t0
            # x .= Taylor1.(q0, order)
            x0 = deepcopy(q0)
            @inbounds tv[1] = t0
            @inbounds xv[:,1] .= q0
            if $V == Val{true}
            @inbounds polynV[:,1] .= deepcopy.(x)
            end
            sign_tstep = copysign(1, tmax-t0)

            # Some auxiliary arrays for root-finding/event detection/Poincaré surface of section evaluation
            g_tupl = g(dx, x, params, t)
            g_tupl_old = g(dx, x, params, t)
            δt = zero(x[1])
            δt_old = zero(x[1])

            x_dx = vcat(x, dx)
            g_dg = vcat(g_tupl[2], g_tupl_old[2])
            x_dx_val = Array{U}(undef, length(x_dx) )
            g_dg_val = vcat(evaluate(g_tupl[2]), evaluate(g_tupl_old[2]))

            tvS = Array{U}(undef, maxsteps+1)
            xvS = similar(xv)
            gvS = similar(tvS)

            # Integration
            nsteps = 1
            nevents = 1 #number of detected events
            while sign_tstep*t0 < sign_tstep*tmax
                δt_old = δt
                δt = neosstep!(f!, t, x, dx, abstol, params, rv) # δt is positive!
                # Below, δt has the proper sign according to the direction of the integration
                δt = sign_tstep * min(δt, sign_tstep*(tmax-t0))
                evaluate!(x, δt, x0) # new initial condition
                g_tupl = g(dx, x, params, t)
                nevents = TaylorIntegration.findroot!(t, x, dx, g_tupl_old, g_tupl, eventorder, tvS, xvS, gvS, t0, 
                                                      δt_old, x_dx, x_dx_val, g_dg, g_dg_val, nrabstol, newtoniter, nevents)
                if $V == Val{true}
                    # Store the Taylor polynomial solution
                    @inbounds polynV[:,nsteps+1] .= deepcopy.(x)
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
    
            if $V == Val{true}
                return TaylorInterpolant(tv[1], view(tv .- tv[1], 1:nsteps), view(transpose(view(polynV, :, 2:nsteps)), 1:nsteps-1, :)),
                       view(tvS, 1:nevents-1), view(transpose(view(xvS, :, 1:nevents-1)), 1:nevents-1, :), view(gvS, 1:nevents-1)
            elseif $V == Val{false}
                return view(tv, 1:nsteps), view(transpose(view(xv, :, 1:nsteps)), 1:nsteps, :), view(tvS, 1:nevents-1), 
                       view(transpose(view(xvS, :, 1:nevents-1)), 1:nevents-1, :), view(gvS, 1:nevents-1)
            end

        end 

    end
end 

@doc raw"""
    lyap_neosinteg(f!, q0::Array{U, 1}, t0::T, tmax::T, order::Int, abstol::T, params = nothing, jacobianfunc! = nothing;
                   maxsteps::Int = 500, parse_eqs::Bool = true) where {T <: Real, U <: Number}

Specialized method of `TaylorIntegration.lyap_taylorinteg` for the calculation of the Lyapunov spectrum of a NEO.

See also [`TaylorIntegration.lyap_taylorinteg`](@ref). 
"""
function lyap_neosinteg(f!, q0::Array{U, 1}, t0::T, tmax::T, order::Int, abstol::T, params = nothing, jacobianfunc! = nothing;
                          maxsteps::Int = 500, parse_eqs::Bool = true) where {T <: Real, U <: Number}

    # Initialize the vector of Taylor1 expansions
    dof = length(q0)
    t = t0 + Taylor1( T, order )
    jt = Matrix{U}(I, dof, dof)
    x0 = vcat(q0, reshape(jt, dof*dof))
    nx0 = length(x0)
    x = Array{Taylor1{U}}(undef, nx0)
    dx = Array{Taylor1{U}}(undef, nx0)
    x .= Taylor1.( x0, order )
    # dx .= zero.(x)
    _dv = Array{TaylorN{Taylor1{U}}}(undef, dof)

    # If user does not provide Jacobian, check number of TaylorN variables and initialize _dv
    if isa(jacobianfunc!, Nothing)
        @assert get_numvars() == dof "`length(q0)` must be equal to number of variables set by `TaylorN`"
        for ind in eachindex(q0)
            _dv[ind] = one(x[1])*TaylorN(Taylor1{U}, ind, order=1)
        end
    end

    # Determine if specialized jetcoeffs! method exists
    parse_eqs, rv = TaylorIntegration._determine_parsing!(parse_eqs, f!, t, view(x, 1:dof), view(dx, 1:dof), params)

    if parse_eqs
        # Re-initialize the Taylor1 expansions
        t = t0 + Taylor1( T, order )
        x .= Taylor1.( x0, order )
        return _lyap_taylorinteg!(f!, t, x, dx, q0, t0, tmax, abstol, jt, _dv, rv, params, jacobianfunc!, maxsteps = maxsteps)
    else
        return _lyap_taylorinteg!(f!, t, x, dx, q0, t0, tmax, abstol, jt, _dv, params, jacobianfunc!, maxsteps = maxsteps)
    end
end

function _lyap_taylorinteg!(f!, t::Taylor1{T}, x::Array{Taylor1{U},1}, dx::Array{Taylor1{U},1}, q0::Array{U,1}, t0::T, tmax::T, 
                            abstol::T, jt::Matrix{U}, _δv::Array{TaylorN{Taylor1{U}},1}, params, jacobianfunc!; maxsteps::Int = 500) where {T <: Real, U <: Number}

    # Allocation
    order = get_order(t)
    tv = Array{T}(undef, maxsteps+1)
    dof = length(q0)
    nx0 = length(x)
    x0 = getcoeff.(x, 0)
    xv = Array{U}(undef, dof, maxsteps+1)
    λ = similar(xv)
    λtsum = similar(q0)
    # q1 = similar(q0)

    # Initial conditions
    @inbounds t[0] = t0
    sign_tstep = copysign(1, tmax-t0)
    t00 = t0
    tspan = zero(T)
    @inbounds tv[1] = t0
    @inbounds for ind in eachindex(q0)
        xv[ind,1] = q0[ind]
        λ[ind,1] = zero(U)
        λtsum[ind] = zero(U)
    end

    # Allocate auxiliary arrays
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

    # Integration
    nsteps = 1
    while sign_tstep*t0 < sign_tstep*tmax
        δt = lyap_neosstep!(f!, t, x, dx, xaux, δx, dδx, jac, abstol, _δv, varsaux, params, jacobianfunc!) # δt is positive!
        # Below, δt has the proper sign according to the direction of the integration
        δt = sign_tstep * min(δt, sign_tstep*(tmax-t0))
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

    return view(tv, 1:nsteps),  view(transpose(xv), 1:nsteps, :),  view(transpose(λ), 1:nsteps, :)
end

function _lyap_taylorinteg!(f!, t::Taylor1{T}, x::Array{Taylor1{U},1}, dx::Array{Taylor1{U},1},
    q0::Array{U,1}, t0::T, tmax::T, abstol::T, jt::Matrix{U}, _δv::Array{TaylorN{Taylor1{U}},1},
    rv::TaylorIntegration.RetAlloc{Taylor1{U}}, params, jacobianfunc!; maxsteps::Int=500) where {T<:Real, U<:Number}

    # Allocation
    order = get_order(t)
    tv = Array{T}(undef, maxsteps+1)
    dof = length(q0)
    x0 = getcoeff.(x, 0)
    xv = Array{U}(undef, dof, maxsteps+1)
    λ = similar(xv)
    λtsum = similar(q0)

    # Initial conditions
    @inbounds t[0] = t0
    sign_tstep = copysign(1, tmax-t0)
    t00 = t0
    tspan = zero(T)
    @inbounds tv[1] = t0
    @inbounds for ind in eachindex(q0)
        xv[ind,1] = q0[ind]
        λ[ind,1] = zero(U)
        λtsum[ind] = zero(U)
    end

    # Allocate auxiliary arrays
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

    # Integration
    nsteps = 1
    while sign_tstep*t0 < sign_tstep*tmax
        δt = lyap_neosstep!(f!, t, x, dx, δx, dδx, jac, abstol, _δv, varsaux, params, rv, jacobianfunc!) # δt is positive!
        # Below, δt has the proper sign according to the direction of the integration
        δt = sign_tstep * min(δt, sign_tstep*(tmax-t0))
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

    return view(tv, 1:nsteps),  view(transpose(xv), 1:nsteps, :),  view(transpose(λ), 1:nsteps, :)
end