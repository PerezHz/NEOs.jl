const V_true = :(Val{true})
const V_false = :(Val{false})
const V_true_false = (V_true, V_false)

@doc raw"""
    neosinteg(f!, g, q0::Array{U, 1}, t0::T, tmax::T, order::Int, abstol::T, Val(true/false), params = nothing; 
              maxsteps::Int = 500, parse_eqs::Bool = true, eventorder::Int = 0, newtoniter::Int = 10, 
              nrabstol::T = eps(T)) where {T <: Real, U <: Number}

Specialized method of `TaylorIntegration.taylorinteg` for the integration of a NEO.

See also [`TaylorIntegration.taylorinteg`](@ref). 
""" neosinteg

for V in V_true_false
    @eval begin

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
                δt = TaylorIntegration.taylorstep!(f!, t, x, dx, xaux, abstol, params) # δt is positive!
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
                jd0 = params[4]
                return TaylorInterpolant(jd0 - JD_J2000, view(tv .- tv[1], 1:nsteps), view(transpose(view(polynV, :, 2:nsteps)), 1:nsteps-1, :)),
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
                δt = TaylorIntegration.taylorstep!(f!, t, x, dx, abstol, params, rv) # δt is positive!
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
                jd0 = params[4]
                return TaylorInterpolant(jd0-JD_J2000, view(tv .- tv[1], 1:nsteps), view(transpose(view(polynV, :, 2:nsteps)), 1:nsteps-1, :)),
                       view(tvS, 1:nevents-1), view(transpose(view(xvS, :, 1:nevents-1)), 1:nevents-1, :), view(gvS, 1:nevents-1)
            elseif $V == Val{false}
                return view(tv, 1:nsteps), view(transpose(view(xv, :, 1:nsteps)), 1:nsteps, :), view(tvS, 1:nevents-1), 
                       view(transpose(view(xvS, :, 1:nevents-1)), 1:nevents-1, :), view(gvS, 1:nevents-1)
            end

        end 

    end
end 