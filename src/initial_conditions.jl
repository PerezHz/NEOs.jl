function initialcond(dq::Vector=zeros(7); quadmath::Bool=false)
    if quadmath
        q0 = Vector{Float128}(undef, 7) #initial condition array
    else
        q0 = Vector{Float64}(undef, 7) #initial condition array
    end
    # # output from JPL solution #199 SPK file for Apophis via SPICE.jl
    # q0[1:3] .= [-0.96330181481543, 0.5100291398744977, 0.16528030015848072]
    # q0[4:6] .= [-0.007118874645519018, -0.0120612341605072, -0.004669513801491289]
    # output from JPL solution #197 SPK file for Apophis via SPICE.jl
    q0[1:3] .= [-0.9633018953468989, 0.5100289806011301, 0.16528024397505386]
    q0[4:6] .= [-0.0071188720570829036, -0.012061235689040535, -0.0046695143453363164]
    q0[7] = 0.0 # A2 Yarkovsky coefficient
    S = typeof(q0[1]+dq[1])
    q = Vector{S}(undef, 7)
    for i in 1:6
        q[i] = q0[i] + dq[i]
    end
    if length(dq) == 7
        q[7] = q0[7] + dq[7]
    elseif length(dq) == 6
        q[7] = q0[7] + 0dq[1]
    end
    return q
end
