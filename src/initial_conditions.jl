function initialcond(dq::Vector=zeros(7))
    q0 = Vector{Float64}(undef, 7) #initial condition array
    # output from JPL solution #199 SPK file for Apophis via SPICE.jl
    q0[1:3] .= [-0.963307185415711, 0.5100200405384891, 0.16527677734731944]
    q0[4:6] .= [-0.007118712801283596, -0.012061319251632138, -0.004669541303113223]
    # # output from JPL solution #197 SPK file for Apophis via SPICE.jl
    # q0[1:3] .= [-0.9633072659452289, 0.5100198812639655, 0.1652767211634811]
    # q0[4:6] .= [-0.007118710212828963, -0.012061320780141318, -0.004669541846949655]
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
