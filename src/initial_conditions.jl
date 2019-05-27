function initialcond()
    # output from JPL Horizons

    q0 = Array{Float64}(undef, 7) #initial condition array

    #99942 apophis
    q0[1:3] .= [-0.9633018148154301, 0.5100291398744978, 0.1652803001584807]
    q0[4:6] .= [-0.007118874645519014, -0.01206123416050721, -0.00466951380149129]

    q0[7] = 0.0 # A2 Yarkovsky coefficient

    return q0
end
