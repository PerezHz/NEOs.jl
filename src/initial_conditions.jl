function initialcond(N)
    # output from JPL Horizons

    q0 = Array{Float64}(undef, 6N+1) #initial condition array

    #sun
    i = 1
    q0[3i-2:3i]         .= [-0.001540802671596496, 0.00438889278161424, 0.001858638798427623]
    q0[3(N+i)-2:3(N+i)] .= [-6.109751994756004e-6, -1.931572052779964e-6, -7.092564632948912e-7]
    #mercury
    i = 2
    q0[3i-2:3i]         .= [0.3046385524098975, -0.2239262872848795, -0.1518487166013703]
    q0[3(N+i)-2:3(N+i)] .= [0.013061958045048, 0.0203012991458782, 0.009489733386707604]
    #venus
    i = 3
    q0[3i-2:3i]         .= [-0.2370720145769874, -0.627040738077614, -0.267324380799051]
    q0[3(N+i)-2:3(N+i)] .= [0.01898821199448838, -0.005631644544425267, -0.003735536566147577]
    #earth
    i = 4
    q0[3i-2:3i]         .= [1.001401996348207, 0.02359641338688291, 0.01018195925304371]
    q0[3(N+i)-2:3(N+i)] .= [-0.0006391556078946668, 0.01572260216573083, 0.006817017673265036]
    #moon
    i = 5
    q0[3i-2:3i]         .= [1.000367432977523, 0.02564117583042711, 0.01116491651187949]
    q0[3(N+i)-2:3(N+i)] .= [-0.001195603426398212, 0.01552432630673204, 0.006676849606937177]
    #mars
    i = 6
    q0[3i-2:3i]         .= [-1.270960531462236, -0.8678508302720224, -0.3639208032155261]
    q0[3(N+i)-2:3(N+i)] .= [0.008879230398263355, -0.009029298818909943, -0.004381331767031147]
    #jupiter
    i = 7
    q0[3i-2:3i]         .= [2.077128040651103, -4.307585374501212, -1.896991978186166]
    q0[3(N+i)-2:3(N+i)] .= [0.006810526029641702, 0.003190231984051103, 0.001201585413837261]
    #saturn
    i = 8
    q0[3i-2:3i]         .= [-8.903241825277506, 2.458572188095974, 1.398754563334379]
    q0[3(N+i)-2:3(N+i)] .= [-0.001977534413001538, -0.004960690616673099, -0.001963830282375063]
    #uranus
    i = 9
    q0[3i-2:3i]         .= [19.82869761020429, -2.880418477382511, -1.542003907582929]
    q0[3(N+i)-2:3(N+i)] .= [0.0006086573593461934, 0.003390327455820016, 0.001476271782483529]
    #neptune
    i = 10
    q0[3i-2:3i]         .= [23.97437761296793, -16.52272563943745, -7.359720948896527]
    q0[3(N+i)-2:3(N+i)] .= [0.001869835257921788, 0.002353940265223347, 0.000916928221659154]
    #ceres
    i = 11
    q0[3i-2:3i]         .= [-1.208270526225632, 1.963641380708032, 1.170481449363465]
    q0[3(N+i)-2:3(N+i)] .= [-0.009348216591255532, -0.005863175034617146, -0.0008596609867992848]
    #apophis
    i = 12
    q0[3i-2:3i]         .= [-0.9633018148154301, 0.5100291398744978, 0.1652803001584807]
    q0[3(N+i)-2:3(N+i)] .= [-0.007118874645519014, -0.01206123416050721, -0.00466951380149129]

    q0[6N+1] = 0.0 # A2 Yarkovsky coefficient

    q1 = q0=[-0.001540802671596496, 0.00438889278161424, 0.001858638798427623, 0.3046385524098975, -0.2239262872848795, -0.1518487166013703, -0.2370720145769874, -0.627040738077614, -0.267324380799051, 1.001401996348207, 0.02359641338688291, 0.01018195925304371, 1.000367432977523, 0.02564117583042711, 0.01116491651187949, -1.270960531462236, -0.8678508302720224, -0.3639208032155261, 2.077128040651103, -4.307585374501212, -1.896991978186166, -8.903241825277506, 2.458572188095974, 1.398754563334379, 19.82869761020429, -2.880418477382511, -1.542003907582929, 23.97437761296793, -16.52272563943745, -7.359720948896527, -1.208270526225632, 1.963641380708032, 1.170481449363465, -0.9633018148154301, 0.5100291398744978, 0.1652803001584807, -6.109751994756004e-6, -1.931572052779964e-6, -7.092564632948912e-7, 0.013061958045048, 0.0203012991458782, 0.009489733386707604, 0.01898821199448838, -0.005631644544425267, -0.003735536566147577, -0.0006391556078946668, 0.01572260216573083, 0.006817017673265036, -0.001195603426398212, 0.01552432630673204, 0.006676849606937177, 0.008879230398263355, -0.009029298818909943, -0.004381331767031147, 0.006810526029641702, 0.003190231984051103, 0.001201585413837261, -0.001977534413001538, -0.004960690616673099, -0.001963830282375063, 0.0006086573593461934, 0.003390327455820016, 0.001476271782483529, 0.001869835257921788, 0.002353940265223347, 0.000916928221659154, -0.009348216591255532, -0.005863175034617146, -0.0008596609867992848, -0.007118874645519014, -0.01206123416050721, -0.00466951380149129, 0.0]
    @show q0 == q1

    return q0
end
