# Nearth-Earth asteroid dynamical model (d=2.0)
# Bodies considered in the model are: the Sun, the eight planets, the Moon and Ceres,
# as well as the asteroid of interest as a test particle with null mass. Dynamical
# effects considered are:
# - post-Newtonian point-mass accelerations between all bodies,
# - figure-effects (oblateness) of the Earth (J2 and J3)
# - J2 effect of the Sun
# - J2 and J3 effect of the Moon
# - Kinematic model for the precession and nutation of the Earth's orientation (IAU 1976/1980 Earth orientation model)
# - Kinematic model for the Moons's orientation (Seidelmann et al., 2006)
# - Non-gravitational accelerations acting upon the asteroid
# are included (Yarkovsky effect) a_nongrav = A2*t_vec*(au/r)^d, where t_vec is the
# unit heliocentric transverse vector, au is 1 astronomical unit, r is the
# asteroid's heliocentric range, A2 is a coefficient (with units of au/day^2),
# and d = 2.0
@taylorize function RNp1BP_pN_A_J23E_J2S_ng_eph!(dq, q, params, t)
    local ss16asteph_t = params[1](t)*one(q[1]) #ss16asteph(t)
    local acceph_t = params[2](t)*one(q[1]) #acc_eph(t)
    local newtonianNb_Potential_t = params[3](t)*one(q[1]) #newtonianNb_Potential(t), massive bodies
    local S = eltype(q[1])
    local N = length(μ) # number of bodies, including NEA
    local _1_to_N = Base.OneTo(N) # iterator over all bodies

    local succ_approx_iter = 1 # number of iterations of post-Newtonian subroutine
    # local j2_body_index = [su, ea]#, mo] # indices of bodies with J2 flattening (note: Earth and Moon also have J3)

    # parameters related to speed of light, c
    local c_p2 = 29979.063823897606 # c^2 = 29979.063823897606 au^2/d^2
    local c_m2 = 3.3356611996764786e-5 # c^-2 = 3.3356611996764786e-5 d^2/au^2

    local zero_q_1 = zero(q[1])

    X = Array{Taylor1{S}}(undef, N)
    Y = Array{Taylor1{S}}(undef, N)
    Z = Array{Taylor1{S}}(undef, N)

    r_p2 = Array{Taylor1{S}}(undef, N)
    r_p1d2 = Array{Taylor1{S}}(undef, N)
    r_p3d2 = Array{Taylor1{S}}(undef, N)
    r_p7d2 = Array{Taylor1{S}}(undef, N)

    newtonianCoeff = Array{Taylor1{S}}(undef, N)

    xi = Array{Taylor1{S}}(undef, N-1)
    yi = Array{Taylor1{S}}(undef, N-1)
    zi = Array{Taylor1{S}}(undef, N-1)
    ui = Array{Taylor1{S}}(undef, N-1)
    vi = Array{Taylor1{S}}(undef, N-1)
    wi = Array{Taylor1{S}}(undef, N-1)
    _3ui = Array{Taylor1{S}}(undef, N-1)
    _3vi = Array{Taylor1{S}}(undef, N-1)
    _3wi = Array{Taylor1{S}}(undef, N-1)
    pNxi = Array{Taylor1{S}}(undef, N-1)
    pNyi = Array{Taylor1{S}}(undef, N-1)
    pNzi = Array{Taylor1{S}}(undef, N-1)

    # #post-Newtonian stuff
    U = Array{Taylor1{S}}(undef, N)
    V = Array{Taylor1{S}}(undef, N)
    W = Array{Taylor1{S}}(undef, N)

    _4U_m_3X = Array{Taylor1{S}}(undef, N)
    _4V_m_3Y = Array{Taylor1{S}}(undef, N)
    _4W_m_3Z = Array{Taylor1{S}}(undef, N)

    UU = Array{Taylor1{S}}(undef, N)
    VV = Array{Taylor1{S}}(undef, N)
    WW = Array{Taylor1{S}}(undef, N)

    postNewtonX = zero_q_1
    postNewtonY = zero_q_1
    postNewtonZ = zero_q_1

    newtonianNb_Potential = Array{Taylor1{S}}(undef, N)
    newtonian1b_Potential = Array{Taylor1{S}}(undef, N)
    newtonianCoeff = Array{Taylor1{S}}(undef, N)

    pntempX = Array{Taylor1{S}}(undef, N-1)
    pntempY = Array{Taylor1{S}}(undef, N-1)
    pntempZ = Array{Taylor1{S}}(undef, N-1)

    pn1 = Array{Taylor1{S}}(undef, N)
    v2 = Array{Taylor1{S}}(undef, N)
    vi_dot_vj = Array{Taylor1{S}}(undef, N)
    pn2 = Array{Taylor1{S}}(undef, N)
    pn3 = Array{Taylor1{S}}(undef, N)

    # J2 acceleration auxiliaries
    t31 = Array{Taylor1{S}}(undef, N)
    t32 = Array{Taylor1{S}}(undef, N)
    t33 = Array{Taylor1{S}}(undef, N)
    r_sin_ϕ = Array{Taylor1{S}}(undef, N)
    F_J2_x = Array{Taylor1{S}}(undef, N)
    F_J2_y = Array{Taylor1{S}}(undef, N)
    F_J2_z = Array{Taylor1{S}}(undef, N)
    F_J2_x1 = Array{Taylor1{S}}(undef, N)
    F_J2_y1 = Array{Taylor1{S}}(undef, N)
    F_J2_z1 = Array{Taylor1{S}}(undef, N)
    F_J2_x2 = Array{Taylor1{S}}(undef, N)
    F_J2_y2 = Array{Taylor1{S}}(undef, N)
    F_J2_z2 = Array{Taylor1{S}}(undef, N)
    temp_accX_i = Array{Taylor1{S}}(undef, N)
    temp_accY_i = Array{Taylor1{S}}(undef, N)
    temp_accZ_i = Array{Taylor1{S}}(undef, N)
    sin_ϕ = Array{Taylor1{S}}(undef, N)
    sin2_ϕ = Array{Taylor1{S}}(undef, N)
    sin3_ϕ = Array{Taylor1{S}}(undef, N)
    sin4_ϕ = Array{Taylor1{S}}(undef, N)
    ϕ = Array{Taylor1{S}}(undef, N)
    cos_ϕ = Array{Taylor1{S}}(undef, N)
    P_2_sin_ϕ = Array{Taylor1{S}}(undef, N)
    ∂P_2_sin_ϕ = Array{Taylor1{S}}(undef, N)
    P_3_sin_ϕ = Array{Taylor1{S}}(undef, N)
    ∂P_3_sin_ϕ = Array{Taylor1{S}}(undef, N)
    m_c_ϕ_∂P_2 = Array{Taylor1{S}}(undef, N)
    m_c_ϕ_∂P_3 = Array{Taylor1{S}}(undef, N)
    Λ2j_div_r4 = Array{Taylor1{S}}(undef, N)
    Λ3j_div_r5 = Array{Taylor1{S}}(undef, N)
    F_J_ξ = Array{Taylor1{S}}(undef, N)
    F_J_η = Array{Taylor1{S}}(undef, N)
    F_J_ζ = Array{Taylor1{S}}(undef, N)
    F_J2_ξ = Array{Taylor1{S}}(undef, N)
    F_J2_η = Array{Taylor1{S}}(undef, N)
    F_J2_ζ = Array{Taylor1{S}}(undef, N)
    F_J3_ξ = Array{Taylor1{S}}(undef, N)
    F_J3_η = Array{Taylor1{S}}(undef, N)
    F_J3_ζ = Array{Taylor1{S}}(undef, N)
    ξx = Array{Taylor1{S}}(undef, N)
    ξy = Array{Taylor1{S}}(undef, N)
    ξz = Array{Taylor1{S}}(undef, N)
    ηx = Array{Taylor1{S}}(undef, N)
    ηy = Array{Taylor1{S}}(undef, N)
    ηz = Array{Taylor1{S}}(undef, N)
    ηx1 = Array{Taylor1{S}}(undef, N)
    ηy1 = Array{Taylor1{S}}(undef, N)
    ηz1 = Array{Taylor1{S}}(undef, N)
    ηx2 = Array{Taylor1{S}}(undef, N)
    ηy2 = Array{Taylor1{S}}(undef, N)
    ηz2 = Array{Taylor1{S}}(undef, N)
    ζx = Array{Taylor1{S}}(undef, N)
    ζy = Array{Taylor1{S}}(undef, N)
    ζz = Array{Taylor1{S}}(undef, N)
    ζx1 = Array{Taylor1{S}}(undef, N)
    ζy1 = Array{Taylor1{S}}(undef, N)
    ζz1 = Array{Taylor1{S}}(undef, N)
    ζx2 = Array{Taylor1{S}}(undef, N)
    ζy2 = Array{Taylor1{S}}(undef, N)
    ζz2 = Array{Taylor1{S}}(undef, N)

    # # extended-body accelerations
    accX = zero_q_1
    accY = zero_q_1
    accZ = zero_q_1

    # # rotations to and from Earth, Sun and Moon pole-oriented frames
    local dsj2k = t-2.451545e6 # J2000.0 = 2.451545e6
    local αs = deg2rad(α_p_sun*one(t))
    local δs = deg2rad(δ_p_sun*one(t))
    local αm = moon_pole_ra(dsj2k)
    local δm = moon_pole_dec(dsj2k)
    local M_ = Array{Taylor1{S}}(undef, 3, 3, N)
    local M_[:,:,ea] = t2c_jpl_de430(dsj2k)
    local M_[:,:,su] = pole_rotation( αs, δs )
    # local M_[:,:,mo] = pole_rotation( αm, δm )

    dq[1] = q[4]
    dq[2] = q[5]
    dq[3] = q[6]

    newtonianNb_Potential[N] = zero_q_1

    for j in j2_body_index
        t31[j] = zero_q_1
        t32[j] = zero_q_1
        t33[j] = zero_q_1
        F_J2_x[j] = zero_q_1
        F_J2_y[j] = zero_q_1
        F_J2_z[j] = zero_q_1
        F_J2_x1[j] = zero_q_1
        F_J2_y1[j] = zero_q_1
        F_J2_z1[j] = zero_q_1
        F_J2_x2[j] = zero_q_1
        F_J2_y2[j] = zero_q_1
        F_J2_z2[j] = zero_q_1
        temp_accX_i[j] = zero_q_1
        temp_accY_i[j] = zero_q_1
        temp_accZ_i[j] = zero_q_1
        sin_ϕ[j] = zero_q_1
        sin2_ϕ[j] = zero_q_1
        sin3_ϕ[j] = zero_q_1
        sin4_ϕ[j] = zero_q_1
        ϕ[j] = zero_q_1
        cos_ϕ[j] = zero_q_1
        P_2_sin_ϕ[j] = zero_q_1
        ∂P_2_sin_ϕ[j] = zero_q_1
        P_3_sin_ϕ[j] = zero_q_1
        ∂P_3_sin_ϕ[j] = zero_q_1
        m_c_ϕ_∂P_2[j] = zero_q_1
        m_c_ϕ_∂P_3[j] = zero_q_1
        Λ2j_div_r4[j] = zero_q_1
        Λ3j_div_r5[j] = zero_q_1
        F_J_ξ[j] = zero_q_1
        F_J_η[j] = zero_q_1
        F_J_ζ[j] = zero_q_1
        F_J2_ξ[j] = zero_q_1
        F_J2_η[j] = zero_q_1
        F_J2_ζ[j] = zero_q_1
        F_J3_ξ[j] = zero_q_1
        F_J3_η[j] = zero_q_1
        F_J3_ζ[j] = zero_q_1
        ξx[j] = zero_q_1
        ξy[j] = zero_q_1
        ξz[j] = zero_q_1
        ηx[j] = zero_q_1
        ηy[j] = zero_q_1
        ηz[j] = zero_q_1
        ηx1[j] = zero_q_1
        ηy1[j] = zero_q_1
        ηz1[j] = zero_q_1
        ηx2[j] = zero_q_1
        ηy2[j] = zero_q_1
        ηz2[j] = zero_q_1
        ζx[j] = zero_q_1
        ζy[j] = zero_q_1
        ζz[j] = zero_q_1
        ζx1[j] = zero_q_1
        ζy1[j] = zero_q_1
        ζz1[j] = zero_q_1
        ζx2[j] = zero_q_1
        ζy2[j] = zero_q_1
        ζz2[j] = zero_q_1
    end #for j in j2_body_index

    #compute point-mass Newtonian accelerations, all bodies
    for i in 1:N-1
        xi[i] = ss16asteph_t[3i-2]
        yi[i] = ss16asteph_t[3i-1]
        zi[i] = ss16asteph_t[3i  ]
        ui[i] = ss16asteph_t[3(N-1+i)-2]
        vi[i] = ss16asteph_t[3(N-1+i)-1]
        wi[i] = ss16asteph_t[3(N-1+i)  ]
        _3ui[i] = 3ui[i]
        _3vi[i] = 3vi[i]
        _3wi[i] = 3wi[i]

        X[i] = xi[i]-q[1]
        Y[i] = yi[i]-q[2]
        Z[i] = zi[i]-q[3]

        U[i] = ui[i]-dq[1]
        V[i] = vi[i]-dq[2]
        W[i] = wi[i]-dq[3]

        _4U_m_3X[i] = (4dq[1]) - _3ui[i]
        _4V_m_3Y[i] = (4dq[2]) - _3vi[i]
        _4W_m_3Z[i] = (4dq[3]) - _3wi[i]

        pn2x = X[i]*_4U_m_3X[i]
        pn2y = Y[i]*_4V_m_3Y[i]
        pn2z = Z[i]*_4W_m_3Z[i]

        UU[i] = ui[i]*dq[1]
        VV[i] = vi[i]*dq[2]
        WW[i] = wi[i]*dq[3]

        vi_dot_vj[i] = ( UU[i]+VV[i] ) + WW[i]

        r_p2[i] = ( (X[i]^2)+(Y[i]^2) ) + (Z[i]^2)

        r_p1d2[i] = sqrt(r_p2[i])
        r_p3d2[i] = r_p2[i]^1.5
        r_p7d2[i] = r_p2[i]^3.5

        newtonianCoeff[i] =  μ[i]/r_p3d2[i]

        pn2[i] = newtonianCoeff[i]*(( pn2x+pn2y ) + pn2z)

        newtonian1b_Potential[i] = μ[i]/r_p1d2[i]
        #newtonianNb_Potential[i] = newtonianNb_Potential_t[i]

        #J2 accelerations, if i-th body is flattened
        if UJ_interaction[i]
            # rotate from inertial frame to extended-body frame
            t31[i] = -X[i]*M_[1,3,i]
            t32[i] = -Y[i]*M_[2,3,i]
            t33[i] = -Z[i]*M_[3,3,i]
            r_sin_ϕ[i] = (t31[i]+t32[i])+t33[i]

            # compute cartesian coordinates of acceleration due to body figure in body frame
            sin_ϕ[i] = r_sin_ϕ[i]/r_p1d2[i] # z/r = latitude of point mass wrt body frame
            ϕ[i] = asin(sin_ϕ[i])
            cos_ϕ[i] = cos(ϕ[i])
            sin2_ϕ[i] = sin_ϕ[i]^2
            sin3_ϕ[i] = sin_ϕ[i]^3
            P_2_sin_ϕ[i] = 1.5sin2_ϕ[i] - 0.5
            ∂P_2_sin_ϕ[i] = 3sin_ϕ[i]
            P_3_sin_ϕ[i] = (-1.5sin_ϕ[i]) + (2.5sin3_ϕ[i])
            ∂P_3_sin_ϕ[i] = -1.5 + 7.5sin2_ϕ[i]
            Λ2j_div_r4[i] = (-Λ2[i])/(r_p2[i]^2)
            Λ3j_div_r5[i] = (-Λ3[i])/(r_p1d2[i]^5)
            m_c_ϕ_∂P_2[i] = (-cos_ϕ[i])*∂P_2_sin_ϕ[i]
            m_c_ϕ_∂P_3[i] = (-cos_ϕ[i])*∂P_3_sin_ϕ[i]
            F_J2_ξ[i] = ( Λ2j_div_r4[i]*(3P_2_sin_ϕ[i]) )
            #F_J2_η[i] = zero_q_1
            F_J2_ζ[i] = Λ2j_div_r4[i]*m_c_ϕ_∂P_2[i]
            F_J3_ξ[i] = ( Λ3j_div_r5[i]*(4P_3_sin_ϕ[i]) )
            #F_J3_η[i] = zero_q_1
            F_J3_ζ[i] = Λ3j_div_r5[i]*m_c_ϕ_∂P_3[i]
            F_J_ξ[i] = F_J2_ξ[i] + F_J3_ξ[i]
            #F_J_η[i] = zero_q_1
            F_J_ζ[i] = F_J2_ζ[i] + F_J3_ζ[i]
            #Compute unit vectors ξ,η,ζ
            ξx[i] = -X[i]/r_p1d2[i]
            ξy[i] = -Y[i]/r_p1d2[i]
            ξz[i] = -Z[i]/r_p1d2[i]
            #Compute η = p x ξ
            ηx1[i] = M_[2,3,i]*ξz[i]
            ηy1[i] = M_[3,3,i]*ξx[i]
            ηz1[i] = M_[1,3,i]*ξy[i]
            ηx2[i] = M_[3,3,i]*ξy[i]
            ηy2[i] = M_[1,3,i]*ξz[i]
            ηz2[i] = M_[2,3,i]*ξx[i]
            ηx[i] = ηx1[i] - ηx2[i]
            ηy[i] = ηy1[i] - ηy2[i]
            ηz[i] = ηz1[i] - ηz2[i]
            #Compute ζ = ξ x η
            ζx1[i] = ξy[i]*ηz[i]
            ζy1[i] = ξz[i]*ηx[i]
            ζz1[i] = ξx[i]*ηy[i]
            ζx2[i] = ξz[i]*ηy[i]
            ζy2[i] = ξx[i]*ηz[i]
            ζz2[i] = ξy[i]*ηx[i]
            ζx[i] = ζx1[i] - ζx2[i]
            ζy[i] = ζy1[i] - ζy2[i]
            ζz[i] = ζz1[i] - ζz2[i]
            # compute cartesian coordinates of acceleration due to body figure in inertial frame
            F_J2_x1[i] = F_J_ξ[i]*ξx[i]
            F_J2_y1[i] = F_J_ξ[i]*ξy[i]
            F_J2_z1[i] = F_J_ξ[i]*ξz[i]
            F_J2_x2[i] = F_J_ζ[i]*ζx[i]
            F_J2_y2[i] = F_J_ζ[i]*ζy[i]
            F_J2_z2[i] = F_J_ζ[i]*ζz[i]
            F_J2_x[i] = F_J2_x1[i] + F_J2_x2[i]
            F_J2_y[i] = F_J2_y1[i] + F_J2_y2[i]
            F_J2_z[i] = F_J2_z1[i] + F_J2_z2[i]

            # # add result to total acceleration on upon j-th body figure due to i-th point mass
            # @show "acc",j,"+μ",i,"Λ2",j
            # temp_accX_j[i] = accX + (μ[i]*F_J2_x[i])
            # accX = temp_accX_j[i]
            # temp_accY_j[i] = accY + (μ[i]*F_J2_y[i])
            # accY = temp_accY_j[i]
            # temp_accZ_j[i] = accZ + (μ[i]*F_J2_z[i])
            # accZ = temp_accZ_j[i]

            # # reaction force on asteroid point mass
            # @show "acc",i,"-μ",j,"Λ2",j
            temp_accX_i[i] = -(μ[i]*F_J2_x[i])
            temp_accY_i[i] = -(μ[i]*F_J2_y[i])
            temp_accZ_i[i] = -(μ[i]*F_J2_z[i])
        end
        v2[i] = ( (ui[i]^2)+(vi[i]^2) ) + (wi[i]^2)
    end #for, i
    v2[N] = ( (q[4]^2)+(q[5]^2) ) + (q[6]^2)
    for i in 1:N-1
        newtonianNb_Potential[N] = newtonianNb_Potential[N] + newtonian1b_Potential[i]
    end
    for i in j2_body_index
        accX = accX + temp_accX_i[i]
        accY = accY + temp_accY_i[i]
        accZ = accZ + temp_accZ_i[i]
    end

    for k in Base.OneTo(succ_approx_iter)
        for i in 1:N-1
            #post-Newtonian corrections to gravitational acceleration
            #Moyer, 1971, page 7 eq. 35
            # acceleration of i-th body
            pNxi[i] = acceph_t[3i-2]
            pNyi[i] = acceph_t[3i-1]
            pNzi[i] = acceph_t[3i  ]

            temp_005a = newtonianNb_Potential_t[i] + (4newtonianNb_Potential[N])
            temp_005b = ( (2v2[i]) - (4vi_dot_vj[i]) ) + v2[N]
            temp_005 = temp_005b - temp_005a

            temp_006a = X[i]*ui[i]
            temp_006b = Y[i]*vi[i]
            temp_006c = Z[i]*wi[i]
            temp_006d = ( temp_006a+temp_006b ) + temp_006c
            # the expression below inside the (...)^2 should have a minus sign in front of the numerator,
            # but upon squaring it is eliminated, so at the end of the day, it is irrelevant ;)
            temp_006e = (temp_006d^2)/r_p2[i]
            temp_006 = 1.5temp_006e

            temp_007a = X[i]*pNxi[i]
            temp_007b = Y[i]*pNyi[i]
            temp_007c = Z[i]*pNzi[i]
            temp_007d = ( temp_007a+temp_007b ) + temp_007c
            temp_007 = 0.5temp_007d

            temp_008 = c_p2 + (temp_005 + temp_007)

            pn1[i] = newtonianCoeff[i]*(temp_006 + temp_008)

            temp_009 = X[i]*pn1[i]
            temp_010 = Y[i]*pn1[i]
            temp_011 = Z[i]*pn1[i]

            pn3[i] = 3.5*newtonian1b_Potential[i]

            temp_013a = pn2[i]*U[i]
            temp_013b = pn3[i]*pNxi[i]
            pntempX[i] = (temp_009 + (temp_013a+temp_013b))

            temp_014a = pn2[i]*V[i]
            temp_014b = pn3[i]*pNyi[i]
            pntempY[i] = (temp_010 + (temp_014a+temp_014b))

            temp_015a = pn2[i]*W[i]
            temp_015b = pn3[i]*pNzi[i]
            pntempZ[i] = (temp_011 + (temp_015a+temp_015b))
        end #for i
        for i in 1:N-1
            postNewtonX = postNewtonX + pntempX[i]
            postNewtonY = postNewtonY + pntempY[i]
            postNewtonZ = postNewtonZ + pntempZ[i]
        end
        postNewtonX = postNewtonX*c_m2
        postNewtonY = postNewtonY*c_m2
        postNewtonZ = postNewtonZ*c_m2
    end #for k in Base.OneTo(succ_approx_iter) # (post-Newtonian iterations)

    #computation of non-gravitational accelerations:
    hx = (Z[1]*V[1])-(Y[1]*W[1])
    hy = (X[1]*W[1])-(Z[1]*U[1])
    hz = (Y[1]*U[1])-(X[1]*V[1])
    r_hs = sqrt(r_p2[1])
    runitx = X[1]/r_hs
    runity = Y[2]/r_hs
    runitz = Z[3]/r_hs

    #cartesian components of transversal unit vector:
    tunitx0 = (hy*runitz) - (hz*runity)
    tunity0 = (hz*runitx) - (hx*runitz)
    tunitz0 = (hx*runity) - (hy*runitx)
    hmag = sqrt( ((tunitx0^2)+(tunity0^2))+(tunitz0^2) )
    tunitx = tunitx0/hmag
    tunity = tunity0/hmag
    tunitz = tunitz0/hmag

    # evaluate non-grav acceleration of NEA (Yarkovsky):
    g_r = r_hs^(-2.0)
    A2_t_g_r = q[7]*g_r

    NGAx = A2_t_g_r*tunitx
    NGAy = A2_t_g_r*tunity
    NGAz = A2_t_g_r*tunitz

    dq[4] = ( postNewtonX + accX ) + NGAx
    dq[5] = ( postNewtonY + accY ) + NGAy
    dq[6] = ( postNewtonZ + accZ ) + NGAz

    dq[7] = zero_q_1

    nothing
end
