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
@taylorize function RNp1BP_pN_A_J23E_J2S_ng_eph!(t, q, dq)
    local ss16asteph_t = ss16asteph(t)
    local acceph_t = acc_eph(t)
    local S = eltype(q[1])
    local N = length(μ) # number of bodies, including NEA
    local _1_to_N = Base.OneTo(N) # iterator over all bodies

    local succ_approx_iter = 1 # number of iterations of post-Newtonian subroutine
    local j2_body_index = [su, ea]#, mo] # indices of bodies with J2 flattening (note: Earth and Moon also have J3)

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

    #Newtonian accelerations
    newtonX = zero_q_1
    newtonY = zero_q_1
    newtonZ = zero_q_1

    xij = Array{Taylor1{S}}(undef, N, N)
    yij = Array{Taylor1{S}}(undef, N, N)
    zij = Array{Taylor1{S}}(undef, N, N)
    rij = Array{Taylor1{S}}(undef, N, N)
    r2ij = Array{Taylor1{S}}(undef, N, N)

    newtonianCoeff = Array{Taylor1{S}}(undef, N)

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

    pntempX = zero_q_1
    pntempY = zero_q_1
    pntempZ = zero_q_1

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
    local αs = deg2rad(268.13*one(t))
    local δs = deg2rad(63.87*one(t))
    local αm = moon_pole_ra(dsj2k)
    local δm = moon_pole_dec(dsj2k)
    local M_ = Array{Taylor1{S}}(undef, 3, 3, N)
    local M_[:,:,ea] = t2c_jpl_de430(dsj2k)
    local M_[:,:,su] = pole_rotation( αs, δs )
    local M_[:,:,mo] = pole_rotation( αm, δm )

    dq[1] = q[4]
    dq[2] = q[5]
    dq[3] = q[6]

    for i in 1:N
        newtonianNb_Potential[i] = zero_q_1
        for j in 1:N
            if j==i
            else
                xij[i,j] = zero_q_1
                yij[i,j] = zero_q_1
                zij[i,j] = zero_q_1
                r2ij[i,j] = zero_q_1
                rij[i,j] = zero_q_1
            end
        end
    end

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
        xi = ss16asteph_t[3i-2]
        yi = ss16asteph_t[3i-1]
        zi = ss16asteph_t[3i  ]
        ui = ss16asteph_t[3(N-1+i)-2]
        vi = ss16asteph_t[3(N-1+i)-1]
        wi = ss16asteph_t[3(N-1+i)  ]

        X[i] = xi-q[1]
        Y[i] = yi-q[2]
        Z[i] = zi-q[3]

        U[i] = ui-dq[1]
        V[i] = vi-dq[2]
        W[i] = wi-dq[3]

        _4U_m_3X[i] = (4dq[1])-(3ui)
        _4V_m_3Y[i] = (4dq[2])-(3vi)
        _4W_m_3Z[i] = (4dq[3])-(3wi)

        pn2x = X[i]*_4U_m_3X[i]
        pn2y = Y[i]*_4V_m_3Y[i]
        pn2z = Z[i]*_4W_m_3Z[i]

        UU[i] = ui*dq[1]
        VV[i] = vi*dq[2]
        WW[i] = wi*dq[3]

        vi_dot_vj[i] = ( UU[i]+VV[i] ) + WW[i]

        r_p2[i] = ( (X[i]^2)+(Y[i]^2) ) + (Z[i]^2)

        r_p1d2[i] = sqrt(r_p2[i])
        r_p3d2[i] = r_p2[i]^1.5
        r_p7d2[i] = r_p2[i]^3.5

        newtonianCoeff[i] =  μ[i]/r_p3d2[i]

        pn2[i] = newtonianCoeff[i]*(( pn2x+pn2y ) + pn2z)

        temp_001 = newtonX + (X[i]*newtonianCoeff[i])
        newtonX = temp_001
        temp_002 = newtonY + (Y[i]*newtonianCoeff[i])
        newtonY = temp_002
        temp_003 = newtonZ + (Z[i]*newtonianCoeff[i])
        newtonZ = temp_003

        newtonian1b_Potential[i] = μ[i]/r_p1d2[i]
        temp_004 = newtonianNb_Potential[N] + newtonian1b_Potential[i]
        newtonianNb_Potential[N] = temp_004

        for j in 1:N-1
            if j==i
            else
                xij[i,j] = xi-ss16asteph_t[3j-2]
                yij[i,j] = yi-ss16asteph_t[3j-1]
                zij[i,j] = zi-ss16asteph_t[3j  ]
                r2ij[i,j] = ((xij[i,j]^2) + (yij[i,j]^2)) + (zij[i,j]^2)
                rij[i,j] = sqrt( r2ij[i,j] )
                temp_004 = newtonianNb_Potential[i] + ( μ[j]/rij[i,j] )
                newtonianNb_Potential[i] = temp_004
            end
        end

        #J2 accelerations, if i-th body is flattened
        if UJ_interaction[N,i]
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
            temp_accX_i[i] = accX - (μ[i]*F_J2_x[i])
            accX = temp_accX_i[i]
            temp_accY_i[i] = accY - (μ[i]*F_J2_y[i])
            accY = temp_accY_i[i]
            temp_accZ_i[i] = accZ - (μ[i]*F_J2_z[i])
            accZ = temp_accZ_i[i]
        end
        v2[i] = ( (ui^2)+(vi^2) ) + (wi^2)
    end #for, i
    v2[N] = ( (q[4]^2)+(q[5]^2) ) + (q[6]^2)

    # postNewtonX = newtonX
    # postNewtonY = newtonY
    # postNewtonZ = newtonZ

    for k in Base.OneTo(succ_approx_iter)
        pntempX = zero_q_1
        pntempY = zero_q_1
        pntempZ = zero_q_1
        for i in 1:N-1
            #post-Newtonian corrections to gravitational acceleration
            #Moyer, 1971, page 7 eq. 35
            # post-Newtonian velocity of i-th body
            ui_ = ss16asteph_t[3(N-1+i)-2]
            vi_ = ss16asteph_t[3(N-1+i)-1]
            wi_ = ss16asteph_t[3(N-1+i)  ]
            # acceleration of i-th body
            pNxi = acceph_t[3i-2]
            pNyi = acceph_t[3i-1]
            pNzi = acceph_t[3i  ]

            temp_005a = newtonianNb_Potential[i] + (4newtonianNb_Potential[N])
            temp_005b = v2[N] + ( (2v2[i]) - (4vi_dot_vj[i]) )
            temp_005 = temp_005b - temp_005a

            temp_006a = X[i]*ui_
            temp_006b = Y[i]*vi_
            temp_006c = Z[i]*wi_
            temp_006d = ( temp_006a+temp_006b ) + temp_006c
            # the expression below inside the (...)^2 should have a minus sign in front of the numerator,
            # but upon squaring it is eliminated, so at the end of the day, it is irrelevant ;)
            temp_006e = (temp_006d^2)/r_p2[i]
            temp_006 = 1.5temp_006e

            temp_007a = X[i]*pNxi
            temp_007b = Y[i]*pNyi
            temp_007c = Z[i]*pNzi
            temp_007d = ( temp_007a+temp_007b ) + temp_007c
            temp_007 = 0.5temp_007d

            temp_008 = c_p2 + (temp_005 + temp_007)

            pn1[i] = newtonianCoeff[i]*(temp_006 + temp_008)

            temp_009 = X[i]*pn1[i]
            temp_010 = Y[i]*pn1[i]
            temp_011 = Z[i]*pn1[i]

            pn3[i] = 3.5*newtonian1b_Potential[i]

            temp_013a = pn2[i]*U[i]
            temp_013b = pn3[i]*pNxi
            temp_013 = pntempX + (temp_009 + (temp_013a+temp_013b))
            pntempX = temp_013

            temp_014a = pn2[i]*V[i]
            temp_014b = pn3[i]*pNyi
            temp_014 = pntempY + (temp_010 + (temp_014a+temp_014b))
            pntempY = temp_014

            temp_015a = pn2[i]*W[i]
            temp_015b = pn3[i]*pNzi
            temp_015 = pntempZ + (temp_011 + (temp_015a+temp_015b))
            pntempZ = temp_015
        end #for i
        postNewtonX = pntempX*c_m2
        postNewtonY = pntempY*c_m2
        postNewtonZ = pntempZ*c_m2
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

# # Nearth-Earth asteroid dynamical model (d=2.25)
# # Bodies considered in the model are: the Sun, the eight planets, the Moon and Ceres,
# # as well as the asteroid of interest as a test particle with null mass. Dynamical
# # effects considered are:
# # - post-Newtonian point-mass accelerations between all bodies,
# # - figure-effects (oblateness) of the Earth (J2, J3 and J4)
# # - J2 effect of the Sun
# # - Kinematic model for the precession and nutation of the Earth's pole (DE430/431 model)
# # - also, a model for non-gravitational accelerations acting upon the asteroid
# # is included (Yarkovsky effect) a_nongrav = A2*t_vec*(au/r)^d, where t_vec is the
# # unit heliocentric transverse vector, au is 1 astronomical unit, r is the
# # asteroid's heliocentric range, A2 is a coefficient (with units of au/day^2),
# # and d = 2.25
# @taylorize function RNp1BP_pN_A_J234E_J2S_ng_d225!(t, q, dq)
#     local S = eltype(q[1])
#     local N = Int((length(q)-1)/6) # number of bodies, including NEA
#     local _1_to_N = Base.OneTo(N) # iterator over all bodies

#     local succ_approx_iter = 1 # number of iterations of post-Newtonian subroutine
#     local su = 1 #Sun's index within `system`
#     local ea = 4 #Earth's index within `system`
#     local j2_body_index = [su, ea] # indices of bodies with J2 flattening (note: Earth also has J3 and J4)

#     # parameters related to speed of light, c
#     local c_p2 = 29979.063823897606 # c^2 = 29979.063823897606 au^2/d^2
#     local c_m2 = 3.3356611996764786e-5 # c^-2 = 3.3356611996764786e-5 d^2/au^2

#     local zero_q_1 = zero(q[1])

#     X = Array{Taylor1{S}}(undef, N, N)
#     Y = Array{Taylor1{S}}(undef, N, N)
#     Z = Array{Taylor1{S}}(undef, N, N)

#     r_p2 = Array{Taylor1{S}}(undef, N, N)
#     r_p3d2 = Array{Taylor1{S}}(undef, N, N)
#     r_p7d2 = Array{Taylor1{S}}(undef, N, N)

#     #Newtonian accelerations
#     newtonX = Array{Taylor1{S}}(undef, N)
#     newtonY = Array{Taylor1{S}}(undef, N)
#     newtonZ = Array{Taylor1{S}}(undef, N)

#     newtonianCoeff = Array{Taylor1{S}}(undef, N, N)

#     #post-Newtonian stuff
#     U = Array{Taylor1{S}}(undef, N, N)
#     V = Array{Taylor1{S}}(undef, N, N)
#     W = Array{Taylor1{S}}(undef, N, N)

#     _4U_m_3X = Array{Taylor1{S}}(undef, N, N)
#     _4V_m_3Y = Array{Taylor1{S}}(undef, N, N)
#     _4W_m_3Z = Array{Taylor1{S}}(undef, N, N)

#     UU = Array{Taylor1{S}}(undef, N, N)
#     VV = Array{Taylor1{S}}(undef, N, N)
#     WW = Array{Taylor1{S}}(undef, N, N)

#     r_p1d2 = Array{Taylor1{S}}(undef, N, N)

#     postNewtonX = Array{Taylor1{S}}(undef, N)
#     postNewtonY = Array{Taylor1{S}}(undef, N)
#     postNewtonZ = Array{Taylor1{S}}(undef, N)

#     newtonianNb_Potential = Array{Taylor1{S}}(undef, N)
#     newtonian1b_Potential = Array{Taylor1{S}}(undef, N, N)
#     newtonianCoeff = Array{Taylor1{S}}(undef, N, N)

#     pntempX = Array{Taylor1{S}}(undef, N)
#     pntempY = Array{Taylor1{S}}(undef, N)
#     pntempZ = Array{Taylor1{S}}(undef, N)

#     pn1 = Array{Taylor1{S}}(undef, N, N)
#     v2 = Array{Taylor1{S}}(undef, N)
#     vi_dot_vj = Array{Taylor1{S}}(undef, N, N)
#     pn2 = Array{Taylor1{S}}(undef, N, N)
#     pn3 = Array{Taylor1{S}}(undef, N, N)

#     # J2 acceleration auxiliaries
#     t11 = Array{Taylor1{S}}(undef, N, N)
#     t12 = Array{Taylor1{S}}(undef, N, N)
#     t13 = Array{Taylor1{S}}(undef, N, N)
#     t21 = Array{Taylor1{S}}(undef, N, N)
#     t22 = Array{Taylor1{S}}(undef, N, N)
#     t23 = Array{Taylor1{S}}(undef, N, N)
#     t31 = Array{Taylor1{S}}(undef, N, N)
#     t32 = Array{Taylor1{S}}(undef, N, N)
#     t33 = Array{Taylor1{S}}(undef, N, N)
#     new_x = Array{Taylor1{S}}(undef, N, N)
#     new_y = Array{Taylor1{S}}(undef, N, N)
#     new_z = Array{Taylor1{S}}(undef, N, N)
#     new_x2 = Array{Taylor1{S}}(undef, N, N)
#     new_y2 = Array{Taylor1{S}}(undef, N, N)
#     ρ_p2_ij = Array{Taylor1{S}}(undef, N, N)
#     z_p2_ij = Array{Taylor1{S}}(undef, N, N)
#     temp_ρ = Array{Taylor1{S}}(undef, N, N)
#     temp_z = Array{Taylor1{S}}(undef, N, N)
#     dum00 = Array{Taylor1{S}}(undef, N, N)
#     dum01 = Array{Taylor1{S}}(undef, N, N)
#     dum = Array{Taylor1{S}}(undef, N, N)
#     dum_ρ = Array{Taylor1{S}}(undef, N, N)
#     dum_z = Array{Taylor1{S}}(undef, N, N)
#     F_J2_bf_x = Array{Taylor1{S}}(undef, N, N)
#     F_J2_bf_y = Array{Taylor1{S}}(undef, N, N)
#     F_J2_bf_z = Array{Taylor1{S}}(undef, N, N)
#     s11 = Array{Taylor1{S}}(undef, N, N)
#     s12 = Array{Taylor1{S}}(undef, N, N)
#     s13 = Array{Taylor1{S}}(undef, N, N)
#     s21 = Array{Taylor1{S}}(undef, N, N)
#     s22 = Array{Taylor1{S}}(undef, N, N)
#     s23 = Array{Taylor1{S}}(undef, N, N)
#     s31 = Array{Taylor1{S}}(undef, N, N)
#     s32 = Array{Taylor1{S}}(undef, N, N)
#     s33 = Array{Taylor1{S}}(undef, N, N)
#     F_J2_x = Array{Taylor1{S}}(undef, N, N)
#     F_J2_y = Array{Taylor1{S}}(undef, N, N)
#     F_J2_z = Array{Taylor1{S}}(undef, N, N)
#     timp_004 = Array{Taylor1{S}}(undef, N, N)
#     timp_005 = Array{Taylor1{S}}(undef, N, N)
#     timp_006 = Array{Taylor1{S}}(undef, N, N)
#     tamp_004 = Array{Taylor1{S}}(undef, N, N)
#     tamp_005 = Array{Taylor1{S}}(undef, N, N)
#     tamp_006 = Array{Taylor1{S}}(undef, N, N)

#     # extended-body accelerations
#     accX = Array{Taylor1{S}}(undef, N)
#     accY = Array{Taylor1{S}}(undef, N)
#     accZ = Array{Taylor1{S}}(undef, N)

#     # rotations to and from Earth and Sun poles
#     local αs = deg2rad(268.13+0t)
#     local δs = deg2rad(63.87+0t)
#     local M_ = Array{Taylor1{S}}(undef, 3, 3, N)
#     local W_ = Array{Taylor1{S}}(undef, 3, 3, N)
#     local M_[:,:,ea] = earth_pole_rotation(t-J2000) # J2000.0 = 2.451545e6
#     local W_[:,:,ea] = inv(M_[:, :, ea])
#     local M_[:,:,su] = pole_rotation( αs, δs )
#     local W_[:,:,su] = inv(M_[:, :, su])

#     for j in _1_to_N
#         newtonX[j] = zero_q_1
#         newtonY[j] = zero_q_1
#         newtonZ[j] = zero_q_1

#         newtonianNb_Potential[j] = zero_q_1

#         accX[j] = zero_q_1
#         accY[j] = zero_q_1
#         accZ[j] = zero_q_1

#         dq[3j-2] = q[3(N+j)-2]
#         dq[3j-1] = q[3(N+j)-1]
#         dq[3j  ] = q[3(N+j)  ]
#     end

#     for j in j2_body_index
#         for i in _1_to_N
#             if i == j
#             else
#                 t11[i,j] = zero_q_1
#                 t12[i,j] = zero_q_1
#                 t13[i,j] = zero_q_1
#                 t21[i,j] = zero_q_1
#                 t22[i,j] = zero_q_1
#                 t23[i,j] = zero_q_1
#                 t31[i,j] = zero_q_1
#                 t32[i,j] = zero_q_1
#                 t33[i,j] = zero_q_1
#                 new_x[i,j] = zero_q_1
#                 new_x2[i,j] = zero_q_1
#                 new_y2[i,j] = zero_q_1
#                 ρ_p2_ij[i,j] = zero_q_1
#                 z_p2_ij[i,j] = zero_q_1
#                 temp_ρ[i,j] = zero_q_1
#                 temp_z[i,j] = zero_q_1
#                 dum00[i,j] = zero_q_1
#                 dum01[i,j] = zero_q_1
#                 dum[i,j] = zero_q_1
#                 dum_ρ[i,j] = zero_q_1
#                 dum_z[i,j] = zero_q_1
#                 F_J2_bf_x[i,j] = zero_q_1
#                 F_J2_bf_y[i,j] = zero_q_1
#                 F_J2_bf_z[i,j] = zero_q_1
#                 s11[i,j] = zero_q_1
#                 s12[i,j] = zero_q_1
#                 s13[i,j] = zero_q_1
#                 s21[i,j] = zero_q_1
#                 s22[i,j] = zero_q_1
#                 s23[i,j] = zero_q_1
#                 s31[i,j] = zero_q_1
#                 s32[i,j] = zero_q_1
#                 s33[i,j] = zero_q_1
#                 F_J2_x[i,j] = zero_q_1
#                 F_J2_y[i,j] = zero_q_1
#                 F_J2_z[i,j] = zero_q_1
#                 timp_004[i,j] = zero_q_1
#                 timp_005[i,j] = zero_q_1
#                 timp_006[i,j] = zero_q_1
#                 tamp_004[i,j] = zero_q_1
#                 tamp_005[i,j] = zero_q_1
#                 tamp_006[i,j] = zero_q_1
#             end #if i == j
#         end #for i in _1_to_N
#     end #for j in j2_body_index

#     #compute point-mass Newtonian accelerations, all bodies
#     for j in _1_to_N
#         for i in _1_to_N
#             # i == j && continue
#             if i == j
#             else
#                 X[i,j] = q[3i-2]-q[3j-2]
#                 Y[i,j] = q[3i-1]-q[3j-1]
#                 Z[i,j] = q[3i]-q[3j]

#                 U[i,j] = dq[3i-2]-dq[3j-2]
#                 V[i,j] = dq[3i-1]-dq[3j-1]
#                 W[i,j] = dq[3i  ]-dq[3j  ]

#                 _4U_m_3X[i,j] = (4dq[3j-2])-(3dq[3i-2])
#                 _4V_m_3Y[i,j] = (4dq[3j-1])-(3dq[3i-1])
#                 _4W_m_3Z[i,j] = (4dq[3j  ])-(3dq[3i  ])

#                 pn2x = X[i,j]*_4U_m_3X[i,j]
#                 pn2y = Y[i,j]*_4V_m_3Y[i,j]
#                 pn2z = Z[i,j]*_4W_m_3Z[i,j]

#                 pn2[i,j] = ( pn2x+pn2y ) + pn2z

#                 UU[i,j] = dq[3i-2]*dq[3j-2]
#                 VV[i,j] = dq[3i-1]*dq[3j-1]
#                 WW[i,j] = dq[3i  ]*dq[3j  ]

#                 vi_dot_vj[i,j] = ( UU[i,j]+VV[i,j] ) + WW[i,j]

#                 r_p2[i,j] = ( (X[i,j]^2)+(Y[i,j]^2) ) + (Z[i,j]^2)

#                 r_p1d2[i,j] = sqrt(r_p2[i,j])
#                 r_p3d2[i,j] = r_p2[i,j]^1.5
#                 r_p7d2[i,j] = r_p2[i,j]^3.5

#                 newtonianCoeff[i,j] =  μ[i]/r_p3d2[i,j]

#                 temp_001 = newtonX[j] + (X[i,j]*newtonianCoeff[i,j])
#                 newtonX[j] = temp_001
#                 temp_002 = newtonY[j] + (Y[i,j]*newtonianCoeff[i,j])
#                 newtonY[j] = temp_002
#                 temp_003 = newtonZ[j] + (Z[i,j]*newtonianCoeff[i,j])
#                 newtonZ[j] = temp_003

#                 newtonian1b_Potential[i,j] = μ[i]/r_p1d2[i, j]

#                 temp_004 = newtonianNb_Potential[j] + newtonian1b_Potential[i, j]
#                 newtonianNb_Potential[j] = temp_004
#             end #if i != j
#         end #for, i
#         v2[j] = ( (dq[3j-2]^2)+(dq[3j-1]^2) ) + (dq[3j]^2)
#     end #for, j

#     for j in _1_to_N
#         postNewtonX[j] = newtonX[j]
#         postNewtonY[j] = newtonY[j]
#         postNewtonZ[j] = newtonZ[j]
#     end

#     for k in Base.OneTo(succ_approx_iter)
#         for j in _1_to_N
#             pntempX[j] = zero_q_1
#             pntempY[j] = zero_q_1
#             pntempZ[j] = zero_q_1
#         end
#         for j in _1_to_N
#             for i in _1_to_N
#                 # i == j && continue
#                 if i == j
#                 else
#                     #post-Newtonian corrections to gravitational acceleration
#                     #Moyer, 1971, page 7 eq. 35
#                     temp_005a = newtonianNb_Potential[i]+(4newtonianNb_Potential[j])
#                     temp_005b = (2v2[i])-(4vi_dot_vj[i,j])
#                     temp_005c = v2[j]+temp_005b
#                     temp_005 = temp_005c-temp_005a
#                     temp_006a = X[i,j]*dq[3i-2]
#                     temp_006b = Y[i,j]*dq[3i-1]
#                     temp_006c = Z[i,j]*dq[3i]
#                     temp_006d = ( temp_006a+temp_006b ) + temp_006c
#                     # the expression below inside the (...)^2 should have a minus sign in front of the numerator,
#                     # but upon squaring it is eliminated, so at the end of the day, it is irrelevant ;)
#                     temp_006e = (temp_006d^2)/r_p2[i,j]
#                     temp_006 = temp_005-(1.5temp_006e)
#                     temp_007a = X[i,j]*postNewtonX[i]
#                     temp_007b = Y[i,j]*postNewtonY[i]
#                     temp_007c = Z[i,j]*postNewtonZ[i]
#                     temp_007d = ( temp_007a+temp_007b ) + temp_007c
#                     temp_007 = temp_006 + (0.5temp_007d)
#                     temp_008 = c_p2+temp_007
#                     pn1[i,j] = newtonianCoeff[i,j]*temp_008

#                     temp_009 = pntempX[j]+(X[i,j]*pn1[i,j])
#                     temp_010 = pntempY[j]+(Y[i,j]*pn1[i,j])
#                     temp_011 = pntempZ[j]+(Z[i,j]*pn1[i,j])

#                     pn3[i,j] = 3.5*newtonian1b_Potential[i,j]

#                     temp_013a = pn2[i,j]*(U[i,j]*newtonianCoeff[i,j])
#                     temp_013b = postNewtonX[i]*pn3[i,j]
#                     temp_013 = temp_009 + (temp_013a+temp_013b)
#                     pntempX[j] = temp_013

#                     temp_014a = pn2[i,j]*(V[i,j]*newtonianCoeff[i,j])
#                     temp_014b = postNewtonY[i]*pn3[i,j]
#                     temp_014 = temp_010 + (temp_014a+temp_014b)
#                     pntempY[j] = temp_014

#                     temp_015a = pn2[i,j]*(W[i,j]*newtonianCoeff[i,j])
#                     temp_015b = postNewtonZ[i]*pn3[i,j]
#                     temp_015 = temp_011 + (temp_015a+temp_015b)
#                     pntempZ[j] = temp_015
#                 end
#             end #for i
#         end #for j
#         for j in _1_to_N
#             postNewtonX[j] = pntempX[j]*c_m2
#             postNewtonY[j] = pntempY[j]*c_m2
#             postNewtonZ[j] = pntempZ[j]*c_m2
#         end
#     end #for k in Base.OneTo(succ_approx_iter) # (post-Newtonian iterations)

#     #J2 accelerations, all flattened bodies
#     for j in j2_body_index
#         for i in _1_to_N
#             if i == j
#             else
#                 # # rotate from inertial frame to extended-body frame
#                 t11[i,j] = X[i,j]*W_[1,1,j]
#                 t21[i,j] = X[i,j]*W_[2,1,j]
#                 t31[i,j] = X[i,j]*W_[3,1,j]
#                 t12[i,j] = Y[i,j]*W_[1,2,j]
#                 t22[i,j] = Y[i,j]*W_[2,2,j]
#                 t32[i,j] = Y[i,j]*W_[3,2,j]
#                 t13[i,j] = Z[i,j]*W_[1,3,j]
#                 t23[i,j] = Z[i,j]*W_[2,3,j]
#                 t33[i,j] = Z[i,j]*W_[3,3,j]
#                 new_x[i,j] = (t11[i,j]+t12[i,j])+t13[i,j]
#                 new_y[i,j] = (t21[i,j]+t22[i,j])+t23[i,j]
#                 new_z[i,j] = (t31[i,j]+t32[i,j])+t33[i,j]

#                 # # compute cartesian components of extended-body acceleration in body frame
#                 new_x2[i,j] = new_x[i,j]^2
#                 new_y2[i,j] = new_y[i,j]^2

#                 ρ_p2_ij[i,j] = new_x2[i,j]+new_y2[i,j]
#                 z_p2_ij[i,j] = new_z[i,j]^2

#                 temp_ρ[i,j] = ρ_p2_ij[i,j] - (4z_p2_ij[i,j])
#                 temp_z[i,j] = (3ρ_p2_ij[i,j]) - (2z_p2_ij[i,j])

#                 dum01[i,j] = Λ2[j]/r_p7d2[i,j]

#                 dum[i,j] = 1.5*dum01[i,j]
#                 dum_ρ[i,j] = dum[i,j]*temp_ρ[i,j]
#                 dum_z[i,j] = dum[i,j]*temp_z[i,j]

#                 F_J2_bf_x[i,j] = dum_ρ[i,j]*new_x[i,j]
#                 F_J2_bf_y[i,j] = dum_ρ[i,j]*new_y[i,j]
#                 F_J2_bf_z[i,j] = dum_z[i,j]*new_z[i,j]

#                 # # rotate components of force from body frame to inertial frame
#                 s11[i,j] = F_J2_bf_x[i,j]*M_[1,1,j]
#                 s21[i,j] = F_J2_bf_x[i,j]*M_[2,1,j]
#                 s31[i,j] = F_J2_bf_x[i,j]*M_[3,1,j]
#                 s12[i,j] = F_J2_bf_y[i,j]*M_[1,2,j]
#                 s22[i,j] = F_J2_bf_y[i,j]*M_[2,2,j]
#                 s32[i,j] = F_J2_bf_y[i,j]*M_[3,2,j]
#                 s13[i,j] = F_J2_bf_z[i,j]*M_[1,3,j]
#                 s23[i,j] = F_J2_bf_z[i,j]*M_[2,3,j]
#                 s33[i,j] = F_J2_bf_z[i,j]*M_[3,3,j]
#                 F_J2_x[i,j] = (s11[i,j]+s12[i,j])+s13[i,j]
#                 F_J2_y[i,j] = (s21[i,j]+s22[i,j])+s23[i,j]
#                 F_J2_z[i,j] = (s31[i,j]+s32[i,j])+s33[i,j]

#                 # # add result to total acceleration upon j-th body figure due to i-th point mass
#                 # @show "acc",j,"+μ",i,"Λ2",j
#                 timp_004[i,j] = accX[j] + (μ[i]*F_J2_x[i,j])
#                 accX[j] = timp_004[i,j]
#                 timp_005[i,j] = accY[j] + (μ[i]*F_J2_y[i,j])
#                 accY[j] = timp_005[i,j]
#                 timp_006[i,j] = accZ[j] + (μ[i]*F_J2_z[i,j])
#                 accZ[j] = timp_006[i,j]

#                 # # reaction force on i-th body
#                 # @show "acc",i,"-μ",j,"Λ2",j
#                 tamp_004[i,j] = accX[i] - (μ[j]*F_J2_x[i,j])
#                 accX[i] = tamp_004[i,j]
#                 tamp_005[i,j] = accY[i] - (μ[j]*F_J2_y[i,j])
#                 accY[i] = tamp_005[i,j]
#                 tamp_006[i,j] = accZ[i] - (μ[j]*F_J2_z[i,j])
#                 accZ[i] = tamp_006[i,j]
#             end #if i == j
#         end #for i in _1_to_N
#     end #for j in j2_body_index

#     #fill the equations of motion for everyone except test particle (Newtonian, post-Newtonian and extended body accelerations)
#     for i in Base.OneTo(N-1)
#         dq[3(N+i)-2] = postNewtonX[i]+accX[i]
#         dq[3(N+i)-1] = postNewtonY[i]+accY[i]
#         dq[3(N+i)  ] = postNewtonZ[i]+accZ[i]
#     end

#     #computation of non-gravitational accelerations:
#     hx = (Y[N,1]*(dq[3N  ]-dq[3]))-(Z[N,1]*(dq[3N-1]-dq[2]))
#     hy = (Z[N,1]*(dq[3N-2]-dq[1]))-(X[N,1]*(dq[3N  ]-dq[3]))
#     hz = (X[N,1]*(dq[3N-1]-dq[2]))-(Y[N,1]*(dq[3N-2]-dq[1]))
#     r_hs = sqrt(r_p2[N,1])
#     runitx = X[N,1]/r_hs
#     runity = Y[N,2]/r_hs
#     runitz = Z[N,3]/r_hs

#     #cartesian components of transversal unit vector:
#     tunitx0 = (hy*runitz)-(hz*runity)
#     tunity0 = (hz*runitx)-(hx*runitz)
#     tunitz0 = (hx*runity)-(hy*runitx)
#     hmag = sqrt( ((tunitx0^2)+(tunity0^2))+(tunitz0^2) )
#     tunitx = tunitx0/hmag
#     tunity = tunity0/hmag
#     tunitz = tunitz0/hmag

#     # evaluate non-grav acceleration of NEA (Yarkovsky):
#     g_r = r_hs^(-2.25)
#     A2_t_g_r = q[6N+1]*g_r

#     NGAx = A2_t_g_r*tunitx
#     NGAy = A2_t_g_r*tunity
#     NGAz = A2_t_g_r*tunitz

#     dq[6N-2] = (postNewtonX[N]+accX[N])+NGAx
#     dq[6N-1] = (postNewtonY[N]+accY[N])+NGAy
#     dq[6N  ] = (postNewtonZ[N]+accZ[N])+NGAz

#     dq[6N+1] = zero_q_1

#     nothing
# end

# # Nearth-Earth asteroid dynamical model (d=2.25 and direct solar radiation pressure)
# # Bodies considered in the model are: the Sun, the eight planets, the Moon and Ceres,
# # as well as the asteroid of interest as a test particle with null mass. Dynamical
# # effects considered are:
# # - post-Newtonian point-mass accelerations between all bodies,
# # - figure-effects (oblateness) of the Earth (J2, J3 and J4)
# # - J2 effect of the Sun
# # - Kinematic model for the precession and nutation of the Earth's pole (DE430/431 model)
# # - also, a model for non-gravitational accelerations acting upon the asteroid
# # is included (Yarkovsky effect) a_nongrav = A2*t_vec*(au/r)^d, where t_vec is the
# # unit heliocentric transverse vector, au is 1 astronomical unit, r is the
# # asteroid's heliocentric range, A2 is a coefficient (with units of au/day^2),
# # and d = 2.25
# @taylorize function RNp1BP_pN_A_J234E_J2S_ng_d225_srp!(t, q, dq)
#     local S = eltype(q[1])
#     local N = Int((length(q)-1)/6) # number of bodies, including NEA
#     local _1_to_N = Base.OneTo(N) # iterator over all bodies

#     local succ_approx_iter = 1 # number of iterations of post-Newtonian subroutine
#     local su = 1 #Sun's index within `system`
#     local ea = 4 #Earth's index within `system`
#     local j2_body_index = [su, ea] # indices of bodies with J2 flattening (note: Earth also has J3 and J4)

#     # parameters related to speed of light, c
#     local c_p2 = 29979.063823897606 # c^2 = 29979.063823897606 au^2/d^2
#     local c_m2 = 3.3356611996764786e-5 # c^-2 = 3.3356611996764786e-5 d^2/au^2

#     local amrat = 3.07e-6 # m^2/kg, area-to-mass ratio (Bennu)

#     local zero_q_1 = zero(q[1])

#     X = Array{Taylor1{S}}(undef, N, N)
#     Y = Array{Taylor1{S}}(undef, N, N)
#     Z = Array{Taylor1{S}}(undef, N, N)

#     r_p2 = Array{Taylor1{S}}(undef, N, N)
#     r_p3d2 = Array{Taylor1{S}}(undef, N, N)
#     r_p7d2 = Array{Taylor1{S}}(undef, N, N)

#     #Newtonian accelerations
#     newtonX = Array{Taylor1{S}}(undef, N)
#     newtonY = Array{Taylor1{S}}(undef, N)
#     newtonZ = Array{Taylor1{S}}(undef, N)

#     newtonianCoeff = Array{Taylor1{S}}(undef, N, N)

#     #post-Newtonian stuff
#     U = Array{Taylor1{S}}(undef, N, N)
#     V = Array{Taylor1{S}}(undef, N, N)
#     W = Array{Taylor1{S}}(undef, N, N)

#     _4U_m_3X = Array{Taylor1{S}}(undef, N, N)
#     _4V_m_3Y = Array{Taylor1{S}}(undef, N, N)
#     _4W_m_3Z = Array{Taylor1{S}}(undef, N, N)

#     UU = Array{Taylor1{S}}(undef, N, N)
#     VV = Array{Taylor1{S}}(undef, N, N)
#     WW = Array{Taylor1{S}}(undef, N, N)

#     r_p1d2 = Array{Taylor1{S}}(undef, N, N)

#     postNewtonX = Array{Taylor1{S}}(undef, N)
#     postNewtonY = Array{Taylor1{S}}(undef, N)
#     postNewtonZ = Array{Taylor1{S}}(undef, N)

#     newtonianNb_Potential = Array{Taylor1{S}}(undef, N)
#     newtonian1b_Potential = Array{Taylor1{S}}(undef, N, N)
#     newtonianCoeff = Array{Taylor1{S}}(undef, N, N)

#     pntempX = Array{Taylor1{S}}(undef, N)
#     pntempY = Array{Taylor1{S}}(undef, N)
#     pntempZ = Array{Taylor1{S}}(undef, N)

#     pn1 = Array{Taylor1{S}}(undef, N, N)
#     v2 = Array{Taylor1{S}}(undef, N)
#     vi_dot_vj = Array{Taylor1{S}}(undef, N, N)
#     pn2 = Array{Taylor1{S}}(undef, N, N)
#     pn3 = Array{Taylor1{S}}(undef, N, N)

#     # J2 acceleration auxiliaries
#     t11 = Array{Taylor1{S}}(undef, N, N)
#     t12 = Array{Taylor1{S}}(undef, N, N)
#     t13 = Array{Taylor1{S}}(undef, N, N)
#     t21 = Array{Taylor1{S}}(undef, N, N)
#     t22 = Array{Taylor1{S}}(undef, N, N)
#     t23 = Array{Taylor1{S}}(undef, N, N)
#     t31 = Array{Taylor1{S}}(undef, N, N)
#     t32 = Array{Taylor1{S}}(undef, N, N)
#     t33 = Array{Taylor1{S}}(undef, N, N)
#     new_x = Array{Taylor1{S}}(undef, N, N)
#     new_y = Array{Taylor1{S}}(undef, N, N)
#     new_z = Array{Taylor1{S}}(undef, N, N)
#     new_x2 = Array{Taylor1{S}}(undef, N, N)
#     new_y2 = Array{Taylor1{S}}(undef, N, N)
#     ρ_p2_ij = Array{Taylor1{S}}(undef, N, N)
#     z_p2_ij = Array{Taylor1{S}}(undef, N, N)
#     temp_ρ = Array{Taylor1{S}}(undef, N, N)
#     temp_z = Array{Taylor1{S}}(undef, N, N)
#     dum00 = Array{Taylor1{S}}(undef, N, N)
#     dum01 = Array{Taylor1{S}}(undef, N, N)
#     dum = Array{Taylor1{S}}(undef, N, N)
#     dum_ρ = Array{Taylor1{S}}(undef, N, N)
#     dum_z = Array{Taylor1{S}}(undef, N, N)
#     F_J2_bf_x = Array{Taylor1{S}}(undef, N, N)
#     F_J2_bf_y = Array{Taylor1{S}}(undef, N, N)
#     F_J2_bf_z = Array{Taylor1{S}}(undef, N, N)
#     s11 = Array{Taylor1{S}}(undef, N, N)
#     s12 = Array{Taylor1{S}}(undef, N, N)
#     s13 = Array{Taylor1{S}}(undef, N, N)
#     s21 = Array{Taylor1{S}}(undef, N, N)
#     s22 = Array{Taylor1{S}}(undef, N, N)
#     s23 = Array{Taylor1{S}}(undef, N, N)
#     s31 = Array{Taylor1{S}}(undef, N, N)
#     s32 = Array{Taylor1{S}}(undef, N, N)
#     s33 = Array{Taylor1{S}}(undef, N, N)
#     F_J2_x = Array{Taylor1{S}}(undef, N, N)
#     F_J2_y = Array{Taylor1{S}}(undef, N, N)
#     F_J2_z = Array{Taylor1{S}}(undef, N, N)
#     timp_004 = Array{Taylor1{S}}(undef, N, N)
#     timp_005 = Array{Taylor1{S}}(undef, N, N)
#     timp_006 = Array{Taylor1{S}}(undef, N, N)
#     tamp_004 = Array{Taylor1{S}}(undef, N, N)
#     tamp_005 = Array{Taylor1{S}}(undef, N, N)
#     tamp_006 = Array{Taylor1{S}}(undef, N, N)

#     # extended-body accelerations
#     accX = Array{Taylor1{S}}(undef, N)
#     accY = Array{Taylor1{S}}(undef, N)
#     accZ = Array{Taylor1{S}}(undef, N)

#     # rotations to and from Earth and Sun poles
#     local αs = deg2rad(268.13+0t)
#     local δs = deg2rad(63.87+0t)
#     local M_ = Array{Taylor1{S}}(undef, 3, 3, N)
#     local W_ = Array{Taylor1{S}}(undef, 3, 3, N)
#     local M_[:,:,ea] = earth_pole_rotation(t-J2000) # J2000.0 = 2.451545e6
#     local W_[:,:,ea] = inv(M_[:, :, ea])
#     local M_[:,:,su] = pole_rotation( αs, δs )
#     local W_[:,:,su] = inv(M_[:, :, su])

#     for j in _1_to_N
#         newtonX[j] = zero_q_1
#         newtonY[j] = zero_q_1
#         newtonZ[j] = zero_q_1

#         newtonianNb_Potential[j] = zero_q_1

#         accX[j] = zero_q_1
#         accY[j] = zero_q_1
#         accZ[j] = zero_q_1

#         dq[3j-2] = q[3(N+j)-2]
#         dq[3j-1] = q[3(N+j)-1]
#         dq[3j  ] = q[3(N+j)  ]
#     end

#     for j in j2_body_index
#         for i in _1_to_N
#             if i == j
#             else
#                 t11[i,j] = zero_q_1
#                 t12[i,j] = zero_q_1
#                 t13[i,j] = zero_q_1
#                 t21[i,j] = zero_q_1
#                 t22[i,j] = zero_q_1
#                 t23[i,j] = zero_q_1
#                 t31[i,j] = zero_q_1
#                 t32[i,j] = zero_q_1
#                 t33[i,j] = zero_q_1
#                 new_x[i,j] = zero_q_1
#                 new_x2[i,j] = zero_q_1
#                 new_y2[i,j] = zero_q_1
#                 ρ_p2_ij[i,j] = zero_q_1
#                 z_p2_ij[i,j] = zero_q_1
#                 temp_ρ[i,j] = zero_q_1
#                 temp_z[i,j] = zero_q_1
#                 dum00[i,j] = zero_q_1
#                 dum01[i,j] = zero_q_1
#                 dum[i,j] = zero_q_1
#                 dum_ρ[i,j] = zero_q_1
#                 dum_z[i,j] = zero_q_1
#                 F_J2_bf_x[i,j] = zero_q_1
#                 F_J2_bf_y[i,j] = zero_q_1
#                 F_J2_bf_z[i,j] = zero_q_1
#                 s11[i,j] = zero_q_1
#                 s12[i,j] = zero_q_1
#                 s13[i,j] = zero_q_1
#                 s21[i,j] = zero_q_1
#                 s22[i,j] = zero_q_1
#                 s23[i,j] = zero_q_1
#                 s31[i,j] = zero_q_1
#                 s32[i,j] = zero_q_1
#                 s33[i,j] = zero_q_1
#                 F_J2_x[i,j] = zero_q_1
#                 F_J2_y[i,j] = zero_q_1
#                 F_J2_z[i,j] = zero_q_1
#                 timp_004[i,j] = zero_q_1
#                 timp_005[i,j] = zero_q_1
#                 timp_006[i,j] = zero_q_1
#                 tamp_004[i,j] = zero_q_1
#                 tamp_005[i,j] = zero_q_1
#                 tamp_006[i,j] = zero_q_1
#             end #if i == j
#         end #for i in _1_to_N
#     end #for j in j2_body_index

#     #compute point-mass Newtonian accelerations, all bodies
#     for j in _1_to_N
#         for i in _1_to_N
#             # i == j && continue
#             if i == j
#             else
#                 X[i,j] = q[3i-2]-q[3j-2]
#                 Y[i,j] = q[3i-1]-q[3j-1]
#                 Z[i,j] = q[3i]-q[3j]

#                 U[i,j] = dq[3i-2]-dq[3j-2]
#                 V[i,j] = dq[3i-1]-dq[3j-1]
#                 W[i,j] = dq[3i  ]-dq[3j  ]

#                 _4U_m_3X[i,j] = (4dq[3j-2])-(3dq[3i-2])
#                 _4V_m_3Y[i,j] = (4dq[3j-1])-(3dq[3i-1])
#                 _4W_m_3Z[i,j] = (4dq[3j  ])-(3dq[3i  ])

#                 pn2x = X[i,j]*_4U_m_3X[i,j]
#                 pn2y = Y[i,j]*_4V_m_3Y[i,j]
#                 pn2z = Z[i,j]*_4W_m_3Z[i,j]

#                 pn2[i,j] = ( pn2x+pn2y ) + pn2z

#                 UU[i,j] = dq[3i-2]*dq[3j-2]
#                 VV[i,j] = dq[3i-1]*dq[3j-1]
#                 WW[i,j] = dq[3i  ]*dq[3j  ]

#                 vi_dot_vj[i,j] = ( UU[i,j]+VV[i,j] ) + WW[i,j]

#                 r_p2[i,j] = ( (X[i,j]^2)+(Y[i,j]^2) ) + (Z[i,j]^2)

#                 r_p1d2[i,j] = sqrt(r_p2[i,j])
#                 r_p3d2[i,j] = r_p2[i,j]^1.5
#                 r_p7d2[i,j] = r_p2[i,j]^3.5

#                 newtonianCoeff[i,j] =  μ[i]/r_p3d2[i,j]

#                 temp_001 = newtonX[j] + (X[i,j]*newtonianCoeff[i,j])
#                 newtonX[j] = temp_001
#                 temp_002 = newtonY[j] + (Y[i,j]*newtonianCoeff[i,j])
#                 newtonY[j] = temp_002
#                 temp_003 = newtonZ[j] + (Z[i,j]*newtonianCoeff[i,j])
#                 newtonZ[j] = temp_003

#                 newtonian1b_Potential[i,j] = μ[i]/r_p1d2[i, j]

#                 temp_004 = newtonianNb_Potential[j] + newtonian1b_Potential[i, j]
#                 newtonianNb_Potential[j] = temp_004
#             end #if i != j
#         end #for, i
#         v2[j] = ( (dq[3j-2]^2)+(dq[3j-1]^2) ) + (dq[3j]^2)
#     end #for, j

#     for j in _1_to_N
#         postNewtonX[j] = newtonX[j]
#         postNewtonY[j] = newtonY[j]
#         postNewtonZ[j] = newtonZ[j]
#     end

#     for k in Base.OneTo(succ_approx_iter)
#         for j in _1_to_N
#             pntempX[j] = zero_q_1
#             pntempY[j] = zero_q_1
#             pntempZ[j] = zero_q_1
#         end
#         for j in _1_to_N
#             for i in _1_to_N
#                 # i == j && continue
#                 if i == j
#                 else
#                     #post-Newtonian corrections to gravitational acceleration
#                     #Moyer, 1971, page 7 eq. 35
#                     temp_005a = newtonianNb_Potential[i]+(4newtonianNb_Potential[j])
#                     temp_005b = (2v2[i])-(4vi_dot_vj[i,j])
#                     temp_005c = v2[j]+temp_005b
#                     temp_005 = temp_005c-temp_005a
#                     temp_006a = X[i,j]*dq[3i-2]
#                     temp_006b = Y[i,j]*dq[3i-1]
#                     temp_006c = Z[i,j]*dq[3i]
#                     temp_006d = ( temp_006a+temp_006b ) + temp_006c
#                     # the expression below inside the (...)^2 should have a minus sign in front of the numerator,
#                     # but upon squaring it is eliminated, so at the end of the day, it is irrelevant ;)
#                     temp_006e = (temp_006d^2)/r_p2[i,j]
#                     temp_006 = temp_005-(1.5temp_006e)
#                     temp_007a = X[i,j]*postNewtonX[i]
#                     temp_007b = Y[i,j]*postNewtonY[i]
#                     temp_007c = Z[i,j]*postNewtonZ[i]
#                     temp_007d = ( temp_007a+temp_007b ) + temp_007c
#                     temp_007 = temp_006 + (0.5temp_007d)
#                     temp_008 = c_p2+temp_007
#                     pn1[i,j] = newtonianCoeff[i,j]*temp_008

#                     temp_009 = pntempX[j]+(X[i,j]*pn1[i,j])
#                     temp_010 = pntempY[j]+(Y[i,j]*pn1[i,j])
#                     temp_011 = pntempZ[j]+(Z[i,j]*pn1[i,j])

#                     pn3[i,j] = 3.5*newtonian1b_Potential[i,j]

#                     temp_013a = pn2[i,j]*(U[i,j]*newtonianCoeff[i,j])
#                     temp_013b = postNewtonX[i]*pn3[i,j]
#                     temp_013 = temp_009 + (temp_013a+temp_013b)
#                     pntempX[j] = temp_013

#                     temp_014a = pn2[i,j]*(V[i,j]*newtonianCoeff[i,j])
#                     temp_014b = postNewtonY[i]*pn3[i,j]
#                     temp_014 = temp_010 + (temp_014a+temp_014b)
#                     pntempY[j] = temp_014

#                     temp_015a = pn2[i,j]*(W[i,j]*newtonianCoeff[i,j])
#                     temp_015b = postNewtonZ[i]*pn3[i,j]
#                     temp_015 = temp_011 + (temp_015a+temp_015b)
#                     pntempZ[j] = temp_015
#                 end
#             end #for i
#         end #for j
#         for j in _1_to_N
#             postNewtonX[j] = pntempX[j]*c_m2
#             postNewtonY[j] = pntempY[j]*c_m2
#             postNewtonZ[j] = pntempZ[j]*c_m2
#         end
#     end #for k in Base.OneTo(succ_approx_iter) # (post-Newtonian iterations)

#     #J2 accelerations, all flattened bodies
#     for j in j2_body_index
#         for i in _1_to_N
#             if i == j
#             else
#                 # # rotate from inertial frame to extended-body frame
#                 t11[i,j] = X[i,j]*W_[1,1,j]
#                 t21[i,j] = X[i,j]*W_[2,1,j]
#                 t31[i,j] = X[i,j]*W_[3,1,j]
#                 t12[i,j] = Y[i,j]*W_[1,2,j]
#                 t22[i,j] = Y[i,j]*W_[2,2,j]
#                 t32[i,j] = Y[i,j]*W_[3,2,j]
#                 t13[i,j] = Z[i,j]*W_[1,3,j]
#                 t23[i,j] = Z[i,j]*W_[2,3,j]
#                 t33[i,j] = Z[i,j]*W_[3,3,j]
#                 new_x[i,j] = (t11[i,j]+t12[i,j])+t13[i,j]
#                 new_y[i,j] = (t21[i,j]+t22[i,j])+t23[i,j]
#                 new_z[i,j] = (t31[i,j]+t32[i,j])+t33[i,j]

#                 # # compute cartesian components of extended-body acceleration in body frame
#                 new_x2[i,j] = new_x[i,j]^2
#                 new_y2[i,j] = new_y[i,j]^2

#                 ρ_p2_ij[i,j] = new_x2[i,j]+new_y2[i,j]
#                 z_p2_ij[i,j] = new_z[i,j]^2

#                 temp_ρ[i,j] = ρ_p2_ij[i,j] - (4z_p2_ij[i,j])
#                 temp_z[i,j] = (3ρ_p2_ij[i,j]) - (2z_p2_ij[i,j])

#                 dum01[i,j] = Λ2[j]/r_p7d2[i,j]

#                 dum[i,j] = 1.5*dum01[i,j]
#                 dum_ρ[i,j] = dum[i,j]*temp_ρ[i,j]
#                 dum_z[i,j] = dum[i,j]*temp_z[i,j]

#                 F_J2_bf_x[i,j] = dum_ρ[i,j]*new_x[i,j]
#                 F_J2_bf_y[i,j] = dum_ρ[i,j]*new_y[i,j]
#                 F_J2_bf_z[i,j] = dum_z[i,j]*new_z[i,j]

#                 # # rotate components of force from body frame to inertial frame
#                 s11[i,j] = F_J2_bf_x[i,j]*M_[1,1,j]
#                 s21[i,j] = F_J2_bf_x[i,j]*M_[2,1,j]
#                 s31[i,j] = F_J2_bf_x[i,j]*M_[3,1,j]
#                 s12[i,j] = F_J2_bf_y[i,j]*M_[1,2,j]
#                 s22[i,j] = F_J2_bf_y[i,j]*M_[2,2,j]
#                 s32[i,j] = F_J2_bf_y[i,j]*M_[3,2,j]
#                 s13[i,j] = F_J2_bf_z[i,j]*M_[1,3,j]
#                 s23[i,j] = F_J2_bf_z[i,j]*M_[2,3,j]
#                 s33[i,j] = F_J2_bf_z[i,j]*M_[3,3,j]
#                 F_J2_x[i,j] = (s11[i,j]+s12[i,j])+s13[i,j]
#                 F_J2_y[i,j] = (s21[i,j]+s22[i,j])+s23[i,j]
#                 F_J2_z[i,j] = (s31[i,j]+s32[i,j])+s33[i,j]

#                 # # add result to total acceleration upon j-th body figure due to i-th point mass
#                 # @show "acc",j,"+μ",i,"Λ2",j
#                 timp_004[i,j] = accX[j] + (μ[i]*F_J2_x[i,j])
#                 accX[j] = timp_004[i,j]
#                 timp_005[i,j] = accY[j] + (μ[i]*F_J2_y[i,j])
#                 accY[j] = timp_005[i,j]
#                 timp_006[i,j] = accZ[j] + (μ[i]*F_J2_z[i,j])
#                 accZ[j] = timp_006[i,j]

#                 # # reaction force on i-th body
#                 # @show "acc",i,"-μ",j,"Λ2",j
#                 tamp_004[i,j] = accX[i] - (μ[j]*F_J2_x[i,j])
#                 accX[i] = tamp_004[i,j]
#                 tamp_005[i,j] = accY[i] - (μ[j]*F_J2_y[i,j])
#                 accY[i] = tamp_005[i,j]
#                 tamp_006[i,j] = accZ[i] - (μ[j]*F_J2_z[i,j])
#                 accZ[i] = tamp_006[i,j]
#             end #if i == j
#         end #for i in _1_to_N
#     end #for j in j2_body_index

#     #fill the equations of motion for everyone except test particle (Newtonian, post-Newtonian and extended body accelerations)
#     for i in Base.OneTo(N-1)
#         dq[3(N+i)-2] = postNewtonX[i]+accX[i]
#         dq[3(N+i)-1] = postNewtonY[i]+accY[i]
#         dq[3(N+i)  ] = postNewtonZ[i]+accZ[i]
#     end

#     #computation of non-gravitational accelerations:
#     hx = (Y[N,1]*(dq[3N  ]-dq[3]))-(Z[N,1]*(dq[3N-1]-dq[2]))
#     hy = (Z[N,1]*(dq[3N-2]-dq[1]))-(X[N,1]*(dq[3N  ]-dq[3]))
#     hz = (X[N,1]*(dq[3N-1]-dq[2]))-(Y[N,1]*(dq[3N-2]-dq[1]))
#     r_hs = sqrt(r_p2[N,1])
#     runitx = X[N,1]/r_hs
#     runity = Y[N,2]/r_hs
#     runitz = Z[N,3]/r_hs

#     #cartesian components of transversal unit vector:
#     tunitx0 = (hy*runitz)-(hz*runity)
#     tunity0 = (hz*runitx)-(hx*runitz)
#     tunitz0 = (hx*runity)-(hy*runitx)
#     hmag = sqrt( ((tunitx0^2)+(tunity0^2))+(tunitz0^2) )
#     tunitx = tunitx0/hmag
#     tunity = tunity0/hmag
#     tunitz = tunitz0/hmag

#     # evaluate non-grav acceleration of asteroid (Yarkovsky):
#     g_r = r_hs^(-2.25)
#     A2_t_g_r = q[6N+1]*g_r

#     NGAx = A2_t_g_r*tunitx
#     NGAy = A2_t_g_r*tunity
#     NGAz = A2_t_g_r*tunitz

#     # evaluate solar radation pressure
#     SRP_coeff = ((S0_sun*amrat)*m2_s3_to_au2_day3)/c_au_per_day
#     SRP_radial = (R_sun/r_hs)^2
#     SRP_fun = SRP_coeff*SRP_radial

#     SRPx = SRP_fun*runitx
#     SRPy = SRP_fun*runity
#     SRPz = SRP_fun*runitz

#     #add non-gravitational accelerations to asteroid acceleration
#     dq[6N-2] = (postNewtonX[N]+accX[N])+(NGAx+SRPx)
#     dq[6N-1] = (postNewtonY[N]+accY[N])+(NGAy+SRPy)
#     dq[6N  ] = (postNewtonZ[N]+accZ[N])+(NGAz+SRPz)

#     dq[6N+1] = zero_q_1

#     nothing
# end
