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
function RNp1BP_pN_A_J23E_J2S_ng_eph!(dq, q, params, t)
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

function TaylorIntegration.jetcoeffs!(::Val{RNp1BP_pN_A_J23E_J2S_ng_eph!}, t::Taylor1{_T}, q::AbstractVector{Taylor1{_S}}, dq::AbstractVector{Taylor1{_S}}, params) where {_T <: Real, _S <: Number}
    order = t.order
    local ss16asteph_t = (params[1])(t) * one(q[1])
    local acceph_t = (params[2])(t) * one(q[1])
    local newtonianNb_Potential_t = (params[3])(t) * one(q[1])
    local S = eltype(q[1])
    local N = length(μ)
    local _1_to_N = Base.OneTo(N)
    local succ_approx_iter = 1
    local c_p2 = 29979.063823897606
    local c_m2 = 3.3356611996764786e-5
    local zero_q_1 = zero(q[1])
    X = Array{Taylor1{S}}(undef, N)
    Y = Array{Taylor1{S}}(undef, N)
    Z = Array{Taylor1{S}}(undef, N)
    r_p2 = Array{Taylor1{S}}(undef, N)
    r_p1d2 = Array{Taylor1{S}}(undef, N)
    r_p3d2 = Array{Taylor1{S}}(undef, N)
    r_p7d2 = Array{Taylor1{S}}(undef, N)
    newtonianCoeff = Array{Taylor1{S}}(undef, N)
    xi = Array{Taylor1{S}}(undef, N - 1)
    yi = Array{Taylor1{S}}(undef, N - 1)
    zi = Array{Taylor1{S}}(undef, N - 1)
    ui = Array{Taylor1{S}}(undef, N - 1)
    vi = Array{Taylor1{S}}(undef, N - 1)
    wi = Array{Taylor1{S}}(undef, N - 1)
    _3ui = Array{Taylor1{S}}(undef, N - 1)
    _3vi = Array{Taylor1{S}}(undef, N - 1)
    _3wi = Array{Taylor1{S}}(undef, N - 1)
    pNxi = Array{Taylor1{S}}(undef, N - 1)
    pNyi = Array{Taylor1{S}}(undef, N - 1)
    pNzi = Array{Taylor1{S}}(undef, N - 1)
    U = Array{Taylor1{S}}(undef, N)
    V = Array{Taylor1{S}}(undef, N)
    W = Array{Taylor1{S}}(undef, N)
    _4U_m_3X = Array{Taylor1{S}}(undef, N)
    _4V_m_3Y = Array{Taylor1{S}}(undef, N)
    _4W_m_3Z = Array{Taylor1{S}}(undef, N)
    UU = Array{Taylor1{S}}(undef, N)
    VV = Array{Taylor1{S}}(undef, N)
    WW = Array{Taylor1{S}}(undef, N)
    postNewtonX = Taylor1(identity(constant_term(zero_q_1)), order)
    postNewtonY = Taylor1(identity(constant_term(zero_q_1)), order)
    postNewtonZ = Taylor1(identity(constant_term(zero_q_1)), order)
    newtonianNb_Potential = Array{Taylor1{S}}(undef, N)
    newtonian1b_Potential = Array{Taylor1{S}}(undef, N)
    pntempX = Array{Taylor1{S}}(undef, N - 1)
    pntempY = Array{Taylor1{S}}(undef, N - 1)
    pntempZ = Array{Taylor1{S}}(undef, N - 1)
    pn1 = Array{Taylor1{S}}(undef, N)
    v2 = Array{Taylor1{S}}(undef, N)
    vi_dot_vj = Array{Taylor1{S}}(undef, N)
    pn2 = Array{Taylor1{S}}(undef, N)
    pn3 = Array{Taylor1{S}}(undef, N)
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
    accX = Taylor1(identity(constant_term(zero_q_1)), order)
    accY = Taylor1(identity(constant_term(zero_q_1)), order)
    accZ = Taylor1(identity(constant_term(zero_q_1)), order)
    local dsj2k = t - 2.451545e6
    local αs = deg2rad(α_p_sun * one(t))
    local δs = deg2rad(δ_p_sun * one(t))
    local αm = moon_pole_ra(dsj2k)
    local δm = moon_pole_dec(dsj2k)
    local M_ = Array{Taylor1{S}}(undef, 3, 3, N)
    local M_[:, :, ea] = t2c_jpl_de430(dsj2k)
    local M_[:, :, su] = pole_rotation(αs, δs)
    dq[1] = Taylor1(identity(constant_term(q[4])), order)
    dq[2] = Taylor1(identity(constant_term(q[5])), order)
    dq[3] = Taylor1(identity(constant_term(q[6])), order)
    newtonianNb_Potential[N] = Taylor1(identity(constant_term(zero_q_1)), order)
    for j = j2_body_index
        t31[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        t32[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        t33[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        F_J2_x[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        F_J2_y[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        F_J2_z[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        F_J2_x1[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        F_J2_y1[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        F_J2_z1[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        F_J2_x2[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        F_J2_y2[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        F_J2_z2[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        temp_accX_i[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        temp_accY_i[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        temp_accZ_i[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        sin_ϕ[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        sin2_ϕ[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        sin3_ϕ[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        sin4_ϕ[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        ϕ[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        cos_ϕ[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        P_2_sin_ϕ[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        ∂P_2_sin_ϕ[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        P_3_sin_ϕ[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        ∂P_3_sin_ϕ[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        m_c_ϕ_∂P_2[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        m_c_ϕ_∂P_3[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        Λ2j_div_r4[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        Λ3j_div_r5[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        F_J_ξ[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        F_J_η[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        F_J_ζ[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        F_J2_ξ[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        F_J2_η[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        F_J2_ζ[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        F_J3_ξ[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        F_J3_η[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        F_J3_ζ[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        ξx[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        ξy[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        ξz[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        ηx[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        ηy[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        ηz[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        ηx1[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        ηy1[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        ηz1[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        ηx2[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        ηy2[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        ηz2[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        ζx[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        ζy[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        ζz[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        ζx1[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        ζy1[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        ζz1[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        ζx2[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        ζy2[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        ζz2[j] = Taylor1(identity(constant_term(zero_q_1)), order)
    end
    tmp784 = Array{Taylor1{_S}}(undef, size(dq))
    tmp784 .= Taylor1(zero(_S), order)
    tmp787 = Array{Taylor1{_S}}(undef, size(dq))
    tmp787 .= Taylor1(zero(_S), order)
    tmp790 = Array{Taylor1{_S}}(undef, size(dq))
    tmp790 .= Taylor1(zero(_S), order)
    pn2x = Array{Taylor1{_S}}(undef, size(X))
    pn2x .= Taylor1(zero(_S), order)
    pn2y = Array{Taylor1{_S}}(undef, size(Y))
    pn2y .= Taylor1(zero(_S), order)
    pn2z = Array{Taylor1{_S}}(undef, size(Z))
    pn2z .= Taylor1(zero(_S), order)
    tmp798 = Array{Taylor1{_S}}(undef, size(UU))
    tmp798 .= Taylor1(zero(_S), order)
    tmp801 = Array{Taylor1{_S}}(undef, size(X))
    tmp801 .= Taylor1(zero(_S), order)
    tmp803 = Array{Taylor1{_S}}(undef, size(Y))
    tmp803 .= Taylor1(zero(_S), order)
    tmp804 = Array{Taylor1{_S}}(undef, size(tmp801))
    tmp804 .= Taylor1(zero(_S), order)
    tmp806 = Array{Taylor1{_S}}(undef, size(Z))
    tmp806 .= Taylor1(zero(_S), order)
    tmp814 = Array{Taylor1{_S}}(undef, size(pn2x))
    tmp814 .= Taylor1(zero(_S), order)
    tmp815 = Array{Taylor1{_S}}(undef, size(tmp814))
    tmp815 .= Taylor1(zero(_S), order)
    tmp818 = Array{Taylor1{_S}}(undef, size(X))
    tmp818 .= Taylor1(zero(_S), order)
    tmp820 = Array{Taylor1{_S}}(undef, size(Y))
    tmp820 .= Taylor1(zero(_S), order)
    tmp822 = Array{Taylor1{_S}}(undef, size(Z))
    tmp822 .= Taylor1(zero(_S), order)
    tmp824 = Array{Taylor1{_S}}(undef, size(t31))
    tmp824 .= Taylor1(zero(_S), order)
    tmp1029 = Array{Taylor1{_S}}(undef, size(sin_ϕ))
    tmp1029 .= Taylor1(zero(_S), order)
    tmp1030 = Array{Taylor1{_S}}(undef, size(ϕ))
    tmp1030 .= Taylor1(zero(_S), order)
    tmp834 = Array{Taylor1{_S}}(undef, size(sin2_ϕ))
    tmp834 .= Taylor1(zero(_S), order)
    tmp840 = Array{Taylor1{_S}}(undef, size(sin_ϕ))
    tmp840 .= Taylor1(zero(_S), order)
    tmp842 = Array{Taylor1{_S}}(undef, size(sin3_ϕ))
    tmp842 .= Taylor1(zero(_S), order)
    tmp846 = Array{Taylor1{_S}}(undef, size(sin2_ϕ))
    tmp846 .= Taylor1(zero(_S), order)
    tmp848 = Array{Taylor1{_S}}(undef, size(Λ2))
    tmp848 .= Taylor1(zero(_S), order)
    tmp850 = Array{Taylor1{_S}}(undef, size(r_p2))
    tmp850 .= Taylor1(zero(_S), order)
    tmp852 = Array{Taylor1{_S}}(undef, size(Λ3))
    tmp852 .= Taylor1(zero(_S), order)
    tmp854 = Array{Taylor1{_S}}(undef, size(r_p1d2))
    tmp854 .= Taylor1(zero(_S), order)
    tmp856 = Array{Taylor1{_S}}(undef, size(cos_ϕ))
    tmp856 .= Taylor1(zero(_S), order)
    tmp858 = Array{Taylor1{_S}}(undef, size(cos_ϕ))
    tmp858 .= Taylor1(zero(_S), order)
    tmp861 = Array{Taylor1{_S}}(undef, size(Λ2j_div_r4))
    tmp861 .= Taylor1(zero(_S), order)
    tmp865 = Array{Taylor1{_S}}(undef, size(Λ3j_div_r5))
    tmp865 .= Taylor1(zero(_S), order)
    tmp870 = Array{Taylor1{_S}}(undef, size(X))
    tmp870 .= Taylor1(zero(_S), order)
    tmp872 = Array{Taylor1{_S}}(undef, size(Y))
    tmp872 .= Taylor1(zero(_S), order)
    tmp874 = Array{Taylor1{_S}}(undef, size(Z))
    tmp874 .= Taylor1(zero(_S), order)
    tmp903 = Array{Taylor1{_S}}(undef, size(μ))
    tmp903 .= Taylor1(zero(_S), order)
    tmp905 = Array{Taylor1{_S}}(undef, size(μ))
    tmp905 .= Taylor1(zero(_S), order)
    tmp907 = Array{Taylor1{_S}}(undef, size(μ))
    tmp907 .= Taylor1(zero(_S), order)
    tmp910 = Array{Taylor1{_S}}(undef, size(ui))
    tmp910 .= Taylor1(zero(_S), order)
    tmp912 = Array{Taylor1{_S}}(undef, size(vi))
    tmp912 .= Taylor1(zero(_S), order)
    tmp913 = Array{Taylor1{_S}}(undef, size(tmp910))
    tmp913 .= Taylor1(zero(_S), order)
    tmp915 = Array{Taylor1{_S}}(undef, size(wi))
    tmp915 .= Taylor1(zero(_S), order)
    for i = 1:N - 1
        xi[i] = Taylor1(identity(constant_term(ss16asteph_t[3i - 2])), order)
        yi[i] = Taylor1(identity(constant_term(ss16asteph_t[3i - 1])), order)
        zi[i] = Taylor1(identity(constant_term(ss16asteph_t[3i])), order)
        ui[i] = Taylor1(identity(constant_term(ss16asteph_t[3 * ((N - 1) + i) - 2])), order)
        vi[i] = Taylor1(identity(constant_term(ss16asteph_t[3 * ((N - 1) + i) - 1])), order)
        wi[i] = Taylor1(identity(constant_term(ss16asteph_t[3 * ((N - 1) + i)])), order)
        _3ui[i] = Taylor1(constant_term(3) * constant_term(ui[i]), order)
        _3vi[i] = Taylor1(constant_term(3) * constant_term(vi[i]), order)
        _3wi[i] = Taylor1(constant_term(3) * constant_term(wi[i]), order)
        X[i] = Taylor1(constant_term(xi[i]) - constant_term(q[1]), order)
        Y[i] = Taylor1(constant_term(yi[i]) - constant_term(q[2]), order)
        Z[i] = Taylor1(constant_term(zi[i]) - constant_term(q[3]), order)
        U[i] = Taylor1(constant_term(ui[i]) - constant_term(dq[1]), order)
        V[i] = Taylor1(constant_term(vi[i]) - constant_term(dq[2]), order)
        W[i] = Taylor1(constant_term(wi[i]) - constant_term(dq[3]), order)
        tmp784[1] = Taylor1(constant_term(4) * constant_term(dq[1]), order)
        _4U_m_3X[i] = Taylor1(constant_term(tmp784[1]) - constant_term(_3ui[i]), order)
        tmp787[2] = Taylor1(constant_term(4) * constant_term(dq[2]), order)
        _4V_m_3Y[i] = Taylor1(constant_term(tmp787[2]) - constant_term(_3vi[i]), order)
        tmp790[3] = Taylor1(constant_term(4) * constant_term(dq[3]), order)
        _4W_m_3Z[i] = Taylor1(constant_term(tmp790[3]) - constant_term(_3wi[i]), order)
        pn2x[i] = Taylor1(constant_term(X[i]) * constant_term(_4U_m_3X[i]), order)
        pn2y[i] = Taylor1(constant_term(Y[i]) * constant_term(_4V_m_3Y[i]), order)
        pn2z[i] = Taylor1(constant_term(Z[i]) * constant_term(_4W_m_3Z[i]), order)
        UU[i] = Taylor1(constant_term(ui[i]) * constant_term(dq[1]), order)
        VV[i] = Taylor1(constant_term(vi[i]) * constant_term(dq[2]), order)
        WW[i] = Taylor1(constant_term(wi[i]) * constant_term(dq[3]), order)
        tmp798[i] = Taylor1(constant_term(UU[i]) + constant_term(VV[i]), order)
        vi_dot_vj[i] = Taylor1(constant_term(tmp798[i]) + constant_term(WW[i]), order)
        tmp801[i] = Taylor1(constant_term(X[i]) ^ constant_term(2), order)
        tmp803[i] = Taylor1(constant_term(Y[i]) ^ constant_term(2), order)
        tmp804[i] = Taylor1(constant_term(tmp801[i]) + constant_term(tmp803[i]), order)
        tmp806[i] = Taylor1(constant_term(Z[i]) ^ constant_term(2), order)
        r_p2[i] = Taylor1(constant_term(tmp804[i]) + constant_term(tmp806[i]), order)
        r_p1d2[i] = Taylor1(sqrt(constant_term(r_p2[i])), order)
        r_p3d2[i] = Taylor1(constant_term(r_p2[i]) ^ constant_term(1.5), order)
        r_p7d2[i] = Taylor1(constant_term(r_p2[i]) ^ constant_term(3.5), order)
        newtonianCoeff[i] = Taylor1(constant_term(μ[i]) / constant_term(r_p3d2[i]), order)
        tmp814[i] = Taylor1(constant_term(pn2x[i]) + constant_term(pn2y[i]), order)
        tmp815[i] = Taylor1(constant_term(tmp814[i]) + constant_term(pn2z[i]), order)
        pn2[i] = Taylor1(constant_term(newtonianCoeff[i]) * constant_term(tmp815[i]), order)
        newtonian1b_Potential[i] = Taylor1(constant_term(μ[i]) / constant_term(r_p1d2[i]), order)
        if UJ_interaction[i]
            tmp818[i] = Taylor1(-(constant_term(X[i])), order)
            t31[i] = Taylor1(constant_term(tmp818[i]) * constant_term(M_[1, 3, i]), order)
            tmp820[i] = Taylor1(-(constant_term(Y[i])), order)
            t32[i] = Taylor1(constant_term(tmp820[i]) * constant_term(M_[2, 3, i]), order)
            tmp822[i] = Taylor1(-(constant_term(Z[i])), order)
            t33[i] = Taylor1(constant_term(tmp822[i]) * constant_term(M_[3, 3, i]), order)
            tmp824[i] = Taylor1(constant_term(t31[i]) + constant_term(t32[i]), order)
            r_sin_ϕ[i] = Taylor1(constant_term(tmp824[i]) + constant_term(t33[i]), order)
            sin_ϕ[i] = Taylor1(constant_term(r_sin_ϕ[i]) / constant_term(r_p1d2[i]), order)
            ϕ[i] = Taylor1(asin(constant_term(sin_ϕ[i])), order)
            tmp1029[i] = Taylor1(sqrt(1 - constant_term(sin_ϕ[i]) ^ 2), order)
            cos_ϕ[i] = Taylor1(cos(constant_term(ϕ[i])), order)
            tmp1030[i] = Taylor1(sin(constant_term(ϕ[i])), order)
            sin2_ϕ[i] = Taylor1(constant_term(sin_ϕ[i]) ^ constant_term(2), order)
            sin3_ϕ[i] = Taylor1(constant_term(sin_ϕ[i]) ^ constant_term(3), order)
            tmp834[i] = Taylor1(constant_term(1.5) * constant_term(sin2_ϕ[i]), order)
            P_2_sin_ϕ[i] = Taylor1(constant_term(tmp834[i]) - constant_term(0.5), order)
            ∂P_2_sin_ϕ[i] = Taylor1(constant_term(3) * constant_term(sin_ϕ[i]), order)
            tmp840[i] = Taylor1(constant_term(-1.5) * constant_term(sin_ϕ[i]), order)
            tmp842[i] = Taylor1(constant_term(2.5) * constant_term(sin3_ϕ[i]), order)
            P_3_sin_ϕ[i] = Taylor1(constant_term(tmp840[i]) + constant_term(tmp842[i]), order)
            tmp846[i] = Taylor1(constant_term(7.5) * constant_term(sin2_ϕ[i]), order)
            ∂P_3_sin_ϕ[i] = Taylor1(constant_term(-1.5) + constant_term(tmp846[i]), order)
            tmp848[i] = Taylor1(-(constant_term(Λ2[i])), order)
            tmp850[i] = Taylor1(constant_term(r_p2[i]) ^ constant_term(2), order)
            Λ2j_div_r4[i] = Taylor1(constant_term(tmp848[i]) / constant_term(tmp850[i]), order)
            tmp852[i] = Taylor1(-(constant_term(Λ3[i])), order)
            tmp854[i] = Taylor1(constant_term(r_p1d2[i]) ^ constant_term(5), order)
            Λ3j_div_r5[i] = Taylor1(constant_term(tmp852[i]) / constant_term(tmp854[i]), order)
            tmp856[i] = Taylor1(-(constant_term(cos_ϕ[i])), order)
            m_c_ϕ_∂P_2[i] = Taylor1(constant_term(tmp856[i]) * constant_term(∂P_2_sin_ϕ[i]), order)
            tmp858[i] = Taylor1(-(constant_term(cos_ϕ[i])), order)
            m_c_ϕ_∂P_3[i] = Taylor1(constant_term(tmp858[i]) * constant_term(∂P_3_sin_ϕ[i]), order)
            tmp861[i] = Taylor1(constant_term(Λ2j_div_r4[i]) * constant_term(3), order)
            F_J2_ξ[i] = Taylor1(constant_term(tmp861[i]) * constant_term(P_2_sin_ϕ[i]), order)
            F_J2_ζ[i] = Taylor1(constant_term(Λ2j_div_r4[i]) * constant_term(m_c_ϕ_∂P_2[i]), order)
            tmp865[i] = Taylor1(constant_term(Λ3j_div_r5[i]) * constant_term(4), order)
            F_J3_ξ[i] = Taylor1(constant_term(tmp865[i]) * constant_term(P_3_sin_ϕ[i]), order)
            F_J3_ζ[i] = Taylor1(constant_term(Λ3j_div_r5[i]) * constant_term(m_c_ϕ_∂P_3[i]), order)
            F_J_ξ[i] = Taylor1(constant_term(F_J2_ξ[i]) + constant_term(F_J3_ξ[i]), order)
            F_J_ζ[i] = Taylor1(constant_term(F_J2_ζ[i]) + constant_term(F_J3_ζ[i]), order)
            tmp870[i] = Taylor1(-(constant_term(X[i])), order)
            ξx[i] = Taylor1(constant_term(tmp870[i]) / constant_term(r_p1d2[i]), order)
            tmp872[i] = Taylor1(-(constant_term(Y[i])), order)
            ξy[i] = Taylor1(constant_term(tmp872[i]) / constant_term(r_p1d2[i]), order)
            tmp874[i] = Taylor1(-(constant_term(Z[i])), order)
            ξz[i] = Taylor1(constant_term(tmp874[i]) / constant_term(r_p1d2[i]), order)
            ηx1[i] = Taylor1(constant_term(M_[2, 3, i]) * constant_term(ξz[i]), order)
            ηy1[i] = Taylor1(constant_term(M_[3, 3, i]) * constant_term(ξx[i]), order)
            ηz1[i] = Taylor1(constant_term(M_[1, 3, i]) * constant_term(ξy[i]), order)
            ηx2[i] = Taylor1(constant_term(M_[3, 3, i]) * constant_term(ξy[i]), order)
            ηy2[i] = Taylor1(constant_term(M_[1, 3, i]) * constant_term(ξz[i]), order)
            ηz2[i] = Taylor1(constant_term(M_[2, 3, i]) * constant_term(ξx[i]), order)
            ηx[i] = Taylor1(constant_term(ηx1[i]) - constant_term(ηx2[i]), order)
            ηy[i] = Taylor1(constant_term(ηy1[i]) - constant_term(ηy2[i]), order)
            ηz[i] = Taylor1(constant_term(ηz1[i]) - constant_term(ηz2[i]), order)
            ζx1[i] = Taylor1(constant_term(ξy[i]) * constant_term(ηz[i]), order)
            ζy1[i] = Taylor1(constant_term(ξz[i]) * constant_term(ηx[i]), order)
            ζz1[i] = Taylor1(constant_term(ξx[i]) * constant_term(ηy[i]), order)
            ζx2[i] = Taylor1(constant_term(ξz[i]) * constant_term(ηy[i]), order)
            ζy2[i] = Taylor1(constant_term(ξx[i]) * constant_term(ηz[i]), order)
            ζz2[i] = Taylor1(constant_term(ξy[i]) * constant_term(ηx[i]), order)
            ζx[i] = Taylor1(constant_term(ζx1[i]) - constant_term(ζx2[i]), order)
            ζy[i] = Taylor1(constant_term(ζy1[i]) - constant_term(ζy2[i]), order)
            ζz[i] = Taylor1(constant_term(ζz1[i]) - constant_term(ζz2[i]), order)
            F_J2_x1[i] = Taylor1(constant_term(F_J_ξ[i]) * constant_term(ξx[i]), order)
            F_J2_y1[i] = Taylor1(constant_term(F_J_ξ[i]) * constant_term(ξy[i]), order)
            F_J2_z1[i] = Taylor1(constant_term(F_J_ξ[i]) * constant_term(ξz[i]), order)
            F_J2_x2[i] = Taylor1(constant_term(F_J_ζ[i]) * constant_term(ζx[i]), order)
            F_J2_y2[i] = Taylor1(constant_term(F_J_ζ[i]) * constant_term(ζy[i]), order)
            F_J2_z2[i] = Taylor1(constant_term(F_J_ζ[i]) * constant_term(ζz[i]), order)
            F_J2_x[i] = Taylor1(constant_term(F_J2_x1[i]) + constant_term(F_J2_x2[i]), order)
            F_J2_y[i] = Taylor1(constant_term(F_J2_y1[i]) + constant_term(F_J2_y2[i]), order)
            F_J2_z[i] = Taylor1(constant_term(F_J2_z1[i]) + constant_term(F_J2_z2[i]), order)
            tmp903[i] = Taylor1(constant_term(μ[i]) * constant_term(F_J2_x[i]), order)
            temp_accX_i[i] = Taylor1(-(constant_term(tmp903[i])), order)
            tmp905[i] = Taylor1(constant_term(μ[i]) * constant_term(F_J2_y[i]), order)
            temp_accY_i[i] = Taylor1(-(constant_term(tmp905[i])), order)
            tmp907[i] = Taylor1(constant_term(μ[i]) * constant_term(F_J2_z[i]), order)
            temp_accZ_i[i] = Taylor1(-(constant_term(tmp907[i])), order)
        end
        tmp910[i] = Taylor1(constant_term(ui[i]) ^ constant_term(2), order)
        tmp912[i] = Taylor1(constant_term(vi[i]) ^ constant_term(2), order)
        tmp913[i] = Taylor1(constant_term(tmp910[i]) + constant_term(tmp912[i]), order)
        tmp915[i] = Taylor1(constant_term(wi[i]) ^ constant_term(2), order)
        v2[i] = Taylor1(constant_term(tmp913[i]) + constant_term(tmp915[i]), order)
    end
    tmp918 = Taylor1(constant_term(q[4]) ^ constant_term(2), order)
    tmp920 = Taylor1(constant_term(q[5]) ^ constant_term(2), order)
    tmp921 = Taylor1(constant_term(tmp918) + constant_term(tmp920), order)
    tmp923 = Taylor1(constant_term(q[6]) ^ constant_term(2), order)
    v2[N] = Taylor1(constant_term(tmp921) + constant_term(tmp923), order)
    for i = 1:N - 1
        newtonianNb_Potential[N] = Taylor1(constant_term(newtonianNb_Potential[N]) + constant_term(newtonian1b_Potential[i]), order)
    end
    for i = j2_body_index
        accX = Taylor1(constant_term(accX) + constant_term(temp_accX_i[i]), order)
        accY = Taylor1(constant_term(accY) + constant_term(temp_accY_i[i]), order)
        accZ = Taylor1(constant_term(accZ) + constant_term(temp_accZ_i[i]), order)
    end
    tmp930 = Array{Taylor1{_S}}(undef, size(newtonianNb_Potential))
    tmp930 .= Taylor1(zero(_S), order)
    temp_005a = Array{Taylor1{_S}}(undef, size(newtonianNb_Potential_t))
    temp_005a .= Taylor1(zero(_S), order)
    tmp933 = Array{Taylor1{_S}}(undef, size(v2))
    tmp933 .= Taylor1(zero(_S), order)
    tmp935 = Array{Taylor1{_S}}(undef, size(vi_dot_vj))
    tmp935 .= Taylor1(zero(_S), order)
    tmp936 = Array{Taylor1{_S}}(undef, size(tmp933))
    tmp936 .= Taylor1(zero(_S), order)
    temp_005b = Array{Taylor1{_S}}(undef, size(tmp936))
    temp_005b .= Taylor1(zero(_S), order)
    temp_005 = Array{Taylor1{_S}}(undef, size(temp_005b))
    temp_005 .= Taylor1(zero(_S), order)
    temp_006a = Array{Taylor1{_S}}(undef, size(X))
    temp_006a .= Taylor1(zero(_S), order)
    temp_006b = Array{Taylor1{_S}}(undef, size(Y))
    temp_006b .= Taylor1(zero(_S), order)
    temp_006c = Array{Taylor1{_S}}(undef, size(Z))
    temp_006c .= Taylor1(zero(_S), order)
    tmp942 = Array{Taylor1{_S}}(undef, size(temp_006a))
    tmp942 .= Taylor1(zero(_S), order)
    temp_006d = Array{Taylor1{_S}}(undef, size(tmp942))
    temp_006d .= Taylor1(zero(_S), order)
    tmp945 = Array{Taylor1{_S}}(undef, size(temp_006d))
    tmp945 .= Taylor1(zero(_S), order)
    temp_006e = Array{Taylor1{_S}}(undef, size(tmp945))
    temp_006e .= Taylor1(zero(_S), order)
    temp_006 = Array{Taylor1{_S}}(undef, size(temp_006e))
    temp_006 .= Taylor1(zero(_S), order)
    temp_007a = Array{Taylor1{_S}}(undef, size(X))
    temp_007a .= Taylor1(zero(_S), order)
    temp_007b = Array{Taylor1{_S}}(undef, size(Y))
    temp_007b .= Taylor1(zero(_S), order)
    temp_007c = Array{Taylor1{_S}}(undef, size(Z))
    temp_007c .= Taylor1(zero(_S), order)
    tmp952 = Array{Taylor1{_S}}(undef, size(temp_007a))
    tmp952 .= Taylor1(zero(_S), order)
    temp_007d = Array{Taylor1{_S}}(undef, size(tmp952))
    temp_007d .= Taylor1(zero(_S), order)
    temp_007 = Array{Taylor1{_S}}(undef, size(temp_007d))
    temp_007 .= Taylor1(zero(_S), order)
    tmp956 = Array{Taylor1{_S}}(undef, size(temp_005))
    tmp956 .= Taylor1(zero(_S), order)
    temp_008 = Array{Taylor1{_S}}(undef, size(tmp956))
    temp_008 .= Taylor1(zero(_S), order)
    tmp958 = Array{Taylor1{_S}}(undef, size(temp_006))
    tmp958 .= Taylor1(zero(_S), order)
    temp_009 = Array{Taylor1{_S}}(undef, size(X))
    temp_009 .= Taylor1(zero(_S), order)
    temp_010 = Array{Taylor1{_S}}(undef, size(Y))
    temp_010 .= Taylor1(zero(_S), order)
    temp_011 = Array{Taylor1{_S}}(undef, size(Z))
    temp_011 .= Taylor1(zero(_S), order)
    temp_013a = Array{Taylor1{_S}}(undef, size(pn2))
    temp_013a .= Taylor1(zero(_S), order)
    temp_013b = Array{Taylor1{_S}}(undef, size(pn3))
    temp_013b .= Taylor1(zero(_S), order)
    tmp967 = Array{Taylor1{_S}}(undef, size(temp_013a))
    tmp967 .= Taylor1(zero(_S), order)
    temp_014a = Array{Taylor1{_S}}(undef, size(pn2))
    temp_014a .= Taylor1(zero(_S), order)
    temp_014b = Array{Taylor1{_S}}(undef, size(pn3))
    temp_014b .= Taylor1(zero(_S), order)
    tmp971 = Array{Taylor1{_S}}(undef, size(temp_014a))
    tmp971 .= Taylor1(zero(_S), order)
    temp_015a = Array{Taylor1{_S}}(undef, size(pn2))
    temp_015a .= Taylor1(zero(_S), order)
    temp_015b = Array{Taylor1{_S}}(undef, size(pn3))
    temp_015b .= Taylor1(zero(_S), order)
    tmp975 = Array{Taylor1{_S}}(undef, size(temp_015a))
    tmp975 .= Taylor1(zero(_S), order)
    for k = Base.OneTo(succ_approx_iter)
        for i = 1:N - 1
            pNxi[i] = Taylor1(identity(constant_term(acceph_t[3i - 2])), order)
            pNyi[i] = Taylor1(identity(constant_term(acceph_t[3i - 1])), order)
            pNzi[i] = Taylor1(identity(constant_term(acceph_t[3i])), order)
            tmp930[N] = Taylor1(constant_term(4) * constant_term(newtonianNb_Potential[N]), order)
            temp_005a[i] = Taylor1(constant_term(newtonianNb_Potential_t[i]) + constant_term(tmp930[N]), order)
            tmp933[i] = Taylor1(constant_term(2) * constant_term(v2[i]), order)
            tmp935[i] = Taylor1(constant_term(4) * constant_term(vi_dot_vj[i]), order)
            tmp936[i] = Taylor1(constant_term(tmp933[i]) - constant_term(tmp935[i]), order)
            temp_005b[i] = Taylor1(constant_term(tmp936[i]) + constant_term(v2[N]), order)
            temp_005[i] = Taylor1(constant_term(temp_005b[i]) - constant_term(temp_005a[i]), order)
            temp_006a[i] = Taylor1(constant_term(X[i]) * constant_term(ui[i]), order)
            temp_006b[i] = Taylor1(constant_term(Y[i]) * constant_term(vi[i]), order)
            temp_006c[i] = Taylor1(constant_term(Z[i]) * constant_term(wi[i]), order)
            tmp942[i] = Taylor1(constant_term(temp_006a[i]) + constant_term(temp_006b[i]), order)
            temp_006d[i] = Taylor1(constant_term(tmp942[i]) + constant_term(temp_006c[i]), order)
            tmp945[i] = Taylor1(constant_term(temp_006d[i]) ^ constant_term(2), order)
            temp_006e[i] = Taylor1(constant_term(tmp945[i]) / constant_term(r_p2[i]), order)
            temp_006[i] = Taylor1(constant_term(1.5) * constant_term(temp_006e[i]), order)
            temp_007a[i] = Taylor1(constant_term(X[i]) * constant_term(pNxi[i]), order)
            temp_007b[i] = Taylor1(constant_term(Y[i]) * constant_term(pNyi[i]), order)
            temp_007c[i] = Taylor1(constant_term(Z[i]) * constant_term(pNzi[i]), order)
            tmp952[i] = Taylor1(constant_term(temp_007a[i]) + constant_term(temp_007b[i]), order)
            temp_007d[i] = Taylor1(constant_term(tmp952[i]) + constant_term(temp_007c[i]), order)
            temp_007[i] = Taylor1(constant_term(0.5) * constant_term(temp_007d[i]), order)
            tmp956[i] = Taylor1(constant_term(temp_005[i]) + constant_term(temp_007[i]), order)
            temp_008[i] = Taylor1(constant_term(c_p2) + constant_term(tmp956[i]), order)
            tmp958[i] = Taylor1(constant_term(temp_006[i]) + constant_term(temp_008[i]), order)
            pn1[i] = Taylor1(constant_term(newtonianCoeff[i]) * constant_term(tmp958[i]), order)
            temp_009[i] = Taylor1(constant_term(X[i]) * constant_term(pn1[i]), order)
            temp_010[i] = Taylor1(constant_term(Y[i]) * constant_term(pn1[i]), order)
            temp_011[i] = Taylor1(constant_term(Z[i]) * constant_term(pn1[i]), order)
            pn3[i] = Taylor1(constant_term(3.5) * constant_term(newtonian1b_Potential[i]), order)
            temp_013a[i] = Taylor1(constant_term(pn2[i]) * constant_term(U[i]), order)
            temp_013b[i] = Taylor1(constant_term(pn3[i]) * constant_term(pNxi[i]), order)
            tmp967[i] = Taylor1(constant_term(temp_013a[i]) + constant_term(temp_013b[i]), order)
            pntempX[i] = Taylor1(constant_term(temp_009[i]) + constant_term(tmp967[i]), order)
            temp_014a[i] = Taylor1(constant_term(pn2[i]) * constant_term(V[i]), order)
            temp_014b[i] = Taylor1(constant_term(pn3[i]) * constant_term(pNyi[i]), order)
            tmp971[i] = Taylor1(constant_term(temp_014a[i]) + constant_term(temp_014b[i]), order)
            pntempY[i] = Taylor1(constant_term(temp_010[i]) + constant_term(tmp971[i]), order)
            temp_015a[i] = Taylor1(constant_term(pn2[i]) * constant_term(W[i]), order)
            temp_015b[i] = Taylor1(constant_term(pn3[i]) * constant_term(pNzi[i]), order)
            tmp975[i] = Taylor1(constant_term(temp_015a[i]) + constant_term(temp_015b[i]), order)
            pntempZ[i] = Taylor1(constant_term(temp_011[i]) + constant_term(tmp975[i]), order)
        end
        for i = 1:N - 1
            postNewtonX = Taylor1(constant_term(postNewtonX) + constant_term(pntempX[i]), order)
            postNewtonY = Taylor1(constant_term(postNewtonY) + constant_term(pntempY[i]), order)
            postNewtonZ = Taylor1(constant_term(postNewtonZ) + constant_term(pntempZ[i]), order)
        end
        postNewtonX = Taylor1(constant_term(postNewtonX) * constant_term(c_m2), order)
        postNewtonY = Taylor1(constant_term(postNewtonY) * constant_term(c_m2), order)
        postNewtonZ = Taylor1(constant_term(postNewtonZ) * constant_term(c_m2), order)
    end
    tmp983 = Taylor1(constant_term(Z[1]) * constant_term(V[1]), order)
    tmp984 = Taylor1(constant_term(Y[1]) * constant_term(W[1]), order)
    hx = Taylor1(constant_term(tmp983) - constant_term(tmp984), order)
    tmp986 = Taylor1(constant_term(X[1]) * constant_term(W[1]), order)
    tmp987 = Taylor1(constant_term(Z[1]) * constant_term(U[1]), order)
    hy = Taylor1(constant_term(tmp986) - constant_term(tmp987), order)
    tmp989 = Taylor1(constant_term(Y[1]) * constant_term(U[1]), order)
    tmp990 = Taylor1(constant_term(X[1]) * constant_term(V[1]), order)
    hz = Taylor1(constant_term(tmp989) - constant_term(tmp990), order)
    r_hs = Taylor1(sqrt(constant_term(r_p2[1])), order)
    runitx = Taylor1(constant_term(X[1]) / constant_term(r_hs), order)
    runity = Taylor1(constant_term(Y[2]) / constant_term(r_hs), order)
    runitz = Taylor1(constant_term(Z[3]) / constant_term(r_hs), order)
    tmp996 = Taylor1(constant_term(hy) * constant_term(runitz), order)
    tmp997 = Taylor1(constant_term(hz) * constant_term(runity), order)
    tunitx0 = Taylor1(constant_term(tmp996) - constant_term(tmp997), order)
    tmp999 = Taylor1(constant_term(hz) * constant_term(runitx), order)
    tmp1000 = Taylor1(constant_term(hx) * constant_term(runitz), order)
    tunity0 = Taylor1(constant_term(tmp999) - constant_term(tmp1000), order)
    tmp1002 = Taylor1(constant_term(hx) * constant_term(runity), order)
    tmp1003 = Taylor1(constant_term(hy) * constant_term(runitx), order)
    tunitz0 = Taylor1(constant_term(tmp1002) - constant_term(tmp1003), order)
    tmp1006 = Taylor1(constant_term(tunitx0) ^ constant_term(2), order)
    tmp1008 = Taylor1(constant_term(tunity0) ^ constant_term(2), order)
    tmp1009 = Taylor1(constant_term(tmp1006) + constant_term(tmp1008), order)
    tmp1011 = Taylor1(constant_term(tunitz0) ^ constant_term(2), order)
    tmp1012 = Taylor1(constant_term(tmp1009) + constant_term(tmp1011), order)
    hmag = Taylor1(sqrt(constant_term(tmp1012)), order)
    tunitx = Taylor1(constant_term(tunitx0) / constant_term(hmag), order)
    tunity = Taylor1(constant_term(tunity0) / constant_term(hmag), order)
    tunitz = Taylor1(constant_term(tunitz0) / constant_term(hmag), order)
    g_r = Taylor1(constant_term(r_hs) ^ constant_term(-2.0), order)
    A2_t_g_r = Taylor1(constant_term(q[7]) * constant_term(g_r), order)
    NGAx = Taylor1(constant_term(A2_t_g_r) * constant_term(tunitx), order)
    NGAy = Taylor1(constant_term(A2_t_g_r) * constant_term(tunity), order)
    NGAz = Taylor1(constant_term(A2_t_g_r) * constant_term(tunitz), order)
    tmp1023 = Taylor1(constant_term(postNewtonX) + constant_term(accX), order)
    dq[4] = Taylor1(constant_term(tmp1023) + constant_term(NGAx), order)
    tmp1025 = Taylor1(constant_term(postNewtonY) + constant_term(accY), order)
    dq[5] = Taylor1(constant_term(tmp1025) + constant_term(NGAy), order)
    tmp1027 = Taylor1(constant_term(postNewtonZ) + constant_term(accZ), order)
    dq[6] = Taylor1(constant_term(tmp1027) + constant_term(NGAz), order)
    dq[7] = Taylor1(identity(constant_term(zero_q_1)), order)
    for __idx = eachindex(q)
        (q[__idx]).coeffs[2] = (dq[__idx]).coeffs[1]
    end
    for ord = 1:order - 1
        ordnext = ord + 1
        TaylorSeries.identity!(postNewtonX, zero_q_1, ord)
        TaylorSeries.identity!(postNewtonY, zero_q_1, ord)
        TaylorSeries.identity!(postNewtonZ, zero_q_1, ord)
        TaylorSeries.identity!(accX, zero_q_1, ord)
        TaylorSeries.identity!(accY, zero_q_1, ord)
        TaylorSeries.identity!(accZ, zero_q_1, ord)
        TaylorSeries.identity!(dq[1], q[4], ord)
        TaylorSeries.identity!(dq[2], q[5], ord)
        TaylorSeries.identity!(dq[3], q[6], ord)
        TaylorSeries.identity!(newtonianNb_Potential[N], zero_q_1, ord)
        for j = j2_body_index
            TaylorSeries.identity!(t31[j], zero_q_1, ord)
            TaylorSeries.identity!(t32[j], zero_q_1, ord)
            TaylorSeries.identity!(t33[j], zero_q_1, ord)
            TaylorSeries.identity!(F_J2_x[j], zero_q_1, ord)
            TaylorSeries.identity!(F_J2_y[j], zero_q_1, ord)
            TaylorSeries.identity!(F_J2_z[j], zero_q_1, ord)
            TaylorSeries.identity!(F_J2_x1[j], zero_q_1, ord)
            TaylorSeries.identity!(F_J2_y1[j], zero_q_1, ord)
            TaylorSeries.identity!(F_J2_z1[j], zero_q_1, ord)
            TaylorSeries.identity!(F_J2_x2[j], zero_q_1, ord)
            TaylorSeries.identity!(F_J2_y2[j], zero_q_1, ord)
            TaylorSeries.identity!(F_J2_z2[j], zero_q_1, ord)
            TaylorSeries.identity!(temp_accX_i[j], zero_q_1, ord)
            TaylorSeries.identity!(temp_accY_i[j], zero_q_1, ord)
            TaylorSeries.identity!(temp_accZ_i[j], zero_q_1, ord)
            TaylorSeries.identity!(sin_ϕ[j], zero_q_1, ord)
            TaylorSeries.identity!(sin2_ϕ[j], zero_q_1, ord)
            TaylorSeries.identity!(sin3_ϕ[j], zero_q_1, ord)
            TaylorSeries.identity!(sin4_ϕ[j], zero_q_1, ord)
            TaylorSeries.identity!(ϕ[j], zero_q_1, ord)
            TaylorSeries.identity!(cos_ϕ[j], zero_q_1, ord)
            TaylorSeries.identity!(P_2_sin_ϕ[j], zero_q_1, ord)
            TaylorSeries.identity!(∂P_2_sin_ϕ[j], zero_q_1, ord)
            TaylorSeries.identity!(P_3_sin_ϕ[j], zero_q_1, ord)
            TaylorSeries.identity!(∂P_3_sin_ϕ[j], zero_q_1, ord)
            TaylorSeries.identity!(m_c_ϕ_∂P_2[j], zero_q_1, ord)
            TaylorSeries.identity!(m_c_ϕ_∂P_3[j], zero_q_1, ord)
            TaylorSeries.identity!(Λ2j_div_r4[j], zero_q_1, ord)
            TaylorSeries.identity!(Λ3j_div_r5[j], zero_q_1, ord)
            TaylorSeries.identity!(F_J_ξ[j], zero_q_1, ord)
            TaylorSeries.identity!(F_J_η[j], zero_q_1, ord)
            TaylorSeries.identity!(F_J_ζ[j], zero_q_1, ord)
            TaylorSeries.identity!(F_J2_ξ[j], zero_q_1, ord)
            TaylorSeries.identity!(F_J2_η[j], zero_q_1, ord)
            TaylorSeries.identity!(F_J2_ζ[j], zero_q_1, ord)
            TaylorSeries.identity!(F_J3_ξ[j], zero_q_1, ord)
            TaylorSeries.identity!(F_J3_η[j], zero_q_1, ord)
            TaylorSeries.identity!(F_J3_ζ[j], zero_q_1, ord)
            TaylorSeries.identity!(ξx[j], zero_q_1, ord)
            TaylorSeries.identity!(ξy[j], zero_q_1, ord)
            TaylorSeries.identity!(ξz[j], zero_q_1, ord)
            TaylorSeries.identity!(ηx[j], zero_q_1, ord)
            TaylorSeries.identity!(ηy[j], zero_q_1, ord)
            TaylorSeries.identity!(ηz[j], zero_q_1, ord)
            TaylorSeries.identity!(ηx1[j], zero_q_1, ord)
            TaylorSeries.identity!(ηy1[j], zero_q_1, ord)
            TaylorSeries.identity!(ηz1[j], zero_q_1, ord)
            TaylorSeries.identity!(ηx2[j], zero_q_1, ord)
            TaylorSeries.identity!(ηy2[j], zero_q_1, ord)
            TaylorSeries.identity!(ηz2[j], zero_q_1, ord)
            TaylorSeries.identity!(ζx[j], zero_q_1, ord)
            TaylorSeries.identity!(ζy[j], zero_q_1, ord)
            TaylorSeries.identity!(ζz[j], zero_q_1, ord)
            TaylorSeries.identity!(ζx1[j], zero_q_1, ord)
            TaylorSeries.identity!(ζy1[j], zero_q_1, ord)
            TaylorSeries.identity!(ζz1[j], zero_q_1, ord)
            TaylorSeries.identity!(ζx2[j], zero_q_1, ord)
            TaylorSeries.identity!(ζy2[j], zero_q_1, ord)
            TaylorSeries.identity!(ζz2[j], zero_q_1, ord)
        end
        for i = 1:N - 1
            TaylorSeries.identity!(xi[i], ss16asteph_t[3i - 2], ord)
            TaylorSeries.identity!(yi[i], ss16asteph_t[3i - 1], ord)
            TaylorSeries.identity!(zi[i], ss16asteph_t[3i], ord)
            TaylorSeries.identity!(ui[i], ss16asteph_t[3 * ((N - 1) + i) - 2], ord)
            TaylorSeries.identity!(vi[i], ss16asteph_t[3 * ((N - 1) + i) - 1], ord)
            TaylorSeries.identity!(wi[i], ss16asteph_t[3 * ((N - 1) + i)], ord)
            TaylorSeries.mul!(_3ui[i], 3, ui[i], ord)
            TaylorSeries.mul!(_3vi[i], 3, vi[i], ord)
            TaylorSeries.mul!(_3wi[i], 3, wi[i], ord)
            TaylorSeries.subst!(X[i], xi[i], q[1], ord)
            TaylorSeries.subst!(Y[i], yi[i], q[2], ord)
            TaylorSeries.subst!(Z[i], zi[i], q[3], ord)
            TaylorSeries.subst!(U[i], ui[i], dq[1], ord)
            TaylorSeries.subst!(V[i], vi[i], dq[2], ord)
            TaylorSeries.subst!(W[i], wi[i], dq[3], ord)
            TaylorSeries.mul!(tmp784[1], 4, dq[1], ord)
            TaylorSeries.subst!(_4U_m_3X[i], tmp784[1], _3ui[i], ord)
            TaylorSeries.mul!(tmp787[2], 4, dq[2], ord)
            TaylorSeries.subst!(_4V_m_3Y[i], tmp787[2], _3vi[i], ord)
            TaylorSeries.mul!(tmp790[3], 4, dq[3], ord)
            TaylorSeries.subst!(_4W_m_3Z[i], tmp790[3], _3wi[i], ord)
            TaylorSeries.mul!(pn2x[i], X[i], _4U_m_3X[i], ord)
            TaylorSeries.mul!(pn2y[i], Y[i], _4V_m_3Y[i], ord)
            TaylorSeries.mul!(pn2z[i], Z[i], _4W_m_3Z[i], ord)
            TaylorSeries.mul!(UU[i], ui[i], dq[1], ord)
            TaylorSeries.mul!(VV[i], vi[i], dq[2], ord)
            TaylorSeries.mul!(WW[i], wi[i], dq[3], ord)
            TaylorSeries.add!(tmp798[i], UU[i], VV[i], ord)
            TaylorSeries.add!(vi_dot_vj[i], tmp798[i], WW[i], ord)
            TaylorSeries.pow!(tmp801[i], X[i], 2, ord)
            TaylorSeries.pow!(tmp803[i], Y[i], 2, ord)
            TaylorSeries.add!(tmp804[i], tmp801[i], tmp803[i], ord)
            TaylorSeries.pow!(tmp806[i], Z[i], 2, ord)
            TaylorSeries.add!(r_p2[i], tmp804[i], tmp806[i], ord)
            TaylorSeries.sqrt!(r_p1d2[i], r_p2[i], ord)
            TaylorSeries.pow!(r_p3d2[i], r_p2[i], 1.5, ord)
            TaylorSeries.pow!(r_p7d2[i], r_p2[i], 3.5, ord)
            TaylorSeries.div!(newtonianCoeff[i], μ[i], r_p3d2[i], ord)
            TaylorSeries.add!(tmp814[i], pn2x[i], pn2y[i], ord)
            TaylorSeries.add!(tmp815[i], tmp814[i], pn2z[i], ord)
            TaylorSeries.mul!(pn2[i], newtonianCoeff[i], tmp815[i], ord)
            TaylorSeries.div!(newtonian1b_Potential[i], μ[i], r_p1d2[i], ord)
            if UJ_interaction[i]
                TaylorSeries.subst!(tmp818[i], X[i], ord)
                TaylorSeries.mul!(t31[i], tmp818[i], M_[1, 3, i], ord)
                TaylorSeries.subst!(tmp820[i], Y[i], ord)
                TaylorSeries.mul!(t32[i], tmp820[i], M_[2, 3, i], ord)
                TaylorSeries.subst!(tmp822[i], Z[i], ord)
                TaylorSeries.mul!(t33[i], tmp822[i], M_[3, 3, i], ord)
                TaylorSeries.add!(tmp824[i], t31[i], t32[i], ord)
                TaylorSeries.add!(r_sin_ϕ[i], tmp824[i], t33[i], ord)
                TaylorSeries.div!(sin_ϕ[i], r_sin_ϕ[i], r_p1d2[i], ord)
                TaylorSeries.asin!(ϕ[i], sin_ϕ[i], tmp1029[i], ord)
                TaylorSeries.sincos!(tmp1030[i], cos_ϕ[i], ϕ[i], ord)
                TaylorSeries.pow!(sin2_ϕ[i], sin_ϕ[i], 2, ord)
                TaylorSeries.pow!(sin3_ϕ[i], sin_ϕ[i], 3, ord)
                TaylorSeries.mul!(tmp834[i], 1.5, sin2_ϕ[i], ord)
                TaylorSeries.subst!(P_2_sin_ϕ[i], tmp834[i], 0.5, ord)
                TaylorSeries.mul!(∂P_2_sin_ϕ[i], 3, sin_ϕ[i], ord)
                TaylorSeries.mul!(tmp840[i], -1.5, sin_ϕ[i], ord)
                TaylorSeries.mul!(tmp842[i], 2.5, sin3_ϕ[i], ord)
                TaylorSeries.add!(P_3_sin_ϕ[i], tmp840[i], tmp842[i], ord)
                TaylorSeries.mul!(tmp846[i], 7.5, sin2_ϕ[i], ord)
                TaylorSeries.add!(∂P_3_sin_ϕ[i], -1.5, tmp846[i], ord)
                TaylorSeries.subst!(tmp848[i], Λ2[i], ord)
                TaylorSeries.pow!(tmp850[i], r_p2[i], 2, ord)
                TaylorSeries.div!(Λ2j_div_r4[i], tmp848[i], tmp850[i], ord)
                TaylorSeries.subst!(tmp852[i], Λ3[i], ord)
                TaylorSeries.pow!(tmp854[i], r_p1d2[i], 5, ord)
                TaylorSeries.div!(Λ3j_div_r5[i], tmp852[i], tmp854[i], ord)
                TaylorSeries.subst!(tmp856[i], cos_ϕ[i], ord)
                TaylorSeries.mul!(m_c_ϕ_∂P_2[i], tmp856[i], ∂P_2_sin_ϕ[i], ord)
                TaylorSeries.subst!(tmp858[i], cos_ϕ[i], ord)
                TaylorSeries.mul!(m_c_ϕ_∂P_3[i], tmp858[i], ∂P_3_sin_ϕ[i], ord)
                TaylorSeries.mul!(tmp861[i], Λ2j_div_r4[i], 3, ord)
                TaylorSeries.mul!(F_J2_ξ[i], tmp861[i], P_2_sin_ϕ[i], ord)
                TaylorSeries.mul!(F_J2_ζ[i], Λ2j_div_r4[i], m_c_ϕ_∂P_2[i], ord)
                TaylorSeries.mul!(tmp865[i], Λ3j_div_r5[i], 4, ord)
                TaylorSeries.mul!(F_J3_ξ[i], tmp865[i], P_3_sin_ϕ[i], ord)
                TaylorSeries.mul!(F_J3_ζ[i], Λ3j_div_r5[i], m_c_ϕ_∂P_3[i], ord)
                TaylorSeries.add!(F_J_ξ[i], F_J2_ξ[i], F_J3_ξ[i], ord)
                TaylorSeries.add!(F_J_ζ[i], F_J2_ζ[i], F_J3_ζ[i], ord)
                TaylorSeries.subst!(tmp870[i], X[i], ord)
                TaylorSeries.div!(ξx[i], tmp870[i], r_p1d2[i], ord)
                TaylorSeries.subst!(tmp872[i], Y[i], ord)
                TaylorSeries.div!(ξy[i], tmp872[i], r_p1d2[i], ord)
                TaylorSeries.subst!(tmp874[i], Z[i], ord)
                TaylorSeries.div!(ξz[i], tmp874[i], r_p1d2[i], ord)
                TaylorSeries.mul!(ηx1[i], M_[2, 3, i], ξz[i], ord)
                TaylorSeries.mul!(ηy1[i], M_[3, 3, i], ξx[i], ord)
                TaylorSeries.mul!(ηz1[i], M_[1, 3, i], ξy[i], ord)
                TaylorSeries.mul!(ηx2[i], M_[3, 3, i], ξy[i], ord)
                TaylorSeries.mul!(ηy2[i], M_[1, 3, i], ξz[i], ord)
                TaylorSeries.mul!(ηz2[i], M_[2, 3, i], ξx[i], ord)
                TaylorSeries.subst!(ηx[i], ηx1[i], ηx2[i], ord)
                TaylorSeries.subst!(ηy[i], ηy1[i], ηy2[i], ord)
                TaylorSeries.subst!(ηz[i], ηz1[i], ηz2[i], ord)
                TaylorSeries.mul!(ζx1[i], ξy[i], ηz[i], ord)
                TaylorSeries.mul!(ζy1[i], ξz[i], ηx[i], ord)
                TaylorSeries.mul!(ζz1[i], ξx[i], ηy[i], ord)
                TaylorSeries.mul!(ζx2[i], ξz[i], ηy[i], ord)
                TaylorSeries.mul!(ζy2[i], ξx[i], ηz[i], ord)
                TaylorSeries.mul!(ζz2[i], ξy[i], ηx[i], ord)
                TaylorSeries.subst!(ζx[i], ζx1[i], ζx2[i], ord)
                TaylorSeries.subst!(ζy[i], ζy1[i], ζy2[i], ord)
                TaylorSeries.subst!(ζz[i], ζz1[i], ζz2[i], ord)
                TaylorSeries.mul!(F_J2_x1[i], F_J_ξ[i], ξx[i], ord)
                TaylorSeries.mul!(F_J2_y1[i], F_J_ξ[i], ξy[i], ord)
                TaylorSeries.mul!(F_J2_z1[i], F_J_ξ[i], ξz[i], ord)
                TaylorSeries.mul!(F_J2_x2[i], F_J_ζ[i], ζx[i], ord)
                TaylorSeries.mul!(F_J2_y2[i], F_J_ζ[i], ζy[i], ord)
                TaylorSeries.mul!(F_J2_z2[i], F_J_ζ[i], ζz[i], ord)
                TaylorSeries.add!(F_J2_x[i], F_J2_x1[i], F_J2_x2[i], ord)
                TaylorSeries.add!(F_J2_y[i], F_J2_y1[i], F_J2_y2[i], ord)
                TaylorSeries.add!(F_J2_z[i], F_J2_z1[i], F_J2_z2[i], ord)
                TaylorSeries.mul!(tmp903[i], μ[i], F_J2_x[i], ord)
                TaylorSeries.subst!(temp_accX_i[i], tmp903[i], ord)
                TaylorSeries.mul!(tmp905[i], μ[i], F_J2_y[i], ord)
                TaylorSeries.subst!(temp_accY_i[i], tmp905[i], ord)
                TaylorSeries.mul!(tmp907[i], μ[i], F_J2_z[i], ord)
                TaylorSeries.subst!(temp_accZ_i[i], tmp907[i], ord)
            end
            TaylorSeries.pow!(tmp910[i], ui[i], 2, ord)
            TaylorSeries.pow!(tmp912[i], vi[i], 2, ord)
            TaylorSeries.add!(tmp913[i], tmp910[i], tmp912[i], ord)
            TaylorSeries.pow!(tmp915[i], wi[i], 2, ord)
            TaylorSeries.add!(v2[i], tmp913[i], tmp915[i], ord)
        end
        TaylorSeries.pow!(tmp918, q[4], 2, ord)
        TaylorSeries.pow!(tmp920, q[5], 2, ord)
        TaylorSeries.add!(tmp921, tmp918, tmp920, ord)
        TaylorSeries.pow!(tmp923, q[6], 2, ord)
        TaylorSeries.add!(v2[N], tmp921, tmp923, ord)
        for i = 1:N - 1
            TaylorSeries.add!(newtonianNb_Potential[N], newtonianNb_Potential[N], newtonian1b_Potential[i], ord)
        end
        for i = j2_body_index
            TaylorSeries.add!(accX, accX, temp_accX_i[i], ord)
            TaylorSeries.add!(accY, accY, temp_accY_i[i], ord)
            TaylorSeries.add!(accZ, accZ, temp_accZ_i[i], ord)
        end
        for k = Base.OneTo(succ_approx_iter)
            for i = 1:N - 1
                TaylorSeries.identity!(pNxi[i], acceph_t[3i - 2], ord)
                TaylorSeries.identity!(pNyi[i], acceph_t[3i - 1], ord)
                TaylorSeries.identity!(pNzi[i], acceph_t[3i], ord)
                TaylorSeries.mul!(tmp930[N], 4, newtonianNb_Potential[N], ord)
                TaylorSeries.add!(temp_005a[i], newtonianNb_Potential_t[i], tmp930[N], ord)
                TaylorSeries.mul!(tmp933[i], 2, v2[i], ord)
                TaylorSeries.mul!(tmp935[i], 4, vi_dot_vj[i], ord)
                TaylorSeries.subst!(tmp936[i], tmp933[i], tmp935[i], ord)
                TaylorSeries.add!(temp_005b[i], tmp936[i], v2[N], ord)
                TaylorSeries.subst!(temp_005[i], temp_005b[i], temp_005a[i], ord)
                TaylorSeries.mul!(temp_006a[i], X[i], ui[i], ord)
                TaylorSeries.mul!(temp_006b[i], Y[i], vi[i], ord)
                TaylorSeries.mul!(temp_006c[i], Z[i], wi[i], ord)
                TaylorSeries.add!(tmp942[i], temp_006a[i], temp_006b[i], ord)
                TaylorSeries.add!(temp_006d[i], tmp942[i], temp_006c[i], ord)
                TaylorSeries.pow!(tmp945[i], temp_006d[i], 2, ord)
                TaylorSeries.div!(temp_006e[i], tmp945[i], r_p2[i], ord)
                TaylorSeries.mul!(temp_006[i], 1.5, temp_006e[i], ord)
                TaylorSeries.mul!(temp_007a[i], X[i], pNxi[i], ord)
                TaylorSeries.mul!(temp_007b[i], Y[i], pNyi[i], ord)
                TaylorSeries.mul!(temp_007c[i], Z[i], pNzi[i], ord)
                TaylorSeries.add!(tmp952[i], temp_007a[i], temp_007b[i], ord)
                TaylorSeries.add!(temp_007d[i], tmp952[i], temp_007c[i], ord)
                TaylorSeries.mul!(temp_007[i], 0.5, temp_007d[i], ord)
                TaylorSeries.add!(tmp956[i], temp_005[i], temp_007[i], ord)
                TaylorSeries.add!(temp_008[i], c_p2, tmp956[i], ord)
                TaylorSeries.add!(tmp958[i], temp_006[i], temp_008[i], ord)
                TaylorSeries.mul!(pn1[i], newtonianCoeff[i], tmp958[i], ord)
                TaylorSeries.mul!(temp_009[i], X[i], pn1[i], ord)
                TaylorSeries.mul!(temp_010[i], Y[i], pn1[i], ord)
                TaylorSeries.mul!(temp_011[i], Z[i], pn1[i], ord)
                TaylorSeries.mul!(pn3[i], 3.5, newtonian1b_Potential[i], ord)
                TaylorSeries.mul!(temp_013a[i], pn2[i], U[i], ord)
                TaylorSeries.mul!(temp_013b[i], pn3[i], pNxi[i], ord)
                TaylorSeries.add!(tmp967[i], temp_013a[i], temp_013b[i], ord)
                TaylorSeries.add!(pntempX[i], temp_009[i], tmp967[i], ord)
                TaylorSeries.mul!(temp_014a[i], pn2[i], V[i], ord)
                TaylorSeries.mul!(temp_014b[i], pn3[i], pNyi[i], ord)
                TaylorSeries.add!(tmp971[i], temp_014a[i], temp_014b[i], ord)
                TaylorSeries.add!(pntempY[i], temp_010[i], tmp971[i], ord)
                TaylorSeries.mul!(temp_015a[i], pn2[i], W[i], ord)
                TaylorSeries.mul!(temp_015b[i], pn3[i], pNzi[i], ord)
                TaylorSeries.add!(tmp975[i], temp_015a[i], temp_015b[i], ord)
                TaylorSeries.add!(pntempZ[i], temp_011[i], tmp975[i], ord)
            end
            for i = 1:N - 1
                TaylorSeries.add!(postNewtonX, postNewtonX, pntempX[i], ord)
                TaylorSeries.add!(postNewtonY, postNewtonY, pntempY[i], ord)
                TaylorSeries.add!(postNewtonZ, postNewtonZ, pntempZ[i], ord)
            end
            TaylorSeries.mul!(postNewtonX, postNewtonX, c_m2, ord)
            TaylorSeries.mul!(postNewtonY, postNewtonY, c_m2, ord)
            TaylorSeries.mul!(postNewtonZ, postNewtonZ, c_m2, ord)
        end
        TaylorSeries.mul!(tmp983, Z[1], V[1], ord)
        TaylorSeries.mul!(tmp984, Y[1], W[1], ord)
        TaylorSeries.subst!(hx, tmp983, tmp984, ord)
        TaylorSeries.mul!(tmp986, X[1], W[1], ord)
        TaylorSeries.mul!(tmp987, Z[1], U[1], ord)
        TaylorSeries.subst!(hy, tmp986, tmp987, ord)
        TaylorSeries.mul!(tmp989, Y[1], U[1], ord)
        TaylorSeries.mul!(tmp990, X[1], V[1], ord)
        TaylorSeries.subst!(hz, tmp989, tmp990, ord)
        TaylorSeries.sqrt!(r_hs, r_p2[1], ord)
        TaylorSeries.div!(runitx, X[1], r_hs, ord)
        TaylorSeries.div!(runity, Y[2], r_hs, ord)
        TaylorSeries.div!(runitz, Z[3], r_hs, ord)
        TaylorSeries.mul!(tmp996, hy, runitz, ord)
        TaylorSeries.mul!(tmp997, hz, runity, ord)
        TaylorSeries.subst!(tunitx0, tmp996, tmp997, ord)
        TaylorSeries.mul!(tmp999, hz, runitx, ord)
        TaylorSeries.mul!(tmp1000, hx, runitz, ord)
        TaylorSeries.subst!(tunity0, tmp999, tmp1000, ord)
        TaylorSeries.mul!(tmp1002, hx, runity, ord)
        TaylorSeries.mul!(tmp1003, hy, runitx, ord)
        TaylorSeries.subst!(tunitz0, tmp1002, tmp1003, ord)
        TaylorSeries.pow!(tmp1006, tunitx0, 2, ord)
        TaylorSeries.pow!(tmp1008, tunity0, 2, ord)
        TaylorSeries.add!(tmp1009, tmp1006, tmp1008, ord)
        TaylorSeries.pow!(tmp1011, tunitz0, 2, ord)
        TaylorSeries.add!(tmp1012, tmp1009, tmp1011, ord)
        TaylorSeries.sqrt!(hmag, tmp1012, ord)
        TaylorSeries.div!(tunitx, tunitx0, hmag, ord)
        TaylorSeries.div!(tunity, tunity0, hmag, ord)
        TaylorSeries.div!(tunitz, tunitz0, hmag, ord)
        TaylorSeries.pow!(g_r, r_hs, -2.0, ord)
        TaylorSeries.mul!(A2_t_g_r, q[7], g_r, ord)
        TaylorSeries.mul!(NGAx, A2_t_g_r, tunitx, ord)
        TaylorSeries.mul!(NGAy, A2_t_g_r, tunity, ord)
        TaylorSeries.mul!(NGAz, A2_t_g_r, tunitz, ord)
        TaylorSeries.add!(tmp1023, postNewtonX, accX, ord)
        TaylorSeries.add!(dq[4], tmp1023, NGAx, ord)
        TaylorSeries.add!(tmp1025, postNewtonY, accY, ord)
        TaylorSeries.add!(dq[5], tmp1025, NGAy, ord)
        TaylorSeries.add!(tmp1027, postNewtonZ, accZ, ord)
        TaylorSeries.add!(dq[6], tmp1027, NGAz, ord)
        TaylorSeries.identity!(dq[7], zero_q_1, ord)
        for __idx = eachindex(q)
            (q[__idx]).coeffs[ordnext + 1] = (dq[__idx]).coeffs[ordnext] / ordnext
        end
    end
    return nothing
end
