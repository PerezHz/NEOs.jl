#auxiliary function to evaluate pl.eph. case: no jet transport
function evaleph(eph::TaylorInterpolant, t::Taylor1, ::Taylor1{T}) where {T<:Real}
    return eph(t)
end

# #auxiliary function to evaluate pl.eph. case: Taylor1 jet transport
function evaleph(eph::TaylorInterpolant, t::Taylor1, q1::Taylor1{Taylor1{T}}) where {T<:Real}
    return map(x->Taylor1( x.coeffs*one(q1[0]) ), eph(t))
end

#auxiliary function to evaluate pl.eph. case: TaylorN jet transport
function evaleph(eph::TaylorInterpolant, t::Taylor1, q1::Taylor1{TaylorN{T}}) where {T<:Real}
    return map(x->Taylor1( x.coeffs*one(q1[0]) ), eph(t))
end

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
    local jd0 = params[4]
    local ss16asteph_t = evaleph(params[1], t, q[1]) # params[2](t)*one(q[1]) #ss16asteph(t)
    local acceph_t = evaleph(params[2], t, q[1]) # params[2](t)*one(q[1]) #acc_eph(t)
    local newtonianNb_Potential_t = evaleph(params[3], t, q[1]) # params[3](t)*one(q[1]) #newtonianNb_Potential(t), massive bodies
    local S = eltype(q[1])
    local N = length(μ) # number of bodies, including NEA
    local _1_to_Nm1 = 1:(N-1) # iterator over all bodies

    # parameters related to speed of light, c
    local c_p2 = 29979.063823897606 # c^2 = 29979.063823897606 au^2/d^2
    local c_m2 = 3.3356611996764786e-5 # c^-2 = 3.3356611996764786e-5 d^2/au^2

    local zero_q_1 = zero(q[1])

    X = Array{Taylor1{S}}(undef, N)
    Y = Array{Taylor1{S}}(undef, N)
    Z = Array{Taylor1{S}}(undef, N)

    r_p2 = Array{Taylor1{S}}(undef, N)
    r_p3d2 = Array{Taylor1{S}}(undef, N)
    r_p7d2 = Array{Taylor1{S}}(undef, N)

    newtonianCoeff = Array{Taylor1{S}}(undef, N)

    ui = Array{Taylor1{S}}(undef, N-1)
    vi = Array{Taylor1{S}}(undef, N-1)
    wi = Array{Taylor1{S}}(undef, N-1)

    #post-Newtonian stuff
    U = Array{Taylor1{S}}(undef, N)
    V = Array{Taylor1{S}}(undef, N)
    W = Array{Taylor1{S}}(undef, N)

    _4U_m_3X = Array{Taylor1{S}}(undef, N)
    _4V_m_3Y = Array{Taylor1{S}}(undef, N)
    _4W_m_3Z = Array{Taylor1{S}}(undef, N)

    UU = Array{Taylor1{S}}(undef, N)
    VV = Array{Taylor1{S}}(undef, N)
    WW = Array{Taylor1{S}}(undef, N)

    r_p1d2 = Array{Taylor1{S}}(undef, N)

    newtonianNb_Potential = Array{Taylor1{S}}(undef, N)
    newtonian1b_Potential = Array{Taylor1{S}}(undef, N)
    newtonianCoeff = Array{Taylor1{S}}(undef, N)
    newton_acc_X = Array{Taylor1{S}}(undef, N)
    newton_acc_Y = Array{Taylor1{S}}(undef, N)
    newton_acc_Z = Array{Taylor1{S}}(undef, N)

    v2 = Array{Taylor1{S}}(undef, N)
    vi_dot_vj = Array{Taylor1{S}}(undef, N)
    pn2 = Array{Taylor1{S}}(undef, N)
    pn3 = Array{Taylor1{S}}(undef, N)
    _4ϕj = Array{Taylor1{S}}(undef, N)
    ϕi_plus_4ϕj = Array{Taylor1{S}}(undef, N)
    sj2_plus_2si2_minus_4vivj = Array{Taylor1{S}}(undef, N)
    ϕs_and_vs = Array{Taylor1{S}}(undef, N)
    U_t_pn2 = Array{Taylor1{S}}(undef, N)
    V_t_pn2 = Array{Taylor1{S}}(undef, N)
    W_t_pn2 = Array{Taylor1{S}}(undef, N)
    pn1t1_7 = Array{Taylor1{S}}(undef, N)

    pntempX = zero_q_1
    pntempY = zero_q_1
    pntempZ = zero_q_1
    pn1 = Array{Taylor1{S}}(undef, N)
    X_t_pn1 = Array{Taylor1{S}}(undef, N)
    Y_t_pn1 = Array{Taylor1{S}}(undef, N)
    Z_t_pn1 = Array{Taylor1{S}}(undef, N)
    pNX_t_pn3 = Array{Taylor1{S}}(undef, N)
    pNY_t_pn3 = Array{Taylor1{S}}(undef, N)
    pNZ_t_pn3 = Array{Taylor1{S}}(undef, N)
    pNX_t_X = Array{Taylor1{S}}(undef, N)
    pNY_t_Y = Array{Taylor1{S}}(undef, N)
    pNZ_t_Z = Array{Taylor1{S}}(undef, N)

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
    # temp_accX_j = Array{Taylor1{S}}(undef, N)
    # temp_accY_j = Array{Taylor1{S}}(undef, N)
    # temp_accZ_j = Array{Taylor1{S}}(undef, N)
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

    # extended-body accelerations
    accX = zero_q_1
    accY = zero_q_1
    accZ = zero_q_1

    # rotations to and from Earth, Sun and Moon pole-oriented frames
    local dsj2k = t+(jd0-2.451545e6) # days since J2000.0 = 2.451545e6
    local M_ = Array{Taylor1{S}}(undef, 3, 3, N)
    local M_[:,:,ea] = t2c_jpl_de430(dsj2k)

    dq[1] = q[4]
    dq[2] = q[5]
    dq[3] = q[6]

    newtonianNb_Potential[N] = zero_q_1

    #compute point-mass Newtonian accelerations, all bodies
    for i in _1_to_Nm1
        ui[i] = ss16asteph_t[3(N-1+i)-2]
        vi[i] = ss16asteph_t[3(N-1+i)-1]
        wi[i] = ss16asteph_t[3(N-1+i)  ]

        X[i] = ss16asteph_t[3i-2]-q[1]
        Y[i] = ss16asteph_t[3i-1]-q[2]
        Z[i] = ss16asteph_t[3i  ]-q[3]

        U[i] = ui[i]-dq[1]
        V[i] = vi[i]-dq[2]
        W[i] = wi[i]-dq[3]

        _4U_m_3X[i] = (4dq[1]) - (3ui[i])
        _4V_m_3Y[i] = (4dq[2]) - (3vi[i])
        _4W_m_3Z[i] = (4dq[3]) - (3wi[i])

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

        newton_acc_X[i] = X[i]*newtonianCoeff[i]
        newton_acc_Y[i] = Y[i]*newtonianCoeff[i]
        newton_acc_Z[i] = Z[i]*newtonianCoeff[i]

        newtonian1b_Potential[i] = μ[i]/r_p1d2[i]
        pn3[i] = 3.5newtonian1b_Potential[i]
        U_t_pn2[i] = pn2[i]*U[i]
        V_t_pn2[i] = pn2[i]*V[i]
        W_t_pn2[i] = pn2[i]*W[i]

        #J2 accelerations, if i-th body is flattened
        if UJ_interaction[i]
            # # rotate from inertial frame to extended-body frame
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
            F_J_ξ[i] = F_J2_ξ[i] #+ F_J3_ξ[i]
            #F_J_η[i] = zero_q_1
            F_J_ζ[i] = F_J2_ζ[i] #+ F_J3_ζ[i]
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
        end # if UJ_interaction[i]
        v2[i] = ( (ui[i]^2)+(vi[i]^2) ) + (wi[i]^2)
    end #for, i
    v2[N] = ( (q[4]^2)+(q[5]^2) ) + (q[6]^2)

    for i in _1_to_Nm1
        temp_004 = newtonian1b_Potential[i] + newtonianNb_Potential[N]
        newtonianNb_Potential[N] = temp_004
        if UJ_interaction[i]
            # # # add result to total acceleration on upon j-th body figure due to i-th point mass
            # # @show "acc",j,"+μ",i,"Λ2",j
            # temp_accX_j[i,j] = accX[j] + (μ[i]*F_J2_x[i,j])
            # accX[j] = temp_accX_j[i,j]
            # temp_accY_j[i,j] = accY[j] + (μ[i]*F_J2_y[i,j])
            # accY[j] = temp_accY_j[i,j]
            # temp_accZ_j[i,j] = accZ[j] + (μ[i]*F_J2_z[i,j])
            # accZ[j] = temp_accZ_j[i,j]

            # # reaction force on i-th body
            # @show "acc",i,"-μ",j,"Λ2",j
            temp_accX_i[i] = accX - (μ[i]*F_J2_x[i])
            accX = temp_accX_i[i]
            temp_accY_i[i] = accY - (μ[i]*F_J2_y[i])
            accY = temp_accY_i[i]
            temp_accZ_i[i] = accZ - (μ[i]*F_J2_z[i])
            accZ = temp_accZ_i[i]
        end
    end

    #post-Newtonian accelerations due to Sun, Moon and planets (Mercury through Neptune)
    #Moyer, 1971, page 7 eq. 35
    # post-Newtonian iterative procedure setup and initialization
    _4ϕj[N] = 4newtonianNb_Potential[N]
    for i in 1:10
        ϕi_plus_4ϕj[i] = newtonianNb_Potential_t[i] + _4ϕj[N]
        sj2_plus_2si2_minus_4vivj[i] = ( (2v2[i]) - (4vi_dot_vj[i]) ) + v2[N]
        ϕs_and_vs[i] = sj2_plus_2si2_minus_4vivj[i] - ϕi_plus_4ϕj[i]
        Xij_t_Ui = X[i]*ui[i]
        Yij_t_Vi = Y[i]*vi[i]
        Zij_t_Wi = Z[i]*wi[i]
        Rij_dot_Vi = ( Xij_t_Ui+Yij_t_Vi ) + Zij_t_Wi
        # the expression below inside the (...)^2 should have a minus sign in front of the numerator,
        # but upon squaring it is eliminated, so at the end of the day, it is irrelevant ;)
        pn1t7 = (Rij_dot_Vi^2)/r_p2[i]
        pn1t2_7 = ϕs_and_vs[i] - (1.5pn1t7)
        pn1t1_7[i] = c_p2 + pn1t2_7

        pNX_t_X[i] = acceph_t[3i-2]*X[i]
        pNY_t_Y[i] = acceph_t[3i-1]*Y[i]
        pNZ_t_Z[i] = acceph_t[3i  ]*Z[i]
        pn1[i] = (  pn1t1_7[i]  +  (0.5*( (pNX_t_X[i]+pNY_t_Y[i]) + pNZ_t_Z[i] ))  )

        X_t_pn1[i] = newton_acc_X[i]*pn1[i]
        Y_t_pn1[i] = newton_acc_Y[i]*pn1[i]
        Z_t_pn1[i] = newton_acc_Z[i]*pn1[i]

        pNX_t_pn3[i] = acceph_t[3i-2]*pn3[i]
        pNY_t_pn3[i] = acceph_t[3i-1]*pn3[i]
        pNZ_t_pn3[i] = acceph_t[3i  ]*pn3[i]
    end #for i
    for i in 1:10
        termpnx = ( X_t_pn1[i] + (U_t_pn2[i]+pNX_t_pn3[i]) )
        sumpnx = pntempX + termpnx
        pntempX = sumpnx
        termpny = ( Y_t_pn1[i] + (V_t_pn2[i]+pNY_t_pn3[i]) )
        sumpny = pntempY + termpny
        pntempY = sumpny
        termpnz = ( Z_t_pn1[i] + (W_t_pn2[i]+pNZ_t_pn3[i]) )
        sumpnz = pntempZ + termpnz
        pntempZ = sumpnz
    end
    # compute Newtonian accelerations due to Pluto and 16 asteroid perturbers
    for i in 11:27
        X_t_pn1[i] = c_p2*newton_acc_X[i]
        Y_t_pn1[i] = c_p2*newton_acc_Y[i]
        Z_t_pn1[i] = c_p2*newton_acc_Z[i]
    end #for i
    for i in 11:27
        termpnx = X_t_pn1[i]
        sumpnx = pntempX + termpnx
        pntempX = sumpnx
        termpny = Y_t_pn1[i]
        sumpny = pntempY + termpny
        pntempY = sumpny
        termpnz = Z_t_pn1[i]
        sumpnz = pntempZ + termpnz
        pntempZ = sumpnz
    end
    postNewtonX = pntempX*c_m2
    postNewtonY = pntempY*c_m2
    postNewtonZ = pntempZ*c_m2

    # compute non-gravitational acceleration
    hx = (Y[1]*W[1])-(Z[1]*V[1])
    hy = (Z[1]*U[1])-(X[1]*W[1])
    hz = (X[1]*V[1])-(Y[1]*U[1])

    #cartesian components of transversal unit vector:
    tunitx0 = (hz*Y[1]) - (hy*Z[1]) # Note: Y[1] = y_Sun - y_Apophis, etc.
    tunity0 = (hx*Z[1]) - (hz*X[1])
    tunitz0 = (hy*X[1]) - (hx*Y[1])
    hmag = sqrt( ((tunitx0^2)+(tunity0^2))+(tunitz0^2) )
    tunitx = tunitx0/hmag
    tunity = tunity0/hmag
    tunitz = tunitz0/hmag

    # evaluate non-grav acceleration of NEA (Yarkovsky):
    g_r = r_p2[1]
    A2_t_g_r = q[7]/g_r

    NGAx = A2_t_g_r*tunitx
    NGAy = A2_t_g_r*tunity
    NGAz = A2_t_g_r*tunitz

    dq[4] = ( postNewtonX + accX ) + NGAx
    dq[5] = ( postNewtonY + accY ) + NGAy
    dq[6] = ( postNewtonZ + accZ ) + NGAz

    dq[7] = zero_q_1

    nothing
end

function RNp1BP_pN_A_J23E_J2S_ng_eph_threads!(dq, q, params, t)
    local jd0 = params[4]
    local ss16asteph_t = evaleph(params[1], t, q[1]) # params[2](t)*one(q[1]) #ss16asteph(t)
    local acceph_t = evaleph(params[2], t, q[1]) # params[2](t)*one(q[1]) #acc_eph(t)
    local newtonianNb_Potential_t = evaleph(params[3], t, q[1]) # params[3](t)*one(q[1]) #newtonianNb_Potential(t), massive bodies
    local S = eltype(q[1])
    local N = length(μ) # number of bodies, including NEA
    local _1_to_Nm1 = 1:(N-1) # iterator over all bodies

    # parameters related to speed of light, c
    local c_p2 = 29979.063823897606 # c^2 = 29979.063823897606 au^2/d^2
    local c_m2 = 3.3356611996764786e-5 # c^-2 = 3.3356611996764786e-5 d^2/au^2

    local zero_q_1 = zero(q[1])

    X = Array{Taylor1{S}}(undef, N)
    Y = Array{Taylor1{S}}(undef, N)
    Z = Array{Taylor1{S}}(undef, N)

    r_p2 = Array{Taylor1{S}}(undef, N)
    r_p3d2 = Array{Taylor1{S}}(undef, N)
    r_p7d2 = Array{Taylor1{S}}(undef, N)

    newtonianCoeff = Array{Taylor1{S}}(undef, N)

    ui = Array{Taylor1{S}}(undef, N-1)
    vi = Array{Taylor1{S}}(undef, N-1)
    wi = Array{Taylor1{S}}(undef, N-1)

    #post-Newtonian stuff
    U = Array{Taylor1{S}}(undef, N)
    V = Array{Taylor1{S}}(undef, N)
    W = Array{Taylor1{S}}(undef, N)

    _4U_m_3X = Array{Taylor1{S}}(undef, N)
    _4V_m_3Y = Array{Taylor1{S}}(undef, N)
    _4W_m_3Z = Array{Taylor1{S}}(undef, N)

    UU = Array{Taylor1{S}}(undef, N)
    VV = Array{Taylor1{S}}(undef, N)
    WW = Array{Taylor1{S}}(undef, N)

    r_p1d2 = Array{Taylor1{S}}(undef, N)

    newtonianNb_Potential = Array{Taylor1{S}}(undef, N)
    newtonian1b_Potential = Array{Taylor1{S}}(undef, N)
    newtonianCoeff = Array{Taylor1{S}}(undef, N)
    newton_acc_X = Array{Taylor1{S}}(undef, N)
    newton_acc_Y = Array{Taylor1{S}}(undef, N)
    newton_acc_Z = Array{Taylor1{S}}(undef, N)

    v2 = Array{Taylor1{S}}(undef, N)
    vi_dot_vj = Array{Taylor1{S}}(undef, N)
    pn2 = Array{Taylor1{S}}(undef, N)
    pn3 = Array{Taylor1{S}}(undef, N)
    _4ϕj = Array{Taylor1{S}}(undef, N)
    ϕi_plus_4ϕj = Array{Taylor1{S}}(undef, N)
    sj2_plus_2si2_minus_4vivj = Array{Taylor1{S}}(undef, N)
    ϕs_and_vs = Array{Taylor1{S}}(undef, N)
    U_t_pn2 = Array{Taylor1{S}}(undef, N)
    V_t_pn2 = Array{Taylor1{S}}(undef, N)
    W_t_pn2 = Array{Taylor1{S}}(undef, N)
    pn1t1_7 = Array{Taylor1{S}}(undef, N)

    pntempX = zero_q_1
    pntempY = zero_q_1
    pntempZ = zero_q_1
    pn1 = Array{Taylor1{S}}(undef, N)
    X_t_pn1 = Array{Taylor1{S}}(undef, N)
    Y_t_pn1 = Array{Taylor1{S}}(undef, N)
    Z_t_pn1 = Array{Taylor1{S}}(undef, N)
    pNX_t_pn3 = Array{Taylor1{S}}(undef, N)
    pNY_t_pn3 = Array{Taylor1{S}}(undef, N)
    pNZ_t_pn3 = Array{Taylor1{S}}(undef, N)
    pNX_t_X = Array{Taylor1{S}}(undef, N)
    pNY_t_Y = Array{Taylor1{S}}(undef, N)
    pNZ_t_Z = Array{Taylor1{S}}(undef, N)

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
    # temp_accX_j = Array{Taylor1{S}}(undef, N)
    # temp_accY_j = Array{Taylor1{S}}(undef, N)
    # temp_accZ_j = Array{Taylor1{S}}(undef, N)
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

    # extended-body accelerations
    accX = zero_q_1
    accY = zero_q_1
    accZ = zero_q_1

    # rotations to and from Earth, Sun and Moon pole-oriented frames
    local dsj2k = t+(jd0-2.451545e6) # days since J2000.0 = 2.451545e6
    local M_ = Array{Taylor1{S}}(undef, 3, 3, N)
    local M_[:,:,ea] = t2c_jpl_de430(dsj2k)

    dq[1] = q[4]
    dq[2] = q[5]
    dq[3] = q[6]

    newtonianNb_Potential[N] = zero_q_1

    #compute point-mass Newtonian accelerations, all bodies
    Threads.@threads for i in _1_to_Nm1
        ui[i] = ss16asteph_t[3(N-1+i)-2]
        vi[i] = ss16asteph_t[3(N-1+i)-1]
        wi[i] = ss16asteph_t[3(N-1+i)  ]

        X[i] = ss16asteph_t[3i-2]-q[1]
        Y[i] = ss16asteph_t[3i-1]-q[2]
        Z[i] = ss16asteph_t[3i  ]-q[3]

        U[i] = ui[i]-dq[1]
        V[i] = vi[i]-dq[2]
        W[i] = wi[i]-dq[3]

        _4U_m_3X[i] = (4dq[1]) - (3ui[i])
        _4V_m_3Y[i] = (4dq[2]) - (3vi[i])
        _4W_m_3Z[i] = (4dq[3]) - (3wi[i])

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

        newton_acc_X[i] = X[i]*newtonianCoeff[i]
        newton_acc_Y[i] = Y[i]*newtonianCoeff[i]
        newton_acc_Z[i] = Z[i]*newtonianCoeff[i]

        newtonian1b_Potential[i] = μ[i]/r_p1d2[i]
        pn3[i] = 3.5newtonian1b_Potential[i]
        U_t_pn2[i] = pn2[i]*U[i]
        V_t_pn2[i] = pn2[i]*V[i]
        W_t_pn2[i] = pn2[i]*W[i]

        #J2 accelerations, if i-th body is flattened
        if UJ_interaction[i]
            # # rotate from inertial frame to extended-body frame
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
            F_J_ξ[i] = F_J2_ξ[i] #+ F_J3_ξ[i]
            #F_J_η[i] = zero_q_1
            F_J_ζ[i] = F_J2_ζ[i] #+ F_J3_ζ[i]
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
        end # if UJ_interaction[i]
        v2[i] = ( (ui[i]^2)+(vi[i]^2) ) + (wi[i]^2)
    end #for, i
    v2[N] = ( (q[4]^2)+(q[5]^2) ) + (q[6]^2)

    for i in _1_to_Nm1
        temp_004 = newtonian1b_Potential[i] + newtonianNb_Potential[N]
        newtonianNb_Potential[N] = temp_004
        if UJ_interaction[i]
            # # # add result to total acceleration on upon j-th body figure due to i-th point mass
            # # @show "acc",j,"+μ",i,"Λ2",j
            # temp_accX_j[i,j] = accX[j] + (μ[i]*F_J2_x[i,j])
            # accX[j] = temp_accX_j[i,j]
            # temp_accY_j[i,j] = accY[j] + (μ[i]*F_J2_y[i,j])
            # accY[j] = temp_accY_j[i,j]
            # temp_accZ_j[i,j] = accZ[j] + (μ[i]*F_J2_z[i,j])
            # accZ[j] = temp_accZ_j[i,j]

            # # reaction force on i-th body
            # @show "acc",i,"-μ",j,"Λ2",j
            temp_accX_i[i] = accX - (μ[i]*F_J2_x[i])
            accX = temp_accX_i[i]
            temp_accY_i[i] = accY - (μ[i]*F_J2_y[i])
            accY = temp_accY_i[i]
            temp_accZ_i[i] = accZ - (μ[i]*F_J2_z[i])
            accZ = temp_accZ_i[i]
        end
    end

    #post-Newtonian accelerations due to Sun, Moon and planets (Mercury through Neptune)
    #Moyer, 1971, page 7 eq. 35
    # post-Newtonian iterative procedure setup and initialization
    _4ϕj[N] = 4newtonianNb_Potential[N]
    Threads.@threads for i in 1:10
        ϕi_plus_4ϕj[i] = newtonianNb_Potential_t[i] + _4ϕj[N]
        sj2_plus_2si2_minus_4vivj[i] = ( (2v2[i]) - (4vi_dot_vj[i]) ) + v2[N]
        ϕs_and_vs[i] = sj2_plus_2si2_minus_4vivj[i] - ϕi_plus_4ϕj[i]
        Xij_t_Ui = X[i]*ui[i]
        Yij_t_Vi = Y[i]*vi[i]
        Zij_t_Wi = Z[i]*wi[i]
        Rij_dot_Vi = ( Xij_t_Ui+Yij_t_Vi ) + Zij_t_Wi
        # the expression below inside the (...)^2 should have a minus sign in front of the numerator,
        # but upon squaring it is eliminated, so at the end of the day, it is irrelevant ;)
        pn1t7 = (Rij_dot_Vi^2)/r_p2[i]
        pn1t2_7 = ϕs_and_vs[i] - (1.5pn1t7)
        pn1t1_7[i] = c_p2 + pn1t2_7

        pNX_t_X[i] = acceph_t[3i-2]*X[i]
        pNY_t_Y[i] = acceph_t[3i-1]*Y[i]
        pNZ_t_Z[i] = acceph_t[3i  ]*Z[i]
        pn1[i] = (  pn1t1_7[i]  +  (0.5*( (pNX_t_X[i]+pNY_t_Y[i]) + pNZ_t_Z[i] ))  )

        X_t_pn1[i] = newton_acc_X[i]*pn1[i]
        Y_t_pn1[i] = newton_acc_Y[i]*pn1[i]
        Z_t_pn1[i] = newton_acc_Z[i]*pn1[i]

        pNX_t_pn3[i] = acceph_t[3i-2]*pn3[i]
        pNY_t_pn3[i] = acceph_t[3i-1]*pn3[i]
        pNZ_t_pn3[i] = acceph_t[3i  ]*pn3[i]
    end #for i
    for i in 1:10
        termpnx = ( X_t_pn1[i] + (U_t_pn2[i]+pNX_t_pn3[i]) )
        sumpnx = pntempX + termpnx
        pntempX = sumpnx
        termpny = ( Y_t_pn1[i] + (V_t_pn2[i]+pNY_t_pn3[i]) )
        sumpny = pntempY + termpny
        pntempY = sumpny
        termpnz = ( Z_t_pn1[i] + (W_t_pn2[i]+pNZ_t_pn3[i]) )
        sumpnz = pntempZ + termpnz
        pntempZ = sumpnz
    end
    # compute Newtonian accelerations due to Pluto and 16 asteroid perturbers
    Threads.@threads for i in 11:27
        X_t_pn1[i] = c_p2*newton_acc_X[i]
        Y_t_pn1[i] = c_p2*newton_acc_Y[i]
        Z_t_pn1[i] = c_p2*newton_acc_Z[i]
    end #for i
    for i in 11:27
        termpnx = X_t_pn1[i]
        sumpnx = pntempX + termpnx
        pntempX = sumpnx
        termpny = Y_t_pn1[i]
        sumpny = pntempY + termpny
        pntempY = sumpny
        termpnz = Z_t_pn1[i]
        sumpnz = pntempZ + termpnz
        pntempZ = sumpnz
    end
    postNewtonX = pntempX*c_m2
    postNewtonY = pntempY*c_m2
    postNewtonZ = pntempZ*c_m2

    # compute non-gravitational acceleration
    hx = (Y[1]*W[1])-(Z[1]*V[1])
    hy = (Z[1]*U[1])-(X[1]*W[1])
    hz = (X[1]*V[1])-(Y[1]*U[1])

    #cartesian components of transversal unit vector:
    tunitx0 = (hz*Y[1]) - (hy*Z[1]) # Note: Y[1] = y_Sun - y_Apophis, etc.
    tunity0 = (hx*Z[1]) - (hz*X[1])
    tunitz0 = (hy*X[1]) - (hx*Y[1])
    hmag = sqrt( ((tunitx0^2)+(tunity0^2))+(tunitz0^2) )
    tunitx = tunitx0/hmag
    tunity = tunity0/hmag
    tunitz = tunitz0/hmag

    # evaluate non-grav acceleration of NEA (Yarkovsky):
    g_r = r_p2[1]
    A2_t_g_r = q[7]/g_r

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
    local jd0 = params[4]
    local ss16asteph_t = evaleph(params[1], t, q[1])
    local acceph_t = evaleph(params[2], t, q[1])
    local newtonianNb_Potential_t = evaleph(params[3], t, q[1])
    local S = eltype(q[1])
    local N = length(μ)
    local _1_to_Nm1 = 1:N - 1
    local c_p2 = 29979.063823897606
    local c_m2 = 3.3356611996764786e-5
    local zero_q_1 = zero(q[1])
    X = Array{Taylor1{S}}(undef, N)
    Y = Array{Taylor1{S}}(undef, N)
    Z = Array{Taylor1{S}}(undef, N)
    r_p2 = Array{Taylor1{S}}(undef, N)
    r_p3d2 = Array{Taylor1{S}}(undef, N)
    r_p7d2 = Array{Taylor1{S}}(undef, N)
    newtonianCoeff = Array{Taylor1{S}}(undef, N)
    ui = Array{Taylor1{S}}(undef, N - 1)
    vi = Array{Taylor1{S}}(undef, N - 1)
    wi = Array{Taylor1{S}}(undef, N - 1)
    U = Array{Taylor1{S}}(undef, N)
    V = Array{Taylor1{S}}(undef, N)
    W = Array{Taylor1{S}}(undef, N)
    _4U_m_3X = Array{Taylor1{S}}(undef, N)
    _4V_m_3Y = Array{Taylor1{S}}(undef, N)
    _4W_m_3Z = Array{Taylor1{S}}(undef, N)
    UU = Array{Taylor1{S}}(undef, N)
    VV = Array{Taylor1{S}}(undef, N)
    WW = Array{Taylor1{S}}(undef, N)
    r_p1d2 = Array{Taylor1{S}}(undef, N)
    newtonianNb_Potential = Array{Taylor1{S}}(undef, N)
    newtonian1b_Potential = Array{Taylor1{S}}(undef, N)
    newton_acc_X = Array{Taylor1{S}}(undef, N)
    newton_acc_Y = Array{Taylor1{S}}(undef, N)
    newton_acc_Z = Array{Taylor1{S}}(undef, N)
    v2 = Array{Taylor1{S}}(undef, N)
    vi_dot_vj = Array{Taylor1{S}}(undef, N)
    pn2 = Array{Taylor1{S}}(undef, N)
    pn3 = Array{Taylor1{S}}(undef, N)
    _4ϕj = Array{Taylor1{S}}(undef, N)
    ϕi_plus_4ϕj = Array{Taylor1{S}}(undef, N)
    sj2_plus_2si2_minus_4vivj = Array{Taylor1{S}}(undef, N)
    ϕs_and_vs = Array{Taylor1{S}}(undef, N)
    U_t_pn2 = Array{Taylor1{S}}(undef, N)
    V_t_pn2 = Array{Taylor1{S}}(undef, N)
    W_t_pn2 = Array{Taylor1{S}}(undef, N)
    pn1t1_7 = Array{Taylor1{S}}(undef, N)
    pntempX = Taylor1(identity(constant_term(zero_q_1)), order)
    pntempY = Taylor1(identity(constant_term(zero_q_1)), order)
    pntempZ = Taylor1(identity(constant_term(zero_q_1)), order)
    pn1 = Array{Taylor1{S}}(undef, N)
    X_t_pn1 = Array{Taylor1{S}}(undef, N)
    Y_t_pn1 = Array{Taylor1{S}}(undef, N)
    Z_t_pn1 = Array{Taylor1{S}}(undef, N)
    pNX_t_pn3 = Array{Taylor1{S}}(undef, N)
    pNY_t_pn3 = Array{Taylor1{S}}(undef, N)
    pNZ_t_pn3 = Array{Taylor1{S}}(undef, N)
    pNX_t_X = Array{Taylor1{S}}(undef, N)
    pNY_t_Y = Array{Taylor1{S}}(undef, N)
    pNZ_t_Z = Array{Taylor1{S}}(undef, N)
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
    local dsj2k = t + (jd0 - 2.451545e6)
    local M_ = Array{Taylor1{S}}(undef, 3, 3, N)
    local M_[:, :, ea] = t2c_jpl_de430(dsj2k)
    dq[1] = Taylor1(identity(constant_term(q[4])), order)
    dq[2] = Taylor1(identity(constant_term(q[5])), order)
    dq[3] = Taylor1(identity(constant_term(q[6])), order)
    newtonianNb_Potential[N] = Taylor1(identity(constant_term(zero_q_1)), order)
    tmp1273 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1273 .= Taylor1(zero(_S), order)
    tmp1275 = Array{Taylor1{_S}}(undef, size(ui))
    tmp1275 .= Taylor1(zero(_S), order)
    tmp1278 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1278 .= Taylor1(zero(_S), order)
    tmp1280 = Array{Taylor1{_S}}(undef, size(vi))
    tmp1280 .= Taylor1(zero(_S), order)
    tmp1283 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1283 .= Taylor1(zero(_S), order)
    tmp1285 = Array{Taylor1{_S}}(undef, size(wi))
    tmp1285 .= Taylor1(zero(_S), order)
    pn2x = Array{Taylor1{_S}}(undef, size(X))
    pn2x .= Taylor1(zero(_S), order)
    pn2y = Array{Taylor1{_S}}(undef, size(Y))
    pn2y .= Taylor1(zero(_S), order)
    pn2z = Array{Taylor1{_S}}(undef, size(Z))
    pn2z .= Taylor1(zero(_S), order)
    tmp1293 = Array{Taylor1{_S}}(undef, size(UU))
    tmp1293 .= Taylor1(zero(_S), order)
    tmp1296 = Array{Taylor1{_S}}(undef, size(X))
    tmp1296 .= Taylor1(zero(_S), order)
    tmp1298 = Array{Taylor1{_S}}(undef, size(Y))
    tmp1298 .= Taylor1(zero(_S), order)
    tmp1299 = Array{Taylor1{_S}}(undef, size(tmp1296))
    tmp1299 .= Taylor1(zero(_S), order)
    tmp1301 = Array{Taylor1{_S}}(undef, size(Z))
    tmp1301 .= Taylor1(zero(_S), order)
    tmp1309 = Array{Taylor1{_S}}(undef, size(pn2x))
    tmp1309 .= Taylor1(zero(_S), order)
    tmp1310 = Array{Taylor1{_S}}(undef, size(tmp1309))
    tmp1310 .= Taylor1(zero(_S), order)
    tmp1321 = Array{Taylor1{_S}}(undef, size(X))
    tmp1321 .= Taylor1(zero(_S), order)
    tmp1323 = Array{Taylor1{_S}}(undef, size(Y))
    tmp1323 .= Taylor1(zero(_S), order)
    tmp1325 = Array{Taylor1{_S}}(undef, size(Z))
    tmp1325 .= Taylor1(zero(_S), order)
    tmp1327 = Array{Taylor1{_S}}(undef, size(t31))
    tmp1327 .= Taylor1(zero(_S), order)
    tmp1521 = Array{Taylor1{_S}}(undef, size(sin_ϕ))
    tmp1521 .= Taylor1(zero(_S), order)
    tmp1522 = Array{Taylor1{_S}}(undef, size(ϕ))
    tmp1522 .= Taylor1(zero(_S), order)
    tmp1337 = Array{Taylor1{_S}}(undef, size(sin2_ϕ))
    tmp1337 .= Taylor1(zero(_S), order)
    tmp1343 = Array{Taylor1{_S}}(undef, size(sin_ϕ))
    tmp1343 .= Taylor1(zero(_S), order)
    tmp1345 = Array{Taylor1{_S}}(undef, size(sin3_ϕ))
    tmp1345 .= Taylor1(zero(_S), order)
    tmp1349 = Array{Taylor1{_S}}(undef, size(sin2_ϕ))
    tmp1349 .= Taylor1(zero(_S), order)
    tmp1351 = Array{Taylor1{_S}}(undef, size(Λ2))
    tmp1351 .= Taylor1(zero(_S), order)
    tmp1353 = Array{Taylor1{_S}}(undef, size(r_p2))
    tmp1353 .= Taylor1(zero(_S), order)
    tmp1355 = Array{Taylor1{_S}}(undef, size(Λ3))
    tmp1355 .= Taylor1(zero(_S), order)
    tmp1357 = Array{Taylor1{_S}}(undef, size(r_p1d2))
    tmp1357 .= Taylor1(zero(_S), order)
    tmp1359 = Array{Taylor1{_S}}(undef, size(cos_ϕ))
    tmp1359 .= Taylor1(zero(_S), order)
    tmp1361 = Array{Taylor1{_S}}(undef, size(cos_ϕ))
    tmp1361 .= Taylor1(zero(_S), order)
    tmp1364 = Array{Taylor1{_S}}(undef, size(Λ2j_div_r4))
    tmp1364 .= Taylor1(zero(_S), order)
    tmp1368 = Array{Taylor1{_S}}(undef, size(Λ3j_div_r5))
    tmp1368 .= Taylor1(zero(_S), order)
    tmp1371 = Array{Taylor1{_S}}(undef, size(X))
    tmp1371 .= Taylor1(zero(_S), order)
    tmp1373 = Array{Taylor1{_S}}(undef, size(Y))
    tmp1373 .= Taylor1(zero(_S), order)
    tmp1375 = Array{Taylor1{_S}}(undef, size(Z))
    tmp1375 .= Taylor1(zero(_S), order)
    tmp1405 = Array{Taylor1{_S}}(undef, size(ui))
    tmp1405 .= Taylor1(zero(_S), order)
    tmp1407 = Array{Taylor1{_S}}(undef, size(vi))
    tmp1407 .= Taylor1(zero(_S), order)
    tmp1408 = Array{Taylor1{_S}}(undef, size(tmp1405))
    tmp1408 .= Taylor1(zero(_S), order)
    tmp1410 = Array{Taylor1{_S}}(undef, size(wi))
    tmp1410 .= Taylor1(zero(_S), order)
    for i = _1_to_Nm1
        ui[i] = Taylor1(identity(constant_term(ss16asteph_t[3 * ((N - 1) + i) - 2])), order)
        vi[i] = Taylor1(identity(constant_term(ss16asteph_t[3 * ((N - 1) + i) - 1])), order)
        wi[i] = Taylor1(identity(constant_term(ss16asteph_t[3 * ((N - 1) + i)])), order)
        X[i] = Taylor1(constant_term(ss16asteph_t[3i - 2]) - constant_term(q[1]), order)
        Y[i] = Taylor1(constant_term(ss16asteph_t[3i - 1]) - constant_term(q[2]), order)
        Z[i] = Taylor1(constant_term(ss16asteph_t[3i]) - constant_term(q[3]), order)
        U[i] = Taylor1(constant_term(ui[i]) - constant_term(dq[1]), order)
        V[i] = Taylor1(constant_term(vi[i]) - constant_term(dq[2]), order)
        W[i] = Taylor1(constant_term(wi[i]) - constant_term(dq[3]), order)
        tmp1273[1] = Taylor1(constant_term(4) * constant_term(dq[1]), order)
        tmp1275[i] = Taylor1(constant_term(3) * constant_term(ui[i]), order)
        _4U_m_3X[i] = Taylor1(constant_term(tmp1273[1]) - constant_term(tmp1275[i]), order)
        tmp1278[2] = Taylor1(constant_term(4) * constant_term(dq[2]), order)
        tmp1280[i] = Taylor1(constant_term(3) * constant_term(vi[i]), order)
        _4V_m_3Y[i] = Taylor1(constant_term(tmp1278[2]) - constant_term(tmp1280[i]), order)
        tmp1283[3] = Taylor1(constant_term(4) * constant_term(dq[3]), order)
        tmp1285[i] = Taylor1(constant_term(3) * constant_term(wi[i]), order)
        _4W_m_3Z[i] = Taylor1(constant_term(tmp1283[3]) - constant_term(tmp1285[i]), order)
        pn2x[i] = Taylor1(constant_term(X[i]) * constant_term(_4U_m_3X[i]), order)
        pn2y[i] = Taylor1(constant_term(Y[i]) * constant_term(_4V_m_3Y[i]), order)
        pn2z[i] = Taylor1(constant_term(Z[i]) * constant_term(_4W_m_3Z[i]), order)
        UU[i] = Taylor1(constant_term(ui[i]) * constant_term(dq[1]), order)
        VV[i] = Taylor1(constant_term(vi[i]) * constant_term(dq[2]), order)
        WW[i] = Taylor1(constant_term(wi[i]) * constant_term(dq[3]), order)
        tmp1293[i] = Taylor1(constant_term(UU[i]) + constant_term(VV[i]), order)
        vi_dot_vj[i] = Taylor1(constant_term(tmp1293[i]) + constant_term(WW[i]), order)
        tmp1296[i] = Taylor1(constant_term(X[i]) ^ float(constant_term(2)), order)
        tmp1298[i] = Taylor1(constant_term(Y[i]) ^ float(constant_term(2)), order)
        tmp1299[i] = Taylor1(constant_term(tmp1296[i]) + constant_term(tmp1298[i]), order)
        tmp1301[i] = Taylor1(constant_term(Z[i]) ^ float(constant_term(2)), order)
        r_p2[i] = Taylor1(constant_term(tmp1299[i]) + constant_term(tmp1301[i]), order)
        r_p1d2[i] = Taylor1(sqrt(constant_term(r_p2[i])), order)
        r_p3d2[i] = Taylor1(constant_term(r_p2[i]) ^ float(constant_term(1.5)), order)
        r_p7d2[i] = Taylor1(constant_term(r_p2[i]) ^ float(constant_term(3.5)), order)
        newtonianCoeff[i] = Taylor1(constant_term(μ[i]) / constant_term(r_p3d2[i]), order)
        tmp1309[i] = Taylor1(constant_term(pn2x[i]) + constant_term(pn2y[i]), order)
        tmp1310[i] = Taylor1(constant_term(tmp1309[i]) + constant_term(pn2z[i]), order)
        pn2[i] = Taylor1(constant_term(newtonianCoeff[i]) * constant_term(tmp1310[i]), order)
        newton_acc_X[i] = Taylor1(constant_term(X[i]) * constant_term(newtonianCoeff[i]), order)
        newton_acc_Y[i] = Taylor1(constant_term(Y[i]) * constant_term(newtonianCoeff[i]), order)
        newton_acc_Z[i] = Taylor1(constant_term(Z[i]) * constant_term(newtonianCoeff[i]), order)
        newtonian1b_Potential[i] = Taylor1(constant_term(μ[i]) / constant_term(r_p1d2[i]), order)
        pn3[i] = Taylor1(constant_term(3.5) * constant_term(newtonian1b_Potential[i]), order)
        U_t_pn2[i] = Taylor1(constant_term(pn2[i]) * constant_term(U[i]), order)
        V_t_pn2[i] = Taylor1(constant_term(pn2[i]) * constant_term(V[i]), order)
        W_t_pn2[i] = Taylor1(constant_term(pn2[i]) * constant_term(W[i]), order)
        if UJ_interaction[i]
            tmp1321[i] = Taylor1(-(constant_term(X[i])), order)
            t31[i] = Taylor1(constant_term(tmp1321[i]) * constant_term(M_[1, 3, i]), order)
            tmp1323[i] = Taylor1(-(constant_term(Y[i])), order)
            t32[i] = Taylor1(constant_term(tmp1323[i]) * constant_term(M_[2, 3, i]), order)
            tmp1325[i] = Taylor1(-(constant_term(Z[i])), order)
            t33[i] = Taylor1(constant_term(tmp1325[i]) * constant_term(M_[3, 3, i]), order)
            tmp1327[i] = Taylor1(constant_term(t31[i]) + constant_term(t32[i]), order)
            r_sin_ϕ[i] = Taylor1(constant_term(tmp1327[i]) + constant_term(t33[i]), order)
            sin_ϕ[i] = Taylor1(constant_term(r_sin_ϕ[i]) / constant_term(r_p1d2[i]), order)
            ϕ[i] = Taylor1(asin(constant_term(sin_ϕ[i])), order)
            tmp1521[i] = Taylor1(sqrt(1 - constant_term(sin_ϕ[i]) ^ 2), order)
            cos_ϕ[i] = Taylor1(cos(constant_term(ϕ[i])), order)
            tmp1522[i] = Taylor1(sin(constant_term(ϕ[i])), order)
            sin2_ϕ[i] = Taylor1(constant_term(sin_ϕ[i]) ^ float(constant_term(2)), order)
            sin3_ϕ[i] = Taylor1(constant_term(sin_ϕ[i]) ^ float(constant_term(3)), order)
            tmp1337[i] = Taylor1(constant_term(1.5) * constant_term(sin2_ϕ[i]), order)
            P_2_sin_ϕ[i] = Taylor1(constant_term(tmp1337[i]) - constant_term(0.5), order)
            ∂P_2_sin_ϕ[i] = Taylor1(constant_term(3) * constant_term(sin_ϕ[i]), order)
            tmp1343[i] = Taylor1(constant_term(-1.5) * constant_term(sin_ϕ[i]), order)
            tmp1345[i] = Taylor1(constant_term(2.5) * constant_term(sin3_ϕ[i]), order)
            P_3_sin_ϕ[i] = Taylor1(constant_term(tmp1343[i]) + constant_term(tmp1345[i]), order)
            tmp1349[i] = Taylor1(constant_term(7.5) * constant_term(sin2_ϕ[i]), order)
            ∂P_3_sin_ϕ[i] = Taylor1(constant_term(-1.5) + constant_term(tmp1349[i]), order)
            tmp1351[i] = Taylor1(-(constant_term(Λ2[i])), order)
            tmp1353[i] = Taylor1(constant_term(r_p2[i]) ^ float(constant_term(2)), order)
            Λ2j_div_r4[i] = Taylor1(constant_term(tmp1351[i]) / constant_term(tmp1353[i]), order)
            tmp1355[i] = Taylor1(-(constant_term(Λ3[i])), order)
            tmp1357[i] = Taylor1(constant_term(r_p1d2[i]) ^ float(constant_term(5)), order)
            Λ3j_div_r5[i] = Taylor1(constant_term(tmp1355[i]) / constant_term(tmp1357[i]), order)
            tmp1359[i] = Taylor1(-(constant_term(cos_ϕ[i])), order)
            m_c_ϕ_∂P_2[i] = Taylor1(constant_term(tmp1359[i]) * constant_term(∂P_2_sin_ϕ[i]), order)
            tmp1361[i] = Taylor1(-(constant_term(cos_ϕ[i])), order)
            m_c_ϕ_∂P_3[i] = Taylor1(constant_term(tmp1361[i]) * constant_term(∂P_3_sin_ϕ[i]), order)
            tmp1364[i] = Taylor1(constant_term(Λ2j_div_r4[i]) * constant_term(3), order)
            F_J2_ξ[i] = Taylor1(constant_term(tmp1364[i]) * constant_term(P_2_sin_ϕ[i]), order)
            F_J2_ζ[i] = Taylor1(constant_term(Λ2j_div_r4[i]) * constant_term(m_c_ϕ_∂P_2[i]), order)
            tmp1368[i] = Taylor1(constant_term(Λ3j_div_r5[i]) * constant_term(4), order)
            F_J3_ξ[i] = Taylor1(constant_term(tmp1368[i]) * constant_term(P_3_sin_ϕ[i]), order)
            F_J3_ζ[i] = Taylor1(constant_term(Λ3j_div_r5[i]) * constant_term(m_c_ϕ_∂P_3[i]), order)
            F_J_ξ[i] = Taylor1(identity(constant_term(F_J2_ξ[i])), order)
            F_J_ζ[i] = Taylor1(identity(constant_term(F_J2_ζ[i])), order)
            tmp1371[i] = Taylor1(-(constant_term(X[i])), order)
            ξx[i] = Taylor1(constant_term(tmp1371[i]) / constant_term(r_p1d2[i]), order)
            tmp1373[i] = Taylor1(-(constant_term(Y[i])), order)
            ξy[i] = Taylor1(constant_term(tmp1373[i]) / constant_term(r_p1d2[i]), order)
            tmp1375[i] = Taylor1(-(constant_term(Z[i])), order)
            ξz[i] = Taylor1(constant_term(tmp1375[i]) / constant_term(r_p1d2[i]), order)
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
        end
        tmp1405[i] = Taylor1(constant_term(ui[i]) ^ float(constant_term(2)), order)
        tmp1407[i] = Taylor1(constant_term(vi[i]) ^ float(constant_term(2)), order)
        tmp1408[i] = Taylor1(constant_term(tmp1405[i]) + constant_term(tmp1407[i]), order)
        tmp1410[i] = Taylor1(constant_term(wi[i]) ^ float(constant_term(2)), order)
        v2[i] = Taylor1(constant_term(tmp1408[i]) + constant_term(tmp1410[i]), order)
    end
    tmp1413 = Taylor1(constant_term(q[4]) ^ float(constant_term(2)), order)
    tmp1415 = Taylor1(constant_term(q[5]) ^ float(constant_term(2)), order)
    tmp1416 = Taylor1(constant_term(tmp1413) + constant_term(tmp1415), order)
    tmp1418 = Taylor1(constant_term(q[6]) ^ float(constant_term(2)), order)
    v2[N] = Taylor1(constant_term(tmp1416) + constant_term(tmp1418), order)
    temp_004 = Array{Taylor1{_S}}(undef, size(newtonian1b_Potential))
    temp_004 .= Taylor1(zero(_S), order)
    tmp1421 = Array{Taylor1{_S}}(undef, size(μ))
    tmp1421 .= Taylor1(zero(_S), order)
    tmp1423 = Array{Taylor1{_S}}(undef, size(μ))
    tmp1423 .= Taylor1(zero(_S), order)
    tmp1425 = Array{Taylor1{_S}}(undef, size(μ))
    tmp1425 .= Taylor1(zero(_S), order)
    for i = _1_to_Nm1
        temp_004[i] = Taylor1(constant_term(newtonian1b_Potential[i]) + constant_term(newtonianNb_Potential[N]), order)
        newtonianNb_Potential[N] = Taylor1(identity(constant_term(temp_004[i])), order)
        if UJ_interaction[i]
            tmp1421[i] = Taylor1(constant_term(μ[i]) * constant_term(F_J2_x[i]), order)
            temp_accX_i[i] = Taylor1(constant_term(accX) - constant_term(tmp1421[i]), order)
            accX = Taylor1(identity(constant_term(temp_accX_i[i])), order)
            tmp1423[i] = Taylor1(constant_term(μ[i]) * constant_term(F_J2_y[i]), order)
            temp_accY_i[i] = Taylor1(constant_term(accY) - constant_term(tmp1423[i]), order)
            accY = Taylor1(identity(constant_term(temp_accY_i[i])), order)
            tmp1425[i] = Taylor1(constant_term(μ[i]) * constant_term(F_J2_z[i]), order)
            temp_accZ_i[i] = Taylor1(constant_term(accZ) - constant_term(tmp1425[i]), order)
            accZ = Taylor1(identity(constant_term(temp_accZ_i[i])), order)
        end
    end
    _4ϕj[N] = Taylor1(constant_term(4) * constant_term(newtonianNb_Potential[N]), order)
    tmp1431 = Array{Taylor1{_S}}(undef, size(v2))
    tmp1431 .= Taylor1(zero(_S), order)
    tmp1433 = Array{Taylor1{_S}}(undef, size(vi_dot_vj))
    tmp1433 .= Taylor1(zero(_S), order)
    tmp1434 = Array{Taylor1{_S}}(undef, size(tmp1431))
    tmp1434 .= Taylor1(zero(_S), order)
    Xij_t_Ui = Array{Taylor1{_S}}(undef, size(X))
    Xij_t_Ui .= Taylor1(zero(_S), order)
    Yij_t_Vi = Array{Taylor1{_S}}(undef, size(Y))
    Yij_t_Vi .= Taylor1(zero(_S), order)
    Zij_t_Wi = Array{Taylor1{_S}}(undef, size(Z))
    Zij_t_Wi .= Taylor1(zero(_S), order)
    tmp1440 = Array{Taylor1{_S}}(undef, size(Xij_t_Ui))
    tmp1440 .= Taylor1(zero(_S), order)
    Rij_dot_Vi = Array{Taylor1{_S}}(undef, size(tmp1440))
    Rij_dot_Vi .= Taylor1(zero(_S), order)
    tmp1443 = Array{Taylor1{_S}}(undef, size(Rij_dot_Vi))
    tmp1443 .= Taylor1(zero(_S), order)
    pn1t7 = Array{Taylor1{_S}}(undef, size(tmp1443))
    pn1t7 .= Taylor1(zero(_S), order)
    tmp1446 = Array{Taylor1{_S}}(undef, size(pn1t7))
    tmp1446 .= Taylor1(zero(_S), order)
    pn1t2_7 = Array{Taylor1{_S}}(undef, size(ϕs_and_vs))
    pn1t2_7 .= Taylor1(zero(_S), order)
    tmp1453 = Array{Taylor1{_S}}(undef, size(pNX_t_X))
    tmp1453 .= Taylor1(zero(_S), order)
    tmp1454 = Array{Taylor1{_S}}(undef, size(tmp1453))
    tmp1454 .= Taylor1(zero(_S), order)
    tmp1455 = Array{Taylor1{_S}}(undef, size(tmp1454))
    tmp1455 .= Taylor1(zero(_S), order)
    for i = 1:10
        ϕi_plus_4ϕj[i] = Taylor1(constant_term(newtonianNb_Potential_t[i]) + constant_term(_4ϕj[N]), order)
        tmp1431[i] = Taylor1(constant_term(2) * constant_term(v2[i]), order)
        tmp1433[i] = Taylor1(constant_term(4) * constant_term(vi_dot_vj[i]), order)
        tmp1434[i] = Taylor1(constant_term(tmp1431[i]) - constant_term(tmp1433[i]), order)
        sj2_plus_2si2_minus_4vivj[i] = Taylor1(constant_term(tmp1434[i]) + constant_term(v2[N]), order)
        ϕs_and_vs[i] = Taylor1(constant_term(sj2_plus_2si2_minus_4vivj[i]) - constant_term(ϕi_plus_4ϕj[i]), order)
        Xij_t_Ui[i] = Taylor1(constant_term(X[i]) * constant_term(ui[i]), order)
        Yij_t_Vi[i] = Taylor1(constant_term(Y[i]) * constant_term(vi[i]), order)
        Zij_t_Wi[i] = Taylor1(constant_term(Z[i]) * constant_term(wi[i]), order)
        tmp1440[i] = Taylor1(constant_term(Xij_t_Ui[i]) + constant_term(Yij_t_Vi[i]), order)
        Rij_dot_Vi[i] = Taylor1(constant_term(tmp1440[i]) + constant_term(Zij_t_Wi[i]), order)
        tmp1443[i] = Taylor1(constant_term(Rij_dot_Vi[i]) ^ float(constant_term(2)), order)
        pn1t7[i] = Taylor1(constant_term(tmp1443[i]) / constant_term(r_p2[i]), order)
        tmp1446[i] = Taylor1(constant_term(1.5) * constant_term(pn1t7[i]), order)
        pn1t2_7[i] = Taylor1(constant_term(ϕs_and_vs[i]) - constant_term(tmp1446[i]), order)
        pn1t1_7[i] = Taylor1(constant_term(c_p2) + constant_term(pn1t2_7[i]), order)
        pNX_t_X[i] = Taylor1(constant_term(acceph_t[3i - 2]) * constant_term(X[i]), order)
        pNY_t_Y[i] = Taylor1(constant_term(acceph_t[3i - 1]) * constant_term(Y[i]), order)
        pNZ_t_Z[i] = Taylor1(constant_term(acceph_t[3i]) * constant_term(Z[i]), order)
        tmp1453[i] = Taylor1(constant_term(pNX_t_X[i]) + constant_term(pNY_t_Y[i]), order)
        tmp1454[i] = Taylor1(constant_term(tmp1453[i]) + constant_term(pNZ_t_Z[i]), order)
        tmp1455[i] = Taylor1(constant_term(0.5) * constant_term(tmp1454[i]), order)
        pn1[i] = Taylor1(constant_term(pn1t1_7[i]) + constant_term(tmp1455[i]), order)
        X_t_pn1[i] = Taylor1(constant_term(newton_acc_X[i]) * constant_term(pn1[i]), order)
        Y_t_pn1[i] = Taylor1(constant_term(newton_acc_Y[i]) * constant_term(pn1[i]), order)
        Z_t_pn1[i] = Taylor1(constant_term(newton_acc_Z[i]) * constant_term(pn1[i]), order)
        pNX_t_pn3[i] = Taylor1(constant_term(acceph_t[3i - 2]) * constant_term(pn3[i]), order)
        pNY_t_pn3[i] = Taylor1(constant_term(acceph_t[3i - 1]) * constant_term(pn3[i]), order)
        pNZ_t_pn3[i] = Taylor1(constant_term(acceph_t[3i]) * constant_term(pn3[i]), order)
    end
    tmp1463 = Array{Taylor1{_S}}(undef, size(U_t_pn2))
    tmp1463 .= Taylor1(zero(_S), order)
    termpnx = Array{Taylor1{_S}}(undef, size(X_t_pn1))
    termpnx .= Taylor1(zero(_S), order)
    sumpnx = Array{Taylor1{_S}}(undef, size(termpnx))
    sumpnx .= Taylor1(zero(_S), order)
    tmp1466 = Array{Taylor1{_S}}(undef, size(V_t_pn2))
    tmp1466 .= Taylor1(zero(_S), order)
    termpny = Array{Taylor1{_S}}(undef, size(Y_t_pn1))
    termpny .= Taylor1(zero(_S), order)
    sumpny = Array{Taylor1{_S}}(undef, size(termpny))
    sumpny .= Taylor1(zero(_S), order)
    tmp1469 = Array{Taylor1{_S}}(undef, size(W_t_pn2))
    tmp1469 .= Taylor1(zero(_S), order)
    termpnz = Array{Taylor1{_S}}(undef, size(Z_t_pn1))
    termpnz .= Taylor1(zero(_S), order)
    sumpnz = Array{Taylor1{_S}}(undef, size(termpnz))
    sumpnz .= Taylor1(zero(_S), order)
    for i = 1:10
        tmp1463[i] = Taylor1(constant_term(U_t_pn2[i]) + constant_term(pNX_t_pn3[i]), order)
        termpnx[i] = Taylor1(constant_term(X_t_pn1[i]) + constant_term(tmp1463[i]), order)
        sumpnx[i] = Taylor1(constant_term(pntempX) + constant_term(termpnx[i]), order)
        pntempX = Taylor1(identity(constant_term(sumpnx[i])), order)
        tmp1466[i] = Taylor1(constant_term(V_t_pn2[i]) + constant_term(pNY_t_pn3[i]), order)
        termpny[i] = Taylor1(constant_term(Y_t_pn1[i]) + constant_term(tmp1466[i]), order)
        sumpny[i] = Taylor1(constant_term(pntempY) + constant_term(termpny[i]), order)
        pntempY = Taylor1(identity(constant_term(sumpny[i])), order)
        tmp1469[i] = Taylor1(constant_term(W_t_pn2[i]) + constant_term(pNZ_t_pn3[i]), order)
        termpnz[i] = Taylor1(constant_term(Z_t_pn1[i]) + constant_term(tmp1469[i]), order)
        sumpnz[i] = Taylor1(constant_term(pntempZ) + constant_term(termpnz[i]), order)
        pntempZ = Taylor1(identity(constant_term(sumpnz[i])), order)
    end
    for i = 11:27
        X_t_pn1[i] = Taylor1(constant_term(c_p2) * constant_term(newton_acc_X[i]), order)
        Y_t_pn1[i] = Taylor1(constant_term(c_p2) * constant_term(newton_acc_Y[i]), order)
        Z_t_pn1[i] = Taylor1(constant_term(c_p2) * constant_term(newton_acc_Z[i]), order)
    end
    for i = 11:27
        termpnx[i] = Taylor1(identity(constant_term(X_t_pn1[i])), order)
        sumpnx[i] = Taylor1(constant_term(pntempX) + constant_term(termpnx[i]), order)
        pntempX = Taylor1(identity(constant_term(sumpnx[i])), order)
        termpny[i] = Taylor1(identity(constant_term(Y_t_pn1[i])), order)
        sumpny[i] = Taylor1(constant_term(pntempY) + constant_term(termpny[i]), order)
        pntempY = Taylor1(identity(constant_term(sumpny[i])), order)
        termpnz[i] = Taylor1(identity(constant_term(Z_t_pn1[i])), order)
        sumpnz[i] = Taylor1(constant_term(pntempZ) + constant_term(termpnz[i]), order)
        pntempZ = Taylor1(identity(constant_term(sumpnz[i])), order)
    end
    postNewtonX = Taylor1(constant_term(pntempX) * constant_term(c_m2), order)
    postNewtonY = Taylor1(constant_term(pntempY) * constant_term(c_m2), order)
    postNewtonZ = Taylor1(constant_term(pntempZ) * constant_term(c_m2), order)
    tmp1481 = Taylor1(constant_term(Y[1]) * constant_term(W[1]), order)
    tmp1482 = Taylor1(constant_term(Z[1]) * constant_term(V[1]), order)
    hx = Taylor1(constant_term(tmp1481) - constant_term(tmp1482), order)
    tmp1484 = Taylor1(constant_term(Z[1]) * constant_term(U[1]), order)
    tmp1485 = Taylor1(constant_term(X[1]) * constant_term(W[1]), order)
    hy = Taylor1(constant_term(tmp1484) - constant_term(tmp1485), order)
    tmp1487 = Taylor1(constant_term(X[1]) * constant_term(V[1]), order)
    tmp1488 = Taylor1(constant_term(Y[1]) * constant_term(U[1]), order)
    hz = Taylor1(constant_term(tmp1487) - constant_term(tmp1488), order)
    tmp1490 = Taylor1(constant_term(hz) * constant_term(Y[1]), order)
    tmp1491 = Taylor1(constant_term(hy) * constant_term(Z[1]), order)
    tunitx0 = Taylor1(constant_term(tmp1490) - constant_term(tmp1491), order)
    tmp1493 = Taylor1(constant_term(hx) * constant_term(Z[1]), order)
    tmp1494 = Taylor1(constant_term(hz) * constant_term(X[1]), order)
    tunity0 = Taylor1(constant_term(tmp1493) - constant_term(tmp1494), order)
    tmp1496 = Taylor1(constant_term(hy) * constant_term(X[1]), order)
    tmp1497 = Taylor1(constant_term(hx) * constant_term(Y[1]), order)
    tunitz0 = Taylor1(constant_term(tmp1496) - constant_term(tmp1497), order)
    tmp1500 = Taylor1(constant_term(tunitx0) ^ float(constant_term(2)), order)
    tmp1502 = Taylor1(constant_term(tunity0) ^ float(constant_term(2)), order)
    tmp1503 = Taylor1(constant_term(tmp1500) + constant_term(tmp1502), order)
    tmp1505 = Taylor1(constant_term(tunitz0) ^ float(constant_term(2)), order)
    tmp1506 = Taylor1(constant_term(tmp1503) + constant_term(tmp1505), order)
    hmag = Taylor1(sqrt(constant_term(tmp1506)), order)
    tunitx = Taylor1(constant_term(tunitx0) / constant_term(hmag), order)
    tunity = Taylor1(constant_term(tunity0) / constant_term(hmag), order)
    tunitz = Taylor1(constant_term(tunitz0) / constant_term(hmag), order)
    g_r = Taylor1(identity(constant_term(r_p2[1])), order)
    A2_t_g_r = Taylor1(constant_term(q[7]) / constant_term(g_r), order)
    NGAx = Taylor1(constant_term(A2_t_g_r) * constant_term(tunitx), order)
    NGAy = Taylor1(constant_term(A2_t_g_r) * constant_term(tunity), order)
    NGAz = Taylor1(constant_term(A2_t_g_r) * constant_term(tunitz), order)
    tmp1515 = Taylor1(constant_term(postNewtonX) + constant_term(accX), order)
    dq[4] = Taylor1(constant_term(tmp1515) + constant_term(NGAx), order)
    tmp1517 = Taylor1(constant_term(postNewtonY) + constant_term(accY), order)
    dq[5] = Taylor1(constant_term(tmp1517) + constant_term(NGAy), order)
    tmp1519 = Taylor1(constant_term(postNewtonZ) + constant_term(accZ), order)
    dq[6] = Taylor1(constant_term(tmp1519) + constant_term(NGAz), order)
    dq[7] = Taylor1(identity(constant_term(zero_q_1)), order)
    for __idx = eachindex(q)
        (q[__idx]).coeffs[2] = (dq[__idx]).coeffs[1]
    end
    for ord = 1:order - 1
        ordnext = ord + 1
        TaylorSeries.identity!(pntempX, zero_q_1, ord)
        TaylorSeries.identity!(pntempY, zero_q_1, ord)
        TaylorSeries.identity!(pntempZ, zero_q_1, ord)
        TaylorSeries.identity!(accX, zero_q_1, ord)
        TaylorSeries.identity!(accY, zero_q_1, ord)
        TaylorSeries.identity!(accZ, zero_q_1, ord)
        TaylorSeries.identity!(dq[1], q[4], ord)
        TaylorSeries.identity!(dq[2], q[5], ord)
        TaylorSeries.identity!(dq[3], q[6], ord)
        TaylorSeries.identity!(newtonianNb_Potential[N], zero_q_1, ord)
        for i = _1_to_Nm1
            TaylorSeries.identity!(ui[i], ss16asteph_t[3 * ((N - 1) + i) - 2], ord)
            TaylorSeries.identity!(vi[i], ss16asteph_t[3 * ((N - 1) + i) - 1], ord)
            TaylorSeries.identity!(wi[i], ss16asteph_t[3 * ((N - 1) + i)], ord)
            TaylorSeries.subst!(X[i], ss16asteph_t[3i - 2], q[1], ord)
            TaylorSeries.subst!(Y[i], ss16asteph_t[3i - 1], q[2], ord)
            TaylorSeries.subst!(Z[i], ss16asteph_t[3i], q[3], ord)
            TaylorSeries.subst!(U[i], ui[i], dq[1], ord)
            TaylorSeries.subst!(V[i], vi[i], dq[2], ord)
            TaylorSeries.subst!(W[i], wi[i], dq[3], ord)
            TaylorSeries.mul!(tmp1273[1], 4, dq[1], ord)
            TaylorSeries.mul!(tmp1275[i], 3, ui[i], ord)
            TaylorSeries.subst!(_4U_m_3X[i], tmp1273[1], tmp1275[i], ord)
            TaylorSeries.mul!(tmp1278[2], 4, dq[2], ord)
            TaylorSeries.mul!(tmp1280[i], 3, vi[i], ord)
            TaylorSeries.subst!(_4V_m_3Y[i], tmp1278[2], tmp1280[i], ord)
            TaylorSeries.mul!(tmp1283[3], 4, dq[3], ord)
            TaylorSeries.mul!(tmp1285[i], 3, wi[i], ord)
            TaylorSeries.subst!(_4W_m_3Z[i], tmp1283[3], tmp1285[i], ord)
            TaylorSeries.mul!(pn2x[i], X[i], _4U_m_3X[i], ord)
            TaylorSeries.mul!(pn2y[i], Y[i], _4V_m_3Y[i], ord)
            TaylorSeries.mul!(pn2z[i], Z[i], _4W_m_3Z[i], ord)
            TaylorSeries.mul!(UU[i], ui[i], dq[1], ord)
            TaylorSeries.mul!(VV[i], vi[i], dq[2], ord)
            TaylorSeries.mul!(WW[i], wi[i], dq[3], ord)
            TaylorSeries.add!(tmp1293[i], UU[i], VV[i], ord)
            TaylorSeries.add!(vi_dot_vj[i], tmp1293[i], WW[i], ord)
            TaylorSeries.pow!(tmp1296[i], X[i], 2, ord)
            TaylorSeries.pow!(tmp1298[i], Y[i], 2, ord)
            TaylorSeries.add!(tmp1299[i], tmp1296[i], tmp1298[i], ord)
            TaylorSeries.pow!(tmp1301[i], Z[i], 2, ord)
            TaylorSeries.add!(r_p2[i], tmp1299[i], tmp1301[i], ord)
            TaylorSeries.sqrt!(r_p1d2[i], r_p2[i], ord)
            TaylorSeries.pow!(r_p3d2[i], r_p2[i], 1.5, ord)
            TaylorSeries.pow!(r_p7d2[i], r_p2[i], 3.5, ord)
            TaylorSeries.div!(newtonianCoeff[i], μ[i], r_p3d2[i], ord)
            TaylorSeries.add!(tmp1309[i], pn2x[i], pn2y[i], ord)
            TaylorSeries.add!(tmp1310[i], tmp1309[i], pn2z[i], ord)
            TaylorSeries.mul!(pn2[i], newtonianCoeff[i], tmp1310[i], ord)
            TaylorSeries.mul!(newton_acc_X[i], X[i], newtonianCoeff[i], ord)
            TaylorSeries.mul!(newton_acc_Y[i], Y[i], newtonianCoeff[i], ord)
            TaylorSeries.mul!(newton_acc_Z[i], Z[i], newtonianCoeff[i], ord)
            TaylorSeries.div!(newtonian1b_Potential[i], μ[i], r_p1d2[i], ord)
            TaylorSeries.mul!(pn3[i], 3.5, newtonian1b_Potential[i], ord)
            TaylorSeries.mul!(U_t_pn2[i], pn2[i], U[i], ord)
            TaylorSeries.mul!(V_t_pn2[i], pn2[i], V[i], ord)
            TaylorSeries.mul!(W_t_pn2[i], pn2[i], W[i], ord)
            if UJ_interaction[i]
                TaylorSeries.subst!(tmp1321[i], X[i], ord)
                TaylorSeries.mul!(t31[i], tmp1321[i], M_[1, 3, i], ord)
                TaylorSeries.subst!(tmp1323[i], Y[i], ord)
                TaylorSeries.mul!(t32[i], tmp1323[i], M_[2, 3, i], ord)
                TaylorSeries.subst!(tmp1325[i], Z[i], ord)
                TaylorSeries.mul!(t33[i], tmp1325[i], M_[3, 3, i], ord)
                TaylorSeries.add!(tmp1327[i], t31[i], t32[i], ord)
                TaylorSeries.add!(r_sin_ϕ[i], tmp1327[i], t33[i], ord)
                TaylorSeries.div!(sin_ϕ[i], r_sin_ϕ[i], r_p1d2[i], ord)
                TaylorSeries.asin!(ϕ[i], sin_ϕ[i], tmp1521[i], ord)
                TaylorSeries.sincos!(tmp1522[i], cos_ϕ[i], ϕ[i], ord)
                TaylorSeries.pow!(sin2_ϕ[i], sin_ϕ[i], 2, ord)
                TaylorSeries.pow!(sin3_ϕ[i], sin_ϕ[i], 3, ord)
                TaylorSeries.mul!(tmp1337[i], 1.5, sin2_ϕ[i], ord)
                TaylorSeries.subst!(P_2_sin_ϕ[i], tmp1337[i], 0.5, ord)
                TaylorSeries.mul!(∂P_2_sin_ϕ[i], 3, sin_ϕ[i], ord)
                TaylorSeries.mul!(tmp1343[i], -1.5, sin_ϕ[i], ord)
                TaylorSeries.mul!(tmp1345[i], 2.5, sin3_ϕ[i], ord)
                TaylorSeries.add!(P_3_sin_ϕ[i], tmp1343[i], tmp1345[i], ord)
                TaylorSeries.mul!(tmp1349[i], 7.5, sin2_ϕ[i], ord)
                TaylorSeries.add!(∂P_3_sin_ϕ[i], -1.5, tmp1349[i], ord)
                TaylorSeries.subst!(tmp1351[i], Λ2[i], ord)
                TaylorSeries.pow!(tmp1353[i], r_p2[i], 2, ord)
                TaylorSeries.div!(Λ2j_div_r4[i], tmp1351[i], tmp1353[i], ord)
                TaylorSeries.subst!(tmp1355[i], Λ3[i], ord)
                TaylorSeries.pow!(tmp1357[i], r_p1d2[i], 5, ord)
                TaylorSeries.div!(Λ3j_div_r5[i], tmp1355[i], tmp1357[i], ord)
                TaylorSeries.subst!(tmp1359[i], cos_ϕ[i], ord)
                TaylorSeries.mul!(m_c_ϕ_∂P_2[i], tmp1359[i], ∂P_2_sin_ϕ[i], ord)
                TaylorSeries.subst!(tmp1361[i], cos_ϕ[i], ord)
                TaylorSeries.mul!(m_c_ϕ_∂P_3[i], tmp1361[i], ∂P_3_sin_ϕ[i], ord)
                TaylorSeries.mul!(tmp1364[i], Λ2j_div_r4[i], 3, ord)
                TaylorSeries.mul!(F_J2_ξ[i], tmp1364[i], P_2_sin_ϕ[i], ord)
                TaylorSeries.mul!(F_J2_ζ[i], Λ2j_div_r4[i], m_c_ϕ_∂P_2[i], ord)
                TaylorSeries.mul!(tmp1368[i], Λ3j_div_r5[i], 4, ord)
                TaylorSeries.mul!(F_J3_ξ[i], tmp1368[i], P_3_sin_ϕ[i], ord)
                TaylorSeries.mul!(F_J3_ζ[i], Λ3j_div_r5[i], m_c_ϕ_∂P_3[i], ord)
                TaylorSeries.identity!(F_J_ξ[i], F_J2_ξ[i], ord)
                TaylorSeries.identity!(F_J_ζ[i], F_J2_ζ[i], ord)
                TaylorSeries.subst!(tmp1371[i], X[i], ord)
                TaylorSeries.div!(ξx[i], tmp1371[i], r_p1d2[i], ord)
                TaylorSeries.subst!(tmp1373[i], Y[i], ord)
                TaylorSeries.div!(ξy[i], tmp1373[i], r_p1d2[i], ord)
                TaylorSeries.subst!(tmp1375[i], Z[i], ord)
                TaylorSeries.div!(ξz[i], tmp1375[i], r_p1d2[i], ord)
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
            end
            TaylorSeries.pow!(tmp1405[i], ui[i], 2, ord)
            TaylorSeries.pow!(tmp1407[i], vi[i], 2, ord)
            TaylorSeries.add!(tmp1408[i], tmp1405[i], tmp1407[i], ord)
            TaylorSeries.pow!(tmp1410[i], wi[i], 2, ord)
            TaylorSeries.add!(v2[i], tmp1408[i], tmp1410[i], ord)
        end
        TaylorSeries.pow!(tmp1413, q[4], 2, ord)
        TaylorSeries.pow!(tmp1415, q[5], 2, ord)
        TaylorSeries.add!(tmp1416, tmp1413, tmp1415, ord)
        TaylorSeries.pow!(tmp1418, q[6], 2, ord)
        TaylorSeries.add!(v2[N], tmp1416, tmp1418, ord)
        for i = _1_to_Nm1
            TaylorSeries.add!(temp_004[i], newtonian1b_Potential[i], newtonianNb_Potential[N], ord)
            TaylorSeries.identity!(newtonianNb_Potential[N], temp_004[i], ord)
            if UJ_interaction[i]
                TaylorSeries.mul!(tmp1421[i], μ[i], F_J2_x[i], ord)
                TaylorSeries.subst!(temp_accX_i[i], accX, tmp1421[i], ord)
                TaylorSeries.identity!(accX, temp_accX_i[i], ord)
                TaylorSeries.mul!(tmp1423[i], μ[i], F_J2_y[i], ord)
                TaylorSeries.subst!(temp_accY_i[i], accY, tmp1423[i], ord)
                TaylorSeries.identity!(accY, temp_accY_i[i], ord)
                TaylorSeries.mul!(tmp1425[i], μ[i], F_J2_z[i], ord)
                TaylorSeries.subst!(temp_accZ_i[i], accZ, tmp1425[i], ord)
                TaylorSeries.identity!(accZ, temp_accZ_i[i], ord)
            end
        end
        TaylorSeries.mul!(_4ϕj[N], 4, newtonianNb_Potential[N], ord)
        for i = 1:10
            TaylorSeries.add!(ϕi_plus_4ϕj[i], newtonianNb_Potential_t[i], _4ϕj[N], ord)
            TaylorSeries.mul!(tmp1431[i], 2, v2[i], ord)
            TaylorSeries.mul!(tmp1433[i], 4, vi_dot_vj[i], ord)
            TaylorSeries.subst!(tmp1434[i], tmp1431[i], tmp1433[i], ord)
            TaylorSeries.add!(sj2_plus_2si2_minus_4vivj[i], tmp1434[i], v2[N], ord)
            TaylorSeries.subst!(ϕs_and_vs[i], sj2_plus_2si2_minus_4vivj[i], ϕi_plus_4ϕj[i], ord)
            TaylorSeries.mul!(Xij_t_Ui[i], X[i], ui[i], ord)
            TaylorSeries.mul!(Yij_t_Vi[i], Y[i], vi[i], ord)
            TaylorSeries.mul!(Zij_t_Wi[i], Z[i], wi[i], ord)
            TaylorSeries.add!(tmp1440[i], Xij_t_Ui[i], Yij_t_Vi[i], ord)
            TaylorSeries.add!(Rij_dot_Vi[i], tmp1440[i], Zij_t_Wi[i], ord)
            TaylorSeries.pow!(tmp1443[i], Rij_dot_Vi[i], 2, ord)
            TaylorSeries.div!(pn1t7[i], tmp1443[i], r_p2[i], ord)
            TaylorSeries.mul!(tmp1446[i], 1.5, pn1t7[i], ord)
            TaylorSeries.subst!(pn1t2_7[i], ϕs_and_vs[i], tmp1446[i], ord)
            TaylorSeries.add!(pn1t1_7[i], c_p2, pn1t2_7[i], ord)
            TaylorSeries.mul!(pNX_t_X[i], acceph_t[3i - 2], X[i], ord)
            TaylorSeries.mul!(pNY_t_Y[i], acceph_t[3i - 1], Y[i], ord)
            TaylorSeries.mul!(pNZ_t_Z[i], acceph_t[3i], Z[i], ord)
            TaylorSeries.add!(tmp1453[i], pNX_t_X[i], pNY_t_Y[i], ord)
            TaylorSeries.add!(tmp1454[i], tmp1453[i], pNZ_t_Z[i], ord)
            TaylorSeries.mul!(tmp1455[i], 0.5, tmp1454[i], ord)
            TaylorSeries.add!(pn1[i], pn1t1_7[i], tmp1455[i], ord)
            TaylorSeries.mul!(X_t_pn1[i], newton_acc_X[i], pn1[i], ord)
            TaylorSeries.mul!(Y_t_pn1[i], newton_acc_Y[i], pn1[i], ord)
            TaylorSeries.mul!(Z_t_pn1[i], newton_acc_Z[i], pn1[i], ord)
            TaylorSeries.mul!(pNX_t_pn3[i], acceph_t[3i - 2], pn3[i], ord)
            TaylorSeries.mul!(pNY_t_pn3[i], acceph_t[3i - 1], pn3[i], ord)
            TaylorSeries.mul!(pNZ_t_pn3[i], acceph_t[3i], pn3[i], ord)
        end
        for i = 1:10
            TaylorSeries.add!(tmp1463[i], U_t_pn2[i], pNX_t_pn3[i], ord)
            TaylorSeries.add!(termpnx[i], X_t_pn1[i], tmp1463[i], ord)
            TaylorSeries.add!(sumpnx[i], pntempX, termpnx[i], ord)
            TaylorSeries.identity!(pntempX, sumpnx[i], ord)
            TaylorSeries.add!(tmp1466[i], V_t_pn2[i], pNY_t_pn3[i], ord)
            TaylorSeries.add!(termpny[i], Y_t_pn1[i], tmp1466[i], ord)
            TaylorSeries.add!(sumpny[i], pntempY, termpny[i], ord)
            TaylorSeries.identity!(pntempY, sumpny[i], ord)
            TaylorSeries.add!(tmp1469[i], W_t_pn2[i], pNZ_t_pn3[i], ord)
            TaylorSeries.add!(termpnz[i], Z_t_pn1[i], tmp1469[i], ord)
            TaylorSeries.add!(sumpnz[i], pntempZ, termpnz[i], ord)
            TaylorSeries.identity!(pntempZ, sumpnz[i], ord)
        end
        for i = 11:27
            TaylorSeries.mul!(X_t_pn1[i], c_p2, newton_acc_X[i], ord)
            TaylorSeries.mul!(Y_t_pn1[i], c_p2, newton_acc_Y[i], ord)
            TaylorSeries.mul!(Z_t_pn1[i], c_p2, newton_acc_Z[i], ord)
        end
        for i = 11:27
            TaylorSeries.identity!(termpnx[i], X_t_pn1[i], ord)
            TaylorSeries.add!(sumpnx[i], pntempX, termpnx[i], ord)
            TaylorSeries.identity!(pntempX, sumpnx[i], ord)
            TaylorSeries.identity!(termpny[i], Y_t_pn1[i], ord)
            TaylorSeries.add!(sumpny[i], pntempY, termpny[i], ord)
            TaylorSeries.identity!(pntempY, sumpny[i], ord)
            TaylorSeries.identity!(termpnz[i], Z_t_pn1[i], ord)
            TaylorSeries.add!(sumpnz[i], pntempZ, termpnz[i], ord)
            TaylorSeries.identity!(pntempZ, sumpnz[i], ord)
        end
        TaylorSeries.mul!(postNewtonX, pntempX, c_m2, ord)
        TaylorSeries.mul!(postNewtonY, pntempY, c_m2, ord)
        TaylorSeries.mul!(postNewtonZ, pntempZ, c_m2, ord)
        TaylorSeries.mul!(tmp1481, Y[1], W[1], ord)
        TaylorSeries.mul!(tmp1482, Z[1], V[1], ord)
        TaylorSeries.subst!(hx, tmp1481, tmp1482, ord)
        TaylorSeries.mul!(tmp1484, Z[1], U[1], ord)
        TaylorSeries.mul!(tmp1485, X[1], W[1], ord)
        TaylorSeries.subst!(hy, tmp1484, tmp1485, ord)
        TaylorSeries.mul!(tmp1487, X[1], V[1], ord)
        TaylorSeries.mul!(tmp1488, Y[1], U[1], ord)
        TaylorSeries.subst!(hz, tmp1487, tmp1488, ord)
        TaylorSeries.mul!(tmp1490, hz, Y[1], ord)
        TaylorSeries.mul!(tmp1491, hy, Z[1], ord)
        TaylorSeries.subst!(tunitx0, tmp1490, tmp1491, ord)
        TaylorSeries.mul!(tmp1493, hx, Z[1], ord)
        TaylorSeries.mul!(tmp1494, hz, X[1], ord)
        TaylorSeries.subst!(tunity0, tmp1493, tmp1494, ord)
        TaylorSeries.mul!(tmp1496, hy, X[1], ord)
        TaylorSeries.mul!(tmp1497, hx, Y[1], ord)
        TaylorSeries.subst!(tunitz0, tmp1496, tmp1497, ord)
        TaylorSeries.pow!(tmp1500, tunitx0, 2, ord)
        TaylorSeries.pow!(tmp1502, tunity0, 2, ord)
        TaylorSeries.add!(tmp1503, tmp1500, tmp1502, ord)
        TaylorSeries.pow!(tmp1505, tunitz0, 2, ord)
        TaylorSeries.add!(tmp1506, tmp1503, tmp1505, ord)
        TaylorSeries.sqrt!(hmag, tmp1506, ord)
        TaylorSeries.div!(tunitx, tunitx0, hmag, ord)
        TaylorSeries.div!(tunity, tunity0, hmag, ord)
        TaylorSeries.div!(tunitz, tunitz0, hmag, ord)
        TaylorSeries.identity!(g_r, r_p2[1], ord)
        TaylorSeries.div!(A2_t_g_r, q[7], g_r, ord)
        TaylorSeries.mul!(NGAx, A2_t_g_r, tunitx, ord)
        TaylorSeries.mul!(NGAy, A2_t_g_r, tunity, ord)
        TaylorSeries.mul!(NGAz, A2_t_g_r, tunitz, ord)
        TaylorSeries.add!(tmp1515, postNewtonX, accX, ord)
        TaylorSeries.add!(dq[4], tmp1515, NGAx, ord)
        TaylorSeries.add!(tmp1517, postNewtonY, accY, ord)
        TaylorSeries.add!(dq[5], tmp1517, NGAy, ord)
        TaylorSeries.add!(tmp1519, postNewtonZ, accZ, ord)
        TaylorSeries.add!(dq[6], tmp1519, NGAz, ord)
        TaylorSeries.identity!(dq[7], zero_q_1, ord)
        for __idx = eachindex(q)
            (q[__idx]).coeffs[ordnext + 1] = (dq[__idx]).coeffs[ordnext] / ordnext
        end
    end
    return nothing
end

function TaylorIntegration.jetcoeffs!(::Val{RNp1BP_pN_A_J23E_J2S_ng_eph_threads!}, t::Taylor1{_T}, q::AbstractVector{Taylor1{_S}}, dq::AbstractVector{Taylor1{_S}}, params) where {_T <: Real, _S <: Number}
    order = t.order
    local jd0 = params[4]
    local ss16asteph_t = evaleph(params[1], t, q[1])
    local acceph_t = evaleph(params[2], t, q[1])
    local newtonianNb_Potential_t = evaleph(params[3], t, q[1])
    local S = eltype(q[1])
    local N = length(μ)
    local _1_to_Nm1 = 1:N - 1
    local c_p2 = 29979.063823897606
    local c_m2 = 3.3356611996764786e-5
    local zero_q_1 = zero(q[1])
    X = Array{Taylor1{S}}(undef, N)
    Y = Array{Taylor1{S}}(undef, N)
    Z = Array{Taylor1{S}}(undef, N)
    r_p2 = Array{Taylor1{S}}(undef, N)
    r_p3d2 = Array{Taylor1{S}}(undef, N)
    r_p7d2 = Array{Taylor1{S}}(undef, N)
    newtonianCoeff = Array{Taylor1{S}}(undef, N)
    ui = Array{Taylor1{S}}(undef, N - 1)
    vi = Array{Taylor1{S}}(undef, N - 1)
    wi = Array{Taylor1{S}}(undef, N - 1)
    U = Array{Taylor1{S}}(undef, N)
    V = Array{Taylor1{S}}(undef, N)
    W = Array{Taylor1{S}}(undef, N)
    _4U_m_3X = Array{Taylor1{S}}(undef, N)
    _4V_m_3Y = Array{Taylor1{S}}(undef, N)
    _4W_m_3Z = Array{Taylor1{S}}(undef, N)
    UU = Array{Taylor1{S}}(undef, N)
    VV = Array{Taylor1{S}}(undef, N)
    WW = Array{Taylor1{S}}(undef, N)
    r_p1d2 = Array{Taylor1{S}}(undef, N)
    newtonianNb_Potential = Array{Taylor1{S}}(undef, N)
    newtonian1b_Potential = Array{Taylor1{S}}(undef, N)
    newton_acc_X = Array{Taylor1{S}}(undef, N)
    newton_acc_Y = Array{Taylor1{S}}(undef, N)
    newton_acc_Z = Array{Taylor1{S}}(undef, N)
    v2 = Array{Taylor1{S}}(undef, N)
    vi_dot_vj = Array{Taylor1{S}}(undef, N)
    pn2 = Array{Taylor1{S}}(undef, N)
    pn3 = Array{Taylor1{S}}(undef, N)
    _4ϕj = Array{Taylor1{S}}(undef, N)
    ϕi_plus_4ϕj = Array{Taylor1{S}}(undef, N)
    sj2_plus_2si2_minus_4vivj = Array{Taylor1{S}}(undef, N)
    ϕs_and_vs = Array{Taylor1{S}}(undef, N)
    U_t_pn2 = Array{Taylor1{S}}(undef, N)
    V_t_pn2 = Array{Taylor1{S}}(undef, N)
    W_t_pn2 = Array{Taylor1{S}}(undef, N)
    pn1t1_7 = Array{Taylor1{S}}(undef, N)
    pntempX = Taylor1(identity(constant_term(zero_q_1)), order)
    pntempY = Taylor1(identity(constant_term(zero_q_1)), order)
    pntempZ = Taylor1(identity(constant_term(zero_q_1)), order)
    pn1 = Array{Taylor1{S}}(undef, N)
    X_t_pn1 = Array{Taylor1{S}}(undef, N)
    Y_t_pn1 = Array{Taylor1{S}}(undef, N)
    Z_t_pn1 = Array{Taylor1{S}}(undef, N)
    pNX_t_pn3 = Array{Taylor1{S}}(undef, N)
    pNY_t_pn3 = Array{Taylor1{S}}(undef, N)
    pNZ_t_pn3 = Array{Taylor1{S}}(undef, N)
    pNX_t_X = Array{Taylor1{S}}(undef, N)
    pNY_t_Y = Array{Taylor1{S}}(undef, N)
    pNZ_t_Z = Array{Taylor1{S}}(undef, N)
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
    local dsj2k = t + (jd0 - 2.451545e6)
    local M_ = Array{Taylor1{S}}(undef, 3, 3, N)
    local M_[:, :, ea] = t2c_jpl_de430(dsj2k)
    dq[1] = Taylor1(identity(constant_term(q[4])), order)
    dq[2] = Taylor1(identity(constant_term(q[5])), order)
    dq[3] = Taylor1(identity(constant_term(q[6])), order)
    newtonianNb_Potential[N] = Taylor1(identity(constant_term(zero_q_1)), order)
    tmp2543 = Array{Taylor1{_S}}(undef, size(dq))
    tmp2543 .= Taylor1(zero(_S), order)
    tmp2545 = Array{Taylor1{_S}}(undef, size(ui))
    tmp2545 .= Taylor1(zero(_S), order)
    tmp2548 = Array{Taylor1{_S}}(undef, size(dq))
    tmp2548 .= Taylor1(zero(_S), order)
    tmp2550 = Array{Taylor1{_S}}(undef, size(vi))
    tmp2550 .= Taylor1(zero(_S), order)
    tmp2553 = Array{Taylor1{_S}}(undef, size(dq))
    tmp2553 .= Taylor1(zero(_S), order)
    tmp2555 = Array{Taylor1{_S}}(undef, size(wi))
    tmp2555 .= Taylor1(zero(_S), order)
    pn2x = Array{Taylor1{_S}}(undef, size(X))
    pn2x .= Taylor1(zero(_S), order)
    pn2y = Array{Taylor1{_S}}(undef, size(Y))
    pn2y .= Taylor1(zero(_S), order)
    pn2z = Array{Taylor1{_S}}(undef, size(Z))
    pn2z .= Taylor1(zero(_S), order)
    tmp2563 = Array{Taylor1{_S}}(undef, size(UU))
    tmp2563 .= Taylor1(zero(_S), order)
    tmp2566 = Array{Taylor1{_S}}(undef, size(X))
    tmp2566 .= Taylor1(zero(_S), order)
    tmp2568 = Array{Taylor1{_S}}(undef, size(Y))
    tmp2568 .= Taylor1(zero(_S), order)
    tmp2569 = Array{Taylor1{_S}}(undef, size(tmp2566))
    tmp2569 .= Taylor1(zero(_S), order)
    tmp2571 = Array{Taylor1{_S}}(undef, size(Z))
    tmp2571 .= Taylor1(zero(_S), order)
    tmp2579 = Array{Taylor1{_S}}(undef, size(pn2x))
    tmp2579 .= Taylor1(zero(_S), order)
    tmp2580 = Array{Taylor1{_S}}(undef, size(tmp2579))
    tmp2580 .= Taylor1(zero(_S), order)
    tmp2591 = Array{Taylor1{_S}}(undef, size(X))
    tmp2591 .= Taylor1(zero(_S), order)
    tmp2593 = Array{Taylor1{_S}}(undef, size(Y))
    tmp2593 .= Taylor1(zero(_S), order)
    tmp2595 = Array{Taylor1{_S}}(undef, size(Z))
    tmp2595 .= Taylor1(zero(_S), order)
    tmp2597 = Array{Taylor1{_S}}(undef, size(t31))
    tmp2597 .= Taylor1(zero(_S), order)
    tmp2791 = Array{Taylor1{_S}}(undef, size(sin_ϕ))
    tmp2791 .= Taylor1(zero(_S), order)
    tmp2792 = Array{Taylor1{_S}}(undef, size(ϕ))
    tmp2792 .= Taylor1(zero(_S), order)
    tmp2607 = Array{Taylor1{_S}}(undef, size(sin2_ϕ))
    tmp2607 .= Taylor1(zero(_S), order)
    tmp2613 = Array{Taylor1{_S}}(undef, size(sin_ϕ))
    tmp2613 .= Taylor1(zero(_S), order)
    tmp2615 = Array{Taylor1{_S}}(undef, size(sin3_ϕ))
    tmp2615 .= Taylor1(zero(_S), order)
    tmp2619 = Array{Taylor1{_S}}(undef, size(sin2_ϕ))
    tmp2619 .= Taylor1(zero(_S), order)
    tmp2621 = Array{Taylor1{_S}}(undef, size(Λ2))
    tmp2621 .= Taylor1(zero(_S), order)
    tmp2623 = Array{Taylor1{_S}}(undef, size(r_p2))
    tmp2623 .= Taylor1(zero(_S), order)
    tmp2625 = Array{Taylor1{_S}}(undef, size(Λ3))
    tmp2625 .= Taylor1(zero(_S), order)
    tmp2627 = Array{Taylor1{_S}}(undef, size(r_p1d2))
    tmp2627 .= Taylor1(zero(_S), order)
    tmp2629 = Array{Taylor1{_S}}(undef, size(cos_ϕ))
    tmp2629 .= Taylor1(zero(_S), order)
    tmp2631 = Array{Taylor1{_S}}(undef, size(cos_ϕ))
    tmp2631 .= Taylor1(zero(_S), order)
    tmp2634 = Array{Taylor1{_S}}(undef, size(Λ2j_div_r4))
    tmp2634 .= Taylor1(zero(_S), order)
    tmp2638 = Array{Taylor1{_S}}(undef, size(Λ3j_div_r5))
    tmp2638 .= Taylor1(zero(_S), order)
    tmp2641 = Array{Taylor1{_S}}(undef, size(X))
    tmp2641 .= Taylor1(zero(_S), order)
    tmp2643 = Array{Taylor1{_S}}(undef, size(Y))
    tmp2643 .= Taylor1(zero(_S), order)
    tmp2645 = Array{Taylor1{_S}}(undef, size(Z))
    tmp2645 .= Taylor1(zero(_S), order)
    tmp2675 = Array{Taylor1{_S}}(undef, size(ui))
    tmp2675 .= Taylor1(zero(_S), order)
    tmp2677 = Array{Taylor1{_S}}(undef, size(vi))
    tmp2677 .= Taylor1(zero(_S), order)
    tmp2678 = Array{Taylor1{_S}}(undef, size(tmp2675))
    tmp2678 .= Taylor1(zero(_S), order)
    tmp2680 = Array{Taylor1{_S}}(undef, size(wi))
    tmp2680 .= Taylor1(zero(_S), order)
    #= REPL[4]:161 =# Threads.@threads for i = _1_to_Nm1
            ui[i] = Taylor1(identity(constant_term(ss16asteph_t[3 * ((N - 1) + i) - 2])), order)
            vi[i] = Taylor1(identity(constant_term(ss16asteph_t[3 * ((N - 1) + i) - 1])), order)
            wi[i] = Taylor1(identity(constant_term(ss16asteph_t[3 * ((N - 1) + i)])), order)
            X[i] = Taylor1(constant_term(ss16asteph_t[3i - 2]) - constant_term(q[1]), order)
            Y[i] = Taylor1(constant_term(ss16asteph_t[3i - 1]) - constant_term(q[2]), order)
            Z[i] = Taylor1(constant_term(ss16asteph_t[3i]) - constant_term(q[3]), order)
            U[i] = Taylor1(constant_term(ui[i]) - constant_term(dq[1]), order)
            V[i] = Taylor1(constant_term(vi[i]) - constant_term(dq[2]), order)
            W[i] = Taylor1(constant_term(wi[i]) - constant_term(dq[3]), order)
            tmp2543[1] = Taylor1(constant_term(4) * constant_term(dq[1]), order)
            tmp2545[i] = Taylor1(constant_term(3) * constant_term(ui[i]), order)
            _4U_m_3X[i] = Taylor1(constant_term(tmp2543[1]) - constant_term(tmp2545[i]), order)
            tmp2548[2] = Taylor1(constant_term(4) * constant_term(dq[2]), order)
            tmp2550[i] = Taylor1(constant_term(3) * constant_term(vi[i]), order)
            _4V_m_3Y[i] = Taylor1(constant_term(tmp2548[2]) - constant_term(tmp2550[i]), order)
            tmp2553[3] = Taylor1(constant_term(4) * constant_term(dq[3]), order)
            tmp2555[i] = Taylor1(constant_term(3) * constant_term(wi[i]), order)
            _4W_m_3Z[i] = Taylor1(constant_term(tmp2553[3]) - constant_term(tmp2555[i]), order)
            pn2x[i] = Taylor1(constant_term(X[i]) * constant_term(_4U_m_3X[i]), order)
            pn2y[i] = Taylor1(constant_term(Y[i]) * constant_term(_4V_m_3Y[i]), order)
            pn2z[i] = Taylor1(constant_term(Z[i]) * constant_term(_4W_m_3Z[i]), order)
            UU[i] = Taylor1(constant_term(ui[i]) * constant_term(dq[1]), order)
            VV[i] = Taylor1(constant_term(vi[i]) * constant_term(dq[2]), order)
            WW[i] = Taylor1(constant_term(wi[i]) * constant_term(dq[3]), order)
            tmp2563[i] = Taylor1(constant_term(UU[i]) + constant_term(VV[i]), order)
            vi_dot_vj[i] = Taylor1(constant_term(tmp2563[i]) + constant_term(WW[i]), order)
            tmp2566[i] = Taylor1(constant_term(X[i]) ^ float(constant_term(2)), order)
            tmp2568[i] = Taylor1(constant_term(Y[i]) ^ float(constant_term(2)), order)
            tmp2569[i] = Taylor1(constant_term(tmp2566[i]) + constant_term(tmp2568[i]), order)
            tmp2571[i] = Taylor1(constant_term(Z[i]) ^ float(constant_term(2)), order)
            r_p2[i] = Taylor1(constant_term(tmp2569[i]) + constant_term(tmp2571[i]), order)
            r_p1d2[i] = Taylor1(sqrt(constant_term(r_p2[i])), order)
            r_p3d2[i] = Taylor1(constant_term(r_p2[i]) ^ float(constant_term(1.5)), order)
            r_p7d2[i] = Taylor1(constant_term(r_p2[i]) ^ float(constant_term(3.5)), order)
            newtonianCoeff[i] = Taylor1(constant_term(μ[i]) / constant_term(r_p3d2[i]), order)
            tmp2579[i] = Taylor1(constant_term(pn2x[i]) + constant_term(pn2y[i]), order)
            tmp2580[i] = Taylor1(constant_term(tmp2579[i]) + constant_term(pn2z[i]), order)
            pn2[i] = Taylor1(constant_term(newtonianCoeff[i]) * constant_term(tmp2580[i]), order)
            newton_acc_X[i] = Taylor1(constant_term(X[i]) * constant_term(newtonianCoeff[i]), order)
            newton_acc_Y[i] = Taylor1(constant_term(Y[i]) * constant_term(newtonianCoeff[i]), order)
            newton_acc_Z[i] = Taylor1(constant_term(Z[i]) * constant_term(newtonianCoeff[i]), order)
            newtonian1b_Potential[i] = Taylor1(constant_term(μ[i]) / constant_term(r_p1d2[i]), order)
            pn3[i] = Taylor1(constant_term(3.5) * constant_term(newtonian1b_Potential[i]), order)
            U_t_pn2[i] = Taylor1(constant_term(pn2[i]) * constant_term(U[i]), order)
            V_t_pn2[i] = Taylor1(constant_term(pn2[i]) * constant_term(V[i]), order)
            W_t_pn2[i] = Taylor1(constant_term(pn2[i]) * constant_term(W[i]), order)
            if UJ_interaction[i]
                tmp2591[i] = Taylor1(-(constant_term(X[i])), order)
                t31[i] = Taylor1(constant_term(tmp2591[i]) * constant_term(M_[1, 3, i]), order)
                tmp2593[i] = Taylor1(-(constant_term(Y[i])), order)
                t32[i] = Taylor1(constant_term(tmp2593[i]) * constant_term(M_[2, 3, i]), order)
                tmp2595[i] = Taylor1(-(constant_term(Z[i])), order)
                t33[i] = Taylor1(constant_term(tmp2595[i]) * constant_term(M_[3, 3, i]), order)
                tmp2597[i] = Taylor1(constant_term(t31[i]) + constant_term(t32[i]), order)
                r_sin_ϕ[i] = Taylor1(constant_term(tmp2597[i]) + constant_term(t33[i]), order)
                sin_ϕ[i] = Taylor1(constant_term(r_sin_ϕ[i]) / constant_term(r_p1d2[i]), order)
                ϕ[i] = Taylor1(asin(constant_term(sin_ϕ[i])), order)
                tmp2791[i] = Taylor1(sqrt(1 - constant_term(sin_ϕ[i]) ^ 2), order)
                cos_ϕ[i] = Taylor1(cos(constant_term(ϕ[i])), order)
                tmp2792[i] = Taylor1(sin(constant_term(ϕ[i])), order)
                sin2_ϕ[i] = Taylor1(constant_term(sin_ϕ[i]) ^ float(constant_term(2)), order)
                sin3_ϕ[i] = Taylor1(constant_term(sin_ϕ[i]) ^ float(constant_term(3)), order)
                tmp2607[i] = Taylor1(constant_term(1.5) * constant_term(sin2_ϕ[i]), order)
                P_2_sin_ϕ[i] = Taylor1(constant_term(tmp2607[i]) - constant_term(0.5), order)
                ∂P_2_sin_ϕ[i] = Taylor1(constant_term(3) * constant_term(sin_ϕ[i]), order)
                tmp2613[i] = Taylor1(constant_term(-1.5) * constant_term(sin_ϕ[i]), order)
                tmp2615[i] = Taylor1(constant_term(2.5) * constant_term(sin3_ϕ[i]), order)
                P_3_sin_ϕ[i] = Taylor1(constant_term(tmp2613[i]) + constant_term(tmp2615[i]), order)
                tmp2619[i] = Taylor1(constant_term(7.5) * constant_term(sin2_ϕ[i]), order)
                ∂P_3_sin_ϕ[i] = Taylor1(constant_term(-1.5) + constant_term(tmp2619[i]), order)
                tmp2621[i] = Taylor1(-(constant_term(Λ2[i])), order)
                tmp2623[i] = Taylor1(constant_term(r_p2[i]) ^ float(constant_term(2)), order)
                Λ2j_div_r4[i] = Taylor1(constant_term(tmp2621[i]) / constant_term(tmp2623[i]), order)
                tmp2625[i] = Taylor1(-(constant_term(Λ3[i])), order)
                tmp2627[i] = Taylor1(constant_term(r_p1d2[i]) ^ float(constant_term(5)), order)
                Λ3j_div_r5[i] = Taylor1(constant_term(tmp2625[i]) / constant_term(tmp2627[i]), order)
                tmp2629[i] = Taylor1(-(constant_term(cos_ϕ[i])), order)
                m_c_ϕ_∂P_2[i] = Taylor1(constant_term(tmp2629[i]) * constant_term(∂P_2_sin_ϕ[i]), order)
                tmp2631[i] = Taylor1(-(constant_term(cos_ϕ[i])), order)
                m_c_ϕ_∂P_3[i] = Taylor1(constant_term(tmp2631[i]) * constant_term(∂P_3_sin_ϕ[i]), order)
                tmp2634[i] = Taylor1(constant_term(Λ2j_div_r4[i]) * constant_term(3), order)
                F_J2_ξ[i] = Taylor1(constant_term(tmp2634[i]) * constant_term(P_2_sin_ϕ[i]), order)
                F_J2_ζ[i] = Taylor1(constant_term(Λ2j_div_r4[i]) * constant_term(m_c_ϕ_∂P_2[i]), order)
                tmp2638[i] = Taylor1(constant_term(Λ3j_div_r5[i]) * constant_term(4), order)
                F_J3_ξ[i] = Taylor1(constant_term(tmp2638[i]) * constant_term(P_3_sin_ϕ[i]), order)
                F_J3_ζ[i] = Taylor1(constant_term(Λ3j_div_r5[i]) * constant_term(m_c_ϕ_∂P_3[i]), order)
                F_J_ξ[i] = Taylor1(identity(constant_term(F_J2_ξ[i])), order)
                F_J_ζ[i] = Taylor1(identity(constant_term(F_J2_ζ[i])), order)
                tmp2641[i] = Taylor1(-(constant_term(X[i])), order)
                ξx[i] = Taylor1(constant_term(tmp2641[i]) / constant_term(r_p1d2[i]), order)
                tmp2643[i] = Taylor1(-(constant_term(Y[i])), order)
                ξy[i] = Taylor1(constant_term(tmp2643[i]) / constant_term(r_p1d2[i]), order)
                tmp2645[i] = Taylor1(-(constant_term(Z[i])), order)
                ξz[i] = Taylor1(constant_term(tmp2645[i]) / constant_term(r_p1d2[i]), order)
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
            end
            tmp2675[i] = Taylor1(constant_term(ui[i]) ^ float(constant_term(2)), order)
            tmp2677[i] = Taylor1(constant_term(vi[i]) ^ float(constant_term(2)), order)
            tmp2678[i] = Taylor1(constant_term(tmp2675[i]) + constant_term(tmp2677[i]), order)
            tmp2680[i] = Taylor1(constant_term(wi[i]) ^ float(constant_term(2)), order)
            v2[i] = Taylor1(constant_term(tmp2678[i]) + constant_term(tmp2680[i]), order)
        end
    tmp2683 = Taylor1(constant_term(q[4]) ^ float(constant_term(2)), order)
    tmp2685 = Taylor1(constant_term(q[5]) ^ float(constant_term(2)), order)
    tmp2686 = Taylor1(constant_term(tmp2683) + constant_term(tmp2685), order)
    tmp2688 = Taylor1(constant_term(q[6]) ^ float(constant_term(2)), order)
    v2[N] = Taylor1(constant_term(tmp2686) + constant_term(tmp2688), order)
    temp_004 = Array{Taylor1{_S}}(undef, size(newtonian1b_Potential))
    temp_004 .= Taylor1(zero(_S), order)
    tmp2691 = Array{Taylor1{_S}}(undef, size(μ))
    tmp2691 .= Taylor1(zero(_S), order)
    tmp2693 = Array{Taylor1{_S}}(undef, size(μ))
    tmp2693 .= Taylor1(zero(_S), order)
    tmp2695 = Array{Taylor1{_S}}(undef, size(μ))
    tmp2695 .= Taylor1(zero(_S), order)
    for i = _1_to_Nm1
        temp_004[i] = Taylor1(constant_term(newtonian1b_Potential[i]) + constant_term(newtonianNb_Potential[N]), order)
        newtonianNb_Potential[N] = Taylor1(identity(constant_term(temp_004[i])), order)
        if UJ_interaction[i]
            tmp2691[i] = Taylor1(constant_term(μ[i]) * constant_term(F_J2_x[i]), order)
            temp_accX_i[i] = Taylor1(constant_term(accX) - constant_term(tmp2691[i]), order)
            accX = Taylor1(identity(constant_term(temp_accX_i[i])), order)
            tmp2693[i] = Taylor1(constant_term(μ[i]) * constant_term(F_J2_y[i]), order)
            temp_accY_i[i] = Taylor1(constant_term(accY) - constant_term(tmp2693[i]), order)
            accY = Taylor1(identity(constant_term(temp_accY_i[i])), order)
            tmp2695[i] = Taylor1(constant_term(μ[i]) * constant_term(F_J2_z[i]), order)
            temp_accZ_i[i] = Taylor1(constant_term(accZ) - constant_term(tmp2695[i]), order)
            accZ = Taylor1(identity(constant_term(temp_accZ_i[i])), order)
        end
    end
    _4ϕj[N] = Taylor1(constant_term(4) * constant_term(newtonianNb_Potential[N]), order)
    tmp2701 = Array{Taylor1{_S}}(undef, size(v2))
    tmp2701 .= Taylor1(zero(_S), order)
    tmp2703 = Array{Taylor1{_S}}(undef, size(vi_dot_vj))
    tmp2703 .= Taylor1(zero(_S), order)
    tmp2704 = Array{Taylor1{_S}}(undef, size(tmp2701))
    tmp2704 .= Taylor1(zero(_S), order)
    Xij_t_Ui = Array{Taylor1{_S}}(undef, size(X))
    Xij_t_Ui .= Taylor1(zero(_S), order)
    Yij_t_Vi = Array{Taylor1{_S}}(undef, size(Y))
    Yij_t_Vi .= Taylor1(zero(_S), order)
    Zij_t_Wi = Array{Taylor1{_S}}(undef, size(Z))
    Zij_t_Wi .= Taylor1(zero(_S), order)
    tmp2710 = Array{Taylor1{_S}}(undef, size(Xij_t_Ui))
    tmp2710 .= Taylor1(zero(_S), order)
    Rij_dot_Vi = Array{Taylor1{_S}}(undef, size(tmp2710))
    Rij_dot_Vi .= Taylor1(zero(_S), order)
    tmp2713 = Array{Taylor1{_S}}(undef, size(Rij_dot_Vi))
    tmp2713 .= Taylor1(zero(_S), order)
    pn1t7 = Array{Taylor1{_S}}(undef, size(tmp2713))
    pn1t7 .= Taylor1(zero(_S), order)
    tmp2716 = Array{Taylor1{_S}}(undef, size(pn1t7))
    tmp2716 .= Taylor1(zero(_S), order)
    pn1t2_7 = Array{Taylor1{_S}}(undef, size(ϕs_and_vs))
    pn1t2_7 .= Taylor1(zero(_S), order)
    tmp2723 = Array{Taylor1{_S}}(undef, size(pNX_t_X))
    tmp2723 .= Taylor1(zero(_S), order)
    tmp2724 = Array{Taylor1{_S}}(undef, size(tmp2723))
    tmp2724 .= Taylor1(zero(_S), order)
    tmp2725 = Array{Taylor1{_S}}(undef, size(tmp2724))
    tmp2725 .= Taylor1(zero(_S), order)
    #= REPL[4]:306 =# Threads.@threads for i = 1:10
            ϕi_plus_4ϕj[i] = Taylor1(constant_term(newtonianNb_Potential_t[i]) + constant_term(_4ϕj[N]), order)
            tmp2701[i] = Taylor1(constant_term(2) * constant_term(v2[i]), order)
            tmp2703[i] = Taylor1(constant_term(4) * constant_term(vi_dot_vj[i]), order)
            tmp2704[i] = Taylor1(constant_term(tmp2701[i]) - constant_term(tmp2703[i]), order)
            sj2_plus_2si2_minus_4vivj[i] = Taylor1(constant_term(tmp2704[i]) + constant_term(v2[N]), order)
            ϕs_and_vs[i] = Taylor1(constant_term(sj2_plus_2si2_minus_4vivj[i]) - constant_term(ϕi_plus_4ϕj[i]), order)
            Xij_t_Ui[i] = Taylor1(constant_term(X[i]) * constant_term(ui[i]), order)
            Yij_t_Vi[i] = Taylor1(constant_term(Y[i]) * constant_term(vi[i]), order)
            Zij_t_Wi[i] = Taylor1(constant_term(Z[i]) * constant_term(wi[i]), order)
            tmp2710[i] = Taylor1(constant_term(Xij_t_Ui[i]) + constant_term(Yij_t_Vi[i]), order)
            Rij_dot_Vi[i] = Taylor1(constant_term(tmp2710[i]) + constant_term(Zij_t_Wi[i]), order)
            tmp2713[i] = Taylor1(constant_term(Rij_dot_Vi[i]) ^ float(constant_term(2)), order)
            pn1t7[i] = Taylor1(constant_term(tmp2713[i]) / constant_term(r_p2[i]), order)
            tmp2716[i] = Taylor1(constant_term(1.5) * constant_term(pn1t7[i]), order)
            pn1t2_7[i] = Taylor1(constant_term(ϕs_and_vs[i]) - constant_term(tmp2716[i]), order)
            pn1t1_7[i] = Taylor1(constant_term(c_p2) + constant_term(pn1t2_7[i]), order)
            pNX_t_X[i] = Taylor1(constant_term(acceph_t[3i - 2]) * constant_term(X[i]), order)
            pNY_t_Y[i] = Taylor1(constant_term(acceph_t[3i - 1]) * constant_term(Y[i]), order)
            pNZ_t_Z[i] = Taylor1(constant_term(acceph_t[3i]) * constant_term(Z[i]), order)
            tmp2723[i] = Taylor1(constant_term(pNX_t_X[i]) + constant_term(pNY_t_Y[i]), order)
            tmp2724[i] = Taylor1(constant_term(tmp2723[i]) + constant_term(pNZ_t_Z[i]), order)
            tmp2725[i] = Taylor1(constant_term(0.5) * constant_term(tmp2724[i]), order)
            pn1[i] = Taylor1(constant_term(pn1t1_7[i]) + constant_term(tmp2725[i]), order)
            X_t_pn1[i] = Taylor1(constant_term(newton_acc_X[i]) * constant_term(pn1[i]), order)
            Y_t_pn1[i] = Taylor1(constant_term(newton_acc_Y[i]) * constant_term(pn1[i]), order)
            Z_t_pn1[i] = Taylor1(constant_term(newton_acc_Z[i]) * constant_term(pn1[i]), order)
            pNX_t_pn3[i] = Taylor1(constant_term(acceph_t[3i - 2]) * constant_term(pn3[i]), order)
            pNY_t_pn3[i] = Taylor1(constant_term(acceph_t[3i - 1]) * constant_term(pn3[i]), order)
            pNZ_t_pn3[i] = Taylor1(constant_term(acceph_t[3i]) * constant_term(pn3[i]), order)
        end
    tmp2733 = Array{Taylor1{_S}}(undef, size(U_t_pn2))
    tmp2733 .= Taylor1(zero(_S), order)
    termpnx = Array{Taylor1{_S}}(undef, size(X_t_pn1))
    termpnx .= Taylor1(zero(_S), order)
    sumpnx = Array{Taylor1{_S}}(undef, size(termpnx))
    sumpnx .= Taylor1(zero(_S), order)
    tmp2736 = Array{Taylor1{_S}}(undef, size(V_t_pn2))
    tmp2736 .= Taylor1(zero(_S), order)
    termpny = Array{Taylor1{_S}}(undef, size(Y_t_pn1))
    termpny .= Taylor1(zero(_S), order)
    sumpny = Array{Taylor1{_S}}(undef, size(termpny))
    sumpny .= Taylor1(zero(_S), order)
    tmp2739 = Array{Taylor1{_S}}(undef, size(W_t_pn2))
    tmp2739 .= Taylor1(zero(_S), order)
    termpnz = Array{Taylor1{_S}}(undef, size(Z_t_pn1))
    termpnz .= Taylor1(zero(_S), order)
    sumpnz = Array{Taylor1{_S}}(undef, size(termpnz))
    sumpnz .= Taylor1(zero(_S), order)
    for i = 1:10
        tmp2733[i] = Taylor1(constant_term(U_t_pn2[i]) + constant_term(pNX_t_pn3[i]), order)
        termpnx[i] = Taylor1(constant_term(X_t_pn1[i]) + constant_term(tmp2733[i]), order)
        sumpnx[i] = Taylor1(constant_term(pntempX) + constant_term(termpnx[i]), order)
        pntempX = Taylor1(identity(constant_term(sumpnx[i])), order)
        tmp2736[i] = Taylor1(constant_term(V_t_pn2[i]) + constant_term(pNY_t_pn3[i]), order)
        termpny[i] = Taylor1(constant_term(Y_t_pn1[i]) + constant_term(tmp2736[i]), order)
        sumpny[i] = Taylor1(constant_term(pntempY) + constant_term(termpny[i]), order)
        pntempY = Taylor1(identity(constant_term(sumpny[i])), order)
        tmp2739[i] = Taylor1(constant_term(W_t_pn2[i]) + constant_term(pNZ_t_pn3[i]), order)
        termpnz[i] = Taylor1(constant_term(Z_t_pn1[i]) + constant_term(tmp2739[i]), order)
        sumpnz[i] = Taylor1(constant_term(pntempZ) + constant_term(termpnz[i]), order)
        pntempZ = Taylor1(identity(constant_term(sumpnz[i])), order)
    end
    #= REPL[4]:345 =# Threads.@threads for i = 11:27
            X_t_pn1[i] = Taylor1(constant_term(c_p2) * constant_term(newton_acc_X[i]), order)
            Y_t_pn1[i] = Taylor1(constant_term(c_p2) * constant_term(newton_acc_Y[i]), order)
            Z_t_pn1[i] = Taylor1(constant_term(c_p2) * constant_term(newton_acc_Z[i]), order)
        end
    for i = 11:27
        termpnx[i] = Taylor1(identity(constant_term(X_t_pn1[i])), order)
        sumpnx[i] = Taylor1(constant_term(pntempX) + constant_term(termpnx[i]), order)
        pntempX = Taylor1(identity(constant_term(sumpnx[i])), order)
        termpny[i] = Taylor1(identity(constant_term(Y_t_pn1[i])), order)
        sumpny[i] = Taylor1(constant_term(pntempY) + constant_term(termpny[i]), order)
        pntempY = Taylor1(identity(constant_term(sumpny[i])), order)
        termpnz[i] = Taylor1(identity(constant_term(Z_t_pn1[i])), order)
        sumpnz[i] = Taylor1(constant_term(pntempZ) + constant_term(termpnz[i]), order)
        pntempZ = Taylor1(identity(constant_term(sumpnz[i])), order)
    end
    postNewtonX = Taylor1(constant_term(pntempX) * constant_term(c_m2), order)
    postNewtonY = Taylor1(constant_term(pntempY) * constant_term(c_m2), order)
    postNewtonZ = Taylor1(constant_term(pntempZ) * constant_term(c_m2), order)
    tmp2751 = Taylor1(constant_term(Y[1]) * constant_term(W[1]), order)
    tmp2752 = Taylor1(constant_term(Z[1]) * constant_term(V[1]), order)
    hx = Taylor1(constant_term(tmp2751) - constant_term(tmp2752), order)
    tmp2754 = Taylor1(constant_term(Z[1]) * constant_term(U[1]), order)
    tmp2755 = Taylor1(constant_term(X[1]) * constant_term(W[1]), order)
    hy = Taylor1(constant_term(tmp2754) - constant_term(tmp2755), order)
    tmp2757 = Taylor1(constant_term(X[1]) * constant_term(V[1]), order)
    tmp2758 = Taylor1(constant_term(Y[1]) * constant_term(U[1]), order)
    hz = Taylor1(constant_term(tmp2757) - constant_term(tmp2758), order)
    tmp2760 = Taylor1(constant_term(hz) * constant_term(Y[1]), order)
    tmp2761 = Taylor1(constant_term(hy) * constant_term(Z[1]), order)
    tunitx0 = Taylor1(constant_term(tmp2760) - constant_term(tmp2761), order)
    tmp2763 = Taylor1(constant_term(hx) * constant_term(Z[1]), order)
    tmp2764 = Taylor1(constant_term(hz) * constant_term(X[1]), order)
    tunity0 = Taylor1(constant_term(tmp2763) - constant_term(tmp2764), order)
    tmp2766 = Taylor1(constant_term(hy) * constant_term(X[1]), order)
    tmp2767 = Taylor1(constant_term(hx) * constant_term(Y[1]), order)
    tunitz0 = Taylor1(constant_term(tmp2766) - constant_term(tmp2767), order)
    tmp2770 = Taylor1(constant_term(tunitx0) ^ float(constant_term(2)), order)
    tmp2772 = Taylor1(constant_term(tunity0) ^ float(constant_term(2)), order)
    tmp2773 = Taylor1(constant_term(tmp2770) + constant_term(tmp2772), order)
    tmp2775 = Taylor1(constant_term(tunitz0) ^ float(constant_term(2)), order)
    tmp2776 = Taylor1(constant_term(tmp2773) + constant_term(tmp2775), order)
    hmag = Taylor1(sqrt(constant_term(tmp2776)), order)
    tunitx = Taylor1(constant_term(tunitx0) / constant_term(hmag), order)
    tunity = Taylor1(constant_term(tunity0) / constant_term(hmag), order)
    tunitz = Taylor1(constant_term(tunitz0) / constant_term(hmag), order)
    g_r = Taylor1(identity(constant_term(r_p2[1])), order)
    A2_t_g_r = Taylor1(constant_term(q[7]) / constant_term(g_r), order)
    NGAx = Taylor1(constant_term(A2_t_g_r) * constant_term(tunitx), order)
    NGAy = Taylor1(constant_term(A2_t_g_r) * constant_term(tunity), order)
    NGAz = Taylor1(constant_term(A2_t_g_r) * constant_term(tunitz), order)
    tmp2785 = Taylor1(constant_term(postNewtonX) + constant_term(accX), order)
    dq[4] = Taylor1(constant_term(tmp2785) + constant_term(NGAx), order)
    tmp2787 = Taylor1(constant_term(postNewtonY) + constant_term(accY), order)
    dq[5] = Taylor1(constant_term(tmp2787) + constant_term(NGAy), order)
    tmp2789 = Taylor1(constant_term(postNewtonZ) + constant_term(accZ), order)
    dq[6] = Taylor1(constant_term(tmp2789) + constant_term(NGAz), order)
    dq[7] = Taylor1(identity(constant_term(zero_q_1)), order)
    for __idx = eachindex(q)
        (q[__idx]).coeffs[2] = (dq[__idx]).coeffs[1]
    end
    for ord = 1:order - 1
        ordnext = ord + 1
        TaylorSeries.identity!(pntempX, zero_q_1, ord)
        TaylorSeries.identity!(pntempY, zero_q_1, ord)
        TaylorSeries.identity!(pntempZ, zero_q_1, ord)
        TaylorSeries.identity!(accX, zero_q_1, ord)
        TaylorSeries.identity!(accY, zero_q_1, ord)
        TaylorSeries.identity!(accZ, zero_q_1, ord)
        TaylorSeries.identity!(dq[1], q[4], ord)
        TaylorSeries.identity!(dq[2], q[5], ord)
        TaylorSeries.identity!(dq[3], q[6], ord)
        TaylorSeries.identity!(newtonianNb_Potential[N], zero_q_1, ord)
        #= REPL[4]:161 =# Threads.@threads for i = _1_to_Nm1
                TaylorSeries.identity!(ui[i], ss16asteph_t[3 * ((N - 1) + i) - 2], ord)
                TaylorSeries.identity!(vi[i], ss16asteph_t[3 * ((N - 1) + i) - 1], ord)
                TaylorSeries.identity!(wi[i], ss16asteph_t[3 * ((N - 1) + i)], ord)
                TaylorSeries.subst!(X[i], ss16asteph_t[3i - 2], q[1], ord)
                TaylorSeries.subst!(Y[i], ss16asteph_t[3i - 1], q[2], ord)
                TaylorSeries.subst!(Z[i], ss16asteph_t[3i], q[3], ord)
                TaylorSeries.subst!(U[i], ui[i], dq[1], ord)
                TaylorSeries.subst!(V[i], vi[i], dq[2], ord)
                TaylorSeries.subst!(W[i], wi[i], dq[3], ord)
                TaylorSeries.mul!(tmp2543[1], 4, dq[1], ord)
                TaylorSeries.mul!(tmp2545[i], 3, ui[i], ord)
                TaylorSeries.subst!(_4U_m_3X[i], tmp2543[1], tmp2545[i], ord)
                TaylorSeries.mul!(tmp2548[2], 4, dq[2], ord)
                TaylorSeries.mul!(tmp2550[i], 3, vi[i], ord)
                TaylorSeries.subst!(_4V_m_3Y[i], tmp2548[2], tmp2550[i], ord)
                TaylorSeries.mul!(tmp2553[3], 4, dq[3], ord)
                TaylorSeries.mul!(tmp2555[i], 3, wi[i], ord)
                TaylorSeries.subst!(_4W_m_3Z[i], tmp2553[3], tmp2555[i], ord)
                TaylorSeries.mul!(pn2x[i], X[i], _4U_m_3X[i], ord)
                TaylorSeries.mul!(pn2y[i], Y[i], _4V_m_3Y[i], ord)
                TaylorSeries.mul!(pn2z[i], Z[i], _4W_m_3Z[i], ord)
                TaylorSeries.mul!(UU[i], ui[i], dq[1], ord)
                TaylorSeries.mul!(VV[i], vi[i], dq[2], ord)
                TaylorSeries.mul!(WW[i], wi[i], dq[3], ord)
                TaylorSeries.add!(tmp2563[i], UU[i], VV[i], ord)
                TaylorSeries.add!(vi_dot_vj[i], tmp2563[i], WW[i], ord)
                TaylorSeries.pow!(tmp2566[i], X[i], 2, ord)
                TaylorSeries.pow!(tmp2568[i], Y[i], 2, ord)
                TaylorSeries.add!(tmp2569[i], tmp2566[i], tmp2568[i], ord)
                TaylorSeries.pow!(tmp2571[i], Z[i], 2, ord)
                TaylorSeries.add!(r_p2[i], tmp2569[i], tmp2571[i], ord)
                TaylorSeries.sqrt!(r_p1d2[i], r_p2[i], ord)
                TaylorSeries.pow!(r_p3d2[i], r_p2[i], 1.5, ord)
                TaylorSeries.pow!(r_p7d2[i], r_p2[i], 3.5, ord)
                TaylorSeries.div!(newtonianCoeff[i], μ[i], r_p3d2[i], ord)
                TaylorSeries.add!(tmp2579[i], pn2x[i], pn2y[i], ord)
                TaylorSeries.add!(tmp2580[i], tmp2579[i], pn2z[i], ord)
                TaylorSeries.mul!(pn2[i], newtonianCoeff[i], tmp2580[i], ord)
                TaylorSeries.mul!(newton_acc_X[i], X[i], newtonianCoeff[i], ord)
                TaylorSeries.mul!(newton_acc_Y[i], Y[i], newtonianCoeff[i], ord)
                TaylorSeries.mul!(newton_acc_Z[i], Z[i], newtonianCoeff[i], ord)
                TaylorSeries.div!(newtonian1b_Potential[i], μ[i], r_p1d2[i], ord)
                TaylorSeries.mul!(pn3[i], 3.5, newtonian1b_Potential[i], ord)
                TaylorSeries.mul!(U_t_pn2[i], pn2[i], U[i], ord)
                TaylorSeries.mul!(V_t_pn2[i], pn2[i], V[i], ord)
                TaylorSeries.mul!(W_t_pn2[i], pn2[i], W[i], ord)
                if UJ_interaction[i]
                    TaylorSeries.subst!(tmp2591[i], X[i], ord)
                    TaylorSeries.mul!(t31[i], tmp2591[i], M_[1, 3, i], ord)
                    TaylorSeries.subst!(tmp2593[i], Y[i], ord)
                    TaylorSeries.mul!(t32[i], tmp2593[i], M_[2, 3, i], ord)
                    TaylorSeries.subst!(tmp2595[i], Z[i], ord)
                    TaylorSeries.mul!(t33[i], tmp2595[i], M_[3, 3, i], ord)
                    TaylorSeries.add!(tmp2597[i], t31[i], t32[i], ord)
                    TaylorSeries.add!(r_sin_ϕ[i], tmp2597[i], t33[i], ord)
                    TaylorSeries.div!(sin_ϕ[i], r_sin_ϕ[i], r_p1d2[i], ord)
                    TaylorSeries.asin!(ϕ[i], sin_ϕ[i], tmp2791[i], ord)
                    TaylorSeries.sincos!(tmp2792[i], cos_ϕ[i], ϕ[i], ord)
                    TaylorSeries.pow!(sin2_ϕ[i], sin_ϕ[i], 2, ord)
                    TaylorSeries.pow!(sin3_ϕ[i], sin_ϕ[i], 3, ord)
                    TaylorSeries.mul!(tmp2607[i], 1.5, sin2_ϕ[i], ord)
                    TaylorSeries.subst!(P_2_sin_ϕ[i], tmp2607[i], 0.5, ord)
                    TaylorSeries.mul!(∂P_2_sin_ϕ[i], 3, sin_ϕ[i], ord)
                    TaylorSeries.mul!(tmp2613[i], -1.5, sin_ϕ[i], ord)
                    TaylorSeries.mul!(tmp2615[i], 2.5, sin3_ϕ[i], ord)
                    TaylorSeries.add!(P_3_sin_ϕ[i], tmp2613[i], tmp2615[i], ord)
                    TaylorSeries.mul!(tmp2619[i], 7.5, sin2_ϕ[i], ord)
                    TaylorSeries.add!(∂P_3_sin_ϕ[i], -1.5, tmp2619[i], ord)
                    TaylorSeries.subst!(tmp2621[i], Λ2[i], ord)
                    TaylorSeries.pow!(tmp2623[i], r_p2[i], 2, ord)
                    TaylorSeries.div!(Λ2j_div_r4[i], tmp2621[i], tmp2623[i], ord)
                    TaylorSeries.subst!(tmp2625[i], Λ3[i], ord)
                    TaylorSeries.pow!(tmp2627[i], r_p1d2[i], 5, ord)
                    TaylorSeries.div!(Λ3j_div_r5[i], tmp2625[i], tmp2627[i], ord)
                    TaylorSeries.subst!(tmp2629[i], cos_ϕ[i], ord)
                    TaylorSeries.mul!(m_c_ϕ_∂P_2[i], tmp2629[i], ∂P_2_sin_ϕ[i], ord)
                    TaylorSeries.subst!(tmp2631[i], cos_ϕ[i], ord)
                    TaylorSeries.mul!(m_c_ϕ_∂P_3[i], tmp2631[i], ∂P_3_sin_ϕ[i], ord)
                    TaylorSeries.mul!(tmp2634[i], Λ2j_div_r4[i], 3, ord)
                    TaylorSeries.mul!(F_J2_ξ[i], tmp2634[i], P_2_sin_ϕ[i], ord)
                    TaylorSeries.mul!(F_J2_ζ[i], Λ2j_div_r4[i], m_c_ϕ_∂P_2[i], ord)
                    TaylorSeries.mul!(tmp2638[i], Λ3j_div_r5[i], 4, ord)
                    TaylorSeries.mul!(F_J3_ξ[i], tmp2638[i], P_3_sin_ϕ[i], ord)
                    TaylorSeries.mul!(F_J3_ζ[i], Λ3j_div_r5[i], m_c_ϕ_∂P_3[i], ord)
                    TaylorSeries.identity!(F_J_ξ[i], F_J2_ξ[i], ord)
                    TaylorSeries.identity!(F_J_ζ[i], F_J2_ζ[i], ord)
                    TaylorSeries.subst!(tmp2641[i], X[i], ord)
                    TaylorSeries.div!(ξx[i], tmp2641[i], r_p1d2[i], ord)
                    TaylorSeries.subst!(tmp2643[i], Y[i], ord)
                    TaylorSeries.div!(ξy[i], tmp2643[i], r_p1d2[i], ord)
                    TaylorSeries.subst!(tmp2645[i], Z[i], ord)
                    TaylorSeries.div!(ξz[i], tmp2645[i], r_p1d2[i], ord)
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
                end
                TaylorSeries.pow!(tmp2675[i], ui[i], 2, ord)
                TaylorSeries.pow!(tmp2677[i], vi[i], 2, ord)
                TaylorSeries.add!(tmp2678[i], tmp2675[i], tmp2677[i], ord)
                TaylorSeries.pow!(tmp2680[i], wi[i], 2, ord)
                TaylorSeries.add!(v2[i], tmp2678[i], tmp2680[i], ord)
            end
        TaylorSeries.pow!(tmp2683, q[4], 2, ord)
        TaylorSeries.pow!(tmp2685, q[5], 2, ord)
        TaylorSeries.add!(tmp2686, tmp2683, tmp2685, ord)
        TaylorSeries.pow!(tmp2688, q[6], 2, ord)
        TaylorSeries.add!(v2[N], tmp2686, tmp2688, ord)
        for i = _1_to_Nm1
            TaylorSeries.add!(temp_004[i], newtonian1b_Potential[i], newtonianNb_Potential[N], ord)
            TaylorSeries.identity!(newtonianNb_Potential[N], temp_004[i], ord)
            if UJ_interaction[i]
                TaylorSeries.mul!(tmp2691[i], μ[i], F_J2_x[i], ord)
                TaylorSeries.subst!(temp_accX_i[i], accX, tmp2691[i], ord)
                TaylorSeries.identity!(accX, temp_accX_i[i], ord)
                TaylorSeries.mul!(tmp2693[i], μ[i], F_J2_y[i], ord)
                TaylorSeries.subst!(temp_accY_i[i], accY, tmp2693[i], ord)
                TaylorSeries.identity!(accY, temp_accY_i[i], ord)
                TaylorSeries.mul!(tmp2695[i], μ[i], F_J2_z[i], ord)
                TaylorSeries.subst!(temp_accZ_i[i], accZ, tmp2695[i], ord)
                TaylorSeries.identity!(accZ, temp_accZ_i[i], ord)
            end
        end
        TaylorSeries.mul!(_4ϕj[N], 4, newtonianNb_Potential[N], ord)
        #= REPL[4]:306 =# Threads.@threads for i = 1:10
                TaylorSeries.add!(ϕi_plus_4ϕj[i], newtonianNb_Potential_t[i], _4ϕj[N], ord)
                TaylorSeries.mul!(tmp2701[i], 2, v2[i], ord)
                TaylorSeries.mul!(tmp2703[i], 4, vi_dot_vj[i], ord)
                TaylorSeries.subst!(tmp2704[i], tmp2701[i], tmp2703[i], ord)
                TaylorSeries.add!(sj2_plus_2si2_minus_4vivj[i], tmp2704[i], v2[N], ord)
                TaylorSeries.subst!(ϕs_and_vs[i], sj2_plus_2si2_minus_4vivj[i], ϕi_plus_4ϕj[i], ord)
                TaylorSeries.mul!(Xij_t_Ui[i], X[i], ui[i], ord)
                TaylorSeries.mul!(Yij_t_Vi[i], Y[i], vi[i], ord)
                TaylorSeries.mul!(Zij_t_Wi[i], Z[i], wi[i], ord)
                TaylorSeries.add!(tmp2710[i], Xij_t_Ui[i], Yij_t_Vi[i], ord)
                TaylorSeries.add!(Rij_dot_Vi[i], tmp2710[i], Zij_t_Wi[i], ord)
                TaylorSeries.pow!(tmp2713[i], Rij_dot_Vi[i], 2, ord)
                TaylorSeries.div!(pn1t7[i], tmp2713[i], r_p2[i], ord)
                TaylorSeries.mul!(tmp2716[i], 1.5, pn1t7[i], ord)
                TaylorSeries.subst!(pn1t2_7[i], ϕs_and_vs[i], tmp2716[i], ord)
                TaylorSeries.add!(pn1t1_7[i], c_p2, pn1t2_7[i], ord)
                TaylorSeries.mul!(pNX_t_X[i], acceph_t[3i - 2], X[i], ord)
                TaylorSeries.mul!(pNY_t_Y[i], acceph_t[3i - 1], Y[i], ord)
                TaylorSeries.mul!(pNZ_t_Z[i], acceph_t[3i], Z[i], ord)
                TaylorSeries.add!(tmp2723[i], pNX_t_X[i], pNY_t_Y[i], ord)
                TaylorSeries.add!(tmp2724[i], tmp2723[i], pNZ_t_Z[i], ord)
                TaylorSeries.mul!(tmp2725[i], 0.5, tmp2724[i], ord)
                TaylorSeries.add!(pn1[i], pn1t1_7[i], tmp2725[i], ord)
                TaylorSeries.mul!(X_t_pn1[i], newton_acc_X[i], pn1[i], ord)
                TaylorSeries.mul!(Y_t_pn1[i], newton_acc_Y[i], pn1[i], ord)
                TaylorSeries.mul!(Z_t_pn1[i], newton_acc_Z[i], pn1[i], ord)
                TaylorSeries.mul!(pNX_t_pn3[i], acceph_t[3i - 2], pn3[i], ord)
                TaylorSeries.mul!(pNY_t_pn3[i], acceph_t[3i - 1], pn3[i], ord)
                TaylorSeries.mul!(pNZ_t_pn3[i], acceph_t[3i], pn3[i], ord)
            end
        for i = 1:10
            TaylorSeries.add!(tmp2733[i], U_t_pn2[i], pNX_t_pn3[i], ord)
            TaylorSeries.add!(termpnx[i], X_t_pn1[i], tmp2733[i], ord)
            TaylorSeries.add!(sumpnx[i], pntempX, termpnx[i], ord)
            TaylorSeries.identity!(pntempX, sumpnx[i], ord)
            TaylorSeries.add!(tmp2736[i], V_t_pn2[i], pNY_t_pn3[i], ord)
            TaylorSeries.add!(termpny[i], Y_t_pn1[i], tmp2736[i], ord)
            TaylorSeries.add!(sumpny[i], pntempY, termpny[i], ord)
            TaylorSeries.identity!(pntempY, sumpny[i], ord)
            TaylorSeries.add!(tmp2739[i], W_t_pn2[i], pNZ_t_pn3[i], ord)
            TaylorSeries.add!(termpnz[i], Z_t_pn1[i], tmp2739[i], ord)
            TaylorSeries.add!(sumpnz[i], pntempZ, termpnz[i], ord)
            TaylorSeries.identity!(pntempZ, sumpnz[i], ord)
        end
        #= REPL[4]:345 =# Threads.@threads for i = 11:27
                TaylorSeries.mul!(X_t_pn1[i], c_p2, newton_acc_X[i], ord)
                TaylorSeries.mul!(Y_t_pn1[i], c_p2, newton_acc_Y[i], ord)
                TaylorSeries.mul!(Z_t_pn1[i], c_p2, newton_acc_Z[i], ord)
            end
        for i = 11:27
            TaylorSeries.identity!(termpnx[i], X_t_pn1[i], ord)
            TaylorSeries.add!(sumpnx[i], pntempX, termpnx[i], ord)
            TaylorSeries.identity!(pntempX, sumpnx[i], ord)
            TaylorSeries.identity!(termpny[i], Y_t_pn1[i], ord)
            TaylorSeries.add!(sumpny[i], pntempY, termpny[i], ord)
            TaylorSeries.identity!(pntempY, sumpny[i], ord)
            TaylorSeries.identity!(termpnz[i], Z_t_pn1[i], ord)
            TaylorSeries.add!(sumpnz[i], pntempZ, termpnz[i], ord)
            TaylorSeries.identity!(pntempZ, sumpnz[i], ord)
        end
        TaylorSeries.mul!(postNewtonX, pntempX, c_m2, ord)
        TaylorSeries.mul!(postNewtonY, pntempY, c_m2, ord)
        TaylorSeries.mul!(postNewtonZ, pntempZ, c_m2, ord)
        TaylorSeries.mul!(tmp2751, Y[1], W[1], ord)
        TaylorSeries.mul!(tmp2752, Z[1], V[1], ord)
        TaylorSeries.subst!(hx, tmp2751, tmp2752, ord)
        TaylorSeries.mul!(tmp2754, Z[1], U[1], ord)
        TaylorSeries.mul!(tmp2755, X[1], W[1], ord)
        TaylorSeries.subst!(hy, tmp2754, tmp2755, ord)
        TaylorSeries.mul!(tmp2757, X[1], V[1], ord)
        TaylorSeries.mul!(tmp2758, Y[1], U[1], ord)
        TaylorSeries.subst!(hz, tmp2757, tmp2758, ord)
        TaylorSeries.mul!(tmp2760, hz, Y[1], ord)
        TaylorSeries.mul!(tmp2761, hy, Z[1], ord)
        TaylorSeries.subst!(tunitx0, tmp2760, tmp2761, ord)
        TaylorSeries.mul!(tmp2763, hx, Z[1], ord)
        TaylorSeries.mul!(tmp2764, hz, X[1], ord)
        TaylorSeries.subst!(tunity0, tmp2763, tmp2764, ord)
        TaylorSeries.mul!(tmp2766, hy, X[1], ord)
        TaylorSeries.mul!(tmp2767, hx, Y[1], ord)
        TaylorSeries.subst!(tunitz0, tmp2766, tmp2767, ord)
        TaylorSeries.pow!(tmp2770, tunitx0, 2, ord)
        TaylorSeries.pow!(tmp2772, tunity0, 2, ord)
        TaylorSeries.add!(tmp2773, tmp2770, tmp2772, ord)
        TaylorSeries.pow!(tmp2775, tunitz0, 2, ord)
        TaylorSeries.add!(tmp2776, tmp2773, tmp2775, ord)
        TaylorSeries.sqrt!(hmag, tmp2776, ord)
        TaylorSeries.div!(tunitx, tunitx0, hmag, ord)
        TaylorSeries.div!(tunity, tunity0, hmag, ord)
        TaylorSeries.div!(tunitz, tunitz0, hmag, ord)
        TaylorSeries.identity!(g_r, r_p2[1], ord)
        TaylorSeries.div!(A2_t_g_r, q[7], g_r, ord)
        TaylorSeries.mul!(NGAx, A2_t_g_r, tunitx, ord)
        TaylorSeries.mul!(NGAy, A2_t_g_r, tunity, ord)
        TaylorSeries.mul!(NGAz, A2_t_g_r, tunitz, ord)
        TaylorSeries.add!(tmp2785, postNewtonX, accX, ord)
        TaylorSeries.add!(dq[4], tmp2785, NGAx, ord)
        TaylorSeries.add!(tmp2787, postNewtonY, accY, ord)
        TaylorSeries.add!(dq[5], tmp2787, NGAy, ord)
        TaylorSeries.add!(tmp2789, postNewtonZ, accZ, ord)
        TaylorSeries.add!(dq[6], tmp2789, NGAz, ord)
        TaylorSeries.identity!(dq[7], zero_q_1, ord)
        for __idx = eachindex(q)
            (q[__idx]).coeffs[ordnext + 1] = (dq[__idx]).coeffs[ordnext] / ordnext
        end
    end
    return nothing
end
