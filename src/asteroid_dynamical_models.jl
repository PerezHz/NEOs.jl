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
@taylorize function RNp1BP_pN_A_J23E_J2S_ng_eph!(dq, q, params, t)
    local ss16asteph_t = evaleph(params[1], t, q[1]) # params[2](t)*one(q[1]) #ss16asteph(t)
    local acceph_t = evaleph(params[2], t, q[1]) # params[2](t)*one(q[1]) #acc_eph(t)
    local newtonianNb_Potential_t = evaleph(params[3], t, q[1]) # params[3](t)*one(q[1]) #newtonianNb_Potential(t), massive bodies
    local S = eltype(q[1])
    local N = length(μ) # number of bodies, including NEA
    local _1_to_Nm1 = Base.OneTo(N-1) # iterator over all bodies

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
    local dsj2k = t-2.451545e6 # J2000.0 = 2.451545e6
    local αs = deg2rad(α_p_sun*one(t))
    local δs = deg2rad(δ_p_sun*one(t))
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

    #post-Newtonian corrections to gravitational acceleration
    #Moyer, 1971, page 7 eq. 35
    # post-Newtonian iterative procedure setup and initialization
    _4ϕj[N] = 4newtonianNb_Potential[N]
    for i in _1_to_Nm1
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
    for i in _1_to_Nm1
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
    postNewtonX = pntempX*c_m2
    postNewtonY = pntempY*c_m2
    postNewtonZ = pntempZ*c_m2

    # compute non-gravitational acceleration
    hx = (Y[1]*W[1])-(Z[1]*V[1])
    hy = (Z[1]*U[1])-(X[1]*W[1])
    hz = (X[1]*V[1])-(Y[1]*U[1])

    #cartesian components of transversal unit vector:
    tunitx0 = (hz*Y[1]) - (hy*Z[1])
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

@taylorize function RNp1BP_pN_A_J23E_J2S_ng_eph_threads!(dq, q, params, t)
    local ss16asteph_t = evaleph(params[1], t, q[1]) # params[2](t)*one(q[1]) #ss16asteph(t)
    local acceph_t = evaleph(params[2], t, q[1]) # params[2](t)*one(q[1]) #acc_eph(t)
    local newtonianNb_Potential_t = evaleph(params[3], t, q[1]) # params[3](t)*one(q[1]) #newtonianNb_Potential(t), massive bodies
    local S = eltype(q[1])
    local N = length(μ) # number of bodies, including NEA
    local _1_to_Nm1 = Base.OneTo(N-1) # iterator over all bodies

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
    local dsj2k = t-2.451545e6 # J2000.0 = 2.451545e6
    local αs = deg2rad(α_p_sun*one(t))
    local δs = deg2rad(δ_p_sun*one(t))
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

    #post-Newtonian corrections to gravitational acceleration
    #Moyer, 1971, page 7 eq. 35
    # post-Newtonian iterative procedure setup and initialization
    _4ϕj[N] = 4newtonianNb_Potential[N]
    Threads.@threads for i in _1_to_Nm1
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
    for i in _1_to_Nm1
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
