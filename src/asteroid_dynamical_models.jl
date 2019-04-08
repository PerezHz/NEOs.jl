# Nearth-Earth asteroid dynamical model (d=2.0)
# Bodies considered in the model are: the Sun, the eight planets, the Moon and Ceres,
# as well as the asteroid of interest as a test particle with null mass. Dynamical
# effects considered are:
# - post-Newtonian point-mass accelerations between all bodies,
# - figure-effects (oblateness) of the Earth (J2, J3 and J4)
# - J2 effect of the Sun
# - Kinematic model for the precession and nutation of the Earth's pole (DE430/431 model)
# - also, a model for non-gravitational accelerations acting upon the asteroid
# is included (Yarkovsky effect) a_nongrav = A2*t_vec*(au/r)^d, where t_vec is the
# unit heliocentric transverse vector, au is 1 astronomical unit, r is the
# asteroid's heliocentric range, A2 is a coefficient (with units of au/day^2),
# and d = 2.0
@taylorize function RNp1BP_pN_A_J234E_J2S_ng!(t, q, dq)
    local S = eltype(q[1])
    local N = Int((length(q)-1)/6) # number of bodies, including NEA
    local _1_to_N = Base.OneTo(N) # iterator over all bodies

    local succ_approx_iter = 1 # number of iterations of post-Newtonian subroutine
    local su = 1 #Sun's index within `system`
    local ea = 4 #Earth's index within `system`
    local j2_body_index = [su, ea] # indices of bodies with J2 flattening (note: Earth also has J3 and J4)

    # parameters related to speed of light, c
    local c_p2 = 29979.063823897606 # c^2 = 29979.063823897606 au^2/d^2
    local c_m2 = 3.3356611996764786e-5 # c^-2 = 3.3356611996764786e-5 d^2/au^2

    local zero_q_1 = zero(q[1])

    X = Array{Taylor1{S}}(undef, N, N)
    Y = Array{Taylor1{S}}(undef, N, N)
    Z = Array{Taylor1{S}}(undef, N, N)

    r_p2 = Array{Taylor1{S}}(undef, N, N)
    r_p3d2 = Array{Taylor1{S}}(undef, N, N)
    r_p7d2 = Array{Taylor1{S}}(undef, N, N)

    #Newtonian accelerations
    newtonX = Array{Taylor1{S}}(undef, N)
    newtonY = Array{Taylor1{S}}(undef, N)
    newtonZ = Array{Taylor1{S}}(undef, N)

    newtonianCoeff = Array{Taylor1{S}}(undef, N, N)

    #post-Newtonian stuff
    U = Array{Taylor1{S}}(undef, N, N)
    V = Array{Taylor1{S}}(undef, N, N)
    W = Array{Taylor1{S}}(undef, N, N)

    _4U_m_3X = Array{Taylor1{S}}(undef, N, N)
    _4V_m_3Y = Array{Taylor1{S}}(undef, N, N)
    _4W_m_3Z = Array{Taylor1{S}}(undef, N, N)

    UU = Array{Taylor1{S}}(undef, N, N)
    VV = Array{Taylor1{S}}(undef, N, N)
    WW = Array{Taylor1{S}}(undef, N, N)

    r_p1d2 = Array{Taylor1{S}}(undef, N, N)

    postNewtonX = Array{Taylor1{S}}(undef, N)
    postNewtonY = Array{Taylor1{S}}(undef, N)
    postNewtonZ = Array{Taylor1{S}}(undef, N)

    newtonianNb_Potential = Array{Taylor1{S}}(undef, N)
    newtonian1b_Potential = Array{Taylor1{S}}(undef, N, N)
    newtonianCoeff = Array{Taylor1{S}}(undef, N, N)

    pntempX = Array{Taylor1{S}}(undef, N)
    pntempY = Array{Taylor1{S}}(undef, N)
    pntempZ = Array{Taylor1{S}}(undef, N)

    pn1 = Array{Taylor1{S}}(undef, N, N)
    v2 = Array{Taylor1{S}}(undef, N)
    vi_dot_vj = Array{Taylor1{S}}(undef, N, N)
    pn2 = Array{Taylor1{S}}(undef, N, N)
    pn3 = Array{Taylor1{S}}(undef, N, N)

    # J2 acceleration auxiliaries
    t11 = Array{Taylor1{S}}(undef, N, N)
    t12 = Array{Taylor1{S}}(undef, N, N)
    t13 = Array{Taylor1{S}}(undef, N, N)
    t21 = Array{Taylor1{S}}(undef, N, N)
    t22 = Array{Taylor1{S}}(undef, N, N)
    t23 = Array{Taylor1{S}}(undef, N, N)
    t31 = Array{Taylor1{S}}(undef, N, N)
    t32 = Array{Taylor1{S}}(undef, N, N)
    t33 = Array{Taylor1{S}}(undef, N, N)
    new_x = Array{Taylor1{S}}(undef, N, N)
    new_y = Array{Taylor1{S}}(undef, N, N)
    new_z = Array{Taylor1{S}}(undef, N, N)
    new_x2 = Array{Taylor1{S}}(undef, N, N)
    new_y2 = Array{Taylor1{S}}(undef, N, N)
    ρ_p2_ij = Array{Taylor1{S}}(undef, N, N)
    z_p2_ij = Array{Taylor1{S}}(undef, N, N)
    temp_ρ = Array{Taylor1{S}}(undef, N, N)
    temp_z = Array{Taylor1{S}}(undef, N, N)
    dum00 = Array{Taylor1{S}}(undef, N, N)
    dum01 = Array{Taylor1{S}}(undef, N, N)
    dum = Array{Taylor1{S}}(undef, N, N)
    dum_ρ = Array{Taylor1{S}}(undef, N, N)
    dum_z = Array{Taylor1{S}}(undef, N, N)
    F_J2_bf_x = Array{Taylor1{S}}(undef, N, N)
    F_J2_bf_y = Array{Taylor1{S}}(undef, N, N)
    F_J2_bf_z = Array{Taylor1{S}}(undef, N, N)
    s11 = Array{Taylor1{S}}(undef, N, N)
    s12 = Array{Taylor1{S}}(undef, N, N)
    s13 = Array{Taylor1{S}}(undef, N, N)
    s21 = Array{Taylor1{S}}(undef, N, N)
    s22 = Array{Taylor1{S}}(undef, N, N)
    s23 = Array{Taylor1{S}}(undef, N, N)
    s31 = Array{Taylor1{S}}(undef, N, N)
    s32 = Array{Taylor1{S}}(undef, N, N)
    s33 = Array{Taylor1{S}}(undef, N, N)
    F_J2_x = Array{Taylor1{S}}(undef, N, N)
    F_J2_y = Array{Taylor1{S}}(undef, N, N)
    F_J2_z = Array{Taylor1{S}}(undef, N, N)
    timp_004 = Array{Taylor1{S}}(undef, N, N)
    timp_005 = Array{Taylor1{S}}(undef, N, N)
    timp_006 = Array{Taylor1{S}}(undef, N, N)
    tamp_004 = Array{Taylor1{S}}(undef, N, N)
    tamp_005 = Array{Taylor1{S}}(undef, N, N)
    tamp_006 = Array{Taylor1{S}}(undef, N, N)

    # extended-body accelerations
    accX = Array{Taylor1{S}}(undef, N)
    accY = Array{Taylor1{S}}(undef, N)
    accZ = Array{Taylor1{S}}(undef, N)

    # rotations to and from Earth and Sun poles
    local αs = deg2rad(268.13+0t)
    local δs = deg2rad(63.87+0t)
    local M_ = Array{Taylor1{S}}(undef, 3, 3, N)
    local W_ = Array{Taylor1{S}}(undef, 3, 3, N)
    local M_[:,:,ea] = earth_pole_rotation(t-J2000) # J2000.0 = 2.451545e6
    local W_[:,:,ea] = inv(M_[:, :, ea])
    local M_[:,:,su] = pole_rotation( αs, δs )
    local W_[:,:,su] = inv(M_[:, :, su])

    for j in _1_to_N
        newtonX[j] = zero_q_1
        newtonY[j] = zero_q_1
        newtonZ[j] = zero_q_1

        newtonianNb_Potential[j] = zero_q_1

        accX[j] = zero_q_1
        accY[j] = zero_q_1
        accZ[j] = zero_q_1

        dq[3j-2] = q[3(N+j)-2]
        dq[3j-1] = q[3(N+j)-1]
        dq[3j  ] = q[3(N+j)  ]
    end

    for j in j2_body_index
        for i in _1_to_N
            if i == j
            else
                t11[i,j] = zero_q_1
                t12[i,j] = zero_q_1
                t13[i,j] = zero_q_1
                t21[i,j] = zero_q_1
                t22[i,j] = zero_q_1
                t23[i,j] = zero_q_1
                t31[i,j] = zero_q_1
                t32[i,j] = zero_q_1
                t33[i,j] = zero_q_1
                new_x[i,j] = zero_q_1
                new_x2[i,j] = zero_q_1
                new_y2[i,j] = zero_q_1
                ρ_p2_ij[i,j] = zero_q_1
                z_p2_ij[i,j] = zero_q_1
                temp_ρ[i,j] = zero_q_1
                temp_z[i,j] = zero_q_1
                dum00[i,j] = zero_q_1
                dum01[i,j] = zero_q_1
                dum[i,j] = zero_q_1
                dum_ρ[i,j] = zero_q_1
                dum_z[i,j] = zero_q_1
                F_J2_bf_x[i,j] = zero_q_1
                F_J2_bf_y[i,j] = zero_q_1
                F_J2_bf_z[i,j] = zero_q_1
                s11[i,j] = zero_q_1
                s12[i,j] = zero_q_1
                s13[i,j] = zero_q_1
                s21[i,j] = zero_q_1
                s22[i,j] = zero_q_1
                s23[i,j] = zero_q_1
                s31[i,j] = zero_q_1
                s32[i,j] = zero_q_1
                s33[i,j] = zero_q_1
                F_J2_x[i,j] = zero_q_1
                F_J2_y[i,j] = zero_q_1
                F_J2_z[i,j] = zero_q_1
                timp_004[i,j] = zero_q_1
                timp_005[i,j] = zero_q_1
                timp_006[i,j] = zero_q_1
                tamp_004[i,j] = zero_q_1
                tamp_005[i,j] = zero_q_1
                tamp_006[i,j] = zero_q_1
            end #if i == j
        end #for i in _1_to_N
    end #for j in j2_body_index

    #compute point-mass Newtonian accelerations, all bodies
    for j in _1_to_N
        for i in _1_to_N
            # i == j && continue
            if i == j
            else
                X[i,j] = q[3i-2]-q[3j-2]
                Y[i,j] = q[3i-1]-q[3j-1]
                Z[i,j] = q[3i]-q[3j]

                U[i,j] = dq[3i-2]-dq[3j-2]
                V[i,j] = dq[3i-1]-dq[3j-1]
                W[i,j] = dq[3i  ]-dq[3j  ]

                _4U_m_3X[i,j] = (4dq[3j-2])-(3dq[3i-2])
                _4V_m_3Y[i,j] = (4dq[3j-1])-(3dq[3i-1])
                _4W_m_3Z[i,j] = (4dq[3j  ])-(3dq[3i  ])

                pn2x = X[i,j]*_4U_m_3X[i,j]
                pn2y = Y[i,j]*_4V_m_3Y[i,j]
                pn2z = Z[i,j]*_4W_m_3Z[i,j]

                pn2[i,j] = ( pn2x+pn2y ) + pn2z

                UU[i,j] = dq[3i-2]*dq[3j-2]
                VV[i,j] = dq[3i-1]*dq[3j-1]
                WW[i,j] = dq[3i  ]*dq[3j  ]

                vi_dot_vj[i,j] = ( UU[i,j]+VV[i,j] ) + WW[i,j]

                r_p2[i,j] = ( (X[i,j]^2)+(Y[i,j]^2) ) + (Z[i,j]^2)

                r_p1d2[i,j] = sqrt(r_p2[i,j])
                r_p3d2[i,j] = r_p2[i,j]^1.5
                r_p7d2[i,j] = r_p2[i,j]^3.5

                newtonianCoeff[i,j] =  μ[i]/r_p3d2[i,j]

                temp_001 = newtonX[j] + (X[i,j]*newtonianCoeff[i,j])
                newtonX[j] = temp_001
                temp_002 = newtonY[j] + (Y[i,j]*newtonianCoeff[i,j])
                newtonY[j] = temp_002
                temp_003 = newtonZ[j] + (Z[i,j]*newtonianCoeff[i,j])
                newtonZ[j] = temp_003

                newtonian1b_Potential[i,j] = μ[i]/r_p1d2[i, j]

                temp_004 = newtonianNb_Potential[j] + newtonian1b_Potential[i, j]
                newtonianNb_Potential[j] = temp_004
            end #if i != j
        end #for, i
        v2[j] = ( (dq[3j-2]^2)+(dq[3j-1]^2) ) + (dq[3j]^2)
    end #for, j

    for j in _1_to_N
        postNewtonX[j] = newtonX[j]
        postNewtonY[j] = newtonY[j]
        postNewtonZ[j] = newtonZ[j]
    end

    for k in Base.OneTo(succ_approx_iter)
        for j in _1_to_N
            pntempX[j] = zero_q_1
            pntempY[j] = zero_q_1
            pntempZ[j] = zero_q_1
        end
        for j in _1_to_N
            for i in _1_to_N
                # i == j && continue
                if i == j
                else
                    #post-Newtonian corrections to gravitational acceleration
                    #Moyer, 1971, page 7 eq. 35
                    temp_005a = newtonianNb_Potential[i]+(4newtonianNb_Potential[j])
                    temp_005b = (2v2[i])-(4vi_dot_vj[i,j])
                    temp_005c = v2[j]+temp_005b
                    temp_005 = temp_005c-temp_005a
                    temp_006a = X[i,j]*dq[3i-2]
                    temp_006b = Y[i,j]*dq[3i-1]
                    temp_006c = Z[i,j]*dq[3i]
                    temp_006d = ( temp_006a+temp_006b ) + temp_006c
                    # the expression below inside the (...)^2 should have a minus sign in front of the numerator,
                    # but upon squaring it is eliminated, so at the end of the day, it is irrelevant ;)
                    temp_006e = (temp_006d^2)/r_p2[i,j]
                    temp_006 = temp_005-(1.5temp_006e)
                    temp_007a = X[i,j]*postNewtonX[i]
                    temp_007b = Y[i,j]*postNewtonY[i]
                    temp_007c = Z[i,j]*postNewtonZ[i]
                    temp_007d = ( temp_007a+temp_007b ) + temp_007c
                    temp_007 = temp_006 + (0.5temp_007d)
                    temp_008 = c_p2+temp_007
                    pn1[i,j] = newtonianCoeff[i,j]*temp_008

                    temp_009 = pntempX[j]+(X[i,j]*pn1[i,j])
                    temp_010 = pntempY[j]+(Y[i,j]*pn1[i,j])
                    temp_011 = pntempZ[j]+(Z[i,j]*pn1[i,j])

                    pn3[i,j] = 3.5*newtonian1b_Potential[i,j]

                    temp_013a = pn2[i,j]*(U[i,j]*newtonianCoeff[i,j])
                    temp_013b = postNewtonX[i]*pn3[i,j]
                    temp_013 = temp_009 + (temp_013a+temp_013b)
                    pntempX[j] = temp_013

                    temp_014a = pn2[i,j]*(V[i,j]*newtonianCoeff[i,j])
                    temp_014b = postNewtonY[i]*pn3[i,j]
                    temp_014 = temp_010 + (temp_014a+temp_014b)
                    pntempY[j] = temp_014

                    temp_015a = pn2[i,j]*(W[i,j]*newtonianCoeff[i,j])
                    temp_015b = postNewtonZ[i]*pn3[i,j]
                    temp_015 = temp_011 + (temp_015a+temp_015b)
                    pntempZ[j] = temp_015
                end
            end #for i
        end #for j
        for j in _1_to_N
            postNewtonX[j] = pntempX[j]*c_m2
            postNewtonY[j] = pntempY[j]*c_m2
            postNewtonZ[j] = pntempZ[j]*c_m2
        end
    end #for k in Base.OneTo(succ_approx_iter) # (post-Newtonian iterations)

    #J2 accelerations, all flattened bodies
    for j in j2_body_index
        for i in _1_to_N
            if i == j
            else
                # # rotate from inertial frame to extended-body frame
                t11[i,j] = X[i,j]*W_[1,1,j]
                t21[i,j] = X[i,j]*W_[2,1,j]
                t31[i,j] = X[i,j]*W_[3,1,j]
                t12[i,j] = Y[i,j]*W_[1,2,j]
                t22[i,j] = Y[i,j]*W_[2,2,j]
                t32[i,j] = Y[i,j]*W_[3,2,j]
                t13[i,j] = Z[i,j]*W_[1,3,j]
                t23[i,j] = Z[i,j]*W_[2,3,j]
                t33[i,j] = Z[i,j]*W_[3,3,j]
                new_x[i,j] = (t11[i,j]+t12[i,j])+t13[i,j]
                new_y[i,j] = (t21[i,j]+t22[i,j])+t23[i,j]
                new_z[i,j] = (t31[i,j]+t32[i,j])+t33[i,j]

                # # compute cartesian components of extended-body acceleration in body frame
                new_x2[i,j] = new_x[i,j]^2
                new_y2[i,j] = new_y[i,j]^2

                ρ_p2_ij[i,j] = new_x2[i,j]+new_y2[i,j]
                z_p2_ij[i,j] = new_z[i,j]^2

                temp_ρ[i,j] = ρ_p2_ij[i,j] - (4z_p2_ij[i,j])
                temp_z[i,j] = (3ρ_p2_ij[i,j]) - (2z_p2_ij[i,j])

                dum01[i,j] = Λ2[j]/r_p7d2[i,j]

                dum[i,j] = 1.5*dum01[i,j]
                dum_ρ[i,j] = dum[i,j]*temp_ρ[i,j]
                dum_z[i,j] = dum[i,j]*temp_z[i,j]

                F_J2_bf_x[i,j] = dum_ρ[i,j]*new_x[i,j]
                F_J2_bf_y[i,j] = dum_ρ[i,j]*new_y[i,j]
                F_J2_bf_z[i,j] = dum_z[i,j]*new_z[i,j]

                # # rotate components of force from body frame to inertial frame
                s11[i,j] = F_J2_bf_x[i,j]*M_[1,1,j]
                s21[i,j] = F_J2_bf_x[i,j]*M_[2,1,j]
                s31[i,j] = F_J2_bf_x[i,j]*M_[3,1,j]
                s12[i,j] = F_J2_bf_y[i,j]*M_[1,2,j]
                s22[i,j] = F_J2_bf_y[i,j]*M_[2,2,j]
                s32[i,j] = F_J2_bf_y[i,j]*M_[3,2,j]
                s13[i,j] = F_J2_bf_z[i,j]*M_[1,3,j]
                s23[i,j] = F_J2_bf_z[i,j]*M_[2,3,j]
                s33[i,j] = F_J2_bf_z[i,j]*M_[3,3,j]
                F_J2_x[i,j] = (s11[i,j]+s12[i,j])+s13[i,j]
                F_J2_y[i,j] = (s21[i,j]+s22[i,j])+s23[i,j]
                F_J2_z[i,j] = (s31[i,j]+s32[i,j])+s33[i,j]

                # # add result to total acceleration upon j-th body figure due to i-th point mass
                # @show "acc",j,"+μ",i,"Λ2",j
                timp_004[i,j] = accX[j] + (μ[i]*F_J2_x[i,j])
                accX[j] = timp_004[i,j]
                timp_005[i,j] = accY[j] + (μ[i]*F_J2_y[i,j])
                accY[j] = timp_005[i,j]
                timp_006[i,j] = accZ[j] + (μ[i]*F_J2_z[i,j])
                accZ[j] = timp_006[i,j]

                # # reaction force on i-th body
                # @show "acc",i,"-μ",j,"Λ2",j
                tamp_004[i,j] = accX[i] - (μ[j]*F_J2_x[i,j])
                accX[i] = tamp_004[i,j]
                tamp_005[i,j] = accY[i] - (μ[j]*F_J2_y[i,j])
                accY[i] = tamp_005[i,j]
                tamp_006[i,j] = accZ[i] - (μ[j]*F_J2_z[i,j])
                accZ[i] = tamp_006[i,j]
            end #if i == j
        end #for i in _1_to_N
    end #for j in j2_body_index

    #fill the equations of motion for everyone except test particle (Newtonian, post-Newtonian and extended body accelerations)
    for i in Base.OneTo(N-1)
        dq[3(N+i)-2] = postNewtonX[i]+accX[i]
        dq[3(N+i)-1] = postNewtonY[i]+accY[i]
        dq[3(N+i)  ] = postNewtonZ[i]+accZ[i]
    end

    #computation of non-gravitational accelerations:
    hx = (Y[N,1]*(dq[3N  ]-dq[3]))-(Z[N,1]*(dq[3N-1]-dq[2]))
    hy = (Z[N,1]*(dq[3N-2]-dq[1]))-(X[N,1]*(dq[3N  ]-dq[3]))
    hz = (X[N,1]*(dq[3N-1]-dq[2]))-(Y[N,1]*(dq[3N-2]-dq[1]))
    r_hs = sqrt(r_p2[N,1])
    runitx = X[N,1]/r_hs
    runity = Y[N,2]/r_hs
    runitz = Z[N,3]/r_hs

    #cartesian components of transversal unit vector:
    tunitx0 = (hy*runitz)-(hz*runity)
    tunity0 = (hz*runitx)-(hx*runitz)
    tunitz0 = (hx*runity)-(hy*runitx)
    hmag = sqrt( ((tunitx0^2)+(tunity0^2))+(tunitz0^2) )
    tunitx = tunitx0/hmag
    tunity = tunity0/hmag
    tunitz = tunitz0/hmag

    # evaluate non-grav acceleration of NEA (Yarkovsky):
    g_r = r_hs^(-2.0)
    A2_t_g_r = q[6N+1]*g_r

    NGAx = A2_t_g_r*tunitx
    NGAy = A2_t_g_r*tunity
    NGAz = A2_t_g_r*tunitz

    dq[6N-2] = (postNewtonX[N]+accX[N])+NGAx
    dq[6N-1] = (postNewtonY[N]+accY[N])+NGAy
    dq[6N  ] = (postNewtonZ[N]+accZ[N])+NGAz

    dq[6N+1] = zero_q_1

    nothing
end

# Nearth-Earth asteroid dynamical model (d=2.25)
# Bodies considered in the model are: the Sun, the eight planets, the Moon and Ceres,
# as well as the asteroid of interest as a test particle with null mass. Dynamical
# effects considered are:
# - post-Newtonian point-mass accelerations between all bodies,
# - figure-effects (oblateness) of the Earth (J2, J3 and J4)
# - J2 effect of the Sun
# - Kinematic model for the precession and nutation of the Earth's pole (DE430/431 model)
# - also, a model for non-gravitational accelerations acting upon the asteroid
# is included (Yarkovsky effect) a_nongrav = A2*t_vec*(au/r)^d, where t_vec is the
# unit heliocentric transverse vector, au is 1 astronomical unit, r is the
# asteroid's heliocentric range, A2 is a coefficient (with units of au/day^2),
# and d = 2.25
@taylorize function RNp1BP_pN_A_J234E_J2S_ng_d225!(t, q, dq)
    local S = eltype(q[1])
    local N = Int((length(q)-1)/6) # number of bodies, including NEA
    local _1_to_N = Base.OneTo(N) # iterator over all bodies

    local succ_approx_iter = 1 # number of iterations of post-Newtonian subroutine
    local su = 1 #Sun's index within `system`
    local ea = 4 #Earth's index within `system`
    local j2_body_index = [su, ea] # indices of bodies with J2 flattening (note: Earth also has J3 and J4)

    # parameters related to speed of light, c
    local c_p2 = 29979.063823897606 # c^2 = 29979.063823897606 au^2/d^2
    local c_m2 = 3.3356611996764786e-5 # c^-2 = 3.3356611996764786e-5 d^2/au^2

    local zero_q_1 = zero(q[1])

    X = Array{Taylor1{S}}(undef, N, N)
    Y = Array{Taylor1{S}}(undef, N, N)
    Z = Array{Taylor1{S}}(undef, N, N)

    r_p2 = Array{Taylor1{S}}(undef, N, N)
    r_p3d2 = Array{Taylor1{S}}(undef, N, N)
    r_p7d2 = Array{Taylor1{S}}(undef, N, N)

    #Newtonian accelerations
    newtonX = Array{Taylor1{S}}(undef, N)
    newtonY = Array{Taylor1{S}}(undef, N)
    newtonZ = Array{Taylor1{S}}(undef, N)

    newtonianCoeff = Array{Taylor1{S}}(undef, N, N)

    #post-Newtonian stuff
    U = Array{Taylor1{S}}(undef, N, N)
    V = Array{Taylor1{S}}(undef, N, N)
    W = Array{Taylor1{S}}(undef, N, N)

    _4U_m_3X = Array{Taylor1{S}}(undef, N, N)
    _4V_m_3Y = Array{Taylor1{S}}(undef, N, N)
    _4W_m_3Z = Array{Taylor1{S}}(undef, N, N)

    UU = Array{Taylor1{S}}(undef, N, N)
    VV = Array{Taylor1{S}}(undef, N, N)
    WW = Array{Taylor1{S}}(undef, N, N)

    r_p1d2 = Array{Taylor1{S}}(undef, N, N)

    postNewtonX = Array{Taylor1{S}}(undef, N)
    postNewtonY = Array{Taylor1{S}}(undef, N)
    postNewtonZ = Array{Taylor1{S}}(undef, N)

    newtonianNb_Potential = Array{Taylor1{S}}(undef, N)
    newtonian1b_Potential = Array{Taylor1{S}}(undef, N, N)
    newtonianCoeff = Array{Taylor1{S}}(undef, N, N)

    pntempX = Array{Taylor1{S}}(undef, N)
    pntempY = Array{Taylor1{S}}(undef, N)
    pntempZ = Array{Taylor1{S}}(undef, N)

    pn1 = Array{Taylor1{S}}(undef, N, N)
    v2 = Array{Taylor1{S}}(undef, N)
    vi_dot_vj = Array{Taylor1{S}}(undef, N, N)
    pn2 = Array{Taylor1{S}}(undef, N, N)
    pn3 = Array{Taylor1{S}}(undef, N, N)

    # J2 acceleration auxiliaries
    t11 = Array{Taylor1{S}}(undef, N, N)
    t12 = Array{Taylor1{S}}(undef, N, N)
    t13 = Array{Taylor1{S}}(undef, N, N)
    t21 = Array{Taylor1{S}}(undef, N, N)
    t22 = Array{Taylor1{S}}(undef, N, N)
    t23 = Array{Taylor1{S}}(undef, N, N)
    t31 = Array{Taylor1{S}}(undef, N, N)
    t32 = Array{Taylor1{S}}(undef, N, N)
    t33 = Array{Taylor1{S}}(undef, N, N)
    new_x = Array{Taylor1{S}}(undef, N, N)
    new_y = Array{Taylor1{S}}(undef, N, N)
    new_z = Array{Taylor1{S}}(undef, N, N)
    new_x2 = Array{Taylor1{S}}(undef, N, N)
    new_y2 = Array{Taylor1{S}}(undef, N, N)
    ρ_p2_ij = Array{Taylor1{S}}(undef, N, N)
    z_p2_ij = Array{Taylor1{S}}(undef, N, N)
    temp_ρ = Array{Taylor1{S}}(undef, N, N)
    temp_z = Array{Taylor1{S}}(undef, N, N)
    dum00 = Array{Taylor1{S}}(undef, N, N)
    dum01 = Array{Taylor1{S}}(undef, N, N)
    dum = Array{Taylor1{S}}(undef, N, N)
    dum_ρ = Array{Taylor1{S}}(undef, N, N)
    dum_z = Array{Taylor1{S}}(undef, N, N)
    F_J2_bf_x = Array{Taylor1{S}}(undef, N, N)
    F_J2_bf_y = Array{Taylor1{S}}(undef, N, N)
    F_J2_bf_z = Array{Taylor1{S}}(undef, N, N)
    s11 = Array{Taylor1{S}}(undef, N, N)
    s12 = Array{Taylor1{S}}(undef, N, N)
    s13 = Array{Taylor1{S}}(undef, N, N)
    s21 = Array{Taylor1{S}}(undef, N, N)
    s22 = Array{Taylor1{S}}(undef, N, N)
    s23 = Array{Taylor1{S}}(undef, N, N)
    s31 = Array{Taylor1{S}}(undef, N, N)
    s32 = Array{Taylor1{S}}(undef, N, N)
    s33 = Array{Taylor1{S}}(undef, N, N)
    F_J2_x = Array{Taylor1{S}}(undef, N, N)
    F_J2_y = Array{Taylor1{S}}(undef, N, N)
    F_J2_z = Array{Taylor1{S}}(undef, N, N)
    timp_004 = Array{Taylor1{S}}(undef, N, N)
    timp_005 = Array{Taylor1{S}}(undef, N, N)
    timp_006 = Array{Taylor1{S}}(undef, N, N)
    tamp_004 = Array{Taylor1{S}}(undef, N, N)
    tamp_005 = Array{Taylor1{S}}(undef, N, N)
    tamp_006 = Array{Taylor1{S}}(undef, N, N)

    # extended-body accelerations
    accX = Array{Taylor1{S}}(undef, N)
    accY = Array{Taylor1{S}}(undef, N)
    accZ = Array{Taylor1{S}}(undef, N)

    # rotations to and from Earth and Sun poles
    local αs = deg2rad(268.13+0t)
    local δs = deg2rad(63.87+0t)
    local M_ = Array{Taylor1{S}}(undef, 3, 3, N)
    local W_ = Array{Taylor1{S}}(undef, 3, 3, N)
    local M_[:,:,ea] = earth_pole_rotation(t-J2000) # J2000.0 = 2.451545e6
    local W_[:,:,ea] = inv(M_[:, :, ea])
    local M_[:,:,su] = pole_rotation( αs, δs )
    local W_[:,:,su] = inv(M_[:, :, su])

    for j in _1_to_N
        newtonX[j] = zero_q_1
        newtonY[j] = zero_q_1
        newtonZ[j] = zero_q_1

        newtonianNb_Potential[j] = zero_q_1

        accX[j] = zero_q_1
        accY[j] = zero_q_1
        accZ[j] = zero_q_1

        dq[3j-2] = q[3(N+j)-2]
        dq[3j-1] = q[3(N+j)-1]
        dq[3j  ] = q[3(N+j)  ]
    end

    for j in j2_body_index
        for i in _1_to_N
            if i == j
            else
                t11[i,j] = zero_q_1
                t12[i,j] = zero_q_1
                t13[i,j] = zero_q_1
                t21[i,j] = zero_q_1
                t22[i,j] = zero_q_1
                t23[i,j] = zero_q_1
                t31[i,j] = zero_q_1
                t32[i,j] = zero_q_1
                t33[i,j] = zero_q_1
                new_x[i,j] = zero_q_1
                new_x2[i,j] = zero_q_1
                new_y2[i,j] = zero_q_1
                ρ_p2_ij[i,j] = zero_q_1
                z_p2_ij[i,j] = zero_q_1
                temp_ρ[i,j] = zero_q_1
                temp_z[i,j] = zero_q_1
                dum00[i,j] = zero_q_1
                dum01[i,j] = zero_q_1
                dum[i,j] = zero_q_1
                dum_ρ[i,j] = zero_q_1
                dum_z[i,j] = zero_q_1
                F_J2_bf_x[i,j] = zero_q_1
                F_J2_bf_y[i,j] = zero_q_1
                F_J2_bf_z[i,j] = zero_q_1
                s11[i,j] = zero_q_1
                s12[i,j] = zero_q_1
                s13[i,j] = zero_q_1
                s21[i,j] = zero_q_1
                s22[i,j] = zero_q_1
                s23[i,j] = zero_q_1
                s31[i,j] = zero_q_1
                s32[i,j] = zero_q_1
                s33[i,j] = zero_q_1
                F_J2_x[i,j] = zero_q_1
                F_J2_y[i,j] = zero_q_1
                F_J2_z[i,j] = zero_q_1
                timp_004[i,j] = zero_q_1
                timp_005[i,j] = zero_q_1
                timp_006[i,j] = zero_q_1
                tamp_004[i,j] = zero_q_1
                tamp_005[i,j] = zero_q_1
                tamp_006[i,j] = zero_q_1
            end #if i == j
        end #for i in _1_to_N
    end #for j in j2_body_index

    #compute point-mass Newtonian accelerations, all bodies
    for j in _1_to_N
        for i in _1_to_N
            # i == j && continue
            if i == j
            else
                X[i,j] = q[3i-2]-q[3j-2]
                Y[i,j] = q[3i-1]-q[3j-1]
                Z[i,j] = q[3i]-q[3j]

                U[i,j] = dq[3i-2]-dq[3j-2]
                V[i,j] = dq[3i-1]-dq[3j-1]
                W[i,j] = dq[3i  ]-dq[3j  ]

                _4U_m_3X[i,j] = (4dq[3j-2])-(3dq[3i-2])
                _4V_m_3Y[i,j] = (4dq[3j-1])-(3dq[3i-1])
                _4W_m_3Z[i,j] = (4dq[3j  ])-(3dq[3i  ])

                pn2x = X[i,j]*_4U_m_3X[i,j]
                pn2y = Y[i,j]*_4V_m_3Y[i,j]
                pn2z = Z[i,j]*_4W_m_3Z[i,j]

                pn2[i,j] = ( pn2x+pn2y ) + pn2z

                UU[i,j] = dq[3i-2]*dq[3j-2]
                VV[i,j] = dq[3i-1]*dq[3j-1]
                WW[i,j] = dq[3i  ]*dq[3j  ]

                vi_dot_vj[i,j] = ( UU[i,j]+VV[i,j] ) + WW[i,j]

                r_p2[i,j] = ( (X[i,j]^2)+(Y[i,j]^2) ) + (Z[i,j]^2)

                r_p1d2[i,j] = sqrt(r_p2[i,j])
                r_p3d2[i,j] = r_p2[i,j]^1.5
                r_p7d2[i,j] = r_p2[i,j]^3.5

                newtonianCoeff[i,j] =  μ[i]/r_p3d2[i,j]

                temp_001 = newtonX[j] + (X[i,j]*newtonianCoeff[i,j])
                newtonX[j] = temp_001
                temp_002 = newtonY[j] + (Y[i,j]*newtonianCoeff[i,j])
                newtonY[j] = temp_002
                temp_003 = newtonZ[j] + (Z[i,j]*newtonianCoeff[i,j])
                newtonZ[j] = temp_003

                newtonian1b_Potential[i,j] = μ[i]/r_p1d2[i, j]

                temp_004 = newtonianNb_Potential[j] + newtonian1b_Potential[i, j]
                newtonianNb_Potential[j] = temp_004
            end #if i != j
        end #for, i
        v2[j] = ( (dq[3j-2]^2)+(dq[3j-1]^2) ) + (dq[3j]^2)
    end #for, j

    for j in _1_to_N
        postNewtonX[j] = newtonX[j]
        postNewtonY[j] = newtonY[j]
        postNewtonZ[j] = newtonZ[j]
    end

    for k in Base.OneTo(succ_approx_iter)
        for j in _1_to_N
            pntempX[j] = zero_q_1
            pntempY[j] = zero_q_1
            pntempZ[j] = zero_q_1
        end
        for j in _1_to_N
            for i in _1_to_N
                # i == j && continue
                if i == j
                else
                    #post-Newtonian corrections to gravitational acceleration
                    #Moyer, 1971, page 7 eq. 35
                    temp_005a = newtonianNb_Potential[i]+(4newtonianNb_Potential[j])
                    temp_005b = (2v2[i])-(4vi_dot_vj[i,j])
                    temp_005c = v2[j]+temp_005b
                    temp_005 = temp_005c-temp_005a
                    temp_006a = X[i,j]*dq[3i-2]
                    temp_006b = Y[i,j]*dq[3i-1]
                    temp_006c = Z[i,j]*dq[3i]
                    temp_006d = ( temp_006a+temp_006b ) + temp_006c
                    # the expression below inside the (...)^2 should have a minus sign in front of the numerator,
                    # but upon squaring it is eliminated, so at the end of the day, it is irrelevant ;)
                    temp_006e = (temp_006d^2)/r_p2[i,j]
                    temp_006 = temp_005-(1.5temp_006e)
                    temp_007a = X[i,j]*postNewtonX[i]
                    temp_007b = Y[i,j]*postNewtonY[i]
                    temp_007c = Z[i,j]*postNewtonZ[i]
                    temp_007d = ( temp_007a+temp_007b ) + temp_007c
                    temp_007 = temp_006 + (0.5temp_007d)
                    temp_008 = c_p2+temp_007
                    pn1[i,j] = newtonianCoeff[i,j]*temp_008

                    temp_009 = pntempX[j]+(X[i,j]*pn1[i,j])
                    temp_010 = pntempY[j]+(Y[i,j]*pn1[i,j])
                    temp_011 = pntempZ[j]+(Z[i,j]*pn1[i,j])

                    pn3[i,j] = 3.5*newtonian1b_Potential[i,j]

                    temp_013a = pn2[i,j]*(U[i,j]*newtonianCoeff[i,j])
                    temp_013b = postNewtonX[i]*pn3[i,j]
                    temp_013 = temp_009 + (temp_013a+temp_013b)
                    pntempX[j] = temp_013

                    temp_014a = pn2[i,j]*(V[i,j]*newtonianCoeff[i,j])
                    temp_014b = postNewtonY[i]*pn3[i,j]
                    temp_014 = temp_010 + (temp_014a+temp_014b)
                    pntempY[j] = temp_014

                    temp_015a = pn2[i,j]*(W[i,j]*newtonianCoeff[i,j])
                    temp_015b = postNewtonZ[i]*pn3[i,j]
                    temp_015 = temp_011 + (temp_015a+temp_015b)
                    pntempZ[j] = temp_015
                end
            end #for i
        end #for j
        for j in _1_to_N
            postNewtonX[j] = pntempX[j]*c_m2
            postNewtonY[j] = pntempY[j]*c_m2
            postNewtonZ[j] = pntempZ[j]*c_m2
        end
    end #for k in Base.OneTo(succ_approx_iter) # (post-Newtonian iterations)

    #J2 accelerations, all flattened bodies
    for j in j2_body_index
        for i in _1_to_N
            if i == j
            else
                # # rotate from inertial frame to extended-body frame
                t11[i,j] = X[i,j]*W_[1,1,j]
                t21[i,j] = X[i,j]*W_[2,1,j]
                t31[i,j] = X[i,j]*W_[3,1,j]
                t12[i,j] = Y[i,j]*W_[1,2,j]
                t22[i,j] = Y[i,j]*W_[2,2,j]
                t32[i,j] = Y[i,j]*W_[3,2,j]
                t13[i,j] = Z[i,j]*W_[1,3,j]
                t23[i,j] = Z[i,j]*W_[2,3,j]
                t33[i,j] = Z[i,j]*W_[3,3,j]
                new_x[i,j] = (t11[i,j]+t12[i,j])+t13[i,j]
                new_y[i,j] = (t21[i,j]+t22[i,j])+t23[i,j]
                new_z[i,j] = (t31[i,j]+t32[i,j])+t33[i,j]

                # # compute cartesian components of extended-body acceleration in body frame
                new_x2[i,j] = new_x[i,j]^2
                new_y2[i,j] = new_y[i,j]^2

                ρ_p2_ij[i,j] = new_x2[i,j]+new_y2[i,j]
                z_p2_ij[i,j] = new_z[i,j]^2

                temp_ρ[i,j] = ρ_p2_ij[i,j] - (4z_p2_ij[i,j])
                temp_z[i,j] = (3ρ_p2_ij[i,j]) - (2z_p2_ij[i,j])

                dum01[i,j] = Λ2[j]/r_p7d2[i,j]

                dum[i,j] = 1.5*dum01[i,j]
                dum_ρ[i,j] = dum[i,j]*temp_ρ[i,j]
                dum_z[i,j] = dum[i,j]*temp_z[i,j]

                F_J2_bf_x[i,j] = dum_ρ[i,j]*new_x[i,j]
                F_J2_bf_y[i,j] = dum_ρ[i,j]*new_y[i,j]
                F_J2_bf_z[i,j] = dum_z[i,j]*new_z[i,j]

                # # rotate components of force from body frame to inertial frame
                s11[i,j] = F_J2_bf_x[i,j]*M_[1,1,j]
                s21[i,j] = F_J2_bf_x[i,j]*M_[2,1,j]
                s31[i,j] = F_J2_bf_x[i,j]*M_[3,1,j]
                s12[i,j] = F_J2_bf_y[i,j]*M_[1,2,j]
                s22[i,j] = F_J2_bf_y[i,j]*M_[2,2,j]
                s32[i,j] = F_J2_bf_y[i,j]*M_[3,2,j]
                s13[i,j] = F_J2_bf_z[i,j]*M_[1,3,j]
                s23[i,j] = F_J2_bf_z[i,j]*M_[2,3,j]
                s33[i,j] = F_J2_bf_z[i,j]*M_[3,3,j]
                F_J2_x[i,j] = (s11[i,j]+s12[i,j])+s13[i,j]
                F_J2_y[i,j] = (s21[i,j]+s22[i,j])+s23[i,j]
                F_J2_z[i,j] = (s31[i,j]+s32[i,j])+s33[i,j]

                # # add result to total acceleration upon j-th body figure due to i-th point mass
                # @show "acc",j,"+μ",i,"Λ2",j
                timp_004[i,j] = accX[j] + (μ[i]*F_J2_x[i,j])
                accX[j] = timp_004[i,j]
                timp_005[i,j] = accY[j] + (μ[i]*F_J2_y[i,j])
                accY[j] = timp_005[i,j]
                timp_006[i,j] = accZ[j] + (μ[i]*F_J2_z[i,j])
                accZ[j] = timp_006[i,j]

                # # reaction force on i-th body
                # @show "acc",i,"-μ",j,"Λ2",j
                tamp_004[i,j] = accX[i] - (μ[j]*F_J2_x[i,j])
                accX[i] = tamp_004[i,j]
                tamp_005[i,j] = accY[i] - (μ[j]*F_J2_y[i,j])
                accY[i] = tamp_005[i,j]
                tamp_006[i,j] = accZ[i] - (μ[j]*F_J2_z[i,j])
                accZ[i] = tamp_006[i,j]
            end #if i == j
        end #for i in _1_to_N
    end #for j in j2_body_index

    #fill the equations of motion for everyone except test particle (Newtonian, post-Newtonian and extended body accelerations)
    for i in Base.OneTo(N-1)
        dq[3(N+i)-2] = postNewtonX[i]+accX[i]
        dq[3(N+i)-1] = postNewtonY[i]+accY[i]
        dq[3(N+i)  ] = postNewtonZ[i]+accZ[i]
    end

    #computation of non-gravitational accelerations:
    hx = (Y[N,1]*(dq[3N  ]-dq[3]))-(Z[N,1]*(dq[3N-1]-dq[2]))
    hy = (Z[N,1]*(dq[3N-2]-dq[1]))-(X[N,1]*(dq[3N  ]-dq[3]))
    hz = (X[N,1]*(dq[3N-1]-dq[2]))-(Y[N,1]*(dq[3N-2]-dq[1]))
    r_hs = sqrt(r_p2[N,1])
    runitx = X[N,1]/r_hs
    runity = Y[N,2]/r_hs
    runitz = Z[N,3]/r_hs

    #cartesian components of transversal unit vector:
    tunitx0 = (hy*runitz)-(hz*runity)
    tunity0 = (hz*runitx)-(hx*runitz)
    tunitz0 = (hx*runity)-(hy*runitx)
    hmag = sqrt( ((tunitx0^2)+(tunity0^2))+(tunitz0^2) )
    tunitx = tunitx0/hmag
    tunity = tunity0/hmag
    tunitz = tunitz0/hmag

    # evaluate non-grav acceleration of NEA (Yarkovsky):
    g_r = r_hs^(-2.25)
    A2_t_g_r = q[6N+1]*g_r

    NGAx = A2_t_g_r*tunitx
    NGAy = A2_t_g_r*tunity
    NGAz = A2_t_g_r*tunitz

    dq[6N-2] = (postNewtonX[N]+accX[N])+NGAx
    dq[6N-1] = (postNewtonY[N]+accY[N])+NGAy
    dq[6N  ] = (postNewtonZ[N]+accZ[N])+NGAz

    dq[6N+1] = zero_q_1

    nothing
end

# Nearth-Earth asteroid dynamical model (d=2.25 and direct solar radiation pressure)
# Bodies considered in the model are: the Sun, the eight planets, the Moon and Ceres,
# as well as the asteroid of interest as a test particle with null mass. Dynamical
# effects considered are:
# - post-Newtonian point-mass accelerations between all bodies,
# - figure-effects (oblateness) of the Earth (J2, J3 and J4)
# - J2 effect of the Sun
# - Kinematic model for the precession and nutation of the Earth's pole (DE430/431 model)
# - also, a model for non-gravitational accelerations acting upon the asteroid
# is included (Yarkovsky effect) a_nongrav = A2*t_vec*(au/r)^d, where t_vec is the
# unit heliocentric transverse vector, au is 1 astronomical unit, r is the
# asteroid's heliocentric range, A2 is a coefficient (with units of au/day^2),
# and d = 2.25
@taylorize function RNp1BP_pN_A_J234E_J2S_ng_d225_srp!(t, q, dq)
    local S = eltype(q[1])
    local N = Int((length(q)-1)/6) # number of bodies, including NEA
    local _1_to_N = Base.OneTo(N) # iterator over all bodies

    local succ_approx_iter = 1 # number of iterations of post-Newtonian subroutine
    local su = 1 #Sun's index within `system`
    local ea = 4 #Earth's index within `system`
    local j2_body_index = [su, ea] # indices of bodies with J2 flattening (note: Earth also has J3 and J4)

    # parameters related to speed of light, c
    local c_p2 = 29979.063823897606 # c^2 = 29979.063823897606 au^2/d^2
    local c_m2 = 3.3356611996764786e-5 # c^-2 = 3.3356611996764786e-5 d^2/au^2

    local amrat = 3.07e-6 # m^2/kg, area-to-mass ratio (Bennu)

    local zero_q_1 = zero(q[1])

    X = Array{Taylor1{S}}(undef, N, N)
    Y = Array{Taylor1{S}}(undef, N, N)
    Z = Array{Taylor1{S}}(undef, N, N)

    r_p2 = Array{Taylor1{S}}(undef, N, N)
    r_p3d2 = Array{Taylor1{S}}(undef, N, N)
    r_p7d2 = Array{Taylor1{S}}(undef, N, N)

    #Newtonian accelerations
    newtonX = Array{Taylor1{S}}(undef, N)
    newtonY = Array{Taylor1{S}}(undef, N)
    newtonZ = Array{Taylor1{S}}(undef, N)

    newtonianCoeff = Array{Taylor1{S}}(undef, N, N)

    #post-Newtonian stuff
    U = Array{Taylor1{S}}(undef, N, N)
    V = Array{Taylor1{S}}(undef, N, N)
    W = Array{Taylor1{S}}(undef, N, N)

    _4U_m_3X = Array{Taylor1{S}}(undef, N, N)
    _4V_m_3Y = Array{Taylor1{S}}(undef, N, N)
    _4W_m_3Z = Array{Taylor1{S}}(undef, N, N)

    UU = Array{Taylor1{S}}(undef, N, N)
    VV = Array{Taylor1{S}}(undef, N, N)
    WW = Array{Taylor1{S}}(undef, N, N)

    r_p1d2 = Array{Taylor1{S}}(undef, N, N)

    postNewtonX = Array{Taylor1{S}}(undef, N)
    postNewtonY = Array{Taylor1{S}}(undef, N)
    postNewtonZ = Array{Taylor1{S}}(undef, N)

    newtonianNb_Potential = Array{Taylor1{S}}(undef, N)
    newtonian1b_Potential = Array{Taylor1{S}}(undef, N, N)
    newtonianCoeff = Array{Taylor1{S}}(undef, N, N)

    pntempX = Array{Taylor1{S}}(undef, N)
    pntempY = Array{Taylor1{S}}(undef, N)
    pntempZ = Array{Taylor1{S}}(undef, N)

    pn1 = Array{Taylor1{S}}(undef, N, N)
    v2 = Array{Taylor1{S}}(undef, N)
    vi_dot_vj = Array{Taylor1{S}}(undef, N, N)
    pn2 = Array{Taylor1{S}}(undef, N, N)
    pn3 = Array{Taylor1{S}}(undef, N, N)

    # J2 acceleration auxiliaries
    t11 = Array{Taylor1{S}}(undef, N, N)
    t12 = Array{Taylor1{S}}(undef, N, N)
    t13 = Array{Taylor1{S}}(undef, N, N)
    t21 = Array{Taylor1{S}}(undef, N, N)
    t22 = Array{Taylor1{S}}(undef, N, N)
    t23 = Array{Taylor1{S}}(undef, N, N)
    t31 = Array{Taylor1{S}}(undef, N, N)
    t32 = Array{Taylor1{S}}(undef, N, N)
    t33 = Array{Taylor1{S}}(undef, N, N)
    new_x = Array{Taylor1{S}}(undef, N, N)
    new_y = Array{Taylor1{S}}(undef, N, N)
    new_z = Array{Taylor1{S}}(undef, N, N)
    new_x2 = Array{Taylor1{S}}(undef, N, N)
    new_y2 = Array{Taylor1{S}}(undef, N, N)
    ρ_p2_ij = Array{Taylor1{S}}(undef, N, N)
    z_p2_ij = Array{Taylor1{S}}(undef, N, N)
    temp_ρ = Array{Taylor1{S}}(undef, N, N)
    temp_z = Array{Taylor1{S}}(undef, N, N)
    dum00 = Array{Taylor1{S}}(undef, N, N)
    dum01 = Array{Taylor1{S}}(undef, N, N)
    dum = Array{Taylor1{S}}(undef, N, N)
    dum_ρ = Array{Taylor1{S}}(undef, N, N)
    dum_z = Array{Taylor1{S}}(undef, N, N)
    F_J2_bf_x = Array{Taylor1{S}}(undef, N, N)
    F_J2_bf_y = Array{Taylor1{S}}(undef, N, N)
    F_J2_bf_z = Array{Taylor1{S}}(undef, N, N)
    s11 = Array{Taylor1{S}}(undef, N, N)
    s12 = Array{Taylor1{S}}(undef, N, N)
    s13 = Array{Taylor1{S}}(undef, N, N)
    s21 = Array{Taylor1{S}}(undef, N, N)
    s22 = Array{Taylor1{S}}(undef, N, N)
    s23 = Array{Taylor1{S}}(undef, N, N)
    s31 = Array{Taylor1{S}}(undef, N, N)
    s32 = Array{Taylor1{S}}(undef, N, N)
    s33 = Array{Taylor1{S}}(undef, N, N)
    F_J2_x = Array{Taylor1{S}}(undef, N, N)
    F_J2_y = Array{Taylor1{S}}(undef, N, N)
    F_J2_z = Array{Taylor1{S}}(undef, N, N)
    timp_004 = Array{Taylor1{S}}(undef, N, N)
    timp_005 = Array{Taylor1{S}}(undef, N, N)
    timp_006 = Array{Taylor1{S}}(undef, N, N)
    tamp_004 = Array{Taylor1{S}}(undef, N, N)
    tamp_005 = Array{Taylor1{S}}(undef, N, N)
    tamp_006 = Array{Taylor1{S}}(undef, N, N)

    # extended-body accelerations
    accX = Array{Taylor1{S}}(undef, N)
    accY = Array{Taylor1{S}}(undef, N)
    accZ = Array{Taylor1{S}}(undef, N)

    # rotations to and from Earth and Sun poles
    local αs = deg2rad(268.13+0t)
    local δs = deg2rad(63.87+0t)
    local M_ = Array{Taylor1{S}}(undef, 3, 3, N)
    local W_ = Array{Taylor1{S}}(undef, 3, 3, N)
    local M_[:,:,ea] = earth_pole_rotation(t-J2000) # J2000.0 = 2.451545e6
    local W_[:,:,ea] = inv(M_[:, :, ea])
    local M_[:,:,su] = pole_rotation( αs, δs )
    local W_[:,:,su] = inv(M_[:, :, su])

    for j in _1_to_N
        newtonX[j] = zero_q_1
        newtonY[j] = zero_q_1
        newtonZ[j] = zero_q_1

        newtonianNb_Potential[j] = zero_q_1

        accX[j] = zero_q_1
        accY[j] = zero_q_1
        accZ[j] = zero_q_1

        dq[3j-2] = q[3(N+j)-2]
        dq[3j-1] = q[3(N+j)-1]
        dq[3j  ] = q[3(N+j)  ]
    end

    for j in j2_body_index
        for i in _1_to_N
            if i == j
            else
                t11[i,j] = zero_q_1
                t12[i,j] = zero_q_1
                t13[i,j] = zero_q_1
                t21[i,j] = zero_q_1
                t22[i,j] = zero_q_1
                t23[i,j] = zero_q_1
                t31[i,j] = zero_q_1
                t32[i,j] = zero_q_1
                t33[i,j] = zero_q_1
                new_x[i,j] = zero_q_1
                new_x2[i,j] = zero_q_1
                new_y2[i,j] = zero_q_1
                ρ_p2_ij[i,j] = zero_q_1
                z_p2_ij[i,j] = zero_q_1
                temp_ρ[i,j] = zero_q_1
                temp_z[i,j] = zero_q_1
                dum00[i,j] = zero_q_1
                dum01[i,j] = zero_q_1
                dum[i,j] = zero_q_1
                dum_ρ[i,j] = zero_q_1
                dum_z[i,j] = zero_q_1
                F_J2_bf_x[i,j] = zero_q_1
                F_J2_bf_y[i,j] = zero_q_1
                F_J2_bf_z[i,j] = zero_q_1
                s11[i,j] = zero_q_1
                s12[i,j] = zero_q_1
                s13[i,j] = zero_q_1
                s21[i,j] = zero_q_1
                s22[i,j] = zero_q_1
                s23[i,j] = zero_q_1
                s31[i,j] = zero_q_1
                s32[i,j] = zero_q_1
                s33[i,j] = zero_q_1
                F_J2_x[i,j] = zero_q_1
                F_J2_y[i,j] = zero_q_1
                F_J2_z[i,j] = zero_q_1
                timp_004[i,j] = zero_q_1
                timp_005[i,j] = zero_q_1
                timp_006[i,j] = zero_q_1
                tamp_004[i,j] = zero_q_1
                tamp_005[i,j] = zero_q_1
                tamp_006[i,j] = zero_q_1
            end #if i == j
        end #for i in _1_to_N
    end #for j in j2_body_index

    #compute point-mass Newtonian accelerations, all bodies
    for j in _1_to_N
        for i in _1_to_N
            # i == j && continue
            if i == j
            else
                X[i,j] = q[3i-2]-q[3j-2]
                Y[i,j] = q[3i-1]-q[3j-1]
                Z[i,j] = q[3i]-q[3j]

                U[i,j] = dq[3i-2]-dq[3j-2]
                V[i,j] = dq[3i-1]-dq[3j-1]
                W[i,j] = dq[3i  ]-dq[3j  ]

                _4U_m_3X[i,j] = (4dq[3j-2])-(3dq[3i-2])
                _4V_m_3Y[i,j] = (4dq[3j-1])-(3dq[3i-1])
                _4W_m_3Z[i,j] = (4dq[3j  ])-(3dq[3i  ])

                pn2x = X[i,j]*_4U_m_3X[i,j]
                pn2y = Y[i,j]*_4V_m_3Y[i,j]
                pn2z = Z[i,j]*_4W_m_3Z[i,j]

                pn2[i,j] = ( pn2x+pn2y ) + pn2z

                UU[i,j] = dq[3i-2]*dq[3j-2]
                VV[i,j] = dq[3i-1]*dq[3j-1]
                WW[i,j] = dq[3i  ]*dq[3j  ]

                vi_dot_vj[i,j] = ( UU[i,j]+VV[i,j] ) + WW[i,j]

                r_p2[i,j] = ( (X[i,j]^2)+(Y[i,j]^2) ) + (Z[i,j]^2)

                r_p1d2[i,j] = sqrt(r_p2[i,j])
                r_p3d2[i,j] = r_p2[i,j]^1.5
                r_p7d2[i,j] = r_p2[i,j]^3.5

                newtonianCoeff[i,j] =  μ[i]/r_p3d2[i,j]

                temp_001 = newtonX[j] + (X[i,j]*newtonianCoeff[i,j])
                newtonX[j] = temp_001
                temp_002 = newtonY[j] + (Y[i,j]*newtonianCoeff[i,j])
                newtonY[j] = temp_002
                temp_003 = newtonZ[j] + (Z[i,j]*newtonianCoeff[i,j])
                newtonZ[j] = temp_003

                newtonian1b_Potential[i,j] = μ[i]/r_p1d2[i, j]

                temp_004 = newtonianNb_Potential[j] + newtonian1b_Potential[i, j]
                newtonianNb_Potential[j] = temp_004
            end #if i != j
        end #for, i
        v2[j] = ( (dq[3j-2]^2)+(dq[3j-1]^2) ) + (dq[3j]^2)
    end #for, j

    for j in _1_to_N
        postNewtonX[j] = newtonX[j]
        postNewtonY[j] = newtonY[j]
        postNewtonZ[j] = newtonZ[j]
    end

    for k in Base.OneTo(succ_approx_iter)
        for j in _1_to_N
            pntempX[j] = zero_q_1
            pntempY[j] = zero_q_1
            pntempZ[j] = zero_q_1
        end
        for j in _1_to_N
            for i in _1_to_N
                # i == j && continue
                if i == j
                else
                    #post-Newtonian corrections to gravitational acceleration
                    #Moyer, 1971, page 7 eq. 35
                    temp_005a = newtonianNb_Potential[i]+(4newtonianNb_Potential[j])
                    temp_005b = (2v2[i])-(4vi_dot_vj[i,j])
                    temp_005c = v2[j]+temp_005b
                    temp_005 = temp_005c-temp_005a
                    temp_006a = X[i,j]*dq[3i-2]
                    temp_006b = Y[i,j]*dq[3i-1]
                    temp_006c = Z[i,j]*dq[3i]
                    temp_006d = ( temp_006a+temp_006b ) + temp_006c
                    # the expression below inside the (...)^2 should have a minus sign in front of the numerator,
                    # but upon squaring it is eliminated, so at the end of the day, it is irrelevant ;)
                    temp_006e = (temp_006d^2)/r_p2[i,j]
                    temp_006 = temp_005-(1.5temp_006e)
                    temp_007a = X[i,j]*postNewtonX[i]
                    temp_007b = Y[i,j]*postNewtonY[i]
                    temp_007c = Z[i,j]*postNewtonZ[i]
                    temp_007d = ( temp_007a+temp_007b ) + temp_007c
                    temp_007 = temp_006 + (0.5temp_007d)
                    temp_008 = c_p2+temp_007
                    pn1[i,j] = newtonianCoeff[i,j]*temp_008

                    temp_009 = pntempX[j]+(X[i,j]*pn1[i,j])
                    temp_010 = pntempY[j]+(Y[i,j]*pn1[i,j])
                    temp_011 = pntempZ[j]+(Z[i,j]*pn1[i,j])

                    pn3[i,j] = 3.5*newtonian1b_Potential[i,j]

                    temp_013a = pn2[i,j]*(U[i,j]*newtonianCoeff[i,j])
                    temp_013b = postNewtonX[i]*pn3[i,j]
                    temp_013 = temp_009 + (temp_013a+temp_013b)
                    pntempX[j] = temp_013

                    temp_014a = pn2[i,j]*(V[i,j]*newtonianCoeff[i,j])
                    temp_014b = postNewtonY[i]*pn3[i,j]
                    temp_014 = temp_010 + (temp_014a+temp_014b)
                    pntempY[j] = temp_014

                    temp_015a = pn2[i,j]*(W[i,j]*newtonianCoeff[i,j])
                    temp_015b = postNewtonZ[i]*pn3[i,j]
                    temp_015 = temp_011 + (temp_015a+temp_015b)
                    pntempZ[j] = temp_015
                end
            end #for i
        end #for j
        for j in _1_to_N
            postNewtonX[j] = pntempX[j]*c_m2
            postNewtonY[j] = pntempY[j]*c_m2
            postNewtonZ[j] = pntempZ[j]*c_m2
        end
    end #for k in Base.OneTo(succ_approx_iter) # (post-Newtonian iterations)

    #J2 accelerations, all flattened bodies
    for j in j2_body_index
        for i in _1_to_N
            if i == j
            else
                # # rotate from inertial frame to extended-body frame
                t11[i,j] = X[i,j]*W_[1,1,j]
                t21[i,j] = X[i,j]*W_[2,1,j]
                t31[i,j] = X[i,j]*W_[3,1,j]
                t12[i,j] = Y[i,j]*W_[1,2,j]
                t22[i,j] = Y[i,j]*W_[2,2,j]
                t32[i,j] = Y[i,j]*W_[3,2,j]
                t13[i,j] = Z[i,j]*W_[1,3,j]
                t23[i,j] = Z[i,j]*W_[2,3,j]
                t33[i,j] = Z[i,j]*W_[3,3,j]
                new_x[i,j] = (t11[i,j]+t12[i,j])+t13[i,j]
                new_y[i,j] = (t21[i,j]+t22[i,j])+t23[i,j]
                new_z[i,j] = (t31[i,j]+t32[i,j])+t33[i,j]

                # # compute cartesian components of extended-body acceleration in body frame
                new_x2[i,j] = new_x[i,j]^2
                new_y2[i,j] = new_y[i,j]^2

                ρ_p2_ij[i,j] = new_x2[i,j]+new_y2[i,j]
                z_p2_ij[i,j] = new_z[i,j]^2

                temp_ρ[i,j] = ρ_p2_ij[i,j] - (4z_p2_ij[i,j])
                temp_z[i,j] = (3ρ_p2_ij[i,j]) - (2z_p2_ij[i,j])

                dum01[i,j] = Λ2[j]/r_p7d2[i,j]

                dum[i,j] = 1.5*dum01[i,j]
                dum_ρ[i,j] = dum[i,j]*temp_ρ[i,j]
                dum_z[i,j] = dum[i,j]*temp_z[i,j]

                F_J2_bf_x[i,j] = dum_ρ[i,j]*new_x[i,j]
                F_J2_bf_y[i,j] = dum_ρ[i,j]*new_y[i,j]
                F_J2_bf_z[i,j] = dum_z[i,j]*new_z[i,j]

                # # rotate components of force from body frame to inertial frame
                s11[i,j] = F_J2_bf_x[i,j]*M_[1,1,j]
                s21[i,j] = F_J2_bf_x[i,j]*M_[2,1,j]
                s31[i,j] = F_J2_bf_x[i,j]*M_[3,1,j]
                s12[i,j] = F_J2_bf_y[i,j]*M_[1,2,j]
                s22[i,j] = F_J2_bf_y[i,j]*M_[2,2,j]
                s32[i,j] = F_J2_bf_y[i,j]*M_[3,2,j]
                s13[i,j] = F_J2_bf_z[i,j]*M_[1,3,j]
                s23[i,j] = F_J2_bf_z[i,j]*M_[2,3,j]
                s33[i,j] = F_J2_bf_z[i,j]*M_[3,3,j]
                F_J2_x[i,j] = (s11[i,j]+s12[i,j])+s13[i,j]
                F_J2_y[i,j] = (s21[i,j]+s22[i,j])+s23[i,j]
                F_J2_z[i,j] = (s31[i,j]+s32[i,j])+s33[i,j]

                # # add result to total acceleration upon j-th body figure due to i-th point mass
                # @show "acc",j,"+μ",i,"Λ2",j
                timp_004[i,j] = accX[j] + (μ[i]*F_J2_x[i,j])
                accX[j] = timp_004[i,j]
                timp_005[i,j] = accY[j] + (μ[i]*F_J2_y[i,j])
                accY[j] = timp_005[i,j]
                timp_006[i,j] = accZ[j] + (μ[i]*F_J2_z[i,j])
                accZ[j] = timp_006[i,j]

                # # reaction force on i-th body
                # @show "acc",i,"-μ",j,"Λ2",j
                tamp_004[i,j] = accX[i] - (μ[j]*F_J2_x[i,j])
                accX[i] = tamp_004[i,j]
                tamp_005[i,j] = accY[i] - (μ[j]*F_J2_y[i,j])
                accY[i] = tamp_005[i,j]
                tamp_006[i,j] = accZ[i] - (μ[j]*F_J2_z[i,j])
                accZ[i] = tamp_006[i,j]
            end #if i == j
        end #for i in _1_to_N
    end #for j in j2_body_index

    #fill the equations of motion for everyone except test particle (Newtonian, post-Newtonian and extended body accelerations)
    for i in Base.OneTo(N-1)
        dq[3(N+i)-2] = postNewtonX[i]+accX[i]
        dq[3(N+i)-1] = postNewtonY[i]+accY[i]
        dq[3(N+i)  ] = postNewtonZ[i]+accZ[i]
    end

    #computation of non-gravitational accelerations:
    hx = (Y[N,1]*(dq[3N  ]-dq[3]))-(Z[N,1]*(dq[3N-1]-dq[2]))
    hy = (Z[N,1]*(dq[3N-2]-dq[1]))-(X[N,1]*(dq[3N  ]-dq[3]))
    hz = (X[N,1]*(dq[3N-1]-dq[2]))-(Y[N,1]*(dq[3N-2]-dq[1]))
    r_hs = sqrt(r_p2[N,1])
    runitx = X[N,1]/r_hs
    runity = Y[N,2]/r_hs
    runitz = Z[N,3]/r_hs

    #cartesian components of transversal unit vector:
    tunitx0 = (hy*runitz)-(hz*runity)
    tunity0 = (hz*runitx)-(hx*runitz)
    tunitz0 = (hx*runity)-(hy*runitx)
    hmag = sqrt( ((tunitx0^2)+(tunity0^2))+(tunitz0^2) )
    tunitx = tunitx0/hmag
    tunity = tunity0/hmag
    tunitz = tunitz0/hmag

    # evaluate non-grav acceleration of asteroid (Yarkovsky):
    g_r = r_hs^(-2.25)
    A2_t_g_r = q[6N+1]*g_r

    NGAx = A2_t_g_r*tunitx
    NGAy = A2_t_g_r*tunity
    NGAz = A2_t_g_r*tunitz

    # evaluate solar radation pressure
    SRP_coeff = ((S0_sun*amrat)*m2_s3_to_au2_day3)/c_au_per_day
    SRP_radial = (R_sun/r_hs)^2
    SRP_fun = SRP_coeff*SRP_radial

    SRPx = SRP_fun*runitx
    SRPy = SRP_fun*runity
    SRPz = SRP_fun*runitz

    #add non-gravitational accelerations to asteroid acceleration
    dq[6N-2] = (postNewtonX[N]+accX[N])+(NGAx+SRPx)
    dq[6N-1] = (postNewtonY[N]+accY[N])+(NGAy+SRPy)
    dq[6N  ] = (postNewtonZ[N]+accZ[N])+(NGAz+SRPz)

    dq[6N+1] = zero_q_1

    nothing
end
