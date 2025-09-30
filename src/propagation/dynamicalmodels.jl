"""
    nongravs!

Asteroid dynamical model specially suited for long-term integrations and impact monitoring.

The model considers the asteroid of interest as a test particle with null mass. Perturbing
bodies included in the model are: the Sun, the eight planets, the Moon, Pluto and the 16
most massive asteroids. Planetary ephemerides are provided by `PlanetaryEphemeris.jl`,
which is based on the JPL DE430 model.

Dynamical effects considered are:

- Post-Newtonian point-mass accelerations between all bodies. See equation (35) in page
    7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract.

- Figure-effects (oblateness) of the Earth (``J_2`` and ``J_3``), ``J_2`` effect of the
    Sun and ``J_2`` and ``J_3`` effect of the Moon. See equation (28) in page 13 of
    https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract and equations (173)
    and (174) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract.

- Kinematic model for the precession and nutation of the Earth's orientation (IAU 1976/1980
    Earth orientation model). See [`PlanetaryEphemeris.c2t_jpl_de430`](@ref).

- Kinematic model for the Moon's orientation (Seidelmann et al., 2006). See equations
    (14)-(15) in page 9 and equations (34)-(35) in page 16 of
    https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract.

- Non-gravitational accelerations model (Marsden et al., 1973). See equations (1)-(5) in
    pages (211)-(212) of https://articles.adsabs.harvard.edu/pdf/1973AJ.....78..211M.

To improve performance, some internal loops are multi-threaded via `@threads`.

For other dynamical models, see [`gravityonly!`](@ref) and [`newtonian!`](@ref).
"""
function nongravs!(dq, q, params, t)
    # Julian date (TDB) of start time
    local jd0 = params.jd0
    # Days since J2000.0 = 2.451545e6
    local dsj2k = t + (jd0 - JD_J2000)
    # Solar system ephemeris at dsj2k
    local ss16asteph_t = params.ssepht
    evaleph!(ss16asteph_t, params.sseph, dsj2k)
    # Accelerations at dsj2k
    local acceph_t = params.accepht
    evaleph!(acceph_t, params.acceph, dsj2k)
    # Newtonian potentials at dsj2k
    local newtonianNb_Potential_t = params.potepht
    evaleph!(newtonianNb_Potential_t, params.poteph, dsj2k)
    # Type of position / velocity components
    local S = eltype(q)
    # Interaction matrix with flattened bodies
    local UJ_interaction = params.UJ_interaction
    # Number of bodies, including NEA
    local N = params.N
    # Number of bodies, except the asteroid
    local Nm1 = N-1
    # Vector of mass parameters GM's
    local μ = params.μ
    # Marsden et al. (1973) radial function constants
    local marsden_α = params.marsden_radial[1]
    local marsden_r₀ = params.marsden_radial[2]
    local marsden_m = -params.marsden_radial[3]
    local marsden_n = params.marsden_radial[4]
    local marsden_k = -params.marsden_radial[5]

    # zero(q[1])
    local zero_q_1 = auxzero(q[1])

    #=
    Point-mass accelerations
    See equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    =#

    # Position of the i-th body - position of the asteroid
    X = Array{S}(undef, N)         # X-axis component
    Y = Array{S}(undef, N)         # Y-axis component
    Z = Array{S}(undef, N)         # Z-axis component

    # Distance between the i-th body and the asteroid
    r_p2 = Array{S}(undef, N)      # r_{i,asteroid}^2
    r_p1d2 = Array{S}(undef, N)    # sqrt(r_p2) <-> r_{i, asteroid}
    r_p3d2 = Array{S}(undef, N)    # r_p2^1.5 <-> r_{i, asteroid}^3
    r_p7d2 = Array{S}(undef, N)    # r_p2^3.5 <-> r_{i, asteroid}^7

    # Newtonian coefficient, i.e., mass parameter / distance^3 -> \mu_i / r_{i, asteroid}^3
    newtonianCoeff = Array{S}(undef, N)

    # Velocity of the i-th body
    ui = Array{S}(undef, N-1)      # X-axis component
    vi = Array{S}(undef, N-1)      # Y-axis component
    wi = Array{S}(undef, N-1)      # Z-axis component

    # Post-Newtonian stuff

    # Velocity of the i-th body - velocity of the asteroid
    U = Array{S}(undef, N)         # X-axis component
    V = Array{S}(undef, N)         # Y-axis component
    W = Array{S}(undef, N)         # Z-axis component

    # 4 * Velocity of the asteroid
    _4dq = Array{S}(undef, 3)  # (X,Y,Z) components

    # 4 * Velocity of the asteroid - 3 * velocity of the i-th body
    _4U_m_3X = Array{S}(undef, N)  # X-axis component
    _4V_m_3Y = Array{S}(undef, N)  # Y-axis component
    _4W_m_3Z = Array{S}(undef, N)  # Z-axis component

    # v_{i,j}v_{asteroid,j} j = x, y, z
    UU = Array{S}(undef, N)        # v_{i,x}v_{asteroid,x}
    VV = Array{S}(undef, N)        # v_{i,y}v_{asteroid,y}
    WW = Array{S}(undef, N)        # v_{i,z}v_{asteroid,z}

    # Newtonian potential of 1 body \mu_i / r_{i, asteroid}
    newtonian1b_Potential = Array{S}(undef, N)
    # Newtonian potential of N bodies
    # \sum_{i\neq l} \frac{\mu_i}{r_{il}} or
    # \sum_{j\neq k} \frac{\mu_j}{r_{jk}}
    newtonianNb_Potential = Array{S}(undef, N)

    # Newtonian coefficient * difference between two positions, i.e.,
    # \mu_i * (\mathbf{r_i} - \mathbf{r_asteroid}) / r_{ij}^3
    newton_acc_X = Array{S}(undef, N)   # X-axis component
    newton_acc_Y = Array{S}(undef, N)   # Y-axis component
    newton_acc_Z = Array{S}(undef, N)   # Z-axis component

    # Combinations of velocities
    v2 = Array{S}(undef, N)             # Velocity magnitude squared ||\mathbf{v}_i||^2
    vi_dot_vj = Array{S}(undef, N)      # <Velocity of the i-th body, velocity of the asteroid>

    # Second term without (\mathbf{v}_i - \mathbf{v}_j)
    pn2 = Array{S}(undef, N)
    # Full second term
    U_t_pn2 = Array{S}(undef, N)        # X-axis component
    V_t_pn2 = Array{S}(undef, N)        # Y-axis component
    W_t_pn2 = Array{S}(undef, N)        # Z-axis component

    # Third term without newtonian accelerations \mathbf{a}_i
    pn3 = Array{S}(undef, N)
    # Full third term of equation (35)
    pNX_t_pn3 = Array{S}(undef, N)      # X-axis component
    pNY_t_pn3 = Array{S}(undef, N)      # Y-axis component
    pNZ_t_pn3 = Array{S}(undef, N)      # Z-axis component

    # First term
    _4ϕj = Array{S}(undef, N)           # 4*\sum term inside {}
    ϕi_plus_4ϕj = Array{S}(undef, N)    # 4*\sum + \sum terms inside {}
    sj2_plus_2si2_minus_4vivj = Array{S}(undef, N)  # \dot{s}_j^2 + 2\dot{s}_i^2 - 4 <, > terms inside {}
    ϕs_and_vs = Array{S}(undef, N)      # -4\sum - \sum + \dot{s}_j^2 + 2\dot{s}_i^2  - 4<, > terms inside {}
    pn1t1_7 = Array{S}(undef, N)        # Everything inside the {} in the first term except for the term with accelerations (last)
    # Last term inside the {}
    pNX_t_X = Array{S}(undef, N)        # X-axis component
    pNY_t_Y = Array{S}(undef, N)        # Y-axis component
    pNZ_t_Z = Array{S}(undef, N)        # Z-axis component
    # Everything inside the {} in the first term
    pn1 = Array{S}(undef, N)
    # Full first term
    X_t_pn1 = Array{S}(undef, N)     # X-axis component
    Y_t_pn1 = Array{S}(undef, N)     # Y-axis component
    Z_t_pn1 = Array{S}(undef, N)     # Z-axis component

    # Temporary post-Newtonian accelerations
    pntempX = zero_q_1               # X-axis component
    pntempY = zero_q_1               # Y-axis component
    pntempZ = zero_q_1               # Z-axis component

    #=
    Extended body accelerations
    See equation (28) in page 13 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    and equations (173) and (174) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    =#

    # J_2 acceleration auxiliaries

    # Auxiliaries to compute body-fixed frame coordinates
    t31 = Array{S}(undef, N)
    t32 = Array{S}(undef, N)
    t33 = Array{S}(undef, N)
    # z-coordinate in body-fixed frame
    r_sin_ϕ = Array{S}(undef, N)

    # Trigonometric functions of latitude ϕ and longitude λ in the body-fixed coordinate system
    # See equations (165)-(168) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    sin_ϕ = Array{S}(undef, N)      # sin(latitude ϕ)
    ϕ = Array{S}(undef, N)          # Latitude ϕ
    cos_ϕ = Array{S}(undef, N)      # cos(latitude ϕ)
    sin2_ϕ = Array{S}(undef, N)     # sin(latitude ϕ)^2
    sin3_ϕ = Array{S}(undef, N)     # sin(latitude ϕ)^3
    sin4_ϕ = Array{S}(undef, N)     # sin(latitude ϕ)^4

    # Acceleration due to zonal harmonics in inertial frame
    F_J2_x = Array{S}(undef, N)
    F_J2_y = Array{S}(undef, N)
    F_J2_z = Array{S}(undef, N)
    # Auxiliaries to compute F_J2_i, i = x, y, z
    F_J2_x1 = Array{S}(undef, N)
    F_J2_y1 = Array{S}(undef, N)
    F_J2_z1 = Array{S}(undef, N)
    F_J2_x2 = Array{S}(undef, N)
    F_J2_y2 = Array{S}(undef, N)
    F_J2_z2 = Array{S}(undef, N)

    # Temporary arrays for the sum of full extended body accelerations
    temp_accX_i = Array{S}(undef, N)
    temp_accY_i = Array{S}(undef, N)
    temp_accZ_i = Array{S}(undef, N)

    # Legendre polynomials
    # See equations (176) and (177) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    P_2_sin_ϕ = Array{S}(undef, N)      # Second Legendre polynomial P_2(sin ϕ)
    ∂P_2_sin_ϕ = Array{S}(undef, N)     # dP_2(sin ϕ)/d(sin ϕ)
    P_3_sin_ϕ = Array{S}(undef, N)      # Third Legendre polynomial P_3(sin ϕ)
    ∂P_3_sin_ϕ = Array{S}(undef, N)     # dP_3(sin ϕ)/d(sin ϕ)
    # -cos ϕ P_n'
    m_c_ϕ_∂P_2 = Array{S}(undef, N)     # -cos ϕ P_2'
    m_c_ϕ_∂P_3 = Array{S}(undef, N)     # -cos ϕ P_3'

    # -J_n * R^n / r^m
    # J_n: n-th zonal harmonic coefficient
    # R: radius of the body
    # r: distance between the body and the asteroid
    Λ2j_div_r4 = Array{S}(undef, N)   # J_2 * R^2 / r^4
    Λ3j_div_r5 = Array{S}(undef, N)   # J_3 * R^3 / r^5

    # Accelerations due to zonal harmonics in body frame

    # Acceleration due to zonal harmonics J_n, n = 2, 3
    F_J_ξ = Array{S}(undef, N)         # ξ-axis component
    F_J_η = Array{S}(undef, N)         # η-axis component
    F_J_ζ = Array{S}(undef, N)         # ζ-axis component
    # Acceleration due to second zonal harmonic J_2
    F_J2_ξ = Array{S}(undef, N)        # ξ-axis component
    F_J2_η = Array{S}(undef, N)        # η-axis component
    F_J2_ζ = Array{S}(undef, N)        # ζ-axis component
    # Acceleration due to third zonal harmonic J_3
    F_J3_ξ = Array{S}(undef, N)        # ξ-axis component
    F_J3_η = Array{S}(undef, N)        # η-axis component
    F_J3_ζ = Array{S}(undef, N)        # ζ-axis component

    # Unit vectors (ξ, η, ζ) in inertial frame

    # ξ vector
    ξx = Array{S}(undef, N)
    ξy = Array{S}(undef, N)
    ξz = Array{S}(undef, N)
    # η vector
    ηx = Array{S}(undef, N)
    ηy = Array{S}(undef, N)
    ηz = Array{S}(undef, N)
    # Auxiliaries to compute η vector
    ηx1 = Array{S}(undef, N)
    ηy1 = Array{S}(undef, N)
    ηz1 = Array{S}(undef, N)
    ηx2 = Array{S}(undef, N)
    ηy2 = Array{S}(undef, N)
    ηz2 = Array{S}(undef, N)
    # ζ vector
    ζx = Array{S}(undef, N)
    ζy = Array{S}(undef, N)
    ζz = Array{S}(undef, N)
    # Auxiliaries to compute ζ vector
    ζx1 = Array{S}(undef, N)
    ζy1 = Array{S}(undef, N)
    ζz1 = Array{S}(undef, N)
    ζx2 = Array{S}(undef, N)
    ζy2 = Array{S}(undef, N)
    ζz2 = Array{S}(undef, N)

    # Full extended-body accelerations
    accX = zero_q_1
    accY = zero_q_1
    accZ = zero_q_1

    # Rotations to and from Earth, Sun and Moon pole-oriented frames
    local M_ = Array{S}(undef, 3, 3, N)

    local M_[:, :, ea] = t2c_jpl_de430(dsj2k, zero_q_1)

    # Fill first 3 elements of dq with velocities
    dq[1] = q[4]
    dq[2] = q[5]
    dq[3] = q[6]

    # Newtonian potential of N bodies
    newtonianNb_Potential[N] = zero_q_1

    #=
    Compute point-mass Newtonian accelerations, all bodies
    See equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    =#
    _4dq[1] = 4dq[1]
    _4dq[2] = 4dq[2]
    _4dq[3] = 4dq[3]
    Threads.@threads for i in 1:Nm1
        # Velocity of the i-th body
        ui[i] = ss16asteph_t[3(N-1+i)-2]    # X-axis component
        vi[i] = ss16asteph_t[3(N-1+i)-1]    # Y-axis component
        wi[i] = ss16asteph_t[3(N-1+i)  ]    # Z-axis component

        # Position of the i-th body - position of the asteroid
        X[i] = ss16asteph_t[3i-2]-q[1]      # X-axis component
        Y[i] = ss16asteph_t[3i-1]-q[2]      # Y-axis component
        Z[i] = ss16asteph_t[3i  ]-q[3]      # Z-axis component

        # Velocity of the i-th body - velocity of the asteroid
        U[i] = ui[i]-dq[1]                  # X-axis component
        V[i] = vi[i]-dq[2]                  # Y-axis component
        W[i] = wi[i]-dq[3]                  # Z-axis component

        # 4 * Velocity of the asteroid - 3 * velocity of the i-th body
        _4U_m_3X[i] = (-3ui[i]) + (_4dq[1])  # X-axis component
        _4V_m_3Y[i] = (-3vi[i]) + (_4dq[2])  # Y-axis component
        _4W_m_3Z[i] = (-3wi[i]) + (_4dq[3])  # Z-axis component

        # Dot product inside the [] in the second term
        pn2x = X[i]*_4U_m_3X[i]
        pn2y = Y[i]*_4V_m_3Y[i]
        pn2z = Z[i]*_4W_m_3Z[i]

        # v_{ij}v_{asteroid} j = x, y, z
        UU[i] = ui[i]*dq[1]
        VV[i] = vi[i]*dq[2]
        WW[i] = wi[i]*dq[3]

        # <Velocity of the i-th body, velocity of the asteroid>
        vi_dot_vj[i] = ( UU[i]+VV[i] ) + WW[i]

        # Distance between the i-th body and the asteroid
        r_p2[i] = ( (X[i]^2)+(Y[i]^2) ) + (Z[i]^2)  # r_{i,asteroid}^2
        r_p1d2[i] = sqrt(r_p2[i])                   # sqrt(r_p2) <-> r_{i,asteroid}
        r_p3d2[i] = r_p2[i]^1.5                     # r_p2^1.5 <-> r_{i, asteroid}^3
        r_p7d2[i] = r_p2[i]^3.5                     # r_p2^3.5 <-> r_{i, asteroid}^7

        # Newtonian coefficient, i.e., mass parameter / distance^3 -> \mu_i / r_{i, asteroid}^3
        newtonianCoeff[i] =  μ[i]/r_p3d2[i]

        # Second term without (\mathbf{v}_i - \mathbf{v}_asteroid)
        pn2[i] = newtonianCoeff[i]*(( pn2x+pn2y ) + pn2z)

        # Newtonian coefficient * difference between two positions, i.e.,
        # \mu_i * (\mathbf{r_i} - \mathbf{r_asteroid}) / r_{ij}^3
        newton_acc_X[i] = X[i]*newtonianCoeff[i]
        newton_acc_Y[i] = Y[i]*newtonianCoeff[i]
        newton_acc_Z[i] = Z[i]*newtonianCoeff[i]

        # Newtonian potential of 1 body \mu_i / r_{i, asteroid}
        newtonian1b_Potential[i] = μ[i]/r_p1d2[i]
        # Third term without newtonian accelerations \mathbf{a}_i
        pn3[i] = 3.5newtonian1b_Potential[i]
        # Full second term
        U_t_pn2[i] = pn2[i]*U[i]    # X-axis component
        V_t_pn2[i] = pn2[i]*V[i]    # Y-axis component
        W_t_pn2[i] = pn2[i]*W[i]    # Z-axis component

        #=
        Extended body accelerations
        See equation (28) in page 13 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
        and equations (173) and (174) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
        =#

        # J_2 accelerations, if i-th body is flattened
        if UJ_interaction[i]
            # Rotate from inertial frame to extended-body frame
            # Here we are rotating only the Z-coordinate
            t31[i] = -X[i]*M_[1,3,i]
            t32[i] = -Y[i]*M_[2,3,i]
            t33[i] = -Z[i]*M_[3,3,i]
            r_sin_ϕ[i] = (t31[i]+t32[i])+t33[i]   # z-coordinate in body-fixed frame

            # Trigonometric functions of latitude ϕ and longitude λ in the body-fixed coordinate system
            # See equations (165)-(168) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract

            sin_ϕ[i] = r_sin_ϕ[i]/r_p1d2[i]       # sin(latitude ϕ)
            ϕ[i] = asin(sin_ϕ[i])                 # Latitude ϕ
            cos_ϕ[i] = cos(ϕ[i])                  # cos(latitude ϕ)
            sin2_ϕ[i] = sin_ϕ[i]^2                # sin(latitude ϕ)^2
            sin3_ϕ[i] = sin_ϕ[i]^3                # sin(latitude ϕ)^3

            # Legendre polynomials

            # See equations (176) and (177) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
            P_2_sin_ϕ[i] = 1.5sin2_ϕ[i] - 0.5               # Second Legendre polynomial P_2(sin ϕ)
            ∂P_2_sin_ϕ[i] = 3sin_ϕ[i]                       # dP_2(sin ϕ)/d(sin ϕ)
            P_3_sin_ϕ[i] = (-1.5sin_ϕ[i]) + (2.5sin3_ϕ[i])  # Third Legendre polynomial P_3(sin ϕ)
            ∂P_3_sin_ϕ[i] = -1.5 + 7.5sin2_ϕ[i]             # dP_3(sin ϕ)/d(sin ϕ)

            # Compute cartesian coordinates of acceleration due to body figure in body frame

            # -J_n * R^n / r^m
            # J_n: n-th zonal harmonic coefficient
            # R: radius of the body
            # r: distance between the body and the asteroid
            Λ2j_div_r4[i] = -(Λ2[i]/(r_p2[i]^2))    # J_2 * R^2 / r^4
            Λ3j_div_r5[i] = -(Λ3[i]/(r_p1d2[i]^5))  # J_3 * R^3 / r^5

            # -cos ϕ P_n'
            m_c_ϕ_∂P_2[i] = (-cos_ϕ[i])*∂P_2_sin_ϕ[i]   # -cos ϕ P_2'
            m_c_ϕ_∂P_3[i] = (-cos_ϕ[i])*∂P_3_sin_ϕ[i]   # -cos ϕ P_3'

            # Acceleration due to second zonal harmonic J_2 in body frame
            F_J2_ξ[i] = ( Λ2j_div_r4[i]*(3P_2_sin_ϕ[i]) )   # ξ-axis component
            # F_J2_η[i] = zero_q_1                          # η-axis component
            F_J2_ζ[i] = Λ2j_div_r4[i]*m_c_ϕ_∂P_2[i]         # ζ-axis component

            # Acceleration due to third zonal harmonic J_3 in body frame
            F_J3_ξ[i] = ( Λ3j_div_r5[i]*(4P_3_sin_ϕ[i]) )   # ξ-axis component
            #F_J3_η[i] = zero_q_1                           # η-axis component
            F_J3_ζ[i] = Λ3j_div_r5[i]*m_c_ϕ_∂P_3[i]         # ζ-axis component

            # Compute accelerations due to zonal harmonics J_n, n = 2, 3 in body frame
            F_J_ξ[i] = F_J2_ξ[i] # + F_J3_ξ[i]              # ξ-axis component
            # F_J_η[i] = zero_q_1                           # η-axis component
            F_J_ζ[i] = F_J2_ζ[i] # + F_J3_ζ[i]              # ζ-axis component

            # Compute unit vectors (ξ, η, ζ) in inertial frame

            # ξ components in inertial frame
            ξx[i] = -X[i]/r_p1d2[i]
            ξy[i] = -Y[i]/r_p1d2[i]
            ξz[i] = -Z[i]/r_p1d2[i]

            # Compute η = p x ξ
            # Auxiliaries
            ηx1[i] = M_[2,3,i]*ξz[i]
            ηy1[i] = M_[3,3,i]*ξx[i]
            ηz1[i] = M_[1,3,i]*ξy[i]
            ηx2[i] = M_[3,3,i]*ξy[i]
            ηy2[i] = M_[1,3,i]*ξz[i]
            ηz2[i] = M_[2,3,i]*ξx[i]
            # η components in inertial frame
            ηx[i] = ηx1[i] - ηx2[i]
            ηy[i] = ηy1[i] - ηy2[i]
            ηz[i] = ηz1[i] - ηz2[i]

            # Compute ζ = ξ x η
            ζx1[i] = ξy[i]*ηz[i]
            ζy1[i] = ξz[i]*ηx[i]
            ζz1[i] = ξx[i]*ηy[i]
            ζx2[i] = ξz[i]*ηy[i]
            ζy2[i] = ξx[i]*ηz[i]
            ζz2[i] = ξy[i]*ηx[i]
            # ζ components in inertial frame
            ζx[i] = ζx1[i] - ζx2[i]
            ζy[i] = ζy1[i] - ζy2[i]
            ζz[i] = ζz1[i] - ζz2[i]

            # Compute cartesian coordinates of acceleration due to body figure in inertial frame
            # Auxiliaries
            F_J2_x1[i] = F_J_ξ[i]*ξx[i]
            F_J2_y1[i] = F_J_ξ[i]*ξy[i]
            F_J2_z1[i] = F_J_ξ[i]*ξz[i]
            F_J2_x2[i] = F_J_ζ[i]*ζx[i]
            F_J2_y2[i] = F_J_ζ[i]*ζy[i]
            F_J2_z2[i] = F_J_ζ[i]*ζz[i]
            # Acceleration due to zonal harmonics in inertial frame
            F_J2_x[i] = F_J2_x1[i] + F_J2_x2[i]
            F_J2_y[i] = F_J2_y1[i] + F_J2_y2[i]
            F_J2_z[i] = F_J2_z1[i] + F_J2_z2[i]
        end
        # Velocity magnitude of the i-th body
        v2[i] = ( (ui[i]^2)+(vi[i]^2) ) + (wi[i]^2)
    end
    # Asteroid velocity magnitude
    v2[N] = ( (q[4]^2)+(q[5]^2) ) + (q[6]^2)

    for i in 1:Nm1
        # Newtonian potential of N bodies
        # \sum_{i\neq l} \frac{\mu_i}{r_{il}} or
        # \sum_{j\neq k} \frac{\mu_j}{r_{jk}}
        temp_004 = newtonian1b_Potential[i] + newtonianNb_Potential[N]
        newtonianNb_Potential[N] = temp_004

        # Extended body accelerations
        # J_n accelerations, if i-th body is flattened
        if UJ_interaction[i]
            # Reaction force on i-th body
            temp_accX_i[i] = accX - (μ[i]*F_J2_x[i])
            accX = temp_accX_i[i]
            temp_accY_i[i] = accY - (μ[i]*F_J2_y[i])
            accY = temp_accY_i[i]
            temp_accZ_i[i] = accZ - (μ[i]*F_J2_z[i])
            accZ = temp_accZ_i[i]
        end
    end

    #=
    Post-Newtonian accelerations due to Sun, Moon and planets (Mercury through Neptune)
    Post-Newtonian iterative procedure setup and initialization
    See equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    =#

    # First term

    # 4*\sum term inside {}
    _4ϕj[N] = 4newtonianNb_Potential[N]
    Threads.@threads for i in 1:10
        # 4*\sum + \sum terms inside {}
        ϕi_plus_4ϕj[i] = newtonianNb_Potential_t[i] + _4ϕj[N]
        # \dot{s}_j^2 + 2\dot{s}_i^2 - 4 <, > terms inside {}
        sj2_plus_2si2_minus_4vivj[i] = ( (2v2[i]) - (4vi_dot_vj[i]) ) + v2[N]
        # -4\sum - \sum + \dot{s}_j^2 + 2\dot{s}_i^2  - 4<, > terms inside {}
        ϕs_and_vs[i] = sj2_plus_2si2_minus_4vivj[i] - ϕi_plus_4ϕj[i]
        # (\mathbf{r}_i - \mathbf{r}_asteroid)\cdot\mathbf{v_i}
        Xij_t_Ui = X[i]*ui[i]
        Yij_t_Vi = Y[i]*vi[i]
        Zij_t_Wi = Z[i]*wi[i]
        Rij_dot_Vi = ( Xij_t_Ui+Yij_t_Vi ) + Zij_t_Wi
        # The expression below inside the (...)^2 should have a minus sign in front of the
        # numerator, but upon squaring it is eliminated, so at the end of the day, it is
        # irrelevant ;)
        # (\mathbf{r}_i - \mathbf{r}_asteroid)\cdot\mathbf{v_i} / r_{i, asteroid}
        pn1t7 = (Rij_dot_Vi^2)/r_p2[i]
        # Everything inside the {} except for the first and last terms
        pn1t2_7 = ϕs_and_vs[i] - (1.5pn1t7)
        # Everything inside the {} except for the last term
        pn1t1_7[i] = c_p2 + pn1t2_7

        # Last term inside the {}
        pNX_t_X[i] = acceph_t[3i-2]*X[i]   # X-axis component
        pNY_t_Y[i] = acceph_t[3i-1]*Y[i]   # Y-axis component
        pNZ_t_Z[i] = acceph_t[3i  ]*Z[i]   # Z-axis component

        # Everything inside the {} in the first term
        pn1[i] = (  pn1t1_7[i]  +  (0.5*( (pNX_t_X[i]+pNY_t_Y[i]) + pNZ_t_Z[i] ))  )

        # Full first term
        X_t_pn1[i] = newton_acc_X[i]*pn1[i]   # X-axis component
        Y_t_pn1[i] = newton_acc_Y[i]*pn1[i]   # Y-axis component
        Z_t_pn1[i] = newton_acc_Z[i]*pn1[i]   # Z-axis component

        # Full third term
        pNX_t_pn3[i] = acceph_t[3i-2]*pn3[i]   # X-axis component
        pNY_t_pn3[i] = acceph_t[3i-1]*pn3[i]   # Y-axis component
        pNZ_t_pn3[i] = acceph_t[3i  ]*pn3[i]   # Z-axis component
    end
    # Temporary post-Newtonian accelerations (planets)
    for i in 1:10
        termpnx = ( X_t_pn1[i] + (U_t_pn2[i]+pNX_t_pn3[i]) )   # X-axis component
        sumpnx = pntempX + termpnx
        pntempX = sumpnx
        termpny = ( Y_t_pn1[i] + (V_t_pn2[i]+pNY_t_pn3[i]) )   # Y-axis component
        sumpny = pntempY + termpny
        pntempY = sumpny
        termpnz = ( Z_t_pn1[i] + (W_t_pn2[i]+pNZ_t_pn3[i]) )   # Z-axis component
        sumpnz = pntempZ + termpnz
        pntempZ = sumpnz
    end
    # Compute Newtonian accelerations due to Pluto and 16 asteroid perturbers
    Threads.@threads for i in 11:Nm1
        # Full first term
        X_t_pn1[i] = c_p2*newton_acc_X[i]
        Y_t_pn1[i] = c_p2*newton_acc_Y[i]
        Z_t_pn1[i] = c_p2*newton_acc_Z[i]
    end
    # Temporary post-Newtonian accelerations (Pluto + 16 asteroid perturbers)
    for i in 11:Nm1
        termpnx = X_t_pn1[i]          # X-axis component
        sumpnx = pntempX + termpnx
        pntempX = sumpnx
        termpny = Y_t_pn1[i]          # Y-axis component
        sumpny = pntempY + termpny
        pntempY = sumpny
        termpnz = Z_t_pn1[i]          # Z-axis component
        sumpnz = pntempZ + termpnz
        pntempZ = sumpnz
    end
     # Post-Newtonian acelerations
    postNewtonX = pntempX*c_m2     # X-axis component
    postNewtonY = pntempY*c_m2     # Y-axis component
    postNewtonZ = pntempZ*c_m2     # Z-axis component

    #=
    Compute non-gravitational acceleration (Marsden et al., 1973)
    =#

    # Angular momentum per unit mass (normal vector)
    h_x = (Y[1]*W[1])-(Z[1]*V[1])   # X-axis component
    h_y = (Z[1]*U[1])-(X[1]*W[1])   # Y-axis component
    h_z = (X[1]*V[1])-(Y[1]*U[1])   # Z-axis component
    # Norm of normal vector
    h_norm = sqrt( ((h_x^2)+(h_y^2))+(h_z^2) )
    # Cartesian components of normal unit vector
    h_x_unit = h_x/h_norm
    h_y_unit = h_y/h_norm
    h_z_unit = h_z/h_norm

    # Cartesian components of transversal vector t = h × (\mathbf{r}_Sun - \mathbf{r}_asteroid)
    t_x = (h_z*Y[1]) - (h_y*Z[1])    # Note: Y[1] = y_Sun - y_asteroid, etc.
    t_y = (h_x*Z[1]) - (h_z*X[1])
    t_z = (h_y*X[1]) - (h_x*Y[1])
    # Norm of transversal vector
    t_norm = sqrt( ((t_x^2)+(t_y^2))+(t_z^2) )
    # Cartesian components of transversal unit vector
    t_x_unit = t_x/t_norm
    t_y_unit = t_y/t_norm
    t_z_unit = t_z/t_norm

    # Cartesian components of radial unit vector
    r_x_unit = -(X[1]/r_p1d2[1])
    r_y_unit = -(Y[1]/r_p1d2[1])
    r_z_unit = -(Z[1]/r_p1d2[1])

    # Marsden et al. (1973) radial function
    g_r_quotient = r_p1d2[1] / marsden_r₀
    g_r_A = marsden_α * (g_r_quotient^marsden_m)
    g_r_B = (1 + (g_r_quotient^marsden_n))^marsden_k
    g_r = g_r_A * g_r_B

    # Evaluate non-grav acceleration
    A2_t_g_r = g_r*q[7]     # Yarkovsky effect
    A1_t_g_r = g_r*q[8]     # Radiation pressure
    A3_t_g_r = g_r*q[9]     # Normal non-gravitational component

    # Non gravitational acceleration: Yarkovsky + radiation pressure + normal component
    NGAx = ((A2_t_g_r*t_x_unit) + (A1_t_g_r*r_x_unit)) + (A3_t_g_r)*h_x_unit
    NGAy = ((A2_t_g_r*t_y_unit) + (A1_t_g_r*r_y_unit)) + (A3_t_g_r)*h_y_unit
    NGAz = ((A2_t_g_r*t_z_unit) + (A1_t_g_r*r_z_unit)) + (A3_t_g_r)*h_z_unit

    # Fill dq[4:6] with accelerations
    # Post-Newton point mass + Extended body + Non-gravitational
    dq[4] = ( postNewtonX + accX ) + NGAx
    dq[5] = ( postNewtonY + accY ) + NGAy
    dq[6] = ( postNewtonZ + accZ ) + NGAz
    # Nongrav acceleration coefficients do not change in time
    dq[7] = zero_q_1
    dq[8] = zero_q_1
    dq[9] = zero_q_1

    nothing
end

"""
    gravityonly!

Asteroid dynamical model specially suited for long-term integrations and multi-arc orbit
determination.

The model considers the asteroid of interest as a test particle with null mass. Perturbing
bodies included in the model are: the Sun, the eight planets, the Moon, Pluto and the 16
most massive asteroids. Planetary ephemerides are provided by `PlanetaryEphemeris.jl`,
which is based on the JPL DE430 model.

Dynamical effects considered are:

- Post-Newtonian point-mass accelerations between all bodies. See equation (35) in page
    7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract.

- Figure-effects (oblateness) of the Earth (``J_2`` and ``J_3``), ``J_2`` effect of the
    Sun and ``J_2`` and ``J_3`` effect of the Moon. See equation (28) in page 13 of
    https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract and equations (173)
    and (174) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract.

- Kinematic model for the precession and nutation of the Earth's orientation (IAU 1976/1980
    Earth orientation model). See [`PlanetaryEphemeris.c2t_jpl_de430`](@ref).

- Kinematic model for the Moon's orientation (Seidelmann et al., 2006). See equations
    (14)-(15) in page 9 and equations (34)-(35) in page 16 of
    https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract.

To improve performance, some internal loops are multi-threaded via `@threads`.

For other dynamical models, see [`nongravs!`](@ref) and [`newtonian!`](@ref).
"""
function gravityonly!(dq, q, params, t)
    # Julian date (TDB) of start time
    local jd0 = params.jd0
    # Days since J2000.0 = 2.451545e6
    local dsj2k = t + (jd0 - JD_J2000)
    # Solar system ephemeris at dsj2k
    local ss16asteph_t = params.ssepht
    evaleph!(ss16asteph_t, params.sseph, dsj2k)
    # Accelerations at dsj2k
    local acceph_t = params.accepht
    evaleph!(acceph_t, params.acceph, dsj2k)
    # Newtonian potentials at dsj2k
    local newtonianNb_Potential_t = params.potepht
    evaleph!(newtonianNb_Potential_t, params.poteph, dsj2k)
    # Type of position / velocity components
    local S = eltype(q)
    # Interaction matrix with flattened bodies
    local UJ_interaction = params.UJ_interaction
    # Number of bodies, including NEA
    local N = params.N
    # Number of bodies, except the asteroid
    local Nm1 = N-1
    # Vector of mass parameters GM's
    local μ = params.μ

    # zero(q[1])
    local zero_q_1 = auxzero(q[1])

    #=
    Point-mass accelerations
    See equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    =#

    # Position of the i-th body - position of the asteroid
    X = Array{S}(undef, N)         # X-axis component
    Y = Array{S}(undef, N)         # Y-axis component
    Z = Array{S}(undef, N)         # Z-axis component

    # Distance between the i-th body and the asteroid
    r_p2 = Array{S}(undef, N)      # r_{i,asteroid}^2
    r_p1d2 = Array{S}(undef, N)    # sqrt(r_p2) <-> r_{i, asteroid}
    r_p3d2 = Array{S}(undef, N)    # r_p2^1.5 <-> r_{i, asteroid}^3
    r_p7d2 = Array{S}(undef, N)    # r_p2^3.5 <-> r_{i, asteroid}^7

    # Newtonian coefficient, i.e., mass parameter / distance^3 -> \mu_i / r_{i, asteroid}^3
    newtonianCoeff = Array{S}(undef, N)

    # Velocity of the i-th body
    ui = Array{S}(undef, N-1)      # X-axis component
    vi = Array{S}(undef, N-1)      # Y-axis component
    wi = Array{S}(undef, N-1)      # Z-axis component

    # Post-Newtonian stuff

    # Velocity of the i-th body - velocity of the asteroid
    U = Array{S}(undef, N)         # X-axis component
    V = Array{S}(undef, N)         # Y-axis component
    W = Array{S}(undef, N)         # Z-axis component

    # 4 * Velocity of the asteroid
    _4dq = Array{S}(undef, 3)  # (X,Y,Z) components

    # 4 * Velocity of the asteroid - 3 * velocity of the i-th body
    _4U_m_3X = Array{S}(undef, N)  # X-axis component
    _4V_m_3Y = Array{S}(undef, N)  # Y-axis component
    _4W_m_3Z = Array{S}(undef, N)  # Z-axis component

    # v_{i,j}v_{asteroid,j} j = x, y, z
    UU = Array{S}(undef, N)        # v_{i,x}v_{asteroid,x}
    VV = Array{S}(undef, N)        # v_{i,y}v_{asteroid,y}
    WW = Array{S}(undef, N)        # v_{i,z}v_{asteroid,z}

    # Newtonian potential of 1 body \mu_i / r_{i, asteroid}
    newtonian1b_Potential = Array{S}(undef, N)
    # Newtonian potential of N bodies
    # \sum_{i\neq l} \frac{\mu_i}{r_{il}} or
    # \sum_{j\neq k} \frac{\mu_j}{r_{jk}}
    newtonianNb_Potential = Array{S}(undef, N)

    # Newtonian coefficient * difference between two positions, i.e.,
    # \mu_i * (\mathbf{r_i} - \mathbf{r_asteroid}) / r_{ij}^3
    newton_acc_X = Array{S}(undef, N)   # X-axis component
    newton_acc_Y = Array{S}(undef, N)   # Y-axis component
    newton_acc_Z = Array{S}(undef, N)   # Z-axis component

    # Combinations of velocities
    v2 = Array{S}(undef, N)             # Velocity magnitude squared ||\mathbf{v}_i||^2
    vi_dot_vj = Array{S}(undef, N)      # <Velocity of the i-th body, velocity of the asteroid>

    # Second term without (\mathbf{v}_i - \mathbf{v}_j)
    pn2 = Array{S}(undef, N)
    # Full second term
    U_t_pn2 = Array{S}(undef, N)        # X-axis component
    V_t_pn2 = Array{S}(undef, N)        # Y-axis component
    W_t_pn2 = Array{S}(undef, N)        # Z-axis component

    # Third term without newtonian accelerations \mathbf{a}_i
    pn3 = Array{S}(undef, N)
    # Full third term of equation (35)
    pNX_t_pn3 = Array{S}(undef, N)      # X-axis component
    pNY_t_pn3 = Array{S}(undef, N)      # Y-axis component
    pNZ_t_pn3 = Array{S}(undef, N)      # Z-axis component

    # First term
    _4ϕj = Array{S}(undef, N)           # 4*\sum term inside {}
    ϕi_plus_4ϕj = Array{S}(undef, N)    # 4*\sum + \sum terms inside {}
    sj2_plus_2si2_minus_4vivj = Array{S}(undef, N)  # \dot{s}_j^2 + 2\dot{s}_i^2 - 4 <, > terms inside {}
    ϕs_and_vs = Array{S}(undef, N)      # -4\sum - \sum + \dot{s}_j^2 + 2\dot{s}_i^2  - 4<, > terms inside {}
    pn1t1_7 = Array{S}(undef, N)        # Everything inside the {} in the first term except for the term with accelerations (last)
    # Last term inside the {}
    pNX_t_X = Array{S}(undef, N)        # X-axis component
    pNY_t_Y = Array{S}(undef, N)        # Y-axis component
    pNZ_t_Z = Array{S}(undef, N)        # Z-axis component
    # Everything inside the {} in the first term
    pn1 = Array{S}(undef, N)
    # Full first term
    X_t_pn1 = Array{S}(undef, N)     # X-axis component
    Y_t_pn1 = Array{S}(undef, N)     # Y-axis component
    Z_t_pn1 = Array{S}(undef, N)     # Z-axis component

    # Temporary post-Newtonian accelerations
    pntempX = zero_q_1               # X-axis component
    pntempY = zero_q_1               # Y-axis component
    pntempZ = zero_q_1               # Z-axis component

    #=
    Extended body accelerations
    See equation (28) in page 13 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    and equations (173) and (174) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    =#

    # J_2 acceleration auxiliaries

    # Auxiliaries to compute body-fixed frame coordinates
    t31 = Array{S}(undef, N)
    t32 = Array{S}(undef, N)
    t33 = Array{S}(undef, N)
    # z-coordinate in body-fixed frame
    r_sin_ϕ = Array{S}(undef, N)

    # Trigonometric functions of latitude ϕ and longitude λ in the body-fixed coordinate system
    # See equations (165)-(168) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    sin_ϕ = Array{S}(undef, N)      # sin(latitude ϕ)
    ϕ = Array{S}(undef, N)          # Latitude ϕ
    cos_ϕ = Array{S}(undef, N)      # cos(latitude ϕ)
    sin2_ϕ = Array{S}(undef, N)     # sin(latitude ϕ)^2
    sin3_ϕ = Array{S}(undef, N)     # sin(latitude ϕ)^3
    sin4_ϕ = Array{S}(undef, N)     # sin(latitude ϕ)^4

    # Acceleration due to zonal harmonics in inertial frame
    F_J2_x = Array{S}(undef, N)
    F_J2_y = Array{S}(undef, N)
    F_J2_z = Array{S}(undef, N)
    # Auxiliaries to compute F_J2_i, i = x, y, z
    F_J2_x1 = Array{S}(undef, N)
    F_J2_y1 = Array{S}(undef, N)
    F_J2_z1 = Array{S}(undef, N)
    F_J2_x2 = Array{S}(undef, N)
    F_J2_y2 = Array{S}(undef, N)
    F_J2_z2 = Array{S}(undef, N)

    # Temporary arrays for the sum of full extended body accelerations
    temp_accX_i = Array{S}(undef, N)
    temp_accY_i = Array{S}(undef, N)
    temp_accZ_i = Array{S}(undef, N)

    # Legendre polynomials
    # See equations (176) and (177) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    P_2_sin_ϕ = Array{S}(undef, N)      # Second Legendre polynomial P_2(sin ϕ)
    ∂P_2_sin_ϕ = Array{S}(undef, N)     # dP_2(sin ϕ)/d(sin ϕ)
    P_3_sin_ϕ = Array{S}(undef, N)      # Third Legendre polynomial P_3(sin ϕ)
    ∂P_3_sin_ϕ = Array{S}(undef, N)     # dP_3(sin ϕ)/d(sin ϕ)
    # -cos ϕ P_n'
    m_c_ϕ_∂P_2 = Array{S}(undef, N)     # -cos ϕ P_2'
    m_c_ϕ_∂P_3 = Array{S}(undef, N)     # -cos ϕ P_3'

    # -J_n * R^n / r^m
    # J_n: n-th zonal harmonic coefficient
    # R: radius of the body
    # r: distance between the body and the asteroid
    Λ2j_div_r4 = Array{S}(undef, N)   # J_2 * R^2 / r^4
    Λ3j_div_r5 = Array{S}(undef, N)   # J_3 * R^3 / r^5

    # Accelerations due to zonal harmonics in body frame

    # Acceleration due to zonal harmonics J_n, n = 2, 3
    F_J_ξ = Array{S}(undef, N)         # ξ-axis component
    F_J_η = Array{S}(undef, N)         # η-axis component
    F_J_ζ = Array{S}(undef, N)         # ζ-axis component
    # Acceleration due to second zonal harmonic J_2
    F_J2_ξ = Array{S}(undef, N)        # ξ-axis component
    F_J2_η = Array{S}(undef, N)        # η-axis component
    F_J2_ζ = Array{S}(undef, N)        # ζ-axis component
    # Acceleration due to third zonal harmonic J_3
    F_J3_ξ = Array{S}(undef, N)        # ξ-axis component
    F_J3_η = Array{S}(undef, N)        # η-axis component
    F_J3_ζ = Array{S}(undef, N)        # ζ-axis component

    # Unit vectors (ξ, η, ζ) in inertial frame

    # ξ vector
    ξx = Array{S}(undef, N)
    ξy = Array{S}(undef, N)
    ξz = Array{S}(undef, N)
    # η vector
    ηx = Array{S}(undef, N)
    ηy = Array{S}(undef, N)
    ηz = Array{S}(undef, N)
    # Auxiliaries to compute η vector
    ηx1 = Array{S}(undef, N)
    ηy1 = Array{S}(undef, N)
    ηz1 = Array{S}(undef, N)
    ηx2 = Array{S}(undef, N)
    ηy2 = Array{S}(undef, N)
    ηz2 = Array{S}(undef, N)
    # ζ vector
    ζx = Array{S}(undef, N)
    ζy = Array{S}(undef, N)
    ζz = Array{S}(undef, N)
    # Auxiliaries to compute ζ vector
    ζx1 = Array{S}(undef, N)
    ζy1 = Array{S}(undef, N)
    ζz1 = Array{S}(undef, N)
    ζx2 = Array{S}(undef, N)
    ζy2 = Array{S}(undef, N)
    ζz2 = Array{S}(undef, N)

    # Full extended-body accelerations
    accX = zero_q_1
    accY = zero_q_1
    accZ = zero_q_1

    # Rotations to and from Earth, Sun and Moon pole-oriented frames
    local M_ = Array{S}(undef, 3, 3, N)

    local M_[:, :, ea] = t2c_jpl_de430(dsj2k, zero_q_1)

    # Fill first 3 elements of dq with velocities
    dq[1] = q[4]
    dq[2] = q[5]
    dq[3] = q[6]

    # Newtonian potential of N bodies
    newtonianNb_Potential[N] = zero_q_1

    #=
    Compute point-mass Newtonian accelerations, all bodies
    See equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    =#
    _4dq[1] = 4dq[1]
    _4dq[2] = 4dq[2]
    _4dq[3] = 4dq[3]
    Threads.@threads for i in 1:Nm1
        # Velocity of the i-th body
        ui[i] = ss16asteph_t[3(N-1+i)-2]    # X-axis component
        vi[i] = ss16asteph_t[3(N-1+i)-1]    # Y-axis component
        wi[i] = ss16asteph_t[3(N-1+i)  ]    # Z-axis component

        # Position of the i-th body - position of the asteroid
        X[i] = ss16asteph_t[3i-2]-q[1]      # X-axis component
        Y[i] = ss16asteph_t[3i-1]-q[2]      # Y-axis component
        Z[i] = ss16asteph_t[3i  ]-q[3]      # Z-axis component

        # Velocity of the i-th body - velocity of the asteroid
        U[i] = ui[i]-dq[1]                  # X-axis component
        V[i] = vi[i]-dq[2]                  # Y-axis component
        W[i] = wi[i]-dq[3]                  # Z-axis component

        # 4 * Velocity of the asteroid - 3 * velocity of the i-th body
        _4U_m_3X[i] = (-3ui[i]) + (_4dq[1])  # X-axis component
        _4V_m_3Y[i] = (-3vi[i]) + (_4dq[2])  # Y-axis component
        _4W_m_3Z[i] = (-3wi[i]) + (_4dq[3])  # Z-axis component

        # Dot product inside the [] in the second term
        pn2x = X[i]*_4U_m_3X[i]
        pn2y = Y[i]*_4V_m_3Y[i]
        pn2z = Z[i]*_4W_m_3Z[i]

        # v_{ij}v_{asteroid} j = x, y, z
        UU[i] = ui[i]*dq[1]
        VV[i] = vi[i]*dq[2]
        WW[i] = wi[i]*dq[3]

        # <Velocity of the i-th body, velocity of the asteroid>
        vi_dot_vj[i] = ( UU[i]+VV[i] ) + WW[i]

        # Distance between the i-th body and the asteroid
        r_p2[i] = ( (X[i]^2)+(Y[i]^2) ) + (Z[i]^2)  # r_{i,asteroid}^2
        r_p1d2[i] = sqrt(r_p2[i])                   # sqrt(r_p2) <-> r_{i,asteroid}
        r_p3d2[i] = r_p2[i]^1.5                     # r_p2^1.5 <-> r_{i, asteroid}^3
        r_p7d2[i] = r_p2[i]^3.5                     # r_p2^3.5 <-> r_{i, asteroid}^7

        # Newtonian coefficient, i.e., mass parameter / distance^3 -> \mu_i / r_{i, asteroid}^3
        newtonianCoeff[i] =  μ[i]/r_p3d2[i]

        # Second term without (\mathbf{v}_i - \mathbf{v}_asteroid)
        pn2[i] = newtonianCoeff[i]*(( pn2x+pn2y ) + pn2z)

        # Newtonian coefficient * difference between two positions, i.e.,
        # \mu_i * (\mathbf{r_i} - \mathbf{r_asteroid}) / r_{ij}^3
        newton_acc_X[i] = X[i]*newtonianCoeff[i]
        newton_acc_Y[i] = Y[i]*newtonianCoeff[i]
        newton_acc_Z[i] = Z[i]*newtonianCoeff[i]

        # Newtonian potential of 1 body \mu_i / r_{i, asteroid}
        newtonian1b_Potential[i] = μ[i]/r_p1d2[i]
        # Third term without newtonian accelerations \mathbf{a}_i
        pn3[i] = 3.5newtonian1b_Potential[i]
        # Full second term
        U_t_pn2[i] = pn2[i]*U[i]    # X-axis component
        V_t_pn2[i] = pn2[i]*V[i]    # Y-axis component
        W_t_pn2[i] = pn2[i]*W[i]    # Z-axis component

        #=
        Extended body accelerations
        See equation (28) in page 13 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
        and equations (173) and (174) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
        =#

        # J_2 accelerations, if i-th body is flattened
        if UJ_interaction[i]
            # Rotate from inertial frame to extended-body frame
            # Here we are rotating only the Z-coordinate
            t31[i] = -X[i]*M_[1,3,i]
            t32[i] = -Y[i]*M_[2,3,i]
            t33[i] = -Z[i]*M_[3,3,i]
            r_sin_ϕ[i] = (t31[i]+t32[i])+t33[i]   # z-coordinate in body-fixed frame

            # Trigonometric functions of latitude ϕ and longitude λ in the body-fixed coordinate system
            # See equations (165)-(168) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract

            sin_ϕ[i] = r_sin_ϕ[i]/r_p1d2[i]       # sin(latitude ϕ)
            ϕ[i] = asin(sin_ϕ[i])                 # Latitude ϕ
            cos_ϕ[i] = cos(ϕ[i])                  # cos(latitude ϕ)
            sin2_ϕ[i] = sin_ϕ[i]^2                # sin(latitude ϕ)^2
            sin3_ϕ[i] = sin_ϕ[i]^3                # sin(latitude ϕ)^3

            # Legendre polynomials

            # See equations (176) and (177) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
            P_2_sin_ϕ[i] = 1.5sin2_ϕ[i] - 0.5               # Second Legendre polynomial P_2(sin ϕ)
            ∂P_2_sin_ϕ[i] = 3sin_ϕ[i]                       # dP_2(sin ϕ)/d(sin ϕ)
            P_3_sin_ϕ[i] = (-1.5sin_ϕ[i]) + (2.5sin3_ϕ[i])  # Third Legendre polynomial P_3(sin ϕ)
            ∂P_3_sin_ϕ[i] = -1.5 + 7.5sin2_ϕ[i]             # dP_3(sin ϕ)/d(sin ϕ)

            # Compute cartesian coordinates of acceleration due to body figure in body frame

            # -J_n * R^n / r^m
            # J_n: n-th zonal harmonic coefficient
            # R: radius of the body
            # r: distance between the body and the asteroid
            Λ2j_div_r4[i] = -(Λ2[i]/(r_p2[i]^2))    # J_2 * R^2 / r^4
            Λ3j_div_r5[i] = -(Λ3[i]/(r_p1d2[i]^5))  # J_3 * R^3 / r^5

            # -cos ϕ P_n'
            m_c_ϕ_∂P_2[i] = (-cos_ϕ[i])*∂P_2_sin_ϕ[i]   # -cos ϕ P_2'
            m_c_ϕ_∂P_3[i] = (-cos_ϕ[i])*∂P_3_sin_ϕ[i]   # -cos ϕ P_3'

            # Acceleration due to second zonal harmonic J_2 in body frame
            F_J2_ξ[i] = ( Λ2j_div_r4[i]*(3P_2_sin_ϕ[i]) )   # ξ-axis component
            # F_J2_η[i] = zero_q_1                          # η-axis component
            F_J2_ζ[i] = Λ2j_div_r4[i]*m_c_ϕ_∂P_2[i]         # ζ-axis component

            # Acceleration due to third zonal harmonic J_3 in body frame
            F_J3_ξ[i] = ( Λ3j_div_r5[i]*(4P_3_sin_ϕ[i]) )   # ξ-axis component
            #F_J3_η[i] = zero_q_1                           # η-axis component
            F_J3_ζ[i] = Λ3j_div_r5[i]*m_c_ϕ_∂P_3[i]         # ζ-axis component

            # Compute accelerations due to zonal harmonics J_n, n = 2, 3 in body frame
            F_J_ξ[i] = F_J2_ξ[i] # + F_J3_ξ[i]              # ξ-axis component
            # F_J_η[i] = zero_q_1                           # η-axis component
            F_J_ζ[i] = F_J2_ζ[i] # + F_J3_ζ[i]              # ζ-axis component

            # Compute unit vectors (ξ, η, ζ) in inertial frame

            # ξ components in inertial frame
            ξx[i] = -X[i]/r_p1d2[i]
            ξy[i] = -Y[i]/r_p1d2[i]
            ξz[i] = -Z[i]/r_p1d2[i]

            # Compute η = p x ξ
            # Auxiliaries
            ηx1[i] = M_[2,3,i]*ξz[i]
            ηy1[i] = M_[3,3,i]*ξx[i]
            ηz1[i] = M_[1,3,i]*ξy[i]
            ηx2[i] = M_[3,3,i]*ξy[i]
            ηy2[i] = M_[1,3,i]*ξz[i]
            ηz2[i] = M_[2,3,i]*ξx[i]
            # η components in inertial frame
            ηx[i] = ηx1[i] - ηx2[i]
            ηy[i] = ηy1[i] - ηy2[i]
            ηz[i] = ηz1[i] - ηz2[i]

            # Compute ζ = ξ x η
            ζx1[i] = ξy[i]*ηz[i]
            ζy1[i] = ξz[i]*ηx[i]
            ζz1[i] = ξx[i]*ηy[i]
            ζx2[i] = ξz[i]*ηy[i]
            ζy2[i] = ξx[i]*ηz[i]
            ζz2[i] = ξy[i]*ηx[i]
            # ζ components in inertial frame
            ζx[i] = ζx1[i] - ζx2[i]
            ζy[i] = ζy1[i] - ζy2[i]
            ζz[i] = ζz1[i] - ζz2[i]

            # Compute cartesian coordinates of acceleration due to body figure in inertial frame
            # Auxiliaries
            F_J2_x1[i] = F_J_ξ[i]*ξx[i]
            F_J2_y1[i] = F_J_ξ[i]*ξy[i]
            F_J2_z1[i] = F_J_ξ[i]*ξz[i]
            F_J2_x2[i] = F_J_ζ[i]*ζx[i]
            F_J2_y2[i] = F_J_ζ[i]*ζy[i]
            F_J2_z2[i] = F_J_ζ[i]*ζz[i]
            # Acceleration due to zonal harmonics in inertial frame
            F_J2_x[i] = F_J2_x1[i] + F_J2_x2[i]
            F_J2_y[i] = F_J2_y1[i] + F_J2_y2[i]
            F_J2_z[i] = F_J2_z1[i] + F_J2_z2[i]
        end
        # Velocity magnitude of the i-th body
        v2[i] = ( (ui[i]^2)+(vi[i]^2) ) + (wi[i]^2)
    end
    # Asteroid velocity magnitude
    v2[N] = ( (q[4]^2)+(q[5]^2) ) + (q[6]^2)

    for i in 1:Nm1
        # Newtonian potential of N bodies
        # \sum_{i\neq l} \frac{\mu_i}{r_{il}} or
        # \sum_{j\neq k} \frac{\mu_j}{r_{jk}}
        temp_004 = newtonian1b_Potential[i] + newtonianNb_Potential[N]
        newtonianNb_Potential[N] = temp_004

        # Extended body accelerations
        # J_n accelerations, if i-th body is flattened
        if UJ_interaction[i]
            # Reaction force on i-th body
            temp_accX_i[i] = accX - (μ[i]*F_J2_x[i])
            accX = temp_accX_i[i]
            temp_accY_i[i] = accY - (μ[i]*F_J2_y[i])
            accY = temp_accY_i[i]
            temp_accZ_i[i] = accZ - (μ[i]*F_J2_z[i])
            accZ = temp_accZ_i[i]
        end
    end

    #=
    Post-Newtonian accelerations due to Sun, Moon and planets (Mercury through Neptune)
    Post-Newtonian iterative procedure setup and initialization
    See equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    =#

    # First term

    # 4*\sum term inside {}
    _4ϕj[N] = 4newtonianNb_Potential[N]
    Threads.@threads for i in 1:10
        # 4*\sum + \sum terms inside {}
        ϕi_plus_4ϕj[i] = newtonianNb_Potential_t[i] + _4ϕj[N]
        # \dot{s}_j^2 + 2\dot{s}_i^2 - 4 <, > terms inside {}
        sj2_plus_2si2_minus_4vivj[i] = ( (2v2[i]) - (4vi_dot_vj[i]) ) + v2[N]
        # -4\sum - \sum + \dot{s}_j^2 + 2\dot{s}_i^2  - 4<, > terms inside {}
        ϕs_and_vs[i] = sj2_plus_2si2_minus_4vivj[i] - ϕi_plus_4ϕj[i]
        # (\mathbf{r}_i - \mathbf{r}_asteroid)\cdot\mathbf{v_i}
        Xij_t_Ui = X[i]*ui[i]
        Yij_t_Vi = Y[i]*vi[i]
        Zij_t_Wi = Z[i]*wi[i]
        Rij_dot_Vi = ( Xij_t_Ui+Yij_t_Vi ) + Zij_t_Wi
        # The expression below inside the (...)^2 should have a minus sign in front of the
        # numerator, but upon squaring it is eliminated, so at the end of the day, it is
        # irrelevant ;)
        # (\mathbf{r}_i - \mathbf{r}_asteroid)\cdot\mathbf{v_i} / r_{i, asteroid}
        pn1t7 = (Rij_dot_Vi^2)/r_p2[i]
        # Everything inside the {} except for the first and last terms
        pn1t2_7 = ϕs_and_vs[i] - (1.5pn1t7)
        # Everything inside the {} except for the last term
        pn1t1_7[i] = c_p2 + pn1t2_7

        # Last term inside the {}
        pNX_t_X[i] = acceph_t[3i-2]*X[i]   # X-axis component
        pNY_t_Y[i] = acceph_t[3i-1]*Y[i]   # Y-axis component
        pNZ_t_Z[i] = acceph_t[3i  ]*Z[i]   # Z-axis component

        # Everything inside the {} in the first term
        pn1[i] = (  pn1t1_7[i]  +  (0.5*( (pNX_t_X[i]+pNY_t_Y[i]) + pNZ_t_Z[i] ))  )

        # Full first term
        X_t_pn1[i] = newton_acc_X[i]*pn1[i]   # X-axis component
        Y_t_pn1[i] = newton_acc_Y[i]*pn1[i]   # Y-axis component
        Z_t_pn1[i] = newton_acc_Z[i]*pn1[i]   # Z-axis component

        # Full third term
        pNX_t_pn3[i] = acceph_t[3i-2]*pn3[i]   # X-axis component
        pNY_t_pn3[i] = acceph_t[3i-1]*pn3[i]   # Y-axis component
        pNZ_t_pn3[i] = acceph_t[3i  ]*pn3[i]   # Z-axis component
    end
    # Temporary post-Newtonian accelerations (planets)
    for i in 1:10
        termpnx = ( X_t_pn1[i] + (U_t_pn2[i]+pNX_t_pn3[i]) )   # X-axis component
        sumpnx = pntempX + termpnx
        pntempX = sumpnx
        termpny = ( Y_t_pn1[i] + (V_t_pn2[i]+pNY_t_pn3[i]) )   # Y-axis component
        sumpny = pntempY + termpny
        pntempY = sumpny
        termpnz = ( Z_t_pn1[i] + (W_t_pn2[i]+pNZ_t_pn3[i]) )   # Z-axis component
        sumpnz = pntempZ + termpnz
        pntempZ = sumpnz
    end
    # Compute Newtonian accelerations due to Pluto and 16 asteroid perturbers
    Threads.@threads for i in 11:Nm1
        # Full first term
        X_t_pn1[i] = c_p2*newton_acc_X[i]
        Y_t_pn1[i] = c_p2*newton_acc_Y[i]
        Z_t_pn1[i] = c_p2*newton_acc_Z[i]
    end
    # Temporary post-Newtonian accelerations (Pluto + 16 asteroid perturbers)
    for i in 11:Nm1
        termpnx = X_t_pn1[i]          # X-axis component
        sumpnx = pntempX + termpnx
        pntempX = sumpnx
        termpny = Y_t_pn1[i]          # Y-axis component
        sumpny = pntempY + termpny
        pntempY = sumpny
        termpnz = Z_t_pn1[i]          # Z-axis component
        sumpnz = pntempZ + termpnz
        pntempZ = sumpnz
    end
     # Post-Newtonian acelerations
    postNewtonX = pntempX*c_m2     # X-axis component
    postNewtonY = pntempY*c_m2     # Y-axis component
    postNewtonZ = pntempZ*c_m2     # Z-axis component

    # Fill dq[4:6] with accelerations
    # Post-Newton point mass + Extended body
    dq[4] = postNewtonX + accX
    dq[5] = postNewtonY + accY
    dq[6] = postNewtonZ + accZ

    nothing
end

"""
    newtonian!

Asteroid dynamical model specially suited for short-term integrations and initial
orbit determination.

The model considers the asteroid of interest as a test particle with null mass.
Perturbing bodies included in the model are: the Sun, the eight planets and the
Moon. Planetary ephemerides are provided by `PlanetaryEphemeris.jl`, which is
based on the JPL DE430 model.

Dynamical effects considered are:

- Newtonian point-mass accelerations between all bodies.

To improve performance, some internal loops are multi-threaded via `@threads`.

For other dynamical models, see [`nongravs!`](@ref) and [`gravityonly!`](@ref).
"""
function newtonian!(dq, q, params, t)
    # Julian date (TDB) of start time
    local jd0 = params.jd0
    # Days since J2000.0 = 2.451545e6
    local dsj2k = t + (jd0 - JD_J2000)
    # Solar system ephemeris at dsj2k
    local ss16asteph_t = params.ssepht
    evaleph!(ss16asteph_t, params.sseph, dsj2k)
    # Type of position / velocity components
    local S = eltype(q)
    # Number of bodies, including NEA
    local N = 10 # Sun, Moon and planets # params.N
    # Number of bodies, except the asteroid
    local Nm1 = N-1
    # Vector of mass parameters GM's
    local μ = params.μ[1:10]

    # zero(q[1])
    local zero_q_1 = auxzero(q[1])

    #=
    Point-mass accelerations
    See equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    =#

    # Position of the i-th body - position of the asteroid
    X = Array{S}(undef, N)         # X-axis component
    Y = Array{S}(undef, N)         # Y-axis component
    Z = Array{S}(undef, N)         # Z-axis component

    # Distance between the i-th body and the asteroid
    r_p2 = Array{S}(undef, N)      # r_{i,asteroid}^2
    r_p3d2 = Array{S}(undef, N)    # r_p2^1.5 <-> r_{i, asteroid}^3

    # Newtonian coefficient, i.e., mass parameter / distance^3 -> \mu_i / r_{i, asteroid}^3
    newtonianCoeff = Array{S}(undef, N)

    # Newtonian coefficient * difference between two positions, i.e.,
    # \mu_i * (\mathbf{r_i} - \mathbf{r_asteroid}) / r_{ij}^3
    newton_acc_X = Array{S}(undef, N)   # X-axis component
    newton_acc_Y = Array{S}(undef, N)   # Y-axis component
    newton_acc_Z = Array{S}(undef, N)   # Z-axis component

    # Temporary arrays for the sum of full extended body accelerations
    temp_accX_i = Array{S}(undef, N)
    temp_accY_i = Array{S}(undef, N)
    temp_accZ_i = Array{S}(undef, N)

    # Full extended-body accelerations
    accX = zero_q_1
    accY = zero_q_1
    accZ = zero_q_1

    # Fill first 3 elements of dq with velocities
    dq[1] = q[4]
    dq[2] = q[5]
    dq[3] = q[6]

    #=
    Compute point-mass Newtonian accelerations, all bodies
    See equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    =#
    Threads.@threads for i in 1:Nm1
        # Position of the i-th body - position of the asteroid
        X[i] = ss16asteph_t[3i-2]-q[1]      # X-axis component
        Y[i] = ss16asteph_t[3i-1]-q[2]      # Y-axis component
        Z[i] = ss16asteph_t[3i  ]-q[3]      # Z-axis component

        # Distance between the i-th body and the asteroid
        r_p2[i] = ( (X[i]^2)+(Y[i]^2) ) + (Z[i]^2)  # r_{i,asteroid}^2
        r_p3d2[i] = r_p2[i]^1.5                     # r_p2^1.5 <-> r_{i, asteroid}^3

        # Newtonian coefficient, i.e., mass parameter / distance^3 -> \mu_i / r_{i, asteroid}^3
        newtonianCoeff[i] =  μ[i]/r_p3d2[i]

        # Newtonian coefficient * difference between two positions, i.e.,
        # \mu_i * (\mathbf{r_i} - \mathbf{r_asteroid}) / r_{ij}^3
        newton_acc_X[i] = X[i]*newtonianCoeff[i]
        newton_acc_Y[i] = Y[i]*newtonianCoeff[i]
        newton_acc_Z[i] = Z[i]*newtonianCoeff[i]
    end

    for i in 1:Nm1
        # Newtonian point-mass accelerations
        temp_accX_i[i] = accX + newton_acc_X[i]
        accX = temp_accX_i[i]
        temp_accY_i[i] = accY + newton_acc_Y[i]
        accY = temp_accY_i[i]
        temp_accZ_i[i] = accZ + newton_acc_Z[i]
        accZ = temp_accZ_i[i]
    end

    # Fill dq[4:6] with accelerations
    # Newtonian point-mass accelerations
    dq[4] = accX
    dq[5] = accY
    dq[6] = accZ

    nothing
end

# Number of degrees of freedom for each dynamical model
dof(::Val{newtonian!}) = 6
dof(::Val{gravityonly!}) = 6
dof(::Val{nongravs!}) = 9

# Number of jet transport variables for each dynamical model
numvars(::Val{newtonian!}, _) = 6
numvars(::Val{gravityonly!}, _) = 6
numvars(::Val{nongravs!}, params::Parameters) = 6 + count(!iszero, params.marsden_scalings)