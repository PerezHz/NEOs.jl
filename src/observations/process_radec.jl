@doc raw"""
    compute_radec(observatory::ObservatoryMPC{T}, t_r_utc::DateTime, niter::Int = 10; eo::Bool = true,
                  xve::Function = earthposvel, xvs::Function = sunposvel, xva::Function = apophisposvel197) where {T <: AbstractFloat}
    compute_radec(obs::RadecMPC{T}, niter::Int = 10; eo::Bool = true, xve::Function = earthposvel, xvs::Function = sunposvel,
                  xva::Function = apophisposvel197) where {T <: AbstractFloat}
    compute_radec(obs::Vector{RadecMPC{T}}, niter::Int=10; eo::Bool=true, xve::Function=earthposvel, xvs::Function=sunposvel,
                  xva::Function=apophisposvel197) where {T <: AbstractFloat}

Compute astrometric right ascension and declination (both in arcsec) for a set of MPC-formatted observations.

# Arguments

- `observatory::ObservatoryMPC{T}`: observation site.
- `t_r_utc::DateTime`: UTC time of astrometric observation.
- `obs::RadecMPC{T}/Vector{RadecMPC{T}}`: observations.
- `niter::Int`: number of light-time solution iterations.
- `eo::Bool`: compute corrections due to Earth orientation, LOD, polar motion.
- `xve::Function`: Earth ephemeris [et seconds since J2000] -> [barycentric position in km and velocity in km/sec].
- `xvs::Function`: Sun ephemeris [et seconds since J2000] -> [barycentric position in km and velocity in km/sec].
- `xva::Function`: asteroid ephemeris [et seconds since J2000] -> [barycentric position in km and velocity in km/sec].
"""
function compute_radec(observatory::ObservatoryMPC{T}, t_r_utc::DateTime, niter::Int = 10; eo::Bool = true,
                       xve::Function = earthposvel, xvs::Function = sunposvel, xva::Function = apophisposvel197) where {T <: AbstractFloat}
    # Transform receiving time from UTC to TDB seconds since J2000
    et_r_secs = datetime2et(t_r_utc)
    # Compute geocentric position/velocity of receiving antenna in inertial frame [km, km/s]
    RV_r = obsposvelECI(observatory, et_r_secs, eo = eo)
    R_r = RV_r[1:3]
    V_r = RV_r[4:6]
    # Earth's barycentric position and velocity at receive time
    rv_e_t_r = xve(et_r_secs)
    r_e_t_r = rv_e_t_r[1:3]
    v_e_t_r = rv_e_t_r[4:6]
    # Receiver barycentric position and velocity at receive time
    r_r_t_r = r_e_t_r + R_r
    v_r_t_r = v_e_t_r + V_r
    # Asteroid barycentric position and velocity at receive time
    rv_a_t_r = xva(et_r_secs)
    r_a_t_r = rv_a_t_r[1:3]
    # Sun barycentric position and velocity at receive time
    rv_s_t_r = xvs(et_r_secs)
    r_s_t_r = rv_s_t_r[1:3]

    # Down-leg iteration
    # τ_D first approximation
    # See equation (1) of https://doi.org/10.1086/116062
    ρ_vec_r = r_a_t_r - r_r_t_r
    ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)
    # -R_b/c, but delay is wrt asteroid Center (Brozovic et al., 2018)
    τ_D = ρ_r/clightkms # (seconds)
    # Bounce time, new estimate
    # See equation (2) of https://doi.org/10.1086/116062
    et_b_secs = et_r_secs - τ_D

    # Allocate memmory for time-delays
    Δτ_D = zero(τ_D)               # Total time delay
    Δτ_rel_D = zero(τ_D)           # Shapiro delay
    # Δτ_corona_D = zero(τ_D)      # Delay due to Solar corona
    Δτ_tropo_D = zero(τ_D)         # Delay due to Earth's troposphere

    for i in 1:niter
        # Asteroid barycentric position (in km) and velocity (in km/s) at bounce time (TDB)
        rv_a_t_b = xva(et_b_secs)
        r_a_t_b = rv_a_t_b[1:3]
        v_a_t_b = rv_a_t_b[4:6]
        # Estimated position of the asteroid's center of mass relative to the recieve point
        # See equation (3) of https://doi.org/10.1086/116062
        ρ_vec_r = r_a_t_b - r_r_t_r
        # Magnitude of ρ_vec_r
        # See equation (4) of https://doi.org/10.1086/116062
        ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)

        # Compute down-leg Shapiro delay
        # NOTE: when using PPN, substitute 2 -> 1+γ in expressions for Shapiro delay,
        # Δτ_rel_[D|U]

        # Earth's position at t_r
        e_D_vec  = r_r_t_r - r_s_t_r
        # Heliocentric distance of Earth at t_r
        e_D = sqrt(e_D_vec[1]^2 + e_D_vec[2]^2 + e_D_vec[3]^2)
        # Barycentric position of Sun at estimated bounce time
        rv_s_t_b = xvs(et_b_secs)
        r_s_t_b = rv_s_t_b[1:3]
        # Heliocentric position of asteroid at t_b
        p_D_vec  = r_a_t_b - r_s_t_b
        # Heliocentric distance of asteroid at t_b
        p_D = sqrt(p_D_vec[1]^2 + p_D_vec[2]^2 + p_D_vec[3]^2)
        # Signal path distance (down-leg)
        q_D = ρ_r

        # Shapiro correction to time-delay
        Δτ_rel_D = shapiro_delay(e_D, p_D, q_D)  # seconds
        # Troposphere correction to time-delay
        # Δτ_tropo_D = tropo_delay(R_r, ρ_vec_r) # seconds
        # Solar corona correction to time-delay
        # Δτ_corona_D = corona_delay(constant_term.(r_a_t_b), r_r_t_r, r_s_t_r, F_tx, station_code) # seconds
        # Total time-delay
        Δτ_D = Δτ_rel_D # + Δτ_tropo_D #+ Δτ_corona_D # seconds

        # New estimate
        p_dot_23 = dot(ρ_vec_r, v_a_t_b)/ρ_r
        # Time delay correction
        Δt_2 = (τ_D - ρ_r/clightkms - Δτ_rel_D)/(1.0-p_dot_23/clightkms)
        # Time delay new estimate
        τ_D = τ_D - Δt_2
        # Bounce time, new estimate
        # See equation (2) of https://doi.org/10.1086/116062
        et_b_secs = et_r_secs - τ_D
    end
    # Asteroid barycentric position (in km) and velocity (in km/s) at bounce time (TDB)
    rv_a_t_b = xva(et_b_secs)
    r_a_t_b = rv_a_t_b[1:3]
    v_a_t_b = rv_a_t_b[4:6]
    # Estimated position of the asteroid's center of mass relative to the recieve point
    # See equation (3) of https://doi.org/10.1086/116062
    ρ_vec_r = r_a_t_b - r_r_t_r
    # Magnitude of ρ_vec_r
    # See equation (4) of https://doi.org/10.1086/116062
    ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)

    # TODO: add aberration and atmospheric refraction corrections

    # Compute gravitational deflection of light
    # See Explanatory Supplement to the Astronomical Almanac (ESAA) 2014 Section 7.4.1.4
    E_H_vec = r_r_t_r -r_s_t_r    # ESAA 2014, equation (7.104)
    U_vec = ρ_vec_r               # r_a_t_b - r_e_t_r, ESAA 2014, equation (7.112)
    U_norm = ρ_r                  # sqrt(U_vec[1]^2 + U_vec[2]^2 + U_vec[3]^2)
    u_vec = U_vec/U_norm
    # Barycentric position and velocity of Sun at converged bounce time
    rv_s_t_b = xvs(et_b_secs)
    r_s_t_b = rv_s_t_b[1:3]
    Q_vec = r_a_t_b - r_s_t_b     # ESAA 2014, equation (7.113)
    q_vec = Q_vec/sqrt(Q_vec[1]^2 + Q_vec[2]^2 + Q_vec[3]^2)
    E_H = sqrt(E_H_vec[1]^2 + E_H_vec[2]^2 + E_H_vec[3]^2)
    e_vec = E_H_vec/E_H
    # See ESAA 2014, equation (7.115)
    g1 = (2μ_DE430[su]/(c_au_per_day^2))/(E_H/au)
    g2 = 1 + dot(q_vec, e_vec)
    # See ESAA 2014, equation (7.116)
    u1_vec = U_norm*(  u_vec + (g1/g2)*( dot(u_vec,q_vec)*e_vec - dot(e_vec,u_vec)*q_vec )  )
    u1_norm = sqrt(u1_vec[1]^2 + u1_vec[2]^2 + u1_vec[3]^2)

    # Compute right ascension, declination angles
    α_rad_ = mod2pi(atan(u1_vec[2], u1_vec[1]))
    α_rad = mod2pi(α_rad_)          # right ascension (rad)
    δ_rad = asin(u1_vec[3]/u1_norm) # declination (rad)

    δ_as = rad2arcsec(δ_rad) # rad -> arcsec + debiasing
    α_as = rad2arcsec(α_rad) # rad -> arcsec + debiasing

    return α_as, δ_as # right ascension, declination both in arcsec
end

function compute_radec(obs::RadecMPC{T}, niter::Int = 10; eo::Bool = true, xve::Function = earthposvel, xvs::Function = sunposvel,
                       xva::Function = apophisposvel197) where {T <: AbstractFloat}
    return compute_radec(obs.observatory, obs.date, niter; eo = eo, xve = xve, xvs = xvs, xva = xva)
end

function compute_radec(obs::Vector{RadecMPC{T}}, niter::Int = 10; eo::Bool = true, xve::Function = earthposvel, xvs::Function = sunposvel,
                       xva::Function = apophisposvel197) where {T <: AbstractFloat}

    # Number of observations
    n_optical_obs = length(obs)
    # UTC time of first astrometric observation
    utc1 = obs[1].date
    # TDB seconds since J2000.0 for first astrometric observation
    et1 = datetime2et(utc1)
    # Asteroid ephemeris at et1
    a1_et1 = xva(et1)[1]
    # Type of asteroid ephemeris
    S = typeof(a1_et1)

    # Right ascension
    vra = Array{S}(undef, n_optical_obs)
    # Declination
    vdec = Array{S}(undef, n_optical_obs)

    # Iterate over the number of observations
    for i in 1:n_optical_obs
        vra[i], vdec[i] = compute_radec(obs[i], niter, eo = eo, xve = xve, xvs = xvs, xva = xva)
    end

    return vra, vdec # arcsec, arcsec
end

# MPC catalogues corresponding to debiasing tables included in https://doi.org/10.1016/j.icarus.2014.07.033
const mpc_catalogue_codes_2014 = ["a", "b", "c", "d", "e", "g", "i", "j", "l", "m", "o", "p", "q", "r",
                                  "u", "v", "w", "L", "N"]

# MPC catalogues corresponding to debiasing tables included in https://doi.org/10.1016/j.icarus.2019.113596
const mpc_catalogue_codes_2018 = ["a", "b", "c", "d", "e", "g", "i", "j", "l", "m", "n", "o", "p", "q",
                                  "r", "t", "u", "v", "w", "L", "N", "Q", "R", "S", "U", "W"]

@doc raw"""
    select_debiasing_table(debias_table::String = "2018")

Return the catalogue codes, truth catalogue, resolution and bias matrix of the corresponding debias table. The possible values
for `debias_table` are:
- `2014` corresponds to https://doi.org/10.1016/j.icarus.2014.07.033,
- `2018` corresponds to https://doi.org/10.1016/j.icarus.2019.113596,
- `hires2018` corresponds to https://doi.org/10.1016/j.icarus.2019.113596.
"""
function select_debiasing_table(debias_table::String = "2018")
    # Debiasing tables are loaded "lazily" via Julia artifacts, according to rules in Artifacts.toml
    if debias_table == "2018"
        debias_path = artifact"debias_2018"
        mpc_catalogue_codes_201X = mpc_catalogue_codes_2018
        # The healpix tesselation resolution of the bias map from https://doi.org/10.1016/j.icarus.2019.113596
        NSIDE = 64
        # In 2018 debias table Gaia DR2 catalogue is regarded as the truth
        truth = "V"
    elseif debias_table == "hires2018"
        debias_path = artifact"debias_hires2018"
        mpc_catalogue_codes_201X = mpc_catalogue_codes_2018
        # The healpix tesselation resolution of the high-resolution bias map from https://doi.org/10.1016/j.icarus.2019.113596
        NSIDE = 256
        # In 2018 debias table Gaia DR2 catalogue is regarded as the truth
        truth = "V"
    elseif debias_table == "2014"
        debias_path = artifact"debias_2014"
        mpc_catalogue_codes_201X = mpc_catalogue_codes_2014
        # The healpix tesselation resolution of the bias map from https://doi.org/10.1016/j.icarus.2014.07.033
        NSIDE= 64
        # In 2014 debias table PPMXL catalogue is regarded as the truth
        truth = "t"
    else
        @error "Unknown bias map: $(debias_table). Possible values are `2014`, `2018` and `hires2018`."
    end

    # Debias table file
    bias_file = joinpath(debias_path, "bias.dat")
    # Read bias matrix
    bias_matrix = readdlm(bias_file, comment_char='!', comments=true)
    # Initialize healpix Resolution variable
    resol = Resolution(NSIDE)
    # Compatibility between bias matrix and resolution
    @assert size(bias_matrix) == (resol.numOfPixels, 4length(mpc_catalogue_codes_201X)) "Bias table file $bias_file dimensions do not match expected parameter NSIDE=$NSIDE and/or number of catalogs in table."

    return mpc_catalogue_codes_201X, truth, resol, bias_matrix
end

@doc raw"""
    debiasing(obs::RadecMPC{T}, mpc_catalogue_codes_201X::Vector{String}, truth::String, resol::Resolution,
              bias_matrix::Matrix{T}) where {T <: AbstractFloat}

Return total debiasing correction in right ascension and declination (both in arcsec).

# Arguments

- `obs::RadecMPC{T}`: optical observation.
- `mpc_catalogue_codes_201X::Vector{String}`: catalogues present in debiasing table.
- `truth::String`: truth catalogue of debiasing table.
- `resol::Resolution`: resolution.
- `bias_matrix::Matrix{T}`: debiasing table.
"""
function debiasing(obs::RadecMPC{T}, mpc_catalogue_codes_201X::Vector{String}, truth::String, resol::Resolution,
                   bias_matrix::Matrix{T}) where {T <: AbstractFloat}

    # Catalogue code
    catcode = obs.catalogue.code

    # If star catalogue is not present in debiasing table, then set corrections equal to zero
    if (catcode ∉ mpc_catalogue_codes_201X) && catcode != "Y"
        # Catalogue exists in mpc_catalogues[]
        if !isunknown(obs.catalogue)
            # Truth catalogue is not present in debiasing table but it does not send a warning
            if catcode != truth
                @warn "Catalogue $(obs.catalogue.name) not found in debiasing table. Setting debiasing corrections equal to zero."
            end
        # Unknown catalogue
        elseif catcode == ""
            @warn "Catalog information not available in observation record. Setting debiasing corrections equal to zero."
        # Catalogue code is not empty but it does not match an MPC catalogue code either
        else
            @warn "Catalog code $catcode does not correspond to MPC catalogue code. Setting debiasing corrections equal to zero."
        end
        α_corr = zero(T)
        δ_corr = zero(T)
    # If star catalogue is present in debiasing table, then compute corrections
    else
        # Get pixel tile index, assuming iso-latitude rings indexing, which is the formatting in `tiles.dat`.
        # Substracting 1 from the returned value of `ang2pixRing` corresponds to 0-based indexing, as in `tiles.dat`;
        # not substracting 1 from the returned value of `ang2pixRing` corresponds to 1-based indexing, as in Julia.
        # Since we use pix_ind to get the corresponding row number in `bias.dat`, it's not necessary to substract 1.
        pix_ind = ang2pixRing(resol, π/2 - obs.δ, obs.α)

        # Handle edge case: in new MPC catalogue nomenclature, "UCAC-5"->"Y"; but in debias tables "UCAC-5"->"W"
        if catcode == "Y"
            cat_ind = findfirst(x -> x == "W", mpc_catalogue_codes_201X)
        else
            cat_ind = findfirst(x -> x == catcode, mpc_catalogue_codes_201X)
        end

        # Read dRA, pmRA, dDEC, pmDEC data from bias.dat
        # dRA: position correction in RA * cos(DEC) at epoch J2000.0 [arcsec]
        # dDEC: position correction in DEC at epoch J2000.0 [arcsec]
        # pmRA: proper motion correction in RA*cos(DEC) [mas/yr]
        # pmDEC: proper motion correction in DEC [mas/yr]
        dRA, dDEC, pmRA, pmDEC = bias_matrix[pix_ind, 4*cat_ind-3:4*cat_ind]
        # Seconds since J2000 (TDB)
        et_secs_i = datetime2et(obs.date)
        # Seconds sinde J2000 (TT)
        tt_secs_i = et_secs_i - ttmtdb(et_secs_i)
        # Years since J2000
        yrs_J2000_tt = tt_secs_i/(daysec*yr)
        # Total debiasing correction in right ascension (arcsec)
        α_corr = dRA + yrs_J2000_tt*pmRA/1_000
        # Total debiasing correction in declination (arcsec)
        δ_corr = dDEC + yrs_J2000_tt*pmDEC/1_000
    end

    return α_corr, δ_corr
end

@doc raw"""
    w8sveres17(obs::RadecMPC{T}) where {T <: AbstractFloat}

Return the statistical weight from Veres et al. (2017) corresponding to `obs`.
"""
function w8sveres17(obs::RadecMPC{T}) where {T <: AbstractFloat}

    obscode = obs.observatory.code
    dt_utc_obs = obs.date
    catalogue = obs.catalogue.code

    # Unit weight (arcseconds)
    w = 1.0
    # Table 2: epoch-dependent astrometric residuals
    if obscode == "703"
        return Date(dt_utc_obs) < Date(2014,1,1) ? w : 0.8w
    elseif obscode == "691"
        return Date(dt_utc_obs) < Date(2003,1,1) ? 0.6w : 0.5w
    elseif obscode == "644"
        return Date(dt_utc_obs) < Date(2003,9,1) ? 0.6w : 0.4w
    # Table 3: most active CCD asteroid observers
    elseif obscode ∈ ("704", "C51", "J75")
        return w
    elseif obscode == "G96"
        return 0.5w
    elseif obscode == "F51"
        return 0.2w
    elseif obscode ∈ ("G45", "608")
        return 0.6w
    elseif obscode == "699"
        return 0.8w
    elseif obscode ∈ ("D29", "E12")
        return 0.75w
    # Table 4:
    elseif obscode ∈ ("645", "673", "H01")
        return 0.3w
    elseif obscode ∈ ("J04", "K92", "K93", "Q63", "Q64", "V37", "W85", "W86", "W87", "K91", "E10", "F65") #Tenerife + Las Cumbres
        return 0.4w
    elseif obscode ∈ ("689", "950", "W84")
        return 0.5w
    #elseif obscode ∈ ("G83", "309") # Applies only to program code assigned to M. Micheli
    #    if catalogue ∈ ("q", "t") # "q"=>"UCAC-4", "t"=>"PPMXL"
    #        return 0.3w
    #    elseif catalogue ∈ ("U", "V") # Gaia-DR1, Gaia-DR2
    #        return 0.2w
    #    end
    elseif obscode ∈ ("Y28",)
        if catalogue ∈ ("t", "U", "V")
            return 0.3w
        else
            return w
        end
    elseif obscode ∈ ("568",)
        if catalogue ∈ ("o", "s") # "o"=>"USNO-B1.0", "s"=>"USNO-B2.0"
            return 0.5w
        elseif catalogue ∈ ("U", "V") # Gaia DR1, DR2
            return 0.1w
        elseif catalogue ∈ ("t",) #"t"=>"PPMXL"
            return 0.2w
        else
            return w
        end
    elseif obscode ∈ ("T09", "T12", "T14") && catalogue ∈ ("U", "V") # Gaia DR1, DR2
        return 0.1w
    elseif catalogue == ""
        return 1.5w
    elseif catalogue != ""
        return w
    else
        return w
    end
end

@doc raw"""
    radec_astrometry(obs::Vector{RadecMPC{T}}, niter::Int = 10; eo::Bool = true, debias_table::String = "2018",
                     xve::Function = earthposvel, xvs::Function = sunposvel, xva::Function = apophisposvel197) where {T <: AbstractFloat}
    radec_astrometry(outfilename::String, opticalobsfile::String, asteph::TaylorInterpolant, ss16asteph::TaylorInterpolant,
                     niter::Int = 5; eo::Bool = true, debias_table::String = "2018")

Compute optical astrometry, i.e. dates of observation plus observed, computed, debiasing corrections and statistical weights for
right ascension and declination (in arcsec).

See also [`compute_radec`](@ref), [`debiasing`](@ref), [`w8sveres17`](@ref) and [`Healpix.ang2pixRing`](@ref).

# Arguments

- `obs::Vector{RadecMPC{T}}`: vector of observations.
- `outfilename::String`: file where to save the results (.jld2).
- `opticalobsfile::String`: file where to retrieve optical observations.
- `asteph::TaylorInterpolant`: NEO's ephemeris.
- `ss16asteph::TaylorInterpolant`: solar system ephemeris.
- `niter::Int`: number of light-time solution iterations.
- `eo::Bool`: compute corrections due to Earth orientation, LOD, polar motion.
- `debias_table::String`: possible values are:
    - `2014` corresponds to https://doi.org/10.1016/j.icarus.2014.07.033,
    - `2018` corresponds to https://doi.org/10.1016/j.icarus.2019.113596,
    - `hires2018` corresponds to https://doi.org/10.1016/j.icarus.2019.113596.
- `xve::Function`: Earth ephemeris [et seconds since J2000] -> [barycentric position in km and velocity in km/sec].
- `xvs::Function`: Sun ephemeris [et seconds since J2000] -> [barycentric position in km and velocity in km/sec].
- `xva::Function`: asteroid ephemeris [et seconds since J2000] -> [barycentric position in km and velocity in km/sec].
"""
function radec_astrometry(obs::Vector{RadecMPC{T}}, niter::Int = 10; eo::Bool = true, debias_table::String = "2018",
                          xve::Function = earthposvel, xvs::Function = sunposvel, xva::Function = apophisposvel197) where {T <: AbstractFloat}

    # Number of observations
    n_optical_obs = length(obs)

    # UTC time of first astrometric observation
    utc1 = obs[1].date
    # TDB seconds since J2000.0 for first astrometric observation
    et1 = datetime2et(utc1)
    # Asteroid ephemeris at et1
    a1_et1 = xva(et1)[1]
    # Type of asteroid ephemeris
    S = typeof(a1_et1)

    # Observed right ascension
    α_obs = Array{Float64}(undef, n_optical_obs)
    # Observed declination
    δ_obs = Array{Float64}(undef, n_optical_obs)

    # Computed right ascension
    α_comp = Array{S}(undef, n_optical_obs)
    # Computed declination
    δ_comp = Array{S}(undef, n_optical_obs)

    # Total debiasing correction in right ascension (arcsec)
    α_corr = Array{Float64}(undef, n_optical_obs)
    # Total debiasing correction in declination (arcsec)
    δ_corr = Array{Float64}(undef, n_optical_obs)

    # Times of observations
    datetime_obs = Array{DateTime}(undef, n_optical_obs)
    # Statistical weights
    w8s = Array{Float64}(undef, n_optical_obs)

    # Select debias table
    mpc_catalogue_codes_201X, truth, resol, bias_matrix = select_debiasing_table(debias_table)

    # Iterate over the observations
    for i in 1:n_optical_obs

        # Time of observation
        datetime_obs[i] = obs[i].date

        # Observed ra/dec
        # Note: ra is multiplied by a metric factor cos(dec) to match the format of debiasing corrections
        δ_obs[i] = rad2arcsec(obs[i].δ)                 # arcsec
        α_obs[i] = rad2arcsec(obs[i].α) * cos(obs[i].δ) # arcsec

        # Computed ra/dec
        α_comp_as, δ_comp_as = compute_radec(obs[i], niter, eo=eo, xve=xve, xvs=xvs, xva=xva)
        # Multiply by metric factor cos(dec)
        α_comp[i] = α_comp_as * cos(obs[i].δ)  # arcsec
        δ_comp[i] = δ_comp_as                  # arcsec

        # Debiasing corrections
        α_corr[i], δ_corr[i] = debiasing(obs[i], mpc_catalogue_codes_201X, truth, resol, bias_matrix)

        # Statistical weights from Veres et al. (2017)
        w8s[i] = w8sveres17(obs[i])

    end

    # Return time of observation, observed ra/dec, computed ra/dec, total debiasing correction in ra/dec and statistical weights
    return datetime_obs, α_obs, δ_obs, α_comp, δ_comp, α_corr, δ_corr, w8s
end

function radec_astrometry(outfilename::String, opticalobsfile::String, asteph::TaylorInterpolant, ss16asteph::TaylorInterpolant,
                          niter::Int = 5; eo::Bool = true, debias_table::String = "2018")

    # Read optical observations
    radec = read_radec_mpc(opticalobsfile)

    # Check that first and last observation times are within interpolation interval
    Δt_0 = datetime2julian(radec[1].date) - JD_J2000 - asteph.t0
    Δt_f = datetime2julian(radec[end].date) - JD_J2000 - asteph.t0
    t_min, t_max = minmax(asteph.t[1], asteph.t[end])
    @assert t_min ≤ Δt_0 ≤ Δt_f ≤ t_max "First and/or last observation times are outside interpolation interval"

    # Number of massive bodies
    Nm1 = (size(ss16asteph.x)[2]-13) ÷ 6
    # Number of bodies, including NEA
    N = Nm1 + 1

    # Change t, x, v units, resp., from days, au, au/day to sec, km, km/sec

    # NEO
    function asteph_et(et)
        return auday2kmsec(asteph(et/daysec)[1:6])
    end
    # Earth
    function earth_et(et)
        return auday2kmsec(ss16asteph(et)[union(3*4-2:3*4,3*(N-1+4)-2:3*(N-1+4))])
    end
    # Sun
    function sun_et(et)
        return auday2kmsec(ss16asteph(et)[union(3*1-2:3*1,3*(N-1+1)-2:3*(N-1+1))])
    end

    # Compute ra/dec astrometry
    datetime_obs, α_obs, δ_obs, α_comp, δ_comp, α_corr, δ_corr, w8s = radec_astrometry(radec, niter; eo = eo, debias_table = debias_table,
                                                                                       xve = earth_et, xvs = sun_et, xva = asteph_et)

    # Save data to file
    println("Saving data to file: $outfilename")
    JLD2.jldopen(outfilename, "w") do file
        # Write variables to jld2 file
        JLD2.write(file,
            "datetime_obs", datetime_obs,
            "α_obs", α_obs,
            "δ_obs", δ_obs,
            "α_comp", α_comp,
            "δ_comp", δ_comp,
            "α_corr", α_corr,
            "δ_corr", δ_corr,
            "w8s", w8s
        )
    end

    return nothing
end

@doc raw"""
    residuals(obs::Vector{RadecMPC{T}}, niter::Int = 10; eo::Bool = true, debias_table::String = "2018",
              xvs::Function = sunposvel, xve::Function = earthposvel, xva::Function = apophisposvel197) where {T <: AbstractFloat}

Compute optical O-C residuals.

See also [`compute_radec`](@ref), [`debiasing`](@ref), [`w8sveres17`](@ref) and [`radec_astrometry`](@ref).

# Arguments

- `obs::Vector{RadecMPC{T}}`: vector of observations.
- `niter::Int`: number of light-time solution iterations.
- `eo::Bool`: compute corrections due to Earth orientation, LOD, polar motion.
- `debias_table::String`: possible values are:
    - `2014` corresponds to https://doi.org/10.1016/j.icarus.2014.07.033,
    - `2018` corresponds to https://doi.org/10.1016/j.icarus.2019.113596,
    - `hires2018` corresponds to https://doi.org/10.1016/j.icarus.2019.113596.
- `xvs::Function`: Sun ephemeris [et seconds since J2000] -> [barycentric position in km and velocity in km/sec].
- `xve::Function`: Earth ephemeris [et seconds since J2000] -> [barycentric position in km and velocity in km/sec].
- `xva::Function`: asteroid ephemeris [et seconds since J2000] -> [barycentric position in km and velocity in km/sec].
"""
function residuals(obs::Vector{RadecMPC{T}}, niter::Int = 10; eo::Bool = true, debias_table::String = "2018",
                   xvs::Function = sunposvel, xve::Function = earthposvel, xva::Function = apophisposvel197) where {T <: AbstractFloat}

    # Optical astrometry (dates + observed + computed + debiasing + weights)
    x_jt = radec_astrometry(obs, niter; eo = eo, debias_table = debias_table, xvs = xvs, xve = xve, xva = xva)
    # Right ascension residuals
    res_α = x_jt[2] .- x_jt[6] .- x_jt[4]
    # Declination residuals
    res_δ = x_jt[3] .- x_jt[7] .- x_jt[5]
    # Total residuals
    res = vcat(res_α, res_δ)
    # Weights
    w = repeat(1 ./ x_jt[end].^2, 2)

    return res, w
end