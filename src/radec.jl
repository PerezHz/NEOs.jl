# Compute astrometric right ascension and declination for an asteroid at
# UTC instant `t_r_utc` from tracking station with code `station_code` from Earth,
# Sun and asteroid ephemerides
# station_code: observing station identifier (MPC nomenclature)
# t_r_utc: UTC time of astrometric observation (DateTime)
# niter: number of light-time solution iterations
# pm: compute polar motion corrections
# lod: compute corrections due to variations in length of day
# eocorr: compute corrections due to Earth orientation parameters
# debias: compute debiasing according to Eggl et al. (2019)
# catalog: Stellar catalog used in astrometric reduction, MPC nomenclature ("a", "b", "c", etc...)
# xve: Earth ephemeris wich takes et seconds since J2000 as input and returns Earth barycentric position in km and velocity in km/second
# xvs: Sun ephemeris wich takes et seconds since J2000 as input and returns Sun barycentric position in km and velocity in km/second
# xva: asteroid ephemeris wich takes et seconds since J2000 as input and returns asteroid barycentric position in km and velocity in km/second
function radec(station_code::Union{Int,String}, t_r_utc::DateTime,
        niter::Int=10; pm::Bool=true, lod::Bool=true, eocorr::Bool=true,
        xve::Function=earth_pv, xvs::Function=sun_pv, xva::Function=apophis_pv_197)
    et_r_secs = str2et(string(t_r_utc))
    # Compute geocentric position/velocity of receiving antenna in inertial frame (au, au/day)
    R_r, V_r = observer_position(station_code, et_r_secs, pm=pm, lod=lod, eocorr=eocorr)
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
    # down-leg iteration
    # τ_D first approximation: Eq. (1) Yeomans et al. (1992)
    ρ_vec_r = r_a_t_r - r_r_t_r
    ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)
    τ_D = ρ_r/clightkms # (seconds) -R_b/c, but delay is wrt asteroid Center (Brozovic et al., 2018)
    # bounce time, new estimate Eq. (2) Yeomans et al. (1992)
    et_b_secs = et_r_secs - τ_D

    Δτ_D = zero(τ_D)
    Δτ_rel_D = zero(τ_D)
    # Δτ_corona_D = zero(τ_D)
    Δτ_tropo_D = zero(τ_D)

    for i in 1:niter
        # asteroid barycentric position (in au) at bounce time (TDB)
        rv_a_t_b = xva(et_b_secs)
        r_a_t_b = rv_a_t_b[1:3]
        v_a_t_b = rv_a_t_b[4:6]
        # Eq. (3) Yeomans et al. (1992)
        ρ_vec_r = r_a_t_b - r_r_t_r
        # Eq. (4) Yeomans et al. (1992)
        ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)
        # compute down-leg Shapiro delay
        # NOTE: when using PPN, substitute 2 -> 1+γ in expressions for Shapiro delay, Δτ_rel_[D|U]
        e_D_vec  = r_r_t_r - r_s_t_r
        e_D = sqrt(e_D_vec[1]^2 + e_D_vec[2]^2 + e_D_vec[3]^2) # heliocentric distance of Earth at t_r
        rv_s_t_b = xvs(et_b_secs) # barycentric position and velocity of Sun at estimated bounce time
        r_s_t_b = rv_s_t_b[1:3]
        p_D_vec  = r_a_t_b - r_s_t_b
        p_D = sqrt(p_D_vec[1]^2 + p_D_vec[2]^2 + p_D_vec[3]^2) # heliocentric distance of asteroid at t_b
        q_D = ρ_r #signal path distance (down-leg)
        # Shapiro correction to time-delay
        Δτ_rel_D = shapiro_delay(e_D, p_D, q_D)
        # # troposphere correction to time-delay
        # Δτ_tropo_D = tropo_delay(R_r, ρ_vec_r) # seconds
        # Δτ_corona_D = corona_delay(constant_term.(r_a_t_b), r_r_t_r, r_s_t_r, F_tx, station_code) # seconds
        Δτ_D = Δτ_rel_D # + Δτ_tropo_D #+ Δτ_corona_D # seconds
        p_dot_23 = dot(ρ_vec_r, v_a_t_b)/ρ_r
        Δt_2 = (τ_D - ρ_r/clightkms - Δτ_rel_D)/(1.0-p_dot_23/clightkms)
        τ_D = τ_D - Δt_2
        et_b_secs = et_r_secs - τ_D
    end
    rv_a_t_b = xva(et_b_secs)
    r_a_t_b = rv_a_t_b[1:3]
    v_a_t_b = rv_a_t_b[4:6]

    ρ_vec_r = r_a_t_b - r_r_t_r
    ρ_r = sqrt(ρ_vec_r[1]^2 + ρ_vec_r[2]^2 + ρ_vec_r[3]^2)

    # TODO: add aberration and atmospheric refraction corrections

    # Compute gravitational deflection of light, ESAA 2014 Section 7.4.1.4
    E_H_vec = r_r_t_r -r_s_t_r # ESAA 2014, Eq. (7.104)
    U_vec = ρ_vec_r #r_a_t_b - r_e_t_r # ESAA 2014, Eq. (7.112)
    U_norm = ρ_r # sqrt(U_vec[1]^2 + U_vec[2]^2 + U_vec[3]^2)
    u_vec = U_vec/U_norm
    rv_s_t_b = xvs(et_b_secs) # barycentric position and velocity of Sun at converged bounce time
    r_s_t_b = rv_s_t_b[1:3]
    Q_vec = r_a_t_b - r_s_t_b # ESAA 2014, Eq. (7.113)
    q_vec = Q_vec/sqrt(Q_vec[1]^2 + Q_vec[2]^2 + Q_vec[3]^2)
    E_H = sqrt(E_H_vec[1]^2 + E_H_vec[2]^2 + E_H_vec[3]^2)
    e_vec = E_H_vec/E_H
    g1 = (2μ[1]/(c_au_per_day^2))/(E_H/au) # ESAA 2014, Eq. (7.115)
    g2 = 1 + dot(q_vec, e_vec)
    # @show g1, g2
    u1_vec = U_norm*(  u_vec + (g1/g2)*( dot(u_vec,q_vec)*e_vec - dot(e_vec,u_vec)*q_vec )  ) # ESAA 2014, Eq. (7.116)

    # Compute aberration of light, ESAA 2014 Section 7.4.1.5
    u1_norm = sqrt(u1_vec[1]^2 + u1_vec[2]^2 + u1_vec[3]^2)
    # @show norm(u1_vec/U_norm), μ[1]/(c_au_per_day^2) # u1_vec/U_norm is a unit vector to order μ[1]/(c_au_per_day^2)
    u_vec_new = u1_vec/u1_norm
    V_vec = v_r_t_r/clightkms
    V_norm = sqrt(V_vec[1]^2 + V_vec[2]^2 + V_vec[3]^2)
    β_m1 = sqrt(1-V_norm^2) # β^-1
    # @show norm(v_r_t_r), V_norm, β_m1
    # @show β_m1
    f1 = dot(u_vec_new, V_vec)
    f2 = 1 + f1/(1+β_m1)
    # u2_vec = ρ_vec_r; u2_norm = ρ_r # uncorrected radial unit vector
    # u2_vec = u1_vec; u2_norm = sqrt(u2_vec[1]^2 + u2_vec[2]^2 + u2_vec[3]^2) # radial unit vector with grav deflection corr.
    # u2_vec = ( β_m1*u1_vec + f2*u1_norm*V_vec )/( 1+f1 ) # ESAA 2014, Eq. (7.118)
    # u2_norm = sqrt(u2_vec[1]^2 + u2_vec[2]^2 + u2_vec[3]^2)
    # u2_vec = u1_vec #/U_norm + V_vec # ESAA 2014, Eq. (7.118)
    # u2_vec = u1_vec + u1_norm*V_vec # ESAA 2014, Eq. (7.119)
    u2_vec = u1_vec
    # u2_vec = u_vec
    u2_norm = sqrt(u2_vec[1]^2 + u2_vec[2]^2 + u2_vec[3]^2)

    # Compute right ascension, declination angles
    α_rad_ = mod2pi(atan(u2_vec[2], u2_vec[1]))
    α_rad = mod2pi(α_rad_) # right ascension (rad)
    δ_rad = asin(u2_vec[3]/u2_norm) # declination (rad)

    δ_as = rad2arcsec(δ_rad) # rad -> arcsec + debiasing
    α_as = rad2arcsec(α_rad) # rad -> arcsec + debiasing

    return α_as, δ_as # right ascension, declination (arcsec, arcsec)
end

# Compute optical astrometric ra/dec ephemeris for a set of observations in a MPC-formatted file
function radec_mpc(astopticalobsfile, niter::Int=10; pm::Bool=true, lod::Bool=true,
        eocorr::Bool=true, xve::Function=earth_pv,
        xvs::Function=sun_pv, xva::Function=apophis_pv_197)

    obs_df = readmp(astopticalobsfile)
    n_optical_obs = nrow(obs_df)

    utc1 = DateTime(obs_df.yr[1], obs_df.month[1], obs_df.day[1]) + Microsecond( round(1e6*86400*obs_df.utc[1]) )
    et1 = str2et(string(utc1))
    a1_et1 = xva(et1)[1]
    S = typeof(a1_et1)

    vra = Array{S}(undef, n_optical_obs)
    vdec = Array{S}(undef, n_optical_obs)

    for i in 1:n_optical_obs
        utc_i = DateTime(obs_df.yr[i], obs_df.month[i], obs_df.day[i]) + Microsecond( round(1e6*86400*obs_df.utc[i]) )
        station_code_i = string(obs_df.obscode[i])
        vra[i], vdec[i] = radec(station_code_i, utc_i, niter, pm=pm, lod=lod,
            eocorr=eocorr, xve=xve, xvs=xvs, xva=xva)
    end

    return vra, vdec # arcsec, arcsec
end

function radec_table(mpcobsfile::String, niter::Int=10; pm::Bool=true,
        lod::Bool=true, eocorr::Bool=true, debias_table::String="2018",
        xve::Function=earth_pv, xvs::Function=sun_pv,
        xva::Function=apophis_pv_197)

    obs_t = table(readmp(mpcobsfile))
    n_optical_obs = length(obs_t)

    utc1 = DateTime(obs_t[1].yr, obs_t[1].month, obs_t[1].day) + Microsecond( round(1e6*86400*obs_t[1].utc) )
    et1 = str2et(string(utc1))
    a1_et1 = xva(et1)[1]
    S = typeof(a1_et1)

    α_obs = Array{Float64}(undef, n_optical_obs)
    δ_obs = Array{Float64}(undef, n_optical_obs)
    α_comp = Array{S}(undef, n_optical_obs)
    δ_comp = Array{S}(undef, n_optical_obs)
    α_corr = Array{Float64}(undef, n_optical_obs)
    δ_corr = Array{Float64}(undef, n_optical_obs)
    datetime_obs = Array{DateTime}(undef, n_optical_obs)

    # Select debiasing table: 2014 corresponds to Farnocchia et al. (2015), 2018 corresponds to Eggl et al. (2020)
    # debiasing tables are loaded "lazily" via Julia artifacts, according to rules in Artifacts.toml
    if debias_table == "2018"
        debias_path = artifact"debias_2018"
        mpc_catalog_codes_201X = mpc_catalog_codes_2018
        NSIDE= 64 #The healpix tesselation resolution of the bias map from Eggl et al. (2020)
        truth = "V" # In 2018 debias table Gaia DR2 catalog is regarded as the truth
    elseif debias_table == "hires2018"
        debias_path = artifact"debias_hires2018"
        mpc_catalog_codes_201X = mpc_catalog_codes_2018
        NSIDE= 256 #The healpix tesselation resolution of the high-resolution bias map from Eggl et al. (2020)
        truth = "V" # In 2018 debias table Gaia DR2 catalog is regarded as the truth
    elseif debias_table == "2014"
        debias_path = artifact"debias_2014"
        mpc_catalog_codes_201X = mpc_catalog_codes_2014
        NSIDE= 64 #The healpix tesselation resolution of the bias map from Farnocchia et al. (2015)
        truth = "t" # In 2014 debias table PPMXL catalog is regarded as the truth
    else
        @error "Unknown bias map: $(debias_table). Possible values are `2014`, `2018` and `hires2018`."
    end
    bias_file = joinpath(debias_path, "bias.dat")

    bias_matrix = readdlm(bias_file, comment_char='!', comments=true)
    resol = Resolution(NSIDE) # initialize healpix Resolution variable
    @assert size(bias_matrix) == (resol.numOfPixels, 4length(mpc_catalog_codes_201X)) "Bias table file $bias_file dimensions do not match expected parameter NSIDE=$NSIDE and/or number of catalogs in table."

    for i in 1:n_optical_obs
        # observed values
        δ_i_deg = sign(obs_t[i].decd)*(abs(obs_t[i].decd) + obs_t[i].decm/60 + obs_t[i].decs/3600) # deg
        # the following if ... block handles the sign of declination, including edge cases in declination such as -00 01
        if obs_t[i].signdec == "+"
            δ_i_deg = +(obs_t[i].decd + obs_t[i].decm/60 + obs_t[i].decs/3600) # deg
        elseif obs_t[i].signdec == "-"
            δ_i_deg = -(obs_t[i].decd + obs_t[i].decm/60 + obs_t[i].decs/3600) # deg
        else
            @warn "Could not parse declination sign: $(obs_t[i].signdec). Setting positive sign."
            δ_i_deg =  (obs_t[i].decd + obs_t[i].decm/60 + obs_t[i].decs/3600) # deg
        end
        δ_i_rad = deg2rad(δ_i_deg) #rad
        α_i_deg = 15(obs_t[i].rah + obs_t[i].ram/60 + obs_t[i].ras/3600) # deg
        α_i_rad = deg2rad(α_i_deg) #rad
        δ_obs[i] = 3600δ_i_deg # arcsec
        α_obs[i] = 3600α_i_deg*cos(δ_i_rad) # arcsec
        # computed values
        datetime_obs[i] = DateTime(obs_t[i].yr, obs_t[i].month, obs_t[i].day) + Microsecond( round(1e6*86400*obs_t[i].utc) )
        # utc_i = datetime_obs[i]
        station_code_i = string(obs_t[i].obscode)
        α_comp_as, δ_comp_as = radec(station_code_i, datetime_obs[i], niter, pm=pm, lod=lod,
            eocorr=eocorr, xve=xve, xvs=xvs, xva=xva)
        α_comp[i] = α_comp_as*cos(δ_i_rad) # arcsec
        δ_comp[i] = δ_comp_as # arcsec
        if obs_t[i].catalog ∉ mpc_catalog_codes_201X
            # Handle case: if star catalog not present in debiasing table, then set corrections equal to zero
            if haskey(mpc_catalog_codes, obs_t[i].catalog)
                if obs_t[i].catalog != truth
                    catalog_not_found = mpc_catalog_codes[obs_t[i].catalog]
                    @warn "Catalog not found in $(debias_table) table: $(catalog_not_found). Setting debiasing corrections equal to zero."
                end
            elseif obs_t[i].catalog == " "
                @warn "Catalog information not available in observation record. Setting debiasing corrections equal to zero."
            else
                @warn "Catalog code $(obs_t[i].catalog) does not correspond to MPC catalog code. Setting debiasing corrections equal to zero."
            end
            α_corr[i] = 0.0
            δ_corr[i] = 0.0
            continue
        else
            # Otherwise, if star catalog is present in debias table, compute corrections
            # get pixel tile index, assuming iso-latitude rings indexing, which is the formatting in tiles.dat
            # substracting 1 from the returned value of `ang2pixRing` corresponds to 0-based indexing, as in tiles.dat
            # not substracting 1 from the returned value of `ang2pixRing` corresponds to 1-based indexing, as in Julia
            # since we use pix_ind to get the corresponding row number in bias.dat, it's not necessary to substract 1
            pix_ind = ang2pixRing(resol, π/2-δ_i_rad, α_i_rad)
            cat_ind = findfirst(x->x==obs_t[i].catalog, mpc_catalog_codes_201X)
            # @show pix_ind, cat_ind
            # read dRA, pmRA, dDEC, pmDEC data from bias.dat
            # dRA, position correction in RA*cos(DEC) at epoch J2000.0 [arcsec];
            # dDEC, position correction in DEC at epoch J2000.0 [arcsec];
            # pmRA, proper motion correction in RA*cos(DEC) [mas/yr];
            # pmDEC, proper motion correction in DEC [mas/yr].
            dRA, dDEC, pmRA, pmDEC = bias_matrix[pix_ind, 4*cat_ind-3:4*cat_ind]
            # @show dRA, dDEC, pmRA, pmDEC
            # utc_i = DateTime(obs_t[i].yr, obs_t[i].month, obs_t[i].day) + Microsecond( round(1e6*86400*obs_t[i].utc) )
            et_secs_i = str2et(string(datetime_obs[i]))
            tt_secs_i = et_secs_i - ttmtdb(et_secs_i)
            yrs_J2000_tt = tt_secs_i/(daysec*yr)
            α_corr[i] = dRA + yrs_J2000_tt*pmRA/1000 # total debiasing correction in right ascension (arcsec)
            δ_corr[i] = dDEC + yrs_J2000_tt*pmDEC/1000 # total debiasing correction in declination (arcsec)
        end
    end

    return insertcolsafter(obs_t, length(columns(obs_t)), :dt_utc_obs => datetime_obs, :α_obs => α_obs, :δ_obs => δ_obs, :α_comp => α_comp, :δ_comp => δ_comp, :α_corr => α_corr, :δ_corr => δ_corr)
end

# Compute ra/dec astrometry debiasing corrections for MPC-formatted data in
# `mpcobsfile`. Star catalog info is extracted from column 72 of each
# observation row. Current debiasing table options are:
# "2014" => Farnocchia et al. (2015)
# "2018" => Eggl et al. (2020)
# "hires2018" => Eggl et al. (2020), high-resolution (NSIDE=256)
# default `table` value is "2018".
function radec_mpc_corr(mpcobsfile::String, debias_table::String="2018")

    obs_df = readmp(mpcobsfile)
    n_optical_obs = nrow(obs_df)

    α_corr_v = Array{Float64}(undef, n_optical_obs)
    δ_corr_v = Array{Float64}(undef, n_optical_obs)

    # Select debiasing table: 2014 corresponds to Farnocchia et al. (2015), 2018 corresponds to Eggl et al. (2020)
    # debiasing tables are loaded "lazily" via Julia artifacts, according to rules in Artifacts.toml
    if debias_table == "2018"
        debias_path = artifact"debias_2018"
        mpc_catalog_codes_201X = mpc_catalog_codes_2018
        NSIDE= 64 #The healpix tesselation resolution of the bias map from Eggl et al. (2020)
        truth = "V" # In 2018 debias table Gaia DR2 catalog is regarded as the truth
    elseif debias_table == "hires2018"
        debias_path = artifact"debias_hires2018"
        mpc_catalog_codes_201X = mpc_catalog_codes_2018
        NSIDE= 256 #The healpix tesselation resolution of the high-resolution bias map from Eggl et al. (2020)
        truth = "V" # In 2018 debias table Gaia DR2 catalog is regarded as the truth
    elseif debias_table == "2014"
        debias_path = artifact"debias_2014"
        mpc_catalog_codes_201X = mpc_catalog_codes_2014
        NSIDE= 64 #The healpix tesselation resolution of the bias map from Farnocchia et al. (2015)
        truth = "t" # In 2014 debias table PPMXL catalog is regarded as the truth
    else
        @error "Unknown bias map: $(debias_table). Possible values are `2014`, `2018` and `hires2018`."
    end
    bias_file = joinpath(debias_path, "bias.dat")

    bias_matrix = readdlm(bias_file, comment_char='!', comments=true)
    resol = Resolution(NSIDE) # initialize healpix Resolution variable
    @assert size(bias_matrix) == (resol.numOfPixels, 4length(mpc_catalog_codes_201X)) "Bias table file $bias_file dimensions do not match expected parameter NSIDE=$NSIDE and/or number of catalogs in table."

    for i in 1:n_optical_obs
        if obs_df.catalog[i] ∉ mpc_catalog_codes_201X
            # Handle case: if star catalog not present in debiasing table, then set corrections equal to zero
            if haskey(mpc_catalog_codes, obs_df.catalog[i]) && obs_df.catalog[i] != truth
                if obs_df.catalog[i] != truth
                    catalog_not_found = mpc_catalog_codes[obs_df.catalog[i]]
                    @warn "Catalog not found in $(debias_table) table: $(catalog_not_found). Setting debiasing corrections equal to zero."
                end
            elseif obs_df.catalog[i] == " "
                @warn "Catalog information not available in observation record. Setting debiasing corrections equal to zero."
            else
                @warn "Catalog code $(obs_df.catalog[i]) does not correspond to MPC catalog code. Setting debiasing corrections equal to zero."
            end
            α_corr_v[i] = 0.0
            δ_corr_v[i] = 0.0
            continue
        else
            # Otherwise, if star catalog is present in debias table, compute corrections
            α_i_deg = 15(obs_df.rah[i] + obs_df.ram[i]/60 + obs_df.ras[i]/3600) # deg
            # the following if ... block handles the sign of declination, including edge cases in declination such as -00 01
            if obs_df.signdec[i] == "+"
                δ_i_deg = +(obs_df.decd[i] + obs_df.decm[i]/60 + obs_df.decs[i]/3600) # deg
            elseif obs_df.signdec[i] == "-"
                δ_i_deg = -(obs_df.decd[i] + obs_df.decm[i]/60 + obs_df.decs[i]/3600) # deg
            else
                @warn "Could not parse declination sign: $(obs_df.signdec[i]). Setting positive sign."
                δ_i_deg =  (obs_df.decd[i] + obs_df.decm[i]/60 + obs_df.decs[i]/3600) # deg
            end
            α_i_rad = deg2rad(α_i_deg) #rad
            δ_i_rad = deg2rad(δ_i_deg) #rad
            # get pixel tile index, assuming iso-latitude rings indexing, which is the formatting in tiles.dat
            # substracting 1 from the returned value of `ang2pixRing` corresponds to 0-based indexing, as in tiles.dat
            # not substracting 1 from the returned value of `ang2pixRing` corresponds to 1-based indexing, as in Julia
            # since we use pix_ind to get the corresponding row number in bias.dat, it's not necessary to substract 1
            pix_ind = ang2pixRing(resol, π/2-δ_i_rad, α_i_rad)
            cat_ind = findfirst(x->x==obs_df.catalog[i], mpc_catalog_codes_201X)
            # @show pix_ind, cat_ind
            # read dRA, pmRA, dDEC, pmDEC data from bias.dat
            # dRA, position correction in RA*cos(DEC) at epoch J2000.0 [arcsec];
            # dDEC, position correction in DEC at epoch J2000.0 [arcsec];
            # pmRA, proper motion correction in RA*cos(DEC) [mas/yr];
            # pmDEC, proper motion correction in DEC [mas/yr].
            dRA, dDEC, pmRA, pmDEC = bias_matrix[pix_ind, 4*cat_ind-3:4*cat_ind]
            # @show dRA, dDEC, pmRA, pmDEC
            utc_i = DateTime(obs_df.yr[i], obs_df.month[i], obs_df.day[i]) + Microsecond( round(1e6*86400*obs_df.utc[i]) )
            et_secs_i = str2et(string(utc_i))
            tt_secs_i = et_secs_i - ttmtdb(et_secs_i)
            yrs_J2000_tt = tt_secs_i/(daysec*yr)
            α_corr_v[i] = ( dRA + yrs_J2000_tt*pmRA/1000 ) / cos(δ_i_rad) # total debiasing correction in right ascension (arcsec)
            δ_corr_v[i] = dDEC + yrs_J2000_tt*pmDEC/1000 # total debiasing correction in declination (arcsec)
        end
    end
    return α_corr_v, δ_corr_v # arcsec, arcsec
end

# MPC minor planet optical observations fixed-width format
# Format is described in: https://minorplanetcenter.net/iau/info/OpticalObs.html
# Format is discussed thoroughly in S. Chesley et al. (2010) Icarus 210, p. 158-181
# Column 72: star catalog that was used in the reduction of the astrometry
# Column 15: measurement technique, or observation type
const mpc_format_mp = (mpnum=(1,5,String),
provdesig=(6,12,String),
discovery=(13,13,String),
publishnote=(14,14,String),
j2000=(15,15,String),
yr=(16,19,Int),
month=(21,22,Int),
day=(24,25,Int),
utc=(26,32,Float64),
rah=(33,34,Int),
ram=(36,37,Int),
ras=(39,44,Float64),
signdec=(45,45,String),
decd=(46,47,Int),
decm=(49,50,Int),
decs=(52,56,Float64),
info1=(57,65,String),
magband=(66,71,String),
catalog=(72,72,String),
info2=(73,77,String),
obscode=(78,80,String)
)

# MPC minor planet optical observations reader
readmp(mpcfile::String) = readfwf(mpcfile, mpc_format_mp)

# `mpc_catalog_nomenclature_2014` corresponds to debiasing tables included in Farnocchia et al. (2015)
mpc_catalog_codes_2014 = ["a", "b", "c", "d", "e", "g", "i", "j", "l", "m",
"o", "p", "q", "r",
"u", "v", "w", "L", "N"]

# `mpc_catalog_nomenclature_2018` corresponds to debiasing tables included in Eggl et al. (2020)
mpc_catalog_codes_2018 = ["a", "b", "c", "d", "e", "g", "i", "j", "l", "m",
    "n", "o", "p", "q", "r",
    "t", "u", "v", "w", "L", "N",
    "Q", "R", "S", "U", "W"
]

### MPC star catalog codes are were retrieved from the link below
### https://minorplanetcenter.net/iau/info/CatalogueCodes.html
mpc_catalog_codes = Dict(
    "a" =>   "USNO-A1.0",
    "b" =>   "USNO-SA1.0",
    "c" =>   "USNO-A2.0",
    "d" =>   "USNO-SA2.0",
    "e" =>   "UCAC-1",
    "f" =>   "Tycho-1",
    "g" =>   "Tycho-2",
    "h" =>   "GSC-1.0",
    "i" =>   "GSC-1.1",
    "j" =>   "GSC-1.2",
    "k" =>   "GSC-2.2",
    "l" =>   "ACT",
    "m" =>   "GSC-ACT",
    "n" =>   "SDSS-DR8",
    "o" =>   "USNO-B1.0",
    "p" =>   "PPM",
    "q" =>   "UCAC-4",
    "r" =>   "UCAC-2",
    "s" =>   "USNO-B2.0",
    "t" =>   "PPMXL",
    "u" =>   "UCAC-3",
    "v" =>   "NOMAD",
    "w" =>   "CMC-14",
    "x" =>   "Hipparcos 2",
    "y" =>   "Hipparcos",
    "z" =>   "GSC (version unspecified)",
    "A" =>   "AC",
    "B" =>   "SAO 1984",
    "C" =>   "SAO",
    "D" =>   "AGK 3",
    "E" =>   "FK4",
    "F" =>   "ACRS",
    "G" =>   "Lick Gaspra Catalogue",
    "H" =>   "Ida93 Catalogue",
    "I" =>   "Perth 70",
    "J" =>   "COSMOS/UKST Southern Sky Catalogue",
    "K" =>   "Yale",
    "L" =>   "2MASS",
    "M" =>   "GSC-2.3",
    "N" =>   "SDSS-DR7",
    "O" =>   "SST-RC1",
    "P" =>   "MPOSC3",
    "Q" =>   "CMC-15",
    "R" =>   "SST-RC4",
    "S" =>   "URAT-1",
    "T" =>   "URAT-2",
    "U" =>   "Gaia-DR1",
    "V" =>   "Gaia-DR2",
    "W" =>   "UCAC-5"
)
