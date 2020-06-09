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

function radec(astopticalobsfile::String,
        niter::Int=10; pm::Bool=true, lod::Bool=true, eocorr::Bool=true,
        xve::Function=earth_pv, xvs::Function=sun_pv, xva::Function=apophis_pv_197)
    # astopticalobsfile = "tholenetal2013_radec_data.dat"
    astopticalobsdata = readdlm(astopticalobsfile, ',', comments=true)

    utc1 = DateTime(astopticalobsdata[1,1]) + Microsecond( round(1e6daysec*astopticalobsdata[1,2]) )
    et1 = str2et(string(utc1))
    a1_et1 = xva(et1)[1]
    S = typeof(a1_et1)

    n_optical_obs = size(astopticalobsdata)[1]

    vra = Array{S}(undef, n_optical_obs)
    vdec = Array{S}(undef, n_optical_obs)

    for i in 1:n_optical_obs
        utc_i = DateTime(astopticalobsdata[i,1]) + Microsecond( round(1e6daysec*astopticalobsdata[i,2]) )
        station_code_i = string(astopticalobsdata[i,17])
        vra[i], vdec[i] = radec(station_code_i, utc_i, niter, pm=pm, lod=lod,
            eocorr=eocorr, xve=xve, xvs=xvs, xva=xva)
    end

    return vra, vdec
end

function radec_mpc_vokr15(niter::Int=10; pm::Bool=true, lod::Bool=true,
        eocorr::Bool=true, xve::Function=earth_pv, xvs::Function=sun_pv,
        xva::Function=apophis_pv_197)

    astopticalobsfile = joinpath(dirname(pathof(Apophis)), "../vokrouhlickyetal2015_mpc.dat")
    vokr15 = readdlm(astopticalobsfile, ',', comments=true)

    utc1 = DateTime(vokr15[1,4], vokr15[1,5], vokr15[1,6]) + Microsecond( round(1e6*86400*vokr15[1,7]) )
    et1 = str2et(string(utc1))
    a1_et1 = xva(et1)[1]
    S = typeof(a1_et1)

    n_optical_obs = size(vokr15)[1]

    vra = Array{S}(undef, n_optical_obs)
    vdec = Array{S}(undef, n_optical_obs)

    for i in 1:n_optical_obs
        utc_i = DateTime(vokr15[i,4], vokr15[i,5], vokr15[i,6]) + Microsecond( round(1e6*86400*vokr15[i,7]) )
        station_code_i = string(vokr15[i,20])
        vra[i], vdec[i] = radec(station_code_i, utc_i, niter, pm=pm, lod=lod,
            eocorr=eocorr, xve=xve, xvs=xvs, xva=xva)
    end

    return vra, vdec
end

# Compute optical astrometric ra/dec ephemeris for a set of observations in a MPC-formatted file
function radec_mpc(astopticalobsfile, niter::Int=10; pm::Bool=true, lod::Bool=true,
        eocorr::Bool=true, xve::Function=earth_pv,
        xvs::Function=sun_pv, xva::Function=apophis_pv_197)

    obs_df = readmp(astopticalobsfile)

    utc1 = DateTime(obs_df.yr[1], obs_df.month[1], obs_df.day[1]) + Microsecond( round(1e6*86400*obs_df.utc[1]) )
    et1 = str2et(string(utc1))
    a1_et1 = xva(et1)[1]
    S = typeof(a1_et1)

    n_optical_obs = nrow(obs_df)

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

# Compute ra/dec debiasing corrections following Eggl et al. (2019)
function radec_mpc_corr(astopticalobsfile::String, table::String="2018")

    obs_df = readmp(astopticalobsfile)
    n_optical_obs = nrow(obs_df)

    α_corr_v = Array{Float64}(undef, n_optical_obs)
    δ_corr_v = Array{Float64}(undef, n_optical_obs)

    # Select debiasing table: 2014 corresponds to Farnocchia et al. (2015), 2018 corresponds to Eggl et al. (2020)
    if table == "2018"
        bias_file = joinpath(dirname(pathof(Apophis)), "../debias/debias_2018/bias.dat")
        mpc_catalog_nomenclature = mpc_catalog_nomenclature_2018
    elseif table == "2014"
        bias_file = joinpath(dirname(pathof(Apophis)), "../debias/debias_2014/bias.dat")
        mpc_catalog_nomenclature = mpc_catalog_nomenclature_2014
    else
        @error "Unknown debias table: $table"
    end

    bias_matrix = readdlm(bias_file, comment_char='!', comments=true)
    NSIDE= 64 #The healpix tesselation resolution of the bias map from Eggl et al. (2019)
    resol = Resolution(NSIDE) # initialize healpix Resolution variable

    for i in 1:n_optical_obs
        α_i_as = 15(obs_df.rah[i] + obs_df.ram[i]/60 + obs_df.ras[i]/3600) # deg
        δ_i_as = obs_df.decd[i] + obs_df.decm[i]/60 + obs_df.decs[i]/3600 # deg
        α_i_rad = deg2rad(α_i_as) #rad
        δ_i_rad = deg2rad(δ_i_as) #rad
        # get pixel tile index, assuming iso-latitude rings indexing, which is the formatting in tiles.dat
        # substracting 1 from the returned value of `ang2pixRing` corresponds to 0-based indexing, as in tiles.dat
        # not substracting 1 from the returned value of `ang2pixRing` corresponds to 1-based indexing, as in Julia
        pix_ind = ang2pixRing(resol, π/2-δ_i_rad, α_i_rad)
        if obs_df.catalog[i] == "t" && table == "2014"
            α_corr_v[i] = 0.0
            δ_corr_v[i] = 0.0
            continue
        end
        cat_ind = mpc_catalog_nomenclature[obs_df.catalog[i]][1]
        @show pix_ind, cat_ind
        utc_i = DateTime(obs_df.yr[i], obs_df.month[i], obs_df.day[i]) + Microsecond( round(1e6*86400*obs_df.utc[i]) )
        # read dRA, pmRA, dDEC, pmDEC data from bias.dat
        # dRA, position correction in RA*cos(DEC) at epoch J2000.0 [arcsec];
        # dDEC, position correction in DEC at epoch J2000.0 [arcsec];
        # pmRA, proper motion correction in RA*cos(DEC) [mas/yr];
        # pmDEC, proper motion correction in DEC [mas/yr].
        dRA, dDEC, pmRA, pmDEC = bias_matrix[pix_ind, 4*cat_ind-3:4*cat_ind]
        @show dRA, dDEC, pmRA, pmDEC
        et_secs_i = str2et(string(utc_i))
        yrs_J2000_tdb = et_secs_i/(daysec*yr)
        α_corr_v[i] = ( dRA + yrs_J2000_tdb*pmRA/1000 ) / cos(δ_i_rad) # total debiasing correction in right ascension (arcsec)
        δ_corr_v[i] = dDEC + yrs_J2000_tdb*pmDEC/1000 # total debiasing correction in declination (arcsec)
    end
    return α_corr_v, δ_corr_v # arcsec, arcsec
end

### FWF reader, due to @aplavin
### https://gist.github.com/aplavin/224a31ea457b6e0ef0f4c1a20bd28850
convert_val(::Type{String}, val::String) = val
convert_val(::Type{Symbol}, val::String) = Symbol(val)
convert_val(typ::Type{<:Integer}, val::String) = parse(typ, val)
convert_val(typ::Type{<:AbstractFloat}, val::String) = parse(typ, replace(val, "D" => "E"))  # tables output from Fortran often have floats as 1D+5 instead of 1E+5
# # usage:
# read_fwf(
#   "file.tab",
#   (
#     colname1=(1, 20, String),
#     colname2=(100, 120, Int),
#     # ...
#   ),
#   skiprows=[1, 2, 3, 7]
# )
function readfwf(io, colspecs; skiprows=[], missingstrings=[])
    cols = Dict(
        k => Vector{Union{typ, Missing}}(undef, 0)
        for (k, (from, to, typ)) in pairs(colspecs)
    )
    for (irow, line) in eachline(io) |> enumerate
        if irow ∈ skiprows continue end
        for (k, (from, to, typ)) in pairs(colspecs)
            s_val = from <= length(line) ? line[from:min(length(line), to)] : ""
            f_val = s_val in missingstrings ? missing : convert_val(typ, s_val)
            push!(cols[k], f_val)
        end
    end
    DataFrame([k => identity.(cols[k]) for k in keys(colspecs)])
end

# MPC minor planet optical observations fixed-width format
# References: https://minorplanetcenter.net/iau/info/OpticalObs.html
# S. Chesley et al. (2010) Icarus 210, p. 158-181
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
utc=(26,31,Float64),
rah=(33,34,Int),
ram=(36,37,Int),
ras=(39,44,Float64),
decd=(45,47,Int),
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

mpc_catalog_flags_2014 = ["a", "b", "c", "d", "e", "g", "i", "j", "l", "m",
"o", "p", "q", "r", "u", "v", "w", "L", "N"]

catalog_names_2014 = ["USNO-A1.0", "USNO-SA1.0", "USNO-A2.0", "USNO-SA2.0",
"UCAC-1", "Tycho-2", "GSC-1.1", "GSC-1.2", "ACT", "GSC-ACT", "USNO-B1.0",
"PPM", "UCAC-4", "UCAC-2", "UCAC-3", "NOMAD", "CMC-14", "2MASS",
"SDSS-DR7"]

mpc_catalog_nomenclature_2014 = Dict{String, Tuple{Int, String}}()
for i in eachindex(mpc_catalog_flags_2014)
    mpc_catalog_nomenclature_2014[mpc_catalog_flags_2014[i]] = (i, catalog_names_2014[i])
end

mpc_catalog_flags_2018 = ["a", "b", "c", "d", "e", "g", "i", "j", "l", "m",
    "n", "o", "p", "q", "r", "t", "u", "v", "w", "L",
    "N", "Q", "R", "S", "U", "W"
]

catalog_names_2018 = ["USNO-A1.0", "USNO-SA1.0", "USNO-A2.0", "USNO-SA2.0", "UCAC-1",
    "Tycho-2", "GSC-1.1", "GSC-1.2", "ACT", "GSC-ACT", "SDSS-DR8", "USNO-B1.0",
    "PPM", "UCAC-4", "UCAC-2", "PPMXL", "UCAC-3", "NOMAD", "CMC-14",
    "2MASS", "SDSS-DR7", "CMC-15", "SST-RC4", "URAT-1", "Gaia-DR1", "UCAC-5"
]

mpc_catalog_nomenclature_2018 = Dict{String, Tuple{Int, String}}()
for i in eachindex(mpc_catalog_flags_2018)
    mpc_catalog_nomenclature_2018[mpc_catalog_flags_2018[i]] = (i, catalog_names_2018[i])
end
