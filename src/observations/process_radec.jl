@doc raw"""
    OpticalResidual{T <: Real, U <: Number}

An astrometric optical observed minus computed residual.

# Fields

- `ξ_α::U`: right ascension residual.
- `ξ_δ::U`: declination residual.
- `w_α::T`: right ascension weight.
- `w_δ::T`: declination weight.
- `relax_factor::T`: relaxation factor.
- `outlier::Bool`: whether the residual is an outlier or not.
"""
@auto_hash_equals struct OpticalResidual{T <: Real, U <: Number}
    ξ_α::U
    ξ_δ::U
    w_α::T
    w_δ::T
    relax_factor::T
    outlier::Bool
    # Inner constructor
    function OpticalResidual{T, U}(ξ_α::U, ξ_δ::U, w_α::T, w_δ::T, relax_factor::T = one(T),
                                   outlier::Bool = false) where {T <: Real, U <:  Number}
        new{T, U}(ξ_α, ξ_δ, w_α, w_δ, relax_factor, outlier)
    end
end
# Outer constructor
function OpticalResidual(ξ_α::U, ξ_δ::U, w_α::T, w_δ::T, relax_factor::T = one(T),
                         outlier::Bool = false) where {T <: Real, U <: Number}
    return OpticalResidual{T, U}(ξ_α, ξ_δ, w_α, w_δ, relax_factor, outlier)
end

# Evaluate methods
function evaluate(res::OpticalResidual{T, TaylorN{T}}, x::Vector{T}) where {T <: Real}
    return OpticalResidual(res.ξ_α(x), res.ξ_δ(x), res.w_α, res.w_δ, res.relax_factor, res.outlier)
end
(res::OpticalResidual{T, TaylorN{T}})(x::Vector{T}) where {T <: Real} = evaluate(res, x)

function evaluate(res::AbstractVector{OpticalResidual{T, TaylorN{T}}}, x::Vector{T}) where {T <: Real}
    res_new = Vector{OpticalResidual{T, T}}(undef, length(res))
    for i in eachindex(res)
        res_new[i] = evaluate(res[i], x)
    end
    return res_new
end
(res::AbstractVector{OpticalResidual{T, TaylorN{T}}})(x::Vector{T}) where {T <: Real} = evaluate(res, x)

# Print method for OpticalResidual
# Examples:
# α: -138.79801 δ: -89.80025
# α: -134.79450 δ: -91.42509 (outlier)
function show(io::IO, x::OpticalResidual{T, U}) where {T <: Real, U <: Number}
    outlier_flag = outlier(x) ? " (outlier)" : ""
    print(io, "α: ", @sprintf("%+.5f", cte(x.ξ_α)), " δ: ",
          @sprintf("%+.5f", cte(x.ξ_δ)), outlier_flag)
end

@doc raw"""
    unfold(ξs::Vector{OpticalResidual{T, U}}) where {T <: Real, U <: Number}

Concatenate right ascension and declination residuals for an orbit fit.
"""
function unfold(ξs::AbstractVector{OpticalResidual{T, U}}) where {T <: Real, U <: Number}
    # Number of non outliers
    L = count(x -> !x.outlier, ξs)
    # Vector of residuals
    res = Vector{U}(undef, 2*L)
    # Vector of weights
    w = Vector{T}(undef, 2*L)
    # Global counter
    k = 1
    # Fill residuals and weights
    for i in eachindex(ξs)
        if !ξs[i].outlier
            # Right ascension
            res[k] = ξs[i].ξ_α
            w[k] = ξs[i].w_α / ξs[i].relax_factor
            # Declination
            res[k+L] = ξs[i].ξ_δ
            w[k+L] = ξs[i].w_δ / ξs[i].relax_factor
            # Update global counter
            k += 1
        end
    end

    return res, w
end

# Functions to get specific fields of a OpticalResidual object
ra(res::OpticalResidual{T, U}) where {T <: Real, U <: Number} = res.ξ_α
dec(res::OpticalResidual{T, U}) where {T <: Real, U <: Number} = res.ξ_δ
weight_ra(res::OpticalResidual{T, U}) where {T <: Real, U <: Number} = res.w_α
weight_dec(res::OpticalResidual{T, U}) where {T <: Real, U <: Number} = res.w_δ
relax_factor(res::OpticalResidual{T, U}) where {T <: Real, U <: Number} = res.relax_factor
outlier(res::OpticalResidual{T, U}) where {T <: Real, U <: Number} = res.outlier

euclid3D(x::Vector{U}) where {U <: Number} = sqrt(x[1]*x[1] + x[2]*x[2] + x[3]*x[3])
dot3D(x::Vector{U}, y::Vector{V}) where {U, V <: Number} = x[1]*y[1] + x[2]*y[2] + x[3]*y[3]

@doc raw"""
    compute_radec(observatory::ObservatoryMPC{T}, t_r_utc::DateTime; kwargs...) where {T <: AbstractFloat}
    compute_radec(obs::RadecMPC{T}; kwargs...) where {T <: AbstractFloat}
    compute_radec(obs::Vector{RadecMPC{T}}; kwargs...) where {T <: AbstractFloat}

Compute astrometric right ascension and declination (both in arcsec) for a set of MPC-formatted
observations. Corrections due to Earth orientation, LOD, polar motion are considered in
computations.

# Arguments

- `observatory::ObservatoryMPC{T}`: observation site.
- `t_r_utc::DateTime`: UTC time of astrometric observation.
- `obs::Vector{RadecMPC{T}}`: optical observation(s).

# Keyword arguments

- `niter::Int = 5`: number of light-time solution iterations.
- `xve::EarthEph = earthposvel`: Earth ephemeris.
- `xvs::SunEph = sunposvel`: Sun ephemeris.
- `xva::AstEph`: asteroid ephemeris.

All ephemeris must take [et seconds since J2000] and return [barycentric position in km
and velocity in km/sec].
"""
function compute_radec(observatory::ObservatoryMPC{T}, t_r_utc::DateTime; niter::Int = 5,
                       xvs::SunEph = sunposvel, xve::EarthEph = earthposvel,
                       xva::AstEph) where {T <: AbstractFloat, SunEph, EarthEph, AstEph}
    # Transform receiving time from UTC to TDB seconds since J2000
    et_r_secs = datetime2et(t_r_utc)
    # Sun barycentric position and velocity at receive time
    rv_s_t_r = xvs(et_r_secs)
    r_s_t_r = rv_s_t_r[1:3]
    # Earth's barycentric position and velocity at receive time
    rv_e_t_r = xve(et_r_secs)
    r_e_t_r = rv_e_t_r[1:3]
    # Asteroid barycentric position and velocity at receive time
    rv_a_t_r = xva(et_r_secs)
    r_a_t_r = rv_a_t_r[1:3]
    # Compute geocentric position/velocity of receiving antenna in inertial frame [km, km/s]
    RV_r = obsposvelECI(observatory, et_r_secs)
    R_r = RV_r[1:3]
    # Receiver barycentric position and velocity at receive time
    r_r_t_r = r_e_t_r + R_r

    # Down-leg iteration
    # τ_D first approximation
    # See equation (1) of https://doi.org/10.1086/116062
    ρ_vec_r = r_a_t_r - r_r_t_r
    ρ_r = euclid3D(ρ_vec_r)
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
        ρ_r = euclid3D(ρ_vec_r)

        # Compute down-leg Shapiro delay
        # NOTE: when using PPN, substitute 2 -> 1+γ in expressions for Shapiro delay,
        # Δτ_rel_[D|U]

        # Earth's position at t_r
        e_D_vec  = r_r_t_r - r_s_t_r
        # Heliocentric distance of Earth at t_r
        e_D = euclid3D(e_D_vec)
        # Barycentric position of Sun at estimated bounce time
        rv_s_t_b = xvs(et_b_secs)
        r_s_t_b = rv_s_t_b[1:3]
        # Heliocentric position of asteroid at t_b
        p_D_vec  = r_a_t_b - r_s_t_b
        # Heliocentric distance of asteroid at t_b
        p_D = euclid3D(p_D_vec)
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
        p_dot_23 = dot3D(ρ_vec_r, v_a_t_b)/ρ_r
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
    ρ_r = euclid3D(ρ_vec_r)

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
    q_vec = Q_vec/ euclid3D(Q_vec)
    E_H = euclid3D(E_H_vec)
    e_vec = E_H_vec/E_H
    # See ESAA 2014, equation (7.115)
    g1 = (2μ_DE430[su]/(c_au_per_day^2))/(E_H/au)
    g2 = 1 + dot3D(q_vec, e_vec)
    # See ESAA 2014, equation (7.116)
    u1_vec = U_norm*(  u_vec + (g1/g2)*( dot3D(u_vec,q_vec)*e_vec - dot3D(e_vec,u_vec)*q_vec )  )
    u1_norm = euclid3D(u1_vec)

    # Compute right ascension, declination angles
    α_rad_ = mod2pi(atan(u1_vec[2], u1_vec[1]))
    α_rad = mod2pi(α_rad_)          # right ascension (rad)
    δ_rad = asin(u1_vec[3]/u1_norm) # declination (rad)

    δ_as = rad2arcsec(δ_rad) # rad -> arcsec + debiasing
    α_as = rad2arcsec(α_rad) # rad -> arcsec + debiasing

    return α_as, δ_as # right ascension, declination both in arcsec
end

function compute_radec(obs::RadecMPC{T}; kwargs...) where {T <: AbstractFloat}
    return compute_radec(obs.observatory, obs.date; kwargs...)
end

function compute_radec(obs::Vector{RadecMPC{T}}; xva::AstEph, kwargs...) where {T <: AbstractFloat, AstEph}

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
        vra[i], vdec[i] = compute_radec(obs[i]; xva = xva, kwargs...)
    end

    return vra, vdec # arcsec, arcsec
end

@doc raw"""
    select_debiasing_table(debias_table::String = "2018")

Return the catalogue codes, truth catalogue, resolution and bias matrix of the corresponding
debias table. The possible values for `debias_table` are:
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
    debiasing(obs::RadecMPC{T}, mpc_catalogue_codes_201X::Vector{String}, truth::String,
              resol::Resolution, bias_matrix::Matrix{T}) where {T <: AbstractFloat}

Return total debiasing correction in right ascension and declination (both in arcsec).

# Arguments

- `obs::RadecMPC{T}`: optical observation.
- `mpc_catalogue_codes_201X::Vector{String}`: catalogues present in debiasing table.
- `truth::String`: truth catalogue of debiasing table.
- `resol::Resolution`: resolution.
- `bias_matrix::Matrix{T}`: debiasing table.
"""
function debiasing(obs::RadecMPC{T}, mpc_catalogue_codes_201X::Vector{String}, truth::String,
                   resol::Resolution, bias_matrix::Matrix{T}) where {T <: AbstractFloat}

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
        tt_secs_i = et_secs_i - ttmtdb(et_secs_i/daysec)
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
    relax_factor(radec::Vector{RadecMPC{T}}) where {T <: AbstractFloat}

Return a relax factor for each element of `radec` quantifying the correlation between
observations taken on the same night by the same observatory.

!!! reference
    See https://doi.org/10.1016/j.icarus.2017.05.021.
"""
function relax_factor(radec::Vector{RadecMPC{T}}) where {T <: AbstractFloat}
    # Convert to DataFrame
    df = DataFrame(radec)
    # Group by observatory and TimeOfDay
    df.TimeOfDay = TimeOfDay.(radec)
    gdf = groupby(df, [:observatory, :TimeOfDay])
    # Number of observations per tracklet
    cdf = combine(gdf, nrow)
    # Count observations in each group
    Nv = cdf[gdf.groups, :nrow]
    # Relaxation factor
    return map(x -> x > 4.0 ? x/4.0 : 1.0, Nv)
end

# Angle difference taking into account the discontinuity in [0, 2π] -> [0, 2π]
# x and y must be in arcsec
function anglediff(x::T, y::S) where {T, S <: Number}
    # Signed difference
    Δ = x - y
    # Absolute difference
    Δ_abs = abs(Δ)
    # Reflection
    if Δ_abs > 648_000 # Half circle in arcsec
        return -sign(cte(Δ)) * (1_296_000 - Δ_abs)
    else
        return Δ
    end
end

@doc raw"""
    residuals(obs::Vector{RadecMPC{T}}; kwargs...) where {T <: AbstractFloat, AstEph}

Compute O-C residuals for optical astrometry. Corrections due to Earth orientation, LOD,
polar motion are computed by default.

See also [`compute_radec`](@ref), [`debiasing`](@ref), [`w8sveres17`](@ref) and
[`Healpix.ang2pixRing`](@ref).

# Arguments

- `obs::Vector{RadecMPC{T}}`: vector of observations.

# Keyword arguments

- `niter::Int = 5`: number of light-time solution iterations.
- `debias_table::String = "2018"`: possible values are:
    - `2014` corresponds to https://doi.org/10.1016/j.icarus.2014.07.033,
    - `2018` corresponds to https://doi.org/10.1016/j.icarus.2019.113596,
    - `hires2018` corresponds to https://doi.org/10.1016/j.icarus.2019.113596.
- `xve::EarthEph = earthposvel`: Earth ephemeris.
- `xvs::SunEph = sunposvel`: Sun ephemeris.
- `xva::AstEph`: asteroid ephemeris.

All ephemeris must take [et seconds since J2000] and return [barycentric position in km
and velocity in km/sec].
"""
function residuals(obs::Vector{RadecMPC{T}}; debias_table::String = "2018", xva::AstEph,
                   kwargs...) where {T <: AbstractFloat, AstEph}
    mpc_catalogue_codes_201X, truth, resol, bias_matrix = select_debiasing_table(debias_table)
    return residuals(obs, mpc_catalogue_codes_201X, truth, resol, bias_matrix; xva, kwargs...)
end
function residuals(obs::Vector{RadecMPC{T}}, mpc_catalogue_codes_201X::Vector{String}, truth::String,
                   resol::Resolution, bias_matrix::Matrix{T}; xva::AstEph, kwargs...) where {T <: AbstractFloat, AstEph}

    # Number of observations
    N_obs = length(obs)

    # UTC time of first astrometric observation
    utc1 = date(obs[1])
    # TDB seconds since J2000.0 for first astrometric observation
    et1 = datetime2et(utc1)
    # Asteroid ephemeris at et1
    a1_et1 = xva(et1)[1]
    # Type of asteroid ephemeris
    S = typeof(a1_et1)

    # Relax factors
    rex = relax_factor(obs)

    # Vector of residuals
    res = Vector{OpticalResidual{T, S}}(undef, N_obs)

    # Iterate over the observations
    Threads.@threads for i in 1:N_obs

        # Observed ra/dec
        α_obs = rad2arcsec(ra(obs[i]))   # arcsec
        δ_obs = rad2arcsec(dec(obs[i]))  # arcsec

        # Computed ra/dec
        α_comp, δ_comp = compute_radec(obs[i]; xva = xva, kwargs...)   # arcsec

        # Debiasing corrections
        α_corr, δ_corr = debiasing(obs[i], mpc_catalogue_codes_201X, truth, resol, bias_matrix)

        # Statistical weights from Veres et al. (2017)
        w8 = w8sveres17(obs[i])

        # O-C residual ra/dec
        # Note: ra is multiplied by a metric factor cos(dec) to match the format of debiasing corrections
        res[i] = OpticalResidual(
            anglediff(α_obs, α_comp) * cos(dec(obs[i])) - α_corr,
            δ_obs - δ_comp - δ_corr,
            1 / w8^2,
            1 / w8^2,
            rex[i],
            false
        )

    end

    return res
end