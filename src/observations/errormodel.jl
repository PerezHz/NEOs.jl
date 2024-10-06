@doc raw"""
    AbstractErrorModel{T <: Real}

Supertype for the optical astrometry error models API.
"""
abstract type AbstractErrorModel{T <: Real} end

# Optical astrometry weighting schemes API
# All weighting schemes `W{T}` must:
# 1. be mutable subtypes of `AbstractWeightingScheme{T}`,
# 2. have a `w8s::Vector{T}` field,
# 3. implement a `W(::AbstractVector{RadecMPC{T}})` constructor,
# 4. override `getid(::W)`,
# 5. override `update!(::W{T}, ::AbstractVector{RadecMPC{T}})`.

@doc raw"""
    AbstractWeightingScheme{T} <: AbstractErrorModel{T}

Supertype for the optical astrometry weighting schemes API.
"""
abstract type AbstractWeightingScheme{T} <: AbstractErrorModel{T} end

# Print method for AbstractWeightingScheme
show(io::IO, w::AbstractWeightingScheme) = print(io, getid(w),
    " weighting scheme with ", length(w.w8s), " observations")

@doc raw"""
    UniformWeights{T} <: AbstractWeightingScheme{T}

Uniform optical astrometry weighting scheme.
"""
mutable struct UniformWeights{T} <: AbstractWeightingScheme{T}
    w8s::Vector{T}
    # Default constructor
    UniformWeights(radec::AbstractVector{RadecMPC{T}}) where {T <: Real} =
        new{T}(ones(T, length(radec)))
end

# Override getid
getid(w::UniformWeights) = "Uniform"

# Override update!
function update!(w::UniformWeights{T},
    radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
    w.w8s = ones(T, length(radec))
    return nothing
end

@doc raw"""
    Veres17{T} <: AbstractWeightingScheme{T}

Veres et al. (2017) optical astrometry weighting scheme.

!!! reference
    See https://doi.org/10.1016/j.icarus.2017.05.021.
"""
mutable struct Veres17{T} <: AbstractWeightingScheme{T}
    w8s::Vector{T}
    # Default constructor
    Veres17(radec::AbstractVector{RadecMPC{T}}) where {T <: Real} =
        new{T}(w8sveres17(radec))
end

# Override getid
getid(w::Veres17) = "Veres et al. (2017)"

# Override update!
function update!(w::Veres17{T},
    radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
    w.w8s = w8sveres17(radec)
    return nothing
end

@doc raw"""
    σsveres17(obs::RadecMPC)

Return the statistical uncertainty of `obs` according to Veres et al. (2017).

!!! reference
    See https://doi.org/10.1016/j.icarus.2017.05.021.
"""
function σsveres17(obs::RadecMPC)

    obscode = obs.observatory.code
    dt_utc_obs = obs.date
    catalogue = obs.catalogue.code

    # Unit weight (arcseconds)
    w = 1.0
    # Table 2: epoch-dependent astrometric residuals
    if obscode == "703"
        return Date(dt_utc_obs) < Date(2014,1,1) ? w : 0.8w
    elseif obscode ∈ ("691", "291") # Spacewatch, Kitt Peak, Arizona
        return Date(dt_utc_obs) < Date(2003,1,1) ? 0.6w : 0.5w
    elseif obscode == "644"
        return Date(dt_utc_obs) < Date(2003,9,1) ? 0.6w : 0.4w
    # Table 3: most active CCD asteroid observers
    elseif obscode ∈ ("704", "C51", "J75")
        return w
    elseif obscode == "G96"
        return 0.5w
    elseif obscode ∈ ("F51", "F52") # Pan-STARRS 1 & 2, Haleakala, Hawaii
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
    elseif obscode ∈ ("J04", "K92", "K93", "Q63", "Q64", "V37", "W85", "W86", "W87",
        "K91", "E10", "F65") # Tenerife + Las Cumbres
        return 0.4w
    elseif obscode ∈ ("689", "950", "W84")
        return 0.5w
    # Applies only to program code assigned to M. Micheli
    #elseif obscode ∈ ("G83", "309")
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
    rexveres17(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}

Return the relax factor for each element of `radec` according to
Veres et al. (2017), which mitigates unresolved systematic errors
in observations taken on the same night by the same observatory.

!!! reference
    See https://doi.org/10.1016/j.icarus.2017.05.021.
"""
function rexveres17(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
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
    return map(N -> N > 4.0 ? sqrt(N/4.0) : 1.0, Nv)
end

@doc raw"""
    w8sveres17(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}

Return the statistical weight of each element of `radec` according
to Veres et al. (2017).

!!! reference
    See https://doi.org/10.1016/j.icarus.2017.05.021.
"""
function w8sveres17(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
    σs = σsveres17.(radec)
    rex = rexveres17(radec)
    return 1 ./ (rex .* σs) .^ 2
end

# Optical astrometry debiasing schemes API
# All debiasing schemes `D{T}` must:
# 1. be mutable subtypes of `AbstractDebiasingScheme{T}`,
# 2. have a `bias::Vector{Tuple{T, T}}` field,
# 3. implement a `D(::AbstractVector{RadecMPC{T}})` constructor,
# 4. override `getid(::D)`,
# 5. override `update!(::D{T}, ::AbstractVector{RadecMPC{T}})`.

@doc raw"""
    AbstractDebiasingScheme{T} <: AbstractErrorModel{T}

Supertype for the optical astrometry debiasing schemes API.
"""
abstract type AbstractDebiasingScheme{T} <: AbstractErrorModel{T} end

# Print method for AbstractDebiasingScheme
show(io::IO, d::AbstractDebiasingScheme) = print(io, getid(d),
    " debiasing scheme with ", length(d.bias), " observations")

@doc raw"""
    Farnocchia15{T} <: AbstractDebiasingScheme{T}

Farnocchia et al. (2015) optical astrometry debiasing scheme.

!!! reference
    See https://doi.org/10.1016/j.icarus.2014.07.033.
"""
mutable struct Farnocchia15{T} <: AbstractDebiasingScheme{T}
    bias::Vector{Tuple{T, T}}
    catcodes::Vector{String}
    truth::String
    resol::Resolution
    table::Matrix{T}
    # Default constructor
    function Farnocchia15(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
        catcodes, truth, resol, table = select_debiasing_table("2014")
        bias = debiasing.(radec, Ref(catcodes), Ref(truth),
            Ref(resol), Ref(table))
        return new{T}(bias, catcodes, truth, resol, table)
    end
end

# Override getid
getid(d::Farnocchia15) = "Farnocchia et al. (2015)"

# Override update!
function update!(d::Farnocchia15{T},
    radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
    d.bias = debiasing.(radec, Ref(d.catcodes), Ref(d.truth),
    Ref(d.resol), Ref(d.table))
    return nothing
end

@doc raw"""
    Eggl20{T} <: AbstractDebiasingScheme{T}

Eggl et al. (2020) optical astrometry debiasing scheme.

!!! reference
    See https://doi.org/10.1016/j.icarus.2019.113596.
"""
mutable struct Eggl20{T} <: AbstractDebiasingScheme{T}
    bias::Vector{Tuple{T, T}}
    hires::Bool
    catcodes::Vector{String}
    truth::String
    resol::Resolution
    table::Matrix{T}
    # Default constructor
    function Eggl20(radec::AbstractVector{RadecMPC{T}},
        hires::Bool = false) where {T <: Real}
        id = hires ? "hires2018" : "2018"
        catcodes, truth, resol, table = select_debiasing_table(id)
        bias = debiasing.(radec, Ref(catcodes), Ref(truth),
            Ref(resol), Ref(table))
        return new{T}(bias, hires, catcodes, truth, resol, table)
    end
end

# Override getid
getid(d::Eggl20) = string("Eggl et al. (2020)", d.hires ?
    " [high resolution]" : "")

# Override update!
function update!(d::Eggl20{T},
    radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
    d.bias = debiasing.(radec, Ref(d.catcodes), Ref(d.truth),
    Ref(d.resol), Ref(d.table))
    return nothing
end

@doc raw"""
    select_debiasing_table(id::String)

Return the catalogue codes, truth catalogue, resolution and bias matrix of the
corresponding debiasing scheme. The possible values for `id` are:
- `2014` corresponds to https://doi.org/10.1016/j.icarus.2014.07.033,
- `2018` corresponds to https://doi.org/10.1016/j.icarus.2019.113596 (default),
- `hires2018` corresponds to https://doi.org/10.1016/j.icarus.2019.113596
    (high resolution).
"""
function select_debiasing_table(id::String = "2018")
    # Debiasing tables are loaded "lazily" via Julia artifacts,
    # according to rules in Artifacts.toml
    if id == "2018"
        debias_path = artifact"debias_2018"
        catcodes = mpc_catalogue_codes_2018
        # The healpix tesselation resolution of the bias map from
        # https://doi.org/10.1016/j.icarus.2019.113596
        NSIDE = 64
        # In 2018 debias table Gaia DR2 catalogue is regarded as the truth
        truth = "V"
    elseif id == "hires2018"
        debias_path = artifact"debias_hires2018"
        catcodes = mpc_catalogue_codes_2018
        # The healpix tesselation resolution of the high-resolution bias map from
        # https://doi.org/10.1016/j.icarus.2019.113596
        NSIDE = 256
        # In 2018 debias table Gaia DR2 catalogue is regarded as the truth
        truth = "V"
    elseif id == "2014"
        debias_path = artifact"debias_2014"
        catcodes = mpc_catalogue_codes_2014
        # The healpix tesselation resolution of the bias map from
        # https://doi.org/10.1016/j.icarus.2014.07.033
        NSIDE = 64
        # In 2014 debias table PPMXL catalogue is regarded as the truth
        truth = "t"
    else
        @error "Unknown bias map: $(id). Possible values are `2014`, `2018` \
            and `hires2018`."
    end

    # Debias table file
    bias_file = joinpath(debias_path, "bias.dat")
    # Read bias matrix
    table = readdlm(bias_file, comment_char = '!', comments = true)
    # Initialize healpix Resolution variable
    resol = Resolution(NSIDE)
    # Compatibility between bias matrix and resolution
    @assert size(table) == (resol.numOfPixels, 4length(catcodes)) "\
        Bias table file $bias_file dimensions do not match expected parameter \
        NSIDE = $NSIDE and/or number of catalogs in table."

    return catcodes, truth, resol, table
end

@doc raw"""
    debiasing(obs::RadecMPC{T}, catcodes::Vector{String}, truth::String,
        resol::Resolution, table::Matrix{T}) where {T <: Real}

Return total debiasing correction [arcsec] in both right ascension and declination.

## Arguments

- `obs::RadecMPC{T}`: optical observation.
- `catcodes::Vector{String}`: catalogues present in debiasing table.
- `truth::String`: truth catalogue of debiasing table.
- `resol::Resolution`: resolution.
- `table::Matrix{T}`: debiasing table.
"""
function debiasing(obs::RadecMPC{T}, catcodes::Vector{String},
    truth::String, resol::Resolution, table::Matrix{T}) where {T <: Real}

    # Catalogue code
    catcode = obs.catalogue.code

    # If star catalogue is not present in debiasing table,
    # then set corrections equal to zero
    if (catcode ∉ catcodes) && catcode != "Y"
        # Catalogue exists in CATALOGUES_MPC[]
        if !isunknown(obs.catalogue)
            # Truth catalogue is not present in debiasing table
            # but it does not send a warning
            if catcode != truth
                @warn "Catalogue $(obs.catalogue.name) not found in debiasing table. \
                Setting debiasing corrections equal to zero."
            end
        # Unknown catalogue
        elseif catcode == ""
            @warn "Catalog information not available in observation record. \
            Setting debiasing corrections equal to zero."
        # Catalogue code is not empty but it does not match an MPC catalogue code either
        else
            @warn "Catalog code $catcode does not correspond to MPC catalogue code. \
            Setting debiasing corrections equal to zero."
        end
        α_corr, δ_corr = zero(T), zero(T)
    # If star catalogue is present in debiasing table, then compute corrections
    else
        # Get pixel tile index, assuming iso-latitude rings indexing, which is the
        # formatting in `tiles.dat`. Substracting 1 from the returned value of
        # `ang2pixRing` corresponds to 0-based indexing, as in `tiles.dat`; not
        # substracting 1 from the returned value of `ang2pixRing` corresponds to
        # 1-based indexing, as in Julia. Since we use pix_ind to get the corresponding
        # row number in `bias.dat`, it's not necessary to substract 1.
        pix_ind = ang2pixRing(resol, π/2 - obs.δ, obs.α)

        # Handle edge case: in new MPC catalogue nomenclature, "UCAC-5" -> "Y";
        # but in debias tables "UCAC-5" -> "W"
        if catcode == "Y"
            cat_ind = findfirst(x -> x == "W", catcodes)
        else
            cat_ind = findfirst(x -> x == catcode, catcodes)
        end

        # Read dRA, pmRA, dDEC, pmDEC data from bias.dat
        # dRA: position correction in RA * cos(DEC) at epoch J2000.0 [arcsec]
        # dDEC: position correction in DEC at epoch J2000.0 [arcsec]
        # pmRA: proper motion correction in RA*cos(DEC) [mas/yr]
        # pmDEC: proper motion correction in DEC [mas/yr]
        dRA, dDEC, pmRA, pmDEC = table[pix_ind, 4*cat_ind-3:4*cat_ind]
        # Seconds since J2000 (TDB)
        et_secs_i = dtutc2et(obs.date)
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

debiasing(obs::RadecMPC{T}, catcodes::RefValue{Vector{String}},
    truth::RefValue{String}, resol::RefValue{Resolution},
    table::RefValue{Matrix{T}}) where {T <: Real} =
    debiasing(obs, catcodes[], truth[], resol[], table[])