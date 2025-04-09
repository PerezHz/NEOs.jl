include("weightingscheme.jl")

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
getid(::Farnocchia15) = "Farnocchia et al. (2015)"

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

@doc raw"""
    ZeroDebiasing{T} <: AbstractDebiasingScheme{T}

Zero optical astrometry debiasing scheme.
"""
mutable struct ZeroDebiasing{T} <: AbstractDebiasingScheme{T}
    bias::Vector{Tuple{T, T}}
    # Default constructor
    function ZeroDebiasing(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
        bias = [(zero(T), zero(T)) for _ in eachindex(radec)]
        return new{T}(bias)
    end
end

# Override getid
getid(::ZeroDebiasing) = "Zero"

# Override update!
function update!(d::ZeroDebiasing{T},
    radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
    d.bias = [(zero(T), zero(T)) for _ in eachindex(radec)]
    return nothing
end

@doc raw"""
    NEOCCDebiasing{T} <: AbstractDebiasingScheme{T}

NEOCC optical astrometry debiasing scheme.

!!! reference
    See https://neo.ssa.esa.int.
"""
mutable struct NEOCCDebiasing{T} <: AbstractDebiasingScheme{T}
    bias::Vector{Tuple{T, T}}
    # Default constructor
    function NEOCCDebiasing(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
        return new{T}(debiasneocc(radec))
    end
end

# Override getid
getid(::NEOCCDebiasing) = "NEOCC"

# Override update!
function update!(d::NEOCCDebiasing{T},
    radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
    d.bias = debiasneocc(radec)
    return nothing
end

@doc raw"""
    debiasneocc(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}

Return the debiasing correction of each element of `radec` according
to the NEOCC.

!!! reference
    See https://neo.ssa.esa.int.
"""
function debiasneocc(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
    # ID
    id = isempty(radec[end].num) ? unpackdesig(radec[end].tmpdesig) : radec[end].num
    id = replace(id, " " => "")
    # HTTP query
    resp = HTTP.get(string(NEOCC_OBS_API_URL, id, ".rwo"))
    # Convert to String
    text = String(resp.body)
    # Parse lines
    lines = split(text, "\n")[8:(7+length(radec))]
    # Parse weights
    debiasra = parse.(Float64, getindex.(lines, Ref(88:93)))
    debiasdec = parse.(Float64, getindex.(lines, Ref(141:146)))

    return @. tuple(debiasra, debiasdec)
end

@doc raw"""
    NEODyS2Debiasing{T} <: AbstractDebiasingScheme{T}

NEODyS-2 optical astrometry debiasing scheme.

!!! reference
    See https://newton.spacedys.com/neodys/.
"""
mutable struct NEODyS2Debiasing{T} <: AbstractDebiasingScheme{T}
    bias::Vector{Tuple{T, T}}
    # Default constructor
    function NEODyS2Debiasing(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
        return new{T}(debiasneodys2(radec))
    end
end

# Override getid
getid(::NEODyS2Debiasing) = "NEODyS-2"

# Override update!
function update!(d::NEODyS2Debiasing{T},
    radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
    d.bias = debiasneodys2(radec)
    return nothing
end

@doc raw"""
    debiasneodys2(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}

Return the debiasing correction of each element of `radec` according
to NEODyS-2.

!!! reference
    See https://newton.spacedys.com/neodys/.
"""
function debiasneodys2(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
    # ID
    id = isempty(radec[end].num) ? unpackdesig(radec[end].tmpdesig) : radec[end].num
    id = replace(id, " " => "")
    # HTTP query
    resp = HTTP.get(string(NEODyS2_OBS_API_URL, id, ".rwo"))
    # Convert to String
    text = String(resp.body)
    # Parse lines
    lines = split(text, "\n")[8:(7+length(radec))]
    # Parse weights
    debiasra = parse.(Float64, getindex.(lines, Ref(88:93)))
    debiasdec = parse.(Float64, getindex.(lines, Ref(141:146)))

    return @. tuple(debiasra, debiasdec)
end