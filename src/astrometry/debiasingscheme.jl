"""
    AbstractDebiasingScheme{T} <: AbstractAstrometryErrorModel{T}

Supertype for the optical astrometry debiasing schemes interface.

Every debiasing scheme `D{T}` must:
- be a mutable subtype of `AbstractDebiasingScheme{T}`,
- implement a `D(::AbstractOpticalVector{T})` constructor,
- override `debias(::D)`,
- override `getid(::D)`,
- override `update!(::D{T}, ::AbstractOpticalVector{T})`.
"""
abstract type AbstractDebiasingScheme{T} <: AbstractAstrometryErrorModel{T} end

# Print method for AbstractDebiasingScheme
show(io::IO, x::AbstractDebiasingScheme) = print(io, getid(x),
    " debiasing scheme with ", length(debias(x)), " observations")

"""
    ZeroDebiasing{T} <: AbstractDebiasingScheme{T}

Zero optical astrometry debiasing scheme.
"""
mutable struct ZeroDebiasing{T} <: AbstractDebiasingScheme{T}
    debias::Vector{NTuple{2, T}}
end

# Constructor
function ZeroDebiasing(optical::AbstractOpticalVector{T}) where {T <: Real}
    return ZeroDebiasing{T}([(zero(T), zero(T)) for _ in eachindex(optical)])
end

# Override debias
debias(x::ZeroDebiasing) = x.debias

# Override getid
getid(::ZeroDebiasing) = "Zero"

# Override update!
function update!(x::ZeroDebiasing{T}, optical::AbstractOpticalVector{T}) where {T <: Real}
    x.debias = [(zero(T), zero(T)) for _ in eachindex(optical)]
    return nothing
end

"""
    SourceDebiasing{T} <: AbstractDebiasingScheme{T}

Source optical astrometry debiasing scheme.
"""
mutable struct SourceDebiasing{T} <: AbstractDebiasingScheme{T}
    debias::Vector{NTuple{2, T}}
end

# Constructors
function SourceDebiasing(optical::AbstractOpticalVector{T}) where {T <: Real}
    return SourceDebiasing{T}(debias.(optical))
end

# Override debias
debias(x::SourceDebiasing) = x.debias

# Override getid
getid(::SourceDebiasing) = "Source"

# Override update!
function update!(x::SourceDebiasing{T}, optical::AbstractOpticalVector{T}) where {T <: Real}
    x.debias = debias.(optical)
    return nothing
end

"""
    Farnocchia15{T} <: AbstractDebiasingScheme{T}

Farnocchia et al. (2015) optical astrometry debiasing scheme.

!!! reference
    See:
    - https://doi.org/10.1016/j.icarus.2014.07.033
"""
mutable struct Farnocchia15{T} <: AbstractDebiasingScheme{T}
    debias::Vector{NTuple{2, T}}
    catcodes::Vector{Char}
    truth::String
    resol::Resolution
    table::Matrix{T}
end

# Constructor
function Farnocchia15(optical::AbstractOpticalVector{T}) where {T <: Real}
    catcodes, truth, resol, table = select_debiasing_table("2014")
    debias = debiasing(optical, catcodes, truth, resol, table)
    return Farnocchia15{T}(debias, catcodes, truth, resol, table)
end

# Override debias
debias(x::Farnocchia15) = x.debias

# Override getid
getid(::Farnocchia15) = "Farnocchia et al. (2015)"

# Override update!
function update!(x::Farnocchia15{T}, optical::AbstractOpticalVector{T}) where {T <: Real}
    x.debias = debiasing(optical, x.catcodes, x.truth, x.resol, x.table)
    return nothing
end

"""
    Eggl20{T} <: AbstractDebiasingScheme{T}

Eggl et al. (2020) optical astrometry debiasing scheme.

!!! reference
    See:
    - https://doi.org/10.1016/j.icarus.2019.113596
"""
mutable struct Eggl20{T} <: AbstractDebiasingScheme{T}
    debias::Vector{NTuple{2, T}}
    hires::Bool
    catcodes::Vector{Char}
    truth::String
    resol::Resolution
    table::Matrix{T}
end

# Default constructor
function Eggl20(optical::AbstractOpticalVector{T}, hires::Bool = false) where {T <: Real}
    id = hires ? "hires2018" : "2018"
    catcodes, truth, resol, table = select_debiasing_table(id)
    debias = debiasing(optical, catcodes, truth, resol, table)
    return Eggl20{T}(debias, hires, catcodes, truth, resol, table)
end

# Override debias
debias(x::Eggl20) = x.debias

# Override getid
getid(x::Eggl20) = string("Eggl et al. (2020)", x.hires ? " [high resolution]" : "")

# Override update!
function update!(x::Eggl20{T}, optical::AbstractOpticalVector{T}) where {T <: Real}
    x.debias = debiasing(optical, x.catcodes, x.truth, x.resol, x.table)
    return nothing
end

# Return the catalogue codes, truth catalogue, resolution and bias matrix of the
# corresponding debiasing scheme. The possible values for `id` are:
# - `2014` corresponds to https://doi.org/10.1016/j.icarus.2014.07.033,
# - `2018` corresponds to https://doi.org/10.1016/j.icarus.2019.113596 (default),
# - `hires2018` corresponds to https://doi.org/10.1016/j.icarus.2019.113596
#     (high resolution).
function select_debiasing_table(id::String = "2018")
    # Debiasing tables are loaded "lazily" via Julia artifacts,
    # according to rules in Artifacts.toml
    if id == "2018"
        debias_path = artifact"debias_2018"
        catcodes = CATALOGUE_MPC_CODES_2018
        # The healpix tesselation resolution of the bias map from
        # https://doi.org/10.1016/j.icarus.2019.113596
        NSIDE = 64
        # In 2018 debias table Gaia DR2 catalogue is regarded as the truth
        truth = "V"
    elseif id == "hires2018"
        debias_path = artifact"debias_hires2018"
        catcodes = CATALOGUE_MPC_CODES_2018
        # The healpix tesselation resolution of the high-resolution bias map from
        # https://doi.org/10.1016/j.icarus.2019.113596
        NSIDE = 256
        # In 2018 debias table Gaia DR2 catalogue is regarded as the truth
        truth = "V"
    elseif id == "2014"
        debias_path = artifact"debias_2014"
        catcodes = CATALOGUE_MPC_CODES_2014
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

# Return the total debiasing correction [arcsec] in both right ascension and declination,
# for each element of a vector of optical astrometry.
#
# # Arguments
#
# - `obs::AbstractOpticalVector{T}`: optical astrometry.
# - `catcodes::Vector{String}`: catalogues present in debiasing table.
# - `truth::String`: truth catalogue of debiasing table.
# - `resol::Resolution`: resolution.
# - `table::Matrix{T}`: debiasing table.
function debiasing(obs::AbstractOpticalVector{T}, catcodes::Vector{Char}, truth::String,
                   resol::Resolution, table::Matrix{T}) where {T <: Real}
    # Allocate memory
    debias = Vector{NTuple{2, T}}(undef, length(obs))
    warncodes = Vector{Int}(undef, length(obs))
    # Fill
    for i in eachindex(obs)
        debias[i], warncodes[i] = debiasing(obs[i], catcodes, truth, resol, table)
    end
    # Print warnings
    I1 = findall(==(1), warncodes)
    N1, C1 = length(I1), unique!(catalogue.(obs[I1]))
    !isempty(I1) && @warn "Catalogues $C1 not found in the debiasing table.\n\
        Setting debiasing corrections equal to zero in $N1 observations."
    I2 = findall(==(2), warncodes)
    N2, C2 = length(I2), unique!(catalogue.(obs[I2]))
    !isempty(I2) && @warn "Catalogues $C2 not available in the observation record.\n\
        Setting debiasing corrections equal to zero in $N2 observations."
    I3 = findall(==(3), warncodes)
    N3, C3 = length(I3), unique!(catalogue.(obs[I3]))
    !isempty(I3) && @warn "Catalogues $C3 do not correspond to an MPC catalogue code.\n\
        Setting debiasing corrections equal to zero in $N3 observations."

    return debias
end

function debiasing(obs::AbstractOpticalAstrometry{T}, catcodes::Vector{Char},
                   truth::String, resol::Resolution, table::Matrix{T}) where {T <: Real}

    # Catalogue code
    catcode = catalogue(obs).code
    # Warning code
    warncode = 0

    # If star catalogue is not present in debiasing table,
    # then set corrections equal to zero
    if (catcode ∉ catcodes) && catcode != "Y"
        # Catalogue exists in CATALOGUES_MPC[]
        if !isunknown(catalogue(obs))
            # Truth catalogue is not present in debiasing table
            # but it does not send a warning
            warncode = catcode == truth ? 0 : 1
        # Unknown catalogue
        elseif catcode == ""
            warncode = 2
        # Catalogue code is not empty but it does not match an MPC catalogue code either
        else
            warncode = 3
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
        pix_ind = ang2pixRing(resol, π/2 - dec(obs), ra(obs))

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
        et_secs_i = dtutc2et(date(obs))
        # Seconds sinde J2000 (TT)
        tt_secs_i = et_secs_i - ttmtdb(et_secs_i/daysec)
        # Years since J2000
        yrs_J2000_tt = tt_secs_i/(daysec*yr)
        # Total debiasing correction in right ascension (arcsec)
        α_corr = dRA + yrs_J2000_tt*pmRA/1_000
        # Total debiasing correction in declination (arcsec)
        δ_corr = dDEC + yrs_J2000_tt*pmDEC/1_000
    end

    return (α_corr, δ_corr), warncode
end