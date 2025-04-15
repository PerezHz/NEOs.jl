@doc raw"""
    AbstractErrorModel{T <: Real}

Supertype for the optical astrometry error models API.
"""
abstract type AbstractErrorModel{T <: Real} end

# Optical astrometry weighting schemes API
# All weighting schemes `W{T}` must:
# 1. be mutable subtypes of `AbstractWeightingScheme{T}`,
# 2. have a `w8s::Vector{Tuple{T, T}}` field,
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
    w8s::Vector{Tuple{T, T}}
    # Default constructor
    function UniformWeights(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
        w8s = [(one(T), one(T)) for _ in eachindex(radec)]
        return new{T}(w8s)
    end
end

# Override getid
getid(::UniformWeights) = "Uniform"

# Override update!
function update!(w::UniformWeights{T},
    radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
    w.w8s = [(one(T), one(T)) for _ in eachindex(radec)]
    return nothing
end

@doc raw"""
    Veres17{T} <: AbstractWeightingScheme{T}

Veres et al. (2017) optical astrometry weighting scheme.

!!! reference
    See https://doi.org/10.1016/j.icarus.2017.05.021.
"""
mutable struct Veres17{T} <: AbstractWeightingScheme{T}
    w8s::Vector{Tuple{T, T}}
    # Default constructor
    function Veres17(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
        return new{T}(w8sveres17(radec))
    end
end

# Override getid
getid(::Veres17) = "Veres et al. (2017)"

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
    return @. tuple(1 / (rex * σs) ^ 2, 1 / (rex * σs) ^ 2)
end

@doc raw"""
    ADESWeights{T} <: AbstractWeightingScheme{T}

ADES optical astrometry weighting scheme.

!!! reference
    See https://minorplanetcenter.net/mpcops/documentation/ades/.
"""
mutable struct ADESWeights{T} <: AbstractWeightingScheme{T}
    w8s::Vector{Tuple{T, T}}
    # Default constructor
    function ADESWeights(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
        return new{T}(w8sades(radec))
    end
end

# Override getid
getid(::ADESWeights) = "ADES"

# Override update!
function update!(w::ADESWeights{T},
    radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
    w.w8s = w8sades(radec)
    return nothing
end

@doc raw"""
    w8sades(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}

Return the statistical weight of each element of `radec` according
to the ADES format. If ADES does not contain the necessary information,
return the weight assigned by Veres et al. (2017).

!!! reference
    See https://minorplanetcenter.net/mpcops/documentation/ades/.
"""
function w8sades(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
    # ID
    id = isempty(radec[end].num) ? unpackdesig(radec[end].tmpdesig) : radec[end].num
    # HTTP parameters
    params = JSON.json(Dict("desigs" => [id], "output_format" => ["ADES_DF"]))
    resp = get(MPC_OBS_API_URL, ("Content-Type" => "application/json",), params)
    # Convert to String
    text = String(resp.body)
    # Parse JSON
    dict = JSON.parse(text)
    # Parse weights
    obs = dict[1]["ADES_DF"]
    w8sra = Vector{T}(undef, length(obs))
    w8sdec = Vector{T}(undef, length(obs))
    for i in eachindex(obs)
        o = obs[i]
        if haskey(o, "rmsra") && !isnothing(o["rmsra"])
            w8sra[i] = parse(T, o["rmsra"])
        else
            w8sra[i] = σsveres17(radec[i])
        end
        if haskey(o, "rmsdec") && !isnothing(o["rmsdec"])
            w8sdec[i] = parse(T, o["rmsdec"])
        else
            w8sdec[i] = σsveres17(radec[i])
        end
    end

    return @. tuple(1 / w8sra^2, 1 / w8sdec^2)
end

@doc raw"""
    NEOCCWeights{T} <: AbstractWeightingScheme{T}

NEOCC optical astrometry weighting scheme.

!!! reference
    See https://neo.ssa.esa.int.
"""
mutable struct NEOCCWeights{T} <: AbstractWeightingScheme{T}
    w8s::Vector{Tuple{T, T}}
    # Default constructor
    function NEOCCWeights(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
        return new{T}(w8sneocc(radec))
    end
end

# Override getid
getid(::NEOCCWeights) = "NEOCC"

# Override update!
function update!(w::NEOCCWeights{T},
    radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
    w.w8s = w8sneocc(radec)
    return nothing
end

@doc raw"""
    w8sneocc(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}

Return the statistical weight of each element of `radec` according
to the NEOCC.

!!! reference
    See https://neo.ssa.esa.int.
"""
function w8sneocc(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
    # ID
    id = isempty(radec[end].num) ? unpackdesig(radec[end].tmpdesig) : radec[end].num
    id = replace(id, " " => "")
    # HTTP query
    resp = get(string(NEOCC_OBS_API_URL, id, ".rwo"))
    # Convert to String
    text = String(resp.body)
    # Parse lines
    lines = split(text, "\n")[8:(7+length(radec))]
    # Parse weights
    w8sra = parse.(Float64, getindex.(lines, Ref(78:82)))
    w8sdec = parse.(Float64, getindex.(lines, Ref(131:135)))

    return @. tuple(1 / w8sra^2, 1 / w8sdec^2)
end

@doc raw"""
    NEODyS2Weights{T} <: AbstractWeightingScheme{T}

NEODyS-2 optical astrometry weighting scheme.

!!! reference
    See https://newton.spacedys.com/neodys/.
"""
mutable struct NEODyS2Weights{T} <: AbstractWeightingScheme{T}
    w8s::Vector{Tuple{T, T}}
    # Default constructor
    function NEODyS2Weights(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
        return new{T}(w8sneodys2(radec))
    end
end

# Override getid
getid(::NEODyS2Weights) = "NEODyS-2"

# Override update!
function update!(w::NEODyS2Weights{T},
    radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
    w.w8s = w8sneodys2(radec)
    return nothing
end

@doc raw"""
    w8sneodys2(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}

Return the statistical weight of each element of `radec` according
to NEODyS-2.

!!! reference
    See https://newton.spacedys.com/neodys/.
"""
function w8sneodys2(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
    # ID
    id = isempty(radec[end].num) ? unpackdesig(radec[end].tmpdesig) : radec[end].num
    id = replace(id, " " => "")
    # HTTP query
    resp = get(string(NEODyS2_OBS_API_URL, id, ".rwo"))
    # Convert to Strings
    text = String(resp.body)
    # Parse lines
    lines = split(text, "\n")[8:(7+length(radec))]
    # Parse weights
    w8sra = parse.(Float64, getindex.(lines, Ref(78:82)))
    w8sdec = parse.(Float64, getindex.(lines, Ref(131:135)))

    return @. tuple(1 / w8sra^2, 1 / w8sdec^2)
end