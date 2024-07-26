@doc raw"""
    ObservatoryMPC{T <: AbstractFloat}

An observatory in MPC format.

# Fields

- `code::String`: observatory's three character identifier.
- `long::T`, `cos::T`,  `sin::T`:
    - For ground-based observatories: longitude [degrees east of Greenwich],
      `ρ*cos(ϕ')` and `ρ*sin(ϕ')`, where `ϕ'` is the geocentric latitude and
      `ρ` is the geocentric distance in earth radii.
    - For space-based observatories: `x`, `y` and `z` geocentric coordinates.
- `name::String`: observatory's name.
- `date::DateTime`: time of observation (for `:satellite`/`:occultation` observatories).
- `type::Symbol`: `:ground`, `:satellite` or `:occultation`.
- `units::Int`: whether position is encrypted in AU (1) or KM (2).

!!! reference
    The format is described in https://minorplanetcenter.net/iau/lists/ObsCodesF.html.
"""
@auto_hash_equals fields = (code,) struct ObservatoryMPC{T <: AbstractFloat}
    code::String
    long::T
    cos::T
    sin::T
    name::String
    date::DateTime
    type::Symbol
    units::Int
    # Inner constructor
    function ObservatoryMPC{T}(code::String, long::T, cos::T, sin::T, name::String, date::DateTime = DateTime(2000, 1, 1),
                               type::Symbol = :ground, units::Int = 0) where {T <: AbstractFloat}
        return new{T}(code, long, cos, sin, name, date, type, units)
    end
end

# Outer constructor
function ObservatoryMPC(code::String, long::T, cos::T, sin::T, name::String, date::DateTime = DateTime(2000, 1, 1),
                        type::Symbol = :ground, units::Int = 0) where {T <: AbstractFloat}
    return ObservatoryMPC{T}(code, long, cos, sin, name, date, type, units)
end

@doc raw"""
    unknownobs()

Return a `ObservatoryMPC` with no code, coordinates or name.
"""
unknownobs() = ObservatoryMPC("", NaN, NaN, NaN, "")

@doc raw"""
    isunknown(m::ObservatoryMPC{T}) where {T <: AbstractFloat}

Check whether `m` equals `unknownobs()`.
"""
function isunknown(m::ObservatoryMPC{T}) where {T <: AbstractFloat}
    return m == unknownobs()
end

isground(m::ObservatoryMPC{T}) where {T <: AbstractFloat} = m.type == :ground
issatellite(m::ObservatoryMPC{T}) where {T <: AbstractFloat} = m.type == :satellite
isoccultation(m::ObservatoryMPC{T}) where {T <: AbstractFloat} = m.type == :occultation

@doc raw"""
    hascoord(m::ObservatoryMPC{T}) where {T <: AbstractFloat}

Check whether `m` has non `NaN` coordinates.
"""
function hascoord(m::ObservatoryMPC{T}) where {T <: AbstractFloat}
    return !isnan(m.long) && !isnan(m.cos) && !isnan(m.sin)
end

# Print method for ObservatoryMPC
# Examples:
# Unknown observatory
# Spitzer Space Telescope [245]
# Greenwich [000] long: 0.0 cos: 0.62411 sin: 0.77873
function show(io::IO, m::ObservatoryMPC{T}) where {T <: AbstractFloat}
    if isunknown(m)
        print(io, "Unknown observatory")
    elseif isground(m)
        print(io, m.name, " [", m.code, "] long: ", m.long, " cos: ", m.cos, " sin: ", m.sin)
    else
        print(io, m.name, " [", m.code, "]")
    end
end

function neoparse(x::RegexMatch, i::Int, ::Type{Float64})
    y = tryparse(Float64, replace(x[i], " " => ""))
    if isnothing(y)
        return NaN
    else
        return y
    end
end

@doc raw"""
    ObservatoryMPC(m::RegexMatch)

Convert a match of `NEOs.OBSERVATORY_MPC_REGEX` to `ObservatoryMPC`.
"""
function ObservatoryMPC(m::RegexMatch)
    # Check that matched regex is correct
    @assert m.regex == OBSERVATORY_MPC_REGEX "Only matches of `NEOs.OBSERVATORY_MPC_REGEX` can be converted to `ObservatoryMPC`."
    # Field types
    types = fieldtypes(ObservatoryMPC{Float64})
    # ObservatoryMPC{Float64} fields
    args = map(i -> neoparse(m, i, types[i]), 1:5)

    if isnan(args[2]) && isnan(args[3]) && isnan(args[4])
        return ObservatoryMPC(args..., DateTime(2000, 1, 1), :satellite, 0)
    else
        return ObservatoryMPC(args...)
    end
end

@doc raw"""
    read_observatories_mpc(s::String)

Return the matches of `NEOs.OBSERVATORY_MPC_REGEX` in `s` as `Vector{ObservatoryMPC{Float64}}`.
`s` can be either a filename or a text.
"""
function read_observatories_mpc(s::String)
    if !contains(s, "\n") && isfile(s)
        # Read MPC formatted file
        s = String(read(s))
    end
    # Remove header
    s = replace(s, OBSERVATORIES_MPC_HEADER => "")
    # Vector of MPC observatories
    obs = Vector{ObservatoryMPC{Float64}}(undef, 0)
    # Iterate over the matches
    for m in eachmatch(OBSERVATORY_MPC_REGEX, s)
        push!(obs, ObservatoryMPC(m))
    end
    # Eliminate repeated entries
    unique!(obs)

    return obs
end

@doc raw"""
    mpc_long_str(x::T) where {T <: AbstractFloat}

Convert `x` to a string according to the `long` field in MPC format.
"""
function mpc_long_str(x::T) where {T <: AbstractFloat}
    # NaN => empty string
    if isnan(x)
        long_s = repeat(" ", 10)
    else
        long_s = @sprintf("%3.5f", x)
        long_s = lpad(long_s, 10)
    end
    return long_s
end

@doc raw"""
    mpc_cos_str(x::T) where {T <: AbstractFloat}

Convert `x` to a string according to the `cos` field in MPC format.
"""
function mpc_cos_str(x::T) where {T <: AbstractFloat}
    # NaN => empty string
    if isnan(x)
        cos_s = repeat(" ", 8)
    else
        cos_s = @sprintf("%1.6f", x)
    end
    return cos_s
end

@doc raw"""
    mpc_sin_str(x::T) where {T <: AbstractFloat}

Convert `x` to a string according to the `sin` field in MPC format.
"""
function mpc_sin_str(x::T) where {T <: AbstractFloat}
    # NaN => empty string
    if isnan(x)
        sin_s = repeat(" ", 9)
    else
        sin_s = @sprintf("%1.6f", abs(x))
        if x ≥ 0
            sin_s = string("+", sin_s)
        else
            sin_s = string("-", sin_s)
        end
    end

    return sin_s
end

# Convert `obs` to a string according to MPC format.
function string(obs::ObservatoryMPC{T}) where {T <: AbstractFloat}
    if isunknown(obs)
        return ""
    else
        # Longitude string
        long_s = mpc_long_str(obs.long)
        # Cosine string
        cos_s = mpc_cos_str(obs.cos)
        # Sine string  string
        sin_s = mpc_sin_str(obs.sin)
        # Join everything
        obs_s = string(obs.code, long_s, cos_s, sin_s, obs.name)
    end

    return obs_s
end

@doc raw"""
    write_observatories_mpc(obs::Vector{ObservatoryMPC{T}}, filename::String) where {T <: AbstractFloat}

Write `obs` to `filename` in MPC format.
"""
function write_observatories_mpc(obs::Vector{ObservatoryMPC{T}}, filename::String) where {T <: AbstractFloat}
    open(filename, "w") do file
        # Header
        write(file, OBSERVATORIES_MPC_HEADER, "\n")
        # Write observatories
        for i in eachindex(obs)
            line = string(obs[i])
            write(file, line, "\n")
        end
    end
end

@doc raw"""
    update_observatories_mpc()

Update the local observatories file.
"""
function update_observatories_mpc()
    # Download and read observatories file
    ObsCodes_path, txt = download_scratch(OBSERVATORIES_MPC_URL, "ObsCodes.txt")
    # Parse observatories
    obs = read_observatories_mpc(txt)
    # Write observatories to local file
    write_observatories_mpc(obs, ObsCodes_path)
    # Update global variable
    global OBSERVATORIES_MPC[] = read_observatories_mpc(ObsCodes_path)

    return nothing
end

@doc raw"""
    search_obs_code(obscode::String)

Return the observatory in `NEOs.OBSERVATORIES_MPC` that matches `obscode`.
"""
function search_obs_code(obscode::String)

    # Find indexes in OBSERVATORIES_MPC that match obscode
    idxs = findall(x -> x.code == obscode, OBSERVATORIES_MPC[])
    L_i = length(idxs)

    # No observatory matches obscode
    if L_i == 0
        observatory = unknownobs()
    # At least one observatory matches obscode
    else
        observatory = OBSERVATORIES_MPC[][idxs[1]]
        # More than one observatory matches obscode
        if L_i > 1
            @warn("""More than one observatory $(OBSERVATORIES_MPC[][idxs]) has code $obscode,
            selecting first: $(observatory.name)""")
        end
    end

    return observatory

end