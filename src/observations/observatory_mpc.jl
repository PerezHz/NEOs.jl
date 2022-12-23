@doc raw"""
    ObservatoryMPC{T <: AbstractFloat} 

An optical observatory in MPC format. The format is described in https://minorplanetcenter.net/iau/lists/ObsCodesF.html.

# Fields 

- `code::String`: observatory's identifier. 
- `long::T`: longitude [degrees east of Greenwich].
- `cos::T`: `ρ*cos(ϕ')`,
- `sin::T`: `ρ*sin(ϕ')`, where `ϕ'` is the geocentric latitude and `ρ` is the geocentric distance in earth radii.
- `name::String`: observatory's name.
"""
struct ObservatoryMPC{T <: AbstractFloat}
    code::String
    long::T
    cos::T
    sin::T
    name::String
    # Inner constructor
    function ObservatoryMPC{T}(code::String, long::T, cos::T, sin::T, name::String) where {T <: AbstractFloat}
        return new{T}(code, long, cos, sin, name)
    end
end

# Outer constructor
function ObservatoryMPC(code::String, long::T, cos::T, sin::T, name::String) where {T <: AbstractFloat}
    return ObservatoryMPC{T}(code, long, cos, sin, name)
end

# Two ObservatoryMPC are equal if ther code, long, cos, sin and name are equal
function hash(a::ObservatoryMPC{T}, h::UInt) where {T <: AbstractFloat}
    return hash((a.code, a.long, a.cos, a.sin, a.name), h)
end

function ==(a::ObservatoryMPC{T}, b::ObservatoryMPC{T}) where {T <: AbstractFloat}
    return hash(a) == hash(b)
end

# Print method for ObservatoryMPC
# Examples: 
# Unknown observatory
# Spitzer Space Telescope [245]
# Greenwich [000] long: 0.0 cos: 0.62411 sin: 0.77873
function show(io::IO, m::ObservatoryMPC{T}) where {T <: AbstractFloat} 
    if isunknown(m)
        print(io, "Unknown observatory")
    elseif hascoord(m)
        print(io, m.name, " [", m.code, "] long: ", m.long, " cos: ", m.cos, " sin: ", m.sin)
    else
        print(io, m.name, " [", m.code, "]")
    end
end

# MPC observatories file url 
const mpc_observatories_url = "https://minorplanetcenter.net/iau/lists/ObsCodes.html"

# Regular expression to parse an observatory in MPC format
const mpc_observatory_regex = Regex(join(
    [
        # Code regex + space (columns 1-4)
        raw"(?P<code>.{3}) ",
        # Longitude regex (columns 5-13)
        raw"(?P<long>.{9})",
        # Cosine regex + space (column 14-21)
        raw"(?P<cos>.{8})",
        # Sine regex (column 22-30)
        raw"(?P<sin>.{9})",
        # Name regex (columns 31-80)
        raw"(?P<name>.*)",
    ]
))

@doc raw"""
    parse_observatory_float(x::String)

Parses `x`; if the result is not a `Float64`, returns `NaN`. 
"""
function parse_observatory_float(x::String)
    # Try parsing x 
    p = tryparse(Float64, x)
    # The result is not a float 
    if isnothing(p)
        return NaN
    # The result is a float 
    else
        return p
    end
end

@doc raw"""
    hascoord(m::ObservatoryMPC{T}) where {T <: AbstractFloat}

Checks whether `m` has non `NaN` coordinates. 
"""
function hascoord(m::ObservatoryMPC{T}) where {T <: AbstractFloat}
    return !isnan(m.long) && !isnan(m.cos) && !isnan(m.sin)
end

@doc raw"""
    unknownobs()

Returns a `ObservatoryMPC` with no code, coordinates or name. 
"""
function unknownobs()
    return ObservatoryMPC("", NaN, NaN, NaN, "")
end

@doc raw"""
    isunknown(m::ObservatoryMPC{T}) where {T <: AbstractFloat}

Checks whether `m` equals `unknownobs()`.
"""
function isunknown(m::ObservatoryMPC{T}) where {T <: AbstractFloat}
    return m == unknownobs()
end 

@doc raw"""
    ObservatoryMPC(m::RegexMatch)

Converts a match of `NEOs.mpc_observatory_regex` to `ObservatoryMPC`.
"""
function ObservatoryMPC(m::RegexMatch)
    
    long_ = parse_observatory_float(string(m["long"]))
    cos_ = parse_observatory_float(string(m["cos"]))
    sin_ = parse_observatory_float(string(m["sin"]))

    return ObservatoryMPC(
        string(m["code"]),
        long_,
        cos_,
        sin_,
        string(m["name"]),
    )

end

@doc raw"""
    read_observatories_mpc(filename::String)

Returns the matches of `NEOs.mpc_observatory_regex` in `filename` as `ObservatoryMPC`.
"""
function read_observatories_mpc(filename::String)
    # Read lines of mpc formatted file (except header)
    lines = readlines(filename)[2:end]
    # Apply regular expressions
    matches = match.(mpc_observatory_regex, lines)
    # Eliminate nothings
    filter!(!isnothing, matches)
    # Convert matches to ObservatoryMPC
    obs = ObservatoryMPC.(matches)
    
    return obs
end

# Header of MPC observatories file 
const mpc_observatories_header = "Code  Long.   cos      sin    Name"

@doc raw"""
    get_raw_html(url::String)

Returns the raw html text of webpage `url`.
"""
function get_raw_html(url::String)
    # Get raw html 
    resp = get(url)
    # Convert to string 
    text = String(resp.body)
    
    return text
end

@doc raw"""
    parse_observatories_mpc(text::String)

Returns de matches of `NEOs.mpc_observatory_regex` in `text` as `ObservatoryMPC`.
"""
function parse_observatories_mpc(text::String)
    # Eliminate observatories file header 
    text = replace(text, mpc_observatories_header => "")
    # Vector of MPC observatories 
    obs = Vector{ObservatoryMPC{Float64}}(undef, 0)
    # Iterate over the matches 
    for m in eachmatch(mpc_observatory_regex, text)
        push!(obs, ObservatoryMPC(m))
    end
    
    return obs
end

# Path to mpc observatories file 
const ObsCodes_path = joinpath(observations_path, "ObsCodes.txt")
# List of mpc observatories
const mpc_observatories = Ref{Vector{ObservatoryMPC{Float64}}}(read_observatories_mpc(ObsCodes_path))

@doc raw"""
    mpc_long_str(x::T) where {T <: AbstractFloat}

Converts `x` to a string according to the `long` field in MPC format.
"""
function mpc_long_str(x::T) where {T <: AbstractFloat}
    # NaN => empty string 
    if isnan(x)
        long_s = repeat(" ", 9)
    else 
        long_s = @sprintf("%3.5f", x)
        long_s = lpad(long_s, 9)
    end
    return long_s
end

@doc raw"""
    mpc_cos_str(x::T) where {T <: AbstractFloat}

Converts `x` to a string according to the `cos` field in MPC format.
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

Converts `x` to a string according to the `sin` field in MPC format.
"""
function mpc_sin_str(x::T) where {T <: AbstractFloat}
    # NaN => empty string
    if isnan(x)
        sin_s = repeat(" ", 9)
    else
        sin_s = @sprintf("%1.6f", abs(x))
        if x ≥ 0
            sin_s = join(["+", sin_s])
        else
            sin_s = join(["-", sin_s])
        end
    end

    return sin_s 
end

@doc raw"""
    mpc_observatory_str(obs::ObservatoryMPC{T}) where {T <: AbstractFloat}

Converts `obs` to a string acoording to MPC format.
"""
function mpc_observatory_str(obs::ObservatoryMPC{T}) where {T <: AbstractFloat}
    # Longitude string 
    long_s = mpc_long_str(obs.long)
    # Cosine string 
    cos_s = mpc_cos_str(obs.cos)
    # Sine string  string 
    sin_s = mpc_sin_str(obs.sin)
    # Join everything
    obs_s = join([
        obs.code,
        " ",
        long_s,
        cos_s,
        sin_s,
        obs.name,
        "\n"
    ])

    return obs_s
end

@doc raw"""
    write_observatories_mpc(obs::Vector{ObservatoryMPC{T}}, filename::String) where {T <: AbstractFloat}

Writes `obs` to `filename` in MPC format. 
"""
function write_observatories_mpc(obs::Vector{ObservatoryMPC{T}}, filename::String) where {T <: AbstractFloat}
    open(filename, "w") do file
        # Header 
        write(file, mpc_observatories_header, "\n")
        # Write observatories 
        for i in eachindex(obs)
            line = mpc_observatory_str(obs[i])
            write(file, line)
        end 
    end
end

@doc raw"""
    update_observatories_mpc()

Updates the local observatories file.
"""
function update_observatories_mpc()
    # Download source file 
    @info "Downloading file $mpc_observatories_url"
    txt = get_raw_html(mpc_observatories_url)
    # Parse observatories 
    obs = parse_observatories_mpc(txt)
    m_before = length(mpc_observatories[])
    m_after = length(obs)
    @info "Found $m_after observatories ($m_before in the previous version of the file)"
    # Write observatories to local file 
    @info "Updating file $ObsCodes_path"
    write_observatories_mpc(obs, ObsCodes_path)
    # Update global variable 
    @info "Updating variable NEOs.mpc_observatories[]"
    global mpc_observatories[] = read_observatories_mpc(ObsCodes_path)

    return 
end

@doc raw"""
    search_obs_code(obscode::String)

Returns the observatory in `NEOs.mpc_observatories` that matches `obscode`.
"""
function search_obs_code(obscode::String)
    
    # Find indexes in mpc_observatories that match obscode
    idxs = findall(x -> x.code == obscode, mpc_observatories[])
    L_i = length(idxs)

    # No observatory matches obscode
    if L_i == 0
        @warn "Unknown observatory code $obscode"
        observatory = unknownobs()
    # At least one observatory matches obscode
    else
        observatory = mpc_observatories[][idxs[1]]
        # More than one observatory matches obscode
        if L_i > 1
            @warn("""More than one observatory $(mpc_observatories[][idxs]) has code $obscode, 
            selecting first: $(observatory.name)""")
        end
    end
    
    return observatory 
    
end