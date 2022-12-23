@doc raw"""
    CatalogueMPC

An astrometric reference catalogue in MPC format. The format is described in https://minorplanetcenter.net/iau/info/CatalogueCodes.html.

# Fields 

- `code::String`: catalogue's identifier. 
- `name::String`: catalogue's name.
"""
struct CatalogueMPC
    code::String
    name::String
    # Inner constructor
    function CatalogueMPC(code::String, name::String)
        return new(code, name)
    end
end

# Two CatalogueMPC are equal if ther code and name are equal
function hash(a::CatalogueMPC, h::UInt)
    return hash((a.code, a.name), h)
end

function ==(a::CatalogueMPC, b::CatalogueMPC)
    return hash(a) == hash(b)
end

# Print method for CatalogueMPC 
# Examples: 
# Unknown observatory
# Spitzer Space Telescope [245]
# Greenwich [000] long: 0.0 cos: 0.62411 sin: 0.77873
function show(io::IO, m::CatalogueMPC)
    if isunknown(m)
        print(io, "Unknown catalogue")
    else
        print(io, m.name, " [", m.code, "]")
    end
end

# MPC catalogues file url 
const mpc_catalogues_url = "https://minorplanetcenter.net/iau/info/CatalogueCodes.html"

# Regular expression to parse a catalogue in MPC format
const mpc_catalogue_regex = r"\s{2}(?P<code>\w{1})\s{4}(?P<name>.*)"

@doc raw"""
    unknowncat()

Returns a `CatalogueMPC` with no code or name. 
"""
function unknowncat()
    return CatalogueMPC("", "")
end

@doc raw"""
    isunknown(m::CatalogueMPC)

Checks whether `m` equals `unknowncat()`.
"""
function isunknown(m::CatalogueMPC)
    return m == unknowncat()
end 

@doc raw"""
    CatalogueMPC(m::RegexMatch)

Converts a match of `NEOs.mpc_catalogue_regex` to `CatalogueMPC`.
"""
function CatalogueMPC(m::RegexMatch)

    return CatalogueMPC(string(m["code"]), string(m["name"]))

end

@doc raw"""
    read_catalogues_mpc(filename::String)

Returns the matches of `NEOs.mpc_catalogue_regex` in `filename` as `CatalogueMPC`.
"""
function read_catalogues_mpc(filename::String)
    # Read lines of mpc formatted file (except header)
    lines = readlines(filename)[2:end]
    # Apply regular expressions
    matches = match.(mpc_catalogue_regex, lines)
    # Eliminate nothings
    filter!(!isnothing, matches)
    # Convert matches to CatalogueMPC
    cats = CatalogueMPC.(matches)
    
    return cats
end

# Header of MPC catalogues file 
const mpc_catalogues_header = "Char   Catalogue"

@doc raw"""
    parse_catalogues_mpc(text::String)

Returns de matches of `NEOs.mpc_catalogue_regex` in `text` as `CatalogueMPC`.
"""
function parse_catalogues_mpc(text::String)
    # Eliminate catalogues file header 
    text = replace(text, mpc_catalogues_header => "")
    # Vector of MPC catalogues  
    cats = Vector{CatalogueMPC}(undef, 0)
    # Iterate over the matches 
    for m in eachmatch(mpc_catalogue_regex, text)
        push!(cats, CatalogueMPC(m))
    end
    
    return cats 
end

# Path to MPC catalogues file 
const CatalogueCodes_path = joinpath(observations_path, "CatalogueCodes.txt")
# List of MPC catalogues 
const mpc_catalogues = Ref{Vector{CatalogueMPC}}(read_catalogues_mpc(CatalogueCodes_path))

@doc raw"""
    mpc_catalogue_str(cat::CatalogueMPC)

Converts `cats` to a string acoording to MPC format.
"""
function mpc_catalogue_str(cat::CatalogueMPC)
    # Code string 
    code_s = join(["  ", cat.code, "    "])
    # Join everything
    cat_s = join([
        code_s,
        cat.name,
        "\n"
    ])

    return cat_s
end

@doc raw"""
    write_catalogues_mpc(cats::Vector{CatalogueMPC}, filename::String)

Writes `cats` to `filename` in MPC format. 
"""
function write_catalogues_mpc(cats::Vector{CatalogueMPC}, filename::String)
    open(filename, "w") do file
        # Header 
        write(file, mpc_catalogues_header, "\n")
        # Write observatories 
        for i in eachindex(cats)
            line = mpc_catalogue_str(cats[i])
            write(file, line)
        end 
    end
end

@doc raw"""
    update_catalogues_mpc()

Updates the local catalogues file.
"""
function update_catalogues_mpc()
    # Download source file 
    @info "Downloading file $mpc_catalogues_url"
    txt = get_raw_html(mpc_catalogues_url)
    # Parse catalogues  
    cats = parse_catalogues_mpc(txt)
    m_before = length(mpc_catalogues[])
    m_after = length(cats)
    @info "Found $m_after catalogues ($m_before in the previous version of the file)"
    # Write catalogues to local file 
    @info "Updating file $CatalogueCodes_path"
    write_catalogues_mpc(cats, CatalogueCodes_path)
    # Update global variable 
    @info "Updating variable NEOs.mpc_catalogues[]"
    global mpc_catalogues[] = read_catalogues_mpc(CatalogueCodes_path)

    return 
end

@doc raw"""
    search_cat_code(catcode::String)

Returns the catalogue in `NEOs.mpc_cataloges` that matches `catcode`.
"""
function search_cat_code(catcode::String)
    
    # Find indexes in mpc_catalogues that match catcode
    idxs = findall(x -> x.code == catcode, mpc_catalogues[])
    L_i = length(idxs)

    # No catalog matches catcode
    if L_i == 0
        @warn "Unknown catalogue code $catcode"
        catalogue = unknowncat()
    # At least one catalogue matches catcode
    else
        catalogue = mpc_catalogues[][idxs[1]]
        # More than one catalogue matches catcode
        if L_i > 1
            @warn("""More than one catalogue $(mpc_catalogues[][idxs]) has code $catcode, 
            selecting first: $(catalogue.name)""")
        end
    end
    
    return catalogue 
    
end