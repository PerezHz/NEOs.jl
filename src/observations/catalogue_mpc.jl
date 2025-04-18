@doc raw"""
    CatalogueMPC

An astrometric reference catalogue in MPC format.

# Fields

- `code::String`: catalogue's single character identifier.
- `name::String`: catalogue's name.

!!! reference
    The format is described in https://minorplanetcenter.net/iau/info/CatalogueCodes.html.
"""
@auto_hash_equals struct CatalogueMPC
    code::String
    name::String
    # Inner constructor
    function CatalogueMPC(code::String, name::String)
        return new(code, name)
    end
end

@doc raw"""
    unknowncat()

Return a `CatalogueMPC` with no code or name.
"""
function unknowncat()
    return CatalogueMPC("", "")
end

@doc raw"""
    isunknown(m::CatalogueMPC)

Check whether `m` equals `unknowncat()`.
"""
function isunknown(m::CatalogueMPC)
    return m == unknowncat()
end

# Print method for CatalogueMPC
# Examples:
# Unknown catalogue
# USNO-A1.0 [a]
function show(io::IO, m::CatalogueMPC)
    if isunknown(m)
        print(io, "Unknown catalogue")
    else
        print(io, m.name, " [", m.code, "]")
    end
end

neoparse(x::RegexMatch, i::Int, ::Type{String}) = String(strip(x[i]))

@doc raw"""
    CatalogueMPC(m::RegexMatch)

Convert a match of `NEOs.CATALOGUE_MPC_REGEX` to `CatalogueMPC`.
"""
function CatalogueMPC(m::RegexMatch)
    # Check that matched regex is correct
    @assert m.regex == CATALOGUE_MPC_REGEX "Only matches of `NEOs.CATALOGUE_MPC_REGEX` can be converted to `CatalogueMPC`."
    # Field types
    types = fieldtypes(CatalogueMPC)
    # CatalogueMPC fields
    args = map(i -> neoparse(m, i, types[i]), 1:length(types))

    return CatalogueMPC(args...)
end

@doc raw"""
    read_catalogues_mpc(s::String)

Return the matches of `NEOs.CATALOGUE_MPC_REGEX` in `s` as `Vector{CatalogueMPC}`.
`s` can be either a filename or a text.
"""
function read_catalogues_mpc(s::String)
    if !contains(s, "\n") && isfile(s)
        # Read MPC formatted file
        s = String(read(s))
    end
    # Remove header
    s = replace(s, CATALOGUES_MPC_HEADER => "")
    # Vector of MPC catalogues
    cats = Vector{CatalogueMPC}(undef, 0)
    # Iterate over the matches
    for m in eachmatch(CATALOGUE_MPC_REGEX, s)
        push!(cats, CatalogueMPC(m))
    end
    # Eliminate repeated entries
    unique!(cats)

    return cats
end

# Convert `cat` to a string according to MPC format.
function string(cat::CatalogueMPC)
    if isunknown(cat)
        return ""
    else
        # Code string
        code_s = string("  ", cat.code, "    ")
        # Join everything
        cat_s = string(code_s, cat.name)

        return cat_s
    end
end

@doc raw"""
    write_catalogues_mpc(cats::Vector{CatalogueMPC}, filename::String)

Write `cats` to `filename` in MPC format.
"""
function write_catalogues_mpc(cats::Vector{CatalogueMPC}, filename::String)
    open(filename, "w") do file
        # Header
        write(file, CATALOGUES_MPC_HEADER, "\n")
        # Write observatories
        for i in eachindex(cats)
            line = string(cats[i])
            write(file, line, "\n")
        end
    end
end

@doc raw"""
    download_scratch(url::String, filename::String; connect_timeout = 180, readtimeout = 180)

Download `url` and save the output to NEOs scratch space as `filename`. Return the local
path and the contents of the file as a `String`.
"""
function download_scratch(url::String, filename::String; connect_timeout = 180, readtimeout = 180)
    # Local file
    path = joinpath(scratch_path[], filename)
    # Get raw html (HTTP.get retries four times by default)
    resp = HTTP.get(url; connect_timeout, readtimeout)
    # Read local file
    txt = String(resp.body)

    return path, txt
end

@doc raw"""
    update_catalogues_mpc()

Update the local catalogues file.
"""
function update_catalogues_mpc()
    # Download and read catalogues file
    CatalogueCodes_path, txt = download_scratch(CATALOGUES_MPC_URL, "CatalogueCodes.txt")
    # Parse catalogues
    cats = read_catalogues_mpc(txt)
    # Write catalogues to local file
    write_catalogues_mpc(cats, CatalogueCodes_path)
    # Update global variable
    global CATALOGUES_MPC[] = read_catalogues_mpc(CatalogueCodes_path)

    return nothing
end

@doc raw"""
    search_cat_code(catcode::String)

Return the catalogue in `NEOs.CATALOGUES_MPC` that matches `catcode`.
"""
function search_cat_code(catcode::String)

    # Find indexes in mpc_catalogues that match catcode
    idxs = findall(x -> x.code == catcode, CATALOGUES_MPC[])
    L_i = length(idxs)

    # No catalog matches catcode
    if L_i == 0
        catalogue = unknowncat()
    # Exactly one catalogue matches catcode
    else
        catalogue = CATALOGUES_MPC[][idxs[1]]
        # More than one catalogue matches catcode
        if L_i > 1
            @warn("""More than one catalogue $(CATALOGUES_MPC[][idxs]) has code $catcode,
            selecting first: $(catalogue.name)""")
        end
    end

    return catalogue

end