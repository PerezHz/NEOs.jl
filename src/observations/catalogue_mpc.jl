@doc raw"""
    CatalogueMPC

An astrometric reference catalogue in MPC format. The format is described in https://minorplanetcenter.net/iau/info/CatalogueCodes.html.

# Fields

- `code::String`: catalogue's identifier.
- `name::String`: catalogue's name.
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

# Regular expression to parse a catalogue in MPC format
const mpc_catalogue_regex = r"\s{2}(?P<code>\w{1})\s{4}(?P<name>.*)"

@doc raw"""
    CatalogueMPC(m::RegexMatch)

Convert a match of `NEOs.mpc_catalogue_regex` to `CatalogueMPC`.
"""
function CatalogueMPC(m::RegexMatch)
    @assert m.regex == mpc_catalogue_regex "Only matches of `NEOs.mpc_catalogue_regex` can be converted to `CatalogueMPC`."
    return CatalogueMPC(string(m["code"]), string(m["name"]))
end

@doc raw"""
    read_catalogues_mpc(filename::String)

Return the matches of `NEOs.mpc_catalogue_regex` in `filename` as `CatalogueMPC`.
"""
function read_catalogues_mpc(filename::String)
    # Check that file exists
    @assert isfile(filename) "Invalid filename"
    # Read lines of mpc formatted file (except header)
    lines = readlines(filename)[2:end]
    # Apply regular expressions
    matches = match.(mpc_catalogue_regex, lines)
    # Eliminate nothings
    filter!(!isnothing, matches)
    # Convert matches to CatalogueMPC
    cats = CatalogueMPC.(matches)
    # Eliminate repeated entries
    unique!(cats)

    return cats
end

@doc raw"""
    get_raw_html(url::String)

Return the raw html text of webpage `url`.
"""
function get_raw_html(url::String)
    # Get raw html
    resp = get(url)
    # Convert to string
    text = String(resp.body)

    return text
end

# Header of MPC catalogues file
const mpc_catalogues_header = "Char   Catalogue"

@doc raw"""
    parse_catalogues_mpc(text::String)

Return de matches of `NEOs.mpc_catalogue_regex` in `text` as `CatalogueMPC`.
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

    # Eliminate repeated entries
    unique!(cats)

    return cats
end

@doc raw"""
    mpc_catalogue_str(cat::CatalogueMPC)

Convert `cat` to a string according to MPC format.
"""
function mpc_catalogue_str(cat::CatalogueMPC)
    if isunknown(cat)
        return ""
    else
        # Code string
        code_s = join(["  ", cat.code, "    "])
        # Join everything
        cat_s = join([code_s, cat.name])

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
        write(file, mpc_catalogues_header, "\n")
        # Write observatories
        for i in eachindex(cats)
            line = mpc_catalogue_str(cats[i])
            write(file, line, "\n")
        end
    end
end

const downloader = Downloads.Downloader()
downloader.easy_hook = (easy, info) -> Downloads.Curl.setopt(easy, Downloads.Curl.CURLOPT_LOW_SPEED_TIME, 60)

@doc raw"""
    download_scratch(url::String, filename::String)

Download `url` and save the output to NEOs scratch space as `filename`. Return the local path and the contents of the file as a `String`.
"""
function download_scratch(url::String, filename::String)
    # Local file
    path = joinpath(scratch_path[], filename)
    # Download function
    f = retry(delays = fill(0.1, 3)) do
        try
            download(url, path; downloader)
        catch
            rethrow()
        end
    end
    # Download source_file
    f()
    # Read local file
    txt = read(path, String)

    return path, txt
end

# List of MPC catalogues
const mpc_catalogues = Ref{Vector{CatalogueMPC}}([unknowncat()])

@doc raw"""
    update_catalogues_mpc()

Update the local catalogues file.
"""
function update_catalogues_mpc()
    # Download and read catalogues file
    CatalogueCodes_path, txt = download_scratch(mpc_catalogues_url, "CatalogueCodes.txt")
    # Parse catalogues
    cats = parse_catalogues_mpc(txt)
    # Write catalogues to local file
    write_catalogues_mpc(cats, CatalogueCodes_path)
    # Update global variable
    global mpc_catalogues[] = read_catalogues_mpc(CatalogueCodes_path)

    return nothing
end

@doc raw"""
    search_cat_code(catcode::String)

Return the catalogue in `NEOs.mpc_catalogues` that matches `catcode`.
"""
function search_cat_code(catcode::String)

    # Find indexes in mpc_catalogues that match catcode
    idxs = findall(x -> x.code == catcode, mpc_catalogues[])
    L_i = length(idxs)

    # No catalog matches catcode
    if L_i == 0
        catalogue = unknowncat()
    # Exactly one catalogue matches catcode
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