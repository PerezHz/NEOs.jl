"""
    CatalogueMPC <: AbstractCatalogue

A star catalogue recognized by the Minor Planet Center.

# Fields

- `code::Char`: single-character catalogue code.
- `value::String`: catalogue short name.
- `meaning::String`: catalogue full name.
- `deprecated::Bool`: deprecation flag.

!!! reference
    The list of catalogue codes is available at:
    - https://minorplanetcenter.net/iau/info/CatalogueCodes.html
    The list of catalogue names and deprecation flags is available at:
    - https://minorplanetcenter.net/mpcops/documentation/valid-ades-values/#astCat
"""
@auto_hash_equals fields = (code,) struct CatalogueMPC <: AbstractCatalogue
    code::Char
    value::String
    meaning::String
    deprecated::Bool
end

isdeprecated(c::CatalogueMPC) = c.deprecated
isunknown(c::CatalogueMPC) = c.code == ' '
unknowncat() = CatalogueMPC(' ', "UNK", "UNKNOWN", true)

# Print method for CatalogueMPC
show(io::IO, c::CatalogueMPC) = isunknown(c) ? print(io, "Unknown catalogue") :
    print(io, c.meaning, " [", c.code, "]")

# Internal constructor
function CatalogueMPC(dict)
    value = astrometryparse(String, dict["Value"])
    code = astrometryparse(Char, CATALOGUE_MPC_NAMES_TO_CODES[value])
    meaning = astrometryparse(String, dict["Meaning"])
    deprecated = astrometryparse(Bool, dict["Deprecated"])

    return CatalogueMPC(code, value, meaning, deprecated)
end

function parse_catalogues_mpc(text::AbstractString)
    # Parse JSON
    dict = JSON.parse(text)
    # Eliminate ZZCAT catalogue, as it does not have a code
    # See https://github.com/IAU-ADES/ADES-Master/issues/57
    filter!(d -> d["Value"] != "ZZCAT", dict)
    # Parse catalogues
    cats = CatalogueMPC.(dict)
    # Eliminate repeated entries
    unique!(cats)
    # Sort by single-character code
    sort!(cats, by = x -> x.code)

    return cats
end

function read_catalogues_mpc(filename::AbstractString)
    # Read file
    text = read(filename, String)
    # Parse catalogues
    cats = parse_catalogues_mpc(text)

    return cats
end

"""
    update_catalogues_mpc()

Update the local catalogues list.
"""
function update_catalogues_mpc()
    # Download and parse catalogues file
    text = fetch_http_text(MPC; mode = 0)
    # Parse catalogues
    cats = parse_catalogues_mpc(text)
    # Write catalogues to local file
    path = joinpath(SCRATCH_PATH[], "astCat_photCat.json")
    open(path, "w") do file
        write(file, text)
    end
    # Update global variable
    global CATALOGUES_MPC[] = cats

    return nothing
end

"""
    search_catalogue_code(code::Char)

Return the catalogue in `NEOs.CATALOGUES_MPC` that matches `code`.
"""
function search_catalogue_code(code::Char)
    i = findfirst(x -> x.code == code, CATALOGUES_MPC[])
    return isnothing(i) ? unknowncat() : CATALOGUES_MPC[][i]
end

"""
    search_catalogue_value(value::AbstractString)

Return the catalogue in `NEOs.CATALOGUES_MPC` that matches `value`.
"""
function search_catalogue_value(value::AbstractString)
    i = findfirst(x -> x.value == value, CATALOGUES_MPC[])
    return isnothing(i) ? unknowncat() : CATALOGUES_MPC[][i]
end