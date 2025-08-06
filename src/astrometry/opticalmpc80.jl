struct MPC80FileReader
    optical::Vector{SubString}
end

function is_mpc80_two_liner(line::AbstractString)
    line[15] == 'S' && return true
    obscode = view(line, 78:80)
    return obscode == "248" || obscode == "270"
end

function MPC80FileReader(text::AbstractString)
    # Parse optical observations
    if Sys.iswindows()
        # Windows uses \r\n as newline instead of \n
        text = replace(text, "\r\n" => '\n')
    end
    N = count(==('\n'), text)
    L = length(text)
    optical = Vector{SubString{String}}(undef, N)
    a = 1
    for i in eachindex(optical)
        b = findnext('\n', text, a)
        optical[i] = view(text, a:b)
        if is_mpc80_two_liner(optical[i])
            c = findnext('\n', text, b+1)
            optical[i] = view(text, a:c)
            a = c + 1
        else
            a = b + 1
        end
        if a > L
            i < N && deleteat!(optical, i+1:N)
            break
        end
    end

    return MPC80FileReader(optical)
end

"""
    OpticalMPC80{T} <: AbstractOpticalAstrometry{T}

An optical astrometric observation in Minor Planet Center 80-column format.

# Fields

- `number::String`: packed minor body number.
- `desig::String`: packed provisional or temporary designation.
- `discovery::Char`: discovery asterisk.
- `note1::Char`: publishable note or a program code.
- `note2::Char`: J2000.0 conversion or observation technique flag.
- `date::DateTime`: date of observation.
- `ra::T`: observed right ascension [rad].
- `dec::T`: observed declination [rad].
- `info1::String`: additional information.
- `mag::T`: observed magnitude.
- `band::Char`: magnitude band (or nuclear/total flag for comets).
- `catalogue::CatalogueMPC`: reference star catalogue.
- `info2::String`: additional information.
- `observatory::ObservatoryMPC{T}`: observatory.
- `source::String`: original record.

!!! reference
    The Minor Planet Center 80-column format is described at:
    - https://minorplanetcenter.net/iau/info/OpticalObs.html
    and discussed thoroughly in pages 158-181 of:
    - https://doi.org/10.1016/j.icarus.2010.06.003
"""
@auto_hash_equals fields = (date, ra, dec, observatory) struct OpticalMPC80{T} <: AbstractOpticalAstrometry{T}
    number::String
    desig::String
    discovery::Char
    note1::Char
    note2::Char
    date::DateTime
    ra::T
    dec::T
    info1::String
    mag::T
    band::Char
    catalogue::CatalogueMPC
    info2::String
    observatory::ObservatoryMPC{T}
    source::String
end

# AbstractAstrometryObservation interface
date(x::OpticalMPC80) = x.date
band(x::OpticalMPC80) = x.band
observatory(x::OpticalMPC80) = x.observatory
catalogue(x::OpticalMPC80) = x.catalogue
rms(::OpticalMPC80{T}) where {T <: Real} = (one(T), one(T))
debias(::OpticalMPC80{T}) where {T <: Real} = (zero(T), zero(T))
corr(::OpticalMPC80{T}) where {T <: Real} = zero(T)

isdiscovery(x::OpticalMPC80) = x.discovery == '*'

# Print method for OpticalMPC80
function show(io::IO, o::OpticalMPC80)
    # If there is no number, use temporary designation
    if isempty(o.number)
        id_str = 'I' ≤ o.desig[1] ≤ 'K' ? unpackdesig(o.desig) : o.desig
    else
        id_str = unpacknum(o.number)
    end

    print(io, id_str, " α: ", @sprintf("%.5f", rad2deg(o.ra)), "° δ: ", @sprintf("%.5f",
        rad2deg(o.dec)), "° t: ", o.date, " obs: ", o.observatory.name)
end

# Parsing

mpc80parse(_, ::Type{String}, x) = String(x)
mpc80parse(_, ::Type{Char}, x) = isempty(x) ? ' ' : first(x)
mpc80parse(_, ::Type{CatalogueMPC}, x) = isempty(x) ? unknowncat() : search_catalogue_code(first(x))
mpc80parse(_, ::Type{ObservatoryMPC{T}}, x) where {T} = search_observatory_code(x)

function mpc80parse(_, ::Type{DateTime}, x)
    date = DateTime(view(x, 1:10), "yyyy mm dd")
    fraction = parse(Float64, view(x, 11:length(x)))
    return date + Microsecond(round(Int, fraction * 8.64e10))
end

function mpc80parse(name, ::Type{T}, x) where {T <: Real}
    if name == :ra
        return hms2rad(x)
    elseif name == :dec
        return dms2rad(x)
    else
        return isempty(x) ? T(NaN) : parse(T, x)
    end
end

OpticalMPC80(r::DataFrameRow) = OpticalMPC80{Float64}(
    r.number, r.desig, r.discovery, r.note1, r.note2, r.date, r.ra, r.dec,
    r.info1, r.mag, r.band, r.catalogue, r.info2, r.observatory, r.source
)

function parse_optical_mpc80(text::AbstractString)
    # File reader
    reader = MPC80FileReader(text)
    L = length(reader.optical)
    # Construct DataFrame
    R = OpticalMPC80{Float64}
    names, types = fieldnames(R), fieldtypes(R)
    df = DataFrame([fill(astrometrydefault(fieldtype(R, name)), L) for name in names],
        collect(names))
    for (i, line) in enumerate(reader.optical)
        for (name, type, idxs) in zip(names, types, MPC80_OPTICAL_COLUMNS)
            x = strip(view(line, idxs))
            df[i, name] = mpc80parse(name, type, x)
        end
    end
    # Two line observations
    mask = findall(istwoliner, df.observatory)
    for i in mask
        obs = df.observatory[i]
        line = view(reader.optical[i], 82:162)
        I1, I2, I3 = isroving(obs) ? (35:44, 46:55, 57:61) : (35:45, 47:57, 59:69)
        if isoccultation(obs) || issatellite(obs)
            frame = line[33] == '1' ? "ICRF_KM" : "ICRF_AU"
        else
            frame = obs.frame
        end
        coords = SVector{3, Float64}(
            parse(Float64, replace(view(line, I1), " " => "")),
            parse(Float64, replace(view(line, I2), " " => "")),
            parse(Float64, replace(view(line, I3), " " => ""))
        )
        df.observatory[i] = ObservatoryMPC(obs; frame, coords)
    end
    # Source string
    df.source = reader.optical
    # Parse observations
    optical = OpticalMPC80.(eachrow(df))
    # Eliminate repeated entries
    unique!(optical)
    # Sort by date
    sort!(optical)

    return optical
end

"""
    fetch_optical_mpc80(id, source)

Return the optical astrometry of minor body `id` in the Minor Planet Center
80-column format. The `source` of the observations can be either `MPC` or
`NEOCP`.

!!! reference
    The Minor Planet Center observations APIs are described at:
    - https://minorplanetcenter.net/mpcops/documentation/observations-api/
    - https://minorplanetcenter.net/mpcops/documentation/neocp-observations-api/
"""
function fetch_optical_mpc80(id::AbstractString, ::Type{MPC})
    # Get and parse HTTP response
    text = fetch_http_text(MPC; mode = 2, id = id, format = "OBS80")
    # Parse JSON
    dict = JSON.parse(text)
    # Parse observations
    optical = parse_optical_mpc80(dict[1]["OBS80"])

    return optical
end

function fetch_optical_mpc80(id::AbstractString, ::Type{NEOCP})
    # Get and parse HTTP response
    text = fetch_http_text(NEOCP; mode = 1, id = id, format = "OBS80")
    # Parse JSON
    dict = JSON.parse(text)
    # Parse observations
    optical = parse_optical_mpc80(dict[1]["OBS80"])

    return optical
end

# Read / write

"""
    read_optical_mpc80(filename)

Read from `filename` a vector of optical astrometry in the Minor Planet
Center 80-column format.
"""
function read_optical_mpc80(filename::AbstractString)
    # Read file
    text = read(filename, String)
    # Parse observations
    optical = parse_optical_mpc80(text)

    return optical
end

"""
    write_optical_mpc80(obs, filename)

Write `obs` to `filename` in the Minor Planet Center 80-column format.
"""
function write_optical_mpc80(obs::AbstractVector{OpticalMPC80{T}},
                             filename::AbstractString) where {T <: Real}
    open(filename, "w") do file
        for i in eachindex(obs)
            write(file, obs[i].source)
        end
    end
end