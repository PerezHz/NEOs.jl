struct RWOFileReader
    source::Symbol
    header::SubString{String}
    optical::Vector{SubString{String}}
    radar::Vector{SubString{String}}
end

get_rwo_optical_shift(x::Symbol) = x == :NEOCC ? 3 : 0

get_rwo_optical_columns(x::Symbol) = x == :NEOCC ?
    NEOCC_OPTICAL_COLUMNS : NEODyS2_OPTICAL_COLUMNS

function is_rwo_two_liner(line::AbstractString, source::Symbol)
    shift = get_rwo_optical_shift(source)
    line[12+shift] == 'S' && return true
    obscode = view(line, 181+shift:183+shift)
    return obscode == "247" || obscode == "270"
end

function RWOFileReader(text::AbstractString)
    # Parse header
    if Sys.iswindows()
        # Windows uses \r\n as newline instead of \n
        text = replace(text, "\r\n" => '\n')
    end
    @assert contains(text, "END_OF_HEADER\n") "Cannot parse text \
        as it does not contain a header"
    S1 = split(text, "END_OF_HEADER\n")
    header, body = S1
    # Optical header
    if contains(body, NEOCC_OPTICAL_HEADER)
        F1, source = true, :NEOCC
        optical_header = NEOCC_OPTICAL_HEADER
    elseif contains(body, NEODyS2_OPTICAL_HEADER)
        F1, source = true, :NEODyS2
        optical_header = NEODyS2_OPTICAL_HEADER
    else
        F1 = false
    end
    # Radar header
    if contains(body, NEOCC_RADAR_HEADER)
        F2, source = true, :NEOCC
        radar_header = NEOCC_RADAR_HEADER
    elseif contains(body, NEODyS2_RADAR_HEADER)
        F2, source = true, :NEODyS2
        radar_header = NEODyS2_RADAR_HEADER
    else
        F2 = false
    end
    # Parse body
    if F1 && F2
        body = first(split(body, optical_header, keepempty = false))
        body1, body2 = split(body, radar_header)
    elseif F1
        body1 = first(split(body, optical_header, keepempty = false))
        body2 = body[2:1]
    elseif F2
        body1 = body[2:1]
        body2 = first(split(body, radar_header, keepempty = false))
    else
        throw(ArgumentError("Cannot parse text as it has no optical nor radar header"))
    end
    # Parse optical observations
    N = count(==('\n'), body1)
    L = length(body1)
    optical = Vector{SubString{String}}(undef, N)
    a = 1
    for i in eachindex(optical)
        b = findnext('\n', body1, a)
        optical[i] = view(body1, a:b)
        if is_rwo_two_liner(optical[i], source)
            c = findnext('\n', body1, b+1)
            optical[i] = view(body1, a:c)
            a = c + 1
        else
            a = b + 1
        end
        if a > L
            i < N && deleteat!(optical, i+1:N)
            break
        end
    end
    # Parse radar observations
    radar = split(body2, '\n', keepempty = false)

    return RWOFileReader(source, header, optical, radar)
end

"""
    OpticalRWO{T} <: AbstractOpticalAstrometry{T}

An optical astrometric observation in the RWO format.

# Fields

- `design::String`: designation of the object.
- `K::Char`: observation type.
- `T::Char`: observation technology.
- `N::Char`: equivalent to the `note1` field of `OpticalMPC80`.
- `date::DateTime`: date of observation.
- `date_accuracy::T`: time accuracy of `date`.
- `ra::T`: observed right ascension [rad].
- `ra_accuracy::T`: accuracy for `ra`.
- `ra_rms::T`: a-priori formal RMS for `ra` [arcsec].
- `ra_flag::Bool`: manual weight flag for `ra`.
- `ra_bias::T`: bias for `ra` [arcsec].
- `ra_resid::T`: residual for `ra`.
- `dec::T`: observed declination [rad].
- `dec_accuracy::T`: accuracy for `dec`.
- `dec_rms::T`: a-priori formal RMS for `dec` [arcsec].
- `dec_flag::Bool`: manual weight flag for `dec`.
- `dec_bias::T`: bias for `dec` [arcsec].
- `dec_resid::T`: residual for `dec`.
- `mag::T`: apparent magnitude.
- `mag_band::Char`: colour band for `mag`.
- `mag_rms::T`: a-priori formal RMS for `mag`.
- `mag_resid::T`: residual for `mag`.
- `catalogue::CatalogueMPC`: star catalogue used for astrometric reduction.
- `observatory::ObservatoryMPC{T}`: observation station.
- `chi::T`: chi-squared value of the observation residuals.
- `sel_A::Int`: astrometry selection flag.
- `sel_M::Int`: photometry selection flag.
- `source::String`: original record.
- `header::String`: file header.

!!! reference
    The RWO format is described at:
    - https://neo.ssa.esa.int/objects/help#obs
    - https://newton.spacedys.com/neodys/index.php?pc=7.1
"""
@auto_hash_equals fields = (date, ra, dec, observatory) struct OpticalRWO{T} <: AbstractOpticalAstrometry{T}
    design::String
    K::Char
    T::Char
    N::Char
    date::DateTime
    date_accuracy::T
    ra::T
    ra_accuracy::T
    ra_rms::T
    ra_flag::Bool
    ra_bias::T
    ra_resid::T
    dec::T
    dec_accuracy::T
    dec_rms::T
    dec_flag::Bool
    dec_bias::T
    dec_resid::T
    mag::T
    mag_band::Char
    mag_rms::T
    mag_resid::T
    catalogue::CatalogueMPC
    observatory::ObservatoryMPC{T}
    chi::T
    sel_A::Int
    sel_M::Int
    source::String
    header::String
end

# AbstractAstrometryObservation interface
date(x::OpticalRWO) = x.date
band(x::OpticalRWO) = x.mag_band
observatory(x::OpticalRWO) = x.observatory
catalogue(x::OpticalRWO) = x.catalogue
function rms(x::OpticalRWO{T}) where {T <: Real}
    (isnan(x.ra_rms) ? one(T) : x.ra_rms, isnan(x.dec_rms) ? one(T) : x.dec_rms)
end
function debias(x::OpticalRWO{T}) where {T <: Real}
    return (isnan(x.ra_bias) ? zero(T) : x.ra_bias, isnan(x.dec_bias) ? zero(T) : x.dec_bias)
end
corr(::OpticalRWO{T}) where {T <: Real} = zero(T)

trackletid(x::OpticalRWO) = ""

isdeprecated(x::OpticalRWO) = x.T == 'X'

# Print method for OpticalRWO
show(io::IO, o::OpticalRWO) = print(io, o.design, " α: ", @sprintf("%.5f",
    rad2deg(o.ra)), "° δ: ", @sprintf("%.5f", rad2deg(o.dec)), "° t: ", o.date,
    " obs: ", o.observatory.name)

# Parsing

rwoparse(_, ::Type{String}, x) = String(x)
rwoparse(_, ::Type{Char}, x) = isempty(x) ? ' ' : first(x)
rwoparse(_, ::Type{CatalogueMPC}, x) = isempty(x) ? unknowncat() : search_catalogue_code(first(x))

function rwoparse(_, ::Type{ObservatoryMPC{T}}, x) where {T <: Real}
    if contains(x, '-')
        return search_observatory_code(JPL_TO_MPC_OBSCODES[x])
    else
        return search_observatory_code(x)
    end
end

function rwoparse(_, ::Type{DateTime}, x)
    if contains(x, ':')
        date = DateTime(x, "yyyy mm dd HH:MM:SS")
    else
        date = DateTime(view(x, 1:10), "yyyy mm dd")
        fraction = parse(Float64, view(x, 11:length(x)))
        return date + Microsecond(round(Int, fraction * 8.64e10))
    end
end

function rwoparse(name, ::Type{T}, x) where {T <: Real}
    if name == :ra
        return hms2rad(x)
    elseif name == :dec
        return dms2rad(x)
    else
        return isempty(x) ? T(NaN) : parse(T, x)
    end
end

function rwoparse(_, ::Type{Bool}, x)
    if x == "F"
        return false
    elseif x == "T"
        return true
    else
        return parse(Bool, x)
    end
end

OpticalRWO(r::DataFrameRow) = OpticalRWO{Float64}(
    r.design, r.K, r.T, r.N, r.date, r.date_accuracy, r.ra, r.ra_accuracy,
    r.ra_rms, r.ra_flag, r.ra_bias, r.ra_resid, r.dec, r.dec_accuracy, r.dec_rms,
    r.dec_flag, r.dec_bias, r.dec_resid, r.mag, r.mag_band, r.mag_rms, r.mag_resid,
    r.catalogue, r.observatory, r.chi, r.sel_A, r.sel_M, r.source, r.header
)

function parse_optical_rwo(text::AbstractString)
    # File reader
    reader = RWOFileReader(text)
    L = length(reader.optical)
    shift = get_rwo_optical_shift(reader.source)
    columns = get_rwo_optical_columns(reader.source)
    # Construct DataFrame
    R = OpticalRWO{Float64}
    names, types = fieldnames(R), fieldtypes(R)
    df = DataFrame([fill(astrometrydefault(fieldtype(R, name)), L) for name in names],
        collect(names))
    for (i, line) in enumerate(reader.optical)
        for (name, type, idxs) in zip(names, types, columns)
            x = strip(view(line, idxs))
            df[i, name] = rwoparse(name, type, x)
        end
    end
    # Two-line observations
    mask = findall(istwoliner, df.observatory)
    for i in mask
        # Skip observations with a 'x' note (e.g. 2010 BO127)
        df.T[i] == 'x' && continue
        obs = df.observatory[i]
        line = view(reader.optical[i], 199+shift:length(reader.optical[i]))
        if isoccultation(obs) || issatellite(obs)
            frame = line[35+shift] == '1' ? "ICRF_KM" : "ICRF_AU"
        else
            frame = obs.frame
        end
        vals = split(line)
        coords = SVector{3, Float64}(
            parse(Float64, vals[end-3]),
            parse(Float64, vals[end-2]),
            parse(Float64, vals[end-1])
        )
        df.observatory[i] = ObservatoryMPC(obs; frame, coords)
    end
    # Source string
    df.source = reader.optical
    # File header
    df.header .= reader.header
    # Parse observations
    optical = OpticalRWO.(eachrow(df))
    # Eliminate repeated entries
    unique!(optical)
    # Sort by date
    sort!(optical)

    return optical
end

"""
    fetch_optical_rwo(id, source)

Return the optical astrometry of minor body `id` in the RWO format.
The `source` of the observations can be either `NEOCC` or `NEODyS2`.

!!! reference
    The NEOCC observations API is described at:
    - https://neo.ssa.esa.int/computer-access
"""
function fetch_optical_rwo(id::AbstractString, ::Type{NEOCC})
    # Get and parse HTTP response
    text = fetch_http_text(NEOCC; id = replace(id, " " => ""))
    # Parse observations
    optical = parse_optical_rwo(text)

    return optical
end

function fetch_optical_rwo(id::AbstractString, ::Type{NEODyS2})
    # Get and parse HTTP response
    text = fetch_http_text(NEODyS2; id = replace(id, " " => ""))
    # Parse observations
    optical = parse_optical_rwo(text)

    return optical
end

# Read / write

"""
    read_optical_rwo(filename)

Read from `filename` a vector of optical astrometry in the RWO format.
"""
function read_optical_rwo(filename::AbstractString)
    # Read file
    text = read(filename, String)
    # Parse observations
    optical = parse_optical_rwo(text)

    return optical
end

function get_rwo_optical_header(x::OpticalRWO)
    @unpack source, K, T, N = x
    if source[12] == K && source[14] == T && source[16] == N
        return NEODyS2_OPTICAL_HEADER
    elseif source[15] == K && source[17] == T && source[19] == N
        return NEOCC_OPTICAL_HEADER
    else
        throw(ArgumentError("Cannot match observation to an optical header"))
    end
end

"""
    write_optical_rwo(obs, filename)

Write `obs` to `filename` in the RWO format.
"""
function write_optical_rwo(obs::AbstractVector{OpticalRWO{T}},
                           filename::AbstractString) where {T <: Real}
    @assert !isempty(obs) "The observation vector is empty"
    header = obs[1].header
    optical_header = get_rwo_optical_header(obs[1])
    open(filename, "w") do file
        # File Header
        write(file, header, "END_OF_HEADER\n")
        # Optical observations header
        write(file, optical_header)
        # Optical observations
        for i in eachindex(obs)
            write(file, obs[i].source)
        end
    end
end