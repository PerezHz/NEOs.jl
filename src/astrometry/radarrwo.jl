"""
    RadarRWO{T} <: AbstractRadarAstrometry{T}

A radar astrometric observation in the RWO format.

# Fields

- `design::String`: designation of the object.
- `K::Char`: observation type.
- `T::Char`: observation technology.
- `N::Char`: note on observational circumstances.
- `date::DateTime`: time of observation in calendar format.
- `measure::T`: radar range [km] or range rate [km/d].
- `accuracy::T`: measurement accuracy.
- `rms::T`: a-priori RMS.
- `flag::Bool`: manual weight flag.
- `bias::T`: bias of `measure`.
- `resid::T`: residual of `measure`.
- `trx::ObservatoryMPC{T}`: transmitting station.
- `rcx::ObservatoryMPC{T}`: receiving station.
- `chi::T`: chi-squared value.
- `S::Bool`: selection flag.
- `source::String`: original record.
- `header::String`: file header.

!!! reference
    The RWO format is described at:
    - https://neo.ssa.esa.int/objects/help#obs
    - https://newton.spacedys.com/neodys/index.php?pc=7.1
"""
@auto_hash_equals fields = (date, measure, K, rcx, trx) struct RadarRWO{T} <: AbstractRadarAstrometry{T}
    design::String
    K::Char
    T::Char
    N::Char
    date::DateTime
    measure::T
    accuracy::T
    rms::T
    flag::Bool
    bias::T
    resid::T
    trx::ObservatoryMPC{T}
    rcx::ObservatoryMPC{T}
    chi::T
    S::Bool
    source::String
    header::String
end

date(x::RadarRWO) = x.date

isdelay(r::RadarRWO) = r.K == 'R'
isdoppler(r::RadarRWO) = r.K == 'V'
ismonostatic(r::RadarRWO) = r.trx == r.rcx

# Print method for RadarRWO
function show(io::IO, r::RadarRWO)
    units = r.K == 'R' ? "km  " : "km/d"
    print(io, r.design, " value: ", @sprintf("%+e", r.measure), " ", units, " sigma: ",
        @sprintf("%.3e", r.rms), " ", units, " t: ", r.date, " rcvr: ", r.rcx.name)
end

# Parsing

RadarRWO(r::DataFrameRow) = RadarRWO{Float64}(
    r.design, r.K, r.T, r.N, r.date, r.measure, r.accuracy, r.rms, r.flag,
    r.bias, r.resid, r.trx, r.rcx, r.chi, r.S, r.source, r.header
)

function parse_radar_rwo(text::String)
    # File reader
    reader = RWOFileReader(text)
    L = length(reader.radar)
    # Construct DataFrame
    R = RadarRWO{Float64}
    names, types = fieldnames(R), fieldtypes(R)
    df = DataFrame([fill(astrometrydefault(fieldtype(R, name)), L) for name in names],
        collect(names))
    for (i, line) in enumerate(reader.radar)
        for (name, type, idxs) in zip(names, types, RWO_RADAR_COLUMNS)
            x = strip(view(line, idxs))
            df[i, name] = rwoparse(name, type, x)
        end
    end
    # Source string
    df.source = reader.radar
    # File header
    df.header .= reader.header
    # Parse observations
    radar = RadarRWO.(eachrow(df))
    # Eliminate repeated entries
    unique!(radar)
    # Sort by date
    sort!(radar)

    return radar
end

"""
    fetch_radar_rwo(id, source)

Return the radar astrometry of minor body `id` in the RWO format.
The `source` of the observations can be either `NEOCC` or `NEODyS2`.

!!! reference
    The NEOCC observations API is described at:
    - https://neo.ssa.esa.int/computer-access
"""
function fetch_radar_rwo(id::AbstractString, ::Type{NEOCC})
    # Get and parse HTTP response
    text = fetch_http_text(NEOCC; id = replace(id, " " => ""))
    # Parse observations
    radar = parse_radar_rwo(text)

    return radar
end

function fetch_radar_rwo(id::AbstractString, ::Type{NEODyS2})
    # Get and parse HTTP response
    text = fetch_http_text(NEODyS2; id = replace(id, " " => ""))
    # Parse observations
    radar = parse_radar_rwo(text)

    return radar
end

# Read / write

"""
    read_radar_rwo(filename)

Read from `filename` a vector of radar astrometry in the RWO format.
"""
function read_radar_rwo(filename::AbstractString)
    # Read file
    text = read(filename, String)
    # Parse observations
    radar = parse_radar_rwo(text)

    return radar
end

"""
    write_radar_rwo(obs, filename)

Write `obs` to `filename` in the RWO format.
"""
function write_radar_rwo(obs::AbstractVector{RadarRWO{T}},
                         filename::AbstractString) where {T <: Real}
    open(filename, "w") do file
        # File Header
        write(file, obs[1].header, "END_OF_HEADER\n")
        # Radar observations header
        write(file, RWO_RADAR_HEADER)
        # Radar observations
        for i in eachindex(obs)
            write(file, obs[i].source, "\n")
        end
    end
end