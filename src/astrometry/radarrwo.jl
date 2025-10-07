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
- `rms::T`: a-priori RMS [same units as `measure`].
- `flag::Bool`: manual weight flag.
- `bias::T`: bias of `measure` [same units as `measure`].
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

# AbstractAstrometryObservation interface
date(x::RadarRWO) = x.date

measure(x::RadarRWO) = isdelay(x) ? range2delay(x.measure) :
    rangerate2doppler(x.measure, frequency(x))

# The RWO format contains no information about the transmitter frequency.
# The function below is based on the 4,591 observations return by the JPL
# small body radar astrometry API by 19 June 2025
# This approach was inspired by:
# # https://github.com/Bill-Gray/find_orb/blob/3fc52c93b3263fdf8cda7a7bb6ee5f4df532f2cd/mpc_obs.cpp#L2316
frequency(x::RadarRWO) = frequency(x.trx, date(x))

function frequency(observatory::ObservatoryMPC{T}, date::DateTime) where {T <: Real}
    code = observatory.code
    year = Dates.year(date)
    # Arecibo uses 2380 MHz almost all the time
    if code == "251"
        # except its first observation of Eros (1975-01-23)
        if year == 1975
            freq = 430.0
        # and its time-delay observations of 2018 EB (2018-10-05 to 2018-10-07)
        elseif Date(2018, 10, 5) < date < Date(2018, 10, 8) && Time(date) > Time(8, 40)
            freq = 8560.0
        else
            freq = 2380.0
        end
    # Haystack, Westford only has two Icarus observations in 1968 both in 7840.0 MHz
    elseif code == "254"
        freq = 7840.0
    # Green Bank has no observations as the transmitter antenna
    # elseif code == "256"
    # Goldstone DSS 13, Fort Irwin uses 7190 MHz
    elseif code == "252"
        # except its Icarus observations in 1968
        if year < 1969
            freq = 2388.0
        else
            freq = 7190.0
        end
    # Goldstone DSS 14, Fort Irwin
    elseif code == "253"
        # observations before 1991
        if year < 1991
            freq = 8495.0
        # observations between 1991 and 1999
        elseif year < 1999
            freq = 8510.0
        # observations after 1999 use 8560 MHz
        else
            # except two 2008 QS11 observations (2008-10-01) and
            # One 523664 observation (2024-07-26)
            if (DateTime(2008, 10, 1, 9) < date < DateTime(2008, 10, 1, 12)) ||
                date == DateTime(2024, 7, 26, 15, 10)
                freq = 8650.0
            # except one 1685 observation (2016-01-19)
            elseif date == DateTime(2016, 1, 19, 4)
                freq = 8660.0
            else
                freq = 8560.0
            end
        end
    # Canberra DSS 34, DSS 35, DSS 36 and DSS 43 accumulate 38 observations
    # between 2019 and 2025 all in 7159.45 MHz
    elseif code in ("272", "263", "264", "265")
        freq = 7159.45
    # Yevpatoriya has no observations as the transmitter antenna
    # elseif code == "255"
    # ATCA DSS 47 has no observations as the transmitter antenna
    # EISCAT Tromso UHF only has one 367943 observation in 929.6 MHz (2013-02-15)
    elseif code == "259"
        freq = 929.6
    # Ceduna 30, UTAS and Usuda DSS 64 have no observations as the transmitter antenna
    # elseif code in ("287", "308")
    # Unknown observatory
    else
        freq = T(NaN)
    end

    return freq
end

observatory(x::RadarRWO) = x.rcx
rms(x::RadarRWO) = isdelay(x) ? range2delay(x.rms) :
    abs(rangerate2doppler(x.rms, frequency(x)))
debias(x::RadarRWO) = isdelay(x) ? range2delay(x.bias) :
    rangerate2doppler(x.bias, frequency(x))

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