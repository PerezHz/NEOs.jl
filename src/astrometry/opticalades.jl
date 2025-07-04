"""
    OpticalADES{T} <: AbstractOpticalAstrometry{T}

An optical astrometric observation in the Minor Planet Center ADES format.

# Fields

- `permid::String`: IAU permanent designation.
- `provid::String`: MPC unpacked provisional designation.
- `obsid::String`: globally unique observation identifier.
- `trkid::String`: globally unique tracklet identifier.
- `mode::String`: mode of instrumentation.
- `stn::ObservatoryMPC{T}`: observing station.
- `sys::String`: coordinate frame for the observatory position.
- `ctr::Int`: origin of the reference frame given by `sys`.
- `pos1/pos2/pos3::T`: position of the observer. Interpretation depends on the
    value of `sys` (see the Extended help of `ObservatoryMPC`).
- `vel1/vel2/vel3::T`: ICRF velocity of space-based observer.
- `poscov11/poscov12/poscov13/poscov22/poscov23/poscov33::T`: Upper triangular
    part of the `(pos1, pos2, po3)` covariance matrix [same units of position
    coordinates].
- `prog::String`: program code assigned by the MPC.
- `obstime::DateTime`: UTC date and time of the observation.
- `rmstime::T`: random uncertainty in `obstime` [sec].
- `ra::T`: observed right ascension [rad].
- `dec::T`: observed declination [rad].
- `rastar/decstar::T`: for occultation observations, the right ascension and
    declination of the occulted star [rad].
- `deltara/deltadec::T`: measured shift in right ascension and declination for
    occultation observations with respect to the star specified by `raStar` and
    `decStar` [rad].
- `rmsra/rmsdec::T`: for `ra`/`dec` and `deltara/deltadec` observations, the random
    component of the uncertainty [arcsec].
- `rmscorr::T`: correlation between `ra` and `dec`.
- `astcat::CatalogueMPC`: star catalog used for the astrometric reduction.
- `mag::T`: apparent magnitude.
- `rmsmag::T`: 1-sigma uncertainty of `mag`.
- `band::String`: passband designation for photometry.
- `photcat::CatalogueMPC`: star catalog used for the photometric reduction.
- `ref::String`: standard reference field used for citations.
- `disc::String`: discovery flag.
- `subfmt::String`: submission format.
- `prectime::Int`: precision of the reported observation time [millionths of a day].
- `precra/precdec::T`: reported precision for archival astrometry [rad].
- `unctime::T`: estimated systematic time error [sec].
- `notes::String`: a set of one-character note flags to communicate observing
    circumstances.
- `remarks::String`: a comment provided by the observer.
- `deprecated::String`: deprecation flag.
- `source::String`: original record.

!!! reference
    The Minor Planet Center ADES format is described at:
    - https://github.com/IAU-ADES/ADES-Master/blob/master/ADES_Description.pdf
"""
@auto_hash_equals fields = (obstime, ra, dec, stn) struct OpticalADES{T} <: AbstractOpticalAstrometry{T}
    # Identification Group Elements
    permid::String
    provid::String
    obsid::String
    trkid::String
    mode::String
    stn::ObservatoryMPC{T}
    # Location Group Elements
    sys::String
    ctr::Int
    pos1::T
    pos2::T
    pos3::T
    vel1::T
    vel2::T
    vel3::T
    poscov11::T
    poscov12::T
    poscov13::T
    poscov22::T
    poscov23::T
    poscov33::T
    prog::String
    obstime::DateTime
    rmstime::T
    # Observation Group Elements
    ra::T
    dec::T
    rastar::T
    decstar::T
    deltara::T
    deltadec::T
    rmsra::T
    rmsdec::T
    rmscorr::T
    astcat::CatalogueMPC
    # Photometry Group Elements
    mag::T
    rmsmag::T
    band::String
    photcat::CatalogueMPC
    ref::String
    disc::String
    subfmt::String
    # Precision Group Elements
    prectime::Int
    precra::T
    precdec::T
    unctime::T
    notes::String
    remarks::String
    deprecated::String
    source::String
end

# AbstractAstrometryObservation interface
date(x::OpticalADES) = x.obstime
band(x::OpticalADES) = first(x.band)
observatory(x::OpticalADES) = x.stn
catalogue(x::OpticalADES) = x.astcat
function rms(x::OpticalADES{T}) where {T <: Real}
    return (isnan(x.rmsra) ? one(T) : x.rmsra, isnan(x.rmsdec) ? one(T) : x.rmsdec)
end
debias(x::OpticalADES{T}) where {T <: Real} = (zero(T), zero(T))

# Print method for OpticalADES
function show(io::IO, o::OpticalADES)
    # If there is no number, use temporary designation
    id_str = isempty(o.permid) ? o.provid : o.permid

    print(io, id_str, " α: ", @sprintf("%.5f", rad2deg(o.ra)), "° δ: ", @sprintf("%.5f",
        rad2deg(o.dec)), "° t: ", o.obstime, " obs: ", o.stn.name)
end

# Parsing

adesparse(_, ::Type{String}, x) = String(x)
function adesparse(_, ::Type{CatalogueMPC}, x)
    if isempty(x) || !(x in keys(CATALOGUE_MPC_NAMES_TO_CODES))
        return unknowncat()
    else
        return search_catalogue_code(CATALOGUE_MPC_NAMES_TO_CODES[x])
    end
end
adesparse(_, ::Type{ObservatoryMPC{T}}, x) where {T <: Real} = search_observatory_code(x)
adesparse(_, ::Type{Int}, x) = isa(x, Number) ? round(Int, x) : parse(Int, x)

function adesparse(_, ::Type{DateTime}, x)
    date = DateTime(view(x, 1:19))
    fraction = length(x) > 20 ? parse(Float64, view(x, 20:length(x)-1)) : zero(Float64)
    return date + Microsecond( round( Int, fraction * 1e6 ) )
end

function adesparse(name, ::Type{T}, x) where {T <: Real}
    isempty(x) && return T(NaN)
    y = parse(T, x)
    if name in (:ra, :dec, :rastar, :decstar)
        return deg2rad(y)
    elseif name in (:deltara, :deltadec, :precra, :precdec)
        return arcsec2rad(y)
    else # name in (:rmsra, :rmsdec, etc.)
        return y
    end
end

OpticalADES(r::DataFrameRow) = OpticalADES{Float64}(
    r.permid, r.provid, r.obsid, r.trkid, r.mode, r.stn, r.sys, r.ctr,
    r.pos1, r.pos2, r.pos3, r.vel1, r.vel2, r.vel3, r.poscov11, r.poscov12,
    r.poscov13, r.poscov22, r.poscov23, r.poscov33, r.prog, r.obstime,
    r.rmstime, r.ra, r.dec, r.rastar, r.decstar, r.deltara, r.deltadec,
    r.rmsra, r.rmsdec, r.rmscorr, r.astcat, r.mag, r.rmsmag, r.band,
    r.photcat, r.ref, r.disc, r.subfmt, r.prectime, r.precra, r.precdec,
    r.unctime, r.notes, r.remarks, r.deprecated, r.source
)

function parse_optical_ades(text::String)
    # Parse XML
    node = parse(LazyNode, text)
    L = length(children(node[2]))
    # Construct DataFrame
    R = OpticalADES{Float64}
    names = fieldnames(R)
    df = DataFrame([fill(astrometrydefault(fieldtype(R, name)), L) for name in names],
        collect(names))
    for (i, line) in enumerate(children(node[2]))
        for child in children(line)
            name = Symbol(lowercase(tag(child)))
            if name in names
                type = fieldtype(R, name)
                x = strip(XML.value(first(child)))
                df[i, name] = adesparse(name, type, x)
            end
        end
    end
    # Satellite observatories
    mask = findall(istwoliner, df.stn)
    for i in mask
        frame = df.sys[i]
        coords = SVector{3, Float64}(df.pos1[i], df.pos2[i], df.pos3[i])
        df.stn[i] = ObservatoryMPC(df.stn[i]; frame, coords)
    end
    # Occultation observations
    mask = findall(==("OCC"), df.mode)
    for i in mask
        df.dec[i] = df.decstar[i] + df.deltadec[i]
        df.ra[i] = df.rastar[i] + df.deltara[i] / cos(df.dec[i])
    end
    # Source string
    df.source = XML.write.(children(node[2]))
    # Parse observations
    optical = OpticalADES.(eachrow(df))

    return optical
end

"""
    fetch_optical_ades(id, source)

Return the optical astrometry of minor body `id` in the Minor Planet Center
ADES format. The `source` of the observations can be either `MPC` or `NEOCP`.

!!! reference
    The Minor Planet Center observations APIs are described at:
    - https://minorplanetcenter.net/mpcops/documentation/observations-api/
    - https://minorplanetcenter.net/mpcops/documentation/neocp-observations-api/
"""
function fetch_optical_ades(id::AbstractString, ::Type{MPC})
    # Get and parse HTTP response
    text = fetch_http_text(MPC; mode = 2, id = id, format = "XML")
    # Parse JSON
    dict = JSON.parse(text)
    # Parse observations
    optical = parse_optical_ades(dict[1]["XML"])

    return optical
end

function fetch_optical_ades(id::AbstractString, ::Type{NEOCP})
    # Get and parse HTTP response
    text = fetch_http_text(NEOCP; mode = 1, id = id, format = "XML")
    # Parse JSON
    dict = JSON.parse(text)
    # Parse observations
    optical = parse_optical_ades(dict[1]["XML"])

    return optical
end

"""
    read_optical_ades(filename)

Read from `filename` a vector of optical astrometry in the Minor Planet
Center ADES format.
"""
function read_optical_ades(filename::AbstractString)
    # Read file
    text = read(filename, String)
    # Parse observations
    optical = parse_optical_ades(text)

    return optical
end

"""
    write_optical_ades(obs, filename)

Write `obs` to `filename` in the Minor Planet Center ADES format.
"""
function write_optical_ades(obs::AbstractVector{OpticalADES{T}},
                           filename::AbstractString) where {T <: Real}
    open(filename, "w") do file
        write(file, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
        write(file, "<ades version=\"2017\">\n")
        # Optical observations
        for i in eachindex(obs)
            for s in eachsplit(obs[i].source, '\n')
                write(file, "  ", s, "\n")
            end
        end
        write(file, "</ades>")
    end
end