# Common interface
include("abstractastrometry.jl")
include("designations.jl")
include("mpec.jl")
include("cataloguempc.jl")
include("observatorympc.jl")
include("magnitudeband.jl")
include("topocentric.jl")
include("sources.jl")
# Optical astrometry
include("opticalmpc80.jl")
include("neocpobject.jl")
include("opticalrwo.jl")
include("opticalades.jl")
include("opticaltracklet.jl")
include("computeradec.jl")
include("opticalresidual.jl")
include("weightingscheme.jl")
include("debiasingscheme.jl")

"""
    read_optical_astrometry(filename; format = :auto)

Read optical astrometry from `filename`, optionally detecting the file format.
Supported formats are `:auto`, `:ades`, `:mpc80`, and `:obs80`.
"""
function read_optical_astrometry(filename::AbstractString; format = :auto)
    fmt = _normalize_optical_astrometry_format(format)
    fmt = fmt === :auto ? _detect_optical_astrometry_format(filename) : fmt
    return _read_optical_astrometry(Val(fmt), filename)
end

function _normalize_optical_astrometry_format(format)
    fmt = lowercase(strip(String(format)))
    if fmt in ("auto", "ades", "xml")
        return fmt == "xml" ? :ades : Symbol(fmt)
    elseif fmt in ("mpc80", "obs80")
        return :mpc80
    else
        throw(ArgumentError(
            "Unknown input format: $format. Use auto, ades, mpc80, or obs80."
        ))
    end
end

function _detect_optical_astrometry_format(filename::AbstractString)
    for line in eachline(filename)
        stripped = strip(line)
        isempty(stripped) && continue
        return startswith(stripped, '<') ? :ades : :mpc80
    end
    throw(ArgumentError("Cannot detect astrometry format from empty file: $filename"))
end

_read_optical_astrometry(::Val{:ades}, filename::AbstractString) =
    read_optical_ades(filename)

_read_optical_astrometry(::Val{:mpc80}, filename::AbstractString) =
    read_optical_mpc80(filename)
# Radar astrometry
include("radarjpl.jl")
include("radarrwo.jl")
include("radarresidual.jl")

"""
    numberofdays(::AbstractVector)

Return the timespan of a vector of dates in days.
"""
function numberofdays(dates::AbstractVector{DateTime})
    t0, tf = extrema(dates)
    return (tf - t0).value / 86_400_000
end
