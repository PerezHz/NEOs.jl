include("catalogue_mpc.jl")
include("observatory_mpc.jl")
include("radec_mpc.jl")
include("radar_jpl.jl")
include("neocp.jl")
include("process_radar.jl")
include("units.jl")
include("jpl_eph.jl")
include("topocentric.jl")
include("tracklet.jl")
include("errormodel/debiasingscheme.jl")
include("process_radec.jl")

@doc raw"""
    numberofdays(::AbstractVector)

Return the timespan of a vector of dates in days.
"""
function numberofdays(dates::AbstractVector{DateTime})
    t0, tf = extrema(dates)
    return (tf - t0).value / 86_400_000
end

function numberofdays(radec::AbstractVector{RadecMPC{T}}) where {T <: Real}
    t0, tf = extrema(date, radec)
    return (tf - t0).value / 86_400_000
end

function numberofdays(tracklets::AbstractVector{Tracklet{T}}) where {T <: Real}
    dates = map(t -> extrema(date, t.radec), tracklets)
    t0, tf = minimum(first, dates), maximum(last, dates)
    return (tf - t0).value / 86_400_000
end