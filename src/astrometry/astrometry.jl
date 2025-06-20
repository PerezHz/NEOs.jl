# Common interface
include("abstractastrometry.jl")
include("utils.jl")
include("cataloguempc.jl")
include("observatorympc.jl")
include("topocentric.jl")
include("sources.jl")
# Optical astrometry
include("opticalmpc80.jl")
include("neocpobject.jl")
include("opticalrwo.jl")
include("opticalades.jl")
include("opticaltracklet.jl")
include("opticalresidual.jl")
include("weightingscheme.jl")
include("debiasingscheme.jl")
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