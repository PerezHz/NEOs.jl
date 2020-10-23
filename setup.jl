@show Threads.nthreads()
import Pkg
Pkg.instantiate()
using Apophis

using EarthOrientation, SatelliteToolbox
#EarthOrientation.update()
EarthOrientation.update(force = true)
get_iers_eop(force_download = true)
