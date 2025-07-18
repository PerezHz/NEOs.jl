# This file is part of the NEOs.jl package; MIT licensed

testfiles = (
    "common.jl",
    "astrometry/common.jl",
    "astrometry/optical.jl",
    "astrometry/radar.jl",
    "propagation.jl",
    "orbitdetermination/orbitdetermination.jl",
    "bplane.jl",
    "aqua.jl",
)

for file in testfiles
    include(file)
end