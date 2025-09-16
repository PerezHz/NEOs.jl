# This file is part of the NEOs.jl package; MIT licensed

testfiles = (
    "common.jl",
    "astrometry.jl",
    "optical.jl",
    "radar.jl",
    "propagation.jl",
    "orbitdetermination.jl",
    "impactmonitoring.jl",
    "aqua.jl",
)

for file in testfiles
    include(file)
end