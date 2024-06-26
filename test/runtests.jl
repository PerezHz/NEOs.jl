# This file is part of the NEOs.jl package; MIT licensed

testfiles = (
    "osculating.jl",
    "observations.jl",
    "propagation.jl",
    "orbitdetermination.jl",
    "bplane.jl",
    "dataframes.jl",
    "aqua.jl",
    )

for file in testfiles
    include(file)
end