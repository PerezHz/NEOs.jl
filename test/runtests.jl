# This file is part of the NEOs.jl package; MIT licensed

testfiles = (
    "osculating.jl",
    "observations.jl",
    "propagation.jl",
    "orbit_determination.jl",
    "dataframes.jl",
    )

for file in testfiles
    include(file)
end