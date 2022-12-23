# This file is part of the NEOs.jl package; MIT licensed

testfiles = (
    "catalogue_mpc.jl",
    "observatory_mpc.jl",
    "radec_mpc.jl",
    )

for file in testfiles
    include(file)
end