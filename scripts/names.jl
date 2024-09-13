using ArgParse, Downloads, Random, NEOs
using NEOs: RadecMPC, numberofdays, issatellite, reduce_tracklets

# Potentially Hazardous Asteroids MPC list
const PHAS_URL = "https://cgi.minorplanetcenter.net/cgi-bin/textversion.cgi?f=lists/PHAs.html"
# Amors MPC list
const AMORS_URL = "https://cgi.minorplanetcenter.net/cgi-bin/textversion.cgi?f=lists/Amors.html"
# Apollos MPC list
const APOLLOS_URL = "https://cgi.minorplanetcenter.net/cgi-bin/textversion.cgi?f=lists/Apollos.html"
# Atens MPC list
const ATENS_URL = "https://cgi.minorplanetcenter.net/cgi-bin/textversion.cgi?f=lists/Atens.html"

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "names.jl"
    # Desciption (for help screen)
    s.description = "Select a list of NEOs for orbitdetermination.jl"

    @add_arg_table! s begin
        "--N", "-N"
            help = "number of NEOs"
            arg_type = Int
            default = 0
        "--output", "-o"
            help = "output file"
            arg_type = String
            default = "names.txt"
    end

    s.epilog = """
        note: if izero(N), then the script will save all compatible NEOs\n
        \n
        examples:\n
        # Save 1000 NEOs to names.txt\n
        julia --project names.jl -N 1000 -o names.txt\n
    """

    return parse_args(s)
end

function isodcompatible(radec::Vector{RadecMPC{T}}) where {T <: Real}
    # Eliminate observations before oficial discovery
    firstobs = findfirst(r -> !isempty(r.discovery), radec)
    isnothing(firstobs) && return false
    radec = radec[firstobs:end]
    # Filter out incompatible observatories
    filter!(radec) do r
        hascoord(r.observatory) && !issatellite(r.observatory)
    end
    length(radec) < 3 && return false
    # There is at least one set of 3 tracklets with a < 15 days timespan
    tracklets = reduce_tracklets(radec)
    for i in 1:length(tracklets)-2
        numberofdays(tracklets[i:i+2]) > 15.0 && continue
        tracklets = tracklets[i:i+2]
        radec = reduce(vcat, getfield.(tracklets, :radec))
        sort!(radec)
        break
    end
    return numberofdays(radec) <= 15.0
end

function main()

    # Parse arguments from commandline
    parsed_args = parse_commandline()
    # Number of NEOs
    N::Int = parsed_args["N"]
    # Output file
    output::String = parsed_args["output"]
    # Print header
    println("NEOs selector for orbitdetermination.jl")
    println("• Number of NEOs: ", N)
    println("• Output file: ", output)

    # Parse PHAs, Atens, Apollos and Amors
    provdesig = Vector{SubString{String}}(undef, 0)

    for url in [PHAS_URL, ATENS_URL, APOLLOS_URL, AMORS_URL]
        # Download MPC file
        Downloads.download(url, "tmp.txt")
        # Parse lines
        lines = readlines("tmp.txt")[3:end]
        filter!(!isempty, lines)
        filter!(l -> isempty(strip(l[1:27])), lines)
        # Parse provisional designations
        provdesig = vcat(provdesig, map(l -> strip(l[28:40]), lines))
        # Delete MPC file
        rm("tmp.txt")
    end
    # Delete repeated NEOs
    unique!(provdesig)
    # We can only process asteroids discovered between 2000 and 2024
    filter!(desig -> "2000" <= desig[1:4] <= "2024", provdesig)

    if 0 < N < length(provdesig)
        # Shuffle the provisional designations list
        shuffle!(provdesig)
    else
        # Consider all compatible NEOs
        N = length(provdesig)
    end

    # Assemble a list with N NEOs
    names = Vector{String}(undef, N)

    i = 1
    while i <= N && !isempty(provdesig)
        neo = String(popfirst!(provdesig))
        radec = fetch_radec_mpc("designation" => neo)
        if isodcompatible(radec)
            names[i] = neo
            i += 1
        end
    end
    # Eliminate #undef elements
    names = names[1:i-1]
    # Sort provisional designations
    sort!(names)

    # Save names list
    open(output, "w") do file
        for i in eachindex(names)
            write(file, names[i], i == N ? "" : "\n")
        end
    end

    println("• Saved ", N, " names to ", output)

    nothing
end

main()