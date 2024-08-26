using ArgParse, Downloads, Random, HORIZONS, NEOs
using NEOs: numberofdays

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
            help = "Number of NEOs"
            arg_type = Int
            default = 100
        "--output", "-o"
            help = "Output file"
            arg_type = String
            default = "names.txt"
    end

    s.epilog = """
        examples:\n
        \n
        julia --project names.jl -N 250 -o names.txt\n
        \n
    """

    return parse_args(s)
end

function main()

    # Parse arguments from commandline
    parsed_args = parse_commandline()
    # Number of NEOs
    N::Int = parsed_args["N"]
    # Output file
    outfile::String = parsed_args["output"]
    # Print header
    println("NEOs selector for orbitdetermination.jl")
    println("• Number of NEOs: ", N)
    println("• Output file: ", outfile)

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
    # Shuffle the provisional designations list
    shuffle!(provdesig)

    # Assemble a list with N NEOs
    names = Vector{String}(undef, N)

    for i in eachindex(names)
        while !isempty(provdesig)
            neo = String(pop!(provdesig))
            radec = fetch_radec_mpc("designation" => neo)
            jplorbit = sbdb("des" => neo)["orbit"]
            if (numberofdays(radec) <= 15.0) &&  (jplorbit["n_obs_used"] == length(radec)) &&
                isnothing(jplorbit["n_del_obs_used"]) && isnothing(jplorbit["n_dop_obs_used"])
                names[i] = neo
                break
            end
        end
    end
    sort!(names)

    # Save names list
    open(outfile, "w") do file
        for i in eachindex(names)
            write(file, names[i], i == N ? "" : "\n")
        end
    end

    println("• Saved ", N, " names to ", outfile)

    nothing
end

main()