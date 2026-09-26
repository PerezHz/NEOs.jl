using ArgParse, HTTP, JSON, Dates, Downloads, ZipArchives

const MAQI_BASE_URL = "https://maqi.astro.umd.edu"
const MAQI_API_BASE_URL = "https://maqi-api.astro.umd.edu"

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "maqi.jl"
    # Desciption (for help screen)
    s.description = "Post a query to MAQI to fetch the result"

    s.epilog = """
        Example:\n
        \n
        julia --project maqi.jl -e **** -p ****\n
        \n
    """

    @add_arg_table! s begin
        "--email", "-e"
            help = "MAQI account email"
            arg_type = String
        "--password", "-p"
            help = "MAQI account password"
            arg_type = String
        "--query", "-q"
            help = "input query"
            arg_type = String
        "--output", "-o"
            help = "output file"
            arg_type = String
        "--poll_interval"
            help = "time between polls in seconds"
            arg_type = Int
            default = 10
        "--max_retries"
            help = "maximum number of polling attempts"
            arg_type = Int
            default = 10
    end

    return parse_args(s)
end

# Authentication
function maqi_login(email::AbstractString, password::AbstractString)
    login_body = JSON.json(Dict("email" => email, "password" => password))
    login_headers = [
        "Content-Type" => "application/json",
        "Connection" => "close"
    ]
    login_resp = HTTP.post(
        "$MAQI_BASE_URL/api/auth/login",
        login_headers,
        login_body;
        retry_non_idempotent = true,
        retry = true
    )
    login_data = JSON.parse(String(login_resp.body))
    login_date = DateTime(login_resp.headers["Date"][1:end-4], RFC1123Format)
    token, expiration_seconds = login_data["idToken"], login_data["expiresIn"]
    expiration_date = login_date + Second(expiration_seconds)
    return token, expiration_date
end

function async_maqi_query!(io::IOBuffer, token::AbstractString, query::AbstractString;
                           poll_interval::Int = 10, max_retries::Int = 10)
    # Query submission
    query_headers = [
        "Content-Type" => "application/json",
        "Authorization" => "Bearer $token",
        "Connection" => "close"
    ]
    query_body = JSON.json(Dict("query" => query))
    query_resp = HTTP.post(
        "$MAQI_API_BASE_URL/api/query",
        query_headers,
        query_body;
        retry_non_idempotent = true,
        retry = true
    )
    query_data = JSON.parse(String(query_resp.body))
    job_id = query_data["job_id"]
    # Poll the job status
    poll_headers = ["Authorization" => "Bearer $token"]
    for i in 1:max_retries
        poll_resp = HTTP.get("$MAQI_API_BASE_URL/api/query/$job_id", poll_headers)
        poll_data = JSON.parse(String(poll_resp.body))["data"]
        status = poll_data["status"]
        if status == "success"
            url = poll_data["result_url"]
            try
                seekstart(io)
                Downloads.download(url, io)
            catch
                @warn "S3 Download connection reset at attempt #$i. \
                    Retrying in $poll_interval seconds..."
                sleep(poll_interval)
                continue
            end
            @info "MAQI job finished with status $status at attempt #$i"
            return nothing
        elseif status in ("partial", "error")
            error("MAQI job failed with status $status: $(poll_data["error_message"])")
        end
        sleep(poll_interval)
    end
    error("Polling timed out after $(max_retries * poll_interval) seconds.")
end

function parse_maqi_result(io::IOBuffer)
    zipreader = ZipReader(take!(io))
    return zip_readentry(zipreader, 1)
end

function main()

    # Parse command line arguments
    parsed_args = parse_commandline()
    email::String = parsed_args["email"]
    password::String = parsed_args["password"]
    query::String = parsed_args["query"]
    output::String = parsed_args["output"]
    poll_interval::Int = parsed_args["poll_interval"]
    max_retries::Int = parsed_args["max_retries"]

    # Authentication
    token, _ = maqi_login(email, password)
    # Post query and fetch result
    io = IOBuffer()
    async_maqi_query!(io, token, query; poll_interval, max_retries)
    # Parse MAQI response
    result = parse_maqi_result(io)

    # Save result to output
    write(output, result)

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end