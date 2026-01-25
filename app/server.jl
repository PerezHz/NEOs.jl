using HTTP, JSON

include(joinpath(dirname(@__DIR__), "scripts", "neocp", "mcmov.jl"))

function handler(req::HTTP.Request)
    try
        if req.method == "POST" && req.target == "/run"
            body = String(req.body)
            input = isempty(body) ? Dict{String, Any}() : JSON.parse(body)
            neos, scout, neoscan, neocp = mcmov(input; write_output = false)
            dict = JSON.json(Dict(
                "ok" => true,
                "neos" => neos,
                "scout" => scout,
                "neoscan" => neoscan,
                "neocp" => neocp
            ))
            return HTTP.Response(200, dict)
        elseif req.method == "GET" && req.target == "/health"
            return HTTP.Response(200, "ok")
        else
            return HTTP.Response(404, "not found")
        end
    catch err
        # Keep errors JSON so clients can handle them
        return HTTP.Response(500, JSON.json(Dict("ok" => false, "error" => sprint(showerror, err))))
    end
end

host = get(ENV, "HOST", "0.0.0.0")
port = parse(Int, get(ENV, "PORT", "8080"))

println("Listening on http://$host:$port")
HTTP.serve(handler, host, port; reuseaddr = true)