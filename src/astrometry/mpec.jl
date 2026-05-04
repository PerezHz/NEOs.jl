"""
    fetch_mpec_information(terms; kwargs...)

Return the information given by the Minor Planet Center MPECs API
for the given search `terms`.

# Keyword arguments

- `issued_before::Union{Nothing, DateTime}`: return only those MPECs
    published before this date (default: `nothing`).
- `issued_after::Union{Nothing, DateTime}`: return only those MPECs
    published after this given date (default: `nothing`).

!!! reference
    The Minor Planet Center MPECs API is described at:
    - https://docs.minorplanetcenter.net/mpc-ops-docs/apis/mpecs/
"""
fetch_mpec_information(s...; kwargs...) = fetch_mpec_information(collect(s); kwargs)

function fetch_mpec_information(terms::AbstractVector{<:AbstractString};
                                issued_before::Union{Nothing, DateTime} = nothing,
                                issued_after::Union{Nothing, DateTime} = nothing)
    # Parse HTTP response as String
    text = fetch_http_text(MPC; mode = -2, terms, issued_before, issued_after)
    # Parse JSON
    dict = JSON.parse(text)

    return dict
end