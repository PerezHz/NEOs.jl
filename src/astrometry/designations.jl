"""
    unpacknum(::AbstractString)

Unpack a Minor Planet Center permanent designation.

See also [`packnum`](@ref).

!!! reference
    The Minor Planet Center unpacking of permanent designations is described at:
    - https://minorplanetcenter.net/iau/info/PackedDes.html#perm
"""
function unpacknum(s::AbstractString)
    # Packed designation corresponds to a comet
    isnumeric(first(s)) && isletter(last(s)) && return strip(s, '0')
    # Packed designation corresponds to a minor planet and...
    # number is less than 100,000
    if all(isnumeric, s)
        return lstrip(s, '0')
    # number is between 100,000 and 619,999
    elseif isletter(s[1])
        num = 100_000 + 10_000 * (isuppercase(s[1]) ? s[1] - 'A' :
            s[1] - 6 - 'A') + parse(Int, s[2:end])
        return string(num)
    # number is above 619,999
    elseif s[1] == '~'
        num = 620_000 + sum((findfirst(c, BASE_62_ENCODING) - 1) * 62^(4-i)
            for (i, c) in enumerate(s[2:5]))
        return string(num)
    else
        return throw(ArgumentError("Invalid MPC packed permanent designation"))
    end
end

"""
    packnum(::AbstractString)

Pack a Minor Planet Center permanent designation.

See also [`unpacknum`](@ref).

!!! reference
    The Minor Planet Center packing of permanent designations is described at:
    - https://minorplanetcenter.net/iau/info/PackedDes.html#perm
"""
function packnum(s::AbstractString)
    # Packed designation corresponds to a comet
    isletter(last(s)) && return return lpad(s, 5, '0')
    # Packed designation corresponds to a minor planet and...
    # number is less than 100,000
    num = parse(Int, s)
    if num < 100_000
        return lpad(num, 5, '0')
    # number is between 100,000 and 619,999
    elseif 100_000 ≤ num < 619_999
        num = num - 100_000
        q = num ÷ 10_000
        s1 = q < 26 ? q + 'A' : q + 'a' - 26
        s2 = lpad(rem(num, 10_000), 4, '0')
        return string(s1, s2)
    # number is above 619,999
    elseif 620_000 ≤ num
        num = num - 620_000
        d = digits(num, base = 62) .+ 1
        encoding = reverse(BASE_62_ENCODING[d])
        return string("~", lpad(encoding, 4, '0'))
    else
        return throw(ArgumentError("Invalid MPC permanent designation"))
    end
end

"""
    unpackdesig(::AbstractString)

Unpack a Minor Planet Center provisional designation.

See also [`packdesig`](@ref).

!!! reference
    The Minor Planet Center unpacking of provisional designations is described at:
    - https://minorplanetcenter.net/iau/info/PackedDes.html#prov
"""
function unpackdesig(s::AbstractString)
    if s[1] == 'I'
        century = 1800
    elseif s[1] == 'J'
        century = 1900
    elseif s[1] == 'K'
        century = 2000
    end
    year = century + parse(Int, s[2:3])
    half_month_letter = s[4]
    second_letter = s[7]
    if all(isnumeric, s[5:6])
        cycle_count = parse(Int, s[5:6])
    elseif isuppercase(s[5])
        cycle_count = 100 + 10 * (s[5] - 'A') + parse(Int, s[6])
    else
        cycle_count = 100 + 10 * (s[5] - 6 - 'A') + parse(Int, s[6])
    end
    cycle_count_s = iszero(cycle_count) ? "" : string(cycle_count)

    return string(year, " ", half_month_letter, second_letter,
        cycle_count_s)
end

"""
    packdesig(::AbstractString)

Pack a Minor Planet Center provisional designation.

See also [`unpackdesig`](@ref).

!!! reference
    The Minor Planet Center packing of provisional designations is described at:
    - https://minorplanetcenter.net/iau/info/PackedDes.html#prov
"""
function packdesig(s::AbstractString)
    year, tail = split(s, " ")
    if year[1:2] == "18"
        packed_year = string("I", year[3:4])
    elseif year[1:2] == "19"
        packed_year = string("J", year[3:4])
    elseif year[1:2] == "20"
        packed_year = string("K", year[3:4])
    end
    half_month_letter, second_letter = tail[1:2]
    number = length(tail) > 2 ? parse(Int, tail[3:end]) : 0
    if number ≤ 99
        cycle_count = lpad(tail[3:end], 2, "0")
    elseif 100 ≤ number ≤ 359
        cycle_count = string(parse(Int, tail[3:4]) - 10 + 'A', tail[end])
    else
        cycle_count = string(parse(Int, tail[3:4]) - 36 + 'a', tail[end])
    end

    return string(packed_year, half_month_letter, cycle_count, second_letter)
end

"""
    fetch_designation_information(ids)

Return the information given by the Minor Planet Center designation
identifier API about the designation(s) `ids`.

!!! reference
    The Minor Planet Center designation identifier API is described at:
    - https://www.minorplanetcenter.net/mpcops/documentation/designation-identifier-api/
"""
fetch_designation_information(s...) = fetch_designation_information(collect(s))

function fetch_designation_information(ids::AbstractVector{<:AbstractString})
    # Parse HTTP response as String
    text = fetch_http_text(MPC; mode = -1, ids)
    # Parse JSON
    dict = JSON.parse(text)

    return dict
end