"""
    AdmissibleRegion{T <: Real}

Subset of topocentric range × range-rate space defined by the following constraints:
- heliocentric energy ≤ `k_gauss^2/(2a_max)`,
- absolute magnitude ≤ `H_max`,
- geocentric energy ≥ `0`.

# Fields

- `date::DateTime`: time of observation [UTC].
- `ra::T`: right ascension [rad].
- `dec::T`: declination [rad].
- `vra::T`: right ascension velocity [rad/day].
- `vdec::T`: declination velocity [rad/day].
- `H_max::T`: maximum absolute magnitude.
- `a_max::T`: maximum semimajor axis [au].
- `ρ_unit/ρ_α/ρ_δ::Vector{T}`: topocentric unit vector and its partials.
- `q::Vector{T}`: heliocentric position of observer.
- `coeffs::Vector{T}`: polynomial coefficients.
- `ρ_domain::Vector{T}`: range domain.
- `v_ρ_domain::Vector{T}`: range-rate domain.
- `Fs::Matrix{T}`: boundary points.
- `observatory::ObservatoryMPC{T}`: observing station.

!!! reference
    See Chapter 8 of:
    - https://doi.org/10.1017/CBO9781139175371
    or
    - https://doi.org/10.1007/s10569-004-6593-5
"""
@auto_hash_equals struct AdmissibleRegion{T <: Real}
    date::DateTime
    ra::T
    dec::T
    vra::T
    vdec::T
    H_max::T
    a_max::T
    ρ_unit::Vector{T}
    ρ_α::Vector{T}
    ρ_δ::Vector{T}
    q::Vector{T}
    coeffs::Vector{T}
    ρ_domain::Vector{T}
    v_ρ_domain::Vector{T}
    Fs::Matrix{T}
    observatory::ObservatoryMPC{T}
end

# Definition of zero AdmissibleRegion{T}
zero(::Type{AdmissibleRegion{T}}) where {T <: Real} = AdmissibleRegion{T}(
    DateTime(2000), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T),
    Vector{T}(undef, 0), Vector{T}(undef, 0), Vector{T}(undef, 0),
    Vector{T}(undef, 0), Vector{T}(undef, 0), Vector{T}(undef, 0),
    Vector{T}(undef, 0), Matrix{T}(undef, 0, 0), unknownobs(T)
)

iszero(x::AdmissibleRegion{T}) where {T <: Real} = x == zero(AdmissibleRegion{T})

date(x::AdmissibleRegion) = x.date
ra(x::AdmissibleRegion) = x.ra
dec(x::AdmissibleRegion) = x.dec
vra(x::AdmissibleRegion) = x.vra
vdec(x::AdmissibleRegion) = x.vdec
observatory(x::AdmissibleRegion) = x.observatory

# Print method for AdmissibleRegion
function show(io::IO, x::AdmissibleRegion)
    v = string(
        @sprintf("%.5f", rad2deg(ra(x))), ", ",
        @sprintf("%.5f", rad2deg(dec(x))), ", ",
        @sprintf("%.5f", rad2deg(vra(x))), ", ",
        @sprintf("%.5f", rad2deg(vdec(x))), "",
    )
    print(io, "AE: [", v, "]", " t: ", date(x), " obs: ", observatory(x).name)
end

# Check whether P is inside A's boundary
for U in (:(AbstractVector{T}), :(NTuple{2, T}))
    @eval begin
        function in(P::$U, A::AdmissibleRegion{T}) where {T <: Real}
            @assert length(P) == 2 "Points in admissible region are of dimension 2"
            if A.ρ_domain[1] <= P[1] <= A.ρ_domain[2]
                y_range = rangerates(A, P[1], :outer)
                if length(y_range) == 1
                    return P[2] == y_range[1]
                else
                    return y_range[1] <= P[2] <= y_range[2]
                end
            else
                return false
            end
        end
    end
end