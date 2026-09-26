"""
    AdmissibleRegion{T <: Real}

Subset of the topocentric range × range-rate space defined
by the following constraints:
- heliocentric energy ≤ `k_gauss^2/(2a_max)`,
- absolute magnitude ≤ `H_max`,
- geocentric energy ≥ `0`.

# Fields

- `date::DateTime`: time of observation [UTC].
- `ra::T`: right ascension [rad].
- `dec::T`: declination [rad].
- `vra::T`: right ascension velocity [rad/day].
- `vdec::T`: declination velocity [rad/day].
- `mag::T`: apparent magnitude.
- `H_max::T`: maximum absolute magnitude.
- `slope::T`: slope parameter.
- `a_max::T`: maximum semimajor axis [au].
- `ρ_unit/ρ_α/ρ_δ::Vector{T}`: topocentric unit vector and its partials.
- `sun::Vector{T}`: barycentric cartesian state vector of the Sun.
- `observer::Vector{T}`: heliocentric cartesian state vector of observer.
- `coeffs::Vector{T}`: polynomial coefficients.
- `ρ_domain::Vector{T}`: range domain.
- `v_ρ_domain::Vector{T}`: range-rate domain.
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
    mag::T
    H_max::T
    slope::T
    a_max::T
    ρ_unit::Vector{T}
    ρ_α::Vector{T}
    ρ_δ::Vector{T}
    sun::Vector{T}
    observer::Vector{T}
    coeffs::Vector{T}
    ρ_domain::Vector{T}
    v_ρ_domain::Vector{T}
    observatory::ObservatoryMPC{T}
end

# Definition of zero AdmissibleRegion{T}
zero(::Type{AdmissibleRegion{T}}) where {T <: Real} = AdmissibleRegion{T}(
    MINDTTDB, zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T),
    Vector{T}(undef, 0), Vector{T}(undef, 0), Vector{T}(undef, 0),
    Vector{T}(undef, 0), Vector{T}(undef, 0), Vector{T}(undef, 0),
    Vector{T}(undef, 0), Vector{T}(undef, 0), unknownobs(T)
)

iszero(x::AdmissibleRegion{T}) where {T <: Real} = x == zero(AdmissibleRegion{T})

# AdmissibleRegion interface
date(x::AdmissibleRegion) = x.date
ra(x::AdmissibleRegion) = x.ra
dec(x::AdmissibleRegion) = x.dec
vra(x::AdmissibleRegion) = x.vra
vdec(x::AdmissibleRegion) = x.vdec
mag(x::AdmissibleRegion) = x.mag
slopeparameter(x::AdmissibleRegion) = x.slope
observatory(x::AdmissibleRegion) = x.observatory
attributable(x::AdmissibleRegion) = [ra(x), dec(x), vra(x), vdec(x), mag(x)]
rangedomain(x::AdmissibleRegion) = x.ρ_domain
rangeratedomain(x::AdmissibleRegion) = x.v_ρ_domain
numberofcomponents(x::AdmissibleRegion) = 1 + length(rangedomain(x)) > 2

# Print methods for AdmissibleRegion
show(io::IO, x::AdmissibleRegion) = print(io, "Admissible region around ",
    date(x), " at ", observatory(x).name)

function show(io::IO, ::MIME"text/plain", x::AdmissibleRegion)
    t = repeat(' ', 4)
    print(io,
        typeof(x), '\n',
        t, rpad("Observatory: ", 21),  observatory(x).name, '\n',
        t, rpad("Date: ", 21),         date(x), '\n',
        t, rpad("Attributable: ", 21), "[",
            @sprintf("%.5f", rad2deg(ra(x))),   ", ",
            @sprintf("%.5f", rad2deg(dec(x))),  ", ",
            @sprintf("%.5f", rad2deg(vra(x))),  ", ",
            @sprintf("%.5f", rad2deg(vdec(x))), ", ",
            @sprintf("%.2f", mag(x)),
        "]",
    )
    return nothing
end

"""
    AdmissibleRegion(::OpticalTracklet, ::Parameters)

Return the admissible region associated to an optical tracklet. For a list of
parameters, see the `Minimization over the MOV` section of [`Parameters`](@ref).
"""
AdmissibleRegion(x::OpticalTracklet, params::Parameters) = AdmissibleRegion(
    date(x), ra(x), dec(x), vra(x), vdec(x), mag(x), observatory(x), params)

function AdmissibleRegion(date::DateTime, α::T, δ::T, v_α::T, v_δ::T,
                          h::T, observatory::ObservatoryMPC{T},
                          params::Parameters{T}) where {T <: Real}
    # Unpack parameters
    @unpack eph_su, eph_ea, H_max, slope, a_max = params
    # Topocentric unit vector and partials
    ρ, ρ_α, ρ_δ = topounitpdv(α, δ)
    # Time of observation [days since J2000 TDB, Julian days UTC]
    t_days, jd_utc = dtutc2days(date), datetime2julian(date)
    # Barycentric cartesian state vector of the Sun
    sun = eph_su(t_days)
    # Heliocentric cartesian state vector of the observer
    observer = eph_ea(t_days) + kmsec2auday(obsposvelECI(observatory, jd_utc)) - sun
    # Admissible region coefficients
    coeffs = arcoeffs(α, δ, v_α, v_δ, ρ, ρ_α, ρ_δ, observer)
    # Range domain
    ρ_domain = _helrangedomain(coeffs, a_max, h, H_max; slope)
    (isempty(ρ_domain) || ρ_domain[1] > ρ_domain[2]) && return zero(AdmissibleRegion{T})
    # Range-rate domain
    v_ρ_domain = _helrangerates(coeffs, a_max, ρ_domain[1])[1:2]
    # Admissible region
    return AdmissibleRegion{T}(date, α, δ, v_α, v_δ, h, H_max, slope, a_max,
        ρ, ρ_α, ρ_δ, sun, observer, coeffs, ρ_domain, v_ρ_domain, observatory)
end