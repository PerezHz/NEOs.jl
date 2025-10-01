"""
    VirtualImpactor{T} <: AbstractVirtualImpactor{T}

A segment of the line of variations that impacts the Earth.

# Fields

- `t::T`: time of impact [days since J2000 TDB].
- `σ::T`: coordinate on the line of variations.
- `ip::T`: impact probability.
- `a::T`: geocentric semimajor axis [au].
- `domain::NTuple{2, T}`: segment of the line of variations.
- `covariance::Matrix{T}`: target plane covariance matrix.
"""
struct VirtualImpactor{T} <: AbstractVirtualImpactor{T}
    t::T
    σ::T
    ip::T
    a::T
    domain::NTuple{2, T}
    covariance::Matrix{T}
end

nominaltime(x::VirtualImpactor) = x.t
date(x::VirtualImpactor) = days2dtutc(nominaltime(x))

sigma(x::VirtualImpactor) = x.σ

impact_probability(x::VirtualImpactor) = x.ip

semimajoraxis(x::VirtualImpactor) = x.a
function vinf(x::VirtualImpactor{T}) where {T <: Real} # km/s
    if ishyperbolic(x)
        return sqrt( PE.μ[ea] / (-semimajoraxis(x)) ) * (au/daysec)
    else
        return zero(T)
    end
end

width(x::VirtualImpactor) = width(x.domain)

isoutlov(x::VirtualImpactor) = iszero(width(x))
ishyperbolic(x::VirtualImpactor) = semimajoraxis(x) < 0
ismarginal(x::VirtualImpactor) = isempty(x.covariance)

"""
    impactenergy(VI, orbit, params; kwargs...)

Return the impact energy [Mt] of a virtual impactor `VI` associated
to an `orbit`.

# Keyword arguments

- `ρ::Real`: density [kg/m³] (default: `2_600`).
- `p::Real`: albedo (default: `0.14`).

!!! reference
    - https://doi.org/10.1006/icar.2002.6910
"""
function impactenergy(VI::VirtualImpactor, orbit::AbstractOrbit, params::Parameters;
                      ρ::Real = 2_600, p::Real = 0.14)
    # Object's mass [kg]
    M = mass(orbit, params, ρ, p)
    # Impact velocity^2 [km^2/s^2]
    V2 = vinf(VI)^2 + EARTH_ESCAPE_VELOCITY^2
    # Impact energy [Mt]
    return 0.5 * M * V2 / 4.184E+9
end

"""
    palermoscale(VI, orbit, params; kwargs)

Return the Palermo Scale of a virtual impactor `VI` associated to an `orbit`.

# Keyword arguments

- `ρ::Real`: density [kg/m³] (default: `2_600`).
- `p::Real`: albedo (default: `0.14`).

!!! reference
    - https://doi.org/10.1006/icar.2002.6910
"""
function palermoscale(VI::VirtualImpactor, orbit::AbstractOrbit, params::Parameters;
                      ρ::Real = 2_600, p::Real = 0.14)
    # Impact probability
    IP = impact_probability(VI)
    # Background impact frequency [yr⁻¹]
    fb = 0.03 * impactenergy(VI, orbit, params; ρ, p)^(-0.8)
    # Time until impact [yr]
    ΔT = (nominaltime(VI) - epoch(orbit)) / yr
    # Palermo scale
    return log10(IP / (fb * ΔT))
end

"""
    torinoscale(VI, orbit, params; kwargs)

Return the Torino Scale of a virtual impactor `VI` associated to an `orbit`.

# Keyword arguments

- `ρ::Real`: density [kg/m³] (default: `2_600`).
- `p::Real`: albedo (default: `0.14`).

!!! reference
    - https://doi.org/10.1016/S0032-0633(00)00006-4
"""
function torinoscale(VI::VirtualImpactor, orbit::AbstractOrbit, params::Parameters;
                     ρ::Real = 2_600, p::Real = 0.14)
    # Impact probability
    IP = impact_probability(VI)
    # Impact energy [Mt]
    E = impactenergy(VI, orbit, params; ρ, p)
    # Logarithms
    logE, logIP = log10(E), log10(IP)
    # Torino scale
    if (logE+1)/3 + (logIP+2)/2 < 0 || logE < 0
        return 0
    elseif logIP < -2 && logE ≥ 0
        if (logE+1)/3 + (logIP+2)/2 ≥ 0 && (logE-2)/3 + (logIP+2)/2 < 0
            return 1
        elseif (logE-2)/3 + (logIP+2)/2 ≥ 0 && (logE-5)/3 + (logIP+2)/2 < 0
            return 2
        elseif (logE-5)/3 + (logIP+2)/2 ≥ 0 && logIP < -2
            return 6
        end
    elseif -2 ≤ logIP && IP < 0.99
        if 0 ≤ logE < 2
            return 3
        elseif logE ≥ 2 && (logE-5)/3 + (logIP+2)/2 < 0
            return 4
        elseif logE < 5 && (logE-5)/3 + (logIP+2)/2 ≥ 0
            return 5
        elseif logE ≥ 5
            return 7
        end
    elseif IP ≥ 0.99
        if 0 ≤ logE < 2
            return 8
        elseif 2 ≤ logE < 5
            return 9
        elseif logE ≥ 5
            return 10
        end
    end
end

impact_probability(a::Real, b::Real) = erf(a/sqrt(2), b/sqrt(2)) / 2

function Splus(x::T, params::NTuple{8, T}) where {T <: Real}
    r, ρ, m_X, m_Y, σ_X, σ_Y, x_C, y_C = params
    C = 1 / sqrt(2 * (1 - ρ^2))
    A = (y_C + sqrt(r^2 - (x - x_C)^2) - m_Y) / σ_Y
    B = ρ * ( (x - m_X) / σ_X)
    return C * (A - B)
end

function Sminus(x::T, params::NTuple{8, T}) where {T <: Real}
    r, ρ, m_X, m_Y, σ_X, σ_Y, x_C, y_C = params
    C = 1 / sqrt(2 * (1 - ρ^2))
    A = (y_C - sqrt(r^2 - (x - x_C)^2) - m_Y) / σ_Y
    B = ρ * ( (x - m_X) / σ_X)
    return C * (A - B)
end

function mangelbounds(params::NTuple{8, T}) where {T <: Real}
    r, m_X, σ_X, x_C = params[1], params[3], params[5], params[7]
    return (x_C - r - m_X) / σ_X, (x_C + r - m_X) / σ_X
end

function mangelintegrand(v::T, params::NTuple{8, T}) where {T <: Real}
    m_X, σ_X = params[3], params[5]
    return exp(-v^2/2) * erf(Sminus(σ_X * v + m_X, params), Splus(σ_X*v + m_X, params))
end

function impact_probability(σ::T, params::NTuple{8, T}) where {T <: Real}
    C = exp(-σ^2/2) / (2 * sqrt(2π))
    a, b = mangelbounds(params)
    return C * quadgk(Base.Fix2(mangelintegrand, params), a, b)[1]
end

function VirtualImpactor(lov::LineOfVariations{D, T}, od::AbstractODProblem{D, T},
                         orbit::AbstractOrbit, params::Parameters{T},
                         σ::T, t::T, domain::NTuple{2, T}) where {D, T <: Real}
    # Set jet transport variables
    Npar = numvars(orbit)
    set_od_order(T, 2, Npar)
    # Scalar initial condition
    jd0 = epoch(lov) + PE.J2000
    q00 = lov(σ)
    # Second order jet transport initial condition
    q0 = q00 + sigmas(orbit) .* get_variables(T, 2)
    # O-C residuals
    res = init_residuals(TaylorN{T}, od, orbit)
    propres!(res, od, q0, jd0, params)
    # Covariance matrix at reference epoch
    Q = nms(res)
    C = notout(res) * TS.hessian(Q)
    Γ = inv(C)
    # Chi parameter
    # χ = sqrt(notoutobs(orbit) * ( cte(Q) - nms(orbit) ))
    # First order jet transport initial condition
    q0 = q00 + sigmas(orbit) .* get_variables(T, 1)
    # Forward propagation
    nyears = (t + 2 + PE.J2000 - jd0) / yr
    fwd, tvS, _, _ = propagate_root(lov.dynamics, q0, jd0, nyears, params)
    # Marginal close approach
    if any(isnan, fwd.t)
        Γ_B = Matrix{T}(undef, 0, 0)
        ip, a = T(NaN), T(NaN)
        return VirtualImpactor{T}(t, σ, ip, a, domain, Γ_B)
    end
    # Close approach
    t_CA = fwd.t0 + tvS[end]
    xae = fwd(t_CA) - params.eph_ea(t_CA)
    # Asteroid's geocentric semimajor axis
    a = semimajoraxis(xae..., PE.μ[ea], zero(T))
    if a < 0
        # Earth's heliocentric state vector
        xes = params.eph_ea(t_CA) - params.eph_su(t_CA)
        # Öpik's coordinates
        tp = bopik(xae, xes)
    else
        # Modified target plane
        tp = mtp(xae)
    end
    # Target plane coordinates
    X, Y, Z = targetplane(tp)
    # Target plane covariance matrix at close approach
    Γ_B = project([X, Y], zeros(6), Γ)
    # Impact probability
    ρ, x_C, y_C = zero(T), zero(T), zero(T)
    m_X, m_Y, r = cte(X), cte(Y), cte(Z)
    σ_X, σ_Y = sqrt.(diag(Γ_B))
    iparams = (r, ρ, m_X, m_Y, σ_X, σ_Y, x_C, y_C)
    if iszero(width(domain))
        ip = impact_probability(σ, iparams)
    else
        ip = impact_probability(domain[1], domain[2])
    end

    return VirtualImpactor{T}(t, σ, ip, cte(a), domain, Γ_B)
end

"""
    virtualimpactors(VAs, ctol)

Return the indices, sigmas, times and distances of all the minima of [`distance`](@ref)
found in a vector of virtual asteroids `VAs` given a convergence tolerance `ctol`.
"""
function virtualimpactors(VAs::Vector{VirtualAsteroid{T}}, ctol::Real) where {T <: Real}
    # Find all the minima of distance over VAs
    is = Vector{Int}(undef, 0)
    σs = Vector{T}(undef, 0)
    ts = Vector{T}(undef, 0)
    ds = Vector{T}(undef, 0)
    for i in eachindex(VAs)
        VA = VAs[i]
        rs = find_zeros(σ -> rvelea(VA, σ, ctol), convergence_domain(VA, ctol),
                        no_pts = max(10, length(VA.CAs)))
        for j in eachindex(rs)
            # Skip non-minima critical points
            concavity(VA, rs[j], ctol) > 0 || continue
            push!(is, i)
            push!(σs, rs[j])
            push!(ts, timeofca(VA, rs[j], ctol))
            push!(ds, distance(VA, rs[j], ctol))
        end
    end
    # Sort by distance
    perm = sortperm(ds)
    permute!(is, perm)
    permute!(σs, perm)
    permute!(ts, perm)
    permute!(ds, perm)
    # Merge duplicated VIs
    diffi = @. abs( (is[1:end-1] - is[2:end]) / is[1:end-1] )
    diffσ = @. abs( (σs[1:end-1] - σs[2:end]) / σs[1:end-1] )
    difft = @. abs( (ts[1:end-1] - ts[2:end]) / ts[1:end-1] )
    idxs = findall(@. diffi > 0 && diffσ < 0.01 && difft < 0.01) .+ 1
    deleteat!(is, idxs)
    deleteat!(σs, idxs)
    deleteat!(ts, idxs)
    deleteat!(ds, idxs)

    return is, σs, ts, ds
end

function impactinglbound(VA::VirtualAsteroid, σ::Real, ctol::Real)
    for (i, CA) in Iterators.reverse(enumerate(VA.CAs))
        distance(VA, lbound(CA), ctol) > 0 && lbound(CA) < σ && return lbound(VA[i])
    end
    return lbound(VA[1])
end

function impactingubound(VA::VirtualAsteroid, σ::Real, ctol::Real)
    for (i, CA) in enumerate(VA.CAs)
        distance(VA, ubound(CA), ctol) > 0 && σ < ubound(CA) && return ubound(VA[i])
    end
    return ubound(VA[end])
end

function impactingdomain(VA::VirtualAsteroid, ctol::Real, σ::Real, d::Real)
    if d < 0
        a, b = convergence_domain(VA, ctol)
        if lbound(VA[1]) < σ
            a = max(a, impactinglbound(VA, σ, ctol))
        end
        if σ < ubound(VA[end])
            b = min(b, impactingubound(VA, σ, ctol))
        end
        a = find_zero(x -> distance(VA, x, ctol), (a, σ), Bisection())
        b = find_zero(x -> distance(VA, x, ctol), (σ, b), Bisection())
        domain = (a, b)
    else
        domain = (σ, σ)
    end
    return domain
end

"""
    virtualimpactors(VAs, ctol, lov, od, orbit, params; kwargs...)

Return the virtual impactors found in a vector of virtual asteroids `VAs` given
a convergence tolerance `ctol`. A line of variations, orbit determination problem,
orbit and parameters are needed to compute the target plane covariance matrix.

# Keyword arguments

- `N::Int`: maximum number of virtual impactors to compute (default: `10`).
"""
function virtualimpactors(VAs::Vector{VirtualAsteroid{T}}, ctol::Real,
                          lov::LineOfVariations{D, T}, od::AbstractODProblem,
                          orbit::AbstractOrbit, params::Parameters;
                          N::Int = 10) where {D, T <: Real}
    # Find all minima of distance over VAs
    is, σs, ts, ds = virtualimpactors(VAs, ctol)
    # Eliminate conditions outside the domain of lov
    mask = @. !($lbound(lov) ≤ σs ≤ $ubound(lov))
    deleteat!(is, mask)
    deleteat!(σs, mask)
    deleteat!(ts, mask)
    deleteat!(ds, mask)
    # Compute virtual impactors
    VIs = Vector{VirtualImpactor{T}}(undef, 0)
    for k in 1:min(N, length(ds))
        i, σ, t, d = is[k], σs[k], ts[k], ds[k]
        domain = impactingdomain(VAs[i], ctol, σ, d)
        VI = VirtualImpactor(lov, od, orbit, params, σ, t, domain)
        if !ismarginal(VI)
            push!(VIs, VI)
        end
    end
    # Sort by time of impact
    sort!(VIs, by = nominaltime)

    return VIs
end

function summary(VIs::AbstractVector{VirtualImpactor{T}}) where {T <: Real}
    s1 = string(
        "Impactor table{$T}\n",
        repeat('-', 56), "\n",
        "Date (UTC)                 Sigma      Impact probability\n",
    )
    s2 = Vector{String}(undef, length(VIs))
    for (i, VI) in enumerate(VIs)
        t = rpad(string(date(VI)), 27)
        σ = rpad(@sprintf("%+.4f", sigma(VI)), 11)
        asterisk = isoutlov(VI) ? " *" : "  "
        ip = rpad(@sprintf("%+.2E", VI.ip) * asterisk, 11)
        s2[i] = string(t, σ, ip, "\n")
    end
    return string(s1, join(s2))
end

impactor_table(x::AbstractVector{VirtualImpactor{T}}) where {T <: Real} = println(summary(x))