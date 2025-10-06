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
@auto_hash_equals struct VirtualImpactor{T} <: AbstractVirtualImpactor{T}
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
    if any(isnan, fwd.t) || isempty(tvS)
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

function impactingbounds(VA::VirtualAsteroid, σ::Real, ctol::Real)
    a, b = convergence_domain(VA, ctol)
    for CA in Iterators.reverse(VA.CAs)
        a ≤ lbound(CA) || break
        if distance(VA, lbound(CA), ctol) > 0 && lbound(CA) < σ
            a = max(a, lbound(CA))
            break
        end
    end
    for CA in VA.CAs
        ubound(CA) ≤ b || break
        if distance(VA, ubound(CA), ctol) > 0 && σ < ubound(CA)
            b = min(b, ubound(CA))
            break
        end
    end
    return a, b
end

# Val{true}: The radial distance function is monotonic over domain
function impactingcondition(::Val{true}, VA::VirtualAsteroid{T},
                            domain::NTuple{2, T}, ctol::T) where {T <: Real}
    a, b = domain
    a, b = max(a, lbound(VA)), min(b, ubound(VA))
    da, db = distance(VA, a, ctol), distance(VA, b, ctol)
    # The radial distance changes sign in [a, b], so there is a root
    if da * db < 0
        σ = find_zero(σ -> distance(VA, σ, ctol), (a, b), Bisection())
        if da < 0 && db > 0
            domain = (a, σ)
        else
            domain = (σ, b)
        end
        σ = (domain[1] + domain[2]) / 2
    # The radial distance is always negative in [a, b]
    elseif da < 0 && db < 0
        σ = sigma(VA)
        domain = (a, b)
    # The radial distance is always positive in [a, b]
    else
        σ, domain = T(NaN), (T(NaN), T(NaN))
    end
    return σ, domain
end

# Val{false}: The radial distance function is not monotonic around σ
function impactingcondition(::Val{false}, VA::VirtualAsteroid{T},
                            domain::NTuple{2, T}, ctol::T) where {T <: Real}
    σ, d = domain
    if d > 0
        domain = (σ, σ)
    else
        a, b = impactingbounds(VA, σ, ctol)
        σa, da = impactingcondition(Val(true), VA, (a, σ), ctol)
        σb, db = impactingcondition(Val(true), VA, (σ, b), ctol)
        σ = (isnan(σa) || isnan(σb)) ? T(NaN) : σ
        domain = (da[1], db[2])
    end
    return σ, domain
end

"""
    virtualimpactors(VAs, ctol)

Return the indices, sigmas, domains, times and distances of all the
impacting conditions found in a vector of virtual asteroids `VAs`
given a convergence tolerance `ctol`.
"""
function virtualimpactors(VAs::Vector{VirtualAsteroid{T}}, ctol::Real) where {T <: Real}
    # Find all the impacting conditions over VAs
    ics = Vector{Tuple{Int, T, Tuple{T, T}, T, T}}(undef, 0)
    for i in eachindex(VAs)
        VA = VAs[i]
        a, b = convergence_domain(VA, ctol)
        rs = find_zeros(σ -> rvelea(VA, σ, ctol), (a, b), no_pts = max(10, length(VA.CAs)))
        # The radial velocity has no roots in [a, b], so the radial distance is monotonic
        if isempty(rs)
            σ, domain = impactingcondition(Val(true), VA, (a, b), ctol)
            if !isnan(σ)
                t, d = timeofca(VA, σ, ctol), distance(VA, σ, ctol)
                push!(ics, (i, σ, domain, t, d))
            end
        # The radial velocity has at least one root in [a, b], so the radial distance is not
        # monotonic
        else
            for j in eachindex(rs)
                σ = rs[j]
                # Skip non-minima critical points
                concavity(VA, σ, ctol) > 0 || continue
                t, d = timeofca(VA, σ, ctol), distance(VA, σ, ctol)
                σ, domain = impactingcondition(Val(false), VA, (σ, d), ctol)
                if !isnan(σ)
                    push!(ics, (i, σ, domain, t, d))
                end
            end
        end
    end
    # Sort by distance
    sort!(ics, by = Base.Fix2(getindex, 5))
    # Merge duplicated impacting conditions
    Nic = length(ics)
    mask = falses(Nic)
    for k in 1:Nic-1
        ia, σa, domaina, ta, _ = ics[k]
        ib, σb, domainb, tb, _ = ics[k+1]
        # Same index and domain
        if ia == ib && domaina == domainb
            i, domain = ia, domaina
            σ = (domain[1] + domain[2]) / 2
            t, d = timeofca(VAs[i], σ, ctol), distance(VAs[i], σ, ctol)
            ics[k] = (i, σ, domain, t, d)
            mask[k+1] = true
        # Different index, similar sigma and time of impact
        elseif abs( (σa - σb) / σa) < 0.01 && abs((ta - tb) / ta) < 0.01
            mask[k+1] = true
        end
    end
    deleteat!(ics, mask)

    return ics
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
    # Find all the impacting conditions over VAs
    ics = virtualimpactors(VAs, ctol)
    Nic = length(ics)
    # Compute virtual impactors
    VIs = Vector{VirtualImpactor{T}}(undef, 0)
    for k in 1:min(N, Nic)
        i, σ, domain, t, _ = ics[k]
        # Condition is outside the domain of lov
        if !(σ in lov)
            σ = (domain[1] + domain[2]) / 2
            σ in lov || continue
            t = timeofca(VAs[i], σ, ctol)
        end
        VI = VirtualImpactor(lov, od, orbit, params, σ, t, domain)
        # Marginal virtual impactor
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