"""
    VirtualImpactor{T <: Real}

A segment of the line of variations (LOV) that impacts the Earth.

# Fields

- `σ::T`: coordinate on the LOV.
- `t::T`: time of impact [days since J2000 TDB].
- `ip::T`: impact probability.
- `domain::NTuple{2, T}`: segment of the LOV.
- `covariance::Matrix{T}`: target plane covariance matrix.
"""
struct VirtualImpactor{T <: Real}
    σ::T
    t::T
    ip::T
    domain::NTuple{2, T}
    covariance::Matrix{T}
end

function show(io::IO, x::VirtualImpactor)
    t = days2dtutc(x.t)
    ip = @sprintf("%.1E", x.ip)
    print(io, "Virtual impactor at ", t, " with probability ", ip)
end

impact_probability(x::VirtualImpactor) = x.ip

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

function VirtualImpactor(lov::LOV{D, T}, od::AbstractODProblem{D, T},
                         orbit::AbstractOrbit, params::Parameters{T},
                         σ::T, t::T, domain::NTuple{2, T}) where {D, T <: Real}
    # Scalar initial condition
    jd0 = epoch(lov)
    q00 = lov(σ)
    # Second order jet transport initial condition
    q0 = q00 + sigmas(orbit) .* get_variables(T, 2)
    # O-C residuals
    res = NEOs.init_residuals(TaylorN{T}, od, orbit)
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
    nyears = (t + 1 + PE.J2000 - jd0) / yr
    fwd, tvS, _, _ = NEOs.propagate_root(lov.dynamics, q0, epoch(lov), nyears, params)
    # Close approach
    t_CA = fwd.t0 + tvS[end]
    xae = fwd(t_CA) - params.eph_ea(t_CA)
    # Asteroid's geocentric semimajor axis
    a = semimajoraxis(xae..., PE.μ[ea], zero(T))
    if a < 0
        # Earth's heliocentric state vector
        xes = params.eph_ea(t_CA) - params.eph_su(t_CA)
        # Öpik's coordinates
        B = bopik(xae, xes)
        X, Y, Z = B.ξ, B.ζ, B.b
    else
        # Modified target plane
        X, Y = mtp(xae)
        Z = 1.0 * one(x)
    end
    # Target plane covariance matrix at close approach
    Γ_B = project([X, Y], zeros(6), Γ)
    # Impact probability
    ρ, x_C, y_C = zero(T), zero(T), zero(T)
    m_X, m_Y, r = cte(X), cte(Y), cte(Z)
    σ_X, σ_Y = sqrt.(diag(Γ_B))
    iparams = (r, ρ, m_X, m_Y, σ_X, σ_Y, x_C, y_C)
    if iszero(width(domain))
        a, b = mangelbounds(iparams)
        ip = ( exp(-σ^2/2) / (2 * sqrt(2π)) ) * quadgk(Base.Fix2(mangelintegrand, iparams), a, b)[1]
    else
        ip = impact_probability(domain[1], domain[2])
    end

    return VirtualImpactor{T}(σ, t, ip, domain, Γ_B)
end

function virtualimpactors(lov::LOV{D, T}, od::AbstractODProblem{D, T},
                          orbit::AbstractOrbit, params::Parameters{T},
                          VAs::Vector{VirtualAsteroid{T}}, ϵ::T;
                          pmin::T = 1E-10) where {D, T <: Real}
    # Find all the minima of distance over VAs
    is = Vector{Int}(undef, 0)
    σs = Vector{T}(undef, 0)
    ts = Vector{T}(undef, 0)
    ds = Vector{T}(undef, 0)
    for i in eachindex(VAs)
        VA = VAs[i]
        rs = find_zeros(σ -> rvelea(VA, σ, ϵ), convergence_domain(VA, ϵ),
                        no_pts = max(10, length(VA.CAs)))
        for j in eachindex(rs)
            (lov.domain[1] ≤ rs[j] ≤ lov.domain[2] && concavity(VA, rs[j], ϵ) > 0) || continue
            push!(is, i)
            push!(σs, rs[j])
            push!(ts, time(VA, rs[j], ϵ))
            push!(ds, distance(VA, rs[j], ϵ))
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
    # Compute virtual impactors
    VIs = Vector{VirtualImpactor{T}}(undef, 0)
    k = 1
    for (i, σ, t, d) in zip(is, σs, ts, ds)
        VA = VAs[i]
        if d < 0
            a, b = convergence_domain(VA, ϵ)
            if lbound(VA[1]) < σ
                a = max(a, lbound(VA[findlast(CA -> distance(CA, lbound(CA)) > 0 && lbound(CA) < σ, VA.CAs)]))
            end
            if σ < ubound(VA[end])
                b = min(b, ubound(VA[findfirst(CA -> distance(CA, ubound(CA)) > 0 && σ < ubound(CA), VA.CAs)]))
            end
            a = find_zero(x -> distance(VA, x, ϵ), (a, σ), Bisection())
            b = find_zero(x -> distance(VA, x, ϵ), (σ, b), Bisection())
            domain = (a, b)
        else
            domain = (σ, σ)
        end
        @time VI = VirtualImpactor(lov, od, orbit, params, σ, t, domain)
        @show VI
        push!(VIs, VI)
        # impact_probability(VIs[end]) < pmin && break
        k == 5 && break
        k += 1
    end
    # Sort by time of impact
    sort!(VIs, by = x -> x.t)

    return VIs
end

function impactor_table(VIs::Vector{VirtualImpactor{T}}) where {T <: Real}
    s1 = string(
        "Impactor table{$T}\n",
        repeat('-', 56), "\n",
        "Date (UTC)                 Sigma      Impact probability\n",
    )
    ts = [rpad(string(days2dtutc(VI.t)), 27) for VI in VIs]
    σs = [rpad(@sprintf("%+.4f", VI.σ), 11) for VI in VIs]
    ips = [rpad(@sprintf("%+.2E", VI.ip), 11) for VI in VIs]
    s2 = Vector{String}(undef, length(VIs))
    for i in eachindex(VIs)
        s2[i] = string(
            ts[i], σs[i], ips[i], "\n"
        )
    end
    println(s1, join(s2))
end