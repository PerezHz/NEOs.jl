"""
    VirtualImpactor{T} <: AbstractVirtualImpactor{T}

A segment of the line of variations that impacts the Earth.

# Fields

- `t::T`: time of impact [days since J2000 TDB].
- `σ::T`: LOV index.
- `ip::T`: impact probability.
- `a::T`: planetocentric semimajor axis [au].
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

# Outer constructorss
function VirtualImpactor(t::T, σ::T, ip::T, a::T, domain::NTuple{2, T},
                         Γ_tp::Matrix{T} = Matrix{T}(undef, 0, 0)) where {T <: Real}
    return VirtualImpactor{T}(t, σ, ip, a, domain, Γ_tp)
end

function VirtualImpactor(RT::ReturnT1{T}, σ::T, domain::NTuple{2, T},
                         ctol::T) where {T <: Real}
    t = timeofca(RT, σ, ctol)
    ip = impact_probability(domain[1], domain[2])
    a = semimajoraxis(RT, σ, ctol)
    return VirtualImpactor(t, σ, ip, a, domain)
end

# Abbreviations
const Shower{T, U} = AbstractVector{Return{T, U}}
const ShowerT1{T} = Shower{T, Taylor1{T}}

# AbstractVirtualImpactor interface
sigma(x::VirtualImpactor) = x.σ
nominaltime(x::VirtualImpactor) = x.t
semimajoraxis(x::VirtualImpactor) = x.a
covariance(x::VirtualImpactor) = x.covariance
impact_probability(x::VirtualImpactor) = x.ip

date(x::VirtualImpactor) = days2dtutc(nominaltime(x))

width(x::VirtualImpactor) = width(x.domain)

isoutlov(x::VirtualImpactor) = iszero(width(x))
ishyperbolic(x::VirtualImpactor) = semimajoraxis(x) < 0

# A marginal virtual impactor is one that was detected
# but could not be verified
ismarginal(x::VirtualImpactor) = isempty(covariance(x))

semiwidth(x::VirtualImpactor{T}) where {T <: Real} =
    ismarginal(x) ? T(NaN) : sqrt(first(eigvals(covariance(x))))

stretching(x::VirtualImpactor{T}) where {T <: Real} =
    ismarginal(x) ? T(NaN) : sqrt(last(eigvals(covariance(x))))

# Impactor table
function summary(VIs::AbstractVector{VirtualImpactor{T}}) where {T <: Real}
    s1 = string(
        "Impactor table{$T}\n",
        repeat('-', 80), "\n",
        "Date (UTC)          Sigma      Semi-width [RE]    Stretching [RE]    IP\n",
    )
    s2 = Vector{String}(undef, length(VIs))
    for (i, VI) in enumerate(VIs)
        t = rpad(Dates.format(date(VI), "yyyy-mm-dd HH:MM"), 20)
        σ = rpad(@sprintf("%+.4f", sigma(VI)), 11)
        w = rpad(@sprintf("%.3f", semiwidth(VI)), 19)
        Λ = rpad(@sprintf("%.2E", stretching(VI)), 19)
        asterisk = isoutlov(VI) ? " *" : "  "
        ip = rpad(@sprintf("%.2E", VI.ip) * asterisk, 11)
        s2[i] = string(t, σ, w, Λ, ip, "\n")
    end
    return string(s1, join(s2))
end

impactor_table(x::AbstractVector{VirtualImpactor{T}}) where {T} = println(summary(x))

# Impact probability computation
impact_probability(a::Real, b::Real) = erf(a/sqrt(2), b/sqrt(2)) / 2

function impact_probability(params::NTuple{6, T}) where {T <: Real}
    w = params[5]
    C = 1 / (2w*sqrt(2π))
    rmin, rmax = milani2005bounds(params)
    if rmin < rmax
        return C * quadgk(Base.Fix2(milani2005integrand, params), rmin, rmax)[1]
    else
        return zero(T)
    end
end

function milani2005bounds(params::NTuple{6, T}) where {T <: Real}
    _X_, _Z_, w = params[1], params[3], params[5]
    rmin = max(-_Z_, _X_ - 8w)
    rmax = min(+_Z_, _X_ + 8w)
    return rmin, rmax
end

function milani2005integrand(u::T, params::NTuple{6, T}) where {T <: Real}
    _X_, _Y_, _Z_, σ, w, Λ = params
    A = (u - _X_) / w
    C = (_Y_ - Λ*σ) / (sqrt(2)*Λ)
    Splus  = C + sqrt(_Z_^2 - u^2) / (sqrt(2)*Λ)
    Sminus = C - sqrt(_Z_^2 - u^2) / (sqrt(2)*Λ)
    return exp(-A^2/2) * erf(Sminus, Splus)
end

# Virtual impactors search
function VirtualImpactor(
        IM::AbstractIMProblem{D, T}, lov::LineOfVariations{D, T},
        VI::VirtualImpactor{T}, params::Parameters{T}
    ) where {D, T <: Real}
    # Unpack
    @unpack orbit, target = IM
    @unpack t, σ, ip, a, domain = VI
    # Set jet transport variables
    Npar = numvars(orbit)
    set_od_order(T, 2, Npar)
    # Scalar initial condition
    jd0 = epoch(lov) + PE.J2000
    q00 = lov(σ)
    # Second order jet transport initial condition
    q0 = q00 + sigmas(orbit) .* get_variables(T, 2)
    # O-C residuals
    res = init_optical_residuals(TaylorN{T}, IM)
    propres!(res, IM, q0, jd0, params)
    # Covariance matrix at reference epoch
    Q = nms(res)
    C = notout(res) * TS.hessian(Q)
    Γ = inv(C)
    # Chi parameter
    # χ = sqrt(notoutobs(orbit) * ( cte(Q) - nms(orbit) ))
    # First order jet transport initial condition
    q0 = q00 + sigmas(orbit) .* get_variables(T, 1)
    # Virtual asteroid
    VA = VirtualAsteroid(epoch(lov), σ, domain, q0)
    # Number of years until impact
    nyears = min(t + 2 + PE.J2000 - jd0, dtutc2days(DateTime(2099, 12, 31))) / yr
    # Close approach
    CAs = closeapproaches(IM, VA, nyears, params)
    if isempty(CAs)
        @warn string(
            "The following virtual impactor was detected but could not be verified\n",
            "Date (UTC)          Sigma      IP\n",
            rpad(Dates.format(days2dtutc(t), "yyyy-mm-dd HH:MM"), 20),
            rpad(@sprintf("%+.4f", σ), 11),
            @sprintf("%.2E", ip),
            width(domain) > 0 ? "" : " *"
        )
        return VI
    end
    CA = CAs[end]
    # Target plane coordinates
    X, Y, Z = CA.x, CA.y, CA.z
    # Target plane covariance matrix at close approach
    Γ_tp = Symmetric(project([X, Y], zeros(Npar), Γ))
    # Eigenpairs of the target plane covariance matrix
    E_tp = eigen(Γ_tp)
    if any(<(0), E_tp.values)
        @warn string(
            "The following virtual impactor was detected but the target plane \
            covariance matrix was not positive definite\n",
            "Date (UTC)          Sigma      IP\n",
            rpad(Dates.format(days2dtutc(t), "yyyy-mm-dd HH:MM"), 20),
            rpad(@sprintf("%+.4f", σ), 11),
            @sprintf("%.2E", ip),
            width(domain) > 0 ? "" : " *"
        )
        return VI
    end
    # Approximate the impact probability using the formula from Milani et al (2005)
    if iszero(width(domain))
        # Semi-width and stretching
        w, Λ = sqrt.(E_tp.values)
        # Angle between the semimajor axis of the TP covariance ellipse and the Y-axis
        α = angle(E_tp.vectors[:, 2], [zero(T), one(T)])
        # Rotate TP coordinates by an angle of -α
        _X_ = cte(X) * cos(-α) - cte(Y) * sin(-α)
        _Y_ = cte(X) * sin(-α) + cte(Y) * cos(-α)
        _Z_ = cte(Z)
        # Compute the 2D integral
        iparams = (_X_, _Y_, _Z_, σ, w, Λ)
        ip = impact_probability(iparams)
    end

    return VirtualImpactor{T}(t, σ, ip, a, domain, Γ_tp)
end

function virtualimpactors(RT::ReturnT1{T}; ctol::Real = T(Inf), σmax::Real = 3.0,
                          no_pts::Int = 100, dmax::Real = zero(T)) where {T <: Real}
    # Allocate memory
    VIs = Vector{VirtualImpactor{T}}(undef, 0)
    # Intersect the convergence domain of RT with (-σmax, σmax)
    a, b = convergence_domain(RT, ctol)
    overlap((a, b), (-σmax, σmax)) || return VIs
    a, b = max(a, -σmax), min(b, σmax)
    # Find the roots of the radial velocity
    rs = find_zeros(σ -> radialvelocity(RT, σ, ctol), (a, b); no_pts)
    # Check if any of the roots of the radial velocity is a local minimum
    # with positive distance
    for σ in rs
        VI = VirtualImpactor(RT, σ, (σ, σ), ctol)
        if 0 ≤ nominaltime(VI) ≤ 36525.0 && 0 < distance(RT, σ, ctol) < dmax &&
            0 < concavity(RT, σ, ctol)
            push!(VIs, VI)
        end
    end
    # Search for virtual impactors between the roots of the radial velocity
    # Since there are no critical points in domain, the radial distance is monotonic
    for domain in zip(Iterators.flatten((a, rs)), Iterators.flatten((rs, b)))
        da, db = distance(RT, domain[1], ctol), distance(RT, domain[2], ctol)
        # The radial distance is always positive in domain
        (da > 0 && db > 0) && continue
        # The radial distance changes sign in domain, so there is a root
        if da * db < 0
            σ = find_zero(σ -> distance(RT, σ, ctol), domain, Bisection())
            domain = (da < 0 && db > 0) ? (domain[1], σ) : (σ, domain[2])
        end
        σ = midpoint(domain)
        VI = VirtualImpactor(RT, σ, domain, ctol)
        if 0 ≤ nominaltime(VI) ≤ 36525.0 && distance(RT, σ, ctol) < dmax
            push!(VIs, VI)
        end
    end
    isempty(VIs) && return VIs
    # Sort by domain
    sort!(VIs, by = Base.Fix2(getfield, :domain))
    # Merge virtual impactors
    newVIs = [VIs[1]]
    for VI in Iterators.drop(VIs, 1)
        domaina, domainb = newVIs[end].domain, VI.domain
        # The intersection of domaina and domainb is not empty
        if overlap(domaina, domainb)
            domain = (min(domaina[1], domainb[1]), max(domaina[2], domainb[2]))
            σ = midpoint(domain)
            newVIs[end] = VirtualImpactor(RT, σ, domain, ctol)
        # The intersection of domaina and domainb is empty
        else
            push!(newVIs, VI)
        end
    end

    return newVIs
end

function virtualimpactors(RTs::ShowerT1{T}; ctol::Real = T(Inf), σmax::Real = 3.0,
                          no_pts::Int = 100, dmax::Real = zero(T)) where {T <: Real}
    # Search for virtual impactors in each return
    _VIs_ = virtualimpactors.(RTs; ctol, σmax, no_pts, dmax)
    ks = reduce(vcat, fill.(eachindex(_VIs_), length.(_VIs_)))
    VIs = reduce(vcat, _VIs_)
    isempty(VIs) && return VIs
    # Sort by domain
    perm = sortperm(VIs, by = Base.Fix2(getfield, :domain))
    permute!(ks, perm)
    permute!(VIs, perm)
    # Merge virtual impactors
    newks = [ks[1]]
    newVIs = [VIs[1]]
    for (k, VI) in zip(Iterators.drop(ks, 1), Iterators.drop(VIs, 1))
        ka, domaina = newks[end], newVIs[end].domain
        kb, domainb = k, VI.domain
        # The intersection of domaina and domainb is not empty
        if overlap(domaina, domainb)
            domain = (min(domaina[1], domainb[1]), max(domaina[2], domainb[2]))
            σ = midpoint(domain)
            flaga, flagb = domaina[1] ≤ σ ≤ domaina[2], domainb[1] ≤ σ ≤ domainb[2]
            if flaga && flagb
                k = nominaltime(newVIs[end]) < nominaltime(VI) ? ka : kb
            else
                k = flaga ? ka : kb
            end
            newVIs[end] = VirtualImpactor(RTs[k], σ, domain, ctol)
        # The intersection of domaina and domainb is empty
        else
            push!(newks, k)
            push!(newVIs, VI)
        end
    end
    # Sort by time of impact
    sort!(newVIs, by = nominaltime)

    return newVIs
end

"""
    virtualimpactors(IM, lov, RTs, params; kwargs...)

Return the virtual impactors, under the impact monitoring problem
`IM` with line of variations `lov`, found in a vector of returns
`RTs`. For a list of parameters, see the `Propagation` section of
[`Parameters`](@ref).

# Keyword arguments

- `ctol::Real`: convergence tolerance (default: `Inf`).
- `no_pts::Int`: number of points used to find the roots of
    `radialvelocity` in each return (default: `100`).
- `dmax::Real`: maximum allowed value of [`distance`](@ref)
    (default: `0.0`).
"""
function virtualimpactors(IM::AbstractIMProblem{D, T}, lov::LineOfVariations{D, T},
                          RTs::ShowerT1{T}, params::Parameters; ctol::Real = T(Inf),
                          no_pts::Int = 100, dmax::Real = zero(T)) where {D, T <: Real}
    # Find all the virtual impactors in RTs
    σmax = ubound(lov)
    VIs = virtualimpactors(RTs; ctol, σmax, no_pts, dmax)
    # Verify each virtual impactor
    newVIs = Vector{VirtualImpactor{T}}(undef, length(VIs))
    for (i, VI) in enumerate(VIs)
        newVIs[i] = VirtualImpactor(IM, lov, VI, params)
    end
    # Sort by time of impact
    sort!(newVIs, by = nominaltime)

    return newVIs
end