"""
    AbstractOsculatingElements{T <: Real, U <: Number}

Supertype for the osculating orbital elements interface.

Every set of osculating elements `x` has a:
- `gm(x)`: gravitational parameter.
- `epoch(x)`: reference epoch.
- `frame(x)`: reference plane, either `:equatorial` or `:ecliptic`.
- `elements(x)`: set of six osculating orbital elements.
- `covariance(x)`: covariance matrix.
"""
abstract type AbstractOsculatingElements{T <: Real, U <: Number} end

numtypes(::AbstractOsculatingElements{T, U}) where {T, U} = T, U

gm(x::AbstractOsculatingElements) = x.gm

epoch(x::AbstractOsculatingElements) = x.epoch
date(x::AbstractOsculatingElements) = julian2datetime(epoch(x) + JD_J2000 - MJD2000)

frame(x::AbstractOsculatingElements) = x.frame
elements(x::AbstractOsculatingElements) = x.elements

covariance(x::AbstractOsculatingElements) = x.covariance
variances(x::AbstractOsculatingElements) = diag(covariance(x))
sigmas(x::AbstractOsculatingElements) = [y < 0 ? (NaN * y) : sqrt(y) for y in variances(x)]

iscircular(x::AbstractOsculatingElements) = iszero(eccentricity(x))
iselliptic(x::AbstractOsculatingElements) = 0 < eccentricity(x) < 1
isparabolic(x::AbstractOsculatingElements) = isone(eccentricity(x))
ishyperbolic(x::AbstractOsculatingElements) = eccentricity(x) > 1
function conicsection(x::AbstractOsculatingElements)
    if iscircular(x)
        return :circular
    elseif iselliptic(x)
        return :elliptic
    elseif isparabolic(x)
        return :parabolic
    else # ishyperbolic(x)
        return :hyperbolic
    end
end

# Print method for AbstractOsculatingElements
function show(io::IO, x::AbstractOsculatingElements)
    O = elementstype(x)
    T, U = numtypes(x)
    e0 = cte.(elements(x))
    σ0 = sigmas(x)
    se0 = [rpad(@sprintf("%+.12E", e0[i]), 25) for i in eachindex(e0)]
    sσ0 = [rpad(@sprintf("%+.12E", σ0[i]), 25) for i in eachindex(σ0)]
    names = elementsnames(x)
    units = elementsunits(x)
    print(io,
        "$O{$T, $U} $(conicsection(x)) osculating elements\n",
        repeat("-", 67), "\n",
        "Mu: $(gm(x)) au³/day²\n",
        "Epoch: $(epoch(x)) MJD ($(date(x)) TDB)\n",
        "Frame: $(frame(x))\n",
        repeat("-", 67), "\n",
        "Variable    Nominal value            Uncertainty              Units\n",
        rpad(names[1], 12), se0[1], sσ0[1], units[1], "\n",
        rpad(names[2], 12), se0[2], sσ0[2], units[2], "\n",
        rpad(names[3], 12), se0[3], sσ0[3], units[3], "\n",
        rpad(names[4], 12), se0[4], sσ0[4], units[4], "\n",
        rpad(names[5], 12), se0[5], sσ0[5], units[5], "\n",
        rpad(names[6], 12), se0[6], sσ0[6], units[6], "\n",
    )
end

function evaldeltas(y::KeplerianElements{T, TaylorN{T}},
                    dx::Vector{T} = zeros(T, get_numvars())) where {T <: Real}
    O = typeof(y).name.wrapper
    return O{T, T}(gm(y), epoch(y), frame(y), elements(y)(dx), covariance(y)(dx))
end

# Rotate state vector `xas` from equatorial plane to the ecliptic
function equatorial2ecliptic(xas::AbstractVector)
    # Rotation matrix (only positions)
    m_eq2ecl = Rx(deg2rad(ϵ0_deg))
    # Rotational matrix (positions + velocities)
    m_xv_eq2ecl = hcat(vcat(m_eq2ecl, zeros(3,3)), vcat(zeros(3,3), m_eq2ecl))
    # Rotated state vector
    return m_xv_eq2ecl * xas
end

@doc raw"""
    yarkp2adot(A2, a, e; kwargs...)

Return the average semimajor axis drift of an orbit with Yarkovsky coefficient `A2`,
semimajor axis `a` and eccentricity `e`.

# Keyword argument

- `μ`: gravitational parameter of the central body (default: `μ_S`).

!!! reference
    See:
    - https://doi.org/10.1016/j.icarus.2013.02.004

# Extended help

The average semimajor axis drift is given by:
```math
\begin{align*}
    \left\langle\dot{a}\right\rangle
    & = \frac{2A_2(1-e^2)}{n}\left(\frac{1 \ \text{AU}}{p}\right)^2 \\
    & = \frac{2A_2}{(1-e^2)\sqrt{a\mu_\odot}}(1 \ \text{AU})^2,
\end{align*}
```
where ``A_2`` is the Yarkovsky parameter, ``\mu_\odot = GM_\odot`` is the Sun's
gravitational parameter, ``e`` is the eccentricity, ``n = \sqrt{\mu/a^3}`` is the
mean motion, ``p = a(1-e^2)`` is the semilatus rectum, and ``a`` is the semimajor axis.
"""
yarkp2adot(A2, a, e; μ = μ_S) = 2A2 / (sqrt(a) * (1 - e^2) * sqrt(μ))
