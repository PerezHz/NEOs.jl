@doc raw"""
    yarkp2adot(A2, a, e; kwargs...)

Return the average semimajor axis drift of an orbit with Yarkovsky coefficient `A2`,
semimajor axis `a` and eccentricity `e`.

See also [`adot2yarkp`](@ref).

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
yarkp2adot(A2::Number, a::Number, e::Number; μ::Real = μ_S) =
    2A2 / (sqrt(a * μ) * (1 - e^2))

@doc raw"""
    adot2yarkp(adot, a, e; kwargs...)

Return the Yarkovsky coefficient of an orbit with average semimajor axis drift `adot`,
semimajor axis `a` and eccentricity `e`.

See also [`yarkp2adot`](@ref).

# Keyword argument

- `μ`: gravitational parameter of the central body (default: `μ_S`).

!!! reference
    See:
    - https://doi.org/10.1016/j.icarus.2013.02.004

# Extended help

The Yarkovsky coefficient is given by:
```math
\begin{align*}
    A_2
    & = \left\langle\dot{a}\right\rangle
    \frac{n}{2(1-e^2)}\left(\frac{p}{1 \ \text{AU}}\right)^2
    & = \left\langle\dot{a}\right\rangle
    \frac{(1-e^2)\sqrt{a\mu_\odot}}{2}\left(\frac{1}{1 \ \text{AU}}\right)^2
\end{align*}
```
where ``\left\langle\dot{a}\right\rangle`` is the average semimajor axis drift,
``\mu_\odot = GM_\odot`` is the Sun's gravitational parameter, ``e`` is the eccentricity,
``n = \sqrt{\mu/a^3}`` is the mean motion, ``p = a(1-e^2)`` is the semilatus rectum,
and ``a`` is the semimajor axis.
"""
adot2yarkp(adot::Number, a::Number, e::Number; μ::Real = μ_S) =
    adot * sqrt(a * μ) * (1 - e^2) / 2

@doc raw"""
    expectedadot(a, e, D, ρ, ϕ, A)

Return the expected average semimajor axis drift [au/Myr] of an asteroid with
semimajor axis `a` [au], eccentricity `e`, diameter `D` [km], density `ρ` [g/cm³],
obliquity `ϕ` [deg] and Bond albedo `A`, taking Bennu as a reference.

See also [`yarkp2adot`](@ref) and [`adot2yarkp`](@ref).

!!! reference
    See:
    - https://doi.org/10.1051/0004-6361/201833153

# Extended help

The expected average semimajor axis drift is given by:
```math
\begin{align*}
    \left\langle\dot{a}\right\rangle_{\text{exp}} =
    \left\langle\dot{a}\right\rangle_{B}\times
    \frac{\sqrt{a_B}(1-e^2)}{\sqrt{a}(1-e^2)}
    \frac{D_B}{D}\frac{\rho_B}{\rho}\frac{\cos\phi_B}{\cos\phi}\frac{1-A}{1-A_B}
\end{align*}
```
where ``\left\langle\dot{a}\right\rangle`` is the average semimajor axis drift,
``a`` is the semimajor axis, ``e`` is the eccentricity, ``D`` is the diameter,
``ρ`` is the density, ``ϕ`` is the obliquity and ``A`` is the Bond albedo; the
``B`` subscript denotes Bennu's phisical parameters.
"""
expectedadot(a::Number, e::Number, D::Number, ρ::Number, ϕ::Number, A::Number) =
    BENNU_ADOT * sqrt(BENNU_SEMIMAJOR_AXIS / a) * ((1 - BENNU_ECCENTRICITY^2) / (1 - e^2)) *
    (BENNU_DIAMETER / D) * (BENNU_DENSITY / ρ) * (cosd(ϕ) / cosd(BENNU_OBLIQUITY)) *
    ((1 - A) / (1 - BENNU_BOND_ALBEDO))