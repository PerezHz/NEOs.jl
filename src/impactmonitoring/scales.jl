"""
    vinf(IM, VI)

Return the hyperbolic excess velocity [km/s], under the impact
monitoring problem `IM`, of a virtual impactor `VI`.
"""
function vinf(IM::IMProblem, VI::VirtualImpactor)
    a = semimajoraxis(VI)
    if isnan(a)
        return a
    elseif ishyperbolic(VI)
        return sqrt( gm(IM) / (-a) ) * (au/daysec)
    else
        return zero(a)
    end
end

"""
    impactenergy(M, V2)

Return the impact energy [Mt] of a virtual impactor with mass `M`
and squared impact velocity `V2`.
"""
impactenergy(M::Real, V2::Real) = 0.5 * M * V2 / 4.184E+9

"""
    impactenergy(IM, VI, params)

Return the impact energy [Mt], under the impact monitoring problem `IM`,
of a virtual impactor `VI`. For a list of parameters, see the `Physical
properties` section of [`Parameters`](@ref).

!!! reference
    - https://doi.org/10.1006/icar.2002.6910
"""
function impactenergy(IM::IMProblem, VI::VirtualImpactor, params::Parameters)
    # Object's mass [kg]
    M = mass(IM, params)
    # Impact velocity^2 [km^2/s^2]
    V2 = vinf(IM, VI)^2 + escapevelocity(IM)^2
    # Impact energy [Mt]
    return impactenergy(M, V2)
end

"""
    palermoscale(E, IP, ΔT)

Return the Palermo Scale of a virtual impactor with energy `E` [Mt],
impact probability `IP` and time until impact `ΔT` [yr].

!!! reference
    - https://doi.org/10.1006/icar.2002.6910
"""
function palermoscale(E::Real, IP::Real, ΔT::Real)
    # Background impact frequency [yr⁻¹]
    fb = 0.03 * E^(-0.8)
    # Palermo scale
    return log10(IP / (fb * ΔT))
end

"""
    palermoscale(IM, VI, params)

Return the Palermo Scale, under the impact monitoring problem `IM`,
of a virtual impactor `VI`. For a list of parameters, see the
`Physical properties` section of [`Parameters`](@ref).

!!! reference
    - https://doi.org/10.1006/icar.2002.6910
"""
function palermoscale(IM::IMProblem, VI::VirtualImpactor, params::Parameters)
    # Impact energy [Mt]
    E = impactenergy(IM, VI, params)
    # Impact probability
    IP = impact_probability(VI)
    # Time until impact [yr]
    ΔT = (nominaltime(VI) - epoch(IM)) / yr
    # Palermo scale
    return palermoscale(E, IP, ΔT)
end

"""
    torinoscale(E, IP)

Return the Torino Scale of a virtual impactor with energy `E` [Mt]
and impact probability `IP`.

!!! reference
    - https://doi.org/10.1016/S0032-0633(00)00006-4
"""
function torinoscale(E::Real, IP::Real)
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

"""
    torinoscale(IM, VI, params)

Return the Torino Scale, under the impact monitoring problem `IM`,
of a virtual impactor `VI`. For a list of parameters, see the
`Physical properties` section of [`Parameters`](@ref).

!!! reference
    - https://doi.org/10.1016/S0032-0633(00)00006-4
"""
function torinoscale(IM::IMProblem, VI::VirtualImpactor, params::Parameters)
    # Impact energy [Mt]
    E = impactenergy(IM, VI, params)
    # Impact probability
    IP = impact_probability(VI)
    # Torino scale
    return torinoscale(E, IP)
end