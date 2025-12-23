# Project `[ρ, v_ρ]` into `A`'s outer boundary.
function boundary_projection(A::AdmissibleRegion{T}, ρ::T, v_ρ::T) where {T <: Real}
    # Project range
    ρ = clamp(ρ,  A.ρ_domain[1], A.ρ_domain[2])
    # Project range-rate
    y_domain = rangerates(A, ρ, :outer)
    if iszero(length(y_domain))
        v_ρ = sum(A.v_ρ_domain) / 2
    elseif isone(length(y_domain))
        v_ρ = y_domain[1]
    else
        v_ρ = clamp(v_ρ, y_domain[1], y_domain[2])
    end

    return ρ, v_ρ
end

"""
    topo2bary(::AdmissibleRegion, ρ, v_ρ)

Convert topocentric range `ρ` and range-rate `v_ρ` to barycentric cartesian coordinates.
The admissible region fixes the line of sight.
"""
function topo2bary(A::AdmissibleRegion, ρ::Number, v_ρ::Number)
    # Barycentric position
    r = A.q[1:3] + ρ * A.ρ_unit + sseph(su, dtutc2days(A.date))[1:3]
    # Barycentric velocity
    v = A.q[4:6] + v_ρ * A.ρ_unit + ρ * A.vra * A.ρ_α + ρ * A.vdec * A.ρ_δ
        + sseph(su, dtutc2days(A.date))[4:6]
    # Barycentric state vector
    return vcat(r, v)
end

"""
    bary2topo(::AdmissibleRegion, q0)

Convert barycentric cartesian coordinates `q0` to topocentric range and range-rate.
The admissible region fixes the line of sight.
"""
function bary2topo(A::AdmissibleRegion, q0::Vector{<:Number})
    # Heliocentric state vector
    r = q0 - sseph(su, dtutc2days(A.date))
    # Topocentric range
    ρ = euclid3D(r - A.q)
    # Topocentric range rate
    v_ρ = dot3D(r[4:6], A.ρ_unit) - dot3D(A.q[4:6], A.ρ_unit) - ρ * A.vra * dot3D(A.ρ_α, A.ρ_unit)
          - ρ * A.vdec * dot3D(A.ρ_δ, A.ρ_unit)

    return ρ, v_ρ
end

"""
    attr2bary(::AdmissibleRegion, a, ::Parameters)

Convert attributable elements `a` to barycentric cartesian coordinates.
The admissible region fixes the reference epoch and the parameters provide
Sun and Earth's ephemerides.
"""
function attr2bary(A::AdmissibleRegion{T}, a::Vector{U},
                   params::Parameters{T}) where {T <: Real, U <: Number}
    # Unfold
    α, δ, v_α, v_δ, ρ, v_ρ = a
    # Admissible region reference epoch
    # Note: we concluded both t and et should not include the relativistic
    # correction -ρ/c for consistency (05/12/2025)
    t = dtutc2days(A.date) # - ρ / c_au_per_day
    # TO DO: `et::TaylorN` is too slow for `mmov` due to
    # SatelliteToolboxTransformations overloads in src/observations/topocentric.jl
    et = dtutc2et(A.date) # - cte(cte(ρ)) / c_au_per_sec
    # Line of sight vectors
    ρ_unit, ρ_α, ρ_δ = topounitpdv(α, δ)
    # Heliocentric position of the observer
    q = params.eph_ea(t) + kmsec2auday(obsposvelECI(A.observatory, et)) - params.eph_su(t)
    # Barycentric position
    r = q[1:3] + ρ * ρ_unit + params.eph_su(t)[1:3]
    # Barycentric velocity
    v = q[4:6] + v_ρ * ρ_unit + ρ * v_α * ρ_α + ρ * v_δ * ρ_δ
        + params.eph_su(t)[4:6]
    # Barycentric state vector
    return vcat(r, v)
end