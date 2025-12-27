# Return the polynomial coefficients for an [`AdmissibleRegion`](@ref).
# See equation (8.8) of https://doi.org/10.1017/CBO9781139175371
function arcoeffs(α::T, δ::T, v_α::T, v_δ::T, ρ::AbstractVector{T}, ρ_α::AbstractVector{T},
                  ρ_δ::AbstractVector{T}, q::AbstractVector{T}) where {T <: Number}
    coeffs = Vector{T}(undef, 6)
    coeffs[1] = dot3D(q[1:3], q[1:3])
    coeffs[2] = 2 * dot3D(q[4:6], ρ)
    coeffs[3] = v_α^2 * cos(δ)^2 + v_δ^2  # Proper motion squared
    coeffs[4] = 2 * v_α * dot3D(q[4:6], ρ_α) + 2 * v_δ * dot3D(q[4:6], ρ_δ)
    coeffs[5] = dot3D(q[4:6], q[4:6])
    coeffs[6] = 2 * dot3D(q[1:3], ρ)
    return coeffs
end

# Return the topocentric line-of-sight unit vector and its
# partial derivatives with respect to `α` and `δ`.
# See between equations (8.5) and (8.6) of https://doi.org/10.1017/CBO9781139175371
function topounitpdv(α::T, δ::T) where {T <: Number}
    sin_α, cos_α = sincos(α)
    sin_δ, cos_δ = sincos(δ)
    sin_α_sin_δ = sin_α * sin_δ
    sin_α_cos_δ = sin_α * cos_δ
    cos_α_sin_δ = cos_α * sin_δ
    cos_α_cos_δ = cos_α * cos_δ
    ρ = [cos_α_cos_δ, sin_α_cos_δ, sin_δ]
    ρ_α = [-sin_α_cos_δ, cos_α_cos_δ, zero(α)]
    ρ_δ = [-cos_α_sin_δ, -sin_α_sin_δ, cos_δ]
    return ρ, ρ_α, ρ_δ
end

"""
    AdmissibleRegion(::OpticalTracklet, ::Parameters)

Return the admissible region associated to an optical tracklet. For a list of
parameters, see the `Minimization over the MOV` section of [`Parameters`](@ref).
"""
AdmissibleRegion(x::OpticalTracklet, params::Parameters) = AdmissibleRegion(
    date(x), ra(x), dec(x), vra(x), vdec(x), mag(x), observatory(x), params)

function AdmissibleRegion(date::DateTime, α::T, δ::T, v_α::T, v_δ::T,
                          mag::T, observatory::ObservatoryMPC{T},
                          params::Parameters{T}) where {T <: Real}
    # Unpack parameters
    @unpack eph_ea, eph_su, H_max, a_max = params
    # Topocentric unit vector and partials
    ρ, ρ_α, ρ_δ = topounitpdv(α, δ)
    # Time of observation [days (et seconds) since J2000]
    t_days, t_et = dtutc2days(date), dtutc2et(date)
    # Heliocentric position of the observer
    q = eph_ea(t_days) + kmsec2auday(obsposvelECI(observatory, t_et)) - eph_su(t_days)
    # Admissible region coefficients
    coeffs = arcoeffs(α, δ, v_α, v_δ, ρ, ρ_α, ρ_δ, q)
    # Maximum range (heliocentric energy constraint)
    ρ_max = _helmaxrange(coeffs, a_max)
    iszero(ρ_max) && return zero(AdmissibleRegion{T})
    # Minimum range
    if isnan(mag)
        # Earth's sphere of influence radius / Earth's physical radius
        ρ_min = R_SI < ρ_max ? R_SI : R_EA
    else
        # Tiny object boundary
        ρ_min = 10^((mag - H_max)/5)
    end
    ρ_min > ρ_max && return zero(AdmissibleRegion{T})
    # Range domain
    ρ_domain = [ρ_min, ρ_max]
    # Range-rate domain
    v_ρ_min, v_ρ_max = _helrangerates(coeffs, a_max, ρ_min)[1:2]
    v_ρ_domain = [v_ρ_min, v_ρ_max]
    # Range rate symmetry level
    v_ρ_mid = _helrangerates(coeffs, a_max, ρ_max)[1]
    # Boundary points
    Fs = Matrix{T}(undef, 3, 2)
    Fs[1, :] .= [ρ_min, v_ρ_min]
    Fs[2, :] .= [ρ_min, v_ρ_max]
    Fs[3, :] .= [ρ_max, v_ρ_mid]

    return AdmissibleRegion{T}(date, α, δ, v_α, v_δ, H_max, a_max,
        ρ, ρ_α, ρ_δ, q, coeffs, ρ_domain, v_ρ_domain, Fs, observatory)
end