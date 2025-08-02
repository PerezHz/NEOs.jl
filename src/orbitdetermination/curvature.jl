"""
    curvature(::AbstractOpticalVector, ::AbstractWeightingScheme)

Return the geodesic curvature, along-track acceleration and their covariance matrix.

!!! reference
    See section 5.1 of:
    - https://doi.org/10.1016/j.icarus.2007.11.033
"""
function curvature(optical::AbstractOpticalVector, w8s::AbstractWeightingScheme)
    # Days of observation [julian days]
    ts = @. dtutc2jdtdb(date(optical))
    # Right ascension [rad]
    αs = ra.(optical)
    # Declination [rad]
    δs = dec.(optical)
    # Weights [rad⁻¹]
    k = 3_600 * 180 / π
    wtα = @. ( k * first(w8s.w8s) )^2
    wtδ = @. ( k * last(w8s.w8s)  )^2

    return curvature(ts, αs, δs, wtα, wtδ)
end

function curvature(ts::AbstractVector{T}, αs::AbstractVector{T}, δs::AbstractVector{T},
                   wtα::AbstractVector{T}, wtδ::AbstractVector{T}) where {T <: Real}
    @assert length(ts) == length(αs) == length(δs) == length(wtα) == length(wtδ) ≥ 3 """
    At least three observations needed for significant curvature computation."""
    # Days of observation [relative to first observation]
    ts = ts .- ts[1]
    # Initial guess for coefficients
    p0 = ones(T, 3)
    # Fit a second order polynomial to ra/dec
    fit_α = curve_fit(polymodel, ts, αs, wtα, p0)
    fit_δ = curve_fit(polymodel, ts, δs, wtδ, p0)
    # Ra/dec covariance matrix
    Γ_α = vcov(fit_α)
    Γ_δ = vcov(fit_δ)
    Γ_αδ = [Γ_α zeros(T, 3, 3); zeros(T, 3, 3) Γ_δ]
    # Mean time of observation
    t_mean = mean(ts)
    # Evaluate ra/dec and its derivatives at mean time
    α = mod2pi(polymodel(t_mean, fit_α.param))
    δ = polymodel(t_mean, fit_δ.param)
    v_α = polymodel(t_mean, diffcoeffs(fit_α.param))
    v_δ = polymodel(t_mean, diffcoeffs(fit_δ.param))
    a_α = polymodel(t_mean, diffcoeffs(diffcoeffs(fit_α.param)))
    a_δ = polymodel(t_mean, diffcoeffs(diffcoeffs(fit_δ.param)))
    # Trigonometric functions
    sin_δ, cos_δ = sincos(δ)
    # Proper motion
    η = sqrt(v_α^2 * cos_δ^2 + v_δ^2)
    # Geodesic curvature
    κ = ((a_δ*v_α - a_α*v_δ) * cos_δ + v_α * (η^2 + v_δ^2) * sin_δ) / η^3
    # Along track acceleration
    v_η = (a_α * v_α * cos_δ^2 + a_δ * v_δ - v_α^2 * v_δ * cos_δ * sin_δ ) / η
    # Angles to curvature jacobian
    J1 = -(-2 * v_α^3 * cos_δ^2 * sin_δ * a_δ + sin_δ * a_δ * v_α * v_δ^2 +
            2 * v_α^2 * cos_δ^2 * sin_δ * a_α * v_δ - sin_δ * a_α * v_δ^3 -
            v_α^5 * cos_δ^3 - 4 * v_α^3 * cos_δ * v_δ^2 + v_α^3 * cos_δ^3 * v_δ^2 -
            2 * v_α * cos_δ * v_δ^4) / η^5
    J2 = -v_α * (sin(2*δ) * (v_α^2 * a_α * cos_δ^2 + 2 * v_δ^2 * a_α - v_α * v_δ * a_δ) +
                 2 * v_α * v_δ^3 * cos(2δ) + 2 * v_α^3 * v_δ * cos_δ^4) / (2 * η^3)
    J3 = (-v_α * cos_δ^3 * (2v_α * a_δ - 3 * v_δ * a_α) + v_δ^2 * (a_δ * cos_δ -
           v_α^2 * cos_δ^2 * sin_δ + 2 * v_δ^2 * sin_δ) ) / η^5
    J4 = - cos_δ * v_δ * (-cos_δ * a_α * v_δ + v_α^3 * sin_δ * cos_δ^2 +
                           2 * v_α * sin_δ * v_δ^2 + cos_δ * a_δ * v_α) / η^3
    J5 = -(cos_δ * (v_α^2 * a_α * cos_δ^2 - 2 * v_δ^2 * a_α + 3 * v_α * v_δ * a_δ) -
                    v_α * v_δ * sin_δ * (v_α^2 * cos_δ^2 - 2 * v_δ^2)) / η^5
    J6 = - v_α * cos_δ^2 * (-a_δ * v_α + v_α^3 * cos_δ * sin_δ + a_α * v_δ) / η^3
    J7 = - v_δ * cos_δ / η^3
    J8 = v_α * cos_δ^2 / η
    J9 = v_α * cos_δ / η^3
    J10 = v_δ / η
    J = [
        0 J1 J3 J5 J7 J9;
        0 J2 J4 J6 J8 J10
    ]
    # Curvature covariance matrix
    Γ_C = J * Γ_αδ * J'

    return [κ, v_η], Γ_C
end