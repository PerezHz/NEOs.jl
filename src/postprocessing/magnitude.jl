# Observed magnitude model
# h = H + 5 * log10(d_BS * d_BO) - 2.5 * log10(qs)
hmodel(t, p) = @. p[1] + t

@doc raw"""
    absolutemagnitude(::AbstractOrbit{T, T}, ::Parameters{T}) where {T <: Real}

Return the absolute magnitude of an orbit, as well as its standard error.
This function uses the linear H and G model for asteroids.

!!! reference
    See https://minorplanetcenter.net/iau/ECS/MPCArchive/1985/MPC_19851227.pdf.
"""
function absolutemagnitude(orbit::AbstractOrbit{T, T}, params::Parameters{T}) where {T <: Real}
    # Extract optical astrometry
    radec = astrometry(orbit.tracklets)
    # Observation times
    ts = @. dtutc2days(date(radec))
    # Observed magnitudes and bands
    hs, bands = @. getfield(radec, :mag), getfield(radec, :band)
    # Convert magnitudes to V band
    hs .+= getindex.(Ref(V_BAND_CORRECTION), bands)
    # Asteroid, Earth and Sun positions
    xa = [orbit(t)[1:3] for t in ts]
    xe = [params.eph_ea(t)[1:3] for t in ts]
    xs = [params.eph_su(t)[1:3] for t in ts]
    # Distances
    d_BS = @. norm(xa - xs)     # Asteroid-Sun
    d_OS = @. norm(xe - xs)     # Earth-Sun
    d_BO = @. norm(xa - xe)     # Asteroid-Earth
    # Phase angles
    αs = @. acos((d_BO^2 + d_BS^2 - d_OS^2) / (2 * d_BO * d_BS))
    # Phase integrals
    ϕ1s = @. exp(-PHASE_INTEGRAL_A1 * tan(αs/2)^PHASE_INTEGRAL_B1)
    ϕ2s = @. exp(-PHASE_INTEGRAL_A2 * tan(αs/2)^PHASE_INTEGRAL_B2)
    qs = (1 - SLOPE_PARAMETER) * ϕ1s + SLOPE_PARAMETER * ϕ2s
    # Absolute magnitude
    tdata = @. 5 * log10(d_BS * d_BO) - 2.5 * log10(qs)
    mask = @. !isnan(tdata) && !isnan(hs)
    H0 = [mean(view(hs, mask))]
    fit = curve_fit(hmodel, view(tdata, mask), view(hs, mask), H0)
    H = fit.param[1]
    dH = stderror(fit)[1]

    return H, dH
end