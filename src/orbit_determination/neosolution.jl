@doc raw"""
    NEOSolution{T <: Real, U <: Number}

The outcome of the orbit determination process for a NEO.

# Fields

- `tracklets::Vector{Tracklet{T}}`: vector of tracklets.
- `bwd/fwd::TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}}`: backward (forward) integration.
- `t_bwd/t_fwd::Vector{U}`: time of Earth close approach.
- `x_bwd/x_fwd::Vector{U}`: state vector at Earth close approach.
- `g_bwd/g_fwd::Vector{U}`: geocentric distance at close approach.
- `res::Vector{OpticalResidual{T, U}}`: vector of optical residuals.
- `fit::LeastSquaresFit{T}`: least squares fit.
- `scalings::Vector{T}`: jet transport scaling factors.
"""
@auto_hash_equals struct NEOSolution{T <: Real, U <: Number}
    tracklets::Vector{Tracklet{T}}
    bwd::TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}}
    t_bwd::Vector{U}
    x_bwd::Matrix{U}
    g_bwd::Vector{U}
    fwd::TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}}
    t_fwd::Vector{U}
    x_fwd::Matrix{U}
    g_fwd::Vector{U}
    res::Vector{OpticalResidual{T, U}}
    fit::LeastSquaresFit{T}
    scalings::Vector{T}
    # Inner constructor
    function NEOSolution{T, U}(
        tracklets::Vector{Tracklet{T}},
        bwd::TaylorInterpolant{T, U, 2, VT, X}, t_bwd::Vector{U}, x_bwd::Matrix{U}, g_bwd::Vector{U},
        fwd::TaylorInterpolant{T, U, 2, VT, X}, t_fwd::Vector{U}, x_fwd::Matrix{U}, g_fwd::Vector{U},
        res::Vector{OpticalResidual{T, U}}, fit::LeastSquaresFit{T}, scalings::Vector{T}
    ) where {T <: Real, U <: Number, VT <: AbstractVector{T}, X <: AbstractMatrix{Taylor1{U}}}
        @assert bwd.t0 == fwd.t0 "Backward and forward propagation initial times must match"
        _bwd_ = TaylorInterpolant(bwd.t0, bwd.t, collect(bwd.x))
        _fwd_ = TaylorInterpolant(fwd.t0, fwd.t, collect(fwd.x))
        new{T, U}(
            tracklets,
            _bwd_, t_bwd, x_bwd, g_bwd,
            _fwd_, t_fwd, x_fwd, g_fwd,
            res, fit, scalings
        )
    end
end
# Outer constructors
function NEOSolution(
    tracklets::Vector{Tracklet{T}},
    bwd::TaylorInterpolant{T, U, 2, VT, X}, t_bwd::Vector{U}, x_bwd::Matrix{U}, g_bwd::Vector{U},
    fwd::TaylorInterpolant{T, U, 2, VT, X}, t_fwd::Vector{U}, x_fwd::Matrix{U}, g_fwd::Vector{U},
    res::Vector{OpticalResidual{T, U}}, fit::LeastSquaresFit{T}, scalings::Vector{T}
) where {T <: Real, U <: Number, VT <: AbstractVector{T}, X <: AbstractMatrix{Taylor1{U}}}
    NEOSolution{T, U}(
        tracklets,
        bwd, t_bwd, x_bwd, g_bwd,
        fwd, t_fwd, x_fwd, g_fwd,
        res, fit, scalings
    )
end

function NEOSolution(
    tracklets::Vector{Tracklet{T}},
    bwd::TaylorInterpolant{T, U, 2, VT, X}, fwd::TaylorInterpolant{T, U, 2, VT, X},
    res::Vector{OpticalResidual{T, U}}, fit::LeastSquaresFit{T}, scalings::Vector{T}
) where {T <: Real, U <: Number, VT <: AbstractVector{T}, X <: AbstractMatrix{Taylor1{U}}}
    # Backward roots
    t_bwd = Vector{U}(undef, 0)
    x_bwd = Matrix{U}(undef, 0, 0)
    g_bwd = Vector{U}(undef, 0)
    # Forward roots
    t_fwd = Vector{U}(undef, 0)
    x_fwd = Matrix{U}(undef, 0, 0)
    g_fwd = Vector{U}(undef, 0)

    NEOSolution{T, U}(
        tracklets,
        bwd, t_bwd, x_bwd, g_bwd,
        fwd, t_fwd, x_fwd, g_fwd,
        res, fit, scalings
    )
end

# Print method for NEOSolution
# Examples:
# NEO solution with 123 residuals
function show(io::IO, x::NEOSolution{T, U}) where {T <: Real, U <: Number}
    print(io, "NEO solution with ", length(x.res), " residuals")
end

# Evaluation in time method
function (sol::NEOSolution{T, U})(t::V = sol.bwd.t0) where {T <: Real, U <: Number, V <: Number}
    if t <= sol.bwd.t0
        return sol.bwd(t)
    else
        return sol.fwd(t)
    end
end

# Evaluation in fit δs method
function evalfit(sol::NEOSolution{T, TaylorN{T}}) where {T <: Real}
    # Fit δs
    δs = sol.fit.x
    # Evaluate backward integration
    new_bwd_x = map(x -> Taylor1(x.coeffs(δs)), sol.bwd.x);
    new_bwd = TaylorInterpolant(sol.bwd.t0, sol.bwd.t, new_bwd_x)
    # Evaluate backward roots
    new_t_bwd = sol.t_bwd(δs)
    new_x_bwd = sol.x_bwd(δs)
    new_g_bwd = sol.g_bwd(δs)
    # Evaluate forward integration
    new_fwd_x = map(x -> Taylor1(x.coeffs(δs)), sol.fwd.x);
    new_fwd = TaylorInterpolant(sol.fwd.t0, sol.fwd.t, new_fwd_x)
    # Evaluate forward roots
    new_t_fwd = sol.t_fwd(δs)
    new_x_fwd = sol.x_fwd(δs)
    new_g_fwd = sol.g_fwd(δs)
    # Evaluate residuals
    new_res = sol.res(δs)

    NEOSolution{T, T}(
        sol.tracklets,
        new_bwd, new_t_bwd, new_x_bwd, new_g_bwd,
        new_fwd, new_t_fwd, new_x_fwd, new_g_fwd,
        new_res, sol.fit, sol.scalings
    )
end

# Definition of zero NEOSolution
function zero(::Type{NEOSolution{T, U}}) where {T <: Real, U <: Number}
    tracklets = Vector{Tracklet{T}}(undef, 0)
    bwd = zero(TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}})
    t_bwd = Vector{U}(undef, 0)
    x_bwd = Matrix{U}(undef, 0, 0)
    g_bwd = Vector{U}(undef, 0)
    fwd = zero(TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}})
    t_fwd = Vector{U}(undef, 0)
    x_fwd = Matrix{U}(undef, 0, 0)
    g_fwd = Vector{U}(undef, 0)
    res = Vector{OpticalResidual{T, U}}(undef, 0)
    fit = LeastSquaresFit{T}(false, Vector{T}(undef, 0), Matrix{T}(undef, 0, 0), :newton)
    scalings = Vector{T}(undef, 0)

    NEOSolution{T, U}(
        tracklets,
        bwd, t_bwd, x_bwd, g_bwd,
        fwd, t_fwd, x_fwd, g_fwd,
        res, fit, scalings
    )
end

iszero(x::NEOSolution{T, U}) where {T <: Real, U <: Number} = x == zero(NEOSolution{T, U})

# Normalized Root Mean Square Error
function nrms(sol::NEOSolution{T, T}) where {T <: Real}
    if iszero(sol)
        return T(Inf)
    else
        return nrms(sol.res)
    end
end
# Normalized Mean Square Error
function nms(sol::NEOSolution{T, T}) where {T <: Real}
    if iszero(sol)
        return T(Inf)
    else
        return nms(sol.res)
    end
end

@doc raw"""
    sigmas(sol::NEOSolution{T, T}) where {T <: Real}

Return `sol`'s initial condition uncertainties in barycentric cartesian coordinates.
"""
sigmas(sol::NEOSolution{T, T}) where {T <: Real} = sqrt.(diag(sol.fit.Γ)) .* sol.scalings

@doc raw"""
    snr(sol::NEOSolution{T, T}) where {T <: Real}

Return `sol`'s initial condition signal-to-noise ratios in barycentric cartesian coordinates.
"""
snr(sol::NEOSolution{T, T}) where {T <: Real} = abs.(sol()) ./ sigmas(sol)

@doc raw"""
    jplcompare(des::String, sol::NEOSolution{T, U}) where {T <: Real, U <: Number}

Return `abs.(sol() - R) ./ sigmas(sol)`, where `R` is JPL's state vector of object
`des` at `sol`'s  initial epoch.
"""
function jplcompare(des::String, sol::NEOSolution{T, U}) where {T <: Real, U <: Number}
    # Load Solar System ephemerides
    loadjpleph()
    # NEOs barycentric state vector
    q1 = cte.(sol())
    # Time of first (last) observation
    t0, tf = date(sol.tracklets[1].radec[1]), date(sol.tracklets[end].radec[end])
    # Download JPL ephemerides
    bsp = smb_spk("DES = $(des);", t0, tf)
    # Object ID
    id = parse(Int, bsp[1:end-4])
    # Load JPL ephemerides
    furnsh(bsp)
    # Delete ephemerides file
    rm(bsp)
    # JPL barycentric state vector
    q2 = kmsec2auday(spkgeo(id, julian2etsecs(sol.bwd.t0 + PE.J2000), "J2000", 0)[1])
    # Absolute difference in sigma units
    return abs.(q1 - q2) ./ sigmas(sol)
end

@doc raw"""
    uncertaintyparameter(radec::Vector{RadecMPC{T}}, sol::NEOSolution{T, T}, params::NEOParameters{T};
                         dynamics::D = newtonian!) where {T <: Real, D}

Return the Minor Planet Center Uncertainty Parameter.

# Arguments

- `radec::Vector{RadecMPC{T}}`: vector of observations.
- `sol::NEOSolution{T, T}:` least squares orbit.
- `params::NEOParameters{T}`: see [`NEOParameters`](@ref).
- `dynamics::D`: dynamical model.

!!! reference
    https://www.minorplanetcenter.net/iau/info/UValue.html
"""
function uncertaintyparameter(radec::Vector{RadecMPC{T}}, sol::NEOSolution{T, T}, params::NEOParameters{T};
                              dynamics::D = newtonian!) where {T <: Real, D}
    # Reduce tracklets by polynomial regression
    tracklets = reduce_tracklets(radec)
    # Epoch [julian days]
    jd0 = sol.bwd.t0 + PE.J2000
    # Barycentric initial conditions
    q0 = sol(sol.bwd.t0)
    # Scaling factors
    scalings = abs.(q0) ./ 10^6
    # Jet transport variables
    dq = scaled_variables("δx", scalings; order = params.varorder)
    # Origin
    x0 = zeros(T, 6)
    # Initial conditions
    q = q0 + dq
    # Propagation and residuals
    bwd, fwd, res = propres(radec, jd0, q, params; dynamics)
    # Orbit fit
    fit = tryls(res, x0, params.niter)
    # Update solution
    _sol_ = NEOSolution(tracklets, bwd, fwd, res, fit, scalings)
    # Osculating keplerian elements
    osc = pv2kep(_sol_() - params.eph_su(sol.bwd.t0); jd = jd0, frame = :ecliptic)
    # Eccentricity
    e = osc.e
    # Gauss gravitational constant [deg]
    k_0 = 180 * k_gauss / π
    # Time of perihelion passage [julian days]
    Tp = osc.tp
    # Orbital period [days]
    P = 2π * sqrt(osc.a^3 / μ_S)
    # Covariance matrix
    Γ = project([Tp, P], fit)
    # Uncertainties
    dTp, dP = sqrt.(diag(Γ))
    # Convert orbital period to years
    P = P/yr
    # In-orbit longitude runoff [arcsec / decade]
    runoff = (dTp * e(fit.x) + 10 * dP / P(fit.x)) * (k_0 / P(fit.x)) * 3_600 * 3
    # Uncertainty parameter
    C = log(648_000) / 9
    U = floor(Int, log(runoff)/C) + 1

    return U
end