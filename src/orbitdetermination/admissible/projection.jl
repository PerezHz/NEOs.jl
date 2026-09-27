# Project `[ρ, v_ρ]` into `A`'s outer boundary.
function boundary_projection(A::AdmissibleRegion, ρ::Number, v_ρ::Number)
    # Number of components
    Nc = numberofcomponents(A)
    # Outer boundary limits
    xmin, xleft = A.ρ_domain[1], A.ρ_domain[2]
    xright, xmax = Nc == 1 ? (Inf * x1, Inf * x2) : A.ρ_domain[3], A.ρ_domain[4]
    ymin, ymax = A.v_ρ_domain
    ymid = (ymin + ymax) / 2
    # ρ is left of the tiny object boundary
    if ρ ≤ xmin
        return xmin, clamp(v_ρ, ymin, ymax)
    # ρ is inside one of the components
    elseif xmin < ρ < xleft || (Nc > 1 && xright < ρ < xmax)
        ys = _helrangerates(A.coeffs, A.a_max, ρ)
        length(ys) < 2 && return ρ, ymid
        ymin, ymax = minmax(ys[1], ys[2])
        ymin ≤ v_ρ ≤ ymax && return ρ, v_ρ
        m = v_ρ > ymid ? :max : :min
        xmin, xmax = xmin < ρ < xleft ? (xmin, xleft) : (xright, xmax)
        x = clamp(ρ, xmin, xmax)
        y, dy, d2y = _helrangerate_derivatives(A.coeffs, A.a_max, x, m)
        for _ in 1:25
            dx = (x - ρ + (y - v_ρ) * dy) / (1 + (y - v_ρ) * d2y + dy^2)
            x = clamp(x - dx, xmin, xmax)
            y, dy, d2y = _helrangerate_derivatives(A.coeffs, A.a_max, x, m)
            abs(dx) < eps(T) && break
        end
        return x, y
    # ρ is right of the first component
    elseif Nc == 1
        return xleft, ymid
    # ρ is right of the second component
    elseif ρ ≥ xmax
        return xmax, ymid
    # ρ is between the two components
    else
        x = abs(ρ - xleft) < abs(ρ - xright) ? xleft : xright
        return x, ymid
    end
end

"""
    topo2bary(::AdmissibleRegion, ρ, v_ρ)

Convert topocentric range `ρ` and range-rate `v_ρ` to barycentric
cartesian coordinates. The admissible region fixes the line of sight.
"""
function topo2bary(A::AdmissibleRegion, ρ::Number, v_ρ::Number)
    # Barycentric position
    r = A.observer[1:3] + ρ * A.ρ_unit + A.sun[1:3]
    # Barycentric velocity
    v = A.observer[4:6] + v_ρ * A.ρ_unit + ρ * A.vra * A.ρ_α +
        ρ * A.vdec * A.ρ_δ + A.sun[4:6]
    # Barycentric state vector
    return vcat(r, v)
end

"""
    bary2topo(::AdmissibleRegion, q0)

Convert barycentric cartesian coordinates `q0` to topocentric range
and range-rate. The admissible region fixes the line of sight.
"""
function bary2topo(A::AdmissibleRegion, q0::AbstractVector)
    # Heliocentric state vector
    r = q0 - A.sun
    # Topocentric range
    ρ = euclid3D(r - A.observer)
    # Topocentric range rate
    v_ρ = dot3D(r[4:6], A.ρ_unit) - dot3D(A.observer[4:6], A.ρ_unit) -
        ρ * A.vra * dot3D(A.ρ_α, A.ρ_unit) - ρ * A.vdec * dot3D(A.ρ_δ, A.ρ_unit)
    return ρ, v_ρ
end

"""
    attr2bary(::AdmissibleRegion, attr)

Convert attributable elements `attr` to barycentric cartesian
coordinates. The admissible region fixes the reference epoch.
"""
function attr2bary(A::AdmissibleRegion, attr::AbstractVector)
    # Unfold
    α, δ, v_α, v_δ, ρ, v_ρ = attr
    # Line of sight vectors
    ρ_unit, ρ_α, ρ_δ = topounitpdv(α, δ)
    # Barycentric position
    r = A.observer[1:3] + ρ * ρ_unit + A.sun[1:3]
    # Barycentric velocity
    v = A.observer[4:6] + v_ρ * ρ_unit + ρ * v_α * ρ_α + ρ * v_δ * ρ_δ + A.sun[4:6]
    # Barycentric state vector
    return vcat(r, v)
end