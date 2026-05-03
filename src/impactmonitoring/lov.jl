"""
    LineOfVariationsBuffer{T <: Real} <: AbstractBuffer

Pre-allocated memory for [`lineofvariations`](@ref).

# Fields

- `resTN::Vector{OpticalResidual{T, TaylorN{T}}}`: buffer for `TaylorN{T}` residuals.
- `resT1::Vector{OpticalResidual{T, Taylor1{T}}}`: buffer for `Taylor1{T}` residuals.
- `bufferTN::PropresBuffer{T, TaylorN{T}, T}`: buffer for `TaylorN{T}` propagations.
- `bufferT1::PropresBuffer{T, Taylor1{T}, T}`: buffer for `Taylor1{T}` propagations.
"""
struct LineOfVariationsBuffer{T <: Real} <: AbstractBuffer
    resTN::Vector{OpticalResidual{T, TaylorN{T}}}
    resT1::Vector{OpticalResidual{T, Taylor1{T}}}
    bufferTN::PropresBuffer{T, TaylorN{T}, T}
    bufferT1::PropresBuffer{T, Taylor1{T}, T}
end

get_order(x::LineOfVariationsBuffer) = get_order(x.bufferT1.prop.cache.x[1][0])

"""
    LineOfVariationsBuffer(IM, lovorder, params)

Return a `LineOfVariationsBuffer` object with pre-allocated
memory for [`lineofvariations`](@ref).

# Arguments

- `IM::IMProblem`: impact monitoring problem.
- `lovorder::Int`: order of Taylor expansions wrt LOV index.
- `params::Parameters`: see the `Propagation` section of [`Parameters`](@ref).
"""
function LineOfVariationsBuffer(IM::AbstractIMProblem{D, T}, lovorder::Int,
                                params::Parameters{T}) where {D, T <: Real}
    # Unpack
    @unpack orbit = IM
    # Refence epoch [julian date TDB]
    jd0 = epoch(orbit) + PE.J2000
    # Set jet transport order
    set_od_order(T, 2, dof(IM))
    # Initial condition
    q00 = orbit()
    q0TN = q00 + sigmas(orbit) .* get_variables(T, 2)
    q0T1 = q00 + sigmas(orbit) .* Taylor1(lovorder)
    # Vectors of residuals
    resTN = init_optical_residuals(TaylorN{T}, IM)
    resT1 = init_optical_residuals(Taylor1{T}, IM)
    # Propagation and residuals buffers
    bufferTN = PropresBuffer(IM, q0TN, jd0, params)
    bufferT1 = PropresBuffer(IM, q0T1, jd0, params)

    return LineOfVariationsBuffer{T}(resTN, resT1, bufferTN, bufferT1)
end

"""
    LineOfVariations{D, T} <: AbstractLineOfVariations{T}

A parametrization of the line of variations (LOV).

# Fields

- `dynamics::D`: dynamical model.
- `epoch::T`: reference epoch [days since J2000 TDB].
- `sun::Vector{T}`: Sun's barycentric cartesian state vector at `epoch` [au, au/day].
- `coord::Symbol`: integration coordinates.
- `domain::NTuple{2, T}`: integrated domain.
- `bwd/fwd::DensePropagation2{T, T}`: backward (forward) integration.
"""
@auto_hash_equals struct LineOfVariations{D, T} <: AbstractLineOfVariations{T}
    dynamics::D
    epoch::T
    coord::Symbol
    sun::Vector{T}
    domain::NTuple{2, T}
    bwd::DensePropagation2{T, T}
    fwd::DensePropagation2{T, T}
end

# Print method for LineOfVariations
function show(io::IO, x::LineOfVariations)
    S = titlecase(string(x.coord))
    D, T = numtypes(x)
    t = date(x)
    domain = x.domain
    print(io, S, " LOV{", D.instance, ", ", T, "} at ", t, " over ", domain)
end

# AbstractLineOfVariations interface
numtypes(::LineOfVariations{D, T}) where {D, T} = D, T

dynamicalmodel(x::LineOfVariations) = x.dynamics

epoch(x::LineOfVariations) = x.epoch
nominaltime(x::LineOfVariations) = x.epoch

sigma(::LineOfVariations{D, T}) where {D, T} = zero(T)

get_order(x::LineOfVariations) = get_order(first(x.bwd.x))

function (x::LineOfVariations)(σ::Number)
    mjd0 = epoch(x) + MJD2000
    y = σ >= 0 ? x.fwd(σ) : x.bwd(σ)
    return lovtransform(mjd0, y, x.sun, Val(x.coord), Val(:cartesian))
end

(x::LineOfVariations)(σ::Number, domain::NTuple{2, <:Number}) =
    x(σ + max(domain[2] - σ, σ - domain[1]) * Taylor1(get_order(x)))

# Return the Taylor expansion of the vector field associated to the
# weak direction of Γ; see https://doi.org/10.1137/23M1551961
function weakfield(Γ::AbstractMatrix{Taylor1{T}}, order::Int) where {T <: Real}
    # Allocate memory
    Lx, Ly = size(Γ)
    λs = Vector{T}(undef, order+1)
    vs = Matrix{T}(undef, Lx, order+1)
    As = Array{T}(undef, Lx, Ly, order+1)
    # 0th-order eigenvalues
    As[:, :, 1] = Symmetric(constant_term.(Γ))
    E0 = eigen(As[:, :, 1])
    λs[1], vs[:, 1] = E0.values[end], E0.vectors[:, end]
    # Coefficients matrix
    E = [0 vs[:, 1]'; vs[:, 1] λs[1]*I-As[:, :, 1]]
    # Main loop
    for k in 1:order
        As[:, :, k+1] = TS.differentiate.(k, Γ)
        y, z = zeros(T, Lx), zero(T)
        for l in 0:k-1
            y += binomial(k, l) * As[:, :, k-l+1] * vs[:, l+1]
            if l >= 1
                y -= binomial(k, l) * vs[:, k-l+1] * λs[l+1]
                if l < k -1
                    z += binomial(k-1, l-1) * vs[:, k-l+1]' * vs[:, l+1]
                end
            end
        end
        Ek = E \ [-z/2; y]
        λs[k+1], vs[:, k+1] = Ek[1], Ek[2:end]
    end
    λ = Taylor1(λs, order)
    v = [Taylor1(vs[i, :], order) for i in axes(vs, 1)]
    F = sqrt(λ) * v
    return F
end

lovtransform(::Real, x::AbstractVector, ::AbstractVector, ::Val{:cartesian},
    ::Val{:cartesian}) = collect(x)

lovtransform(::Real, x::AbstractVector, sun::AbstractVector, ::Val{:cartesian},
    ::Val{:equinoctial}) = collect(cartesian2equinoctial(equatorial2ecliptic(x - sun)))

lovtransform(::Real, x::AbstractVector, sun::AbstractVector, ::Val{:equinoctial},
    ::Val{:cartesian}) = collect(ecliptic2equatorial(equinoctial2cartesian(x)) + sun)

lovtransform(::Real, x::AbstractVector, sun::AbstractVector, ::Val{:cartesian},
    ::Val{:attributable}) = collect(cartesian2attributable(equatorial2ecliptic(x - sun)))

lovtransform(::Real, x::AbstractVector, sun::AbstractVector, ::Val{:attributable},
    ::Val{:cartesian}) = collect(ecliptic2equatorial(attributable2cartesian(x)) + sun)

lovtransform(t::Real, x::AbstractVector, sun::AbstractVector, ::Val{:cartesian},
    ::Val{:keplerian}) = collect(cartesian2keplerian(equatorial2ecliptic(x - sun), t))

lovtransform(t::Real, x::AbstractVector, sun::AbstractVector, ::Val{:keplerian},
    ::Val{:cartesian}) = collect(ecliptic2equatorial(keplerian2cartesian(x, t)) + sun)

# Return the Taylor expansion of the vector field associated to the
# weak direction of an orbit
function weakfield(
        IM::AbstractIMProblem{D, T}, coord00::Vector{T}, coord::Symbol,
        buffer::LineOfVariationsBuffer{T}, params::Parameters{T}
    ) where {D, T <: Real}
    # Unpack
    @unpack orbit = IM
    @unpack eph_su = params
    @unpack resTN, resT1, bufferTN, bufferT1 = buffer
    # Refence epoch [days since J2000 / JDTDB / MJDTDB]
    t0 = epoch(orbit)
    jd0 = t0 + PE.J2000
    mjd0 = t0 + MJD2000
    # Order with respect to LOV index
    order = get_order(buffer)
    # Sun's state vector at jd0
    sun = eph_su(t0)
    # Jet transpot initial condition
    scalings = fill(1e-8, 6)
    if get_numvars() == 9
        scalings = vcat(scalings, params.marsden_scalings...)
    end
    car00 = lovtransform(mjd0, coord00, sun, Val(coord), Val(:cartesian))
    carTN = car00 + scalings .* get_variables(T, 2)
    # TaylorN propagation and residuals
    propres!(resTN, IM, carTN, jd0, params; buffer = bufferTN)
    # Covariance matrix in residuals space
    QTN = nms(resTN)
    CTN_res = notout(resTN) * TS.hessian(QTN)
    ΓTN_res = inv(Symmetric(CTN_res))
    # Covariance matrix in coordinate space
    coordTN = lovtransform(mjd0, carTN, sun, Val(:cartesian), Val(coord))
    Γ_coord = Symmetric(project(coordTN, zeros(T, get_numvars()), ΓTN_res))
    # Weak direction eigenpair
    E_coord = eigen(Γ_coord)
    k1, v1 = sqrt(E_coord.values[end]), E_coord.vectors[:, end]
    # Perturbation along the weak direction
    coordT1 = coord00 + k1 * v1 * Taylor1(order)
    carT1 = lovtransform(mjd0, coordT1, sun, Val(coord), Val(:cartesian))
    # Note: the above transformation is only consistent to first order
    for i in eachindex(carT1)
        carT1[i][2:end] .= zero(T)
    end
    # Taylor1 propagation and residuals
    propres!(resT1, IM, carT1, jd0, params; buffer = bufferT1)
    # Covariance matrix in residuals space
    QT1 = nms(resT1)
    dQT1, d2QT1 = zero(QT1), zero(QT1)
    TS.differentiate!(dQT1, QT1)
    TS.differentiate!(d2QT1, dQT1)
    CT1_res = notout(resT1) * d2QT1
    ΓT1_res = inv(CT1_res)
    # Covariance matrix in coordinate space
    coordT1 = lovtransform(mjd0, carT1, sun, Val(:cartesian), Val(coord))
    ΓT1_coord = Symmetric(project(coordT1, ΓT1_res))
    # Vector field associated to the weak direction
    F = weakfield(ΓT1_coord, order)

    return F
end

"""
    lov!

Definition of the line of variations as a differential equation.

!!! reference
    See equation (10.2) in section 10.1:
    - https://doi.org/10.1017/CBO9781139175371
"""
function lov!(dq, q, params, t)
    local F = weakfield(params[1], cte(q), params[2], params[3], params[4])
    dq[1] = F[1]
    dq[2] = F[2]
    dq[3] = F[3]
    dq[4] = F[4]
    dq[5] = F[5]
    dq[6] = F[6]
    return nothing
end

# Methods of `TaylorIntegration._allocate_jetcoeffs!` and `TaylorIntegration.jetcoeffs!`
# generated by @taylorize for lov!
function TaylorIntegration._allocate_jetcoeffs!(
        ::Val{lov!}, t::Taylor1{_T}, q::AbstractArray{Taylor1{_S}, _N},
        dq::AbstractArray{Taylor1{_S}, _N}, params
    ) where {_T <: Real, _S <: Number, _N}
    order = get_order(t)
    local F = weakfield(params[1], cte(q), params[2], params[3], params[4])
    dq[1] = Taylor1(identity(constant_term(F[1])), order)
    dq[2] = Taylor1(identity(constant_term(F[2])), order)
    dq[3] = Taylor1(identity(constant_term(F[3])), order)
    dq[4] = Taylor1(identity(constant_term(F[4])), order)
    dq[5] = Taylor1(identity(constant_term(F[5])), order)
    dq[6] = Taylor1(identity(constant_term(F[6])), order)
    return TaylorIntegration.RetAlloc{Taylor1{_S}}(
        Taylor1{_S}[], [Array{Taylor1{_S}, 1}(undef, 0)],
        [Array{Taylor1{_S}, 2}(undef, 0, 0)],
        [Array{Taylor1{_S}, 3}(undef, 0, 0, 0)],
        [Array{Taylor1{_S}, 4}(undef, 0, 0, 0, 0)]
    )
end

function TaylorIntegration.jetcoeffs!(
        ::Val{lov!}, t::Taylor1{_T}, q::AbstractArray{Taylor1{_S}, _N},
        dq::AbstractArray{Taylor1{_S}, _N}, params,
        __ralloc::TaylorIntegration.RetAlloc{Taylor1{_S}}
    ) where {_T <: Real, _S <: Number, _N}
    order = get_order(t)
    local F = weakfield(params[1], cte(q), params[2], params[3], params[4])
    for ord = 0:order - 1
        ordnext = ord + 1
        TaylorSeries.identity!(dq[1], F[1], ord)
        TaylorSeries.identity!(dq[2], F[2], ord)
        TaylorSeries.identity!(dq[3], F[3], ord)
        TaylorSeries.identity!(dq[4], F[4], ord)
        TaylorSeries.identity!(dq[5], F[5], ord)
        TaylorSeries.identity!(dq[6], F[6], ord)
        for __idx = eachindex(q)
            TaylorIntegration.solcoeff!(q[__idx], dq[__idx], ordnext)
        end
    end
    return nothing
end

"""
    lineofvariations(IM, params; kwargs...)

Return a parametrization of the line of variations associated to an
impact monitoring problem `IM`. For a list of parameters, see the
`Propagation` section of [`Parameters`](@ref).

# Keyword arguments

- `coord::Symbol`: coordinates for the LOV; accepted values are: `:cartesian`
    (default), `:equinoctial`, `:keplerian` and `:attributable`.
- `σmax::Real`: maximum (absolute) value of the LOV index (default: `3.0`).
- `lovorder::Int`: order of Taylor expansions wrt LOV index (default: `12`).
- `lovtol::Real`: absolute tolerance used to integrate the LOV (default: `1E-20`).
- `lovsteps::Int`: maximum number of steps for the integration (default: `100`).
- `lovparse::Bool`: whether to use the specialized method of `jetcoeffs` or not
    (default: `true`).
"""
function lineofvariations(IM::AbstractIMProblem{D, T}, params::Parameters{T};
                          coord::Symbol = :cartesian, σmax::Real = 3.0,
                          lovorder::Int = 12, lovtol::Real = 1E-20,
                          lovsteps::Int = 100, lovparse::Bool = true) where {D, T <: Real}
    # Unpack
    @unpack orbit = IM
    @unpack eph_su = params
    # Refence epoch [days since J2000 / MJDTDB]
    t0 = epoch(orbit)
    mjd0 = t0 + MJD2000
    # Sun's state vector at t0 [au, au/day]
    sun = eph_su(t0)
    # Initial condition
    q00 = orbit()
    # Line of variations buffer
    buffer = LineOfVariationsBuffer(IM, lovorder, params)
    # Taylor expansion of the line of variations
    coord00 = lovtransform(mjd0, q00, sun, Val(:cartesian), Val(coord))
    lovparams = (IM, coord, buffer, params)
    _bwd_ = taylorinteg(lov!, coord00, zero(T), -σmax, lovorder, lovtol, lovparams;
                        maxsteps = lovsteps, parse_eqs = lovparse, dense = true)
    _fwd_ = taylorinteg(lov!, coord00, zero(T), σmax, lovorder, lovtol, lovparams;
                        maxsteps = lovsteps, parse_eqs = lovparse, dense = true)
    bwd = TaylorInterpolant{T, T, 2}(zero(T), _bwd_.t, _bwd_.p)
    fwd = TaylorInterpolant{T, T, 2}(zero(T), _fwd_.t, _fwd_.p)
    # The positive sigma direction is given by the semimajor axis
    q0T1 = lovtransform(mjd0, fwd.x[1, :], sun, Val(coord), Val(:cartesian)) - sun
    a = semimajoraxis(q0T1..., μ_S, 0.0)
    if differentiate(1, a) < 0
        bwd, fwd = flipsign(fwd), flipsign(bwd)
    end
    domain = (last(bwd.t), last(fwd.t))

    return LineOfVariations{D, T}(dynamicalmodel(IM), t0, coord, sun, domain, bwd, fwd)
end