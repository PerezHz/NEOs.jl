"""
    LineOfVariations{D, T} <: AbstractLineOfVariations{T}

A parametrization of the line of variations (LOV).

# Fields

- `dynamics::D`: dynamical model.
- `epoch::T`: reference epoch [days since J2000 TDB].
- `domain::NTuple{2, T}`: integrated domain.
- `bwd/fwd::DensePropagation2{T, T}`: backward (forward) integration.
"""
@auto_hash_equals struct LineOfVariations{D, T} <: AbstractLineOfVariations{T}
    dynamics::D
    epoch::T
    domain::NTuple{2, T}
    bwd::DensePropagation2{T, T}
    fwd::DensePropagation2{T, T}
end

# Print method for LineOfVariations
function show(io::IO, x::LineOfVariations)
    D, T = numtypes(x)
    t = date(x)
    domain = x.domain
    print(io, "LineOfVariations{", D.instance, ", ", T, "} at ", t, " over ", domain)
end

# AbstractLineOfVariations interface
numtypes(::LineOfVariations{D, T}) where {D, T} = D, T

dynamicalmodel(x::LineOfVariations) = x.dynamics

epoch(x::LineOfVariations) = x.epoch
nominaltime(x::LineOfVariations) = x.epoch

sigma(::LineOfVariations{D, T}) where {D, T} = zero(T)

get_order(x::LineOfVariations) = get_order(first(x.bwd.x))

(x::LineOfVariations)(σ::Number) = σ >= 0 ? x.fwd(σ) : x.bwd(σ)

(x::LineOfVariations)(σ::Number, domain::NTuple{2, <:Number}) =
    x(σ + max(domain[2] - σ, σ - domain[1]) * Taylor1(get_order(x)))

# Return the Taylor expansion of the covariance matrix of an initial condition,
# with respect to the normalized direction of the eigenvector associated to the
# largest eigenvalue
function lovcovariance(
        resTN::Vector{OpticalResidual{T, TaylorN{T}}},
        resT1::Vector{OpticalResidual{T, Taylor1{T}}},
        IM::AbstractIMProblem{D, T}, q00::Vector{T},
        jd0::T, params::Parameters{T};
        bufferTN::Union{Nothing, PropresBuffer{T, TaylorN{T}, T}} = nothing,
        bufferT1::Union{Nothing, PropresBuffer{T, Taylor1{T}, T}} = nothing,
        order::Int = 12
    ) where {D, T <: Real}
    # Set jet transport variables
    Npar = dof(IM)
    set_od_order(T, 2, Npar)
    # Jet transpot initial condition
    scalings = fill(1e-8, 6)
    if Npar == 9
        scalings = vcat(scalings, params.marsden_scalings...)
    end
    dq = scalings .* get_variables(T, 2)
    q0TN = q00 + dq
    # Propagation and residuals
    propres!(resTN, IM, q0TN, jd0, params; buffer = bufferTN)
    # Covariance matrix in residuals space
    QTN = nms(resTN)
    C = notout(resTN) * TS.hessian(QTN)
    Γ_ξ = inv(Symmetric(C))
    # Covariance matrix in cartesian space
    J = Matrix(TS.jacobian(dq))
    Γ_c = Symmetric(J * Γ_ξ * J')
    # Greatest eigenpair
    E = eigen(Γ_c)
    k1, v1 = sqrt(E.values[end]), E.vectors[:, end]
    # Line of variations initial condition
    q0T1 = q00 + k1 * v1 * Taylor1(order)
    # Propagation and residuals
    propres!(resT1, IM, q0T1, jd0, params; buffer = bufferT1)
    # Covariance matrix in residuals space
    QT1 = nms(resT1)
    dQT1, d2QT1 = zero(QT1), zero(QT1)
    TS.differentiate!(dQT1, QT1)
    TS.differentiate!(d2QT1, dQT1)
    C = notout(resT1) * [d2QT1;;]
    Γ_ξ = inv(Symmetric(C))
    # Covariance matrix in cartesian space
    J = zero.(q0T1)
    for i in eachindex(J)
        TS.differentiate!(J[i], q0T1[i] - q00[i])
    end
    Γ_c = Symmetric(J * Γ_ξ * J')

    return Γ_c
end

# Return the Taylor expansions of the largest eigenpair of Γ
# See https://doi.org/10.1137/23M1551961
function machseries(Γ::AbstractMatrix{Taylor1{T}}, order::Int) where {T <: Real}
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
    λ, v = Taylor1(λs, order), [Taylor1(vs[i, :], order) for i in axes(vs, 1)]

    return λ, v
end

"""
    lov!

Definition of the line of variations as a differential equation.

!!! reference
    See equation (10.2) in section 10.1:
    - https://doi.org/10.1017/CBO9781139175371
"""
function lov!(dq, q, params, t)
    local resTN, resT1, IM, jd0, _params, bufferTN, bufferT1, order = params
    local Γ = lovcovariance(resTN, resT1, IM, cte(q), jd0, _params;
                            bufferTN, bufferT1, order)
    local λ, v = machseries(Γ, order)
    sqrt_λ = sqrt(λ)
    dq[1] = sqrt_λ * v[1]
    dq[2] = sqrt_λ * v[2]
    dq[3] = sqrt_λ * v[3]
    dq[4] = sqrt_λ * v[4]
    dq[5] = sqrt_λ * v[5]
    dq[6] = sqrt_λ * v[6]

    nothing
end

# Methods of `TaylorIntegration._allocate_jetcoeffs!` and `TaylorIntegration.jetcoeffs!`
# generated by @taylorize for lov!
function TaylorIntegration._allocate_jetcoeffs!(
        ::Val{lov!}, t::Taylor1{_T}, q::AbstractArray{Taylor1{_S}, _N},
        dq::AbstractArray{Taylor1{_S}, _N}, params
    ) where {_T <: Real, _S <: Number, _N}
    order = t.order
    local resTN, resT1, IM, jd0, _params, bufferTN, bufferT1, _ = params
    local Γ = lovcovariance(resTN, resT1, IM, cte(q), jd0, _params;
                            bufferTN, bufferT1, order)
    local (λ, v) = machseries(Γ, order)
    sqrt_λ = Taylor1(sqrt(constant_term(λ)), order)
    dq[1] = Taylor1(constant_term(sqrt_λ) * constant_term(v[1]), order)
    dq[2] = Taylor1(constant_term(sqrt_λ) * constant_term(v[2]), order)
    dq[3] = Taylor1(constant_term(sqrt_λ) * constant_term(v[3]), order)
    dq[4] = Taylor1(constant_term(sqrt_λ) * constant_term(v[4]), order)
    dq[5] = Taylor1(constant_term(sqrt_λ) * constant_term(v[5]), order)
    dq[6] = Taylor1(constant_term(sqrt_λ) * constant_term(v[6]), order)
    return TaylorIntegration.RetAlloc{Taylor1{_S}}(
        [sqrt_λ], [Array{Taylor1{_S}, 1}(undef, 0)], [Array{Taylor1{_S}, 2}(undef, 0, 0)],
        [Array{Taylor1{_S}, 3}(undef, 0, 0, 0)], [Array{Taylor1{_S}, 4}(undef, 0, 0, 0, 0)]
    )
end

function TaylorIntegration.jetcoeffs!(
        ::Val{lov!}, t::Taylor1{_T}, q::AbstractArray{Taylor1{_S}, _N},
        dq::AbstractArray{Taylor1{_S}, _N}, params,
        __ralloc::TaylorIntegration.RetAlloc{Taylor1{_S}}
    ) where {_T <: Real, _S <: Number, _N}
    order = t.order
    sqrt_λ = __ralloc.v0[1]
    local resTN, resT1, IM, jd0, _params, bufferTN, bufferT1, _ = params
    local Γ = lovcovariance(resTN, resT1, IM, cte(q), jd0, _params;
                            bufferTN, bufferT1, order)
    local (λ, v) = machseries(Γ, order)
    for ord = 0:order - 1
        ordnext = ord + 1
        TaylorSeries.sqrt!(sqrt_λ, λ, ord)
        TaylorSeries.mul!(dq[1], sqrt_λ, v[1], ord)
        TaylorSeries.mul!(dq[2], sqrt_λ, v[2], ord)
        TaylorSeries.mul!(dq[3], sqrt_λ, v[3], ord)
        TaylorSeries.mul!(dq[4], sqrt_λ, v[4], ord)
        TaylorSeries.mul!(dq[5], sqrt_λ, v[5], ord)
        TaylorSeries.mul!(dq[6], sqrt_λ, v[6], ord)
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

- `σmax::Real`: maximum (absolute) value of the LOV index (default: `3.0`).
- `lovorder::Int`: order of Taylor expansions wrt LOV index (default: `12`).
- `lovtol::Real`: absolute tolerance used to integrate the LOV (default: `1E-20`).
"""
function lineofvariations(IM::AbstractIMProblem{D, T}, params::Parameters{T};
                          σmax::Real = 3.0, lovorder::Int = 12,
                          lovtol::Real = 1E-20) where {D, T <: Real}
    # Unpack
    @unpack orbit = IM
    @unpack maxsteps, parse_eqs = params
    # Number of degrees of freedom
    Npar = dof(IM)
    # Jet transpot initial condition
    set_od_order(T, 2, Npar)
    # Refence epoch [julian date TDB]
    jd0 = epoch(orbit) + PE.J2000
    # Initial condition
    q00 = orbit()
    q0TN = q00 + sigmas(orbit) .* get_variables(T, 2)
    q0T1 = q00 + sigmas(orbit) .* Taylor1(lovorder)
    # Line of variations parameters
    resTN = init_optical_residuals(TaylorN{T}, IM)
    resT1 = init_optical_residuals(Taylor1{T}, IM)
    bufferTN = PropresBuffer(IM, q0TN, jd0, params)
    bufferT1 = PropresBuffer(IM, q0T1, jd0, params)
    lovparams = (resTN, resT1, IM, jd0, params, bufferTN, bufferT1, lovorder)
    # Taylor expansion of the line of variations
    _bwd_ = taylorinteg(lov!, q00, zero(T), -σmax, lovorder, lovtol, lovparams;
                        maxsteps, parse_eqs, dense = true)
    _fwd_ = taylorinteg(lov!, q00, zero(T), σmax, lovorder, lovtol, lovparams;
                        maxsteps, parse_eqs, dense = true)
    bwd = TaylorInterpolant{T, T, 2}(zero(T), _bwd_.t, _bwd_.p)
    fwd = TaylorInterpolant{T, T, 2}(zero(T), _fwd_.t, _fwd_.p)
    domain = (last(bwd.t), last(fwd.t))

    return LineOfVariations{D, T}(dynamicalmodel(IM), epoch(orbit), domain, bwd, fwd)
end