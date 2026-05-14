"""
    LineOfVariationsBuffer{T <: Real} <: AbstractBuffer

Pre-allocated memory for [`lineofvariations`](@ref).

# Fields

- `t0::T`: reference epoch [TDB days since J2000].
- `sun::Vector{T}`: Sun barycentric cartesian state vector at `t0` [au, au/day].
- `scalings::Vector{T}`: covariance matrix scaling factors.
- `resTN::Vector{OpticalResidual{T, TaylorN{T}}}`: buffer for `TaylorN{T}` residuals.
- `bufferTN::PropresBuffer{T, TaylorN{T}, T}`: buffer for `TaylorN{T}` propagations.
"""
struct LineOfVariationsBuffer{T <: Real} <: AbstractBuffer
    t0::T
    sun::Vector{T}
    scalings::Vector{T}
    resTN::Vector{OpticalResidual{T, TaylorN{T}}}
    bufferTN::PropresBuffer{T, TaylorN{T}, T}
end

get_order(x::LineOfVariationsBuffer) = get_order(x.bufferTN.prop.cache.x[1][0])

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
    @unpack eph_su = params
    # Set jet transport order
    Ndof = dof(IM)
    set_od_order(T, lovorder, Ndof)
    # Refence epoch [julian date TDB]
    t0 = epoch(orbit)
    jd0 = t0 + PE.J2000
    # Sun's state vector at jd0
    sun = eph_su(t0)
    # Covariance matrix scaling factors
    scalings = fill(1E-8, 6)
    if Ndof == 9
        scalings = vcat(scalings, params.marsden_scalings...)
    end
    # Initial condition
    q00 = orbit()
    q0TN = q00 + sigmas(orbit) .* get_variables(T, lovorder)
    # Vectors of residuals
    resTN = init_optical_residuals(TaylorN{T}, IM)
    # Propagation and residuals buffers
    bufferTN = PropresBuffer(IM, q0TN, jd0, params)

    return LineOfVariationsBuffer{T}(t0, sun, scalings, resTN, bufferTN)
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

# Hessian matrix without evaluation
# TO DO: move this function to TaylorSeries
function hessianmatrix(f::TaylorN{T}) where {T <: Number}
    numVars = get_numvars()
    hess = Matrix{TaylorN{T}}(undef, numVars, numVars)
    ntup = zeros(Int, numVars)
    for j in axes(hess, 2)
        for i in axes(hess, 1)
            ntup .= 0
            ntup[i] += 1
            ntup[j] += 1
            hess[i, j] = differentiate(f, tuple(ntup...))
        end
    end
    return hess
end

# Coordinate transformations
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

# Return the jet transport expansion of the covariance matrix
# of an orbit in the required coordinates
function covariance(
        IM::AbstractIMProblem{D, T}, coord00::Vector{T}, coord::Symbol,
        buffer::LineOfVariationsBuffer{T}, params::Parameters{T}
    ) where {D, T <: Real}
    # Unpack
    @unpack t0, sun, scalings, resTN, bufferTN = buffer
    # Refence epoch [MJDTDB, JDTDB]
    mjd0 = t0 + MJD2000
    jd0 = t0 + PE.J2000
    # Order with respect to LOV index
    order = get_order(buffer)
    # Jet transpot initial condition
    car00 = lovtransform(mjd0, coord00, sun, Val(coord), Val(:cartesian))
    carTN = car00 + scalings .* get_variables(T, order)
    # TaylorN propagation and residuals
    propres!(resTN, IM, carTN, jd0, params; buffer = bufferTN)
    # Covariance matrix in residuals space
    QTN = nms(resTN)
    CTN_res = notout(resTN) * hessianmatrix(QTN)
    ΓTN_res = inv(Symmetric(CTN_res))
    # Covariance matrix in coordinate space
    coordTN = lovtransform(mjd0, carTN, sun, Val(:cartesian), Val(coord))
    ΓTN_coord = Symmetric(project(coordTN, ΓTN_res))
    return ΓTN_coord
end

# Return the Taylor expansion of the vector field associated to the
# weak direction of Γ; see https://doi.org/10.1137/23M1551961
function weakfield(Γ::AbstractMatrix{Taylor1{T}}, order::Int) where {T <: Real}
    # 0th-order eigenvalues
    Γ0 = Symmetric(constant_term.(Γ))
    E0 = eigen(Γ0)
    λ = Taylor1(E0.values[end], order)
    v = Taylor1.(E0.vectors[:, end], order)
    # Coefficients matrix
    E = [0 E0.vectors[:, end]'; E0.vectors[:, end] E0.values[end]*I-Γ0]
    # Main loop
    for k in 1:order
        y, z = zeros(T, length(v)), zero(T)
        for l in 0:k-1
            y .+= getcoeff.(Γ, k-l) * getcoeff.(v, l)
            if l >= 1
                y .-= getcoeff.(v, k-l) * λ[l]
                if l < k - 1
                    z += l * getcoeff.(v, k-l)' * getcoeff.(v, l)
                end
            end
        end
        Ek = E \ [-z/(2k); y]
        λ[k] = Ek[1]
        for i in eachindex(v)
            v[i][k] = Ek[i+1]
        end
    end
    F = sqrt(λ) * v
    return F
end

"""
    lov!

Definition of the line of variations as a differential equation.

!!! reference
    See equation (10.2) in section 10.1:
    - https://doi.org/10.1017/CBO9781139175371
"""
function lov! end

# Methods of `TaylorIntegration._allocate_jetcoeffs!` and `TaylorIntegration.jetcoeffs!`
# generated by @taylorize for lov!
function TaylorIntegration._allocate_jetcoeffs!(
        ::Val{lov!}, t::Taylor1{_T}, q::AbstractArray{Taylor1{_S}, _N},
        dq::AbstractArray{Taylor1{_S}, _N}, params
    ) where {_T <: Real, _S <: Number, _N}
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
    local IM, coord, buffer, _params_ = params
    @unpack t0, sun, scalings = buffer
    local mjd0 = t0 + MJD2000
    local ΓTN_coord = covariance(IM, cte(q), coord, buffer, _params_)
    for ord = 0:order - 1
        ordnext = ord + 1
        carT1 = lovtransform(mjd0, reduceorder.(q, ord), sun, Val(coord), Val(:cartesian))
        dx = @. (carT1 - cte(carT1)) / scalings
        ΓT1_coord = ΓTN_coord(dx)
        F = weakfield(ΓT1_coord, ord)
        for i in eachindex(F)
            TaylorSeries.identity!(dq[i], F[i], ord)
            TaylorIntegration.solcoeff!(q[i], dq[i], ordnext)
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
- `lovorder::Int`: order of Taylor expansions wrt LOV index (default: `4`).
- `lovtol::Real`: absolute tolerance used to integrate the LOV (default: `1E-8`).
- `lovsteps::Int`: maximum number of steps for the integration (default: `100`).
"""
function lineofvariations(IM::AbstractIMProblem{D, T}, params::Parameters{T};
                          coord::Symbol = :cartesian, σmax::Real = 3.0,
                          lovorder::Int = 4, lovtol::Real = 1E-8,
                          lovsteps::Int = 100) where {D, T <: Real}
    # Line of variations buffer
    buffer = LineOfVariationsBuffer(IM, lovorder, params)
    # Unpack
    @unpack orbit = IM
    @unpack t0, sun = buffer
    # Refence epoch [days since J2000 / MJDTDB]
    mjd0 = t0 + MJD2000
    # Initial condition
    q00 = orbit()
    coord00 = lovtransform(mjd0, q00, sun, Val(:cartesian), Val(coord))
    # TaylorIntegration cache
    lovparams = (IM, coord, buffer, params)
    cache = init_cache(Val(true), zero(T), coord00, lovsteps, lovorder, lov!,
                       lovparams; parse_eqs = true)
    # Taylor expansion of the line of variations
    _bwd_ = taylorinteg!(Val(true), lov!, coord00, zero(T), -σmax, lovtol,
                        cache, lovparams; maxsteps = lovsteps)
    bwd = TaylorInterpolant{T, T, 2}(zero(T), _bwd_.t, _bwd_.p)
    _fwd_ = taylorinteg!(Val(true), lov!, coord00, zero(T), σmax, lovtol,
                        cache, lovparams; maxsteps = lovsteps)
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