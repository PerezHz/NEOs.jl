module NEOsRecipesBaseExt

using RecipesBase
using PlanetaryEphemeris: TaylorInterpolant
using NEOs: OpticalResidual, AdmissibleRegion, LeastSquaresOrbit, cte, ra, dec, arboundary

@recipe function f(res::AbstractVector{OpticalResidual{T, U}}) where {T <: Real, U <: Number}
    seriestype --> :scatter
    return cte.(ra.(res)), cte.(dec.(res))
end

@recipe function f(A::AdmissibleRegion{T}; boundary::Symbol = :outer,
                   N::Int = 100, ρscale::Symbol = :linear) where {T <: Real}
    seriestype --> :path
    tmax = boundary == :outer ? 3 : 2
    ps = map(t -> arboundary(A, t, boundary, ρscale), LinRange(0, tmax, N))
    xs, ys = first.(ps), last.(ps)

    return xs, ys
end

@recipe function f(sol::U, t0::T, tf::T; N::Int = 100,
        projection::Symbol = :xyz) where {T <: Real,
        U <: Union{TaylorInterpolant, LeastSquaresOrbit}}
    seriestype --> :path
    ts = LinRange(t0, tf, N)
    rvs = Matrix{T}(undef, 6, N)
    for i in eachindex(ts)
        rvs[:, i] .= cte.(sol(ts[i])[1:6])
    end
    xs, ys, zs = rvs[1, :], rvs[2, :], rvs[3, :]

    if projection == :x
        return xs
    elseif projection == :y
        return ys
    elseif projection == :z
        return zs
    elseif projection == :xy
        return xs, ys
    elseif projection == :xz
        return xs, zs
    elseif projection == :yz
        return ys, zs
    elseif projection == :xyz
        return xs, ys, zs
    end
end

end