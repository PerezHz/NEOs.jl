module NEOsRecipesBaseExt

using RecipesBase
using NEOs: OpticalResidual, AdmissibleRegion, cte, ra, dec, boundary

@recipe function f(res::AbstractVector{OpticalResidual{T, U}}) where {T <: Real, U <: Number}
    seriestype --> :scatter
    return cte.(ra.(res)), cte.(dec.(res))
end

@recipe function f(A::AdmissibleRegion{T}; N::Int = 100, ρscale::Symbol = :linear) where {T <: Real}
    seriestype --> :path
    ps = map(t -> boundary(A, t), LinRange(0, 3, N));
    xs, ys = first.(ps), last.(ps)
    if ρscale == :log
        xs .= log10.(xs)
    end
    return xs, ys
end

end