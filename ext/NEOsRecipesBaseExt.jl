module NEOsRecipesBaseExt

using RecipesBase
using TaylorIntegration: TaylorSolution
using NEOs: OpticalResidual, AdmissibleRegion, AbstractOrbit, cte, ra, dec,
      arboundary, body2observer

@recipe function f(res::AbstractVector{<:OpticalResidual})
    @series begin
        seriestype := :scatter
        α, δ = @. cte(ra(res)), cte(dec(res))
        return α, δ
    end
end

@recipe function f(A::AdmissibleRegion;
                   N = 100, ρscale = :linear,
                   Hs = [], Hcolor = :magenta, Hlinewidth = 1,
                   outer = true, outercolor = :red, outerlinewidth = 2,
                   inner = true, innercolor = :lime, innerlinewidth = 2)
    @assert isa(N, Int) && N > 0 "Number of points must be an integer greater than zero"
    @assert isa(ρscale, Symbol) && ρscale in (:linear, :log) "Possible values for \
        `ρscale` are: `:linear` and `:log`"
    # Outer boundary
    if outer
        @series begin
            label := ""
            color := outercolor
            seriestype := :path
            linewidth := outerlinewidth
            ts = LinRange(0, 3, N)
            ps = arboundary.(Ref(A), ts, Ref(:outer), Ref(ρscale))
            return first.(ps), last.(ps)
        end
    end
    # Inner boundary
    if inner
        @series begin
            label := ""
            z_order := 1
            color := innercolor
            seriestype := :path
            linewidth := innerlinewidth
            ts = LinRange(0, 2, N)
            ps = arboundary.(Ref(A), ts, Ref(:inner), Ref(ρscale))
            return first.(ps), last.(ps)
        end
    end
    # Shooting star limit
    if !isempty(Hs)
        @series begin
            label := ""
            z_order := 1
            color := Hcolor
            seriestype := :vline
            linewidth := Hlinewidth
            ρs = body2observer.(Ref(A), Hs)
            if ρscale === :linear
                return ρs
            else
                return log10.(ρs)
            end
        end
    end
end

@recipe function f(sol::Union{TaylorSolution, AbstractOrbit},
                   t0::Real, tf::Real; N = 100, projection = :xyz)
    @assert isa(N, Int) && N > 0 "Number of points must be an integer greater than zero"
    @assert isa(projection, Symbol) && projection in (:x, :y, :z, :xy, :xz,
        :yz, :xyz) "Possible values for `projection` are: `:x`, `:y`, `:z`, `:xy`, \
        `:xz`, `:yz` and `:xyz`"
    @series begin
        seriestype := :path
        ts = LinRange(t0, tf, N)
        rvs = Matrix{typeof(t0)}(undef, 6, N)
        for i in eachindex(ts)
            rvs[:, i] .= cte.(sol(ts[i])[1:6])
        end
        xs, ys, zs = rvs[1, :], rvs[2, :], rvs[3, :]
        if projection === :x
            return xs
        elseif projection === :y
            return ys
        elseif projection === :z
            return zs
        elseif projection === :xy
            return xs, ys
        elseif projection === :xz
            return xs, zs
        elseif projection === :yz
            return ys, zs
        elseif projection === :xyz
            return xs, ys, zs
        end
    end
end

end
