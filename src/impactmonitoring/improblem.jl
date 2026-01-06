"""
    IMProblem{D, T, O <: AbstractOrbit{D, T, T}} <: AbstractIMProblem{D, T}

An impact monitoring problem.

# Fields

- `orbit::O`: reference orbit.
- `target::ImpactTarget{T}`: target celestial body.
"""
struct IMProblem{D, T, O <: AbstractOrbit{D, T, T}} <: AbstractIMProblem{D, T}
    orbit::O
    target::ImpactTarget{T}
end

# Outer constructor
IMProblem(orbit::AbstractOrbit{D, T, T}, target::ImpactTarget{T}) where {D, T <: Real} =
    IMProblem{D, T, typeof(orbit)}(orbit, target)

# Print method for IMProblem
function show(io::IO, x::IMProblem)
    t = repeat(' ', 4)
    print(io,
        "Impact monitoring problem\n",
        t, rpad("Orbit:", 10), string(x.orbit), "\n",
        t, rpad("Target:", 10), string(x.target), "\n",
    )
end