# Initial position and velocity for Apophis via SPICE.jl
# See https://github.com/JuliaAstro/SPICE.jl

# Output from JPL solution #197 SPK file for Apophis via SPICE.jl
const x0_JPL_s197 = [-0.9633018953468989, 0.5100289806011301, 0.16528024397505386,
    -0.0071188720570829036, -0.012061235689040535, -0.0046695143453363164]

# Output from JPL solution #199 SPK file for Apophis via SPICE.jl
const x0_JPL_s199 = [-0.96330181481543, 0.5100291398744977, 0.16528030015848072,
    -0.007118874645519018, -0.0120612341605072, -0.004669513801491289]

@doc raw"""
    initialcond(dq::Vector=zeros(7); quadmath::Bool=false)

Returns the vector of initial conditions (3 positions + 3 velocities + ``A_2``
Yarkovsky coefficient) for Apophis.

# Arguments

- `dq`: `eltype(dq)` determines the type of the output. `dq` must be of size 6 or 7. 
- `quadmath`: Whether to use quadruple precision (`true`) or not.
"""
function initialcond(dq::Vector=zeros(7); quadmath::Bool=false)
    # Allocate memory for the initial conditions vector
    if quadmath
        # Use quadruple precision
        q0 = Vector{Float128}(undef, 7) 
    else
        # Use double precision
        q0 = Vector{Float64}(undef, 7) 
    end
    # Initial position and velocity (JPL #197 solution)
    # q0 .= x0_JPL_s199
    q0[1:6] .= x0_JPL_s197
    # A2 Yarkovsky coefficient
    q0[7] = 0.0 

    # Type of output
    S = typeof(q0[1]+dq[1])
    # Vector of initial conditions
    q = Vector{S}(undef, 7)

    # Convert q0 to type S

    # Fill positions + velocities
    for i in 1:6
        q[i] = q0[i] + dq[i]
    end
    # Fill Yarkovsky parameter 
    if length(dq) == 7
        q[7] = q0[7] + dq[7]
    elseif length(dq) == 6
        q[7] = q0[7] + 0dq[1]
    end

    return q
end
