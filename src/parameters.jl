@doc raw"""
    Parameters([::Parameters{T};] kwargs...) where {T <: Real}

A collection of the most important parameters in `NEOs.jl` functions.

# Common

-  `verbose::Bool`: whether to print output or not (default: `true`).

# Propagation

- `maxsteps::Int`: maximum number of steps for the integration (default: `500`).
- `order::Int`: order of Taylor expansions wrt time (default: 25).
- `abstol::T`: absolute tolerance used to compute propagation timestep
    (default: `1e-20`).
- `parse_eqs::Bool`: whether to use the specialized method of `jetcoeffs` or not
    (default: `true`).
- `bwdoffset/fwdoffset::T`: days to propagate beyond first (bwd) / last (fwd) observation
    (default: `0.5`).
- `coeffstol::T`: maximum size of the coefficients (default: `10.0`).

The following parameters are only used when propagating the Marsden et al. (1973)
nongravitational accelerations model:

- `marsden_coeffs::NTuple{3, T}`: Yarkovsky effect (A2), solar radiation pressure (A1),
    and normal component (A3) coefficients (default: `(0.0, 0.0, 0.0)`).
- `marsden_scalings::NTuple{3, T}`: scaling factors for the jet transport perturbation
    to the coefficients above (default: `(0.0, 0.0, 0.0)`).
- `marsden_radial::NTuple{5, T}`: radial function constants; i.e. coefficient
    `α`, normalizing distance `r₀` [au], and exponents `m`, `n` and `k` (default:
    `(1.0, 1.0, 2.0, 0.0, 0.0)`)

# Gauss Method

- `gaussorder::Int`: order of the jet transport perturbation (default: `2`).
- `safegauss::Bool`: if true, use the tracklets to search for Gauss triplets;
    otherwise, use the observations directly (default: `true`).
- `gaussntrip::Int`: number of Gauss triplets to try (default: `1`).
- `gaussτmin::T`: minimum desired total time span for `gaussmetric`
    in days (default: `1`).
- `gaussτmax::T`: maximum desired total time span for `gaussmetric`
    in days (default: `30`).
- `refscale::Symbol`: horizontal scale for ADAM refinement of Gauss preliminary
    orbit (default: `:log`).

# Minimization over the MOV

- `H_max::T`: maximum absolute magnitude (default: `34.5`).
- `a_max::T`: maximum semimajor axis (default: `100.0`).
- `adamiter::Int`: maximum number of iterations for `ADAM` optimizer (default: `200`).
- `adammode::Bool`: whether to perform ADAM iterations with all the observations
    (default: `true`).
- `adamQtol::T`: target function relative tolerance (default: `0.001`).
- `mmovproject::Bool`: whether to project the orbits onto the admissible region
    (default: `true`).
- `tsaorder::Int`: order of the jet transport perturbation (default: `6`).

# Least Squares

- `lspenalty::Bool`: whether to add a penalty to the orbit's eccentricity to
    the target function (default: `false`).
- `lsiter::Int`: maximum number of iterations for `leastsquares` (default: `10`).
- `lsQtol::T`: target function absolute tolerance (default: `1E-3`).
- `lsMtol::T`: Mahalanobis distance tolerance (default: `1E-3`).

# Jet Transport Least Squares

- `jtlsiter::Int`: maximum number of iterations for `jtls` (default: `5`).
- `jtlsorder::Int`: order of the jet transport perturbation in `jtls` (default: `6`).
- `significance::T`: chi-square significance level (default: `0.99`).
- `jtlsmask::Bool`: whether to use `isjtlsfit` to skip bad-conditioned
    preliminary orbits in `jtls` (default: `true`).
- `jtlsproject::Bool`: whether to project the orbits onto the admissible region
    (default: `false`).

# Outlier Rejection

- `outrej::Bool`: whether to perform outlier rejection during least squares
    iterations (default: `false`).
- `χ2_rec::T`: recovery threshold (default: `7.0`).
- `χ2_rej::T`: rejection threshold (default: `8.0`).
- `fudge::T`: rejection fudge term coefficient (default: `400.0`).
- `max_per::T`: maximum allowed rejection percentage (default: `10.0`).

# Physical properties

- `slope::T`: slope parameter (default: `0.15`).
- `albedo::T`: albedo (default: `0.14`).
- `density::T`: density in kg/m³ (default: `2_600.0`).
"""
@with_kw struct Parameters{T <: Real}
    # Common
    verbose::Bool = true
    # Propagation
    maxsteps::Int = 500
    order::Int = 25
    abstol::T = 1e-20
    parse_eqs::Bool = true
    bwdoffset::T = 0.5
    fwdoffset::T = 0.5
    coeffstol::T = 10.0
    marsden_coeffs::NTuple{3, T} = (0.0, 0.0, 0.0)
    marsden_scalings::NTuple{3, T} = (0.0, 0.0, 0.0)
    marsden_radial::NTuple{5, T} = (1.0, 1.0, 2.0, 0.0, 0.0)
    # Sun (earth) ephemeris
    eph_su::DensePropagation2{T, T} = selecteph(sseph, su)
    eph_ea::DensePropagation2{T, T} = selecteph(sseph, ea)
    # Gauss' Method
    gaussorder::Int = 2
    safegauss::Bool = true
    gaussntrip::Int = 1
    gaussτmin::T = 1.0
    gaussτmax::T = 30.0
    refscale::Symbol = :log
    # Minimization over the MOV
    H_max::T = 34.5
    a_max::T = 100.0
    adamiter::Int = 200
    adammode::Bool = true
    adamQtol::T = 0.001
    mmovproject::Bool = true
    tsaorder::Int = 6
    # Least Squares
    lspenalty::Bool = false
    lsiter::Int = 10
    lsQtol::T = 1E-3
    lsMtol::T = 1E-3
    # Jet Transport Least Squares
    jtlsiter::Int = 5
    jtlsorder::Int = 6
    significance::T = 0.99
    jtlsmask::Bool = true
    jtlsproject::Bool = false
    # Outlier Rejection
    outrej::Bool = false
    χ2_rec::T = 7.0
    χ2_rej::T = 8.0
    fudge::T = 400.0
    max_per::T = 10.0
    # Physical properties
    slope::T = 0.15
    albedo::T = 0.14
    density::T = 2_600.0
end
