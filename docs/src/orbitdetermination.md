# [Orbit determination](@id OrbitDetermination)

NEOs uses the procedure described in [Ramírez2025](@cite) to determine the orbit of a Near-Earth Object, either from scratch (initial orbit determination) or starting from a pre-existing orbit (orbit refinement). Below, we illustrate a typical orbit determination workflow by computing the orbit of 2024 XA1, an asteroid which hit Earth on December 3rd, 2024.

## Orbit determination problem

First, we declare the orbit determination problem that we want to solve, which consists of three basic elements [Milani2010](@cite):
- the dynamical model,
- the observations and,
- the error model.
[Selecting the dynamical model](@ref DynamicalModel) and [fetching the optical astrometry](@ref OpticalAstrometry) of a particular object have already been discussed in previous sections of this documentation. On the other hand, in NEOs an error model is composed of a weighting scheme and a debiasing scheme.

Currently, NEOs supports three weighting schemes:
- `UniformWeights`: all observations have the same 1 arcsec weight.
- `SourceWeights`: the weights are extracted from the observations' format.
- `Veres17`: assign weights according to [Veres2017](@cite);
and four debiasing schemes:
- `ZeroDebiasing`: all debiasing factors are set to zero.
- `SourceDebiasing`: the debiasing factors are extracted from the observations' format.
- `Farnocchia15`: assign debiasing factors according to [Farnocchia2015](@cite).
- `Eggl20`: assign debiasing factors according to [Eggl2020](@cite).

!!! note
    Most weighting and debiasing schemes do not define a custom value for each observatory. In that case, the default weight is 1 arcsec and the default debasing factor is 0 arcsec.

Thus, to declare the orbit determination problem use the following sintax:
```@example 2024XA1
using NEOs

optical = fetch_optical_ades("2024 XA1", MPC)

od = ODProblem(newtonian!, optical[2:5];
               weights = Veres17, debias = Eggl20)
```

## Initial orbit determination

Next, to compute a preliminary orbit of 2024 XA1, we use the `initialorbitdetermination` function:
```@example 2024XA1
params = Parameters(
    maxsteps = 100, order = 15, abstol = 1E-12,
    coeffstol = Inf, bwdoffset = 0.042, fwdoffset = 0.042,
    tsaorder = 2, adamiter = 500, adamQtol = 1e-5,
    jtlsorder = 2, jtlsmask = false, jtlsiter = 20,
    lsiter = 10, significance = 0.99, outrej = false,
)
orbit0 = initialorbitdetermination(od, params)
nothing # hide
```
!!! tip "Parameters"
    Some of the most important parameters used in `initialorbitdetermination` are:
    - `bwdoffset/fwdoffset`: how much time [in days] to integrate beyond the first/last observation.
    - `tsaorder/gaussorder/jtlsorder`: the degree of the jet transport expansions used in different sections of the orbit determination routine.
    - `lsiter/jtlsiter`: the maximum number of iterations for normal and jet transport least squares.
    - `significance`: chi-square significance level, used to decide whether an orbit is acceptable.

`orbit0` is an instance of `LeastSquaresOrbit`, a struct containing all the information from the orbit determination process. For instance, from `orbit0` we can compute common quality metrics, e.g.:
```@example 2024XA1
# Normalized mean square residual
nrms(orbit0)
# Signal to noise ratio
snr(orbit0)
# Minor Planet Center's Uncertainty Parameter
uncertaintyparameter(orbit0, params)
```

!!! note
    We can propagate `orbit0` ten hours into the future
    ```@example 2024XA1
    using PlanetaryEphemeris, Dates

    jd0 = epoch(orbit0) + J2000
    nyears_fwd = (datetime2julian(DateTime(2024, 12, 3, 16, 15)) - jd0) / yr
    q0 = orbit0(epoch(orbit0))

    fwd = NEOs.propagate(newtonian!, q0, jd0, nyears_fwd, params)
    ```
    to visualize the asteroid's orbit approaching Earth
    ```@example 2024XA1
    using Plots

    t0, tf = epoch(orbit0), epoch(orbit0) + 10/24

    plot(params.eph_ea, t0, tf, xlabel = "x [au]", ylabel = "y [au]",
        label = "Earth", color = :green, linewidth = 2, N = 1_000,
        projection = :xy, legend = :bottomleft)
    plot!(fwd, t0, tf, label = "2024 XA1 (4 obs)", color = :deepskyblue,
        linewidth = 2, N = 1_000, projection = :xy)
    ```

## Orbit refinement

Now, to refine the orbit of 2024 XA1, we simply include the remaining observations into the orbit determination problem:
```@example 2024XA1
# Update orbit determination problem and parameters
NEOs.update!(od, optical)
params = Parameters(params; jtlsorder = 4, outrej = true, χ2_rec = sqrt(9.21),
                    χ2_rej = sqrt(10), fudge = 100.0, max_per = 34.0)
# Orbit refinement
orbit1 = orbitdetermination(od, orbit0, params)
nothing # hide
```

!!! tip "Parameters"
    The `outrej` parameter turns on/off the outlier rejection, which is based on the algorithm described by [Carpino2003](@cite). In that context, `χ2_rec` and `χ2_rej` are the recovery and rejection thresholds, respectively; `fudge` corresponds to the coefficient of the fudge term and `max_per` indicates the maximum percentage of observations that can be rejected.

!!! note
    We can compare `orbit0` with `orbit1`
    ```@example 2024XA1
    jd1 = epoch(orbit1) + J2000
    q1 = orbit1(epoch(orbit1))

    fwd1 = NEOs.propagate(newtonian!, q1, jd1, nyears_fwd, params)

    t0, tf = epoch(orbit1), epoch(orbit1) + 10/24

    plot!(fwd1, t0, tf, label = "2024 XA1 (79 obs)", color = :red,
        linewidth = 2, N = 1_000, projection = :xy)
    ```

## References
```@bibliography
Pages = ["orbitdetermination.md"]
Canonical = false
```