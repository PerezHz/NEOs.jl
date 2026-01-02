# [Propagation](@id Propagation)

NEOs uses an adaptive-step Taylor method to propagate an initial condition in time; for details, see the [`TaylorIntegration.jl`](https://github.com/PerezHz/TaylorIntegration.jl) documentation. Below, we explain how to propagate an initial condition using Apophis as an example.

## Dynamical model

First, we need to choose a dynamical model. Currently, NEOs supports three different
models: `newtonian!`, `gravityonly!` and `nongravs!`. All of the former:
- treat the object of interest as a test particle with null mass,
- include the contribution from the Sun, the eight planets and the Moon,
- use multi-threading to improve the performance of some internal loops.

In the following table we summarize the dynamical effects included in each model.
For futher details, check the documentation of each function.

| Effect                                | `newtonian!` | `gravityonly!` | `nongravs!` |
| :------------------------------------ | :----------: | :------------: | :---------: |
| Newtonian accelerations               |      ✅      |      ✅      |      ✅      |
| Contribution from Ceres               |      ❌      |      ✅      |      ✅      |
| Post-Newtonian accelerations          |      ❌      |      ✅      |      ✅      |
| Figure-effects (oblateness)           |      ❌      |      ✅      |      ✅      |
| IAU 1976/1980 Earth orientation model |      ❌      |      ✅      |      ✅      |
| Moon's orientation model              |      ❌      |      ✅      |      ✅      |
| Non-gravitational accelerations       |      ❌      |      ❌      |      ✅      |

For now, let's use `newtonian!` as the dynamical model.
```@example Apophis
using NEOs

dynamics = newtonian!
```

## Initial condition

Next, we need an initial condition, i.e.:
- An initial time [julian date TDB],
- An initial barycentric cartesian state vector [au, au/day].

For instance, consider JPL's solution #220 for Apophis at May 5th, 2025:
```@example Apophis
using Dates

jd0 = datetime2julian(DateTime(2025, 5, 5)) # TDB
q00 = [1.250262833433083E-01, 8.819538568253468E-01, 3.308678718381778E-01,   # au
      -1.635141471921021E-02, 5.302965133786512E-03, 1.558030070573754E-03]   # au/day
```
## Propagate

Finally, propagate `q00` ten years into the future (well beyond the April 2029
close approach)
```@example Apophis
params = Parameters(maxsteps = 1_000, order = 25, abstol = 1E-20)
nyears_fwd = 10.0
fwd = propagate(dynamics, q00, jd0, nyears_fwd, params)
```

!!! tip "Parameters"
    Hola

!!! note
    With NEOs plotting recipes we can visualize Apophis' orbit and compare it
    with Earth's.
    ````@example Apophis
    using Plots

    t0 = dtutc2days(DateTime(2025, 5, 5))
    tf = dtutc2days(DateTime(2029, 4, 13))

    scatter(params.eph_su, t0, tf, xlabel = "x [au]", ylabel = "y [au]",
        label = "Sun", color = :gold, markersize = 8, markershape = :star5,
        N = 10, projection = :xy, aspect_ratio = 1)
    plot!(params.eph_ea, t0, tf, label = "Earth", color = :green, linewidth = 2,
        N = 1_000, projection = :xy)
    plot!(fwd, t0, tf, label = "Apophis", color = :deepskyblue, linewidth = 2,
        N = 1_000, projection = :xy)
    ````

## Jet transport

NEOs can exploit jet transport techniques to propagate a neighbourhood
of initial conditions around the nominal orbit. For instance, we can parametrize
the $3\sigma$ uncertainty box as follows:
```@example Apophis
using TaylorSeries

# JPL's solution #220 1-sigma uncertainties
ss = [6.56679380E-09, 1.15256858E-09, 3.02629146E-09,     # au
      4.18459496E-11, 8.86601155E-11, 5.21822853E-11]     # au/day
dq = set_variables("dx", numvars = 6, order = 2)
q0 = q00 .+ 3ss .* dq
```
and propagate it using the same syntax as before
```@example Apophis
fwdJT = propagate(newtonian!, q0, jd0, nyears_fwd, params)
```

Now, we have a polynomial approximation to the evolution of initial
conditions close to the nominal orbit, which we can sample a la
Monte Carlo
```@example Apophis
using Distributions

lb, up = fill(-3.0, 6), fill(3.0, 6)
dist = Product(Uniform.(lb, up))
dx = [rand(dist) for _ in 1:100]
```
and then evaluate in time, in order to visualize the growth of the
uncertainty region after the 2029 close approach.
```@example Apophis
t0 = dtutc2days(DateTime(2028, 4, 13))
tf = dtutc2days(DateTime(2030, 4, 14))
ts = LinRange(t0, tf, 25)

scatter(params.eph_su, t0, tf, xlabel = "x [au]", ylabel = "y [au]",
    label = "Sun", color = :gold, markersize = 8, markershape = :star5,
    N = 10, projection = :xy, aspect_ratio = 1)
plot!(params.eph_ea, t0, tf, label = "Earth", color = :green, linewidth = 2,
    N = 1_000, projection = :xy)
plot!(fwd, t0, tf, label = "Apophis", color = :deepskyblue, linewidth = 2,
    N = 1_000, projection = :xy)
for t in ts
    ps = map(d -> fwdJT(t)(d)[1:2], dx)
    scatter!(first.(ps), last.(ps), color = :blue, markersize = 3,
        markerstrokewidth = 0, label = "")
end
plot!()
```