# NEOCP

NEOs can be used to sample the manifold of variations of an object listed in the [NEO Confirmation Page (NEOCP)](https://www.minorplanetcenter.net/iau/NEO/toconfirm_tabular.html) in a similar way to how services like [NEOScan (NEODyS)](https://newton.spacedys.com/neodys/NEOScan/) and [Scout (JPL)](https://cneos.jpl.nasa.gov/scout/#/) do. In particular, the `mcmov.jl` script located in the `scripts/neocp` folder of the NEOs repository, implements the algorithm described by [Spoto2018](@cite).

To run the aforementioned script, first activate the corresponding enviorment:
```julia
pkg> activate scripts/neocp
pkg> instantiate
```
Next, use the following command to print a help message explaining all the argments that the script can receive:
```
julia --project=scripts/neocp scripts/neocp/mcmov.jl -h
```

!!! tip "Paralellism"
    The `mcmov.jl` script can leverage both distributed and multi-threaded paralellism to speed up the computations. We highly recommend at least using multiple threads to run the script.

Below, we illustrate the results obtained with the `mcmov.jl` script for four objects that were listed in the NEOCP on January 17th, 2025: A11xtrt, JKu348, TF26AC0 and X88481. These results can be reproduced by substituting the corresponding designation in the `{}` of the following command:
```
julia -t 5 --project=scripts/neocp scripts/neocp/mcmov.jl -i {} -o {}.txt -s log --Nx 100 --Ny 100 --maxchi 5.0 --refine --scout --neoscan --neocp
```

![Manifold of variations sampling for NEOCP object A11xtrt](figures/A11xtrt.png)

![Manifold of variations sampling for NEOCP object JKu348](figures/JKu348.png)

![Manifold of variations sampling for NEOCP object TF26AC0](figures/TF26AC0.png)

![Manifold of variations sampling for NEOCP object X88481](figures/X88481.png)

!!! note
    In the figures above, the main differences between NEOs, NEOCP, NEOScan and Scout are due to the sampling strategy. On one hand, both NEOs and NEOScan follow [Spoto2018](@cite) and use a refined square grid. On the other hand, NEOCP uses a Markov Chain Monte Carlo algorithm to sample the manifold of variations, while Scout follows the systematic ranging approach described by [Farnocchia2015b](@cite).

## References
```@bibliography
Pages = ["neocp.md"]
Canonical = false
```