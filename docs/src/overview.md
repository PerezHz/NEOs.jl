# Overview

This documentation and the code in NEOs is divided into four blocks, each one representing a different aspect of the study of Near-Earth Objects:
- The [astrometry](@ref Astrometry) section is dedicated to fetch, read and write optical and radar observations of Near-Earth objects in a variety of formats.
- The [propagation](@ref Propagation) section is meant to integrate in time an initial condition using various dynamical models.
- The **orbit determination** section is used to determine the orbit of an asteroid, either from scratch or starting from a pre-existing orbit.
- The **impact monitoring** section deals with the computation of the probability that a given object collides with a planet.

!!! warning "Impact monitoring section"
    The impact monitoring code is currently being developed, so the corresponding section in this documentation is still pending.

## Parameters

NEOs exports a structure called `Parameters` which contains the most important and repeated arguments of the package's functions. To set the parameters in given session, simply run
```julia
using NEOs

# We can use the default values...
params = Parameters()
# or customize some parameters
params = Parameters(params; maxsteps = 1_000, order = 15, abstol = 1E-12)
```
To see a list of all the available parameters, check the documentation of `Parameters`.
```@docs
Parameters
```

!!! tip "Parameters"
    Across the documentation, the green boxes will indicate which parameters are relevant for a particular block of code.
