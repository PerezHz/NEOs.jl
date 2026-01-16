# NEOs.jl

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5152449.svg)]
(https://doi.org/10.5281/zenodo.5152449)
[![CI](https://github.com/PerezHz/NEOs.jl/actions/workflows/CI.yml/badge.svg?branch=main)]
(https://github.com/PerezHz/NEOs.jl/actions/workflows/CI.yml)
[![Coverage Status](https://coveralls.io/repos/github/PerezHz/NEOs.jl/badge.svg?branch=main)]
(https://coveralls.io/github/PerezHz/NEOs.jl?branch=main)
[![codecov](https://codecov.io/gh/PerezHz/NEOs.jl/branch/main/graph/badge.svg?token=F1IY79YP3J)]
(https://codecov.io/gh/PerezHz/NEOs.jl)
![v0.20.0](https://img.shields.io/badge/version-v0.20.0-green.svg)

A [Julia](http://julialang.org) package for high-accuracy propagation, orbit determination
and impact monitoring of Near-Earth Objects. It exploits jet transport techniques via
[TaylorSeries.jl](https://github.com/JuliaDiff/TaylorSeries.jl) and [TaylorIntegration.jl]
(https://github.com/PerezHz/TaylorIntegration.jl).

### Authors

- [Jorge A. Pérez-Hernández](https://github.com/PerezHz),
    Minor Planet Center, Harvard & Smithsonian Center for Astrophysics
- [Luis Benet](http://www.cicc.unam.mx/~benet/),
    Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México (UNAM)
- [Luis E. Ramírez Montoya](https://github.com/LuEdRaMo),
    Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México (UNAM)

### License

`NEOs` is licensed under the [MIT "Expat" license](https://github.com/PerezHz/NEOs.jl/blob/main/LICENSE.md).

### Installation

NEOs is a [registered package](https://julialang.org/packages/), and is simply installed by running

```julia
pkg> add NEOs
```

!!! compat "Julia 1.10 or later is required"
    NEOs requires Julia 1.10 or later.

!!! info "Tested Julia versions"
    NEOs is currently tested on Julia 1.10 and 1.12.

### Acknowledgments

We acknowledge financial support from UNAM-PAPIIT grants IG100819 and IG-101122, as well as
computing resources provided by LANCAD-UNAM-DGTIC-284.

### Papers using NEOs

```@bibliography
Pages = []
Canonical = false

Perez-Hernandez2022
Ramírez2025
```