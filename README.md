# NEOs.jl

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5152449.svg)](https://doi.org/10.5281/zenodo.5152449)
[![CI](https://github.com/PerezHz/NEOs.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/PerezHz/NEOs.jl/actions/workflows/CI.yml)
[![Coverage Status](https://coveralls.io/repos/github/PerezHz/NEOs.jl/badge.svg?branch=main)](https://coveralls.io/github/PerezHz/NEOs.jl?branch=main)
[![codecov](https://codecov.io/gh/PerezHz/NEOs.jl/branch/main/graph/badge.svg?token=F1IY79YP3J)](https://codecov.io/gh/PerezHz/NEOs.jl)

`NEOs.jl` is a Julia package for high-accuracy orbit determination and propagation of
Near-Earth Objects. `NEOs.jl` exploits jet transport techniques via
[TaylorIntegration.jl](https://github.com/PerezHz/TaylorIntegration.jl).

## Authors

- [Jorge A. Pérez](https://github.com/PerezHz),
Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México (UNAM)
- [Luis Benet](http://www.cicc.unam.mx/~benet/),
Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México (UNAM)
- [Luis Eduardo Ramírez Montoya](https://github.com/LuEdRaMo),
Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México (UNAM)

## Installation

The current version of this package may be installed in Julia pkg manager via:
```
] add NEOs
```

## Usage

The `pha` directory contains the `apophis.jl` script which performs an
orbit determination for asteroid (99942) Apophis from optical and radar astrometry. In order
to run this script, the environment corresponding to the `Project.toml` contained in the
`pha` directory has to be active and instantiated. This can be done, for example, by running
the following command from the `pha` directory:

``julia -e `import Pkg; Pkg.activate(); Pkg.instantiate()` # run this from the `pha` directory ``

Once the `pha` environment is active, this script may be called from the `pha` directory
with the default settings as:

`julia --project apophis.jl`

The `--help` option can be passed to see a list of the customizable settings

`julia --project apophis.jl --help`

`NEOs.propagate` also supports multi-threading:

`julia -t <number-of-threads> --project apophis.jl --help`

## Acknowledgments

We acknowledge financial support from UNAM-PAPIIT grants IG100819 and IG-101122, as well as
computing resources provided by LANCAD-UNAM-DGTIC-284.

## References

- Pérez-Hernández, J.A., Benet, L. Non-zero Yarkovsky acceleration for near-Earth asteroid
    (99942) Apophis. Commun Earth Environ 3, 10 (2022). https://doi.org/10.1038/s43247-021-00337-x
- Pérez-Hernández, Jorge A., & Benet, Luis. (2023).
    [PerezHz/TaylorIntegration.jl](https://github.com/PerezHzTaylorIntegration.jl):
    v0.14.2 (Version v0.14.2). Zenodo. https://doi.org/10.5281/zenodo.8104080
- Ramírez-Montoya, L.E., Pérez-Hernández, J.A. & Benet, L. Jet transport applications to
    the preliminary orbit determination problem. Celest Mech Dyn Astron 137, 16 (2025).
    https://doi.org/10.1007/s10569-025-10246-2