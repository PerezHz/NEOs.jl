# NEOs.jl

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5152449.svg)](https://doi.org/10.5281/zenodo.5152449)

`NEOs.jl` is a Near-Earth Object orbital propagator and
fitter in Julia. `NEOs.jl` exploits jet transport techniques via
[TaylorIntegration.jl](https://github.com/PerezHzTaylorIntegration.jl).

## Authors

- [Jorge A. Pérez](https://www.linkedin.com/in/perezhz),
Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México (UNAM)
- [Luis Benet](http://www.cicc.unam.mx/~benet/),
Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México (UNAM)

## Installation

The current development version of this package may be installed in Julia via:
```
] add NEOs
```

## Usage

The `apophis.jl` file in the `pha` directory contains an example script. This
script may be called as:

`julia --project=@. apophis.jl --help`

`NEOs.propagate` also supports multi-threading:

`julia -t <number-of-threads> --project=@. apophis.jl --help`

## Acknowledgments

We acknowledge financial support from UNAM-PAPIIT grant IG100819 and computing
resources provided by LANCAD-UNAM-DGTIC-284.

## References

- Pérez-Hernández, Jorge A., & Benet, Luis. (2020, October 13).
    [PerezHz/TaylorIntegration.jl](https://github.com/PerezHzTaylorIntegration.jl):
    v0.8.4 (Version v0.8.4). Zenodo. http://doi.org/10.5281/zenodo.4086166
