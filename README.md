# PhD thesis

This is code associated to my phd thesis. I will also add
the thesis itself upon its publication.

## Usage

In the root directory of the repository, run

```sh
julia --project -e 'using Pkg; Pkg.instantiate()'
```

to install the dependencies. Then, use

```sh
julia --project find-indistinguishable-cm-configurations.jl
```

to scan for indistinguishable Calogero-Moser configurations
[as in this paper](https://aip.scitation.org/doi/abs/10.1063/1.4705269).

Use

```sh
julia --project root-system-equidistant-hyperplanes.jl
```

to compute the sizes and multiplicities of equidistant hyperplane
orbits under the automorphism group.

In both cases, one can edit the last line of the script to specify
exactly which root system the computation should be about.
