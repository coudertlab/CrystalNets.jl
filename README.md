# CrystalNets

[![Build Status](https://github.com/coudertlab/CrystalNets.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/coudertlab/CrystalNets.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://coudertlab.github.io/CrystalNets.jl/dev)

A Julia package for determining the topology of a net.

You can use this package through our website interface without any installation required: https://progs.coudert.name/topology !

You can also use it by [calling it from Python](https://molsim.info/CrystalNets.jl/dev/python).

The installation follows the usual procedure. Start by downloading and installing [Julia](https://julialang.org/), version 1.6 at least for `CrystalNets.jl`. This package was optimized with Julia version 1.8 so performance and latency will be better on the more recent versions of Julia. Then, either

- open the Julia REPL and enter the package manager by typing `]`, then install `CrystalNets.jl` by entering:
  ```julia
  pkg> add CrystalNets
  ```
- alternatively, you can do it from a shell by executing:
  ```bash
  julia -e 'import Pkg; Pkg.add("CrystalNets")'
  ```

Please read [the documentation](https://coudertlab.github.io/CrystalNets.jl/dev) for more
information on the use of CrystalNets.jl and [alternative installation as an executable](https://molsim.info/CrystalNets.jl/dev/#Full-installation) to reduce latency.

The companion article is published in SciPost Chemistry: [doi: 10.21468/SciPostChem.1.2.005](https://doi.org/10.21468/SciPostChem.1.2.005).
