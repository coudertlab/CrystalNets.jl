# CrystalNets

<!---
[![Build Status](https://travis-ci.com/Liozou/CrystalNets.jl.svg?branch=master)](https://travis-ci.com/Liozou/CrystalNets.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/Liozou/CrystalNets.jl?svg=true)](https://ci.appveyor.com/project/Liozou/CrystalNets-jl)
-->

A Julia package for determining the topology of a net.
Currently only works for 3-periodic nets.

## <a name="installation"></a>Installation

Install julia and add this module either from a console with:

```bash
julia -e 'import Pkg; Pkg.add(url="https://github.com/coudertlab/CrystalNets.jl")'
```

or from the julia REPL by typing

```
] add https://github.com/coudertlab/CrystalNets.jl
```

## Use as a module

CrystalNets can read any chemical file parseable by [Chemfiles](https://chemfiles.org/chemfiles/latest/formats.html#list-of-supported-formats). In order to be correctly interpreted, the file should only contain information about the crystal (or crystalline framework), which must be a connected 3-periodic structure.

Parsing can be done with the dedicated function `parse_chemfile`. A parsed file will be stored as an intermediate object of type `Crystal`, as in the following session:

```julia
julia> using CrystalNets

julia> crystal = parse_chemfile("/path/to/diamond.cif");
```

`Crystal` objects retain information on the input geometry, such as the cell and the positions. These are irrelevant to the elucidation of the topological structure however, where only the underlying net matters. A `Crystal` must thus be converted to a net, using the type `CrystalNet`:

```julia
julia> net = CrystalNet(crystal);
```

Finally, the topological genome can be extracted with the dedicated `topological_genome` function, which yields a [`PeriodicGraph`](https://github.com/Liozou/PeriodicGraphs.jl) that uniquely represents the net. Its string version is the dense genome that is characteristic of the net: for instance, with the diamond net:

```julia
julia> string(topological_genome(net))
"3 1 2 0 0 0 1 2 0 0 1 1 2 0 1 0 1 2 1 0 0"
```

This string can be fed to the `reckognize_topology` function in order to determine the RCSR name of the topology:

```julia
julia> reckognize_topology("3 1 2 0 0 0 1 2 0 0 1 1 2 0 1 0 1 2 1 0 0")
"dia"
```

All these steps can be summarised by executing

```julia
julia> reckognize_topology(topological_genome(CrystalNet(parse_chemfile(FILE))))
```

where `FILE` is the path to the crystal in a reckognizable chemical file format.

## Use as an executable

### Without additional installation

Once installed following the [Installation](#installation) section, the following script will yield the path to the local installation of the module:

```julia
julia> using CrystalNets

julia> path = pathof(CrystalNets)
"/the/path/to/CrystalNets/src/CrystalNets.jl"

julia> normpath(path, "..", "..")
"/the/path/to/CrystalNets/"

```

Once you have this path, execute the CrystalNets.jl file using the CrystalNet directory as project from a shell with julia. For instance:

```
$ julia --project=/the/path/to/CrystalNets/ /the/path/to/CrystalNets/src/CrystalNets.jl
Missing a CRYSTAL_FILE.
See --help for usage.

$ julia --project=/the/path/to/CrystalNets/ /the/path/to/CrystalNets/src/CrystalNets.jl -c mof /path/to/HKUST-1.cif
Clustering of vertices represented represented at /path/to/clusters.pdb
tbo
```

### Full installation

To create a proper executable (which will execute faster than the previous option), follow the previous [Installation](#installation) section and then run the following julia script after changing the `INSTALLATION_PATH` variable to the location where the CrystalNets executable will be installed. Note that this requires the [PackageCompiler](https://github.com/JuliaLang/PackageCompiler.jl/) module, which can be installed by replacing "CrystalNets" with "PackageCompiler" in the [Installation](#installation) section.

```julia
INSTALLATION_PATH = "/fill/with/installation/path"

using PackageCompiler
using CrystalNets

root = normpath(pathof(CrystalNets), "..", "..")

create_app(root, INSTALLATION_PATH; precompile_statements_file=abspath(root, "src", "precompile.jl"))
```

The executable will be located in the "bin" subdirectory of the specified `INSTALLATION_PATH`, under the name "CrystalNets".

The executable can then simply be used on a chemical file:

```
$ CrystalNets /path/to/diamond.cif
```

will yield

```
dia
```

Run `CrystalNets --help` for the full list of options.
