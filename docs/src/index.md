# CrystalNets.jl

`CrystalNets.jl` is a Julia package for automatic detection and identification of
net topologies underlying crystalline materials.
Its inputs can be chemical files in any format recognized by [chemfiles](https://chemfiles.org/).

To use it directly without any installation, simply use the website interface: [https://progs.coudert.name/topology](https://progs.coudert.name/topology)

To use it through Python, check the [Python interface](@ref) tutorial.

## Package installation

The installation follows the usual procedure. Start by downloading and installing [Julia](https://julialang.org/), v1.10 or more recent. Then, either

- open the Julia REPL and enter the package manager by typing `]`, then install `CrystalNets.jl` by entering:
  ```julia
  pkg> add CrystalNets
  ```
- alternatively, you can do it from a shell by executing:
  ```bash
  julia -e 'import Pkg; Pkg.add("CrystalNets")'
  ```

To use the package, open a REPL and enter

```julia
julia> using CrystalNets
```


## Quick usage as a module

To determine the topology of a structure stored in a file at location `path`, simply call

```julia
julia> determine_topology(path)
```

#### Known nets

If recognized, this yields the name of the net. For example:

```julia
julia> determine_topology("/path/to/diamond.cif")
Export of input is enabled: saving file at /tmp/input_diamond_0.vtf
Export of subnet_Auto is enabled: saving file at /tmp/subnet_Auto_diamond_0.vtf
dia
```

By default, the parsed input and the extracted underlying nets are exported as .vtf files
(see [Visualization](@ref)). To toggle on or off automatic exports, use
[`CrystalNets.toggle_export`](@ref), and similarly with [`CrystalNets.toggle_warning`](@ref)
for warnings.

#### Unnamed nets

If the net is not recognized as part of either the RCSR, EPINET or as a known zeolite, its topological genome is displayed preceded by an "UNKNOWN" mention:

```julia
julia> determine_topology("/path/to/new/material.cif")
UNKNOWN 2 1 2 -2 0 1 2 0 0 1 2 0 1 2 2 1 0
```

In both known and unknown cases, the result is an [`InterpenetratedTopologyResult`](@ref).

#### Interpenetrating substructures

If the file contains multiple interpenetrating substructures, each substructure and its catenation multiplicity can be extracted from the [`InterpenetratedTopologyResult`](@ref).

For example:

```julia
julia> x = determine_topology("/path/to/intertwinned/structures.cif")
2 interpenetrated substructures:
⋅ Subnet 1 → pcu
⋅ Subnet 2 → srs
```

#### Using options

[`Options`](@ref CrystalNets.Options) can be added as individual keyword arguments to the call. For instance:

```julia
julia> path_to_mof5 = joinpath(dirname(dirname(pathof(CrystalNets))), "test", "cif", "MOF-5.cif");

julia> determine_topology(path_to_mof5; structure=StructureType.MOF,
                                        clusterings=[Clustering.PE,Clustering.Standard],
                                        split_O_vertex=false)
PE: cab
Standard: fff
```

## Installation as an executable

To obtain an executable, `CrystalNets.jl` can be statically compiled.
To do so, run the following julia script after changing the `INSTALLATION_PATH` variable to the location where the `CrystalNets.jl` executable will be installed.
Note that this requires the [PackageCompiler](https://github.com/JuliaLang/PackageCompiler.jl/) module.

```julia
const INSTALLATION_PATH = "/fill/with/installation/path"

using PackageCompiler
using CrystalNets

const root = dirname(dirname(pathof(CrystalNets)))

create_app(root, INSTALLATION_PATH;
           precompile_statements_file=abspath(root, "src", "precompile_statements.jl"),
           filter_stdlibs=true)
```

Compilation can take a long time, depending on your hardware and the version of Julia.

The executable will be located in the "bin" subdirectory of the specified `INSTALLATION_PATH`, under the name "CrystalNets".

The executable can then simply be used on a chemical file:

```bash
$ CrystalNets /path/to/diamond.cif
dia
```

Run `CrystalNets --help` for the list of options available to the executable.

!!! tip
    For casual usage, using [the website](https://progs.coudert.name/topology) is the most convenient option, unless the nets you study are too big.

    For intensive workloads with many structures to identify, it is best to use `CrystalNets.jl` as a Julia module through the
    [`determine_topology_dataset`](@ref) function. The module is also the best option to perform more advanced analyses on the net in Julia, or to use the [`Options`](@ref) unavailable to the executable or the website.
