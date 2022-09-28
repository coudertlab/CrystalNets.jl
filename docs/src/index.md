# CrystalNets.jl

`CrystalNets.jl` is a a Julia package for automatic detection and identification of
net topologies underlying crystalline materials.
Its inputs can be chemical files in any format recognized by [chemfiles](https://chemfiles.org/).

## Package installation

The installation follows the usual procedure. Start by downloading and installing [Julia](https://julialang.org/), v1.6 at least for `CrystalNets.jl`. This package was optimized with Julia v1.8 so performance and latency will be better on the more recent versions of Julia. Then, either

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
or proceed with the [Full installation](@ref) to obtain an executable.


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

#### Unknown nets

If the net is not recognized, its topological genome is displayed preceded by an "UNKNOWN"
mention, or "unstable" if the net is unstable:

```julia
julia> determine_topology("/path/to/new/material.cif")
UNKNOWN 2 1 2 -2 0 1 2 0 0 1 2 0 1 2 2 1 0

julia> determine_topology("/path/to/unstable/net.cif")
unstable 1 1 1 1 1 2 0 2 2 1
```

In both known and unknown cases, the result is a [`TopologyResult`](@ref).

#### Interpenetrating substructures

If the file contains multiple interpenetrating substructures, the result is a
`Vector{Tuple{Vector{Int}, TopologyResult}}`, where each entry is a tuple
`(vmap, result)` with:

- `vmap`: the list of vertices of the initial graph that were kept for this substructure.
  The initial graph is the one exported in .vtf as `input`. See also
  [`parse_chemfile`](@ref) and [`CrystalNets.Crystal`](@ref) for manipulations on the initial graph.
- `result`: the [`TopologyResult`](@ref) for this substructure.
For example:

```julia
julia> determine_topology("/path/to/intertwinned/structures.cif")
2-element Vector{Tuple{Vector{Int64}, TopologyResult}}:
 ([2, 3, 4, 6], pcu)
 ([1, 5, 7, 8], srs)
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

## Full installation

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

Compilation can take between fifteen minutes and an hour, depending on your hardware and the version of Julia.

The executable will be located in the "bin" subdirectory of the specified `INSTALLATION_PATH`, under the name "CrystalNets".

The executable can then simply be used on a chemical file:

```bash
$ CrystalNets /path/to/diamond.cif
dia
```

Run `CrystalNets --help` for the list of options available to the executable.

!!! tip
    In terms of performance, the compiled executable is the best option if you only want to identify a few structures from time to time. For intensive workloads with many structures to identify, it is best to use `CrystalNets.jl` as a Julia module through the
    [`determine_topology_dataset`](@ref) and [`guess_topology_dataset`](@ref) functions. The module is also the best option to perform more advanced analyses on the net in Julia, or to use the [`Options`](@ref) unavailable to the executable.
