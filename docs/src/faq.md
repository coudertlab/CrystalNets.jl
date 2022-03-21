# Troubleshooting

This page covers a few issues that may occur when using CrystalNets.jl. If your problem is not mentioned here, you can [open an issue](https://github.com/coudertlab/CrystalNets.jl/issues/new).

## How can I silence warnings or remove exports?

- From the module: see [`CrystalNets.toggle_warning`](@ref) and [`CrystalNets.toggle_export`](@ref)
- From the executable: use the `--no-warn` and `--no-export` flags.

For exports, each `export_...` keyword argument to [`Options`](@ref CrystalNets.Options) can be individually set to `false` or `""` to silence this particular export. See also the paragraph on [Export options](@ref).

## How can I check that the detected topology corresponds to my input?

This is the focus of the [Visualization](@ref) tutorial.

## How can I choose a particular clustering algorithm?

Pass the appropriate option from [`Clustering`](@ref) to the `clusterings`
keyword argument of [`Options`](@ref CrystalNets.Options), or use the `-c` flag from the executable.
For example, to compute the topology of UiO-66 with the Points-of-Extension (PE) clustering algorithm, do

```@meta
DocTestSetup = quote
    using CrystalNets
    import CrystalNets: Options, Clustering, Bonding, StructureType
    const PeriodicGraphs = CrystalNets.PeriodicGraphs
    using .PeriodicGraphs

    CrystalNets.toggle_export(false)
    CrystalNets.toggle_warning(false)
end
```

```jldoctest uio66
julia> path_to_uio66 = joinpath(dirname(dirname(pathof(CrystalNets))), "test", "cif", "UiO-66.cif");

julia> determine_topology(path_to_uio66; structure=StructureType.MOF, clusterings=[Clustering.PE])
PE: ubt
```

The `clusterings` keyword argument accepts a list of [`Clustering`](@ref)s you can check multiple clusterings at once while factoring useless computations. For example:

```jldoctest uio66
julia> determine_topology(path_to_uio66; structure=StructureType.MOF, clusterings=[Clustering.Standard, Clustering.Auto])
Standard: xbi
AllNodes, SingleNodes: fcu
```

Note that `Auto` is equivalent to both `AllNodes` and `SingleNodes` when the [`StructureType`](@ref) is set to `MOF`.

See [below](@ref ToposProStandard) for the reason why the `Standard` topology is **xbi** and not ToposPro's "3,4,8T15".

## The topology has a name given by ToposPro but CrystalNets.jl yields "UNKNOWN ..."

In CrystalNets.jl, a net is identified if the graph of the crystal is [isomorphic](https://en.wikipedia.org/wiki/Graph_isomorphism) to the net. This correspondence is exact and the algorithm can only fail in case of unstable nets, which are reported as such. To do so, CrystalNets.jl actually solves the more complex [graph canonicalization problem](https://en.wikipedia.org/wiki/Graph_canonization) which consists in finding a "genome" (a sequence of numbers) provably unique for each net and such that two isomorphic nets have the same genome. Each genome is then associated with a name. See [this article](https://doi.org/10/dbg89q) for more information on the implementation in the program Systre, from which CrystalNets.jl is derived.

ToposPro does not exactly solve the periodic graph isomorphism problem: instead, it identifies the graph of a crystal as a known net if both share a number of properties (coordination sequences up to a certain point, point symbols and vertex symbols, see [this article](https://doi.org/gd4tgf)). This is usually an excellent approximation, but it is mathematically unsound as there may exist two non-isormorphic nets sharing these three properties. In particular, this means that there cannot be a unique topological genome defined for these "nets" recognized by ToposPro. As a consequence, they cannot be used by CrystalNets.jl.

## [Why is the topology computed with Standard different from ToposPro's standard?](@id ToposProStandard)

While both algorithms usually align, the result may not be the same in all cases because they are not defined in the same way.

Most often, the difference will come from either:

- an oxygen atom having three (or more) bonds becomes a vertex for ToposPro but is removed by CrystalNets.jl. To solve this, use `split_O_vertex=false` in the [`Options`](@ref CrystalNets.Options).
- a paddle-wheel pattern is grouped into a single cluster by CrystalNets.jl but not by ToposPro. To solve this, use `detect_paddlewheels=false` in the [`Options`](@ref CrystalNets.Options).

## How can I do a database topology analysis with CrystalNets.jl?

The built-in way to do this consists in using the [`determine_topology_dataset`](@ref) function, or [`guess_topology_dataset`](@ref) in some cases.
These functions expect the path of a directory containing CIF files within (possibly in subdirectories).
