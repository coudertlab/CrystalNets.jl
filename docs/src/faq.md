# Troubleshooting

This page covers a few issues that may occur when using CrystalNets.jl. If your problem is not mentioned here, you can [open an issue](https://github.com/coudertlab/CrystalNets.jl/issues/new).

## Can I use this without programming, or if I do not know Julia?

Sure! If you only want to identify the topology of a structure, the easiest way is to use the website interface at [https://progs.coudert.name/topology](https://progs.coudert.name/topology)

If you want to use the package in a programmatic way without learning Julia, it is possible to use it through Python. To do so, please check the dedicated [Python interface](@ref) tutorial.

## How can I check that the detected topology corresponds to my input?

This is the focus of the [Visualization](@ref) tutorial.

## How can I silence warnings or remove exports?

- From the module: see [`CrystalNets.toggle_warning`](@ref), [`CrystalNets.toggle_error`](@ref) and [`CrystalNets.toggle_export`](@ref)
- From the executable: use the `--no-warn`, `--no-error` and `--no-export` flags.

For exports, each `export_...` keyword argument to [`Options`](@ref CrystalNets.Options) can be individually set to `false` or `""` to silence this particular export. See also the paragraph on [export options](@ref exports).

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

## CrystalNets.jl incorrectly detects bonds or reports a periodic structure as "0-dimensional"

You can use the [Visualization](@ref) tutorial to identify which bonds are incorrectly
detected or missing.

Check that the input file is clean: if an atom is represented in multiple possible
positions for instance (as is common in CIF files), CrystalNets.jl could mistake them as
multiple atoms, which may break some heuristics of the bond-detection algorithm. Solvent
residues may also be incorrectly bonded in some circumstances and should be removed. The
`ignore_atoms` keyword argument in the [`Options`] may be useful in this regard.

You may want to provide a different `structure` keyword argument in the [`Options`](@ref)
taken among the possible instances of [`StructureType`](@ref). These can modify the
default heuristics for bond-guessing: for example, using `structure=StructureType.MOF`
gives metals a larger radius for the purpose of guessing bonds.

There are several customizable heuristics for bond-guessing available among the
[`Options`](@ref). Of particular interest are `cutoff_coeff`, `ignore_homoatomic_bonds` and
several of the options in the "Miscellaneous" section.

If nothing works, the last solution consists in providing an input file with explicit bonds
set, and use the `bonding=Input` keyword argument to [`Options`](@ref).

## The topology has a name given by ToposPro but CrystalNets.jl yields "UNKNOWN ..."

In CrystalNets.jl, a net is identified if the graph of the crystal is [isomorphic](https://en.wikipedia.org/wiki/Graph_isomorphism) to the net. This correspondence is exact and the algorithm can only fail in case of unstable nets, which are reported as such. To do so, CrystalNets.jl actually solves the more complex [graph canonicalization problem](https://en.wikipedia.org/wiki/Graph_canonization) which consists in finding a "genome" (a sequence of numbers) provably unique for each net and such that two isomorphic nets have the same genome. Each genome is then associated with a name. See [this article](https://doi.org/10/dbg89q) for more information on the implementation in the program Systre, from which CrystalNets.jl is derived.

ToposPro does not exactly solve the periodic graph isomorphism problem: instead, it identifies the graph of a crystal as a known net if both share a number of properties (coordination sequences up to a certain point, point symbols and vertex symbols, see [this article](https://doi.org/gd4tgf)). This is usually an excellent approximation, but it is mathematically unsound as there may exist two non-isormorphic nets sharing these three properties. In particular, this means that there cannot be a unique topological genome defined for these "nets" recognized by ToposPro. As a consequence, they cannot be used by CrystalNets.jl.

## [Why is the topology computed with Standard different from ToposPro's standard?](@id ToposProStandard)

While both algorithms usually align, the result may not be the same in all cases because they are not defined in the same way.

Most often, the difference will come from either:

- an oxygen atom having three (or more) bonds becomes a vertex for ToposPro but is removed by CrystalNets.jl. To solve this, use `split_O_vertex=false` in the [`Options`](@ref CrystalNets.Options).
- a paddle-wheel pattern is grouped into a single cluster by CrystalNets.jl but not by ToposPro. To solve this, use `detect_paddlewheels=false` in the [`Options`](@ref CrystalNets.Options).

## How can I do a database topology analysis with CrystalNets.jl?

The built-in way to do this consists in using the [`determine_topology_dataset`](@ref) function.
This function expects the path of a directory containing CIF files within (possibly in subdirectories).

## How can I directly access the genome of my structure instead of its name?

The result `x` of [`determine_topology`](@ref) is an [`InterpenetratedTopologyResult`](@ref). Its `length` gives the number of interpenetrated substructures. Each of its values, for instance `x[1]`, is a tuple `(topo, n)` meaning that the substructure is an `n`-fold catenated net of topology `topo`. `topo` itself is a [`TopologyResult`](@ref), which stores the result of a topology computation for possibly several clusterings. The [`TopologicalGenome`](@ref) associated to a given clustering can be extracted by indexing the [`TopologyResult`](@ref), for instance `t = topo[Clustering.SingleNodes]` (or simply `t = topo[:SingleNodes]`).

For example:

```jldoctest im19faq
julia> path_to_im19 = joinpath(dirname(dirname(pathof(CrystalNets))), "test", "cif", "IM-19.cif");

julia> result = determine_topology(path_to_im19; structure=StructureType.MOF)
AllNodes: rna
SingleNodes: bpq

julia> typeof(result)
InterpenetratedTopologyResult

julia> length(result)
1

julia> topo, n = only(result);

julia> n # catenation multiplicity
1

julia> topo
AllNodes: rna
SingleNodes: bpq

julia> typeof(topo)
TopologyResult

julia> genome_allnodes = topo[Clustering.AllNodes]
rna

julia> typeof(genome_allnodes)
TopologicalGenome
```

In case where all clusterings lead to the same genome, it can simply be accessed
by calling `first(topo)`.

Having obtained a [`TopologicalGenome`](@ref), the topological genome itself can accessed
by converting it to a `PeriodicGraph`:

```jldoctest im19faq
julia> genome = PeriodicGraph(genome_allnodes)
PeriodicGraph3D(6, PeriodicEdge3D[(1, 2, (0,0,0)), (1, 3, (0,0,0)), (1, 4, (0,0,0)), (1, 4, (0,0,1)), (1, 5, (0,0,0)), (1, 6, (0,0,0)), (2, 4, (0,0,1)), (2, 6, (-1,0,0)), (3, 4, (0,0,1)), (3, 5, (0,-1,0)), (4, 5, (0,0,0)), (4, 6, (0,0,0))])
```

In case of error during topology identification, the returned `genome` is a `PeriodicGraph{0}`.

The string representation of the genome is simply `string(genome)`:

``` im19faq
julia> string(genome)
"3 1 2 0 0 0 1 3 0 0 0 1 4 0 0 0 1 4 0 0 1 1 5 0 0 0 1 6 0 0 0 2 4 0 0 1 2 6 -1 0 0 3 4 0 0 1 3 5 0 -1 0 4 5 0 0 0 4 6 0 0 0"
```
