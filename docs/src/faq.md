# Troubleshooting

This page covers a few issues that may occur when using CrystalNets.jl. If your problem is not mentioned here, you can [open an issue](https://github.com/coudertlab/CrystalNets.jl/issues/new).

## Can I use this without programming, or if I do not know Julia?

Sure! If you only want to identify the topology of a structure, the easiest way is to use the website interface at [https://progs.coudert.name/topology](https://progs.coudert.name/topology)

If you want to use the package in a programmatic way through Python, please check the dedicated [Python interface](@ref) tutorial.

## How can I check that the detected topology corresponds to my input?

This is the focus of the [Visualization](@ref) tutorial.
You may also want to check the [Vertex mapping](@ref mapping) tutorial.

On the website, a visualization panel automatically opens once the topology is computed, which enables checking the correspondence between the topology and the input.

## How can I silence warnings or remove exports?

- From the module: see [`CrystalNets.toggle_warning`](@ref), [`CrystalNets.toggle_error`](@ref) and [`CrystalNets.toggle_export`](@ref)
- From the executable: use the `--no-warn`, `--no-error` and `--no-export` flags.

For exports, each `export_...` keyword argument to [`Options`](@ref) can be individually set to `false` or `""` to silence this particular export. See also the paragraph on [export options](@ref exports).

## How can I choose a particular clustering algorithm?

Pass the appropriate option from [`Clustering`](@ref) to the `clusterings`
keyword argument of [`Options`](@ref), or use the `-c` flag from the executable.
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
`ignore_atoms` keyword argument in the [`Options`](@ref) may be useful in this regard.

You may want to provide a different `structure` keyword argument in the [`Options`](@ref)
taken among the possible instances of [`StructureType`](@ref). These can modify the
default heuristics for bond-guessing: for example, using `structure=StructureType.MOF`
gives metals a larger radius for the purpose of guessing bonds.

There are several customizable heuristics for bond-guessing available among the
[`Options`](@ref). Of particular interest are `cutoff_coeff`, `ignore_homoatomic_bonds` and
several of the options in the "Miscellaneous" section.

If nothing works, the last solution consists in providing an input file with explicit bonds
set, and use the `bonding=Bonding.Input` keyword argument to [`Options`](@ref).

## What is the meaning of a topology name such as "**cds**", "MFI", "sqc2423" or "3-xcoknootigbe" ?

The names of the topology come from different possible sources:
1. The Reticular Chemistry Structural Resource, [RCSR](http://rcsr.net/nets): these are three-character names with bold lowercase letters, possibly followed by one-letter extensions preceded by hyphens. For instance: **pcu**, **bcu-x**, **mtn-e-a**.
2. The International Zeolite Association Structure Commission, [IZA-SC](https://europe.iza-structure.org/IZA-SC/ftc_table.php): these are three-character long names with uppercase letters, possibly preceded by a hyphen or an asterisk. For instance: FAU, LTA, -ITV.
3. The Euclidean Patterns In Non-Euclidean Tilings, [EPINET](https://epinet.anu.edu.au/systre_net_searches/new) project: these are names that start by the three characters "sqc" followed by a number. For instance: sqc2, sqc800, sqc14645.
4. If the topology does not appear in any of the previous databases, CrystalNets.jl attributes a custom name made of an integer, a hyphen, and a sequence of twelve lowercase letters. The integer is the dimensionality of the net, and the sequence of letters is the hash of the topological genome. For instance: 2-rgavcufvusao, 3-bytqigwigfea, 3-wbiudmokgaek.

Note that the custom CrystalNets.jl hash is not guaranteed to be unique: since it is computed by reducing ("hashing") the topological genome into twelve letters, there could be two topologies with different genomes but the same reduced name. Comparing the topological genomes (the sequence of numbers after the name) is the only guaranteed way to make sure that two topologies are different or the same. In practice, the probability of a hash collision is extremely low and we checked that no two topologies shared the same name among more than 700 000 unique topologies, taken from different sources including the Materials Project and hypothetical zeolite databases.

When a topology is present in both the RCSR and IZA-SC, both names are given. Otherwise, the only displayed name is that given by the first database in the previous order of the sources. For example, the zeolite net **fau** in the RCSR is called FAU by the IZA-SC, sqc13519 by EPINET and 3-cqujdekuisif by the custom naming scheme of CrystalNets.jl, but it is only printed as "**fau**, FAU".

## The topology has a name given by ToposPro but CrystalNets.jl yields another name

In CrystalNets.jl, a net is identified if the graph of the crystal is [isomorphic](https://en.wikipedia.org/wiki/Graph_isomorphism) to the net. This correspondence is exact. To do so, CrystalNets.jl actually solves the more complex [graph canonicalization problem](https://en.wikipedia.org/wiki/Graph_canonization) which consists in finding a "genome" (a sequence of numbers) provably unique for each net and such that two isomorphic nets have the same genome. Each genome is then associated with a name. See [this article](https://doi.org/10/dbg89q) for more information on the implementation in the program Systre, from which CrystalNets.jl is derived.

ToposPro does not exactly solve the periodic graph isomorphism problem: instead, it identifies the graph of a crystal as a known net if both share a number of properties (coordination sequences up to a certain point, point symbols and vertex symbols, see [this article](https://doi.org/gd4tgf)). This is usually an excellent approximation, but it is mathematically unsound as there may exist two non-isormorphic nets sharing these three properties. In particular, this means that there cannot be a unique topological genome defined for these "nets" recognized by ToposPro. As a consequence, they cannot be used by CrystalNets.jl.

## [Why is the topology computed with Standard different from ToposPro's standard?](@id ToposProStandard)

While both algorithms usually align, the result may not be the same in all cases because they are not defined in the same way.

Most often, the difference will come from either:

- an oxygen atom having three (or more) bonds becomes a vertex for ToposPro but is removed by CrystalNets.jl in MOFs. To use ToposPro's default behaviour, set `split_O_vertex=false` in the [`Options`](@ref).
- a paddle-wheel pattern is grouped into a single cluster by CrystalNets.jl but not by ToposPro. To use ToposPro's default behaviour, use `detect_paddlewheels=false` in the [`Options`](@ref).

## How can I do a database topology analysis with CrystalNets.jl?

The built-in way to do this consists in using the [`determine_topology_dataset`](@ref) function.
This function expects the path of a directory containing CIF files within (possibly in subdirectories).

## [How can I directly access the genome of my structure instead of its name?](@id genomefromname)

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

```jldoctest im19faq
julia> string(genome)
"3 1 2 0 0 0 1 3 0 0 0 1 4 0 0 0 1 4 0 0 1 1 5 0 0 0 1 6 0 0 0 2 4 0 0 1 2 6 -1 0 0 3 4 0 0 1 3 5 0 -1 0 4 5 0 0 0 4 6 0 0 0"
```

In case the topology is unique, or unique for the given clustering, you can use the shortcut function [`one_topology`](@ref).

## Can I identify which input atom maps to which vertex of the returned genome?

Yes, that is the purpose of the `track_mapping` option, detailed in the
[Vertex mapping](@ref mapping) tutorial.
