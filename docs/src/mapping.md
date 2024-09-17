# [Mapping of vertices](@id mapping)

## Principle

CrystalNets.jl extracts the net underlying the input structure, and can identify its name
when it exists. Sometimes, it is useful to retrieve the mapping between the input structure
and the vertices of the net. This is possible with CrystalNets.jl, using the
`track_mapping` option (see full list of [`Options`](@ref))

There are multiple ways of using it. One consists in preparing an empty list and passing it
to the `track_mapping` option. The list will contain the mapping at the end of the
computation.

In the case of a call to [`determine_topology`](@ref), the net can be extracted from the
[`InterpenetratedTopologyResult`](@ref) returned by the function, while the number of the
input atoms is that corresponding to the "name" field of the .vtf file exported with the
`export_input` option. For example:

```@meta
DocTestSetup = quote
    using CrystalNets
    import CrystalNets: Options, Clustering, Bonding, StructureType
    const PeriodicGraphs = CrystalNets.PeriodicGraphs
    using .PeriodicGraphs
    using PeriodicGraphs.Graphs

    CrystalNets.toggle_export(false)
    CrystalNets.toggle_warning(false)
end
```

```jldoctest PIVCIJ
julia> path_to_PIVCIJ = joinpath(dirname(dirname(pathof(CrystalNets))), "test", "cif", "PIVCIJ.cif");

julia> mapping = Int[]; # pass an empty list to track_mapping

julia> topology = determine_topology(path_to_PIVCIJ; structure=StructureType.MOF, clusterings=[Clustering.SingleNodes], track_mapping=mapping, export_input=true)
Export of input is enabled: saving file at /tmp/input_PIVCIJ.vtf
SingleNodes: fsc

julia> fsc = last(only(first(only(topology))))
fsc

julia> fsc_graph = PeriodicGraph(fsc)
PeriodicGraph3D(2, PeriodicEdge3D[(1, 1, (0,0,1)), (1, 2, (-1,0,0)), (1, 2, (-1,1,0)), (1, 2, (0,0,0)), (1, 2, (0,1,0))])

julia> degree(fsc_graph) # fsc has two vertices, of degree respectively 6 and 4
2-element Vector{Int64}:
 6
 4

julia> first(Iterators.drop(eachline("/tmp/input_PIVCIJ.vtf"), 4))
"atom 0 type Co name 1 resid 0 atomicnumber 27"

julia> mapping[1] # atom 1, Co, is mapped to the 6-coordinated vertex
1

julia> first(Iterators.drop(eachline("/tmp/input_PIVCIJ.vtf"), 11))
"atom 7 type C name 8 resid 0 atomicnumber 6"

julia> mapping[8] # atom 8, C, is part of the 4-coordinated ligand
2

julia> first(Iterators.drop(eachline("/tmp/input_PIVCIJ.vtf"), 13))
"atom 9 type H name 10 resid 0 atomicnumber 1"

julia> mapping[10] # atom 10, H, has been removed early and has no mapping
0
```

```@setup
rm("/tmp/input_PIVCIJ.vtf");
```

!!! note
    The operation that goes from the [`InterpenetratedTopologyResult`](@ref), returned by
    [`determine_topology`](@ref), to the actual genome stored as a `PeriodicGraph`, is
    explained in [a dedicated FAQ question](@ref genomefromname)

When using the lower-level [`topological_genome`](@ref) function, the call will include an
[`Options`](@ref), either given explicitly (for `PeriodicGraph` input) or within the passed
object if it is a [`CrystalNet`](@ref) or an [`UnderlyingNets`](@ref).
In such cases, the mapping can be fetched from the `track_mapping` field of the given
[`Options`](@ref) after the computation, as long as it was initialized with either an empty
list, like in the previous example, or simply the boolean `true`.
For example :

```jldoctest PIVCIJ
julia> crystal_PIVCIJ = parse_chemfile(path_to_PIVCIJ; clusterings=[Clustering.SingleNodes], track_mapping=true);

julia> topological_genome(CrystalNet(crystal_PIVCIJ))
fsc

julia> crystal_PIVCIJ.options.track_mapping == mapping
true
```

## Restriction

!!! warning
    In order to have a one-to-one correspondence between the input structure and the net,
    there must be exactly one net that corresponds to the input. Using `track_mapping` will
    thus fail when there is only one `track_mapping` field for either multiple interpenetrated substructures, or multiple clusterings.

For example, the following command errors:

```jldoctest PIVCIJ
julia> determine_topology(path_to_PIVCIJ; structure=StructureType.MOF, track_mapping=mapping)
ERROR: ArgumentError: Cannot keep a single mapping track for multiple sub-nets. Please use keep_single_track=true only on single components with a single clustering.
[...]
```

The reason is that the `structure=StructureType.MOF` argument makes the default
`Clustering.Auto` value equivalent to both `Clustering.SingleNodes` and
`Clustering.AllNodes`, so the previous call is actually evaluating two different
clusterings, hence the error.

### Workaround: the `keep_single_track=false` option

A quick workaround consists in adding the `keep_single_track=false` option. The result of
the mapping will then be printed after the computation, but it cannot be retrieved
computationally:

```jldoctest PIVCIJ
julia> mapping2 = Int[];

julia> determine_topology(path_to_PIVCIJ; structure=StructureType.MOF, track_mapping=mapping2, keep_single_track=false)
Mapping for CrystalNets.Clustering._Clustering[CrystalNets.Clustering.AllNodes][1, 1, 1, 1, 1, 0, 0, 2, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
Mapping for CrystalNets.Clustering._Clustering[CrystalNets.Clustering.SingleNodes][1, 1, 1, 1, 1, 0, 0, 2, 2, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 2, 2, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
AllNodes: sqc27
SingleNodes: fsc

julia> mapping2 # untouched, cannot be used to retrieve the mapping
Int64[]
```

### Subnet separation

A more proper solution to this issue consists in creating separate [`CrystalNet`](@ref)
instances, one for each substructure and for each target clustering, and thus with their
own separate `track_mapping` field.
The [`UnderlyingNets`](@ref) constructor can be used for this purpose, starting from a
[`CrystalNets.Crystal`](@ref).
For example:

```@meta
DocTestSetup = quote
    using CrystalNets
    import CrystalNets: Options, Clustering, Bonding, StructureType
    const PeriodicGraphs = CrystalNets.PeriodicGraphs
    using .PeriodicGraphs
    using PeriodicGraphs.Graphs

    CrystalNets.toggle_export(false)
    CrystalNets.toggle_warning(false)
end
```

```jldoctest WEBZEK
julia> path_to_WEBZEK = joinpath(dirname(dirname(pathof(CrystalNets))), "test", "cif", "MOFs", "WEBZEK.cif");

julia> crystal_webzek = parse_chemfile(path_to_WEBZEK; structure=StructureType.MOF);

julia> unets_webzek = UnderlyingNets(crystal_webzek);

julia> topological_genome(unets_webzek) # there are 2 interpenetrated subnets, with 2 clusterings each
2 interpenetrated substructures:
⋅ Subnet 1 → AllNodes,SingleNodes: bcu
⋅ Subnet 2 → AllNodes,SingleNodes: dia

julia> webzekA, webzekB = first.(unets_webzek.D3) # only keep 3D nets here
2-element Vector{Vector{CrystalNet3D}}:
 [CrystalNet3D{Rational{Int32}} of WEBZEK_1 with 8 vertices and 32 edges (clustering: AllNodes), CrystalNet3D{Rational{Int32}} of WEBZEK_1 with 8 vertices and 32 edges (clustering: SingleNodes)]
 [CrystalNet3D{Rational{Int32}} of WEBZEK_2 with 4 vertices and 8 edges (clustering: AllNodes), CrystalNet3D{Rational{Int32}} of WEBZEK_2 with 4 vertices and 8 edges (clustering: SingleNodes)]
```

Say we focus on a particular net, for instance the `AllNodes` clustering of substructure B.
After selecting it, its [`Options`](@ref) can be modified by wrapping it in
the [`CrystalNet`](@ref) constructor and using the appropriate keyword arguments:

```jldoctest WEBZEK
julia> subnet = webzekB[1]
CrystalNet3D{Rational{Int32}} of WEBZEK_2 with 4 vertices and 8 edges (clustering: AllNodes)

julia> newmapping = Int[];

julia> newsubnet = CrystalNet(subnet; track_mapping=newmapping);

julia> topological_genome(newsubnet)
dia

julia> newmapping
4-element Vector{Int64}:
 1
 2
 2
 1
```

Note that the call `CrystalNet(subnet; track_mapping=newmapping)` does not modify `subnet`:
it simply creates a copy of `subnet` with a modified `track_mapping` field. As a
consequence, `subnet.options.track_mapping` cannot be used to track the result of the
computation of `topological_genome(newsubnet)`:

```jldoctest WEBZEK
julia> newmapping === newsubnet.options.track_mapping != subnet.options.track_mapping
true
```
