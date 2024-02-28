## Computation options

"""
    Bonding

Selection mode for the detection of bonds. The choices are:
- `Input`: use the input bonds. Fail if those are not specified.
- `Guess`: guess bonds using a variant of chemfiles / VMD algorithm.
- `Auto`: if the input specifies bonds, use them unless they look suspicious (too small or
  or too large according to a heuristic). Otherwise, fall back to `Guess`.
- `NoBond`: do not guess or use any bond. This cannot be used to determine topology.
"""
module Bonding
    @enum _Bonding begin
        Input
        Guess
        Auto
        NoBond
    end
    """See help for [`Bonding`](@ref)"""
    Input, Guess, Auto, NoBond
end
import .Bonding


module StructureType
    @enum _StructureType begin
        Auto
        MOF
        Cluster
        Zeolite
        Guess
    end
    """See help for [`StructureType`](@ref)"""
    Auto, MOF, Cluster, Zeolite, Guess
    export MOF, Cluster, Zeolite
end
import .StructureType
import .StructureType: _StructureType

module Clustering
    @enum _Clustering begin
        Auto         = 1
        Input        = 2
        EachVertex   = 3
        SingleNodes  = 4
        AllNodes     = 5
        Standard     = 6
        PE           = 7
        PEM          = 8
    end
    """See help for [`Clustering`](@ref)"""
    Auto, Input, EachVertex, SingleNodes, AllNodes, Standard, PE, PEM
    export EachVertex, SingleNodes, AllNodes, Standard, PE, PEM
end
import .Clustering
import .Clustering: _Clustering

function clustering_from_num(x::Integer)
    if x == 1
        return Clustering.Auto
    elseif x == 2
        return Clustering.Input
    elseif x == 3
        return Clustering.EachVertex
    elseif x == 4
        return Clustering.SingleNodes
    elseif x == 5
        return Clustering.AllNodes
    elseif x == 6
        return Clustering.Standard
    elseif x == 7
        return Clustering.PE
    elseif x == 8
        return Clustering.PEM
    end
    throw(ArgumentError(lazy"No clustering from number $x"))
end

function clustering_from_symb(x::Symbol)
    if x === :Auto
        return Clustering.Auto
    elseif x === :Input
        return Clustering.Input
    elseif x === :EachVertex
        return Clustering.EachVertex
    elseif x === :SingleNodes
        return Clustering.SingleNodes
    elseif x === :AllNodes
        return Clustering.AllNodes
    elseif x === :Standard
        return Clustering.Standard
    elseif x === :PE
        return Clustering.PE
    elseif x === :PEM
        return Clustering.PEM
    end
    throw(ArgumentError(lazy"No clustering from symbol $x"))
end

function Base.parse(::Type{_Clustering}, s::AbstractString)
    if s == "Auto"
        return Clustering.Auto
    elseif s == "Input"
        return Clustering.Input
    elseif s == "EachVertex"
        return Clustering.EachVertex
    elseif s == "SingleNodes"
        return Clustering.SingleNodes
    elseif s == "AllNodes"
        return Clustering.AllNodes
    elseif s == "Standard"
        return Clustering.Standard
    elseif s == "PE"
        return Clustering.PE
    elseif s == "PEM"
        return Clustering.PEM
    end
    throw(ArgumentError(lazy"No clustering from string \"$s\""))
end

"""
    StructureType

Selection mode for the crystal structure. This choice impacts the bond detection algorithm
as well as the clustering algorithm used.

The choices are:
- `Auto`: No specific structure information. Use Van der Waals radii for bond detection and
  `Input` as [`Clustering`](@ref), or `EachVertex` if the input does not provide residues.
- `MOF`: Use Van der Waals radii for non-metallic atoms and larger radii for metals. Detect
  organic and inorganic clusters and subdivide them according to `AllNodes` and
  `SingleNodes` to identify underlying nets.
- `Cluster`: similar to MOF but metallic atoms are not given a wider radius.
- `Zeolite`: Same as `Auto` but use larger radii for metals (and metalloids) and attempt
  to enforce that each O atom has exactly two neighbours and that they are not O atoms.
- `Guess`: try to identify the clusters as in `Cluster`. If it fails, fall back to `Auto`.
"""
StructureType

"""
    Clustering

The clustering algorithm used to group atoms into vertices.

This choice only affects the creation of a `UnderlyingNets` from a `Crystal`, not the
`Crystal` itself, and in particular not the bond detection algorithm.

The basic choices are:
- `Auto`: determined using the [`StructureType`](@ref).
- `Input`: use the input residues as vertices. Fail if some atom does not belong to a
  residue.
- `EachVertex`: each atom is its own vertex. Vertices with degree 2 or lower are
  iteratively collapsed into edges until all vertices have degree 3 or more.

The next clustering options are designed for MOFs but may target other kinds of frameworks.
In all cases, the clusters are refinements on top of already-defined clusters, such as the
organic and inorganic SBUs defined by the `MOF` structure.
Except for `AllNodes`, infinite clusters (such as the inorganic clusters in a rod MOF) are
split into new finite clusters using heuristics.
- `SingleNodes`: each already-defined cluster is mapped to a vertex.
- `AllNodes`: keep points of extension for organic clusters.
- `Standard`: make each metallic atom its own vertex and do not bond those together if they
  share a common non-metallic neighbour.
- `PE`: stands for Points of Extension. Keep points of extension for organic clusters,
  remove metallic centres and bond their surrounding points of extension.
- `PEM`: stands for Points of Extension and Metals. Keep points of extension for organic
  clusters and each metal centre as a separate vertex.
"""
Clustering

"""
    ClusterKinds(sbus, toclassify=Int[])

Description of the different kinds of SBUs there should be when making clusters.

`sbus` should be a list of set of symbols, each set containing the different elements
acceptable in this SBU (an empty set designates all remaining elements). All elements of
the same category of the periodic table can be grouped together by putting the name of the
category.
For example, `ClusterKinds([[:Au, :halogen, :nonmetal], [:metal, :metalloid], []])` means
that there are three kinds of SBUs:
- the first kind can only hold halogens, non-metals and Au atoms
- the second kind can only hold metalloids and metals (except Au)
- the third kind can hold all the other elements.

The list of possible categories is: :actinide, :noble (for noble gas), :halogen,
:lanthanide, :metal, :metalloid and :nonmetal.

`toclassify` contains the list of SBUs which are not actual SBUs but only groups of
atoms waiting to be merged to a neighboring SBU. The neighboring SBU is chosen
by order in the `sbus` list.

The cluster kinds used by default are
`CrystalNets.ClusterKinds([[:metal, :actinide, :lanthanide], [:C,], [:P, :S],
                           [:nonmetal, :metalloid, :halogen], [:noble]], [3, 4])`.
This means that all atoms that are either metals, actinides or lanthanides are assigned to
class 1 and all C atoms in SBUs of class 2.
Afterwards, each group of adjacent P or S atoms is assigned either class 1 if any of its
neighbor is of class 1, or class 2 otherwise if any of its neighbor is of class 2.
If no such neighbor exist, it is assigned to class 1.
Finally, each group of adjacent nonmetals, metalloids and halogens is assigned class 1 or 2
following the same rule as for P and S atoms.

At the end of the procedure, all atoms are thus given a class between `1` and `length(sbus)`
which is not in `toclassify`. See also [`find_sbus!`](@ref) for the implementation of this
procedure.

To determine which SBU kind corresponds to a given atom, use `getindex`:
```jldoctest
julia> sbu_kinds = CrystalNets.ClusterKinds([[:nonmetal, :halogen], [:metal, :F]]);

julia> sbu_kinds[:O] # nonmetal
1

julia> sbu_kinds[:Au] # metal
2

julia> sbu_kinds[:F] # specifically F
2

julia> sbu_kinds[:Ne] # no given SBU kind
0
```
If no empty set has been explicitly added to `sbus` and an element falls outside of the
included categories, the returned SBU kind is 0.

An exception is made for nonmetals which are part of an aromatic heterocycle: those will be
treated separately and put in the SBU of the corresponding carbons.
"""
struct ClusterKinds
    dict::Dict{Symbol,Int}
    default::Int
    tomerge::Vector{Int}
    len::Int

    function ClusterKinds(sbus, tomerge=Int[])
        dict = Dict{Symbol,Int}()
        default = 0
        for (i, sbu) in enumerate(sbus)
            if isempty(sbu)
                if default != 0
                    throw(ArgumentError("ClusterKinds cannot contain multiple empty sets."))
                end
                default = i
                continue
            end
            for x in sbu
                old = get!(dict, x, i)
                if old != i
                    throw(ArgumentError(lazy"Element \"$x\" appears in at least two SBU kinds ($old and $i)."))
                end
            end
        end
        return new(dict, default, tomerge, length(sbus))
    end
end

function Base.getindex(sbus::ClusterKinds, i::Int)
    get(sbus.dict, element_categories[i], sbus.default)
end
function Base.getindex(sbus::ClusterKinds, x::Symbol)
    ret = get(sbus.dict, x, 0)
    ret == 0 || return ret
    num = get(atomic_numbers, x, 0)
    num == 0 && throw(MissingAtomInformation(lazy"unknown atom name: $x."))
    return sbus[num]
end
function getmetal(sbus::ClusterKinds)
    ret = get(sbus.dict, :metal, nothing)
    ret isa Int && return ret
    ret = get(sbus.dict, :Fe, nothing)
    ret isa Int && return ret
    ret = get(sbus.dict, :Cu, nothing)
    ret isa Int && return ret
    for (sym, at) in atomic_numbers
        ismetal[at] || continue
        ret = get(sbus.dict, sym, nothing)
        ret isa Int && return ret
    end
    return 0
end

Base.length(sbus::ClusterKinds) = sbus.len
Base.:(==)(k1::ClusterKinds, k2::ClusterKinds) = k1.dict == k2.dict && k1.default == k2.default && k1.tomerge == k2.tomerge

const default_sbus = ClusterKinds([
    [:metal, :actinide, :lanthanide], [:C, :Pc, :Ss], [:P, :S],
    [:nonmetal, :metalloid, :halogen], [:noble]
], [3, 4])
# Note: Pc (respectively Ss) is the internal designation of an organic P (resp. S) atom


function ifbooltempdirorempty(x)::String
    if x isa Bool
        x ? tempdir() : ""
    else
        x
    end
end

"""
    Options

Different options, passed as keyword arguments.

## Basic options
- `name`: a name for the structure.
- `bonding`: one of the [`Bonding`](@ref) options. Default is `Bonding.Auto`.
- `structure`: one of the [`StructureType`](@ref) options. Default is `StructureType.Auto`.
- `clusterings`: a list of [`Clustering`](@ref) options. Default is `[Clustering.Auto]`.

## Exports
For each export option, the accepted values are either a string, indicating the path to
the directory in which to store the export, or a boolean, specifying whether or not to do
the export. If the value is `true`, a path will be automatically determined. An empty
string is equivalent to `false`.
- `export_input`: the parsed structure, as a .vtf
- `export_trimmed`: the parsed structure after iteratively removing all atoms having only one
  neighbour, as a .vtf
- `export_attributions`: the attribution of vertices into SBUs, as a .pdb. Only relevant for
  the `MOF` [`StructureType`](@ref).
- `export_clusters`: the clustering of vertices, as a .vtf
- `export_net`: the overall extracted net on which the topology is computed, as a .vtf.
- `export_subnets`: each connected component of the overall net as a separate .vtf file.
  These subnets are defined after grouping vertices according to their [`Clustering`](@ref).

## Other options
- `ignore_atoms`: set of atom symbols to ignore (for instance `[:C,:H]` will
  remove carbohydrate solvent residues).
- `ignore_types`: disregard atom types to compute the topology, making pcu and
  pcu-b identical for example (default is true)
- `cutoff_coeff`: coefficient used to detect bonds. Default is 0.75, higher
  values will include bonds that were considered too long before.
- `skip_minimize`: assume that the cell is already the unit cell (default is false).
- `dimensions`: the set of crystal net dimensions to consider. For instance, putting
  `Set(3)` will ensure that only 3-dimensional nets are considered.
  Default is `Set([1,2,3])`.
- `cluster_kinds`: a [`ClusterKinds`](@ref). Default separates organic and inorganic SBUs.
- `ignore_homoatomic_bonds`: a `Set{Symbol}` such that all X-X bonds of the net are removed
  if X is an atom whose type is in `ignore_homoatomic_bonds`.
- `max_polyhedron_radius`: an integer specifying the maximum number of bonds between two
  corners of the coordination polyhedron built for the [`Clustering.PE`](@ref Clustering) option.
  Default is 4.

## Miscellaneous options
These boolean options have a default value that may be determined by [`Bonding`](@ref),
[`StructureType`](@ref) and [`Clustering`](@ref). They can be directly overriden here.
- `bond_adjacent_sbus`: bond together SBUs which are only separated by a single C atom.
- `authorize_pruning`: remove colliding atoms in the input. Default is true.
- `wider_metallic_bonds`: for bond detections, metals have a radius equal to 1.5× their Van
  der Waals radius. Default is false, unless [`StructureType`](@ref) is `MOF` or `Zeolite`.
- `ignore_homometallic_bonds`: remove all bonds between two metal atoms of the same kind.
- `reduce_homometallic_bonds`: when guessing bonds, do not bond two metallic atoms of the
  same type if they are up to third neighbours anyway.
  Default is false, unless [`StructureType`](@ref) is `MOF`.
- `ignore_metal_cluster_bonds`: do not bond two metallic clusters together if they share at
  least one non-metallic neighbour. Default is false.
- `ignore_low_occupancy`: atoms with occupancy lower than 0.5 are ignored. Default is false.
- `detect_paddlewheels`: detect paddle-wheel pattern and group them into an inorganic vertex.
  Default is true.
- `detect_organiccycles`: detect organic cycles and collapse all belonging C atoms into a new
  vertex. Default is true.
- `detect_pe`: detect organic points-of-extension (organic atoms bonded to another SBU) and
  transform them into vertices. Default is true.
- `cluster_simple_pe`: cluster adjacent points-of-extension if they are not part of a cycle.
  Default is true.
- `separate_metals`: separate each metal atom into its own vertex (instead of grouping them
  to form metallic clusters if they are adjacent or bonded by an oxygen). Default is false,
  unless [`Clustering`](@ref) is `Standard` or `PEM`.
- `premerge_metalbonds`: when a periodic metallic SBU is detected, cluster together bonded
  metal atoms of the same kind before splitting the SBU.
- `split_O_vertex`: if a vertex is composed of a single O, remove it and bond together all of
  its neighbors. Default is true.
- `unify_sbu_decomposition`: apply the same rule to decompose both periodic and finite SBUs.
  Default is false.
- `force_warn`: force printing warning and information even during `..._dataset` function
  calls. Default is false.

## Internal fields
These fields are for internal use and should not be modified by the user:
- `dryrun`: store information on possible options to try (for `guess_topology`).
- `_pos`: the positions of the centre of the clusters collapsed into vertices.
- `error`: store the first error that occured when building the net.
- `throw_error`: if set, throw the error instead of storing it in the `error` field.
- `track_mapping`: if set to `true`, this field will contain, at the end of the computation,
  the mapping from the vertices of the input to the vertices of the returned genome.
  In the case of `determine_topology...` calls, the "input" is the file exported through
  the `export_input` option. Default is false.
"""
struct Options
    name::String # used for exports

    # Input options
    bonding::Bonding._Bonding
    cutoff_coeff::Float64
    structure::StructureType._StructureType
    authorize_pruning::Bool
    wider_metallic_bonds::Bool
    ignore_atoms::Set{Symbol}
    ignore_homoatomic_bonds::Set{Symbol}
    ignore_homometallic_bonds::Bool
    reduce_homometallic_bonds::Bool
    ignore_low_occupancy::Bool
    export_input::String
    export_trimmed::String
    force_warn::Bool

    # Clustering options
    clusterings::Vector{_Clustering}
    bond_adjacent_sbus::Bool
    ignore_metal_cluster_bonds::Union{Nothing,Bool}
    cluster_kinds::ClusterKinds
    detect_paddlewheels::Bool
    detect_organiccycles::Bool
    detect_pe::Bool
    cluster_simple_pe::Bool
    split_O_vertex::Bool
    unify_sbu_decomposition::Bool
    separate_metals::Union{Nothing,Bool}
    premerge_metalbonds::Bool
    max_polyhedron_radius::Int
    export_attributions::String
    export_clusters::String

    # Topology computation options
    skip_minimize::Bool
    dimensions::Set{Int}
    ignore_types::Bool
    export_net::String
    export_subnets::String

    # Internal
    _pos::Vector{SVector{3,Float64}}
    dryrun::Union{Nothing,Dict{Symbol,Union{Nothing,Set{Symbol}}}}
    error::String
    throw_error::Bool
    track_mapping::Union{Nothing,Vector{Int}}

    function Options(; name="unnamed",
                       bonding=Bonding.Auto,
                       structure=StructureType.Auto,
                       cutoff_coeff=0.75,
                       wider_metallic_bonds=nothing,
                       authorize_pruning=true,
                       ignore_atoms=Set{Symbol}(),
                       ignore_homoatomic_bonds=Set{Symbol}(),
                       ignore_homometallic_bonds=false,
                       reduce_homometallic_bonds=nothing,
                       ignore_low_occupancy=false,
                       export_input=DOEXPORT[],
                       export_trimmed=false,
                       force_warn=false,
                       clusterings=_Clustering[Clustering.Auto],
                       bond_adjacent_sbus=false,
                       ignore_metal_cluster_bonds=nothing,
                       cluster_kinds=default_sbus,
                       detect_paddlewheels=true,
                       detect_organiccycles=true,
                       detect_pe=true,
                       cluster_simple_pe=true,
                       split_O_vertex=true,
                       unify_sbu_decomposition=false,
                       separate_metals=nothing,
                       premerge_metalbonds=true,
                       max_polyhedron_radius=4,
                       export_attributions="",
                       export_clusters="",
                       skip_minimize=false,
                       dimensions=Set{Int}([1,2,3]),
                       ignore_types=true,
                       export_net=false,
                       export_subnets=DOEXPORT[],
                       _pos=SVector{3,Float64}[],
                       dryrun=nothing,
                       error="",
                       throw_error=false,
                       track_mapping=nothing,
                    )

        _export_input = ifbooltempdirorempty(export_input)
        _export_trimmed = ifbooltempdirorempty(export_trimmed)
        _export_attributions = ifbooltempdirorempty(export_attributions)
        _export_clusters = ifbooltempdirorempty(export_clusters)
        _export_net = ifbooltempdirorempty(export_net)
        _export_subnets = ifbooltempdirorempty(export_subnets)

        _reduce_homometallic_bonds = if reduce_homometallic_bonds === nothing
            structure == StructureType.MOF
        else
            _reduce_homometallic_bonds
        end

        _wider_metallic_bonds = if wider_metallic_bonds === nothing
            structure == StructureType.MOF || structure == StructureType.Zeolite
        else
            wider_metallic_bonds
        end

        new(
            name,
            bonding,
            cutoff_coeff,
            structure,
            _wider_metallic_bonds,
            authorize_pruning,
            Set{Symbol}(ignore_atoms),
            Set{Symbol}(ignore_homoatomic_bonds),
            ignore_homometallic_bonds,
            _reduce_homometallic_bonds,
            ignore_low_occupancy,
            _export_input,
            _export_trimmed,
            force_warn,
            clusterings,
            bond_adjacent_sbus,
            ignore_metal_cluster_bonds,
            cluster_kinds,
            detect_paddlewheels,
            detect_organiccycles,
            detect_pe,
            cluster_simple_pe,
            split_O_vertex,
            unify_sbu_decomposition,
            separate_metals,
            premerge_metalbonds,
            max_polyhedron_radius,
            _export_attributions,
            _export_clusters,
            skip_minimize,
            Set{Int}(dimensions),
            ignore_types,
            _export_net,
            _export_subnets,
            _pos,
            dryrun,
            error,
            throw_error,
            track_mapping,
        )
    end
end

function Options(options::Options; kwargs...)
    isempty(kwargs) && return options
    base = Dict{Symbol,Any}([x => getfield(options, x) for x in fieldnames(Options)])
    for (kwarg, val) in kwargs
        T = fieldtype(Options, kwarg)
        val = if isconcretetype(T) && !(T <: Enum)
            if T <: Set
                union(base[kwarg], T(val))
            elseif T === String
                if kwarg === :error
                    ifelse(isempty(options.error), val::String, options.error)
                else
                    val::Union{String,Bool}
                end
            else
                T(val)
            end
        else
            val
        end
        base[kwarg] = val
    end
    return Options(; base...)
end

function permute_mapping(options::Options, vmap)
    options.track_mapping isa Vector{Int} || return options
    mapping::Vector{Int} = copy(options.track_mapping)
    if isempty(mapping)
        append!(mapping, vmap)
    else
        for (i, x) in enumerate(mapping)
            x != 0 && (mapping[i] = vmap[x])
        end
    end
    return Options(options; track_mapping=mapping)
end

function rev_permute_mapping(options::Options, rev_vmap, len=length(rev_vmap))
    options.track_mapping isa Vector{Int} || return options
    vmap = zeros(Int, len)
    for (i, x) in enumerate(rev_vmap)
        vmap[x] = i
    end
    return permute_mapping(options, vmap)
end

function (Base.:(==))(opt1::Options, opt2::Options)
    opt1 === opt2 && return true
    for i in 1:nfields(opt1)
        getfield(opt1, i) == getfield(opt2, i) || return false
    end
    return true
end
