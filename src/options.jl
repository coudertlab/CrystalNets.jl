## Computation options

"""
    BondingMode

Selection mode for the detection of bonds. The choices are:
- `Input`: use the input bonds. Fail if those are not specified.
- `Guess`: guess bonds using a variant of chemfiles / VMD algorithm.
- `Auto`: if the input specifies bonds, use them unless they look suspicious (too small or
  or too large according to a heuristic). Otherwise, fall back to `Guess`.
"""
module BondingMode
    @enum _BondingMode begin
        Input
        Guess
        Auto
    end
    """See help for [`BondingMode`](@ref)"""
    Input, Guess, Auto
end
import .BondingMode


"""
    StructureType

Selection mode for the clustering of vertices. The choices are:
- `Input`: use the input residues as clusters. Fail if some atom does not belong to a
  residue.
- `EachVertex`: each vertex is its own cluster.
- `Auto`: attempt `Input` and fall back to `EachVertex` if the input does not provide
  adequate residues.
- `MOF`: discard the input residues and consider the input as a MOF. Identify organic and
  inorganic clusters using heuristics based on the atom types.
- `Cluster`: similar to MOF but metallic atoms are not given a wider radius.
- `Zeolite`: Same as `EachVertex` but attempt to enforce that each O atom has exactly two
  neighbours and that they are not O atoms.
- `Guess`: try to identify the clusters as in `Cluster`. If it fails, fall back to `Auto`.
"""
module StructureType
    @enum _StructureType begin
        Input
        EachVertex
        MOF
        MOFWiderOrganicSBUs
        MOFMetalloidIsMetal
        Cluster
        Zeolite
        Guess
        Auto
    end
    """See help for [`StructureType`](@ref)"""
    Input, EachVertex, MOF, Cluster, Zeolite, Guess, Auto
    """Internal clustering modes, similar to MOF but with different heuristics"""
    MOFWiderOrganicSBUs, MOFMetalloidIsMetal
end
import .StructureType
import .StructureType: _StructureType


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
`CrystalNets.ClusterKinds([[:metal, :actinide, :lanthanide], [:C, :halogen],
                           [:P], [:nonmetal, :metalloid], [:noble]], [3, 4])`.
This means that all atoms that are either metals, actinides or lanthanides are assigned to
class 1 and all halogens and C atoms in SBUs of class 2.
Afterwards, each group of adjacent P atoms is assigned either class 1 if any of its
neighbor is of class 1, or class 2 otherwise if any of its neighbor is of class 2.
If no such neighbor exist, it is assigned to class 1.
Finally, each group of adjacent nonmetals and metalloids is assigned class 1 or 2 following
the same rule as for P atoms.

At the end of the procedure, all atoms are thus given a class between `1` and `length(sbus)`
which is not in `toclassify`. See also [`find_sbus`](@ref) for the implementation of this
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
                    throw(ArgumentError("Element \"$x\" appears in at least two SBU kinds ($old and $i)."))
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
    num == 0 && throw(MissingAtomInformation("unknown atom name: $x."))
    return sbus[num]
end

Base.length(sbus::ClusterKinds) = sbus.len

const default_sbus = ClusterKinds([
    [:metal, :actinide, :lanthanide], [:C, :halogen], [:P], [:nonmetal, :metalloid], [:noble]
], [3, 4])



function ifbooltempdirorempty(x)::String
    if x isa Bool
        x ? tempdir() : ""
    else
        x
    end
end

"""
    Options

Different options, passed as keyword arguments:
- name: a name for the structure
- export_input: path to the directory in which to store the .vtf representing
            the parsed structure. Empty string if none.
- export_attributions: path to the directory in which to store the .pdb
            representing the attribution of vertices into SBUs. Empty if none.
- export_clusters: path to the directory in which to store the .vtf representing
            the clustering of vertices. Empty string if none.
- export_net: path to the directory in which to store the .vtf representing the
            extracted net on which the topology is computed. Empty string if none.
- bonding_mode: one of the [`BondingMode`](@ref) options.
- cutoff_coeff: coefficient used to detect bonds. Default is 0.75, higher
            values will include bonds that were considered to long before.
- ignore_atoms: set of atom symbols to ignore (for instance `[:C,:H]` will
            remove carbohydrate solvent residues).
- bond_adjacent_sbus: bond together SBUs which are only separated by a single C atom.
- skip_minimize: assume that the cell is already the unit cell (default is false).
- dimensions: the set of crystal net dimensions to consider. For instance, putting
            `Set(3)` will ensure that only 3-dimensional nets are considered.
            Default is `Set([1,2,3])`.
- ignore_types: disregard atom types to compute the topology, making pcu and
            pcu-b identical for example (default is true)
"""
struct Options
    name::String # used for exports

    # Input options
    bonding_mode::BondingMode._BondingMode
    cutoff_coeff::Float64
    structure::StructureType._StructureType
    authorize_pruning::Bool
    wider_metallic_bonds::Bool
    ignore_atoms::Set{Symbol}
    ignore_homoatomic_bonds::Set{Symbol}
    ignore_homometallic_bonds::Bool
    ignore_low_occupancy::Bool
    export_input::String

    dryrun::Union{Nothing,Dict{Symbol,Union{Nothing,Set{Symbol}}}}

    # Clustering options
    bond_adjacent_sbus::Bool
    cluster_kinds::ClusterKinds
    detect_paddlewheels::Bool
    detect_heterocycles::Bool
    split_O_vertex::Bool
    unify_sbu_decomposition::Bool
    export_attributions::String
    export_clusters::String

    # Topology computation options
    skip_minimize::Bool
    dimensions::Set{Int}
    ignore_types::Bool
    export_net::String

    # Internal
    _pos::Vector{SVector{3,Float64}}

    function Options(; name="unnamed",
                       bonding_mode=BondingMode.Auto,
                       structure=StructureType.Auto,
                       cutoff_coeff=0.75,
                       wider_metallic_bonds=nothing,
                       authorize_pruning=true,
                       ignore_atoms=Set{Symbol}(),
                       ignore_homoatomic_bonds=Set{Symbol}(),
                       ignore_homometallic_bonds=nothing,
                       ignore_low_occupancy=false,
                       export_input=DOEXPORT[],
                       dryrun=nothing,
                       bond_adjacent_sbus=false,
                       cluster_kinds=default_sbus,
                       detect_paddlewheels=true,
                       detect_heterocycles=true,
                       split_O_vertex=true,
                       unify_sbu_decomposition=false,
                       export_attributions="",
                       export_clusters="",
                       skip_minimize=false,
                       dimensions=Set{Int}([1,2,3]),
                       ignore_types=true,
                       export_net=DOEXPORT[],
                       _pos=SVector{3,Float64}[],
                    )

        _export_input = ifbooltempdirorempty(export_input)
        _export_attributions = ifbooltempdirorempty(export_attributions)
        _export_clusters = ifbooltempdirorempty(export_clusters)
        _export_net = ifbooltempdirorempty(export_net)

        _ignore_homometallic_bonds = if ignore_homometallic_bonds === nothing
            structure ∈ (StructureType.MOF, StructureType.MOFMetalloidIsMetal,
                               StructureType.MOFWiderOrganicSBUs)
        else
            ignore_homometallic_bonds
        end

        _wider_metallic_bonds = if wider_metallic_bonds === nothing
            structure ∈ (StructureType.MOF, StructureType.MOFMetalloidIsMetal,
                               StructureType.MOFWiderOrganicSBUs, StructureType.Zeolite)
        else
            wider_metallic_bonds
        end

        new(
            name,
            bonding_mode,
            cutoff_coeff,
            structure,
            _wider_metallic_bonds,
            authorize_pruning,
            Set{Symbol}(ignore_atoms),
            Set{Symbol}(ignore_homoatomic_bonds),
            _ignore_homometallic_bonds,
            ignore_low_occupancy,
            _export_input,
            dryrun,
            bond_adjacent_sbus,
            cluster_kinds,
            detect_paddlewheels,
            detect_heterocycles,
            split_O_vertex,
            unify_sbu_decomposition,
            _export_attributions,
            _export_clusters,
            skip_minimize,
            Set{Int}(dimensions),
            ignore_types,
            _export_net,
            _pos,
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
                val::Union{String,Bool}
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
