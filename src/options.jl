## Computation options

"""
    BondingMode

Selection mode for the detection of bonds. The choices are:
-   `Input`: use the input bonds. Fail if those are not specified.
-   `Guess`: guess bonds using a variant of chemfiles / VMD algorithm.
-   `Auto`: if the input specifies bonds, use them unless they look suspicious (too or too
    large according to a heuristic). Otherwise, fall back to `Guess`.
"""
module BondingMode
    @enum _BondingMode begin
        Input
        Guess
        Auto
    end
    """See help for `BondingMode`"""
    Input, Guess, Auto
end
import .BondingMode


"""
    ClusteringMode

Selection mode for the clustering of vertices. The choices are:
-   `Input`: use the input residues as clusters. Fail if some atom does
    not belong to a residue.
-   `EachVertex`: each vertex is its own cluster.
-   `MOF`: discard the input residues and consider the input as a MOF. Identify
    organic and inorganic clusters using a simple heuristic based on the atom types.
-   `Auto`: attempt `Input` and fall back to `EachVertex` if the input does not
    provide adequate residues.
-   `Guess`: try to identify the clusters as in `MOF`.
    If it fails, fall back to `Auto`.
"""
module ClusteringMode
    @enum _ClusteringMode begin
        Input
        EachVertex
        MOF
        MOFWiderOrganicSBUs
        MOFMetalloidIsMetal
        Guess
        Auto
    end
    """See help for `ClusteringMode`"""
    Input, EachVertex, MOF, Guess, Auto
    """Internal clustering modes, similar to MOF but with different heuristics"""
    MOFWiderOrganicSBUs, MOFMetalloidIsMetal
end
import .ClusteringMode
import .ClusteringMode: _ClusteringMode


"""
    SBUKinds(sbus, toclassify=Set{Int}())

Description of the different kinds of SBUs there should be when making clusters.
`sbus` should be a list of set of symbols, each set containing the different
elements acceptable in this SBU (an empty set designates all remaining elements).
All elements of the same category of the periodic table can be grouped together
by putting the name of the category. For example
`SBUKinds([[:Au, :halogen, :nonmetal], [:metal, :metalloid], []])` means that
there are three kinds of SBUs:
- the thirst kind can only hold halogens, non-metals and Au atoms
- the second kind can only hold metalloids and metals (except Au)
- the third kind can hold all the other elements.

The list of possible categories is: :actinide, :noble (for noble gas), :halogen,
:lanthanide, :metal, :metalloid and :nonmetal.

`tomerge` contains the list of SBUs which are not actual SBUs but only groups of
atoms waiting to be merged to a neighboring SBU. The neighboring SBU is chosen
by order in the `sbus` list.
For example, `SBUKinds([[:metal], [], [:nonmetal], [:C, :H]], Set{Int}([4]))`
means that all atoms except carbohydrate chains will be grouped into SBUs
corresponding to their nature (SBU of kind 1 for metals, of kind 3 for
non-metals, of kind 2 for the rest), then each carbohydrate chain will be
merged into a neighboring SBU, of kind 1 if any, otherwise of kind 2 if any,
otherwise of kind 3.

Note that, since the entire graph is connected, false SBUs in `tomerge` will
always eventually be merged with another SBU except if the entire graph
comprises only one SBU (in which case there will be an error).

To determine which SBU kind corresponds to a given atom, use `getindex`:
```julia
julia> sbu_kinds = SBUKindsSBUKinds([[:nonmetal, :halogen], [:metal, :F]]);

julia> sbu_kinds[:O]
1

julia> sbu_kinds[:Au]
2

julia> sbu_kinds[:F]
2

julia> sbu_kinds[:Ne] # corresponds to no given SBU kind
0
```
If no empty set has been explicitly added to `sbus` and an element falls outside
of the included categories, the returned SBU kind is 0.

An exception is made for nonmetals which are part of an aromatic carbon cycle: those will
be treated separately and put in the SBU of the corresponding carbons.
"""
struct SBUKinds
    dict::Dict{Symbol,Int}
    default::Int
    tomerge::Set{Int}
    len::Int

    function SBUKinds(sbus, tomerge=Set{Int}())
        dict = Dict{Symbol,Int}()
        default = 0
        for (i, sbu) in enumerate(sbus)
            if isempty(sbu)
                if default != 0
                    throw(ArgumentError("SBUKinds cannot contain multiple empty sets."))
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
        return new(dict, default, Set{Int}(tomerge), length(sbus))
    end
end

function Base.getindex(sbus::SBUKinds, i::Int)
    get(sbus.dict, element_categories[i], sbus.default)
end
function Base.getindex(sbus::SBUKinds, x::Symbol)
    ret = get(sbus.dict, x, 0)
    ret == 0 || return ret
    num = get(atomic_numbers, x, 0)
    num == 0 && throw(MissingAtomInformation("unknown atom name: $x."))
    return sbus[num]
end

Base.length(sbus::SBUKinds) = sbus.len
false_sbus(sbus::SBUKinds) = sbus.tomerge

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
- bonding_mode: one of the [@BondingMode] options, see above.
- cutoff_coeff: coefficient used to detect bonds. Default is 0.75, higher
            values will include bonds that were considered to long before.
- ignore_atoms: set of atom symbols to ignore (for instance [:C,:H] will
            remove carbohydrate solvent residues).
- bond_adjacent_sbus: bond together SBUs which are only separated by a single C atom.
- skip_minimize: assume that the cell is already the unit cell (default is false).
- dimensions: the set of crystal net dimensions to consider. For instance, putting
            Set(3) will ensure that only 3-dimensional nets are considered.
            Default is empty, meaning that all nets are considered.
- ignore_types: disregard atom types to compute the topology, making pcu and
            pcu-b identical for example (default is true)
"""
struct Options
    name::String # used for exports

    # Input options
    bonding_mode::BondingMode._BondingMode
    cutoff_coeff::Float64
    clustering_mode::ClusteringMode._ClusteringMode
    authorize_pruning::Bool
    ignore_atoms::Set{Symbol}
    ignore_homoatomic_bonds::Set{Symbol}
    ignore_homometallic_bonds::Bool
    ignore_low_occupancy::Bool
    export_input::String

    dryrun::Union{Nothing,Dict{Symbol,Union{Nothing,Set{Symbol}}}}

    # Clustering options
    bond_adjacent_sbus::Bool
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
                       clustering_mode=ClusteringMode.Auto,
                       cutoff_coeff=0.75,
                       authorize_pruning=true,
                       ignore_atoms=Set{Symbol}(),
                       ignore_homoatomic_bonds=Set{Symbol}(),
                       ignore_homometallic_bonds=clustering_mode ∈ (ClusteringMode.MOF, ClusteringMode.MOFMetalloidIsMetal, ClusteringMode.MOFWiderOrganicSBUs),
                       ignore_low_occupancy=false,
                       export_input=DOEXPORT[],
                       dryrun=nothing,
                       bond_adjacent_sbus=clustering_mode ∈ (ClusteringMode.MOF, ClusteringMode.MOFMetalloidIsMetal, ClusteringMode.MOFWiderOrganicSBUs),
                       export_attributions="",
                       export_clusters="",
                       skip_minimize=false,
                       dimensions=Set{Int}(),
                       ignore_types=true,
                       export_net=DOEXPORT[],
                       _pos=SVector{3,Float64}[],
                    )

        _export_input = ifbooltempdirorempty(export_input)
        _export_attributions = ifbooltempdirorempty(export_attributions)
        _export_clusters = ifbooltempdirorempty(export_clusters)
        _export_net = ifbooltempdirorempty(export_net)

        new(
            name,
            bonding_mode,
            cutoff_coeff,
            clustering_mode,
            authorize_pruning,
            Set{Symbol}(ignore_atoms),
            Set{Symbol}(ignore_homoatomic_bonds),
            ignore_homometallic_bonds,
            ignore_low_occupancy,
            _export_input,
            dryrun,
            bond_adjacent_sbus,
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
