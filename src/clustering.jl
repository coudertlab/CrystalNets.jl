using Statistics: mean
import Graphs: SimpleEdge

const element_categories = String[ # populated using PeriodicTable.jl
    "nonmetal", "gas", "metal", "metal", "metalloid", "nonmetal", "nonmetal",
    "nonmetal", "halogen", "gas", "metal", "metal", "metal", "metalloid",
    "nonmetal", "nonmetal", "halogen", "gas", "metal", "metal", "metal",
    "metal", "metal", "metal", "metal", "metal", "metal", "metal", "metal",
    "metal", "metal", "metalloid", "metalloid", "nonmetal", "halogen", "gas",
    "metal", "metal", "metal", "metal", "metal", "metal", "metal", "metal",
    "metal", "metal", "metal", "metal", "metal", "metal", "metalloid",
    "metalloid", "halogen", "gas", "metal", "metal", "lanthanide",
    "lanthanide", "lanthanide", "lanthanide", "lanthanide", "lanthanide",
    "lanthanide", "lanthanide", "lanthanide", "lanthanide", "lanthanide",
    "lanthanide", "lanthanide", "lanthanide", "lanthanide", "metal", "metal",
    "metal", "metal", "metal", "metal", "metal", "metal", "metal", "metal",
    "metal", "metal", "metal", "halogen", "gas", "metal", "metal", "actinide",
    "actinide", "actinide", "actinide", "actinide", "actinide", "actinide",
    "actinide", "actinide", "actinide", "actinide", "actinide", "actinide",
    "actinide", "actinide", "metal", "metal", "metal", "metal", "metal",
    "metal", "metal", "metal", "metal", "metal", "metal", "metal", "metal",
    "halogen", "gas", "metal"]


"""
    regroup_sbus(graph::PeriodicGraphs.PeriodicGraph3D, classes::AbstractVector{<:Integer})

Given a classification of vertices into classes, separate the vertices into clusters
of contiguous vertices belonging to the same class.
Class 0 is treated separately: each atom belonging to class 0 is considered to be actually of
class 1, but all its neighbors of classes different than 1 are considered themselves contiguous.
This is specifically used to merge several inorganic SBUs into one if they share a common
neighboring atom.
"""
function regroup_sbus(graph::PeriodicGraph3D, classes::AbstractVector{<:Integer})
    n = length(classes)
    sbus = Vector{PeriodicVertex3D}[]
    sbu_classes = Int[]
    attributions = zeros(Int, n)
    offsets = zeros(SVector{3,Int}, n)
    for i in 1:n
        iszero(attributions[i]) || continue
        class = max(1, classes[i])
        push!(sbus, [PeriodicVertex3D(i)])
        push!(sbu_classes, class)
        attr = length(sbus)
        attributions[i] = attr
        Q = PeriodicVertex3D[PeriodicVertex3D(i)]
        for u in Q
            @assert attributions[u.v] == attr || (class != 1 && classes[u.v] == 0)
            for neigh in neighbors(graph, u.v)
                x = neigh.v
                if classes[x] == class || (class == 1 && classes[x] == 0)
                    attrx = attributions[x]
                    if iszero(attrx)
                        ofs = neigh.ofs .+ u.ofs
                        offsets[x] = ofs
                        attributions[x] = attr
                        y = PeriodicVertex3D(x, ofs)
                        push!(sbus[end], y)
                        push!(Q, y)
                    else
                        @assert attrx == attr
                    end
                elseif class != 1 && classes[x] == 0
                    ofs = neigh.ofs .+ u.ofs
                    push!(Q, PeriodicVertex3D(x, ofs))
                end
            end
        end
    end

    return Clusters(sbus, sbu_classes, attributions, offsets)
end


struct MissingAtomInformation <: Exception
    msg::String
end
function Base.showerror(io::IO, e::MissingAtomInformation)
    print(io, "CrystalNets clustering error: ", e.msg)
end


"""
    find_sbus_naive(crystal)

Naive automatic clustering functions for MOFs
Simply assign to an organic clusters carbons and hydrogens, and to inorganic
clusters any atom in (Cd, Co, Cu, N, Ni, O, Sc, W, Zn, Zr). Fail for other atoms.
"""
function find_sbus_naive(crystal)
    n = nv(crystal.graph)
    classes = Vector{Int}(undef, n)
    for i in 1:n
        typ = crystal.types[i]
        if typ === :C || typ === :H
            classes[i] = 1
        elseif typ ∈ (:O, :Cd, :Co, :Cu, :Ni, :Sc, :W, :Zn, :Zr, :N)
            classes[i] = 2
        else
            throw(MissingAtomInformation("unknown atom type: $typ."))
        end
    end
    return regroup_sbus(crystal.graph, classes)
end



"""
    find_sbus(crystal)

This is an automatic clustering function for MOFs.
Reckognize SBUs using a simple heuristic based on the atom types.
"""
function find_sbus(crystal)
    n = nv(crystal.graph)
    classes = Vector{Int}(undef, n)
    unclassified = Int[]
    for i in 1:n
        atom_name = crystal.types[i]
        if atom_name === :C
            classes[i] = 1 # Class 1 contains organic SBUs
        else
            atom = get(atomic_numbers, atom_name, 0)
            if atom == 0
                if atom_name === Symbol("")
                    throw(MissingAtomInformation("""
                    the input is a periodic graph with no atom information, it cannot be reckognized as a MOF.
                    """))
                elseif atom_name === Symbol("_")
                    throw(MissingAtomInformation("""
                    the input file format does not contain enough information on the atoms to distinguish the organic and inorganic SBUs.
                    Please use a file format containing at least the atom types and in order to use MOF reckognition.
                    """))
                else
                    throw(MissingAtomInformation("unknown atom name: $atom_name"))
                end
            end
            category = element_categories[atom]
            if category == "metal" || category == "actinide" || category == "lanthanide"
                classes[i] = 2 # Class 2 contains inorganic SBUs
            elseif category == "nonmetal" || category == "metalloid"
                classes[i] = 3 # Class 3 is temporary and contains unclassified elements
                push!(unclassified, i)
            elseif category == "halogen"
                classes[i] = 1 # Special case for halogens, which should never be both
                # part of inorganic SBUs and more than monovalent.
            else
                throw(MissingAtomInformation("unhandled atom type: $category (for atom $atom_name)"))
            end
        end
    end
    #= We now classify elements of class 3 according to the following rules:
      - All elements of class 3 that are connected to each other are replaced by a big
        virtual element of class 3 whose neighbours are the sum of the neighbours of each
        of its constituents.
        From this point, each element of class 3 only has neighbours of class 1 or 2.
      - If it has a neighbour of class 2, it becomes of class 2 (inorganic SBUs)
      - Otherwise, it only has neighbours of class 1 and then it becomes of class 1.
    =#
    @assert issorted(unclassified)
    rev_unclassified = zeros(Int, n)
    m = length(unclassified)
    for i in 1:m
        rev_unclassified[unclassified[i]] = i
    end
    edges = SimpleEdge{Int}[]
    for i in 1:m
        for neigh in neighbors(crystal.graph, unclassified[i])
            j = rev_unclassified[neigh.v]
            if !iszero(j)
                # j is the index of unclassified that corresponds to the neighbour of i
                j >= i && break # We will encounter and handle the symmetric case later
                push!(edges, SimpleEdge(i,j))
            end
        end
    end
    _graph = SimpleGraph(edges)
    if nv(_graph) < m
        add_vertices!(_graph, m - nv(_graph))
    end
    clusters = connected_components(_graph)
    for cluster in clusters
        new_class = 1
        for j in cluster
            for neigh in neighbors(crystal.graph, unclassified[j])
                k = classes[neigh.v]
                if k == 2
                    new_class = 2
                    break
                end
            end
            new_class == 2 && break
        end
        for j in cluster
            @assert classes[unclassified[j]] == 3
            classes[unclassified[j]] = new_class
        end
    end
    if crystal.options.cluster_adjacent_sbus
        # As a last pass, we set a class of 0 for each atom of an organic SBU (class 1)
        # that is bonded to an inorganic SBU (class 2).
        for i in 1:n
            classes[i] == 1 || continue
            for x in neighbors(crystal.graph, i)
                if classes[x.v] == 2
                    classes[i] = 0
                    break
                end
            end
        end
    end
    return regroup_sbus(crystal.graph, classes)
end


"""
    coalesce_sbus(crystal::Crystal)

Return the new crystal corresponding to the input where each cluster has been
transformed into a new vertex.
"""
function coalesce_sbus(crystal::Crystal{Clusters})
    clusters = crystal.clusters
    n = length(clusters.sbus)
    pos = Vector{SVector{3,Float64}}(undef, n)
    types = Vector{Symbol}(undef, n)
    for (i, sbu) in enumerate(clusters.sbus)
        pos[i] = mean(crystal.pos[x.v] .+ x.ofs for x in sbu)
        types[i] = length(sbu) == 1 ? crystal.types[only(sbu).v] : Symbol(clusters.classes[i]) #Symbol(join(sort!([crystal.types[x.v] for x in sbu])))
    end
    edgs = PeriodicEdge3D[]
    for e in edges(crystal.graph)
        s, d, of = src(e), dst(e), PeriodicGraphs.ofs(e)
        atts = clusters.attributions[s]
        attd = clusters.attributions[d]
        atts == attd && continue
        push!(edgs, PeriodicEdge3D(atts, attd, of .+ clusters.offsets[s] .- clusters.offsets[d]))
    end
    if isempty(edgs)
        throw(MissingAtomInformation("Coalescence of SBUs into new nodes leads to an edgeless graph: the clustering is probably wrong and the structure is not connected."))
    end
    return Crystal{Nothing}(crystal.cell, types, pos, PeriodicGraph3D(n, edgs),
                            crystal.options)
end
