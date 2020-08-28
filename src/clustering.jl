using Statistics: mean
import LightGraphs: SimpleEdge

const elements = Dict{Symbol,String}( # populated using PeriodicTable.jl
    :H => "diatomic nonmetal",
    :He => "noble gas",
    :Li => "alkali metal",
    :Be => "alkaline earth metal",
    :B => "metalloid",
    :C => "polyatomic nonmetal",
    :N => "diatomic nonmetal",
    :O => "diatomic nonmetal",
    :F => "diatomic nonmetal",
    :Ne => "noble gas",
    :Na => "alkali metal",
    :Mg => "alkaline earth metal",
    :Al => "post-transition metal",
    :Si => "metalloid",
    :P => "polyatomic nonmetal",
    :S => "polyatomic nonmetal",
    :Cl => "diatomic nonmetal",
    :Ar => "noble gas",
    :K => "alkali metal",
    :Ca => "alkaline earth metal",
    :Sc => "transition metal",
    :Ti => "transition metal",
    :V => "transition metal",
    :Cr => "transition metal",
    :Mn => "transition metal",
    :Fe => "transition metal",
    :Co => "transition metal",
    :Ni => "transition metal",
    :Cu => "transition metal",
    :Zn => "transition metal",
    :Ga => "post-transition metal",
    :Ge => "metalloid",
    :As => "metalloid",
    :Se => "polyatomic nonmetal",
    :Br => "diatomic nonmetal",
    :Kr => "noble gas",
    :Rb => "alkali metal",
    :Sr => "alkaline earth metal",
    :Y => "transition metal",
    :Zr => "transition metal",
    :Nb => "transition metal",
    :Mo => "transition metal",
    :Tc => "transition metal",
    :Ru => "transition metal",
    :Rh => "transition metal",
    :Pd => "transition metal",
    :Ag => "transition metal",
    :Cd => "transition metal",
    :In => "post-transition metal",
    :Sn => "post-transition metal",
    :Sb => "metalloid",
    :Te => "metalloid",
    :I => "diatomic nonmetal",
    :Xe => "noble gas",
    :Cs => "alkali metal",
    :Ba => "alkaline earth metal",
    :La => "lanthanide",
    :Ce => "lanthanide",
    :Pr => "lanthanide",
    :Nd => "lanthanide",
    :Pm => "lanthanide",
    :Sm => "lanthanide",
    :Eu => "lanthanide",
    :Gd => "lanthanide",
    :Tb => "lanthanide",
    :Dy => "lanthanide",
    :Ho => "lanthanide",
    :Er => "lanthanide",
    :Tm => "lanthanide",
    :Yb => "lanthanide",
    :Lu => "lanthanide",
    :Hf => "transition metal",
    :Ta => "transition metal",
    :W => "transition metal",
    :Re => "transition metal",
    :Os => "transition metal",
    :Ir => "transition metal",
    :Pt => "transition metal",
    :Au => "transition metal",
    :Hg => "transition metal",
    :Tl => "post-transition metal",
    :Pb => "post-transition metal",
    :Bi => "post-transition metal",
    :Po => "post-transition metal",
    :At => "metalloid",
    :Rn => "noble gas",
    :Fr => "alkali metal",
    :Ra => "alkaline earth metal",
    :Ac => "actinide",
    :Th => "actinide",
    :Pa => "actinide",
    :U => "actinide",
    :Np => "actinide",
    :Pu => "actinide",
    :Am => "actinide",
    :Cm => "actinide",
    :Bk => "actinide",
    :Cf => "actinide",
    :Es => "actinide",
    :Fm => "actinide",
    :Md => "actinide",
    :No => "actinide",
    :Lr => "actinide",
    :Rf => "transition metal",
    :Db => "transition metal",
    :Sg => "transition metal",
    :Bh => "transition metal",
    :Hs => "transition metal",
    :Mt => "unknown, probably transition metal",
    :Ds => "unknown, probably transition metal",
    :Rg => "unknown, probably transition metal",
    :Cn => "transition metal",
    :Nh => "unknown, probably transition metal",
    :Fl => "post-transition metal",
    :Mc => "unknown, probably post-transition metal",
    :Lv => "unknown, probably post-transition metal",
    :Ts => "unknown, probably metalloid",
    :Og => "unknown, predicted to be noble gas",
    :Uue => "unknown, but predicted to be an alkali metal"
)


"""
    regroup_sbus(graph::PeriodicGraphs.PeriodicGraph3D, classes::AbstractVector{<:Integer})

Given a classification of vertices into classes, separate the vertices into clusters
of contiguous vertices belonging to the same class.
"""
function regroup_sbus(graph::PeriodicGraph3D, classes::AbstractVector{<:Integer})
    n = length(classes)
    sbus = Vector{PeriodicVertex3D}[]
    sbu_classes = Int[]
    attributions = zeros(Int, n)
    offsets = zeros(SVector{3,Int}, n)
    for i in 1:n
        iszero(attributions[i]) || continue
        class = classes[i]
        push!(sbus, [PeriodicVertex3D(i)])
        push!(sbu_classes, class)
        attr = length(sbus)
        attributions[i] = attr
        Q = PeriodicVertex3D[PeriodicVertex3D(i)]
        for u in Q
            @assert attributions[u.v] == attr
            for neigh in neighbors(graph, u.v)
                x = neigh.v
                if classes[x] == class
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
                end
            end
        end
    end
    return Clusters(sbus, sbu_classes, attributions, offsets)
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
            error("Unknown atom type")
        end
    end
    return regroup_sbus(crystal.graph, classes)
end


struct MissingAtomInformation <: Exception
    msg::String
end


"""
    find_sbus(crystal)

This is an automatic clustering function for MOFs.
Reckognize SBUs using a simple heuristic based on the atom types.
"""
function find_sbus(crystal)
    n = nv(crystal.graph)
    classes = Vector{Int}(undef, n)
    inv_classes = SVector{3,Vector{Int}}(Int[], Int[], Int[])
    for i in 1:n
        atom_name = crystal.types[i]
        if atom_name === :C
            classes[i] = 1 # Class 1 contains organic SBUs
            push!(inv_classes[1], i)
        else
            if !haskey(elements, atom_name)
                if atom_name === Symbol("")
                    throw(MissingAtomInformation("""
                    The input is a periodic graph with no atom information, it cannot be reckognized as a MOF.
                    """))
                elseif atom_name === Symbol("_")
                    throw(MissingAtomInformation("""
                    The input file format does not contain enough information on the atoms to distinguish the organic and inorganic SBUs.
                    Please use a file format containing at least the atom types and in order to use MOF reckognition.
                    """))
                else
                    throw(MissingAtomInformation("Unknown atom name: $atom_name"))
                end
            end
            category = elements[atom_name]
            last_category = last(split(category, ' '))
            if last_category == "metal" || category == "actinide" || category == "lanthanide"
                classes[i] = 2 # Class 2 contains inorganic SBUs
                push!(inv_classes[2], i)
            elseif last_category == "nonmetal" || category == "metalloid"
                classes[i] = 3 # Class 3 is temporary and contains unclassified elements
                push!(inv_classes[3], i)
            else
                throw(MissingAtomInformation("Unhandled atom type: $atom_name"))
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
    unclassified = inv_classes[3]
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
    global clas = classes
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
    pos = Matrix{Float64}(undef, 3, n)
    types = Vector{Symbol}(undef, n)
    for (i, sbu) in enumerate(clusters.sbus)
        pos[:,i] = mean(crystal.pos[:,x.v] .+ x.ofs for x in sbu)
        types[i] = Symbol(clusters.classes[i]) #Symbol(join(sort!([crystal.types[x.v] for x in sbu])))
    end
    edgs = PeriodicEdge3D[]
    for e in edges(crystal.graph)
        s, d, of = src(e), dst(e), PeriodicGraphs.ofs(e)
        atts = clusters.attributions[s]
        attd = clusters.attributions[d]
        atts == attd && continue
        push!(edgs, PeriodicEdge3D(atts, attd, of .+ clusters.offsets[s] .- clusters.offsets[d]))
    end
    return Crystal{Nothing}(crystal.cell, types, nothing, pos, PeriodicGraph3D(n, edgs))
end
