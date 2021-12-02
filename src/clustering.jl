using Statistics: mean
import Graphs: SimpleEdge

abstract type ClusteringError end
function Base.showerror(io::IO, e::ClusteringError)
    print(io, "CrystalNets clustering error: ", e.msg)
end

"""
    regroup_sbus(graph::PeriodicGraphs.PeriodicGraph3D, classes::AbstractVector{<:Integer})

Given a classification of vertices into classes, separate the vertices into clusters
of contiguous vertices belonging to the same class.

Class 0 is treated separately: each atom belonging to class 0 is considered to be actually of
class 2, but all its neighbors of classes different than 2 are considered themselves contiguous.
This is specifically used to merge several inorganic SBUs into one if they share a common
neighboring atom.
"""
function regroup_sbus(graph::PeriodicGraph3D, classes::AbstractVector{<:Integer})
    n = length(classes)
    sbus = Vector{PeriodicVertex3D}[]
    sbu_classes = Int[]
    attributions = zeros(Int, n)
    offsets = zeros(SVector{3,Int}, n)
    periodicsbus = Int[0]
    for i in 1:n
        iszero(attributions[i]) || continue
        class = classes[i] == 0 ? 2 : classes[i]
        push!(sbus, [PeriodicVertex3D(i)])
        push!(sbu_classes, class)
        attr = length(sbus)
        attributions[i] = attr
        Q = PeriodicVertex3D[PeriodicVertex3D(i)]
        encounteredclass0 = Set{Int}()
        for u in Q
            @assert attributions[u.v] == attr || (class != 2 && classes[u.v] == 0)
            for neigh in neighbors(graph, u.v)
                x = neigh.v
                ofs = neigh.ofs .+ u.ofs
                if classes[x] == class || (class == 2 && classes[x] == 0)
                    attrx = attributions[x]
                    if iszero(attrx)
                        offsets[x] = ofs
                        attributions[x] = attr
                        y = PeriodicVertex3D(x, ofs)
                        push!(sbus[end], y)
                        push!(Q, y)
                    else
                        @assert attrx == attr
                        if ofs != offsets[x] && last(periodicsbus) != attr
                            push!(periodicsbus, attr)
                        end
                    end
                elseif class != 2 && classes[x] == 0 && x ∉ encounteredclass0
                    push!(encounteredclass0, x)
                    push!(Q, PeriodicVertex3D(x, ofs))
                end
            end
        end
    end

    popfirst!(periodicsbus)
    return Clusters(sbus, sbu_classes, attributions, offsets), Set(periodicsbus)
end



struct MissingAtomInformation <: ClusteringError
    msg::String
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
            classes[i] = 2
        elseif typ ∈ (:O, :Cd, :Co, :Cu, :Ni, :Sc, :W, :Zn, :Zr, :N)
            classes[i] = 1
        else
            throw(MissingAtomInformation("unknown atom type: $typ."))
        end
    end
    return first(regroup_sbus(crystal.graph, classes))
end


function delete_target_from_list!(l, v)
    len = length(l)
    for (j, x) in enumerate(l)
        if x.v == v
            if j != len
                l[j] = l[end]
            end
            resize!(l, len - 1)
            break
        end
    end
    nothing
end

function fix_sbu_connectivity!(sbus, graph, oldperiodic, potentialnewperiodic)
    for i_sbu in potentialnewperiodic
        sbu = sbus.sbus[i_sbu]
        toexplore = [sbu[1]]
        explored = Set{Int}(sbu[1].v)
        while !isempty(toexplore)
            u = pop!(toexplore)
            for x in neighbors(graph, u.v)
                sbus.attributions[x.v] == i_sbu || continue
                ofs = u.ofs .+ x.ofs
                if ofs != sbus.offsets[x.v]
                    @show x.v, x.ofs, u.v, u.ofs
                    throw(InvalidSBU("At least one SBU is periodic itself: cannot coalesce SBUs into new vertices."))
                end
                if x.v ∉ explored
                    push!(explored, x.v)
                    push!(toexplore, PeriodicVertex(x.v, ofs))
                end
            end
        end
        @assert length(explored) == length(sbu)
    end
    periodicsbus = Set{Int}()
    for j_initial in oldperiodic
        sbu = [x.v for x in sbus.sbus[j_initial]]
        empty!(sbus.sbus[j_initial])
        j = j_initial
        toexplore = [PeriodicVertex3D(sbu[1])]
        explored = Set{Int}(sbu[1])
        while true
            seeninthissbu = Dict{Int, SVector{3,Int}}()
            while !isempty(toexplore)
                u = pop!(toexplore)
                push!(sbus.sbus[j], u)
                for x in neighbors(graph, u.v)
                    sbus.attributions[x.v] == j_initial || continue
                    sbus.attributions[x.v] = j
                    if x.v ∈ explored
                        if seeninthissbu[x.v] != x.ofs .+ u.ofs
                            push!(periodicsbus, j)
                        end
                    else
                        ofs = x.ofs .+ u.ofs
                        seeninthissbu[x.v] = ofs
                        sbus.offsets[x.v] = ofs
                        push!(explored, x.v)
                        push!(toexplore, PeriodicVertex(x.v, ofs))
                    end
                end
            end
            if length(explored) < length(sbu)
                push!(sbus.sbus, PeriodicVertex3D[])
                j = length(sbus.sbus)
                push!(sbus.classes, sbus.classes[j_initial])
                push!(toexplore, PeriodicVertex3D(sbu[findfirst(∉(explored), sbu)]))
            end
        end
    end
    return periodicsbus
end

function lazy_avg_dst(pos, sbus, mat)
    n = length(sbus.sbus)
    stored_avg = Vector{SVector{3,Float64}}(undef, n)
    alreadycomputed = falses(n)
    function avg_dst(target, c)
        position = if alreadycomputed[c]
            stored_avg[c]
        else
            p = mean(pos[x.v] .+ x.ofs for x in sbus.sbus[c])
            stored_avg[c] = p .- floor.(p)
        end
        return periodic_distance(position, pos[target], mat)
    end
    return avg_dst
end


const default_sbus = SBUKinds([
    [:metal, :actinide, :lanthanide], [:C, :halogen], [:nonmetal, :metalloid],
], Set{Int}(3))

"""
    find_sbus(crystal)

This is an automatic clustering function for MOFs.
Reckognize SBUs using a simple heuristic based on the atom types.
"""
function find_sbus(crystal, kinds=default_sbus)
    n = nv(crystal.graph)
    classes = Vector{Int}(undef, n)
    unclassified = Int[]
    for i in 1:n
        atom_name = crystal.types[i]
        class = kinds[atom_name]
        if class == 0 # atom not identified
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
                throw(MissingAtomInformation("unhandled atom type: $category (for atom $atom_name)."))
            end
        end
        if class ∈ kinds.tomerge
            push!(unclassified, i)
        end
        classes[i] = class
    end
    #= We now merge `unclassified` elements according to the following rules for
        the special case where kinds == default_sbus
      - All elements of class 3 that are connected to each other are replaced by
        a big virtual element of class 3 whose neighbours are the sum of the
        neighbours of each of its constituents.
    From this point, each element of class 3 only has neighbours of class 1 or 2.
      - If it has a neighbour of class 1, it becomes of class 1 (inorganic SBUs)
      - Otherwise, it only has neighbours of class 2 and then it becomes of
        class 2.
    =#
    @assert issorted(unclassified)
    rev_unclassified = zeros(Int, n)
    m = length(unclassified)
    for i in 1:m
        rev_unclassified[unclassified[i]] = i
    end
    edges = SimpleEdge{Int}[]
    for _i in 1:m
        i = unclassified[_i]
        classi = classes[i]
        for neigh in neighbors(crystal.graph, i)
            j = rev_unclassified[neigh.v]
            (j == 0 || classes[neigh.v] != classi) && continue
            # j is the index of unclassified that corresponds to the neighbour of i
            j >= i && break # We will encounter and handle the symmetric case later
            push!(edges, SimpleEdge(_i, j))
        end
    end
    _graph = SimpleGraph(edges)
    nv(_graph) < m && add_vertices!(_graph, m - nv(_graph))
    # _graph contains the subgraph of unclassified elements with only the bonds
    # between unclassified of the same class
    clusters = connected_components(_graph)

    oldlength = 0
    while oldlength != length(clusters)
        oldlength = length(clusters)
        newclusters = Vector{Int}[]
        for cluster in clusters
            new_class = length(kinds) + 1
            for j in cluster
                for neigh in neighbors(crystal.graph, unclassified[j])
                    k = classes[neigh.v]
                    (k ∈ kinds.tomerge || k > new_class) && continue
                    new_class = k
                    k == 1 && break
                end
                new_class == 1 && break
            end
            if new_class == length(kinds) + 1
                push!(newclusters, cluster)
            else
                for j in cluster
                    @assert classes[unclassified[j]] ∈ kinds.tomerge
                    classes[unclassified[j]] = new_class
                end
            end
        end
        clusters = newclusters
    end
    # oldlength may be != 0, this means that there are dangling unclassified clusters.
    # These will be eliminated as 0-dimensional residues later.
    if crystal.options.cluster_adjacent_sbus && kinds == default_sbus
        # As a last pass, we set a class of 0 for each atom of an organic SBU (class 2)
        # that is bonded to an inorganic SBU (class 1).
        organicclass = kinds[:C]
        inorganicclass = kinds[:metal]
        for i in 1:n
            classes[i] == organicclass || continue
            for x in neighbors(crystal.graph, i)
                if classes[x.v] == inorganicclass
                    classes[i] = 0
                    break
                end
            end
        end
    end

    sbus, periodicsbus = regroup_sbus(crystal.graph, classes)

    while !isempty(periodicsbus)
        # In case of wrong clustering, one or more created SBUs may actually be
        # periodic, making them infinite and meaningless.
        # In this case, we attempt to trim these by reassigning the outermost
        # atoms of these SBUs to other ones until they are no longer periodic.
        # This procedure fails if any other SBU has itself become periodic in
        # the process.
        compositions = [sort([crystal.types[x.v] for x in sbu])
                        for sbu in sbus.sbus]
                encounteredcompositions = Dict{Vector{Symbol},Int}()
        categories = Vector{Int}[]
        category_of = Vector{Int}(undef, length(sbus.sbus))
        for (i, compo) in enumerate(compositions)
            n = length(categories) + 1
            category = get!(encounteredcompositions, compo, n)
            category == n && push!(categories, Int[])
            push!(categories[category], i)
            category_of[i] = category
        end
        # Each element of `categories` is a list of sbu index such that all
        # sbus in the list have the same composition (and are hence considered
        # equivalent)

        reclassifiable = Vector{Int}[]
        reclassifiable_parenttype = Vector{Symbol}[]
        targeted = Set{Int}[]
        # reclassifiable[i] is a set of atoms, each corresponding to an element
        # of a periodic SBU which is the neighbor of an atom X of the i-th SBU.
        # The element type of X is accordingly stored in
        # reclassifiable_parenttype[i] while the involved periodic SBU is
        # stored in targeted[i]
        # If the i-th category is that of a periodic SBU, the i-th element of
        # these three lists is empty
        for (i, sbu) in enumerate(sbus.sbus)
            if i ∈ periodicsbus
                push!(reclassifiable, Int[])
                push!(reclassifiable_parenttype, Symbol[])
                push!(targeted, Set{Int}())
                continue
            end
            _reclassifiable = Int[]
            push!(reclassifiable, _reclassifiable)
            _reclassifiable_parenttype = Symbol[]
            push!(reclassifiable_parenttype, _reclassifiable_parenttype)
            _targeted = Set{Int}()
            push!(targeted, _targeted)
            for x in sbu
                typx = crystal.types[x.v]
                for y in neighbors(crystal.graph, x.v)
                    attr = sbus.attributions[y.v]
                    if attr ∈ periodicsbus
                        push!(_reclassifiable, y.v)
                        push!(_reclassifiable_parenttype, typx)
                        push!(_targeted, attr)
                    end
                end
            end
        end

        # avg_dst(target, c) is the average distance between atom `target` and
        # the atoms in SBU c.
        avg_dst = lazy_avg_dst(crystal.pos, sbus, crystal.cell.mat)

        length_targeted = unique!(sort!(length.(targeted)))
        @assert length_targeted[1] == 0
        popfirst!(length_targeted)
        actually_targeted = Set{Int}()
        modified_sbus = Set{Int}()
        for l in length_targeted
            involvedsbus = Set{Int}() # TODO: use BitSet instead of Set{Int} at multiple places
            _involvedcategories = Set{Int}()
            for i in 1:length(sbus.sbus)
                if length(targeted[i]) == l
                    push!(involvedsbus, i)
                    push!(_involvedcategories, category_of[i])
                end
            end

            involvedcategories::Vector{Int} = collect(_involvedcategories)
            involved_parenttypes = [Set{Symbol}() for _ in involvedcategories]
            for (i, cat) in enumerate(involvedcategories)
                for j in cat
                    j ∈ involvedsbus || continue
                    union!(involved_parenttypes[i], reclassifiable_parenttype[j])
                end
            end

            while !isempty(involvedcategories)
                new_involvedcategories = Int[]
                for (i, i_cat) in enumerate(involvedcategories)
                    thisparenttype = pop!(involved_parenttypes[i])
                    if !isempty(involved_parenttypes)
                        push!(new_involvedcategories, i_cat)
                    end

                    cat = categories[i_cat]
                    toremove = Int[]
                    for i_sbu in cat
                        for (j, parenttype) in enumerate(reclassifiable_parenttype[i_sbu])
                            parenttype === thisparenttype || continue
                            target = reclassifiable[i_sbu][j]
                            push!(toremove, target)
                            if sbus.attributions[target] ∈ periodicsbus # not already dealt with
                                contenders = Dict{Int, Tuple{Int, PeriodicVertex3D}}() # SBUs neighboring this atom
                                for x in neighbors(crystal.graph, target)
                                    attr = sbus.attributions[x.v]
                                    if crystal.types[x.v] === parenttype && attr ∈ cat
                                        count, _ = get(contenders, attr, (0, PeriodicVertex3D(0)))
                                        contenders[attr] = (count + 1, x)
                                    end
                                end
                                μ = maximum(first.(values(contenders)))
                                candidates = [(attr, x) for (attr, (c, x)) in contenders if c == μ]
                                j_sbu, neigh = if length(candidates) > 1
                                    dsts = [avg_dst(target, c[1]) for c in candidates]
                                    candidates[argmin(dsts)]
                                else
                                    only(candidates)
                                end
                                target_sbu = sbus.attributions[target]
                                sbus.attributions[target] = j_sbu
                                push!(actually_targeted, target_sbu)
                                delete_target_from_list!(sbus.sbus[target_sbu], target)
                                newofs = sbus.offsets[neigh.v] .- neigh.ofs
                                sbus.offsets[target] = newofs
                                push!(sbus.sbus[j_sbu], PeriodicVertex(target, newofs))
                                push!(modified_sbus, j_sbu)
                            end
                        end
                    end
                    length(actually_targeted) == length(periodicsbus) && break
                end
                length(actually_targeted) == length(periodicsbus) && break
                involvedcategories = new_involvedcategories
            end

            if length(actually_targeted) == length(periodicsbus)
                # prevperiodicsbus = [Set([x.v for x in sbus.sbus[i]]) for i in periodicsbus]
                periodicsbus = fix_sbu_connectivity!(sbus, crystal.graph, periodicsbus, modified_sbus)
                # errorflag = true
                # for i in 1:length(periodicsbus)
                #     newsbu = Set(periodicsbus[i])
                #     x = first(newsbu).v
                #     xnotfound = true
                #     for j in 1:length(prevperiodicsbus)
                #         if x ∈ prevperiodicsbus[j]
                #             if newsbu ⊆ prevperiodicsbus[j]
                #                 xnotfound = false
                #                 if !issetequal(newsbu, prevperiodicsbus)
                #                     errorflag = false # at least one periodic SBU is smaller
                #                 end
                #             end
                #             # If newsbu ⊈ prevperiodicsbus[j] then xnotfound is
                #             # set, so there will be an error, as should be.
                #             break
                #         end
                #     end
                #     if xnotfound
                #         errorflag = true
                #         break
                #     end
                # end
                break
            end
        end

        # perm = sortperm(compositions; lt = _x -> begin
        #     compo1, compo2 = _x
        #     hasmetal1 = any(x -> ismetal[atomic_numbers[x]], compo1)
        #     hasmetal2 = any(x -> ismetal[atomic_numbers[x]], compo2)
        #     if hasmetal1 != hasmetal2
        #         # If compo1 has metal and not compo2, place compo1 in priority
        #         return hasmetal1
        #     end
        #     # Otherwise, place heterogeneous clusters in priority
        #     if length(compo1) != length(compo2)
        #         return length(compo1) > length(compo2)
        #     end
        #     # Otherwise, fallback to default
        #     return isless(compo1, compo2)
        # end)


        # @show periodicsbus
        # for i in 1:length(periodicsbus)
        #     if !periodicsbus[i]
        #         sbu = getfield.(sbus.sbus[i], :v)
        #         tps = Ref(crystal.types)
        #         @show getindex.(tps, sbu), sbu
        #     end
        # end
        # println()

    end
    return sbus
end


struct InvalidSBU <: ClusteringError
    msg::String
end

"""
    coalesce_sbus(crystal::Crystal, clusters::Clusters)

Return the new crystal corresponding to the input where each cluster has been
transformed into a new vertex.
"""
function coalesce_sbus(crystal::Crystal, mode::_ClusteringMode=crystal.options.clustering_mode, _attempt=1)
    clusters, clustering = find_clusters(crystal, mode)
    if clustering == ClusteringMode.EachVertex
        return Crystal{Nothing}(crystal)
    end
    periodicsbuflag = false
    edgs = PeriodicEdge3D[]
    for e in edges(crystal.graph)
        s, d, of = src(e), dst(e), PeriodicGraphs.ofs(e)
        atts = clusters.attributions[s]
        attd = clusters.attributions[d]
        newofs = of .+ clusters.offsets[s] .- clusters.offsets[d]
        if atts == attd
            iszero(newofs) && continue
            periodicsbuflag = true
            break
        end
        push!(edgs, PeriodicEdge3D(atts, attd, newofs))
    end
    if periodicsbuflag || isempty(edgs)
        if _attempt == 1 && clustering == ClusteringMode.MOF
            return coalesce_sbus(crystal, ClusteringMode.MOFWiderOrganicSBUs, 2)
        end
        if _attempt == 2
           return coalesce_sbus(crystal, ClusteringMode.MOFMetalloidIsMetal, 3)
        end
        if periodicsbuflag
            throw(InvalidSBU("At least one SBU is periodic itself: cannot coalesce SBUs into new vertices."))
        else
            throw(InvalidSBU("Coalescence of SBUs into new nodes leads to an edgeless graph: the clustering is probably wrong and the structure is not connected."))
        end
    end
    n = length(clusters.sbus)
    graph = PeriodicGraph3D(n, edgs)
    if mode == ClusteringMode.Guess && nv(graph) == 1
        return coalesce_sbus(crystal, ClusteringMode.Auto)
    end
    if ne(graph) == 0
        throw(EmptyGraphException())
    end
    if !isempty(crystal.options.export_clusters)
        path = tmpexportname(crystal.options.export_clusters, "clusters_", crystal.options.name, ".pdb")
        export_clusters(Crystal{Clusters}(crystal, clusters), path)
        println("Clustering of vertices represented represented at ", replace(path, ('\\'=>'/')))
    end
    pos = Vector{SVector{3,Float64}}(undef, n)
    types = Vector{Symbol}(undef, n)
    for (i, sbu) in enumerate(clusters.sbus)
        pos[i] = mean(crystal.pos[x.v] .+ x.ofs for x in sbu)
        types[i] = length(sbu) == 1 ? crystal.types[only(sbu).v] : Symbol(clusters.classes[i]) #Symbol(join(sort!([crystal.types[x.v] for x in sbu])))
    end
    return Crystal{Nothing}(crystal.cell, types, pos, graph, crystal.options)
end
