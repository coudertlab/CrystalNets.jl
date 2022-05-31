using Statistics: mean
import Graphs: SimpleEdge

abstract type ClusteringError end
function Base.showerror(io::IO, e::ClusteringError)
    print(io, "CrystalNets clustering error: ", e.msg)
end

"""
    regroup_sbus(graph::PeriodicGraph3D, classes::AbstractVector{<:Integer},
                 isolate=Int[])

Given a classification of vertices into classes, separate the vertices into clusters
of contiguous vertices belonging to the same class.

`isolate` is a list where each atom is separated from the rest of its class. Once all such
atoms of its class are isolated, we look for the connected components of non-isolated atoms
in that class. If such a component has only one neighbours which is an isolated atom, it
is added to the vertex of the isolated atom.
"""
function regroup_sbus(graph::PeriodicGraph3D, classes::AbstractVector{<:Integer},
                      isolate=Int[])
    n = length(classes)
    sbus = Vector{PeriodicVertex3D}[]
    sbu_classes = Int[]
    attributions = zeros(Int, n)
    offsets = zeros(SVector{3,Int}, n)
    periodicsbus = Int[0]
    isolate_map = zeros(Int, n)
    for (i, j) in enumerate(isolate)
        isolate_map[j] = i
    end
    isolate_sbus = [[[PeriodicVertex3D(i)]] for i in isolate]
    isolate_ofs = [[zero(SVector{3,Int})] for _ in isolate]
    for i in 1:n
        iszero(attributions[i]) || continue
        isolate_map[i] == 0 || continue
        class = classes[i]
        thissbu = [PeriodicVertex3D(i)]
        attr = length(sbus) + 1
        attributions[i] = attr
        Q = PeriodicVertex3D[PeriodicVertex3D(i)]
        isolate_neighbors = Set{PeriodicVertex3D}()
        for u in Q
            @toggleassert attributions[u.v] == attr
            for neigh in neighbors(graph, u.v)
                x = neigh.v
                classes[x] == class || continue
                ofs = neigh.ofs .+ u.ofs
                if isolate_map[x] != 0
                    push!(isolate_neighbors, PeriodicVertex3D(isolate_map[x], ofs))
                    continue
                end
                attrx = attributions[x]
                if iszero(attrx)
                    offsets[x] = ofs
                    attributions[x] = attr
                    y = PeriodicVertex3D(x, ofs)
                    push!(thissbu, y)
                    push!(Q, y)
                else
                    @toggleassert attrx == attr
                    if ofs != offsets[x] && last(periodicsbus) != attr
                        push!(periodicsbus, attr)
                    end
                end
            end
        end
        if length(isolate_neighbors) == 1 && last(periodicsbus) != attr
            isolate_neighbor = first(isolate_neighbors)
            push!(isolate_sbus[isolate_neighbor.v], thissbu)
            push!(isolate_ofs[isolate_neighbor.v], .- isolate_neighbor.ofs)
            for x in thissbu
                attributions[x.v] = 0
            end
        else
            push!(sbus, thissbu)
            push!(sbu_classes, class)
        end
    end

    for (i, sbulist, ofslist) in zip(isolate, isolate_sbus, isolate_ofs)
        newsbu = Vector{PeriodicVertex3D}(undef, sum(length, sbulist))
        @toggleassert length(newsbu) ≥ 1
        push!(sbus, newsbu)
        i_newsbu = length(sbus)
        push!(sbu_classes, classes[i])
        counter = 0
        for (mergesbu, ofs) in zip(sbulist, ofslist)
            for x in mergesbu
                newofs = x.ofs .+ ofs
                offsets[x.v] = newofs
                counter += 1
                newsbu[counter] = PeriodicVertex3D(x.v, newofs)
                attributions[x.v] = i_newsbu
                @toggleassert classes[x.v] == classes[i]
            end
        end
    end

    popfirst!(periodicsbus)
    periodic = falses(length(sbus))
    periodic[periodicsbus] .= true
    clusters = Clusters(sbus, sbu_classes, attributions, offsets, periodic)
    return clusters, Set(periodicsbus)
end



struct MissingAtomInformation <: ClusteringError
    msg::String
end

function _trim_monovalent!(graph)
    flag = true
    vmap = collect(1:nv(graph))
    while flag
        flag = false
        toremove = Int[]
        for i in vertices(graph)
            if degree(graph, i) ≤ 1
                push!(toremove, i)
            end
        end
        if !isempty(toremove)
            flag = true
            vmap = vmap[rem_vertices!(graph, toremove)]
        end
    end
    return vmap
end

"""
    trim_monovalent(crystal)

Repeatedly remove monovalent atoms from the crystal until none is left.
"""
function trim_monovalent(crystal::Crystal{T}) where T
    graph = deepcopy(crystal.pge.g)
    vmap = _trim_monovalent!(graph)
    types = crystal.types[vmap]
    pge = PeriodicGraphEmbedding3D(graph, crystal.pge.pos[vmap], crystal.pge.cell)
    if T === Nothing
        return Crystal{Nothing}(pge, types, crystal.options)
    else
        return Crystal{Clusters}(pge, types, crystal.clusters[vmap], crystal.options)
    end
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


function is_paddlewheel_candidate!(memoized, sbus, i, types, periodicsbus)
    memo = memoized[i]
    if memo isa Tuple{Symbol,PeriodicVertex3D}
        return memo
    end
    default_return = (Symbol(""), PeriodicVertex3D(0))
    if i ∈ periodicsbus
        memoized[i] = default_return
        return default_return
    end
    sbu = sbus[i]
    4 ≤ length(sbu) ≤ 6 || return default_return
    typs = [types[x.v] for x in sbu]
    metal = PeriodicVertex3D(0)
    nonmetalcounter = IdDict{Symbol,Int}()
    singleexception = false
    for (j, typ) in enumerate(typs)
        typ === :C && return default_return
        eltyp = element_categories[atomic_numbers[typ]]
        if ((eltyp === :nonmetal) | (eltyp === :halogen))
            nonmetalcounter[typ] = get(nonmetalcounter, typ, 0) + 1
        elseif ismetal[atomic_numbers[typ]]
            if metal.v != 0
                return default_return
            end
            metal = sbu[j]
        else
            singleexception && return default_return
            singleexception = true
        end
    end
    n = length(nonmetalcounter)
    (n == 0 || n ≥ 3) && return default_return
    nonmetal = if n == 2
        singleexception && return default_return
        singleexception = true
        x1, x2 = nonmetalcounter
        x1[2] > x2[2] ? x1[1] : x1[2] == x2[2] ? Symbol("") : x2[1]
    else
        first(nonmetalcounter)[1]
    end
    length(sbu) == 6 && !singleexception && return default_return
    ret = metal.v == 0 ? (default_return) : (nonmetal, metal)
    memoized[i] = ret
    return ret
end

# Do not use, tend to give worse results
function bond_carboxylic_acid!(graph, types)
    for u in vertices(graph)
        if types[u] === :C
            neighs = neighbors(graph, u)
            length(neighs) == 3 || continue
            O1 = PeriodicVertex3D(0); O2 = PeriodicVertex3D(0)
            for x in neighs
                typ = types[x.v]
                if typ === :O
                    if O1.v == 0 || O2.v != 0
                        O1 = x
                    else O2.v == 0
                        O2 = x
                    end
                end
            end
            (O1.v | O2.v) == 0 && continue
            add_edge!(graph, O1.v, PeriodicVertex3D(O2.v, O2.ofs .- O1.ofs))
        end
    end
    nothing
end

"""
    regroup_paddlewheel!(graph, clusters::Clusters, types, periodicsbus)

Identify paddle-wheel patterns made of two opposite SBUs and regroup them into
one.
"""
function regroup_paddlewheel!(graph, clusters::Clusters, types, periodicsbus)
    n = length(clusters.sbus)
    replacements = collect(1:n)
    anypaddlewheel = false
    memoized = Union{Missing,Tuple{Symbol,PeriodicVertex3D}}[missing for _ in 1:n]
    opposite_sbus = zeros(Int, n)
    opposite_ofs = Vector{SVector{3,Int}}(undef, n)
    # alt_sbus = zeros(Int, n)
    # alt_ofs = Vector{SVector{3,Int}}(undef, n)
    metals = Vector{Tuple{PeriodicVertex3D,PeriodicVertex3D}}(undef, n)
    for (i, sbu) in enumerate(clusters.sbus)
        nonmetal, metal = is_paddlewheel_candidate!(memoized, clusters.sbus, i, types, periodicsbus)
        if metal.v == 0
            opposite_sbus[i] = -1
            continue
        end
        thissbu = Set{Int}([y.v for y in sbu])
        opposite_metal = PeriodicVertex3D(0)
        # alt_metal = PeriodicVertex3D(0)
        class_sbu = clusters.classes[i]
        opposite_sbu = 0
        # alt_sbu = 0
        ofs_diff = zero(SVector{3,Int})
        # ofs_diff_alt = zero(SVector{3,Int})
        num_contact = 0
        # num_contact_alt = 0
        for u in sbu
            types[u.v] === nonmetal || nonmetal === Symbol("") || continue
            invalid_paddlewheel = false
            contact = false
            # contact_alt = false
            for x in neighbors(graph, u.v)
                clusters.attributions[x.v] == i && continue
                if types[x.v] !== :C
                    invalid_paddlewheel = true
                    break
                end
                for y in neighbors(graph, x.v)
                    y.v ∈ thissbu && continue
                    types[y.v] !== nonmetal && nonmetal !== Symbol("") && continue
                    attr = clusters.attributions[y.v]
                    class_candidate = clusters.classes[attr]
                    class_candidate == class_sbu || continue
                    if attr == opposite_sbu
                        if u.ofs .+ x.ofs .+ y.ofs .- clusters.offsets[y.v] != ofs_diff
                            invalid_paddlewheel = true
                            break
                        end
                        contact = true
                        continue
                    end
                    # if attr == alt_sbu
                    #     if u.ofs .+ x.ofs .+ y.ofs .- clusters.offsets[y.v] != ofs_diff
                    _, opposite_metal = is_paddlewheel_candidate!(memoized, clusters.sbus, attr, types, periodicsbus)
                    if opposite_metal.v != 0
                        same_metal = types[metal.v] == types[opposite_metal.v]
                        if same_metal
                            if opposite_sbu != 0
                                invalid_paddlewheel = true
                                break
                            end
                            ofs_diff = u.ofs .+ x.ofs .+ y.ofs .- clusters.offsets[y.v]
                            opposite_sbu = attr
                            contact = true
                        end
                    end
                end
                invalid_paddlewheel && break
            end
            num_contact += contact
            if invalid_paddlewheel
                opposite_sbu = 0
                break
            end
        end
        if opposite_sbu == 0 || num_contact ≤ 1
            opposite_sbus[i] = -1
            continue
        end
        if opposite_sbus[opposite_sbu] != i && opposite_sbus[opposite_sbu] != 0
            opposite_sbus[opposite_sbu] = -1
            opposite_sbus[i] = -1
        else
            opposite_sbus[i] = opposite_sbu
            opposite_ofs[i] = ofs_diff
            metals[i] = (metal, opposite_metal)
            anypaddlewheel = true
        end
    end

    if anypaddlewheel
        anypaddlewheel = false
        for (i, (opposite_sbu, ofs_diff, (metal, opposite_metal))) in
                enumerate(zip(opposite_sbus, opposite_ofs, metals))
            opposite_sbu ≤ 0 && continue
            opp = opposite_sbus[opposite_sbu]
            if opp != i
                @toggleassert opp == -1
                continue
            end
            anypaddlewheel = true
            opposite_sbus[opposite_sbu] = 0
            sbu = clusters.sbus[i]
            for u in clusters.sbus[opposite_sbu]
                newofs = u.ofs .+ ofs_diff
                push!(sbu, PeriodicVertex(u.v, newofs))
                clusters.offsets[u.v] = newofs
                clusters.attributions[u.v] = i
            end
            empty!(clusters.sbus[opposite_sbu])
            add_edge!(graph, metal.v, PeriodicVertex(opposite_metal.v, opposite_metal.ofs .+ ofs_diff .- metal.ofs))
            @toggleassert clusters.classes[opposite_sbu] == clusters.classes[i]
            clusters.classes[opposite_sbu] = 0
            replacements[opposite_sbu] = i
        end
    end

    if anypaddlewheel
        k = 0
        removed = Int[]
        for (i, j) in enumerate(replacements)
            replacements[i] -= k
            if j != i
                k += 1
                push!(removed, i)
            end
        end
        for rem in removed
            @toggleassert isempty(clusters.sbus[rem])
            @toggleassert clusters.classes[rem] == 0
        end
        deleteat!(clusters.classes, removed)
        deleteat!(clusters.sbus, removed)
        setremoved = Set{Int}(removed)
        for (i, attr) in enumerate(clusters.attributions)
            @toggleassert attr ∉ setremoved
            clusters.attributions[i] = replacements[attr]
        end
        newperiodicsbus = [replacements[x] for x in periodicsbus]
        empty!(periodicsbus)
        union!(periodicsbus, newperiodicsbus)
    end

    nothing
end


"""
    split_sbu!(sbus, graph, i_sbu, classes)

Split SBU number `i_sbu` into new SBUs according to the updated
`classes`. The first argument `sbus` is modified in-place.
Return the list of newly-created periodic SBUs, if any.
"""
function split_sbu!(sbus, graph, i_sbu, classes)
    sbu = [x.v for x in sbus.sbus[i_sbu]]
    @toggleassert allunique(sbu)
    empty!(sbus.sbus[i_sbu])
    periodicsbus = Set{Int}()
    j = i_sbu
    toexplore = [PeriodicVertex3D(sbu[1])]
    explored = Set{Int}()
    while true
        _u = only(toexplore)
        @toggleassert iszero(_u.ofs)
        class = classes[_u.v]
        push!(explored, _u.v)
        sbus.offsets[_u.v] = _u.ofs
        sbus.attributions[_u.v] = j
        seeninthissbu = Dict{Int, SVector{3,Int}}(_u.v => _u.ofs)
        while !isempty(toexplore)
            u = pop!(toexplore)
            push!(sbus.sbus[j], u)
            for x in neighbors(graph, u.v)
                classes[x.v] == class || continue
                attrx = sbus.attributions[x.v]
                ofs = x.ofs .+ u.ofs
                wasnotexplored = true
                if attrx == i_sbu || attrx == j
                    seen = get(seeninthissbu, x.v, nothing)
                    if seen isa SVector{3,Int}
                        wasnotexplored = false
                        seen == ofs || push!(periodicsbus, j)
                    end
                end
                attrx == i_sbu || continue
                sbus.attributions[x.v] = j
                if wasnotexplored
                    seeninthissbu[x.v] = ofs
                    sbus.offsets[x.v] = ofs
                    push!(explored, x.v)
                    push!(toexplore, PeriodicVertex(x.v, ofs))
                end
            end
        end
        length(explored) == length(sbu) && break
        push!(sbus.sbus, PeriodicVertex3D[])
        j = length(sbus.sbus)
        push!(sbus.classes, sbus.classes[i_sbu])
        push!(toexplore, PeriodicVertex3D(sbu[findfirst(∉(explored), sbu)]))
    end

    return periodicsbus
end

"""
    reclassify!(sbus, newperiodicsbus, newclass, graph, types, classof, i_sbu)

Reclassify the atoms of `sbus.sbus[i_sbu])` according to the following algorithm:
- Let's call "target atom" any atom of type `typ` where `classof[typ] == deg` and either
  `deg == 0` or `deg > 0` and the degree of the atom is `deg`.
- Assign a new SBU for each target atom (one new per atom).
- Look at the connected components of atoms in the SBU which are not target atoms.
  For each connected component that is finite (aperiodic) and has only one neighbor
  which is a target atom, put that component in the same SBU as the neighbor.
"""
function reclassify!(sbus, newperiodicsbus, newclass, graph, types, classof, i_sbu)
    thissbu = Set{Int}(x.v for x in sbus.sbus[i_sbu])
    targets = Dict{Int,Int}() # New SBU number of each target atom
    newsbus = Vector{PeriodicVertex3D}[]
    handled = Set{Int}()
    for v in thissbu
        deg = get(classof, types[v], -1)
        deg == -1 && continue
        if deg == 0 || degree(graph, v) == deg
            push!(newsbus, [PeriodicVertex3D(v)])
            targets[v] = length(newsbus)
            push!(handled, v)
        end
    end
    if isempty(handled)
        push!(newperiodicsbus, i_sbu)
        return false
    end
    isempty(newsbus) && return false
    delete!(thissbu, keys(targets))
    max_inclassof_sbu_counter = length(newsbus)
    periodicsbus = Int[]
    neighbor_configurations = Dict{Set{Pair{Int,SVector{3,Int}}},Int}()
    iterlist = collect(keys(targets))
    push!(iterlist, 0)
    for t_init in iterlist
        neighs = t_init == 0 ? [PeriodicVertex3D(x) for x in thissbu] :
                               neighbors(graph, t_init)
        for u_init in neighs
            v_init = u_init.v
            v_init ∈ thissbu || continue
            v_init ∈ handled && continue
            seeninthissbu = Dict{Int, SVector{3,Int}}(v_init => u_init.ofs)
            Q::Vector{PeriodicVertex3D} = [u_init]
            periodicsbuflag = false
            newsbu_neighbors = Set{Pair{Int,SVector{3,Int}}}()
            for u in Q
                for x in neighbors(graph, u.v)
                    sbus.attributions[x.v] == i_sbu || continue
                    seenofs = get(seeninthissbu, x.v, nothing)
                    if seenofs isa SVector{3,Int}
                        if !periodicsbuflag && u.ofs .+ x.ofs != seenofs
                            periodicsbuflag = true
                        end
                        continue
                    end
                    targetsbu = get(targets, x.v, 0)
                    newofs = u.ofs .+ x.ofs
                    if targetsbu != 0
                        push!(newsbu_neighbors, (targetsbu => .-newofs))
                    else
                        seeninthissbu[x.v] = newofs
                        push!(Q, PeriodicVertex3D(x.v, u.ofs .+ x.ofs))
                    end
                end
            end
            newsbu = [PeriodicVertex3D(x, o) for (x, o) in seeninthissbu]
            union!(handled, keys(seeninthissbu))
            setdiff!(thissbu, keys(seeninthissbu))
            if periodicsbuflag
                push!(periodicsbus, length(sbus.sbus) + length(newsbus) + 1)
            elseif length(newsbu_neighbors) == 1
                targetsbu, ofs = first(newsbu_neighbors)
                append!(newsbus[targetsbu], [PeriodicVertex3D(x.v, x.ofs .+ ofs) for x in newsbu])
                continue
            end
            newsbunumber = length(newsbus)+1
            dest = if periodicsbuflag
                newsbunumber
            else
                if length(newsbu) == 1 && types[first(newsbu).v] === :O
                    # Special case since those special SBUs are handled separately
                    newsbunumber
                else
                    get!(neighbor_configurations, newsbu_neighbors, length(newsbus)+1)
                end
            end
            if dest == length(newsbus)+1
                push!(newsbus, newsbu)
            else
                append!(newsbus[dest], newsbu)
            end
        end
    end

    curr_class = sbus.classes[i_sbu]
    for i in 1:length(newsbus)
        sbu = newsbus[i]
        this_new_class = i > max_inclassof_sbu_counter ? newclass : curr_class
        idx = if i < length(newsbus)
            push!(sbus.sbus, sbu)
            push!(sbus.classes, this_new_class)
            length(sbus.sbus)
        else
            sbus.sbus[i_sbu] = sbu
            sbus.classes[i_sbu] = this_new_class
            i_sbu
        end

        for x in sbu
            sbus.attributions[x.v] = idx
            sbus.offsets[x.v] = x.ofs
        end
    end

    if !isempty(periodicsbus)
        if periodicsbus[end] == length(sbus.sbus) + 1
            periodicsbus[end] = i_sbu
        end
        union!(newperiodicsbus, periodicsbus)
    end
    return length(newsbus) > max_inclassof_sbu_counter
end

"""
    add_to_newclass!(classes, graph, sbus, new_class, v, types, noneighborof)

Set the class of `v` to `new_class`. Then, grow the newly created class by adding connected
components of the SBU of `v` such that the new class does not become periodic and does not
contain any vertex that is a neighbor of a vertex whose type is in `noneighborof`.

If `types === nothing`, disregard the condition on `noneighborof`.
"""
function add_to_newclass!(classes, graph, sbus, new_class, v, types, noneighborof)
    classes[v] = new_class
    current_attribution = sbus.attributions[v]
    global_encountered = Set{Int}()
    init_offset = zero(SVector{3,Int})
    typv = types isa Vector{Symbol} ? types[v] : Symbol("")
    for u_init in neighbors(graph, v)
        sbus.attributions[u_init.v] == current_attribution || continue
        u_init.v ∈ global_encountered && continue
        if types isa Vector{Symbol}
            typ = types[u_init.v]
            typ !== typv && typ ∈ noneighborof && continue
        end
        encountered = Set{Int}([v, u_init.v])
        forbidden = Set{Int}()
        #authorized = Set{Int}()
        periodicflag = false
        Q = [u_init]
        while !isempty(Q)
            u = pop!(Q)
            for x in neighbors(graph, u.v)
                sbus.attributions[x.v] == current_attribution || continue
                if !periodicflag && x.ofs .+ u.ofs != init_offset
                    periodicflag = true
                end
                if x.v == v
                    push!(global_encountered, u.v)
                    #if types isa Vector{Symbol} &&
                    #push!(authorized, u.v)
                else
                    condition = x.v ∉ encountered
                    if condition && types isa Vector{Symbol}
                        typ = types[x.v]
                        typcondition = typ ∈ noneighborof
                        if typcondition
                            push!(forbidden, u.v, x.v)
                            push!(encountered, x.v)
                            condition = false
                        end
                    end
                    if condition
                        push!(encountered, x.v)
                        push!(Q, x)
                    end
                end
            end
        end
        if !periodicflag
            #for x in union!(authorized, setdiff!(encountered, forbidden))
            for x in setdiff!(encountered, forbidden)
                classes[x] = new_class
            end
        end
    end
end

function add_to_merge_or_newclass!(classes, mergeto, graph, sbus, periodicsbus, new_class, v)
    otherattr = -1
    ofs = zero(SVector{3,Int})
    for u in neighbors(graph, v)
        attr = sbus.attributions[u.v]
        if attr ∉ periodicsbus
            if otherattr == -1
                otherattr = attr
                ofs = sbus.offsets[u.v] .- u.ofs
            elseif otherattr != attr || ofs != sbus.offsets[u.v] .- u.ofs
                otherattr = -2
                break
            end
        end
    end
    if otherattr ≥ 0
        mergeto[v] = (ofs, otherattr)
        return false
    end
    add_to_newclass!(classes, graph, sbus, new_class, v, nothing, Set{Symbol}())
    return true
end


function small_cycles_around!(handled, newcycleidx, graph, pos, mat, i, u_init, classes, acceptedclasses, Cclass)
    pos_init = pos[u_init.v] .+ u_init.ofs
    posu = [pos_init]
    init_vec = mat * (pos_init .- pos[i])
    vec = [init_vec]
    init_vec = .- init_vec
    angles = [180.0]
    diffangles = [0.0]
    visited = [PeriodicVertex3D(i)]
    visited_set = Set{PeriodicVertex3D}(visited)
    incycles = Set{Int}()
    parent = [PeriodicVertex3D(i)]
    toexplore = [u_init]
    while !isempty(toexplore)
        u = pop!(toexplore)
        last_visited = last(visited)
        last_parent = pop!(parent)
        while last_parent != last_visited
            delete!(visited_set, pop!(visited))
            last_visited = last(visited)
        end
        if u ∈ visited_set
            for (j, x) in enumerate(visited)
                if x == u
                    viewcycle = @view visited[j:end]
                    # countC = count(x -> classes[x.v] == Cclass, viewcycle)
                    # if 2*countC ≥ length(viewcycle)
                        union!(incycles, x.v for x in viewcycle)
                        last_x = visited[j]
                        l0 = get!(handled, directedge(last_x, visited[end]), Int[])
                        push!(l0, newcycleidx)
                        for k in j+1:length(visited)
                            new_x = visited[k]
                            l = get!(handled, directedge(last_x, new_x), Int[])
                            push!(l, newcycleidx)
                            last_x = new_x
                        end
                    # end
                    break
                end
            end
            continue
        end
        last_posu = pop!(posu)
        last_angle = pop!(angles)
        last_diffangle = pop!(diffangles)
        last_vec = pop!(vec)
        push!(visited, u)
        push!(visited_set, u)
        for x in neighbors(graph, u)
            (x == last_visited || degree(graph, x) == 1) && continue
            classes[x.v] ∈ acceptedclasses || continue
            new_vec = mat * (pos[x.v] .+ x.ofs .- last_posu)
            α = angle(last_vec, .-new_vec)
            90 < α < 145 || continue
            #β = dihedral(last_prev_vec, last_vec, new_vec)
            #β < 20 || β > 160 || continue
            γ = Float64(angle(init_vec, mat * (pos[x.v] .+ x.ofs .- pos_init)))
            if last_parent == PeriodicVertex3D(i)
                last_diffangle = (180 - γ)/2
                last_angle = γ + last_diffangle
            end
            γ < last_angle || continue
            diffangle = last_angle - γ
            abs(diffangle - last_diffangle) < last_diffangle/4 || continue
            push!(toexplore, x)
            push!(parent, u)
            push!(posu, pos[x.v] .+ x.ofs)
            push!(angles, γ)
            push!(diffangles, diffangle)
            push!(vec, new_vec)
        end
    end
    # !isempty(incycles) && @show incycles
    return incycles
end

function detect_organiccycles(classes, graph, pos, mat, Cclass, modifiables)
    handled = Dict{PeriodicEdge3D,Vector{Int}}()
    cycles = Set{Int}[]
    acceptedclasses = union(modifiables, Cclass)
    k = 1
    for e in edges(graph)
        if classes[e.src] == Cclass
            i, u = (e.src, e.dst)
        elseif classes[e.dst.v] == Cclass
            i, u = e.dst.v, PeriodicVertex3D(e.src, .- e.dst.ofs)
        else
            continue
        end
        degree(graph, i) ≤ 1 && continue
        (classes[u.v] ∈ acceptedclasses && degree(graph, u.v) > 1) || continue
        haskey(handled, directedge(e)) && continue
        incycle = small_cycles_around!(handled, k, graph, pos, mat, i, u, classes, acceptedclasses, Cclass)
        if !isempty(incycle)
            push!(cycles, incycle)
            k += 1
        end
    end
    return cycles
end

"""
    group_cycle(organiccycle, types, graph)

Return a list of Vector{PeriodicVertex3D} where each sublist consists in atoms belonging
to the same cycle, and which should thus belong to the same vertex eventually.
"""
function group_cycle(organiccycle, types, graph)
    _same_SBU = Dict{Int,Vector{Int}}() # For each atom, the list of cycles it belongs to
    for (i, cycle) in enumerate(organiccycle)
        for x in cycle
            # types[x] === :C || continue
            push!(get!(_same_SBU, x, Int[]), i)
        end
    end
    m = length(organiccycle)
    union_find_cycles = [i for i in 1:m]
    rev_union_find_cycles = [Int[i] for i in 1:m]
    for v in values(_same_SBU)
        length(v) == 1 && continue
        minsbu = minimum(union_find_cycles[vj] for vj in v)
        for w in v
            k = union_find_cycles[w]
            k == minsbu && continue
            union_find_cycles[w] = minsbu
            for j2 in rev_union_find_cycles[k]
                union_find_cycles[j2] = minsbu
            end
            append!(rev_union_find_cycles[minsbu], rev_union_find_cycles[k])
            empty!(rev_union_find_cycles[k])
        end
    end
    groups = Vector{PeriodicVertex3D}[]
    for cycles in rev_union_find_cycles
        isempty(cycles) && continue
        group_v = Set{Int}()
        for i_cycle in cycles
            union!(group_v, organiccycle[i_cycle])
        end
        i = pop!(group_v)
        Q = [PeriodicVertex3D(i)]
        group = PeriodicVertex3D[]
        while !isempty(Q)
            u = pop!(Q)
            # if types[u.v] === :C
                push!(group, u)
            # end
            for x in neighbors(graph, u.v)
                x.v ∈ group_v || continue
                delete!(group_v, x.v)
                pushfirst!(Q, PeriodicVertex3D(x.v, x.ofs .+ u.ofs))
            end
        end
        push!(groups, group)
    end
    return groups
end

function identify_metallic_type(t::Symbol, kinds::ClusterKinds, metalkind::Int=getmetal(kinds))
    t, at = representative_atom(t)
    ismetal[at] && return true
    k = kinds[t]
    k == metalkind && return true
    k ∈ kinds.tomerge && return missing
    return false
end


function remove_metal_cluster_bonds!(graph, types, opts)
    opts.ignore_metal_cluster_bonds === true || return false
    n = length(types)
    metallic = falses(n)
    metallic_list = Int[]
    kinds = opts.cluster_kinds
    metalkind = getmetal(kinds)
    for (i, t) in enumerate(types)
        _metal_identified = identify_metallic_type(t, kinds, metalkind)
        metal_identified = _metal_identified isa Bool ? _metal_identified::Bool : false
        if metal_identified
            metallic[i] = true
            push!(metallic_list, i)
        end
    end
    nonmetallic_neighbours = Vector{Dict{Int,Set{SVector{3,Int}}}}(undef, n)
    for j in metallic_list
        any(x -> metallic[x.v], neighbors(graph, j)) || continue
        nonmetallic_neighs = Dict{Int,Set{SVector{3,Int}}}()
        for x in neighbors(graph, j)
            if !metallic[x.v]
                push!(get!(nonmetallic_neighs, x.v, Set{SVector{3,Int}}()), x.ofs)
            end
        end
        nonmetallic_neighbours[j] = nonmetallic_neighs
    end
    edges_toremove = PeriodicEdge3D[]
    for e in edges(graph)
        s, (d, o) = e
        if metallic[s] && metallic[d]
            for k in intersect(keys(nonmetallic_neighbours[s]), keys(nonmetallic_neighbours[d]))
                if !isempty(intersect(nonmetallic_neighbours[s][k],
                                      nonmetallic_neighbours[d][k] .+ Ref(o)))
                    push!(edges_toremove, e)
                    break
                end
            end
        end
    end
    for e in edges_toremove
        rem_edge!(graph, e)
    end
    return !isempty(edges_toremove)
end




"""
    find_sbus(crystal, kinds=default_sbus)

Recognize SBUs using heuristics based on the atom types corresponding to the `AllNodes`
clustering algorithm.
"""
function find_sbus(crystal, kinds=default_sbus)
    separate_metals = crystal.options.separate_metals::Bool
    n = nv(crystal.pge.g)
    classes = Vector{Int}(undef, n)
    for i in 1:n
        atom_name = crystal.types[i]
        class = kinds[atom_name]
        if class == 0 # atom not identified
            if atom_name === Symbol("")
                throw(MissingAtomInformation(lazy"""
                the input is a periodic graph with no atom information, it cannot be recognized as a MOF.
                """))
            elseif atom_name === Symbol("_")
                throw(MissingAtomInformation(lazy"""
                the input file format does not contain enough information on the atoms to distinguish the organic and inorganic SBUs.
                Please use a file format containing at least the atom types and in order to use MOF recognition.
                """))
            else
                throw(MissingAtomInformation(lazy"unhandled atom type: $(element_categories[atomic_numbers[atom_name]]) (for atom $atom_name)."))
            end
        end
        classes[i] = class
    end

    temporary_classes = copy(kinds.tomerge)
    temporary_classes_set = Set{Int}(temporary_classes)
    nclasses::Int = length(kinds)
    last_class = nclasses
    aromaticcycleclassmin = -1
    aromaticcycleclassmax = -2
    graph = crystal.pge.g
    organickind = kinds[:C]

    if crystal.options.detect_organiccycles
        aromaticcycleclassmin = last_class + 1
        organiccycle = detect_organiccycles(classes, graph, crystal.pge.pos,
                                    crystal.pge.cell.mat, organickind, temporary_classes_set)
        same_SBU = group_cycle(organiccycle, crystal.types, graph)
        if !isempty(same_SBU)
            graph = PeriodicGraph(crystal.pge.g) # make a modificable copy
            for cycle in same_SBU
                # Each atom in an aromatic cycle is bonded to the other atoms of the
                # same cycle. The newly-formed clique is assigned a new class.
                last_class += 1
                aromaticcycleclass = last_class
                for (i, x) in enumerate(cycle)
                    classes[x.v] = last_class
                    for j in (i+1):length(cycle)
                        y = cycle[j]
                        add_edge!(graph, x.v, PeriodicVertex3D(y.v, y.ofs .- x.ofs))
                    end
                end
            end
        end
        aromaticcycleclassmax = last_class
    end

    rev_nclasses = [0 for _ in 1:last_class]
    for (i, k) in enumerate(temporary_classes)
        rev_nclasses[k] = i
    end
    unclassified = [Int[] for _ in temporary_classes]
    for i in 1:n
        k = rev_nclasses[classes[i]]
        if k != 0
            push!(unclassified[k], i)
        end
    end

    #= We now merge `unclassified` elements according to the following rules
    (exemplified for the special case where kinds == default_sbus)
      - All elements of class 3 (resp. 4) that are connected to each other are replaced by
        a big virtual element of class 3 (resp. 4) whose neighbours are the sum of the
        neighbours of each of its constituents.
      - For each virtual element of class 3, if it has a neighbour of class 1, it becomes
        of class 1 (inorganic SBU). Otherwise, if it has a neighbour of class 2, it becomes
        of class 2. Otherwise (all its neighbours are in class 4), it becomes of class 1.
      - Same for each virtual element of class 4: if it has no neighbour at all or at least
        one of class 1, it becomes of class 1, otherwise it becomes of class 2.
    =#

    for in_tmp_class in unclassified
        tmp_class = popfirst!(temporary_classes)
        # in_tmp_class contains the atoms whose current class is tmp_class
        delete!(temporary_classes_set, tmp_class)
        for x in in_tmp_class
            classes[x] != tmp_class && continue # already handled earlier
            # we now collect the virtual element
            visited = Set{Int}(x)
            current_class = last_class + 1 # the class of the virtual element, will be updated
            Q = Int[x]
            for u in Q
                for x in neighbors(graph, u)
                    v = x.v
                    class = classes[v]
                    if class == tmp_class
                        if v ∉ visited
                            push!(Q, v)
                            push!(visited, v)
                        end
                    elseif class ∉ temporary_classes
                        current_class = min(current_class, ifelse(class > nclasses, organickind, class))
                    end
                end
            end
            current_class == last_class + 1 && (current_class = 1) # no classified neighbour
            for u in visited
                classes[u] = current_class
                if crystal.types[u] === :P && current_class == organickind
                    crystal.types[u] = :Pc # marks an organic P
                elseif crystal.types[u] == :S && current_class == organickind
                    crystal.types[u] = :Ss # marks an organic S
                end
            end
        end
    end

    if crystal.options.detect_pe
        separate_simple_pe = !crystal.options.cluster_simple_pe
        first_pe = true
        initial_last_class = last_class
        for (i, class) in enumerate(classes)
            class == organickind || continue
            for x in neighbors(graph, i)
                classx = classes[x.v]
                if classx ≤ initial_last_class && classx != organickind && !(aromaticcycleclassmin ≤ classx ≤ aromaticcycleclassmax)
                    last_class += first_pe + separate_simple_pe
                    first_pe = false
                    classes[i] = last_class
                    break
                end
            end
        end
    end

    metalkind = getmetal(kinds)
    isolate = Int[]
    if separate_metals
        for (i, t) in enumerate(crystal.types)
            if ismetal[atomic_numbers[t]] || (t === :P && classes[i] == metalkind)
                push!(isolate, i)
            end
        end
    end

    # bond_carboxylic_acid!(graph, crystal.types) # gives worse results
    sbus, periodicsbus = regroup_sbus(graph, classes, isolate)
    if crystal.options.detect_paddlewheels
        regroup_paddlewheel!(graph, sbus, crystal.types, periodicsbus)
    end


    while !isempty(periodicsbus)
        mergeto = Dict{Int,Tuple{SVector{3,Int},Int}}() # List of atoms to merge to the neighboring SBU
        compositions = [[crystal.types[x.v] for x in sbu] for sbu in sbus.sbus]
        unique_compositions = [unique!(sort(x)) for x in compositions]
        elements_for_composition = Vector{Dict{Symbol,Int}}(undef, length(sbus.sbus))
        list_periodicsbus = collect(periodicsbus)
        newperiodicsbus = Set{Int}()
        for i_sbu in list_periodicsbus
            sbu = sbus.sbus[i_sbu]
            uniquecompo = unique_compositions[i_sbu]
            if length(uniquecompo) == 1
                numneighbors = zeros(Int32, length(sbu))
                for (i, x) in enumerate(sbu)
                    for u in neighbors(graph, x.v)
                        numneighbors[i] += sbus.attributions[u.v] != i_sbu
                    end
                end
                _m, _M = extrema(numneighbors)
                if _m != _M && _M != 1
                    for (x, num) in zip(sbu, numneighbors)
                        num == _M || continue
                        last_class += add_to_merge_or_newclass!(classes, mergeto, graph, sbus, periodicsbus, last_class+1, x.v)
                    end
                else # always the same number of neighbors in different SBUs
                    if length(sbu) == 1
                        # Since this SBU is periodic, it consists in a single
                        # atom bonded to one of its replicates. This should not happen.
                        throw(InvalidSBU("Irreducible periodic SBU consisting of a single atom bonded to one of its replicates."))
                    end
                    # Abandon: atomize the SBU.
                    for x in sbu
                        last_class += 1
                        classes[x.v] = last_class
                    end
                end

                for (v, (ofs, attr)) in mergeto
                    sbus.attributions[v] = attr
                    delete_target_from_list!(sbu, v)
                    push!(sbus.sbus[attr], PeriodicVertex(v, ofs))
                    sbus.offsets[v] = ofs
                end
                empty!(mergeto)

                union!(newperiodicsbus, split_sbu!(sbus, graph, i_sbu, classes))

            else # multiple types
                classof::Dict{Symbol,Int} = if isassigned(elements_for_composition, i_sbu)
                    @inbounds elements_for_composition[i_sbu]
                else
                    _ats = [atomic_numbers[k] for k in uniquecompo]
                    metallics = [k for (at, k) in zip(_ats, uniquecompo) if ismetal[at] || k === :P || k === :S]
                    if isempty(metallics)
                        metallics = [k for (at, k) in zip(_ats, uniquecompo) if ismetalloid[at]]
                    end
                    _singulars = if isempty(metallics)
                        # pick the least represented to make a new class
                        hist = Dict{Symbol,Int}([typ => 0 for typ in uniquecompo])
                        for x in sbu
                            hist[crystal.types[x.v]] += 1
                        end
                        minimum_hits = minimum(values(hist))
                        Set{Symbol}([k for (k,v) in hist if v == minimum_hits])
                    else
                        Set{Symbol}(metallics)
                    end
                    #=
                    max_degree = 0
                    max_degree_types = Set{Symbol}()
                    min_degree = 16
                    min_degree_per_type = Dict{Symbol,Int}([t => 16 for t in _singulars])
                    for x in sbu
                        _typ = crystal.types[x.v]
                        if _typ ∈ _singulars
                            deg = degree(graph, x.v)
                            if deg ≥ max_degree
                                if deg == max_degree
                                    push!(max_degree_types, _typ)
                                else
                                    max_degree = deg
                                    empty!(max_degree_types)
                                    push!(max_degree_types, _typ)
                                end
                            end
                            _min_degree = min_degree_per_type[_typ]
                            if deg < _min_degree
                                min_degree_per_type[_typ] = deg
                                min_degree = min(min_degree, deg)
                            end
                        end
                    end
                    _classof::Dict{Symbol,Int} = if min_degree == max_degree
                        Dict{Symbol,Int}([(t => 0) for t in _singulars])
                    else
                        anyeq = any(t -> min_degree_per_type[t] == max_degree,
                                      max_degree_types)
                        __classof = Vector{Pair{Symbol,Int}}(undef, length(_singulars))
                        for (i, t) in enumerate(_singulars)
                            __classof[i] = (t => if t ∈ max_degree_types
                                max_degree == min_degree_per_type[t] ? 0 : anyeq ?
                                    -1 : max_degree
                            else
                                -1
                            end)
                        end
                        Dict{Symbol,Int}(__classof)
                    end=#
                    _classof::Dict{Symbol,Int} = Dict{Symbol,Int}([(t => 0) for t in _singulars])
                    if crystal.options.unify_sbu_decomposition
                        toadd = Int[]
                        for (i, compo) in enumerate(unique_compositions)
                            if compo == uniquecompo
                                elements_for_composition[i] = _classof
                                push!(toadd, i)
                            end
                        end
                        append!(list_periodicsbus, setdiff(toadd, periodicsbus))
                    end
                    # @show _classof
                    _classof
                end

                last_class += reclassify!(sbus, newperiodicsbus, last_class, graph,
                                         crystal.types, classof, i_sbu)

            end
        end

        periodicsbus = newperiodicsbus
    end

    return sbus
end

function _split_this_sbu!(toremove, graph, k, types, stopiftype, sbus)
    if sbus !== nothing
        for x in neighbors(graph, k)
            othersbu = sbus[x.v]
            length(othersbu) == 1 && types[othersbu[1].v] === stopiftype && return nothing
        end
    end
    push!(toremove, k)
    neighs = reverse(neighbors(graph, k))
    n = length(neighs)
    for (i, x) in enumerate(neighs)
        for j in (i+1):n
            y = neighs[j]
            add_edge!(graph, x.v, PeriodicVertex(y.v, y.ofs .- x.ofs))
        end
    end
end

function split_special_sbu!(graph, sbus, subgraph, types, splitO)
    toremove = Int[]
    uniqueCs = Int[]
    for (k, sbu) in enumerate(sbus)
        uniquetypes = :O
        for x in sbu
            typ = types[x.v]
            if (typ === :C) | (typ === :N)
                uniquetypes = :C
            elseif typ !== :O
                uniquetypes = Symbol("")
                break
            end
        end
        uniquetypes === Symbol("") && continue
        if uniquetypes === :O && splitO
            _split_this_sbu!(toremove, graph, k, types, uniquetypes, sbus)
        elseif uniquetypes == :C && length(sbu) == 1
            push!(uniqueCs, k)
        end
    end
    for k in uniqueCs
        # eliminate the case of points of extension:
        any(x -> types[x.v] === :C || types[x.v] === :N, neighbors(subgraph, sbus[k][1].v)) && continue

        neighs = neighbors(graph, k)
        otherwiseconnected = falses(length(neighs))
        for (i, x) in enumerate(neighs)
            done = otherwiseconnected[i]
            for j in (i+1):length(neighs)
                done && otherwiseconnected[j] && continue
                y = neighs[j]
                if has_edge(graph, x.v, PeriodicVertex(y.v, y.ofs .- x.ofs))
                    otherwiseconnected[i] = true
                    otherwiseconnected[j] = true
                    done = true
                end
            end
        end
        if all(otherwiseconnected)
            push!(toremove, k)
        else
            _split_this_sbu!(toremove, graph, k, types, types[only(sbus[k]).v], sbus)
        end
    end
    return rem_vertices!(graph, toremove)
end

function split_O_vertices(c)
    toremove = Int[]
    graph = deepcopy(c.pge.g)
    for (k, t) in enumerate(c.types)
        (t === :O && degree(graph, k) > 2) || continue
        _split_this_sbu!(toremove, graph, k, c.types, :O, nothing)
    end
    vmap = rem_vertices!(graph, toremove)
    pge = PeriodicGraphEmbedding3D(graph, c.pge.pos[vmap], c.pge.cell)
    return Crystal{Nothing}(pge, c.types[vmap], Options(c.options; _pos=pge.pos))
end


struct InvalidSBU <: ClusteringError
    msg::String
end


function identify_clustering(c::Crystal{T}, structure::_StructureType, clustering::_Clustering) where T
    if clustering == Clustering.Auto
        if structure == StructureType.Auto
            if T === Clusters
                return structure, clustering, c.clusters
            else
                return identify_clustering(c, structure, Clustering.EachVertex)
            end
        elseif structure == StructureType.Zeolite
            return identify_clustering(c, structure, Clustering.EachVertex)
        elseif structure == StructureType.Guess
            separate_metals = c.options.separate_metals isa Bool ? c.options.separate_metals::Bool : false
            return structure, clustering, separate_metals
        end
        @toggleassert structure == StructureType.MOF || structure == StructureType.Cluster
        return identify_clustering(c, structure, Clustering.AllNodes)
    elseif clustering == Clustering.EachVertex
        return structure, clustering, Clusters(length(c.types))
    elseif clustering == Clustering.Input
        if T === Clusters
            return structure, clustering, c.clusters
        else
            throw(ArgumentError("Cannot use input residues as clusters: the input does not have residues."))
        end
    end
    if clustering == Clustering.SingleNodes
        clustering = Clustering.AllNodes
    elseif clustering == Clustering.PE || clustering == Clustering.Standard
        clustering = Clustering.PEM
    end
    separate_metals2 = c.options.separate_metals isa Bool ? c.options.separate_metals::Bool : clustering == Clustering.PEM
    return structure, clustering, separate_metals2
end


function order_atomtype(sym)
    sym === :C && return 119 # C is put first to identify organic clusters
    sym === :P && return 1   # Put S and O before the other non-metals, but after
    sym === :S && return 0   # everything else
    anum = atomic_numbers[sym]
    elcat = element_categories[anum]
    if elcat === :nonmetal || elcat === :noble || elcat === :halogen
        return anum - 120 # Put non-metals at last
    end
    return anum
end

function _collapse_clusters(crystal::Crystal{Nothing}, clusters::Clusters, onlynv::Bool=false, bypassexport=false)
    structure = crystal.options.structure
    clustering = only(crystal.options.clusterings)
    if clustering == Clustering.EachVertex || structure == StructureType.Zeolite
        if crystal.options.split_O_vertex
            return split_O_vertices(crystal)
        end
        return onlynv ? crystal : Crystal{Nothing}(crystal; _pos=crystal.pge.pos)
    end
    edgs = PeriodicEdge3D[]
    for s in vertices(crystal.pge.g)
        atts = clusters.attributions[s]
        neigh0 = neighbors(crystal.pge.g, s)
        @toggleassert length(neigh0) ≥ 2
        typisnotC = crystal.types[s] !== :C
        for x in neigh0
            d = x.v
            if crystal.options.bond_adjacent_sbus && crystal.types[d] === :C && typisnotC
                neighs = neighbors(crystal.pge.g, d)
                if length(neighs) > 2
                    for y in neighs
                        atty = clusters.attributions[y.v]
                        crystal.types[y.v] === :C && continue
                        newofs = x.ofs .+ y.ofs .+ clusters.offsets[s] .- clusters.offsets[y.v]
                        atts == atty && iszero(newofs) && continue
                        push!(edgs, PeriodicEdge3D(atts, atty, newofs))
                    end
                end
            end
            if d > s || (d == s && x.ofs > zero(SVector{3,Int}))
                attd = clusters.attributions[d]
                newofs = x.ofs .+ clusters.offsets[s] .- clusters.offsets[d]
                if atts == attd && d != s
                    @toggleassert iszero(newofs)
                    continue
                end
                push!(edgs, PeriodicEdge3D(atts, attd, newofs))
            end
        end
    end
    n = length(clusters.sbus)
    graph = PeriodicGraph3D(n, edgs)
    if ne(graph) == 0
        if !isempty(crystal.options.export_trimmed)
            @ifwarn @warn "Could not export trimmed crystal: the resulting structure is empty."
        end
        return Crystal{Nothing}(PeriodicGraphEmbedding{3,Float64}(crystal.pge.cell), Symbol[],
                                 Options(crystal.options; _pos=SVector{3,Float64}[]))
    end
    if !onlynv && !bypassexport
        export_default(crystal, "trimmed", crystal.options.name, crystal.options.export_trimmed)
        if !isempty(crystal.options.export_attributions)
            path = tmpexportname(crystal.options.export_attributions, "attribution_", crystal.options.name, ".pdb")
            export_attributions(Crystal{Clusters}(crystal, clusters), path)
            println("Attributions of atoms into SBUs represented represented at ", replace(path, ('\\'=>'/')))
        end
    end

    sbus = clusters.sbus[split_special_sbu!(graph, clusters.sbus, crystal.pge.g, crystal.types, crystal.options.split_O_vertex)]
    n = length(sbus)
    pos = Vector{SVector{3,Float64}}(undef, n)
    types = Vector{Symbol}(undef, n)
    onlynv && return Crystal{Nothing}(crystal.pge.cell, types, pos, graph, crystal.options)
    for (i, sbu) in enumerate(sbus)
        pos[i] = mean(crystal.pge.pos[x.v] .+ x.ofs for x in sbu)
        name = sort!([crystal.types[x.v] for x in sbu])
        push!(name, Symbol(""))
        newname = Tuple{Int,String}[]
        counter = 1
        for j in 2:length(name)
            if name[j] == name[j-1]
                counter += 1
            else
                sym = name[j-1]
                str_sym = counter == 1 ? string(sym) : string(sym, counter)
                anum = order_atomtype(sym)
                push!(newname, (anum, str_sym))
                counter = 1
            end
        end
        sort!(newname; by=first, rev=true)
        types[i] = length(sbu) == 1 ? crystal.types[first(sbu).v] : Symbol(join(last.(newname)))
    end
    ret = Crystal{Nothing}(crystal.pge.cell, types, pos, graph, Options(crystal.options; _pos=pos))
    return ret
end

function _find_clusters(c::Crystal{T}, guess::Bool, separate_metals::Bool)::Tuple{Bool,Clusters} where T
    if c.options.separate_metals === nothing
        c = Crystal{T}(c; separate_metals)
    end
    if guess
        clusters::Clusters = try
            find_sbus(c, c.options.cluster_kinds)
        catch e
            if !(e isa ClusteringError)
                rethrow()
            end
            return true, T === Clusters ? c.clusters : Clusters(length(c.types))
        end
        newc = Crystal{Nothing}(c; clusterings=[Clustering.Auto], export_attributions=false, export_input=false)
        if nv(_collapse_clusters(newc, clusters, true).pge.g) <= 1
            return true, T === Clusters ? c.clusters : Clusters(length(c.types))
        end
        return false, clusters
    end
    return false, find_sbus(c, c.options.cluster_kinds)
end

function find_clusters(_c::Crystal{T}) where T
    c = trim_monovalent(_c)
    _structure = c.options.structure
    clusterings = c.options.clusterings
    ret = Vector{Tuple{Crystal{Nothing},Union{Int,Clusters}}}(undef, length(clusterings))
    encountered = Dict{Tuple{_StructureType,_Clustering,Bool},Int}()
    for (i, _clustering) in enumerate(clusterings)
        identified = identify_clustering(c, _structure, _clustering)
        structure = identified[1]
        maybeclusters = identified[3]
        idx = maybeclusters isa Bool ? identified::Tuple{_StructureType,_Clustering,Bool} :
                                       (structure, identified[2], true)
        clusters::Union{Int,Clusters} = get!(encountered, idx, i)
        if clusters == i
            if maybeclusters isa Clusters
                clusters = maybeclusters
            else
                separate_metals = maybeclusters
                guess = structure == StructureType.Guess && identified[2] == Clustering.Auto
                _guess, clusters = _find_clusters(c, guess, separate_metals)
                if guess && !_guess
                    structure = StructureType.Cluster
                end
            end
        end
        ret[i] = (Crystal{Nothing}(c; structure, clusterings=[_clustering]), clusters)
    end
    return ret
end


"""
    collapse_clusters(crystal::Crystal)

Return the list of crystals corresponding to the input where each cluster has been
transformed into a new vertex, for each targeted clustering.
"""
function collapse_clusters(crystal::Crystal)
    crystalclusters::Vector{Tuple{Crystal{Nothing},Union{Int,Clusters}}} = find_clusters(crystal)
    ret = Vector{Crystal{Nothing}}(undef, length(crystalclusters))
    alreadyexported = false
    for (i, (cryst, clust)) in enumerate(crystalclusters)
        if clust isa Int
            ret[i] = Crystal{Nothing}(ret[clust]; clusterings=cryst.options.clusterings)
        else
            ret[i] = _collapse_clusters(cryst, clust, false, alreadyexported)
            alreadyexported = true
        end
    end
    return ret
end



function update_new_edgs!(new_edgs, atoms, keep_edges, toremove, rev_vmap)
    n = length(atoms)
    for i in 1:n
        x = atoms[i]
        toremove[x.v] && continue
        for j in (i+1):n
            y = atoms[j]
            toremove[y.v] && continue
            val = keep_edges[i,j]
            x.v == y.v && y.ofs == x.ofs && !val && continue
            e = directedge(PeriodicEdge3D(rev_vmap[x.v], rev_vmap[y.v], y.ofs .- x.ofs))
            if get!(new_edgs, e, val) && !val
                new_edgs[e] = false
            end
        end
    end
    nothing
end

function edges_of_convex_hull(atoms, num_targets, ref_mat, pos, toremove, visited_coplanar4,
                              previous_keep_edges=trues(0,0))
    n = length(atoms)
    _m = size(previous_keep_edges)[1]
    m = min(num_targets, _m)
    m == n && return previous_keep_edges
    keep_edges = trues(n, n)
    if m == _m
        keep_edges[1:m,1:m] .= previous_keep_edges
    else
        keep_edges[1:m,1:m] .= (@view previous_keep_edges[1:m,1:m])
    end
    for i in (m+1):n
        keep_edges[i,i] = false
    end
    n ≤ 3 && return keep_edges
    poss = [pos[x.v] .+ x.ofs for x in atoms]

    for i in 1:num_targets
        ref = poss[i]
        for j1 in 1:n
            keep_edges[i,j1] || continue
            vec1 = poss[j1] .- ref
            j2array = [j for j in 1:n if keep_edges[i,j] && j != j1]
            for j2 in j2array
                vec2 = poss[j2] .- ref
                j3array = [j for j in j2array if keep_edges[i,j] && j != j2]
                for j3 in j3array
                    ordered4 = sort(SVector{4,Int}((i, j1, j2, j3)))
                    ordered4 ∈ visited_coplanar4 && continue
                    vec3 = poss[j3] .- ref
                    β = dihedral(vec1, vec2, vec3)
                    β < 7 || β > 173 || continue
                    push!(visited_coplanar4, ordered4)
                    any(x -> toremove[atoms[x].v], (i, j1, j2, j3)) && continue
                    a, b, c, d = [i,j1,j2,j3][sortperm(poss[[i,j1,j2,j3]])]
                    if poss[d][1] ≤ 0.9*poss[a][1] + 0.1*poss[b][1]
                        a, b, c, d = [i,j1,j2,j3][sortperm(poss[[i,j1,j2,j3]]; by=x->x[2])]
                    end
                    flagloop = 1
                    keepac = false
                    keepbd = false
                    keepsmallest = false
                    while flagloop > 0
                        flagloop = - flagloop
                        ba = poss[a] .- poss[b]
                        bc = poss[c] .- poss[b]
                        bd = poss[d] .- poss[b]
                        αc = angle(ba, bc)
                        αd = angle(ba, bd)
                        if αd > αc
                            c, d = d, c
                            bc, bd = bd, bc
                        end
                        cd = poss[d] .- poss[c]
                        cb = .- bc
                        γmid = angle(cb, cd)
                        ca = poss[a] .- poss[c]
                        γ1 = angle(cb, ca)
                        γ2 = angle(ca, cd)
                        keepsmallest |= ((γmid < γ1) | (γmid < γ2))
                        if γmid < γ1 - 30 || γmid < γ2 - 30
                            if αd > αc
                                c, d = d, c
                            end
                            if flagloop == -6
                                a, b, c, d = d, a, b, c
                                flagloop = 0
                            elseif flagloop == -4
                                b, c, d = b, d, c
                            elseif mod(-flagloop, 2) == 1
                                b, c, d = d, c, b
                            else
                                b, c, d = d, b, c
                            end
                            flagloop = -flagloop + 1
                            continue
                        end
                        Δac = norm(ref_mat * ca)
                        Δbd = norm(ref_mat * bd)
                        # Δab = norm(ref_mat * ba)
                        # Δbc = norm(ref_mat * bc)
                        # Δcd = norm(ref_mat * cd)
                        # Δda = norm(ref_mat * (poss[d] .- poss[a]))
                        # minΔ = min(Δab, Δbc, Δcd, Δda)*1.01
                        # if keepsmallest
                            if Δac < 0.8Δbd
                                keepac = true
                            elseif Δbd < 0.8Δac
                                keepbd = true
                            end
                        # else
                    end
                    keepac || (keep_edges[a, c] = keep_edges[c, a] = false)
                    keepbd || (keep_edges[b, d] = keep_edges[d, b] = false)
                    n == 4 && return keep_edges
                end
            end
        end
    end

    n ≥ 5 || return keep_edges
    for i in 1:num_targets
        ref = poss[i]
        for k in (i+1):n
            keep_edges[i,k] || continue
            toremove[atoms[k].v] && continue
            vec = poss[k] .- ref
            j1array = [j for j in 1:n if keep_edges[i,j] && j != k]
            breakflag = false
            for j1 in j1array
                oneremovedflag1 = toremove[atoms[j1].v]
                vec1 = poss[j1] .- ref
                j2array = [j for j in j1array if keep_edges[i,j] && j != j1 &&
                            sort(SVector{4,Int}(i, j1, j, k)) ∉ visited_coplanar4]
                for j2 in j2array
                    oneremovedflag = oneremovedflag1
                    if toremove[atoms[j2].v]
                        oneremovedflag && continue
                        oneremovedflag = true
                    end
                    vec2 = poss[j2] .- ref
                    for j3 in j2array
                        j3 == j2 && continue
                        keep_edges[i,j3] || continue
                        oneremovedflag && toremove[atoms[j3].v] && continue
                        sort(SVector{4,Int}((i, j1, j2, j3))) ∈ visited_coplanar4 && continue
                        sort(SVector{4,Int}((i, j1, j3, k))) ∈ visited_coplanar4 && continue
                        sort(SVector{4,Int}((i, j2, j3, k))) ∈ visited_coplanar4 && continue
                        mat = inv(hcat(vec1, vec2, poss[j3] .- ref))
                        coeffs = mat * vec
                        breakflag = all(≥(-3e-16), coeffs)
                        keep_edges[i,k] = keep_edges[k,i] = !breakflag
                        breakflag && break
                    end
                    breakflag && break
                end
                breakflag && break
            end
        end
    end
    return keep_edges
end


function regroup_toremove(cryst, tomerge, toremove_list, msg)
    if isempty(toremove_list)
        ret = trimmed_crystal(Crystal{Nothing}(cryst.pge, cryst.types, Options(cryst.options; _pos=cryst.pge.pos)))
        export_default(ret, lazy"clusters_$msg", cryst.options.name, cryst.options.export_clusters)
        return ret
    end
    toremove_all = reduce(vcat, toremove_list)
    previous_graph = PeriodicGraph(cryst.pge.g)
    for v in tomerge
        neighs = neighbors(previous_graph, v)
        for (i, x) in enumerate(neighs)
            for j in (i+1):length(neighs)
                y = neighs[j]
                add_edge!(previous_graph, x.v, PeriodicVertex(y.v, y.ofs .- x.ofs))
            end
        end
    end
    graph = PeriodicGraph(cryst.pge.g)
    vmap = rem_vertices!(graph, toremove_all)
    rev_vmap = zeros(Int, nv(cryst.pge.g))
    for (i,j) in enumerate(vmap)
        rev_vmap[j] = i
    end

    n = length(cryst.types)
    toremove = falses(n)
    toremove[toremove_all] .= true
    mat = Float64.(cryst.pge.cell.mat)

    new_edgs = Dict{PeriodicEdge3D,Bool}()
    for Q in toremove_list
        @toggleassert all(x -> rev_vmap[x] == 0, Q)
        handled = Set{PeriodicEdge3D}()
        for u in Q
            neighs = PeriodicVertex3D[]
            metalneighs = PeriodicVertex3D[]
            for x in neighbors(previous_graph, u)
                if toremove[x.v]
                    push!(metalneighs, x)
                else
                    push!(neighs, x)
                end
            end
            isempty(neighs) && continue
            num_targets = length(neighs)
            # @show neighs
            visited_coplanar4 = Set{SVector{4,Int}}()
            keep_edges = edges_of_convex_hull(neighs, num_targets, mat, cryst.pge.pos, toremove, visited_coplanar4)
            update_new_edgs!(new_edgs, neighs, keep_edges, toremove, rev_vmap)
            #=
            for metal in metalneighs
                e = PeriodicEdge3D(u, metal)
                reverse(e) ∈ handled && continue
                push!(handled, e)
                Q2 = [(PeriodicVertex(u, .- metal.ofs), metal, 1)]
                encounteredofss = Dict{Int,SVector{3,Int}}((x.v => x.ofs) for x in neighs)
                encountered = Set{PeriodicVertex3D}(neighs)
                copyneighs = copy(neighs)
                for (parent, x, dist) in Q2
                    neighs2 = PeriodicVertex3D[]
                    metalneighs2 = PeriodicVertex3D[]
                    encounteredflags = Union{Nothing,SVector{3,Int}}[]
                    for y in neighbors(cryst.pge.g, x.v)
                        y == parent && continue
                        newvertex = PeriodicVertex3D(y.v, y.ofs .+ x.ofs)
                        newvertex ∈ encountered && continue
                        push!(encountered, newvertex)
                        if toremove[y.v]
                            push!(encounteredflags, get(encounteredofss, newvertex.v, nothing))
                            push!(metalneighs2, newvertex)
                            encounteredofss[newvertex.v] = newvertex.ofs
                        else
                            push!(neighs2, newvertex)
                        end
                    end
                    if !isempty(neighs2)
                        append!(copyneighs, neighs2, metalneighs2)
                    else
                        for (metal2, flag) in zip(metalneighs2, encounteredflags)
                            if !(flag isa Nothing) # to not lose the dimension
                                ofs = metal2.ofs .- flag
                                for _n in neighs
                                    _v = rev_vmap[_n.v]
                                    e = directedge(PeriodicEdge3D(_v, _v, ofs))
                                    new_edgs[e] = true
                                end
                                continue
                            end
                            metal2 == PeriodicVertex3D(u) && continue
                            e = PeriodicEdge3D(u, metal2)
                            #reverse(e) ∈ handled && continue
                            push!(handled, e)
                            if dist < cryst.options.max_polyhedron_radius
                                push!(Q2, (PeriodicVertex(metal.v, metal.ofs .- metal2.ofs), metal2, dist+1))
                            end
                        end
                    end
                end
                # u == 1 && @show metal, copyneighs
                keep_edges = edges_of_convex_hull(copyneighs, num_targets, mat, cryst.pge.pos, toremove,
                                                  copy(visited_coplanar4), keep_edges)
                update_new_edgs!(new_edgs, copyneighs, keep_edges, toremove, rev_vmap)
            end
            =#
        end
    end
    for (e, _b) in new_edgs
        _b && add_edge!(graph, e)
    end

    pge = PeriodicGraphEmbedding(graph, cryst.pge.pos[vmap], cryst.pge.cell)
    #remove_triangles!(graph, pge.pos, nothing, Float64.(cryst.pge.cell.mat), new_edgs)
    types = cryst.types[vmap]


    ret = trimmed_crystal(Crystal{Nothing}(pge, types, Options(cryst.options; _pos=pge.pos)))
    export_default(ret, lazy"clusters_$msg", cryst.options.name, cryst.options.export_clusters)
    return ret
end


"""
    pem_to_pe(cryst::Crystal{Nothing})

Convert `PEM` result to `PE` by removing all metallic sbus.
"""
function pem_to_pe(cryst::Crystal{Nothing})
    n = length(cryst.types)
    visited = falses(n)
    metallic = falses(n)
    metallic_list = Int[]
    kinds = cryst.options.cluster_kinds
    metalkind = getmetal(kinds)
    for (i,t) in enumerate(cryst.types)
        visited[i] && continue
        visited[i] = true
        _metal_identified = identify_metallic_type(t, kinds, metalkind)
        metal_identified = _metal_identified isa Bool ? _metal_identified : begin
            # in this case, it depends on whether there exists a metallic neighbor
            flag = false
            for x in neighbors(cryst.pge.g, i)
                if metallic[x.v]
                    flag = true
                    break
                end
                if !visited[x.v]
                    _met = identify_metallic_type(cryst.types[x.v], kinds, metalkind)
                    if _met isa Bool
                        visited[x.v] = true
                        metallic[x.v] = _met
                        if _met
                            flag = true
                            push!(metallic_list, x.v)
                            break
                        end
                    end
                end
            end
            flag
        end
        if metal_identified
            push!(metallic_list, i)
            metallic[i] = true
        end
    end
    if isempty(metallic_list)
        export_default(cryst, "clusters_PE", cryst.options.name, cryst.options.export_clusters)
        return cryst
    end

    toremove_list = Vector{Int}[]
    tomerge = Int[]
    handled = falses(n)
    metal_surrounded = BitVector(undef, n)
    for i in metallic_list
        handled[i] && continue
        handled[i] = true
        # seen = Dict{Int,SVector{3,Int}}(i => zero(SVector{3,Int}))
        neighs = Set{PeriodicVertex3D}()
        # seenofs = SVector{3,Int}[]
        metal_surrounded .= false
        metal_surrounded_list = Int[]
        Q = [PeriodicVertex3D(i)]
        for u in Q
            surrounded = true
            for x in neighbors(cryst.pge.g, u.v)
                ofs = u.ofs .+ x.ofs
                if !metallic[x.v]
                    push!(neighs, PeriodicVertex(x.v, ofs))
                    surrounded = false
                    continue
                end
                # if handled[x.v]
                #     if ofs != seen[x.v]
                #         push!(seenofs, ofs .- seen[x.v])
                #     end
                # else
                if !handled[x.v]
                    handled[x.v] = true
                    # seen[x.v] = ofs
                    push!(Q, PeriodicVertex(x.v, ofs))
                end
            end
            metal_surrounded[u.v] = surrounded
            surrounded && push!(metal_surrounded_list, u.v)
        end
        isempty(neighs) && continue
        surroundedflag = false
        for v in metal_surrounded_list
            for x in neighbors(cryst.pge.g, v)
                if metal_surrounded[x.v]
                    surroundedflag = true
                    break
                end
            end
            surroundedflag && break
        end
        surroundedflag && continue
        append!(tomerge, metal_surrounded_list)
        push!(toremove_list, [y.v for y in Q])
    end

    return regroup_toremove(cryst, tomerge, toremove_list, "PE")
end


function regroup_vmap(cryst, vmap, isolate, msg)
    clusters, periodicsbus = regroup_sbus(cryst.pge.g, vmap, isolate)
    if !isempty(periodicsbus)
        # Abandon. This can happen for example if a single-C SBU has been split previously
        trimmed_c = trimmed_crystal(cryst)
        export_default(trimmed_c, lazy"clusters_$msg", trimmed_c.options.name, trimmed_c.options.export_clusters)
        return trimmed_c
    end
    edgs = PeriodicEdge3D[]
    for s in vertices(cryst.pge.g)
        atts = clusters.attributions[s]
        for x in neighbors(cryst.pge.g, s)
            d = x.v
            attd = clusters.attributions[d]
            newofs = x.ofs .+ clusters.offsets[s] .- clusters.offsets[d]
            if atts == attd && d != s
                @toggleassert iszero(newofs)
                continue
            end
            push!(edgs, PeriodicEdge3D(atts, attd, newofs))
        end
    end
    graph = PeriodicGraph3D(edgs)

    m = length(clusters.sbus)
    @toggleassert nv(graph) == m
    pos = Vector{SVector{3,Float64}}(undef, m)
    types = Vector{Symbol}(undef, m)
    for (i, sbu) in enumerate(clusters.sbus)
        compo = Vector{Vector{Tuple{Symbol,Int}}}(undef, length(sbu))
        for (k, x) in enumerate(sbu)
            styp = string(cryst.types[x.v])
            thiscompo = Tuple{Symbol,Int}[]
            insymb = true
            currsymb = Symbol("")
            last_j = 1
            for j in 2:length(styp)
                c = styp[j]
                if insymb
                    if isuppercase(c)
                        push!(thiscompo, (Symbol(styp[last_j:j-1]), 1))
                        last_j = j
                    elseif isnumeric(c)
                        insymb = false
                        currsymb = Symbol(styp[last_j:j-1])
                        last_j = j
                    end
                elseif isletter(c)
                    @toggleassert isuppercase(c)
                    insymb = true
                    push!(thiscompo, (currsymb, parse(Int, styp[last_j:j-1])))
                    last_j = j
                end
            end
            if insymb
                push!(thiscompo, (Symbol(styp[last_j:end]), 1))
            else
                @toggleassert isnumeric(last(styp))
                push!(thiscompo, (currsymb, parse(Int, styp[last_j:end])))
            end
            compo[k] = thiscompo
        end
        lengths = [sum(last, comp) for comp in compo]
        total_length = sum(lengths)
        #pos[i] = cryst.pge.pos[first(sbu).v] .+ first(sbu).ofs
        pos[i] = sum((cryst.pge.pos[x.v] .+ x.ofs) .* (lengths[k]/total_length) for (k,x) in enumerate(sbu))

        name = reduce(vcat, compo)
        sort!(name; by=first)
        push!(name, (Symbol(""), 0))
        newname = Tuple{Int,String}[]
        counter = name[1][2]
        for j in 2:length(name)
            if name[j][1] == name[j-1][1]
                counter += name[j][2]
            else
                sym = name[j-1][1]
                str_sym = counter == 1 ? string(sym) : string(sym, counter)
                anum = order_atomtype(sym)
                push!(newname, (anum, str_sym))
                counter = name[j][2]
            end
        end
        sort!(newname; by=first, rev=true)
        types[i] = Symbol(join(last.(newname)))
    end
    ret = trimmed_crystal(Crystal{Nothing}(cryst.pge.cell, types, pos, graph, Options(cryst.options; _pos=pos)))
    export_default(ret, lazy"clusters_$msg", cryst.options.name, cryst.options.export_clusters)
    return ret
end

"""
    allnodes_to_singlenodes(cryst::Crystal)

Convert `AllNodes` result to `SingleNodes` by collapsing all points of extension
clusters bonded together into a new organic cluster.
"""
function allnodes_to_singlenodes(cryst::Crystal{Nothing})
    n = length(cryst.types)
    organics = trues(n)
    for (i,t) in enumerate(cryst.types)
        if t !== :C
            s = string(t)
            el = length(s) ≥ 2 && islowercase(s[2]) ? Symbol(s[1:2]) : Symbol(s[1])
            if ismetal[atomic_numbers[el]] || el === :P || el === :S
                organics[i] = false
            end
        end
    end
    clustername = only(cryst.options.clusterings) == Clustering.Standard ? "standard" : "singlenodes"
    if !any(organics)
        c = trimmed_crystal(cryst)
        export_default(c, lazy"clusters_$clustername", c.options.name, c.options.export_clusters)
        return c
    end
    vmap = zeros(Int, n)
    counter = 0
    inorganic = Int[]
    for i in 1:n
        vmap[i] == 0 || continue
        counter += 1
        if !organics[i]
            vmap[i] = counter
            push!(inorganic, i)
            continue
        end
        encountered = Dict{Int,SVector{3,Int}}(i => zero(SVector{3,Int}))
        Q = [PeriodicVertex3D(i)]
        periodic = false
        for u in Q
            for x in neighbors(cryst.pge.g, u.v)
                if organics[x.v]
                    @toggleassert vmap[x.v] == 0
                    ofs = u.ofs .+ x.ofs
                    if haskey(encountered, x.v)
                        if !periodic && encountered[x.v] != ofs
                            periodic = true
                        end
                    else
                        encountered[x.v] = ofs
                        push!(Q, PeriodicVertex3D(x.v, ofs))
                    end
                end
            end
        end
        for (j, u) in enumerate(Q)
            # if the organic SBU is periodic, do not coalesce its vertices
            vmap[u.v] = ifelse(periodic, counter+j-1, counter)
        end
        counter += length(Q)*periodic
    end

    regroup_vmap(cryst, vmap, inorganic, clustername)
end


# """
#     allnodes_to_standard(cryst::Crystal)

# Convert `AllNodes` result to `Standard` by collapsing all points of extension clusters
# bonded together into a new organic cluster, and removing bonds between metallic clusters.

# To correspond to the actual `Standard` representation, the `separate_metals` option must
# have been set when computing the `AllNodes` representation.
# """
# function allnodes_to_standard(cryst::Crystal)
#     return allnodes_to_singlenodes(cryst)
# end
