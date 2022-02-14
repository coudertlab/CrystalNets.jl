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
    graph = deepcopy(crystal.graph)
    vmap = _trim_monovalent!(graph)
    types = crystal.types[vmap]
    pos = crystal.pos[vmap]
    if T === Nothing
        return Crystal{Nothing}(crystal.cell, types, pos, graph, crystal.options)
    else
        clusters = crystal.clusters[vmap]
        return Crystal{Clusters}(crystal.cell, types, clusters, pos, graph, crystal.options)
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
    function regroup_paddlewheel!(graph, types, clusters)

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


function small_cycles_around(graph, pos, mat, i, u_init, classes, acceptedclasses)
    pos_init = pos[u_init.v]
    posu = [pos_init]
    init_vec = mat * (pos_init .+ u_init.ofs .- pos[i])
    vec = [init_vec]
    prev_vec = [init_vec]
    init_vec = .- init_vec
    angles = [180.0]
    diffangles = [0.0]
    visited = Int[i]
    visited_set = Set{Int}(visited)
    incycles = Set{Int}()
    parent = [i]
    toexplore = [u_init]
    offsets = [zero(SVector{3,Int})]
    while !isempty(toexplore)
        u = pop!(toexplore)
        last_visited = last(visited)
        last_parent = pop!(parent)
        while last_parent != last_visited
            delete!(visited_set, pop!(visited))
            last_visited = last(visited)
        end
        if u.v ∈ visited_set
            for (j, x) in enumerate(visited)
                if x == u.v
                    union!(incycles, visited[j:end])
                    break
                end
            end
            continue
        end
        last_posu = pop!(posu)
        last_angle = pop!(angles)
        last_diffangle = pop!(diffangles)
        last_vec = pop!(vec)
        last_prev_vec = pop!(prev_vec)
        last_offset = pop!(offsets)
        push!(visited, u.v)
        push!(visited_set, u.v)
        for x in neighbors(graph, u.v)
            (x.v == last_visited || degree(graph, x.v) == 1) && continue
            classes[x.v] ∈ acceptedclasses || continue
            new_vec = mat * (pos[x.v] .+ x.ofs .- last_posu)
            α = angle(last_vec, .-new_vec)
            100 < α < 145 || continue
            #β = dihedral(last_prev_vec, last_vec, new_vec)
            #β < 20 || β > 160 || continue
            ofs = last_offset .+ x.ofs
            γ = angle(init_vec, mat * (pos[x.v] .+ ofs .- pos_init))
            if last_parent == i
                last_diffangle = (180 - γ)/2
                last_angle = γ + last_diffangle
            end
            γ < last_angle || continue
            diffangle = last_angle - γ
            abs(diffangle - last_diffangle) < last_diffangle/5 || continue
            push!(toexplore, PeriodicVertex3D(x.v, ofs))
            push!(parent, u.v)
            push!(posu, pos[x.v])
            push!(angles, γ)
            push!(diffangles, diffangle)
            push!(prev_vec, last_vec)
            push!(vec, new_vec)
            push!(offsets, ofs)
        end
    end
    # !isempty(incycles) && @show incycles
    return incycles
end

"""
    in_small_cycles_around(graph, pos, mat, i, classes, acceptedclasses)

Return the set of atoms belonging to a small cycle to which also belongs atom `i`.
This cycle must only contain atoms whose class is in `acceptedclasses`
"""
function in_small_cycles_around(graph, pos, mat, i, classes, acceptedclasses)
    neighs = neighbors(graph, i)
    _, state = iterate(neighs)
    incycle = Set{Int}()
    for u in Base.rest(neighs, state)
        degree(graph, u.v) == 1 && continue
        union!(incycle, small_cycles_around(graph, pos, mat, i, u, classes, acceptedclasses))
    end
    return incycle
end
 

function detect_organiccycles(classes, graph, pos, mat, Cclass, modifiables)
    handled = Set{Int}()
    cycles = Set{Int}[]
    acceptedclasses = union(modifiables, Cclass)
    for (i, class) in enumerate(classes)
        class == Cclass || continue
        (i ∈ handled || degree(graph, i) ≤ 1) && continue
        skipflag = true
        for x in neighbors(graph, i)
            if classes[x.v] ∈ acceptedclasses && degree(graph, x.v) > 1
                skipflag = false
                break
            end
        end
        skipflag && continue
        incycle = in_small_cycles_around(graph, pos, mat, i, classes, acceptedclasses)
        if !isempty(incycle)
            union!(handled, incycle)
            push!(cycles, incycle)
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
        sort!(v)
        minsbu = minimum([union_find_cycles[vj] for vj in v])
        for j in 1:length(v)
            k = union_find_cycles[j]
            (k == minsbu || isempty(rev_union_find_cycles[k])) && continue
            union_find_cycles[j] = k
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



"""
    find_sbus(crystal, kinds=default_sbus)

Recognize SBUs using heuristics based on the atom types corresponding to the `Intermediate`
clustering algorithm.
"""
function find_sbus(crystal, kinds=default_sbus, clustering=crystal.options.clustering)
    separate_metals = crystal.options.separate_metals isa Bool ? _crystal.options.separate_metals :
        ((clustering == Clustering.Standard) | (clustering == Clustering.PEM))
    n = nv(crystal.graph)
    classes = Vector{Int}(undef, n)
    for i in 1:n
        atom_name = crystal.types[i]
        class = kinds[atom_name]
        if class == 0 # atom not identified
            if atom_name === Symbol("")
                throw(MissingAtomInformation("""
                the input is a periodic graph with no atom information, it cannot be recognized as a MOF.
                """))
            elseif atom_name === Symbol("_")
                throw(MissingAtomInformation("""
                the input file format does not contain enough information on the atoms to distinguish the organic and inorganic SBUs.
                Please use a file format containing at least the atom types and in order to use MOF recognition.
                """))
            else
                throw(MissingAtomInformation("unhandled atom type: $(element_categories[atomic_numbers[atom_name]]) (for atom $atom_name)."))
            end
        end
        classes[i] = class
    end

    temporary_classes = copy(kinds.tomerge)
    temporary_classes_set = Set{Int}(temporary_classes)
    nclasses::Int = length(kinds)
    last_class = nclasses
    graph = crystal.graph
    organickind = kinds[:C]

    if crystal.options.detect_organiccycles
        organiccycle = detect_organiccycles(classes, graph, crystal.pos,
                                    crystal.cell.mat, organickind, temporary_classes_set)
        same_SBU = group_cycle(organiccycle, crystal.types, graph)
        if !isempty(same_SBU)
            graph = PeriodicGraph(crystal.graph) # make a modificable copy
            for cycle in same_SBU
                # Each atom in an aromatic cycle is bonded to the other atoms of the
                # same cycle. The newly-formed clique is assigned a new class.
                last_class += 1
                for (i, x) in enumerate(cycle)
                    classes[x.v] = last_class
                    for j in (i+1):length(cycle)
                        y = cycle[j]
                        add_edge!(graph, x.v, PeriodicVertex3D(y.v, y.ofs .- x.ofs))
                    end
                end
            end
        end
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
                    crystal.types[u] = :Pc # mark an organic P
                end
            end
        end
    end

    if crystal.options.detect_points_of_extension
        for (i, class) in enumerate(classes)
            class == organickind || continue
            for x in neighbors(graph, i)
                if classes[x.v] != organickind
                    last_class += 1
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
                    metallics = [k for (at, k) in zip(_ats, uniquecompo) if ismetal[at] || k === :P]
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
        # export_default(Crystal{Nothing}(crystal.cell, Symbol.(classes), crystal.pos, graph, crystal.options))
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
    graph = deepcopy(c.graph)
    for (k, t) in enumerate(c.types)
        (t === :O && degree(graph, k) > 2) || continue
        _split_this_sbu!(toremove, graph, k, c.types, :O, nothing)
    end
    vmap = rem_vertices!(graph, toremove)
    pos = c.pos[vmap]
    return Crystal{Nothing}(c.cell, c.types[vmap], pos, graph, Options(c.options; _pos=pos))
end


struct InvalidSBU <: ClusteringError
    msg::String
end


find_clusters(c::Crystal, structure=c.options.structure) = find_clusters(c, structure, c.options.clustering)
function find_clusters(c::Crystal{T}, structure::_StructureType, clustering::_Clustering)::Tuple{Clusters,_StructureType,_Clustering} where T
    if clustering == Clustering.Auto
        if structure == StructureType.Auto
            if T === Clusters
                return c.clusters, structure, clustering
            else
                return find_clusters(c, structure, Clustering.EachVertex)
            end
        elseif structure == StructureType.Zeolite
            return find_clusters(c, structure, Clustering.EachVertex)
        elseif structure == StructureType.Guess
            crystal = Crystal{Nothing}(c)
            try
                return find_clusters(crystal, StructureType.Cluster, clustering)
            catch e
                if !(e isa ClusteringError)
                    rethrow()
                end
            end
            return find_clusters(c, StructureType.Auto, clustering)
        elseif structure == StructureType.MOF || structure == StructureType.Cluster
            return find_clusters(c, structure, Clustering.Intermediate)
        end
    elseif clustering == Clustering.EachVertex
        return Clusters(length(c.types)), structure, clustering
    elseif clustering == Clustering.Input
        if T === Nothing
            throw(ArgumentError("Cannot use input residues as clusters: the input does not have residues."))
        else
            return c.clusters, structure, clustering
        end
    else
        @toggleassert clustering == Clustering.Intermediate ||
                      clustering == Clustering.AllNodes     ||
                      clustering == Clustering.SingleNodes  ||
                      clustering == Clustering.Standard     ||
                      clustering == Clustering.PEM
        if clustering == Clustering.AllNodes || clustering == Clustering.SingleNodes
            clustering = Clustering.Intermediate
        end
        clusters = find_sbus(c, c.options.cluster_kinds, clustering)
        return clusters, structure, clustering
    end
end


function order_atomtype(sym)
    sym === :C && return 119 # C is put first to identify organic clusters
    sym === :P && return 1 # Put P before the other non-metals, but after everything else
    anum = atomic_numbers[sym]
    elcat = element_categories[anum]
    if elcat === :nonmetal || elcat === :noble || elcat === :halogen
        return -anum # Put non-metals at last
    end
    return anum
end

"""
    collapse_clusters(crystal::Crystal)

Return the new crystal corresponding to the input where each cluster has been
transformed into a new vertex.
"""
function collapse_clusters(c::Crystal, _structure::_StructureType=c.options.structure,
                                       _clustering::_Clustering=c.options.clustering)
    crystal = trim_monovalent(c)
    clusters, structure, clustering = find_clusters(crystal, _structure, _clustering)
    if clustering == Clustering.EachVertex || structure == StructureType.Zeolite
        if crystal.options.split_O_vertex
            return split_O_vertices(crystal), structure, clustering
        end
        return Crystal{Clusters}(crystal, clusters; _pos=crystal.pos), structure, clustering
    end
    edgs = PeriodicEdge3D[]
    for s in vertices(crystal.graph)
        atts = clusters.attributions[s]
        neigh0 = neighbors(crystal.graph, s)
        @toggleassert length(neigh0) ≥ 2
        typisnotC = crystal.types[s] !== :C
        for x in neigh0
            d = x.v
            if crystal.options.bond_adjacent_sbus && crystal.types[d] === :C && typisnotC
                neighs = neighbors(crystal.graph, d)
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
    export_default(crystal, "trimmed", crystal.options.name,
                   crystal.options.export_input)
    n = length(clusters.sbus)
    graph = PeriodicGraph3D(n, edgs)
    if _structure == StructureType.Guess && nv(graph) == 1
        return collapse_clusters(crystal, StructureType.Auto)
    end
    if ne(graph) == 0
        return (Crystal{Clusters}(crystal.cell, Symbol[], clusters, SVector{3,Float64}[], graph,
                                Options(crystal.options; _pos=SVector{3,Float64}[])),
                structure, clustering)
    end
    if !isempty(crystal.options.export_attributions)
        path = tmpexportname(crystal.options.export_attributions, "attribution_", crystal.options.name, ".pdb")
        export_attributions(Crystal{Clusters}(crystal, clusters), path)
        println("Attributions of atoms into SBUs represented represented at ", replace(path, ('\\'=>'/')))
    end

    sbus = clusters.sbus[split_special_sbu!(graph, clusters.sbus, crystal.graph, crystal.types, crystal.options.split_O_vertex)]
    n = length(sbus)
    pos = Vector{SVector{3,Float64}}(undef, n)
    types = Vector{Symbol}(undef, n)
    for (i, sbu) in enumerate(sbus)
        pos[i] = mean(crystal.pos[x.v] .+ x.ofs for x in sbu)
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
    ret = Crystal{Clusters}(crystal.cell, types, clusters, pos, graph, Options(crystal.options; _pos=pos))
    if crystal.options.clustering == Clustering.Intermediate ||
            (clustering != Clustering.Intermediate && 
            crystal.options.clustering != Clustering.Standard)
        export_default(ret, "clusters", crystal.options.name, crystal.options.export_clusters)
    end
    return ret, structure, clustering
end

function _convex_hull(neighs, keep_edges, toremove, rev_vmap)
    new_edgs = PeriodicEdge3D[]
    for (i, x) in enumerate(neighs)
        toremove[x.v] && continue
        for j in (i+1):length(neighs)
            keep_edges[i,j] || continue
            y = neighs[j]
            toremove[y.v] && continue
            push!(new_edgs, PeriodicEdge3D(rev_vmap[x.v], rev_vmap[y.v], y.ofs .- x.ofs))
        end
    end
    return new_edgs
end

function convex_hull(neighs, pos, toremove, rev_vmap)
    n = length(neighs)
    keep_edges = trues(n, n)
    n ≤ 3 && return _convex_hull(neighs, keep_edges, toremove, rev_vmap)
    poss = [pos[x.v] .+ x.ofs for x in neighs]
    visited_coplanar4 = Set{SVector{4,Int}}()
    convex_check_3D = n ≥ 5
    circ_noref = collect(2:n)
    for (i, x) in enumerate(neighs)
        circ_noref[ifelse(i == 1, 1, i-1)] = i-1
        i == n || (circ_noref[i] = i+1)
        ref = poss[i]
        circ_noj1 = circ_noref[2:end]
        for (_j1, j1) in enumerate(circ_noref)
            _j1 == 1 || (circ_noj1[_j1-1] = circ_noref[_j1-1])
            _j1 == n-1 || (circ_noj1[_j1] = circ_noref[_j1+1])
            vec1 = poss[j1] .- ref
            circ_noj2 = circ_noj1[2:end]
            for (_j2, j2) in enumerate(circ_noj1)
                _j2 == 1 || (circ_noj2[_j2-1] = circ_noj1[_j2-1])
                _j2 == n-2 || (circ_noj2[_j2] = circ_noj1[_j2+1])
                vec2 = poss[j2] .- ref
                for j3 in circ_noj2
                    ordered4 = sort(SVector{4,Int}((i, j1, j2, j3)))
                    ordered4 ∈ visited_coplanar4 && continue
                    vec3 = poss[j3] .- ref
                    β = dihedral(vec1, vec2, vec3)
                    if β < 10 || β > 170
                        push!(visited_coplanar4, ordered4)
                        couples = [(i,j1), (i,j2), (i,j3), (j1,j2), (j1,j3), (j2,j3)]
                        lens = [norm(vec1), norm(vec2), norm(vec3), norm(poss[j2] .- poss[j1]),
                                norm(poss[j3] .- poss[j1]), norm(poss[j3] .- poss[j2])]
                        I = sortperm(lens)
                        k1, k2 = couples[pop!(I)]
                        k1, k2 = minmax(k1, k2)
                        keep_edges[k1, k2] = false
                        k3, k4 = couples[pop!(I)]
                        k3, k4 = minmax(k3, k4)
                        keep_edges[k3, k4] = false
                        n == 4 && return _convex_hull(neighs, keep_edges, toremove, rev_vmap)
                    elseif convex_check_3D
                        otherneighs = [k for k in circ_noj2 if k > i && k != j3 && keep_edges[i,k]]
                        isempty(otherneighs) && continue
                        mat = inv(hcat(vec1, vec2, vec3))
                        for k in otherneighs
                            coeffs = mat * (poss[k] .- ref)
                            keep_edges[i,k] = !all(≥(-3e-16), coeffs)
                        end
                    end
                end
            end
        end
    end
    return _convex_hull(neighs, keep_edges, toremove, rev_vmap)
end

function regroup_toremove(cryst, toremove_list, bond_neighbors, msg)
    graph = PeriodicGraph(cryst.graph)
    vmap = rem_vertices!(graph, toremove_list)
    rev_vmap = zeros(Int, nv(cryst.graph))
    for (i,j) in enumerate(vmap)
        rev_vmap[j] = i
    end

    n = length(cryst.types)
    toremove = falses(n)
    toremove[toremove_list] .= true

    new_edgs = PeriodicEdge3D[]
    for (u, b) in zip(toremove_list, bond_neighbors)
        b || continue
        @toggleassert rev_vmap[u] == 0
        neighs = neighbors(cryst.graph, u)
        append!(new_edgs, convex_hull(neighs, cryst.pos, toremove, rev_vmap))
    end
    for e in new_edgs
        add_edge!(graph, e)
    end

    pos = cryst.pos[vmap]
    #remove_triangles!(graph, pos, nothing, Float64.(cryst.cell.mat), new_edgs)
    types = cryst.types[vmap]
    

    ret = trimmed_crystal(Crystal{Nothing}(cryst.cell, types, pos, graph, Options(cryst.options; _pos=pos)))
    export_default(ret, "clusters_$msg", cryst.options.name, cryst.options.export_clusters)
    return ret
end


function identify_metallic_type(t, kinds, metalkind=getmetal(kinds))
    at = representative_atom(t)
    ismetal[at] && return true
    k = kinds[at]
    k == metalkind && return true
    k ∈ kinds.tomerge && return missing
    return false
end

"""
    intermediate_to_allnodes(cryst::Crystal{Nothing})

Convert `Intermediate` result to `AllNodes` by removing all periodic metallic sbus.
"""
function intermediate_to_allnodes(cryst::Crystal{Nothing})
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
            for x in neighbors(cryst.graph, i)
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
        export_default(cryst, "clusters_allnodes", cryst.options.name, cryst.options.export_clusters)
        return cryst
    end

    toremove_list = Int[]
    bond_neighbors = Bool[]
    handled = falses(n)
    for i in metallic_list
        handled[i] && continue
        handled[i] = true
        remove_flag = false
        bond_neighbor_flag = true
        seen = Dict{Int,SVector{3,Int}}(i => zero(SVector{3,Int}))
        neighs = Set{PeriodicVertex3D}()
        Q = [PeriodicVertex3D(i)]
        for u in Q
            for x in neighbors(cryst.graph, u.v)
                ofs = u.ofs .+ x.ofs
                if !metallic[x.v]
                    push!(neighs, PeriodicVertex(x.v, ofs))
                    continue
                end
                if handled[x.v]
                    if !remove_flag && ofs != seen[x.v]
                        remove_flag = true
                    end
                else
                    handled[x.v] = true
                    seen[x.v] = ofs
                    push!(Q, PeriodicVertex(x.v, ofs))
                end
            end
        end
        if !remove_flag && !isempty(neighs)
            a = first(neighs)
            P = [a]
            seenneighs = Set{PeriodicVertex3D}(P)
            for u in P
                for x in neighbors(cryst.graph, u.v)
                    y = PeriodicVertex(x.v, u.ofs .+ x.ofs)
                    y ∈ seenneighs && continue
                    if y ∈ neighs
                        push!(seenneighs, y)
                        push!(P, y)
                    end
                end
            end
            bond_neighbor_flag = remove_flag = length(neighs) == length(seenneighs)
        end
        if remove_flag && !isempty(neighs)
            append!(toremove_list, u.v for u in Q)
            append!(bond_neighbors, bond_neighbor_flag for _ in Q)
        end
    end

    return regroup_toremove(cryst, toremove_list, bond_neighbors, "allnodes")
end


function regroup_vmap(cryst, vmap, isolate, msg)
    clusters, periodicsbus = regroup_sbus(cryst.graph, vmap, isolate)
    if !isempty(periodicsbus)
        # Abandon. This can happen for example if a single-C SBU has been split previously
        export_default(cryst, "clusters_$msg", cryst.options.name, cryst.options.export_clusters)
        return cryst
    end
    edgs = PeriodicEdge3D[]
    for s in vertices(cryst.graph)
        atts = clusters.attributions[s]
        for x in neighbors(cryst.graph, s)
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
        #pos[i] = cryst.pos[first(sbu).v] .+ first(sbu).ofs
        pos[i] = sum((cryst.pos[x.v] .+ x.ofs) .* (lengths[k]/total_length) for (k,x) in enumerate(sbu))

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
    ret = trimmed_crystal(Crystal{Nothing}(cryst.cell, types, pos, graph, Options(cryst.options; _pos=pos)))
    export_default(ret, "clusters_$msg", cryst.options.name, cryst.options.export_clusters)
    return ret
end

"""
    intermediate_to_singlenodes(cryst::Crystal)

Convert `Intermediate` result to `SingleNodes` by collapsing all points of extension
clusters bonded together into a new organic cluster.
"""
function intermediate_to_singlenodes(cryst::Crystal)
    n = length(cryst.types)
    organics = trues(n)
    for (i,t) in enumerate(cryst.types)
        if t !== :C
            s = string(t)
            el = length(s) ≥ 2 && islowercase(s[2]) ? Symbol(s[1:2]) : Symbol(s[1])
            if ismetal[atomic_numbers[el]] || el === :P
                organics[i] = false
            end
        end
    end
    if !any(organics)
        export_default(cryst, "clusters_singlenodes", cryst.options.name, cryst.options.export_clusters)
        return cryst
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
            for x in neighbors(cryst.graph, u.v)
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

    regroup_vmap(cryst, vmap, inorganic, cryst.options.clustering == Clustering.Standard ? "standard" : "singlenodes")
end
