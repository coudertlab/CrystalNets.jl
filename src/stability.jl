## Functions related to unstable nets (some vertices have the same equilibrium placement)


"""
    shrink_collisions(net::CrystalNet, collisions)

Remove all colliding vertices and replace them by one new vertex per `CollisionNode`, whose
neighbours are that of the vertices within.
"""
function shrink_collisions(net::CrystalNet{D,T}, collisions::Vector{UnitRange{Int}}) where {D,T}
    first_colliding = first(first(collisions))
    collision_vmap = collect(1:length(net.types))
    for (i, node) in enumerate(collisions)
        @toggleassert all(==(net.pos[first(node)]), @view net.pos[node])
        collision_vmap[node] .= i
    end
    edgs = PeriodicEdge{D}[]
    for e in edges(net.graph)
        push!(edgs, collision_vmap[e.src], collision_vmap[e.dst.v], e.dst.ofs)
    end
    newgraph = PeriodicGraph{D}(edgs)

    m = length(collisions)
    n = first_colliding + m - 1
    @toggleassert nv(newgraph) == n

    newpos = net.pos[1:first_colliding]
    append!(newtypes, net.pos[first(collisions[i])] for i in 2:m)
    newtypes = net.types[1:first_colliding-1]
    append!(newtypes, :* for i in 1:m) # type of collision nodes will be Symbol('*')

    newnet = CrystalNet{D,T}(net.cell, newtypes, newpos, newgraph, net.options)

    return newnet, collision_vmap
end


function _order_collision(subgraph, subnodes)
    m = nv(subgraph)
    if m == 2 # 2
        return [1, 2]
    end

    if m == 3
        length(subnodes) == 1 && return sortperm(degree(subgraph))
        @toggleassert length(subnodes) == 2 # 1-2 or 2-1
        singled = length(subnodes[1]) == 1 ? subnodes[1][1] : subnodes[2][1]
        v = get(neighbors(subgraph, singled), 1, singled == 1 ? 2 : 1)
        remaining = 6 - singled - v
        return length(subnodes[1]) == 1 ? [singled, v, remaining] : [v, remaining, singled]
    end

    # m == 4
    degs = degree(subgraph)
    @toggleassert all(<(4), degs)
    J = sortperm(degs)
    num = length(subnodes)

    if num == 1
        if degs[J[1]] == degs[J[4]]
            if degs[1] == 1 # 1-1-1-1
                neigh1 = only(neighbors(subgraph, 1))
                other1 = ifelse(neigh1 == 2, 3, 2)
                other2 = 10 - neigh1 - other1 - 1
                @toggleassert only(neighbors(subgraph, other1)) == other2
                return [1, neigh1, other1, other2]
            end
            if degs[1] == 2 # 2-2-2-2
                neighs1 = neighbors(subgraph, 1)
                neighl = neighs1[1]
                neighr = neighs1[2]
                neigho = 10 - 1 - neighl - neighr
                @toggleassert issetequal(neighbors(subgraph, neigho), [neighl, neighr])
                return [neighl, 1, neighr, neigho]
            end
        end
        if degs[J] == [1, 1, 2, 2]
            neigha = only(neighbors(subgraph, J[1]))
            neighb = only(neighbors(subgraph, J[2]))
            @toggleassert issetequal(neighbors(subgraph, neigha), [J[1], neighb])
            return [J[1], J[2], neigha, neighb]
        end
        return J
    end

    if num == 2
        if length(subnodes[1]) == 2 # 2-2
            a, b = subnodes[1]
            rem_edge!(subgraph, a, b)
            c, d = subnodes[2]
            rem_edge!(subgraph, c, d)
            if degree(subgraph, b) < degree(subgraph, a)
                a, b = b, a
            end
            if (degree(subgraph, b) == 1 && only(neighbors(subgraph, b)) == c) ||
                degree(subgraph, d) < degree(subgraph, c)
                c, d = d, c
            end
            return [a, b, c, d]
        end
        # 1-3 or 3-1
        i1 = length(subnodes[1]) == 1 ? 1 : 2
        (t, others) = only(subnodes[i1]), subnodes[3-i1]
        u, v, w = others[sortperm(degs[others])]
        @toggleassert issetequal([u,v,w,t], 1:4)
        nts = neighbors(subgraph, t)
        dt = length(nts)
        if dt == 1
            nt = only(nts)
            if degs[v] == degs[w]
                if degs[w] == 1
                    @toggleassert degs[u] == 1
                    if nt == w
                        w, v = v, w
                    end
                    if nt == v
                        u, v = v, u
                    end
                elseif nt == v
                    w, v = v, w
                end
            end
        elseif dt == 2
            if degs[u] == degs[v] 
                if degs[u] == 1
                    @toggleassert degs[w] == degs[v]+1 == 2
                    @toggleassert w ∈ nts
                    if v ∈ nts
                        u, v = v, u
                    end
                else
                    @toggleassert degs[u] == 2
                    if degs[v] == 2
                        @toggleassert degs[w] == 2
                        ww = 10 - t - nts[1] - nts[2]
                        if ww == u
                            u, w = w, u
                        elseif ww == v
                            v, w = w, v
                        end
                    end
                end
            end
        end
        return i1 == 1 ? [t, u, v, w] : [u, v, w, t]
    end

    @toggleassert num == 3 # 1-1-2 or 1-2-1 or 2-1-1
    n2 = length(subnodes[1]) == 2 ? 1 : length(subnodes[2]) == 2 ? 2 : 3
    xi, xj = subnodes[n2]
    rem_edge!(subgraph, xi, xj)
    xk = only(subnodes[n2 == 1 ? 2 : 1])
    xl = only(subnodes[n2 == 1 ? 3 : 5-n2])
    rem_edge!(subgraph, xk, xl)
    neighk = neighbors(subgraph, xk)
    if length(neighk) == 1
        if only(neighk) == xi
            xi, xj = xj, xi
        end
    else
        neighl = neighbors(subgraph, xl)
        if length(neighl) == 1 && only(neighl) == xi
            xi, xj = xj, xi
        end
    end
    return n2 == 1 ? [xi, xj, only(subnodes[2]), only(subnodes[3])] :
           n2 == 2 ? [only(subnodes[1]), xi, xj, only(subnodes[3])] :
                     [only(subnodes[1]), only(subnodes[2]), xi, xj]
end

"""
    order_collision(graph::PeriodicGraph, rnge)

Given collision node (in the form of the corresponding list of colliding vertices), find
an ordering of them which is independent of the current ordering of these vertices and of
vertices which are neither in the collision node nor any of its neighbours.
Return an empty list if it fails.

This function assumes that no vertex in the node has a neighbour in another collision node
and that there are no two representatives of the same vertex that are neighbour to some
vertices in the range.
"""
function order_collision(graph::PeriodicGraph{D}, colliding) where D
    m = length(colliding)
    neigh_per_v = [Int[] for _ in 1:m]
    subgraph = SimpleGraph(m)
    for (i, (nlist, u)) in enumerate(zip(neigh_per_v, colliding))
        for x in neighbors(graph, u)
            j = findfirst(==(x.v), colliding)
            if j isa Int
                add_edge!(subgraph, i, j)
            else
                push!(nlist, x.v)
            end
        end
        sort!(nlist)
    end
    I = sortperm(neigh_per_v)
    subnodes = Vector{Int}[[I[1]]]
    last_list = neigh_per_v[I[1]]
    for _i in 2:length(I)
        i = I[_i]
        nlist = neigh_per_v[i]
        if nlist == last_list
            push!(last(subnodes), i)
        else
            push!(subnodes, [i])
            last_list = nlist
        end
    end
    length(subnodes) == m && return I

    m > 4 && return Int[]
    return _order_collision(subgraph, subnodes)
end

"""
    expand_collisions(collisions, graph::PeriodicGraph, initial_graph, vmap, collision_vmap)

Expand each collision node into the appropriate number of vertices so that the resulting
graph is isomorphic to the initial one, in a manner that only depends on the current graph.
Return the resulting graph.

`vmap` is the map of vertices between `initial_graph` (with collapsed collision nodes) and
`graph`
"""
function expand_collisions(collisions::Vector{UnitRange{Int}}, graph::PeriodicGraph{D}, initial_graph, vmap, collision_vmap) where D
    n = nv(initial_graph)
    m = length(collisions)
    num_withcolliding = length(vmap)
    num_nocolliding = num_withcolliding - m

    @toggleassert num_nocolliding + sum(x->length(x), collisions) == length(collision_vmap) == n
    @toggleassert nv(graph) == num_withcolliding

    rev_vmap = Vector{Int}(undef, num_withcolliding)
    for (i, j) in enumerate(vmap)
        rev_vmap[j] = i
    end

    prev_ranges = Vector{Int}(undef, m+1)
    prev_ranges[1] = num_nocolliding
    for (i,node) in enumerate(collisions)
        @toggleassert first(node) == prev_ranges[i] + 1
        prev_ranges[i+1] = last(node)
    end

    colliding_nodes = [rev_vmap[collision_vmap[prev_ranges[i]]] for i in 2:(m+1)]
    I = sortperm(colliding_nodes)
    colliding_nodes = colliding_nodes[I]
    collisions = collisions[I] # The order of collisions is now independent of the initial representation

    _new_ranges = num_nocolliding .+ cumsum(length(node) for node in collisions)
    pushfirst!(_new_ranges, num_nocolliding)
    new_ranges = [(1+_new_ranges[i]):new_ranges[i+1] for i in 1:m]

    # We now push all the collision nodes to the end
    perm_nocolliding = Vector{Int}(undef, num_nocolliding)
    perm_withcolliding = Vector{Int}(undef, num_withcolliding)
    perm = Vector{Int}(undef, n)
    current_idx_colliding = 1
    for i in 1:nv(graph)
        if i == colliding_nodes[current_idx_colliding]
            perm_withcolliding[num_nocolliding+current_idx_colliding] = i
            newpos = new_ranges[current_idx_colliding]
            oldpos = (1+prev_ranges[current_idx_colliding]):prevprev_rangessums[current_idx_colliding+1]
            perm[newpos] = oldpos

            current_idx_colliding += 1
        else
            idx = i + 1 - current_idx_colliding
            perm_withcolliding[idx] = perm_nocolliding[idx] = i
            perm[idx] = vmap[i]
        end
    end
    newgraph = graph[perm_nocolliding] # subgraph of the non-colliding nodes keeping the order]
    graph = graph[perm_withcolliding]
    initial_graph = initial_graph[perm]

    # At this point, the graph has all the collision nodes stacked at its end, but kept in
    # the same order as that given by the topological genome of the graph.

    add_vertices!(newgraph, n - num_nocolliding)
    for i in 1:m
        rnge = new_ranges[i]
        for x in neighbors(graph, i) # x is a (non-colliding) neighbour of the collision node
            for j in neighbors(initial_graph, x.v)
                if j.v ∈ rnge # j.v is the number of the vertex in the collision node whose neighbour is x
                    add_edge!(newgraph, PeriodicEdge{D}(j.v, x))
                end
            end
        end
    end

    # newgraph contains all vertices from initial_graph
    # The non-colliding ones are ordered like in graph, the others are stacked in the order
    # of their collision node in graph.
    # What remains is the ordering of the vertices within each collision node.

    newperm = collect(1:n)
    for rnge in new_ranges
        neworder = order_collision(newgraph, rnge)
        isempty(neworder) && error("Internal error: invalid unstable node detected at expansion")
        newperm[rnge] = (rnge[1] - 1) .+ neworder
    end
    retgraph = newgraph[newperm]
    newperm[1:num_nocolliding] = vmap[perm_nocolliding]

    return newperm, retgraph
end


function subgraph_identifier(graph, node, vmap=nothing)
    neighbors_of_node_set = Set{Int}
    for u in node
        for x in neighbors(graph, u)
            push!(neighbors_of_node_set, x.v)
        end
    end
    setdiff!(neighbors_of_node_set, node)
    neighbors_of_node = collect(neighbors_of_node_set)
    sort!(neighbors_of_node; by=(vmap isa Vector{Int} ? (x->vmap[x]) : identity))

    newmap = Vector{Int}(undef, nv(graph))
    for (i, u) in enumerate(node)
        newmap[u] = i
    end
    l = length(node)
    for (i, u) in enumerate(neighbors_of_node)
        newmap[u] = l + i
    end

    edgs = PeriodicEdge{0}[]
    for u in enumerate(node)
        v = newmap[u]
        for x in neighbors(graph, u)
            push!(edgs, PeriodicEdge{0}(v, newmap[x.v], ()))
        end
    end

    newgraph = PeriodicGraph{0}(edgs)
    n = nv(newgraph)
    @toggleassert n == length(neighbors_of_node) + length(node)
    neworder = order_collision(newgraph, 1:length(node))
    isempty(neworder) && return ""
    perm = collect(1:n)
    perm[1:length(node)] = neworder
    return string(newgraph[perm])
end


function collect_collisions(net)
    @toggleassert issorted(net.pos)
    @toggleassert all(pos -> all(x -> 0 ≤ x < 1, pos), net.pos)
    collision_sites = Vector{Int}[]
    collision_vertices = Int[]
    flag_collision = false
    for i in 2:length(net.pos)
        if net.pos[i-1] == net.pos[i]
            if flag_collision
                push!(collision_sites[end], i)
                push!(collision_vertices, i)
            else
                push!(collision_sites, [i-1, i])
                push!(collision_vertices, i-1, i)
                flag_collision = true
            end
        else
            flag_collision = false
        end
    end
    return collision_sites, collision_vertices
end

"""
    collision_nodes(c::CrystalNet)

Check that the net is stable, i.e. that no two vertices have the same equilibrium placement.

A net is still considered stable if the collisions in equilibrium placement cannot lead to
different topological genomes. In practice, this happens when:
A) there is no edge between two collision sites and
B) there is no edge between a collision site and two representatives of the same vertex and
C) for each collision site either:
   - the site is made of at most 4 vertices, or
   - no 2 vertices on the site share the same exact set of neighbours out of the site.

In this case, return the list of corresponding `CollisionNode`, the list being empty if the
net is truly stable. Otherwise, return `nothing`.

This function reorders the vertices in the net so that the vertices in a `CollisionNode`
appear after the others and are contiguous when they belong to the same `CollisionNode`.
"""
function collision_nodes(net::CrystalNet{D}) where D
    collision_sites, collision_vertices = collect_collisions(net)

    isempty(collision_sites) && return net, UnitRange{Int}[] # no collision at all

    # Check that conditions A, B and C are respected
    collision_vertices_set = BitSet(collision_vertices)
    for site in collision_sites
        known_neighbors = Dict{Int,SVector{D,Int}}()
        known_nlist = Set{Vector{PeriodicVertex{D}}}()
        uniquenlist = true
        for u in site
            this_nlist = PeriodicVertex{D}[]
            for x in neighbors(net.graph, u)
                x.v ∈ site && continue
                x.v ∈ collision_vertices_set && return net, nothing # condition A
                get!(known_neighbors, x.v, x.ofs) == x.ofs || return net, nothing # condition B
                if uniquenlist
                    push!(this_nlist, x)
                end
            end
            if uniquenlist
                if this_nlist ∈ known_nlist
                    uniquenlist = false
                else
                    push!(known_nlist, this_nlist)
                end
            end
        end
        !uniquenlist && length(site) > 4 && return net, nothing # condition C
    end

    # Reorder the vertices of the net so that those belonging to a collision node are
    # contiguous, all collision nodes are sorted according to length and they appear after
    # non-colliding vertices.
    sort!(collision_sites; by=length)
    n = length(net.types)
    vmap = Vector{Int}(undef, n)
    idx = 1
    for i in 1:n
        if i ∉ collision_vertices_set
            vmap[idx] = i
            idx += 1
        end
    end
    @toggleassert idx == n - length(collision_vertices) + 1
    collisions = Vector{UnitRange{Int}}(undef, length(collision_sites))
    for (k, node) in enumerate(collision_sites)
        collisions[k] = idx:(idx+length(node)-1)
        for x in node
            vmap[idx] = x
            idx += 1
        end
    end
    @toggleassert idx == n+1

    return CrystalNet(net.cell, net.types[vmap], net.pos[vmap], net.graph[vmap], net.options),
           collisions

    # for site in collision_site
    #     adjacent_nodes = Vector{Int}[]
    #     handled = BitSet()
    #     for simplenode in site
    #         simplenode ∈ handled && continue
    #         push!(handled, simplenode)
    #         Q = [simplenode]
    #         ref_ext_neighbours = Set{PeriodicVertex{D}}()
    #         internal_degree = Int[]
    #         for u in Q
    #             ext_neighbours = PeriodicVertex{D}[]
    #             push!(internal_degree, 0)
    #             for x in neighbors(net.graph, u)
    #                 v = x.v
    #                 if iszero(x.ofs)
    #                     if v ∈ Q
    #                         internal_degree[end] += 1
    #                         continue
    #                     end
    #                     if v ∈ site
    #                         push!(Q, v)
    #                         push!(handled, v)
    #                         internal_degree[end] += 1
    #                         continue
    #                     end
    #                 end
    #                 if v ∈ collision_vertices
    #                     return nothing # edge between two collision sites
    #                 else
    #                     push!(ext_neighbours, x)
    #                 end
    #             end
    #         end
    #         @toggleassert !isempty(ext_neighbours) # otherwise, non-periodic input?
    #         if isempty(ref_ext_neighbours)
    #             length(site) > 1 && (ref_ext_neighbours = Set(ext_neighbours))
    #         else
    #             issetequal(ext_neighbours, ref_ext_neighbours) || return nothing # all vertices do not share the same external neighbours
    #         end
    #         push!(adjacent_nodes, Q)
    #     end
    # end

    #=
    each `vertices` field must be sorted
    The list of collision nodes must be sorted by length and the corresponding vertices
    appear in that order in the net, after all the non-colliding vertices
    =#
end