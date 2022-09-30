## Functions related to unstable nets (some vertices have the same equilibrium placement)

"""
    CollisionNode

Store the structure of a collision node through the subgraph `g` extracted with only the
edges bond to the vertices in the node.

The `num` field corresponds to the number of vertices in `g` that collide in the initial
graph, while the `neighs` field stores the indices of their neighbors out of the collision
site.

The colliding vertices are the `num` first vertices of `g`, the others are the neighbors.
In particular, `nv(g) == num + length(neighs)` and `g[(num+1):(nv(g))]` has no edge.
"""
struct CollisionNode
    g::PeriodicGraph{0}
    num::Int
    neighs::Vector{Int}
end
(Base.:(==))(c1::CollisionNode, c2::CollisionNode) = c1.num == c2.num && c1.g == c2.g && c1.neighs == c2.neighs
Base.hash(c::CollisionNode, h::UInt) = hash(c.g, hash(c.num, hash(c.neighs, h)))
Base.length(c::CollisionNode) = c.num

"""
    shrink_collisions(net::CrystalNet, collisions)

Remove all colliding vertices and replace them by one new vertex per `CollisionNode`, whose
neighbours are that of the vertices within.
"""
function shrink_collisions(net::CrystalNet{D,T}, collision_ranges::Vector{UnitRange{Int}}) where {D,T}
    first_colliding = first(first(collision_ranges))
    m = length(collision_ranges)
    n = first_colliding + m - 1

    collision_vmap = collect(1:length(net.types))
    for (i, node) in enumerate(collision_ranges)
        @toggleassert all(==(net.pge.pos[first(node)]), @view net.pge.pos[node])
        collision_vmap[node] .= first_colliding + i - 1
    end
    edgs = PeriodicEdge{D}[]
    for e in edges(net.pge.g)
        src = collision_vmap[e.src]
        dst = collision_vmap[e.dst.v]
        if src < first_colliding || dst < first_colliding
            push!(edgs, PeriodicEdge{D}(src, dst, e.dst.ofs))
        end
    end
    newgraph = PeriodicGraph{D}(edgs)

    @toggleassert nv(newgraph) == n

    newpos = net.pge.pos[1:first_colliding]
    append!(newpos, net.pge.pos[first(collision_ranges[i])] for i in 2:m)
    newtypes = net.types[1:first_colliding-1]
    append!(newtypes, :* for i in 1:m) # type of collision nodes will be Symbol('*')

    return CrystalNet{D,T}(net.pge.cell, newtypes, newpos, newgraph, net.options)
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
        neighls = neighbors(subgraph, xl)
        if length(neighls) == 1 && neighls[1] == xi
            xi, xj = xj, xi
        end
    end
    return n2 == 1 ? [xi, xj, only(subnodes[2]), only(subnodes[3])] :
           n2 == 2 ? [only(subnodes[1]), xi, xj, only(subnodes[3])] :
                     [only(subnodes[1]), only(subnodes[2]), xi, xj]
end

"""
    order_collision(graph::PeriodicGraph, colliding)

Given collision nodes (in the form of the corresponding list of colliding vertices), find
an ordering of them which is independent of the current ordering of these vertices and of
vertices which are neither in the collision node nor any of its neighbours.
Return an empty list if it fails.

This function assumes that no vertex in the node has a neighbour in another collision node
and that there are no two representatives of the same vertex that are neighbour to some
vertices in the range.
"""
function order_collision(graph::PeriodicGraph{D}, colliding) where D
    m = length(colliding)
    neigh_per_v = [Int[] for _ in 1:m] # for each colliding vertex, its non-colliding neighbors
    # subgraph is the graph extracted by only keeping the bonds between the colliding vertices
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
    # subnodes is a list of sublists of colliding vertex indices, such that all indices in
    # a sublist correspond to vertices sharing the same non-colliding neighbors.
    # subnodes is sorted by the lexicographical order of the corersponding lists of neighbors.
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
    unsorted_node(graph::PeriodicGraph, node::UnitRange{Int})

Create a CollisionNode corresponding to the input `node` in `graph`.

Colliding vertices are put first in the CollisionNode (see the definition of a
[`CollisionNode`](@ref)) but otherwise the order of vertices is kept the same as in the
initial graph.
"""
function unsorted_node(graph::PeriodicGraph, node::UnitRange{Int})
    neighbors_of_node_set = BitSet()
    for u in node
        union!(neighbors_of_node_set, x.v for x in neighbors(graph, u))
    end
    setdiff!(neighbors_of_node_set, node)
    neighbors_of_node = sort!(collect(neighbors_of_node_set))

    # newmap[u] is the index of vertex u in the returned CollisionNode
    newmap = Vector{Int}(undef, nv(graph))
    l = length(node)
    for (i, u) in enumerate(node)
        newmap[u] = i
    end
    for (i, u) in enumerate(neighbors_of_node)
        newmap[u] = l + i
    end

    edgs = PeriodicEdge{0}[]
    for u in node
        v = newmap[u]
        append!(edgs, PeriodicEdge{0}(v, newmap[x.v], ()) for x in neighbors(graph, u))
    end

    return CollisionNode(PeriodicGraph{0}(edgs), length(node), neighbors_of_node)
end

"""
    CollisionNode(c::CollisionNode, rev_vmap)

Return a collision node representing the same colliding vertices than in `c`, but ordered
according the order of their neighbors, given by `rev_vmap`.
In other words, if `rev_vmap[i] == 1` and vertex `i` is a neighbor of vertices in `c`,
those vertices will be ordered first (compared to other vertices in the same bracket).

If `rev_vmap === nothing`, assume `∀i, rev_vmap[i] == i`.
"""
function CollisionNode(c::CollisionNode, rev_vmap)
    sorted_neighbors = rev_vmap isa Vector{Int} ? rev_vmap[c.neighs] : c.neighs
    I = sortperm(sorted_neighbors)

    newmap = collect(1:nv(c.g))
    l = length(c)
    for (j, i) in enumerate(I)
        newmap[l+j] = l+i
    end
    sorted_neighbors = sorted_neighbors[I]

    newgraph = c.g[newmap]
    n = nv(newgraph)
    @toggleassert n == length(sorted_neighbors) + l
    neworder = order_collision(newgraph, 1:l)
    @toggleassert !isempty(neworder)
    perm = collect(1:nv(newgraph))
    perm[1:l] = neworder

    return CollisionNode(newgraph[perm], l, sorted_neighbors)
end

"""
    CollisionNode(graph::PeriodicGraph, node::UnitRange{Int}, vmap=nothing)

Return a canonical `CollisionNode` identifying the collision node.
This depends on the order of its neighbours, which is given by their order in the graph or
in `vmap` if provided.
"""
function CollisionNode(graph::PeriodicGraph, node::UnitRange{Int}, rev_vmap=nothing)
    return CollisionNode(unsorted_node(graph, node), rev_vmap)
end


"""
    collect_collisions(net::CrystalNet)

Return a tuple `(collision_sites, collision_vertices)` where
- `collision_sites` is the list a ranges of vertices that share a same position in `net`.
- `collision_vertices` is the list of all such vertices.
"""
function collect_collisions(net::CrystalNet)
    issorted(net.pge.pos) || error("Internal error: unsorted positions while collecting collisions")
    positions = net.pge.pos
    @toggleassert all(pos -> all(x -> 0 ≤ x < 1, pos), positions)
    collision_sites = UnitRange{Int}[] # each is a list of vertices at the same position
    collision_vertices = Int[] # list of vertices on a collision site
    start_collision = 0 # if not 0, first index of a node in a collision
    n = length(positions)
    for i in 2:n
        if positions[i-1] == positions[i]
            if start_collision != 0
                push!(collision_vertices, i)
            else
                push!(collision_vertices, i-1, i)
                start_collision = i-1
            end
        elseif start_collision != 0
            push!(collision_sites, start_collision:(i-1))
            start_collision = 0
        end
    end
    start_collision == 0 || push!(collision_sites, start_collision:n)
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
   α) the site is made of at most 4 vertices, or
   β) no 2 vertices on the site share the same exact set of neighbours out of the site.

In this case, return the `CollisionNodeList` with the corresponding `CollisionNode`s, the
list being empty if the net is truly stable. Otherwise, return `nothing`.

Also return an updated net where the vertices in a `CollisionNode` are collapsed into a new
vertex, appearing after the non-colliding vertices.
"""
function collision_nodes(net::CrystalNet{D,T}) where {D,T}
    collision_sites, collision_vertices = collect_collisions(net)

    isempty(collision_sites) && return net, CollisionNode[] # no collision at all

    # Check that conditions A, B and C are respected
    collision_vertices_set = BitSet(collision_vertices)
    for site in collision_sites
        known_neighbors = Dict{Int,SVector{D,Int}}()
        known_nlist = Set{Vector{PeriodicVertex{D}}}()
        uniquenlist = true # set while condition C)β) holds
        for u in site
            this_nlist = PeriodicVertex{D}[]
            for x in neighbors(net.pge.g, u)
                get!(known_neighbors, x.v, x.ofs) == x.ofs || return net, nothing # condition B
                x.v ∈ site && continue
                x.v ∈ collision_vertices_set && return net, nothing # condition A
                uniquenlist && push!(this_nlist, x)
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
    collision_ranges = Vector{UnitRange{Int}}(undef, length(collision_sites))
    for (k, node) in enumerate(collision_sites)
        collision_ranges[k] = idx:(idx+length(node)-1)
        for x in node
            vmap[idx] = x
            idx += 1
        end
    end
    @toggleassert idx == n+1

    newgraph = net.pge.g[vmap]
    newpos = net.pge.pos[vmap]
    ofs = newpos[1]
    offsets = Vector{SVector{D,Int}}(undef, n)
    for (i, x) in enumerate(newpos)
        offsets[i] = floor.(Int, x .- ofs)
        newpos[i] = x .- ofs .- offsets[i]
    end
    offset_representatives!(newgraph, .-offsets)
    # net.pge.pos[1] should always be [0,0,0]

    @toggleassert iszero(first(net.pge.pos))

    newnet = CrystalNet{D,T}(net.pge.cell, net.types[vmap], newpos, newgraph, net.options)
    newnodes = [CollisionNode(newnet.pge.g, node) for node in collision_ranges]
    return shrink_collisions(newnet, collision_ranges), newnodes
end




function collision_utils(collisions::Vector{CollisionNode}, num_withcolliding::Int, rev_vmap=nothing)
    m = length(collisions)
    num_nocolliding = num_withcolliding - m
    n = num_nocolliding + sum(length, collisions)
    new_collisions = collisions

    if rev_vmap isa Vector{Int}
        colliding_nodes = rev_vmap[(num_nocolliding+1):num_withcolliding]
        new_collisions = [CollisionNode(collisions[i], rev_vmap) for i in sortperm(colliding_nodes)]
        # The order of collisions and is made independent of the initial representation
    end

    return n, m, num_withcolliding, num_nocolliding, new_collisions
end

function collision_utils(collisions::Vector{CollisionNode}, vmap::Vector{Int})
    num_withcolliding = length(vmap)
    rev_vmap = Vector{Int}(undef, num_withcolliding)
    for (i, j) in enumerate(vmap)
        rev_vmap[j] = i
    end
    return collision_utils(collisions, num_withcolliding, rev_vmap)
end

"""
    expand_collisions(collisions::Vector{CollisionNode}, graph::PeriodicGraph, vmap)

Expand each collision node into the appropriate number of vertices so that the resulting
graph is isomorphic to the initial one, in a manner that only depends on the current graph.
Return the resulting graph.

`vmap` is the map of vertices between `initial_graph` (with collapsed collision nodes) and
`graph`
"""
function expand_collisions(collisions::Vector{CollisionNode}, graph::PeriodicGraph{D}, vmap) where D
    n, m, num_withcolliding, num_nocolliding, collisions = collision_utils(collisions, vmap)
    # newtypes = Symbol[j > num_nocolliding ? :O : :C for j in vmap]
    # export_vtf("/tmp/cgraph.vtf", CrystalNet3D(Cell(), newtypes, graph, Options()), 5)

    # collisions are now sorted according to vmap, both in the list and for each subgraph

    @toggleassert nv(graph) == num_withcolliding
    @toggleassert length(vmap) == nv(graph)

    # We now push all the collision nodes to the end
    perm = Vector{Int}(undef, num_withcolliding)
    current_idx_colliding = 1
    for (i, j) in enumerate(vmap)
        if j > num_nocolliding # vertex i in the current graph is a collision node
            perm[num_nocolliding + current_idx_colliding] = i
            current_idx_colliding += 1
        else
            idx = i + 1 - current_idx_colliding
            perm[idx] = i
        end
    end
    graph = graph[perm]
    newgraph = graph[1:num_nocolliding] # subgraph of the non-colliding nodes keeping the order
    newvmap = Vector{Int}(undef, n)
    newvmap[1:num_nocolliding] = vmap[1:num_nocolliding]

    # At this point, the graph has all the collision nodes stacked at its end, but kept in
    # the same order as that given by the topological genome of the graph.

    add_vertices!(newgraph, n - num_nocolliding)
    offsetcounter = num_nocolliding
    for (i, node) in enumerate(collisions)
        k = length(node)
        # bonds internal to the collision site
        for j in 1:k
            for x in neighbors(node.g, j)
                x.v ≤ k || continue
                add_edge!(newgraph, offsetcounter + j, PeriodicVertex{D}(offsetcounter + x.v))
            end
        end
        # other bonds
        for (j, x) in enumerate(neighbors(graph, num_nocolliding + i))
            # x is a (non-colliding) neighbour of the collision node
            for y in neighbors(node.g, length(node) + j)
                add_edge!(newgraph, offsetcounter + y.v, x)
            end
        end
        offsetcounter += length(node)
    end

    # newgraph contains all vertices from initial_graph
    # The non-colliding ones are ordered like in graph, the others are stacked in the order
    # of their collision node in graph, each colliding subgraph sorted with order_collision.

    return newgraph
end
