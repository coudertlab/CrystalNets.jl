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
CollisionNode(n) = CollisionNode(PeriodicGraph{0}(), n, Int[]) # fake collision node

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

struct EmptyList{T} <: AbstractVector{T} end
Base.size(::EmptyList) = (0,)

"""
    CollisionList <: AbstractVector{CollisionNode}

List of [`CollisionNode`](@ref).

Also contains a `priorities` field that give a priority to each non-colliding vertex of
the net
"""
struct CollisionList <: AbstractVector{CollisionNode}
    list::Vector{CollisionNode}
end

Base.size(cl::CollisionList) = (length(cl.list),)
Base.getindex(cl::CollisionList, i::Int) = cl.list[i]

function CollisionList(g::PeriodicGraph, collision_ranges::Vector{UnitRange{Int}})
    return CollisionList([CollisionNode(unsorted_node(g, node), nothing) for node in collision_ranges])
end

function CollisionList(collisions::CollisionList, vmap::Vector{Int}, kept_collisions::Vector{Int})
    return CollisionList([CollisionNode(collisions.list[i], vmap) for i in kept_collisions])
end

"""
    shrink_collisions(net::CrystalNet, collision_ranges)

Remove all colliding vertices and replace them by one new vertex per collision range, whose
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
        # if src < first_colliding || dst < first_colliding
        if src != dst || !iszero(e.dst.ofs) # TODO: check that this is correct
            push!(edgs, PeriodicEdge{D}(src, dst, e.dst.ofs))
        end
    end
    newgraph = PeriodicGraph{D}(edgs)

    @toggleassert nv(newgraph) == n

    newpos = net.pge.pos[1:first_colliding]
    append!(newpos, net.pge.pos[first(collision_ranges[i])] for i in 2:m)
    newtypes = net.types[1:first_colliding-1]
    append!(newtypes, :* for i in 1:m) # type of collision nodes will be Symbol("*")
    opts = permute_mapping!(net.options, collision_vmap)

    return CrystalNet{D,T}(net.pge.cell, newtypes, newpos, newgraph, opts)
end

macro makepriority(n)
    esc(:(return (J, $n%Int8 .+ degs)))
end
macro makepriority(n, m)
    esc(:(return (J, $n%Int8 .+ (degs .== $m%Int8))))
end

"""
    get_priority(subgraph::SimpleGraph)

Return `priority`, a `Vector{Int}` such that all non-colliding neighbors of vertex `i` in
`subgraph` (before permutation) are attributed `priority[i]`.
Also return a permutation of vertices of `subgraph` that makes priority non-decreasing.
"""
function get_priority(subgraph::SimpleGraph)
    m = nv(subgraph)
    degs = Int8.(degree(subgraph))
    J = sortperm(degs)
    m == 2 && @makepriority 0 # 00 or 11: 0-1
    sortdeg = degs[J]
    if m == 3
        first(sortdeg) == last(sortdeg) && @makepriority 2 2 # 000 or 222: 2-3
        @makepriority (3 + sortdeg[1]) # 011 or 112: 3-6
    end
    @toggleassert m == 4
    first(sortdeg) == last(sortdeg) && @makepriority 7 # 0000 or 1111 or 2222 or 3333: 7-10
    last(sortdeg) == 1 && @makepriority 11 # 0011: 11-12
    if last(sortdeg) == 2
        sortdeg[3] == 1 && @makepriority 13 # 0112: 13-15
        sortdeg[1] == 0 && @makepriority 16 2 # 0222: 16-17
        @makepriority 17 # 1122: 18-19
    end
    @toggleassert last(sortdeg) == 3
    sortdeg[2] == 1 && @makepriority 20 3 # 1113: 20-21
    sortdeg[1] == 1 && @makepriority 21 # 1223: 22-24
    @toggleassert sortdeg[1] == 2
    @makepriority 23 # 2233: 25-26
end

"""
    _order_collision!((priority, perm)::Tuple{Vector{Int8},Vector{Int}}, subnodes::Vector{Vector{Int}})

Modify `perm` to contain a permutation of the vertices of `subgraph` uniquely determined by
the connectivity of the vertex within `subgraph` first, and by its neighbors next.

`subgraph` should be the graph of bonds between vertices of the colliding node. It is not
explicitly passed as argument, but should satisfy `perm, priority == get_priority(subgraph)`
`subnodes[i]` is a list of colliding vertex indices sharing the same set of non-colliding
neighbors. `subnodes` itself is sorted, see `_extract_subnodes`.
"""
function _order_collision!((perm, priority)::Tuple{Vector{Int},Vector{Int8}}, subnodes::Vector{Vector{Int}})
    m = length(priority)
    equivalents = UnitRange{Int}[]
    last_priority = priority[perm[1]]
    new = true
    for i in 2:m
        this_priority = priority[perm[i]]
        if this_priority == last_priority
            if new
                push!(equivalents, (i-1):i)
                new = false
            else
                equivalents[end] = (first(equivalents[end])):i
            end
        else
            new = true
        end
    end
    (isempty(equivalents) || length(subnodes) == 1) && return nothing

    neighbor_priority = Vector{Int}(undef, m)
    for (k, node) in enumerate(subnodes)
        for j in node
            neighbor_priority[j] = k
        end
    end
    for eq in equivalents
        sort!(@view(perm[eq]); by=Base.Fix1(getindex, neighbor_priority))
    end
    nothing
end


function _extract_subnodes(neigh_per_v::Vector{Vector{Int}})
    I = sortperm(neigh_per_v; by=x->(length(x), x))
    # subnodes is a list of sublists of colliding vertex indices, such that all indices in
    # a sublist correspond to vertices sharing the same non-colliding neighbors.
    # subnodes is itself sorted by the lexicographical order of the lists of neighbors.
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
    subnodes
end

"""
    order_collision(graph::PeriodicGraph, colliding)

Given collision nodes (in the form of the corresponding list of colliding vertices), find
an ordering of them which is independent of the current ordering of these vertices and of
vertices which are neither in the collision node nor any of its neighbours.
Return a this ordering and a priority list for each colliding vertex, or two empty lists if
it fails.

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
    subnodes = _extract_subnodes(neigh_per_v)
    m > 4 && return Int[], Int8[]
    priorities = get_priority(subgraph)
    _order_collision!(priorities, subnodes)
    return priorities
end

function order_collision(c::CollisionNode)
    m = c.num
    subgraph = SimpleGraph(collect(Edge(u, v) for (u, (v, _)) in edges(c.g) if u ≤ m && v ≤ m))
    add_vertices!(subgraph, m - nv(subgraph))
    neigh_per_v = [Int[] for _ in 1:m]
    for (i, nlist) in enumerate(neigh_per_v)
        for (x, _) in neighbors(c.g, i)
            x > m && push!(nlist, x)
        end
    end
    subnodes = _extract_subnodes(neigh_per_v)
    m > 4 && return Int[], Int8[]
    priorities = get_priority(subgraph)
    _order_collision!(priorities, subnodes)
    return priorities
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
    nv(c.g) == 0 && return c
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
    neworder = first(order_collision(newgraph, 1:l))
    @toggleassert !isempty(neworder)
    perm = collect(1:nv(newgraph))
    perm[1:l] = neworder

    return CollisionNode(newgraph[perm], l, sorted_neighbors)
end

# """
#     CollisionNode(graph::PeriodicGraph, node::UnitRange{Int}, vmap=nothing)

# Return a canonical `CollisionNode` identifying the collision node.
# This depends on the order of its neighbours, which is given by their order in the graph or
# in `vmap` if provided.
# """
# function CollisionNode(graph::PeriodicGraph, node::UnitRange{Int}, rev_vmap=nothing)
#     return CollisionNode(unsorted_node(graph, node), rev_vmap)
# end


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
C) for each collision site, the site is made of at most 4 vertices

In this case, return the `CollisionNodeList` with the corresponding `CollisionNode`s, the
list being empty if the net is truly stable. Otherwise, return the list of ranges of
colliding vertices.

Also return an updated net where the vertices in a `CollisionNode` are collapsed into a new
vertex, appearing after the non-colliding vertices.

Finally, return a net equivalent to the initial one, with its vertices rotated so that
colliding ones are contiguous, and after non-colliding ones. The "list of ranges" mentioned
before refers to the vertices of that net.
"""
function collision_nodes(net::CrystalNet{D,T}) where {D,T}
    collision_sites, collision_vertices = collect_collisions(net)

    isempty(collision_sites) && return CollisionList(net.pge.g, UnitRange{Int}[]), net, net # no collision at all

    # Check that conditions A, B and C are respected
    collision_vertices_set = BitSet(collision_vertices)
    unstable = false
    for site in collision_sites
        if length(site) > 4
            unstable = true # condition C
        end
        known_neighbors = Dict{Int,SVector{D,Int}}()
        for u in site
            for x in neighbors(net.pge.g, u)
                if get!(known_neighbors, x.v, x.ofs) != x.ofs
                    unstable = true # condition B
                end
                x.v ∈ site && continue
                if !unstable && x.v ∈ collision_vertices_set
                    unstable = true # condition A
                end
            end
        end
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

    opts = rev_permute_mapping!(net.options, vmap, length(net.types))
    newnet = CrystalNet{D,T}(net.pge.cell, net.types[vmap], newpos, newgraph, opts)
    # newnodes = [CollisionNode(newnet.pge.g, node) for node in collision_ranges]
    newnodes = unstable ? collision_ranges : CollisionList(newnet.pge.g, collision_ranges)
    return newnodes, shrink_collisions(newnet, collision_ranges), newnet
end


"""
    collision_utils(collisions::CollisionList, vmap::Vector{Int})

Return:
- n: number of vertices in the full graph (with expanded collision nodes).
- num_withcolliding: number of vertices, including vertices standing for collision nodes.
- num_nocolliding: number of non-colliding vertices.
- new_collisions: a [`CollisionList`](@ref) corresponding to `collisions` in the graph
  resulting from applying `vmap` to its vertices. This new list contains collision nodes in
  the same order as they appear in the graph resulting from `vmap`.
"""
function collision_utils(collisions::CollisionList, vmap::Vector{Int})
    num_withcolliding = length(vmap)
    rev_vmap = invperm(vmap)
    m = length(collisions)
    num_nocolliding = num_withcolliding - m
    n = num_nocolliding + sum(length, collisions)

    colliding_nodes = rev_vmap[(num_nocolliding+1):num_withcolliding]
    new_collisions = [CollisionNode(collisions[i], rev_vmap) for i in sortperm(colliding_nodes)]
    # The order of collisions and is made independent of the initial representation

    return n, num_withcolliding, num_nocolliding, new_collisions
end

"""
    expand_collisions(collisions::Vector{CollisionNode}, graph::PeriodicGraph, vmap)

Expand each collision node into the appropriate number of vertices so that the resulting
graph is isomorphic to the initial one, in a manner that only depends on the current graph.
Return the resulting graph.

`vmap` is the map of vertices between `initial_graph` (with collapsed collision nodes) and
`graph`

Also return the permutations of nodes of `graph` prior to expansion.
"""
function expand_collisions(collisions::CollisionList, graph::PeriodicGraph{D}, vmap) where D
    n, num_withcolliding, num_nocolliding, collisions = collision_utils(collisions, vmap)
    # newtypes = Symbol[j > num_nocolliding ? :O : :C for j in vmap]
    # export_vtf("/tmp/cgraph.vtf", CrystalNet3D(Cell(), newtypes, graph, Options()), 5)

    # collisions are now sorted according to vmap, both in the list and for each subgraph

    @toggleassert length(vmap) == nv(graph) == num_withcolliding

    # We now push all the collision nodes to the end
    perm = Vector{Int}(undef, num_withcolliding)
    current_idx_colliding = 1
    for (i, j) in enumerate(vmap)
        if j > num_nocolliding # vertex i in the current graph is a collision node
            perm[num_nocolliding + current_idx_colliding] = i
            current_idx_colliding += 1
        else
            perm[i+1-current_idx_colliding] = i
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

    return newgraph, perm
end

"""
    build_collision_ranges(collisions::CollisionList, n_shrunk)

Given the `collisions` list and the number `n_shrunk` of nodes in the shrunk net, return
the list of ranges of indices such that the `i`-th such range `r` contains the indices of
the colliding vertices in the expanded net that form the `i`-th collision node.
"""
function build_collision_ranges(collisions::CollisionList, n_shrunk)
    collision_ranges = Vector{UnitRange{Int}}(undef, length(collisions))
    counter = n_shrunk - length(collisions)
    for (i, collision_node) in enumerate(collisions)
        next_counter = counter + length(collision_node)
        collision_ranges[i] = (counter+1):next_counter
        counter = next_counter
    end
    collision_ranges
end


function _edge_split(x)
    if x isa Tuple
        x
    else
        s, (d, ofs) = x
        s, d, ofs
    end
end

function _edge_make(x, y, t::SVector{D,T}) where {D,T}
    T <: Rational ? (x, y, t) : PeriodicEdge{D}(x, y, t)
end

# swaps vertices a and a+1 in the list of edges
function vertex_swap!(edgs::Vector{T}, a) where T
    i = 1
    m = length(edgs)
    while i ≤ m
        s, d, ofs = _edge_split(edgs[i])
        if s == a
            p = i
            sp = s; dp = d; ofsp = ofs
            while true
                edgs[p] = _edge_make(a+1, dp - (dp==(a+1)) + (dp==a), ofsp)
                p += 1
                p > m && break
                sp, dp, ofsp = _edge_split(edgs[p])
                sp != a && break
            end
            q = p
            while sp == a+1
                edgs[p] = _edge_make(a, dp - (dp==(a+1)) + (dp==a), ofsp)
                p += 1
                p > m && break
                sp, dp, ofsp = _edge_split(edgs[p])
            end
            if p != q
                reverse!(edgs, i, q-1)
                reverse!(edgs, q, p-1)
                reverse!(edgs, i, p-1)
            end
            i = p
        elseif s == a+1
            edgs[i] = _edge_make(a, d - (d==(a+1)) + (d==a), ofs)
            i += 1
        elseif d == a
            k = i
            sk = s; dk = d; ofsk = ofs
            while true
                edgs[k] = _edge_make(s, a+1, ofsk)
                k += 1
                k > m && break
                sk, dk, ofsk = _edge_split(edgs[k])
                (sk != s || dk != a) && break
            end
            j = k
            while sk == s && dk == a+1
                edgs[k] = _edge_make(s, a, ofsk)
                k += 1
                k > m && break
                sk, dk, ofsk = _edge_split(edgs[k])
            end
            if j != k
                reverse!(edgs, i, j-1)
                reverse!(edgs, j, k-1)
                reverse!(edgs, i, k-1)
            end
            i = k
        elseif d == a+1
            edgs[i] = _edge_make(s, a, ofs)
            i += 1
        else
            i += 1
        end
    end
    nothing
end


function candidate_key_unstable(net::CrystalNet{D,T}, shrunk_candidate, u_s, basis, collision_ranges) where {D,T}
    n = nv(net.pge.g)
    newpos_s, offsets_s, vmap_s = shrunk_candidate

    newpos = Vector{SVector{D,T}}(undef, n) # positions of the kept representatives
    offsets = Vector{SVector{D,Int32}}(undef, n) # offsets of the new representatives w.r.t. the original one, in the original basis
    vmap = Vector{Int}(undef, n) # bijection from the old to the new node number
    rev_vmap = Vector{Int}(undef, n)
    i = 1
    first_collision_m1 = first(first(collision_ranges)) - 1
    for (j, v) in enumerate(vmap_s)
        if v ≤ first_collision_m1
            vmap[i] = v
            rev_vmap[v] = i
            newpos[i] = newpos_s[j]
            offsets[i] = offsets_s[j]
            i += 1
        else
            rnge = collision_ranges[v - first_collision_m1]
            newim1 = i + length(rnge) - 1
            vmap[i:newim1] .= rnge
            newpos[i:newim1] .= (newpos_s[j],)
            offsets[i:newim1] .= (offsets_s[j],)
            rev_vmap[rnge] .= i:newim1
            i = newim1 + 1
        end
    end

    u = u_s ≤ first_collision_m1 ? u_s : first(collision_ranges[u_s - first_collision_m1])
    origin = net.pge.pos[u]
    edgs = Tuple{Int,Int,SVector{D,T}}[]
    bigbasis = T == Rational{BigInt} ? basis : widen(T).(basis)
    mat = T == Rational{BigInt} ? inv(bigbasis) : T.(inv(bigbasis))
    for t in 1:n # t is the node being processed
        neighs = neighbors(net.pge.g, vmap[t])
        ofst = offsets[t]
        for x in neighs
            coordinate = mat*(net.pge.pos[x.v] .+ x.ofs .- origin .+ ofst)
            idx = rev_vmap[x.v]
            realofs = coordinate .- newpos[idx]
            push!(edgs, (t, idx, realofs))
        end
    end
    @toggleassert allunique(edgs)
    return vmap, edgs
end

function topological_key_unstable(net::CrystalNet{D,T}, collision_ranges, shrunk_net, candidates) where {D,T}
    @assert issorted(first.(collision_ranges))
    n = length(net.pge)
    @assert n == last(last(collision_ranges))
    first_collision = first(first(collision_ranges))
    @assert sum(length.(collision_ranges)) == n - first_collision + 1
    permutations_d = [[i] for i in 1:(first_collision-1)]
    for cr in collision_ranges
        append!(permutations_d, collect(cr) for _ in 1:length(cr))
    end

    dummy_edges = [(n + 1, 0, zero(SVector{D,T}))]
    minimal_edgs = copy(dummy_edges)
    vmap = collect(1:length(net.pge))
    # TODO: add constraints on equivalent vertices

    force_consecutive = Dict(first(rnge) => length(rnge) for rnge in collision_ranges)
    for rnge in collision_ranges
        for x in Iterators.drop(rnge, 1)
            force_consecutive[x] = -first(rnge)
        end
    end

    for (v, basis) in candidates
        shrunk_candidate = candidate_key(shrunk_net, v, basis, dummy_edges, Val(true))
        newvmap, edgs = candidate_key_unstable(net, shrunk_candidate, v, basis, collision_ranges)
        # @show newvmap
        @toggleassert !isempty(newvmap)

        sort!(edgs)
        if edgs < minimal_edgs
            minimal_edgs = copy(edgs)
            vmap = copy(newvmap)
            # minimal_basis = basis
        end

        new_positions = [(j = findfirst(==(first(range)), newvmap); j:(j+length(range)-1)) for range in collision_ranges]
        @toggleassert [newvmap[x] for x in new_positions] == collision_ranges
        swaps = ContiguousPlainChangesIterator(new_positions)

        # display(edgs)
        for swap in swaps
            # @show swap
            vertex_swap!(edgs, swap)
            sort!(edgs)
            # display(edgs)
            newvmap[swap], newvmap[swap+1] = newvmap[swap+1], newvmap[swap]
            if edgs < minimal_edgs
                minimal_edgs = copy(edgs)
                vmap = copy(newvmap)
                # minimal_basis = basis
            end
        end
        yield()
    end

    # rev_vmap = Vector{Int}(undef, n)
    # for (i, j) in enumerate(minimal_vmap)
    #     rev_vmap[j] = i
    # end

    newbasis, newedges = findbasis(minimal_edgs)
    graph = PeriodicGraph{D}(n, newedges)
    # return graph, collisions, minimal_vmaps

    # vmap = first(minimal_vmaps)
    @toggleassert let mapped = net.pge.g[vmap]; degree(graph) == degree(mapped) && quotient_graph(graph) == quotient_graph(mapped) end

    # tmpnet = CrystalNet{D,T}(PeriodicGraphEmbedding{D,T}(graph, net.pge.pos[minimal_vmap], net.pge.cell), net.types[minimal_vmap], net.options)
    # export_vtf("/tmp/tmpnet.vtf", tmpnet, 3)

    if !isnothing(net.options.track_mapping)
        map = rev_permute_mapping!(net.options, vmap).track_mapping
        if !net.options.keep_single_track
            _clust = first(net.options.clusterings)
            println("Mapping for ", length(net.options.clusterings) == 1 ? _clust : net.options.clusterings, map)
        end
    end

    return graph
end
