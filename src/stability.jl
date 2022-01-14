## Functions related to unstable nets (some vertices have the same equilibrium placement)

"""
    CollisionRootNode

A `CollisionRootNode` is a connected set of colliding vertices in a `CollisionSubNode`.
It is characterized by the subgraph induced by those of vertices.

The topological algorithm used can reliably identify a `CollisionRootNode` either if it
is a clique, or by the list of degrees of the vertices if it has 4 or less vertices.
"""
struct CollisionRootNode
    clique::Bool
    d::Vector{Int}
    vertices::Vector{Int} # sorted
    SimpleCollisionNode(vertices, clique, d) = new(vertices, clique, sort!(d))
end

function Base.isless(x::CollisionRootNode, y::CollisionRootNode)
    length(x.vertices) != length(y.vertices) && return length(x.vertices) < length(y.vertices)
    x.clique != y.clique && return y.clique
    return x.d < y.d
end
function Base.isequal(x::CollisionRootNode, y::CollisionRootNode)
    length(x.vertices) == length(y.vertices) && x.clique == y.clique && x.d == y.d
end


"""
    CollisionSubNode

In a `CollisionNode`, each group of vertices sharing the same neighbours out of the
`CollisionNode` are clustered in a `CollisionSubNode`.

The `CollisionSubNode` itself is split in `CollisionRootNode`, corresponding to its
connected components according to the subgraph restricted to the `CollisionSubNode`.
"""
struct CollisionSubNode
    rootnodes::Vector{CollisionRootNode}
    parent::Vector{Int}
    vertices::Vector{Int} # sorted
end


"""
    CollisionNode

On an equilibrium placement, if `n` vertices collide, they are clustered into a
`CollisionNode`. It is itself split into `CollisionSubNode`s containing vertices that share
the same neighbours out of the `CollisionNode`.
"""
struct CollisionNode
    subnodes::Vector{CollisionSubNode}
    vertrange::UnitRange{Int}
end


"""
    shrink_collisions(net::CrystalNet, collisions)

Remove all colliding vertices and replace them by one new vertex per `CollisionNode`, whose
neighbours are that of the vertices within.
"""
function shrink_collisions(net::CrystalNet{D,T}, collisions) where {D,T}
    first_colliding = first(first(collisions).vertrange)
    collision_vmap = collect(1:length(net.types))
    for (i, node) in enumerate(collisions)
        @toggleassert all(==(net.pos[first(node.vertrange)]), @view net.pos[node.vertrange])
        collision_vmap[node.vertrange] .= i
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
    append!(newtypes, net.pos[first(collisions[i].vertrange)] for i in 2:m)
    newtypes = net.types[1:first_colliding-1]
    append!(newtypes, :* for i in 1:m) # type of collision nodes will be Symbol('*')

    newnet = CrystalNet{D,T}(net.cell, newtypes, newpos, newgraph, net.options)

    return newnet, collision_vmap
end


"""
    expand_collision!(node::CollisionNode, g::PeriodicGraph, vmap)

Return the list of vertices in the `CollisionNode` sorted in a way that only depends on `g`

"""
function expand_collision!(node::CollisionNode, g::PeriodicGraph, ext_vmap, int_vmap)
    subnodes = nodes.subnodes
    for sub in subnodes
        for (i,p) in enumerate(sub.parent)
            sub.parent[i] = vmap[p]
        end
    end
    sort!(subnodes; by=(x->x.parent))

end


function expand_collisions!(collisions, graph::PeriodicGraph{D}, initial_graph, vmap, collision_vmap) where D
    n = nv(initial_graph)
    m = length(collisions)
    num_withcolliding = length(vmap)
    num_nocolliding = num_withcolliding - m

    @toggleassert num_nocolliding + sum(x->length(x.vertrange), collisions) == length(collision_vmap) == n

    rev_vmap = Vector{Int}(undef, num_withcolliding)
    for (i, j) in enumerate(vmap)
        rev_vmap[j] = i
    end


    prev_ranges = Vector{Int}(undef, m+1)
    prev_ranges[1] = num_nocolliding
    for (i,node) in enumerate(collisions)
        @toggleassert first(node.vertrange) == prev_ranges[i] + 1
        prev_ranges[i+1] = last(node.vertrange)
    end

    colliding_nodes = [rev_vmap[collision_vmap[prev_ranges[i]]] for i in 2:(m+1)]
    I = sortperm(colliding_nodes)
    colliding_nodes = colliding_nodes[I]
    collisions = collisions[I] # The order of collisions is now independent of the initial representation

    new_ranges = num_nocolliding .+ cumsum(length(node.vertrange) for node in collisions)
    pushfirst!(new_ranges, num_nocolliding)

    perm_nocolliding = Vector{Int}(undef, num_nocolliding)
    perm_withcolliding = Vector{Int}(undef, num_withcolliding)
    perm = Vector{Int}(undef, n)
    current_idx_colliding = 1
    for i in 1:nv(graph)
        if i == colliding_nodes[current_idx_colliding]
            perm_withcolliding[num_nocolliding+current_idx_colliding] = i
            newpos = (1+new_ranges[current_idx_colliding]):new_ranges[current_idx_colliding+1]
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
    rev_vmap = perm_withcolliding[rev_vmap]
    ext_vmap = perm[1:num_nocolliding]
    initial_graph = initial_graph[perm]

    # At this point, the graph has all the collision nodes stacked at its end, but kept in
    # the same order as that given by the topological genome of the graph.

    add_vertices!(newgraph, n - num_nocolliding)
    for i in 1:m
        rnge = (1+new_ranges[i]):(new_ranges[i+1])
        for x in neighbors(graph, i)
            for j in neighbors(initial_graph, x.v)
                if j.v ∈ rnge
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
    for (i, node) in collisions
        neworder = expand_collision!(node, newgraph)
        newperm[(1+new_ranges[i]):(new_ranges[i+1])] = (1+new_ranges[i]) .+ neworder
    end
    retgraph = newgraph[newperm]
    newperm[1:num_nocolliding] = vmap[perm_nocolliding]

    return newperm, retgraph
end


"""
    collision_nodes!(c::CrystalNet)

Check that the net is stable, i.e. that no two vertices have the same equilibrium placement.

A net is still considered stable if the collisions in equilibrium placement cannot lead to
different topological genomes. In practice, this happens when:
A) there is no edge between two collision sites and
B) there is no edge between a collision site and two representatives of the same vertex and
C) for each collision site, for all colliding vertices on the site sharing the exact same
   neighbours out of the site, either:
    - there are 4 or less vertices in the collision site, or
    - there is no edge between any pair of vertices in the collision site, or
    - there is an edge between each pair of vertices (the subgraph is a clique).

In this case, return the list of corresponding `CollisionNode`, the list being empty if the
net is truly stable. Otherwise, return `nothing`.

This function reorders the vertices in the net so that the vertices in a `CollisionNode`
appear after the others and are contiguous when they belong to the same `CollisionNode`.
"""
function collision_nodes!(net::CrystalNet{D}) where D
    @toggleassert issorted(net.pos)
    collision_sites = Vector{Int}[]
    _collision_vertices = Int[]
    flag_collision = false
    for i in 2:length(net.pos)
        if net.pos[i-1] == net.pos[i]
            if flag_collision
                push!(collision_sites[end], i)
                push!(_collision_vertices, i)
            else
                push!(collision_sites, [i-1, i])
                push!(_collision_vertices, i-1, i)
                flag_collision = true
            end
        else
            flag_collision = false
        end
    end

    ret = CollisionNode[]
    isempty(collision_sites) && return ret # no collision at all

    collision_vertices = BitSet(_collision_vertices)

    for site in collision_site
        adjacent_nodes = Vector{Int}[]
        handled = BitSet()
        for simplenode in site
            simplenode ∈ handled && continue
            push!(handled, simplenode)
            Q = [simplenode]
            ref_ext_neighbours = Set{PeriodicVertex{D}}()
            internal_degree = Int[]
            for u in Q
                ext_neighbours = PeriodicVertex{D}[]
                push!(internal_degree, 0)
                for x in neighbors(net.graph, u)
                    v = x.v
                    if iszero(x.ofs)
                        if v ∈ Q
                            internal_degree[end] += 1
                            continue
                        end
                        if v ∈ site
                            push!(Q, v)
                            push!(handled, v)
                            internal_degree[end] += 1
                            continue
                        end
                    end
                    if v ∈ collision_vertices
                        return nothing # edge between two collision sites
                    else
                        push!(ext_neighbours, x)
                    end
                end
            end
            @toggleassert !isempty(ext_neighbours)
            if isempty(ref_ext_neighbours)
                length(site) > 1 && (ref_ext_neighbours = Set(ext_neighbours))
            else
                issetequal(ext_neighbours, ref_ext_neighbours) || return nothing # all vertices do not share the same external neighbours
            end
            push!(adjacent_nodes, Q)
        end
    end


    #=
    each `vertices` field must be sorted
    The list of collision nodes must be sorted by length and the corresponding vertices
    appear in that order in the net, after all the non-colliding vertices
    =#
end