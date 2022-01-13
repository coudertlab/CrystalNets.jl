## Functions related to unstable nets (some vertices have the same equilibrium placement)

"""
    SimpleCollisionNode

A `SimpleCollisionNode` is a connected set of colliding vertices. It is characterized by
the subgraph induced by those of vertices.

The topological algorithm used can reliably identify a `SimpleCollisionNode` either if it
is a clique, or by the list of degrees of the vertices if it has 4 or less vertices.
"""
struct SimpleCollisionNode
    num::Int16
    full::Bool
    d::Vector{Int16}
    SimpleCollisionNode(num, full, d) = new(num, full, sort!(d))
end

function Base.isless(x::SimpleCollisionNode, y::SimpleCollisionNode)
    x.num != y.num && return x.num < y.num
    x.full != y.full && return y.full
    return x.d < y.d
end
function Base.isequal(x::SimpleCollisionNode, y::SimpleCollisionNode)
    x.num == y.num && x.full == y.full && x.d == y.d
end


"""
    CollisionNode

On an equilibrium placement, if `n` vertices collide, they are clustered into
a `CollisionNode`. It is itself split into its connected components when considering the
subgraph with only the edges between vertices in the `CollisionNode`. Each component is a
`SimpleCollisionNode`
"""
struct CollisionNode
    nodes::Vector{SimpleCollisionNode}
    CollisionNode(l) = new(sort!(l))
end


"""
    collision_nodes(c::CrystalNet)

Check that the net is stable, i.e. that no two vertices have the same equilibrium placement.

A net is still considered stable if the collisions in equilibrium placement cannot lead to
different topological genomes. This happens if there is no edge between two collision sites
and if, for each collision site, all colliding vertices share the exact same neighbours out
of the collision site and either:
- there are 4 or less vertices in the collision site, or
- there is no edge between any pair of vertices in the collision site, or
- there is an edge between each pair of vertices (the subgraph is a clique).

In this case, return the list of corresponding `CollisionNode`, the list being empty if the
net is truly stable. Otherwise, return `nothing`.
"""
function collision_nodes(net::CrystalNet{D}) where D
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
end