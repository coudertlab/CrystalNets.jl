## Functions related to unstable nets (some vertices have the same equilibrium placement)


struct CollisionNode
    rnge::UnitRange{Int}
    uniquecsequences::Vector{Vector{Int}} # unique coordination sequences
    uniques::Vector{Int} # position of the first unique sequences in the sorted list of coordination sequences
    subranges::Vector{UnitRange{Int}} # indices of rnge; the first subrange starts at 1
end

function get_CollisionNode(g::PeriodicGraph, node::UnitRange)
    csequences = [coordination_sequence(g, i, 6) for i in node]
    I = sortperm(csequences)
    csequences = csequences[I]
    subranges = UnitRange{Int}[]
    next_start = 1
    last_csequence = csequences[1]
    uniques = Int[1]
    for i in 2:length(csequences)
        csequence = csequences[i]
        if csequence != last_csequence
            last_csequence = csequence
            push!(subranges, next_start:(i-1))
            next_start = i
            push!(uniques, i)
        end
    end
    keepat!(csequences, uniques)
    push!(subranges, next_start:length(node))
    filter!(>(1)∘length, subranges)
    # TODO: remove subranges containing only equivalent nodes
    CollisionNode(node, csequences, uniques, subranges), I
end

function possibly_equivalent_nodes(node1::CollisionNode, node2::CollisionNode)
    length(node1.rnge) == length(node2.rnge) &&
        node1.uniquecsequences == node2.uniquecsequences &&
        node1.uniques == node2.uniques
        node1.subranges == node2.subranges
end

"""
    CollisionList <: AbstractVector{CollisionNode}

List of [`CollisionNode`](@ref).
"""
struct CollisionList <: AbstractVector{CollisionNode}
    list::Vector{CollisionNode}
end

Base.size(cl::CollisionList) = (length(cl.list),)
Base.getindex(cl::CollisionList, i::Int) = cl.list[i]

function get_CollisionList(g::PeriodicGraph, collision_ranges::Vector{UnitRange{Int}})
    vmap = collect(1:nv(g))
    list = Vector{CollisionNode}(undef, length(collision_ranges))
    for (i, node) in enumerate(collision_ranges)
        cnode, I = get_CollisionNode(g, node)
        list[i] = cnode
        vmap[node] .= (first(node)-1) .+ I
    end
    CollisionList(list), vmap
end
CollisionList() = CollisionList(CollisionNode[])

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

Check whether the net is stable, i.e. that no two vertices have the same equilibrium
placement.

Return `(shrunk_net, (equiv_net, collisions))` where:
- `equiv_net` is a net equivalent to the initial one, with its vertices rotated so that
  colliding ones are contiguous, and after non-colliding ones.
- `collisions` is a list of ranges corresponding to the indices of the vertices that
  collide to the same node in `equiv_net`.
  For example, `collisions == [4:5, 6:9]` indicates that the vertices 4 and 5 collide to a
  node, and the vertices 6 to 9 collide to another node.
  By definition of `equiv_net`, these ranges are always contiguous and the last range ends
  with the last vertex of the net. `collisions` is empty iff the net is stable.
- `shrunk_net` is the net where all the vertices that collide to a single node are
  collapsed into a single new vertex.
  For example, if `collisions = [4:5, 6:9]`, the `shrunk_net` will only have 5 vertices:
  vertices 1, 2 and 3 correspond vertices 1, 2 and 3 of the initial net, vertex 4
  corresponds to the colliding vertices 4 and 5 of the initial net, and vertex 5
  corresponds to the colliding vertices 6 to 9 of the initial net.
"""
function collision_nodes(net::CrystalNet{D,T}) where {D,T}
    collision_sites, collision_vertices = collect_collisions(net)

    isempty(collision_sites) && return net, (net, CollisionList()) # no collision at all

    collision_vertices_set = BitSet(collision_vertices)

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

    clist, cvmap = get_CollisionList(newgraph, collision_ranges)
    vmap = vmap[cvmap]
    opts = rev_permute_mapping!(net.options, vmap, length(net.types))
    newnet = CrystalNet{D,T}(net.pge.cell, net.types[vmap], newpos, newgraph[cvmap], opts)
    return shrink_collisions(newnet, collision_ranges), (newnet, clist)
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

function minute_collision_ranges(collisions)
    ret = UnitRange{Int}[]
    for node in collisions
        cofs = first(node.rnge) - 1
        append!(ret, rnge .+ cofs for rnge in node.subranges)
    end
    ret
end

function topological_key_unstable(net::CrystalNet{D,T}, collisions::CollisionList, shrunk_net, candidates) where {D,T}
    collision_ranges = [node.rnge for node in collisions]
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

    permutationranges = minute_collision_ranges(collisions)

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

        new_positions = [(j = findfirst(==(first(range)), newvmap); j:(j+length(range)-1)) for range in permutationranges]
        @toggleassert [newvmap[x] for x in new_positions] == permutationranges
        swaps = ContiguousPlainChangesIterator(new_positions, false)

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
