## Cell minimization and determination of the minimal net and crystal

"""
    possible_translations(c::Union{CrystalNet,CrystalNet})

Return a list of tuples `(nz, i_max_den, max_den, t)` where
- `t` is a translation mapping at the origin vertex to another one in the unit cell.
- `max_den` is the maximum denominator in the `D` coefficients of `t`.
- `i_max_den` is the index.
- `nz` is the number of zeros in `t`.

The list is guaranteed to contain all the possible valid translations but may contain some
invalid translations.

See also: [`find_all_valid_translations`](@ref), `PeriodicGraphEmbeddings.check_valid_symmetry`
"""
function possible_translations(c::CrystalNet{D,T}) where {D,T}
    ts = Tuple{Int, Int, Int, SVector{D,T}}[]
    sortedpos = copy(c.pge.pos)
    origin = popfirst!(sortedpos)
    @toggleassert iszero(origin)
    for t in sortedpos
        @toggleassert t == back_to_unit.(t)
        max_den, i_max_den = findmax(denominator.(t))
        numerator(t[i_max_den]) == 1 || continue
        nz = count(iszero, t)
        push!(ts, (nz, i_max_den, max_den, t))
    end
    return sort!(ts; by=(x->(x[1], x[2], x[3], norm(x[4]))))
end


"""
    find_all_valid_translations(c::Union{Crystal,CrystalNet{D}}, collisions::CollisionList) where D

Return a `D`-tuple of list of tuples `(i_max_den, max_den, t)` (see
[`possible_translations`](@ref) for interpretation) where the `n`-th list contains all
valid translations of the net having exactly `n-1` zeros.

A translation is valid if it maps exactly each vertex to a vertex and each edge to an edge.

See also: [`possible_translations`](@ref), `PeriodicGraphEmbeddings.check_valid_symmetry`
"""
function find_all_valid_translations(shrunk_net::CrystalNet{D,T}, collisions::CollisionList) where {D,T}
    ret = NTuple{D, Vector{Tuple{Int, Int, SVector{D,T}}}}(ntuple(_->[], Val(D)))
    check_symmetry = CheckSymmetryWithCollisions(collisions)
    for (nz, i_max_den, max_den, t) in possible_translations(shrunk_net)
        vmap = check_symmetry(shrunk_net.pge, t, nothing, shrunk_net.types)
        if vmap isa Vector{PeriodicVertex{D}}
            push!(ret[nz+1], (i_max_den, max_den, t))
        end
    end
    return ret
end

"""
    minimal_volume_matrix(translations::NTuple{D}) where {D}

Given the output of [`find_all_valid_translations`](@ref), compute the transformation that
allows reducing the net to its minimal cell.
"""
function minimal_volume_matrix end

function minimal_volume_matrix(translations::Tuple{Vector{Tuple{Int,Int,SVector{1,T}}}}) where T
    nz0 = translations[1]
    denmax = 1
    imax = 0
    for j in 1:length(nz0)
        _, den, _ = nz0[j]
        if den > denmax
            imax = j
            denmax = den
        end
    end
    if imax == 0
        push!(nz0, (1, 1, [1]))
        imax = length(nz0)
    end
    _nz0 = [nz0[imax]]
    empty!(nz0)
    append!(nz0, _nz0)

    n = length(nz0)
    detmin = one(T)
    best = 0
    @inbounds for i in 1:n
        d = nz0[i][3][1]
        d == 0 && continue
        if d < detmin
            detmin = d
            best = i
        end
    end
    @toggleassert !iszero(best)
    ret = hcat(nz0[best][3])
    if ret[1] < 0
        ret = hcat(.-nz0[best][3])
    end
    return ret
end

function minimal_volume_matrix(translations::NTuple{2, Vector{Tuple{Int,Int,SVector{2,T}}}}) where T
    nz0, nz1 = translations
    denmax = [1, 1]
    imax = [0, 0]
    for j in 1:length(nz1)
        i, den, _ = nz1[j]
        if den > denmax[i]
            imax[i] = j
            denmax[i] = den
        end
    end
    for j in 1:2
        if imax[j] == 0
            push!(nz1, (j, 1, [j==1, j==2]))
            imax[j] = length(nz1)
        end
    end
    _nz1 = [nz1[i] for i in imax]
    empty!(nz1)
    append!(nz1, _nz1)

    # TODO optimize
    all = vcat(nz0, nz1)
    n = length(all)
    detmin = one(T)
    best = (0, 0)
    @inbounds for i in 1:n-1
        for j in i+1:n
            d = abs(det(hcat(all[i][3], all[j][3])))
            d == 0 && continue
            if d < detmin
                detmin = d
                best = (i, j)
            end
        end
    end
    @toggleassert !any(iszero.(best))
    i, j = best
    ret = hcat(all[i][3], all[j][3])
    if det(ret) < 0
        ret = hcat(all[i][3], .-all[j][3])
    end
    return ret
end

function minimal_volume_matrix(translations::NTuple{3, Vector{Tuple{Int,Int,SVector{3,T}}}}) where T
    nz0, nz1, nz2 = translations

    denmax = [1, 1, 1]
    imax = [0, 0, 0]
    for j in 1:length(nz2)
        i, den, _ = nz2[j]
        if den > denmax[i]
            imax[i] = j
            denmax[i] = den
        end
    end
    for j in 1:3
        if imax[j] == 0
            push!(nz2, (j, 1, [j==1, j==2, j==3]))
            imax[j] = length(nz2)
        end
    end
    _nz2 = [nz2[i] for i in imax]
    empty!(nz2)
    append!(nz2, _nz2)

    # TODO optimize
    all = vcat(nz0, nz1, nz2)
    n = length(all)
    detmin = one(T)
    best = (0, 0, 0)
    @inbounds for i in 1:n-2
        for j in i+1:n-1
            for k in j+1:n
                d = abs(det(hcat(all[i][3], all[j][3], all[k][3])))
                d == 0 && continue
                if d < detmin
                    detmin = d
                    best = (i, j, k)
                end
            end
        end
    end
    @toggleassert !any(iszero.(best))
    i, j, k = best
    ret = hcat(all[i][3], all[j][3], all[k][3])
    if det(ret) < 0
        ret = hcat(all[i][3], all[j][3], .-all[k][3])
    end
    return ret
end


"""
    reduce_with_matrix(c::CrystalNet, mat, collisions)

Given the net and the output of `minimal_volume_matrix` computed on the valid translations
of the net, return the new net representing the initial net in the computed unit cell.
"""
function reduce_with_matrix(c::CrystalNet{D,Rational{T}}, mat, collisions) where {D,T}
    if D == 3
        cell = Cell(c.pge.cell, c.pge.cell.mat * mat)
    else
        _mat = SizedMatrix{3,3,BigFloat}(LinearAlgebra.I)
        _mat[1:D,1:D] .= mat
        cell = Cell(c.pge.cell, c.pge.cell.mat * _mat)
    end

    imat = T.(inv(mat)) # The inverse should only have integer coefficients
    poscol = (imat,) .* c.pge.pos # position in the new reference, will be in [0,1)
    n = length(poscol)
    offset = Vector{SVector{D,Int}}(undef, n) # offsets in the new reference
    for (i, pos) in enumerate(poscol)
        ofs = floor.(Int, pos)
        offset[i] = ofs
        poscol[i] = pos .- ofs
    end
    I_sort = sort(1:n; by=i->(poscol[i], hash_position(offset[i])))
    _i = popfirst!(I_sort)
    @toggleassert iszero(offset[_i])

    # I_kept is the list of index of the vertices kept in the extracted subnet, in the
    # order in which they will appear in that subnet. This order corresponds to that which
    # sorts the vertices by position, then puts unstable nodes at the end.
    I_kept = Int[_i]
    last_sortedcol = poscol[_i]

    for i in I_sort
        x = poscol[i]
        if x != last_sortedcol
            push!(I_kept, i)
            last_sortedcol = x
        end
    end

    kept_collisions = Int[]
    if !isempty(collisions)
        m = n - length(collisions)
        reorder = Vector{Int}(undef, length(I_kept))
        idx = 0
        for (j, i) in enumerate(I_kept)
            if i > m
                reorder[end-idx] = j
                idx += 1
                push!(kept_collisions, i-m)
            else
                reorder[j-idx] = j
            end
        end
        I_kept = I_kept[reorder]
        reverse!(@view I_kept[end-idx+1:end])
    end

    sortedcol = SVector{D,Rational{T}}[SVector{D,Rational{T}}(poscol[i]) for i in I_kept]
    # @toggleassert allunique(sortedcol)

    local vmap::Vector{Int}
    if isempty(collisions) # implies issorted(sortedcol)
        vmap = Vector{Int}(undef, n)
        # @toggleassert issorted(sortedcol)
        for (i, pos) in enumerate(poscol)
            j = searchsortedfirst(sortedcol, pos)
            @toggleassert j <= length(sortedcol) && sortedcol[j] == pos
            vmap[i] = j
        end
    else
        rev_dict = Dict{SVector{D,Rational{T}},Int}(pos => j for (j, pos) in enumerate(sortedcol))
        vmap = [rev_dict[pos] for pos in poscol]
    end

    edges = PeriodicEdge{D}[]
    for i in 1:length(I_kept)
        img = I_kept[i]
        ofs_i = offset[img]
        for neigh in neighbors(c.pge.g, img)
            ofs_x = offset[neigh.v]
            push!(edges, (i, vmap[neigh.v], ofs_x - ofs_i .+ imat*neigh.ofs))
        end
    end

    newcollisions = CollisionList(collisions, vmap, kept_collisions)
    for newnode in newcollisions
        if !allunique(newnode.neighs)
            # contravenes rule B of collision_nodes(::CrystalNet))
            return c, CollisionList(UnitRange{Int}[])
        end
    end

    graph = PeriodicGraph{D}(edges)
    opts = permute_mapping!(c.options, vmap)
    return CrystalNet{D,Rational{T}}(cell, c.types[I_kept], sortedcol, graph, opts), newcollisions
end


"""
    minimize(net::CrystalNet, [collisions])

Return a CrystalNet representing the same net as the input, but in a unit cell.
If `collisions` is given, also return the corresponding collisions after minimization.

The computed unit cell may depend on the representation of the input, i.e. it is not
topologicallly invariant.

If the input `collisions` is not a [`CollisionList`](@ref) but the returned one is, this
signals that minimization failed because of unknown collisions.
"""
function minimize(net::CrystalNet, collisions)
    collisions isa CollisionList || return minimize_unstable(net, collisions)
    collisions::CollisionList
    translations = find_all_valid_translations(net, collisions)
    all(isempty.(translations)) && return net, collisions
    mat = minimal_volume_matrix(translations)
    _net, newcollisions = reduce_with_matrix(net, mat, collisions)
    @toggleassert all(isempty.(find_all_valid_translations(_net, newcollisions)))
    return _net, newcollisions
end
minimize(net::CrystalNet) = minimize(net, CollisionList(net.pge.g, UnitRange{Int}[]))[1]


function minimize_unstable(shrunk_net::CrystalNet, collisions)
    while true
        x = find_first_valid_translation_unstable(shrunk_net, collisions)
        x isa Nothing && break
        shrunk_net, collisions = x
    end
    shrunk_net, collisions
end


function find_transformation_matrix(t::SVector{D,T}) where {D,T}
    _transformation = MMatrix{D,D,T,D*D}(LinearAlgebra.I)
    _transformation[:,1] .= t
    transformation = SMatrix{D,D,T,D*D}(_transformation)
    if D == 2 && issingular(transformation)
        _transformation[:,1] .= (1, 0)
        _transformation[:,2] .= t
        transformation = SMatrix{D,D,T,D*D}(_transformation)
    end
    if D == 3 && issingular(transformation)
        _transformation[:,1] .= (1, 0, 0)
        _transformation[:,2] .= t
        transformation = SMatrix{D,D,T,D*D}(_transformation)
        if issingular(transformation)
            _transformation[:,2] .= (0, 1, 0)
            _transformation[:,3] .= t
            transformation = SMatrix{D,D,T,D*D}(_transformation)
        end
    end
    transformation
end


"""
    orbits_pvmap(shrunk_pvmap::Vector{PeriodicVertex{D}}, m) where D

Given a list `shrunk_pvmap` that maps each vertex `i` of the graph to a new vertex
`shrunk_pvmap[i]` obtained from a valid translation of the graph, return a list of `m`
sublists `subgraphlists` such that:
- `subgraphlists[1]` is mapped to `subgraphlists[2]`, which is mapped to `subgraphlists[3]`,
  etc., and `subgraphlists[end]` is mapped to `subgraphlists[1]` plus an offset (which is
  `length(subgraphlists)-1` times the valid translation).
- `subgraphlists` forms a partition of the vertices of the graph: each vertex number
  appears in exactly one sublist.
- within each element of `subgraphlists`, the collision nodes appear after the rest.
- all vertices in `subgraphlists[1]` have a zero offset.
"""
function orbits_pvmap(shrunk_pvmap::Vector{PeriodicVertex{D}}, m) where D
    orbits = Vector{PeriodicVertex{D}}[]
    n = length(shrunk_pvmap)
    visited = falses(n)
    for i in 1:n
        visited[i] && continue
        visited[i] = true
        tot_ofs = zero(SVector{D,Int})
        orbit = [PeriodicVertex{D}(i)]
        j, ofs = shrunk_pvmap[i]
        k = 1
        while j != i
            @toggleassert !visited[j]
            visited[j] = true
            tot_ofs += ofs
            if k == m
                k = 0
                tot_ofs = zero(SVector{D,Int})
                push!(orbits, orbit)
                orbit = PeriodicVertex{D}[]
            end
            push!(orbit, PeriodicVertex{D}(j, tot_ofs))
            j, ofs = shrunk_pvmap[j]
            k += 1
        end
        @toggleassert k == m
        push!(orbits, orbit)
    end
    @toggleassert m == length(first(orbits))

    ret = Vector{Vector{PeriodicVertex{D}}}(undef, m)
    for i in 1:m
        ret[i] = [orbit[i] for orbit in orbits]
    end
    ret
end

"""
    find_ref_edges(net::CrystalNet{D}, collision_ranges, nodes, shrunk_pge, inv_transformation, refdict, collision_offsets) where D

Find the list of edges of the subnet obtained from `net` by only keeping the vertices
corresponding to the expansion of the `nodes` into the initial vertices they represent
through `collision_ranges`, obtained from a valid translation.

`shrunk_pge` is the `PeriodicGraphEmbedding` of `net` after grouping colliding vertices
into "collision nodes" (one node per set of colliding vertices). The `nodes` refer to such
collision nodes, and the correspondence between collision nodes and the initial vertices
is given by `collision_ranges` or, equivalently, through `collision_offsets`.

`inv_transformation` is the inverse transformation matrix, such that, given `pos` the
position of a vertex in `net`, `inv_transformation*pos` is the position of that same vertex
in the unit cell given the transformation obtained from the valid translation.
`refdict` maps each such position (after the transformation) to the vertex number of the
node in this new unit cell.

Note: all edges, both direct and indirect, are returned.
"""
function find_ref_edges(net::CrystalNet{D}, collision_ranges, nodes, shrunk_pge, inv_transformation, refdict, collision_offsets, vmap) where D
    refofss = [floor.(Int, inv_transformation * shrunk_pge[node]) for node in nodes]
    refedges = PeriodicEdge{D}[]
    first_collision_m1 = first(first(collision_ranges)) - 1

    i = 0
    for (_i, (_u, shrunk_ofs)) in enumerate(nodes)
        refofs = refofss[_i]
        for u in (_u > first_collision_m1 ? collision_ranges[_u - first_collision_m1] : _u:_u)
            v = vmap[u]
            i += 1
            for x in neighbors(net.pge.g, PeriodicVertex{D}(v, shrunk_ofs))
                pos = inv_transformation * net.pge[x]
                ofs = floor.(Int, pos)
                j = refdict[pos .- ofs] + collision_offsets[x.v]
                i == j && ofs == refofs && (empty!(refedges); return refedges) # forbidden edge: abort
                push!(refedges, PeriodicEdge{D}(i, j, ofs - refofs))
            end
        end
    end
    sort!(refedges)
    # display(refedges)
    @toggleassert begin # assert that each edge is unique
        flag = true
        for i in 1:(length(refedges)-1)
            if refedges[i] == refedges[i+1]
                flag = false
                break
            end
        end
        flag
    end
    refedges
end


"""
    reduce_unstable_net(shrunk_net::CrystalNet{D}, net, collision_ranges, shrunk_pvmap, transformation, collision_offsets) where D

Reduce the `net` with a translation symmetry, given as the corresponding `shrunk_pvmap`
that maps each vertex of the `shrunk_net` to its image, the cell matrix `transformation`,
and the `collision_offsets` such that vertices `u` and `v` of `net` are related by the
translation symmetry if and only if:
- they are in the same orbit and
- `collision_offsets[u] == collision_offsets[v]`.

Note that for all vertices `u` which are not in a collision node, `collision_offsets[u] == 1`.

Return `(new_shrunk_net, (new_net, new_collision_ranges))` which mirror the input
`(shrunk_net, (net, collision_ranges))`.
"""
function reduce_unstable_net(shrunk_net::CrystalNet{D}, net, collision_ranges, shrunk_pvmap, transformation, collision_offsets) where D
    inv_transformation = inv(transformation)
    subgraphlists = orbits_pvmap(shrunk_pvmap, Int(det(inv_transformation)))
    _subgraphlist_head = first(subgraphlists)
    @toggleassert all(iszero, last.(_subgraphlist_head))
    subgraphlist_head = first.(_subgraphlist_head)

    rev_shrunkmap = zeros(Int, length(shrunk_net.pge))
    for list in subgraphlists, (i, (v, _)) in enumerate(list)
        rev_shrunkmap[v] = i
    end

    first_collision_subgraphlist = 1 + (length(shrunk_net.pge) - length(collision_ranges))*length(subgraphlist_head)÷length(shrunk_net.pge)
    first_collision_m1 = first_collision_m1 = first(first(collision_ranges)) - 1

    rev_map = zeros(Int, length(net.pge))
    for list in subgraphlists
        counter = 0
        for (v, _) in list
            if v ≤ first_collision_m1
                counter += 1
                rev_map[v] = counter
            else
                rnge = collision_ranges[v-first_collision_m1]
                for u in rnge
                    colloffs = collision_offsets[u]
                    rev_map[u] = counter + colloffs
                end
                counter += length(rnge)
            end
        end
    end


    new_mat = Cell(net.pge.cell.mat * Nmatrix_to_3D(transformation))
    new_shrunk_pos = [inv_transformation * shrunk_net.pge.pos[v] for v in subgraphlist_head]
    new_shrunk_ofs = [floor.(Int, x) for x in new_shrunk_pos]
    for (i, ofs) in enumerate(new_shrunk_ofs); new_shrunk_pos[i] -= ofs; end

    shrunk_newedges = PeriodicEdge{D}[]
    newedges = PeriodicEdge{D}[]
    for (i, v) in enumerate(subgraphlist_head)
        ofsi = new_shrunk_ofs[i]
        for shrunkx in neighbors(shrunk_net.pge.g, v)
            posshrunkx = inv_transformation * shrunk_net.pge[shrunkx]
            newofs = floor.(Int, posshrunkx) - ofsi
            push!(shrunk_newedges, PeriodicEdge{D}(i, rev_shrunkmap[shrunkx.v], newofs))
        end
        for u in (v ≤ first_collision_m1 ? (v,) : collision_ranges[v-first_collision_m1]), x in neighbors(net.pge.g, u)
            posx = inv_transformation * net.pge[x]
            newofs = floor.(Int, posx) - ofsi
            push!(newedges, PeriodicEdge{D}(rev_map[u], rev_map[x.v], newofs))
        end
    end

    subcollisionindices = first_collision_subgraphlist:length(subgraphlist_head)
    virtualvmap = subgraphlist_head[1:first_collision_subgraphlist-1]
    new_collision_ranges = Vector{UnitRange{Int}}(undef, length(subcollisionindices))
    counter = first_collision_subgraphlist
    for (i, j) in enumerate(subcollisionindices)
        jv = subgraphlist_head[j]
        rnge = collision_ranges[jv - first_collision_m1]
        append!(virtualvmap, rnge)
        nextcounter = counter + length(rnge)
        new_collision_ranges[i] = counter:(nextcounter-1)
        counter = nextcounter
    end

    new_shrunk_pge = PeriodicGraphEmbedding{D}(PeriodicGraph{D}(shrunk_newedges), new_shrunk_pos, new_mat)
    new_shrunk_net = CrystalNet{D}(new_shrunk_pge, shrunk_net.types[subgraphlist_head], shrunk_net.options)

    new_pos = [inv_transformation * net.pge.pos[v] for v in virtualvmap]
    new_ofs = [floor.(Int, x) for x in new_pos]
    for (i, ofs) in enumerate(new_ofs); new_pos[i] -= ofs; end
    new_pge = PeriodicGraphEmbedding{D}(PeriodicGraph{D}(newedges), new_pos, new_mat)
    new_options = permute_mapping!(net.options, virtualvmap)
    new_net = CrystalNet{D}(new_pge, net.types[virtualvmap], new_options)
    return (new_shrunk_net, (new_net, new_collision_ranges))
end


function find_first_valid_translation_unstable(shrunk_net::CrystalNet{D,T}, collisions) where {D,T}
    check_symmetry = CheckSymmetryWithCollisions(collisions, false)
    translations_to_check = Tuple{SVector{D,T},Vector{PeriodicVertex{D}}}[]
    for pt in possible_translations(shrunk_net)
        t = last(pt)
        pvmap = check_symmetry(shrunk_net.pge, t, nothing, shrunk_net.types)
        if pvmap isa Vector{PeriodicVertex{D}}
            push!(translations_to_check, (t, pvmap))
        end
    end
    isempty(translations_to_check) && return nothing
    net, collision_ranges = collisions
    n = length(net.pge)
    first_collision_m1 = first(first(collision_ranges)) - 1
    collision_offsets = zeros(Int, n)
    collision_lengths = ones(Int, n)
    for rnge in collision_ranges
        for (i, x) in enumerate(rnge)
            collision_offsets[x] = i-1
            collision_lengths[x] = length(rnge)
        end
    end

    for (t, shrunk_pvmap) in translations_to_check

        ## First, compute the transformation matrix corresponding to this valid translation
        transformation = find_transformation_matrix(t)
        inv_transformation = inv(transformation)

        ## Identify how this transformation divides the net into subnets related by the translation
        subgraphlists = orbits_pvmap(shrunk_pvmap, Int(det(inv_transformation)))
        subgraphlist_head, subgraphlists_tail = Iterators.peel(subgraphlists)

        ## Take one such subnet and compute its corresponding periodic graph
        refdict = Dict{SVector{D,T},Int}() # map each vertex position in the new cell to the index of its first node
        shrunk_refdict = Dict{SVector{D,T},Int}() # same as refdict but without the collisions
        refpos = Vector{SVector{D,T}}(undef, length(subgraphlist_head)) # inverse of shrunk_refdict
        first_collision_subgraphlist = 0 # index of the first collision node in the new reduced graph
        counter_refdict = 1
        for (i, node) in enumerate(subgraphlist_head)
            oldpos = shrunk_net.pge[node]
            npos = inv_transformation * oldpos
            newpos = npos .- floor.(Int, npos)
            refdict[newpos] = counter_refdict
            shrunk_refdict[newpos] = i
            refpos[i] = newpos
            if first_collision_subgraphlist == 0 && node.v > first_collision_m1
                first_collision_subgraphlist = i
            end
            counter_refdict += collision_lengths[node.v]
        end
        @toggleassert length(refdict) == length(refpos) # assert that all positions are unique

        vmap = collect(1:n)
        refedges = find_ref_edges(net, collision_ranges, subgraphlist_head, shrunk_net.pge, inv_transformation, refdict, collision_offsets, vmap)
        cpci = ContiguousPlainChangesIterator(collect(Iterators.drop(collision_ranges, 1)))
        valid = false
        for swap in cpci
            valid = true
            # refedges will be the edges of the periodic graph obtained on the subnet
            refedges = find_ref_edges(net, collision_ranges, subgraphlist_head, shrunk_net.pge, inv_transformation, refdict, collision_offsets, vmap)
            for subgraphlist in subgraphlists_tail
                newedges = find_ref_edges(net, collision_ranges, subgraphlist, shrunk_net.pge, inv_transformation, refdict, collision_offsets, vmap)
                if newedges != refedges
                    valid = false
                    break
                end
            end
            valid && break

            collision_offsets[swap], collision_offsets[swap+1] = collision_offsets[swap+1], collision_offsets[swap]
            vmap[swap], vmap[swap+1] = vmap[swap+1], vmap[swap]
        end

        if valid # found a valid translation!
            return reduce_unstable_net(shrunk_net, net, collision_ranges, shrunk_pvmap, transformation, collision_offsets)
        end
    end
    return nothing
end

# Variants for a Crystal


function find_first_valid_translations(c::Crystal{Nothing})
    origin, rest = Iterators.peel(c.pge.pos)
    for pos in rest
        t = pos .- origin
        t = Float64.(rationalize.(Int8, t, 0.1))
        t -= floor.(Int, t)
        all(iszero.(t)) && continue
        vmap = check_valid_symmetry(c.pge, t, nothing, c.types)
        vmap isa Vector{PeriodicVertex3D} && return t
    end
    return nothing
end

function max_nearest(c::Crystal{Nothing}, trans)
    max = -Inf
    i_max = 0
    for (i, pos) in enumerate(c.pge.pos)
        t = pos .+ trans
        t -= floor.(Int, t)
        j = find_nearest(c.pge.pos, t)
        x = norm(c.pge.pos[j] .- t)
        if x > max
            max = x
            i_max = i
        end
    end
    return i_max, max
end

function dist2(x::SVector{D,T}, y::SVector{D,T}) where {D,T}
    r2 = zero(T)
    for j in 1:D
        r2 += (x[j] - y[j])^2
    end
    r2
end

function find_nearest(l::Vector{SVector{D,T}}, pos::SVector{D,T}) where {D,T}
    minr2 = Inf
    mini = 0
    for (i, x) in enumerate(l)
        r2 = dist2(pos, x)
        if r2 < minr2
            minr2 = r2
            mini = i
        end
    end
    return mini
end

function reduce_with_matrix(c::Crystal{Nothing}, _mat)
    D = 3
    mat = _mat isa SMatrix{3,3,BigFloat,9} ? _mat : SMatrix{3,3,BigFloat,9}(_mat)
    cell = Cell(c.pge.cell, c.pge.cell.mat * mat)

    imat = round.(Int, inv(mat)) # The inverse should only have integer coefficients
    poscol = (imat,) .* c.pge.pos # position in the new reference, will be in [0,1)
    n = length(poscol)
    offset = Vector{SVector{D,Int}}(undef, n) # offsets in the new reference
    for (i, pos) in enumerate(poscol)
        ofs = floor.(Int, pos)
        offset[i] = ofs
        poscol[i] = pos .- ofs
    end
    I_sort = sort(1:n; by=i->(poscol[i], hash_position(offset[i])))
    _i = popfirst!(I_sort)
    @toggleassert iszero(offset[_i])

    # I_kept is the list of index of the vertices kept in the extracted subnet, in the
    # order in which they will appear in that subnet. This order corresponds to that which
    # sorts the vertices by position, then puts unstable nodes at the end.
    I_kept = Int[_i]
    last_sortedcol = poscol[_i]
    for i in I_sort
        x = poscol[i]
        if norm(x - last_sortedcol) > 0.05
            push!(I_kept, i)
            last_sortedcol = x
        end
    end

    sortedcol = SVector{D,Float64}[SVector{D,Float64}(poscol[i]) for i in I_kept]
    # @toggleassert allunique(sortedcol)
    vmap = Vector{Int}(undef, n)
    # @toggleassert issorted(sortedcol)
    for (i, pos) in enumerate(poscol)
        j = find_nearest(sortedcol, pos)
        vmap[i] = j
    end

    edges = PeriodicEdge{D}[]
    for i in 1:length(I_kept)
        img = I_kept[i]
        ofs_i = offset[img]
        for neigh in neighbors(c.pge.g, img)
            ofs_x = offset[neigh.v]
            push!(edges, (i, vmap[neigh.v], ofs_x - ofs_i .+ imat*neigh.ofs))
        end
    end

    graph = PeriodicGraph{D}(edges)
    opts = permute_mapping!(c.options, vmap)
    return Crystal{Nothing}(PeriodicGraphEmbedding{D,Float64}(graph, sortedcol, cell), c.types[I_kept], opts)
end

"""
    minimize(net::Crystal{Nothing})

Return a Crystal representing the same crystal as the input, but in a unit cell.
"""
function minimize(c::Crystal{Nothing})
    trans = find_first_valid_translations(c)
    while !isnothing(trans)
        mat = if iszero(trans[1])
            if iszero(trans[2])
                @assert !iszero(trans[3])
                [[1, 0, 0];; [0, 1, 0];; trans]
            else
                [[1, 0, 0];; trans;; [0, 0, 1]]
            end
        else
            [trans;; [0, 1, 0];; [0, 0, 1]]
        end
        c = reduce_with_matrix(c, mat)
        trans = find_first_valid_translations(c)
    end
    return c
end

