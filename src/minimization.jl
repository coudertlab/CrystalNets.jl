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
    find_all_valid_translations(net::CrystalNet{D}) where D

Return a `D`-tuple of list of tuples `(i_max_den, max_den, t)` (see
[`possible_translations`](@ref) for interpretation) where the `n`-th list contains all
valid translations of the net having exactly `n-1` zeros.

A translation is valid if it maps exactly each vertex to a vertex and each edge to an edge.

See also: [`possible_translations`](@ref), `PeriodicGraphEmbeddings.check_valid_symmetry`
"""
function find_all_valid_translations(net::CrystalNet{D,T}) where {D,T}
    ret = NTuple{D, Vector{Tuple{Int, Int, SVector{D,T}}}}(ntuple(_->[], Val(D)))
    for (nz, i_max_den, max_den, t) in possible_translations(net)
        vmap = check_valid_symmetry(net.pge, t, nothing, net.types)
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
    reduce_with_matrix(c::CrystalNet, mat)

Given the net and the output of `minimal_volume_matrix` computed on the valid translations
of the net, return the new net representing the initial net in the computed unit cell.
"""
function reduce_with_matrix(c::CrystalNet{D,Rational{T}}, mat) where {D,T}
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

    sortedcol = SVector{D,Rational{T}}[SVector{D,Rational{T}}(poscol[i]) for i in I_kept]
    # @toggleassert allunique(sortedcol)

    vmap = Vector{Int}(undef, n)
    # @toggleassert issorted(sortedcol)
    for (i, pos) in enumerate(poscol)
        j = searchsortedfirst(sortedcol, pos)
        @toggleassert j <= length(sortedcol) && sortedcol[j] == pos
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
    return CrystalNet{D,Rational{T}}(cell, c.types[I_kept], sortedcol, graph, opts)
end


"""
    minimize(net::CrystalNet, [collisionsetup])

Return a CrystalNet representing the same net as the input, but in a unit cell.
If `collisionsetup` is given, also return the corresponding collisions after minimization.

The computed unit cell may depend on the representation of the input, i.e. it is not
topologicallly invariant.
"""
function minimize(net::CrystalNet, collisionsetup)
    isempty(collisionsetup[2]) || return minimize_unstable(net, collisionsetup)
    translations = find_all_valid_translations(net)
    all(isempty.(translations)) && return net, collisionsetup
    mat = minimal_volume_matrix(translations)
    _net = reduce_with_matrix(net, mat)
    @toggleassert all(isempty.(find_all_valid_translations(_net)))
    return _net, collisionsetup
end
minimize(net::CrystalNet) = minimize(net, (net, CollisionList()))[1]


function minimize_unstable(shrunk_net::CrystalNet, collisionsetup)
    while true
        x = find_first_valid_translation_unstable(shrunk_net, collisionsetup)
        x isa Nothing && break
        shrunk_net, collisionsetup = x
    end
    shrunk_net, collisionsetup
end


"""
    find_transformation_matrix(t::SVector{D,T}) where {D,T}

Given a valid translation `t` of a periodic graph, return the smallest unit cell that can be
extracted from that translation.

The returned matrix only has non-negative coefficients.
"""
function find_transformation_matrix(t::SVector{D,T}) where {D,T}
    images_set = Set((t - floor.(Int, t),))
    current_image = t
    while true
        current_image += t
        all(isinteger, current_image) && break
        push!(images_set, current_image - floor.(Int, current_image))
    end
    images = collect(images_set)
    n = length(images)
    periodicity = 1 + length(images)
    ret = if D == 1
        SMatrix{1,1,T,1}(minimum(x[1] for x in images))
    elseif D == 2
        push!(images, SA[one(T), zero(T)], SA[zero(T), one(T)])
        n += 2
        transformation2 = zero(T == Rational{BigInt} ? SizedMatrix{2,2,T,2,Matrix{T}} : MMatrix{2,2,T,4})
        mindet2 = (typemax(Int)-1)//1
        best2 = (0,0)
        for i in 1:(n-1)
            transformation2[:,1] = images[i]
            for j in (i+1):n
                transformation2[:,2] = images[j]
                d = abs(det(transformation2))
                if 0 < d < mindet2
                    mindet2 = d
                    best2 = (i,j)
                end
            end
        end
        SMatrix{2,2,T,4}(images[best2[1]]..., images[best2[2]]...)
    else
        @toggleassert D == 3
        oT = one(T); zT = zero(T)
        push!(images, SA[oT, zT, zT], SA[zT, oT, zT], SA[zT, zT, oT])
        n += 3
        transformation3 = zero(T == Rational{BigInt} ? SizedMatrix{3,3,T,2,Matrix{T}} : MMatrix{3,3,T,9})
        mindet3 = (typemax(Int)-1)//1
        best3 = (0,0,0)
        for i in 1:(n-1)
            transformation3[:,1] = images[i]
            for j in (i+1):n
                transformation3[:,2] = images[j]
                for k in (j+1):n
                    transformation3[:,3] = images[k]
                    d = abs(det(transformation3))
                    if 0 < d < mindet3
                        mindet3 = d
                        best3 = (i,j,k)
                    end
                end
            end
        end
        SMatrix{3,3,T,9}(images[best3[1]]..., images[best3[2]]..., images[best3[3]]...)
    end
    @toggleassert abs(LinearAlgebra.det(ret)) == 1//periodicity
    ret
end


"""
    orbits_pvmap(shrunk_pvmap::Vector{PeriodicVertex{D}}) where D

Given a list `shrunk_pvmap` that maps each vertex `i` of the graph to a new vertex
`shrunk_pvmap[i]` obtained from a valid translation of the graph, return a list of
sublists `subgraphlists` such that:
- `subgraphlists[1]` is mapped to `subgraphlists[2]`, which is mapped to `subgraphlists[3]`,
  etc., and `subgraphlists[end]` is mapped to `subgraphlists[1]` plus an offset (which is
  `length(subgraphlists)-1` times the valid translation).
- `subgraphlists` forms a partition of the vertices of the graph: each vertex number
  appears in exactly one sublist.
- within each element of `subgraphlists`, the collision nodes appear after the rest.
- all vertices in `subgraphlists[1]` have a zero offset.
"""
function orbits_pvmap(shrunk_pvmap::Vector{PeriodicVertex{D}}) where D
    orbits = Vector{PeriodicVertex{D}}[]
    n = length(shrunk_pvmap)
    visited = falses(n)
    for i in 1:n
        visited[i] && continue
        visited[i] = true
        tot_ofs = zero(SVector{D,Int})
        orbit = [PeriodicVertex{D}(i)]
        j, ofs = shrunk_pvmap[i]
        while j != i
            @toggleassert !visited[j]
            visited[j] = true
            tot_ofs += ofs
            push!(orbit, PeriodicVertex{D}(j, tot_ofs))
            j, ofs = shrunk_pvmap[j]
        end
        push!(orbits, orbit)
    end
    @toggleassert allequal(length, orbits)

    m = length(first(orbits))
    ret = Vector{Vector{PeriodicVertex{D}}}(undef, m)
    for i in 1:m
        ret[i] = [orbit[i] for orbit in orbits]
    end
    ret
end


"""
    reduce_unstable_net(shrunk_net::CrystalNet{D}, net, collisions, shrunk_pvmap, transformation, collision_offsets) where D

Reduce the `net` with a translation symmetry, given as the corresponding `shrunk_pvmap`
that maps each vertex of the `shrunk_net` to its image, the cell matrix `transformation`,
and the `collision_offsets` such that vertices `u` and `v` of `net` are related by the
translation symmetry if and only if:
- they are in the same orbit and
- `collision_offsets[u] == collision_offsets[v]`.

Note that for all vertices `u` which are not in a collision node, `collision_offsets[u] == 1`.

Return `(new_shrunk_net, (new_net, new_collisions))` which mirror the input
`(shrunk_net, (net, collisions))`.
"""
function reduce_unstable_net(shrunk_net::CrystalNet{D}, net, collisions, shrunk_pvmap, transformation, collision_offsets) where D
    inv_transformation = Int.(inv(transformation))
    periodicity = abs(det(inv_transformation))
    subgraphlists = orbits_pvmap(shrunk_pvmap)
    @toggleassert periodicity == length(subgraphlists)
    _subgraphlist_head = first(subgraphlists)
    @toggleassert all(iszero, last.(_subgraphlist_head))
    subgraphlist_head = first.(_subgraphlist_head)

    rev_shrunkmap = zeros(Int, length(shrunk_net.pge))
    for list in subgraphlists, (i, (v, _)) in enumerate(list)
        rev_shrunkmap[v] = i
    end

    first_collision_subgraphlist = 1 + (length(shrunk_net.pge) - length(collisions))*length(subgraphlist_head)÷length(shrunk_net.pge)
    first_collision_m1 = first_collision_m1 = first(first(collisions).rnge) - 1

    rev_map = zeros(Int, length(net.pge))
    for list in subgraphlists
        counter = 0
        for (v, _) in list
            if v ≤ first_collision_m1
                counter += 1
                rev_map[v] = counter
            else
                rnge = collisions[v-first_collision_m1].rnge
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
        for u in (v ≤ first_collision_m1 ? (v,) : collisions[v-first_collision_m1].rnge), x in neighbors(net.pge.g, u)
            posx = inv_transformation * net.pge[x]
            newofs = floor.(Int, posx) - ofsi
            push!(newedges, PeriodicEdge{D}(rev_map[u], rev_map[x.v], newofs))
        end
    end

    subcollisionindices = first_collision_subgraphlist:length(subgraphlist_head)
    virtualvmap = subgraphlist_head[1:first_collision_subgraphlist-1]
    new_collisions = Vector{CollisionNode}(undef, length(subcollisionindices))
    counter = first_collision_subgraphlist
    for (i, j) in enumerate(subcollisionindices)
        jv = subgraphlist_head[j]
        node = collisions[jv - first_collision_m1]
        append!(virtualvmap, node.rnge)
        nextcounter = counter + length(node.rnge)
        new_collisions[i] = CollisionNode(counter:(nextcounter-1), node.uniquecsequences, node.uniques, node.subranges)
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
    return (new_shrunk_net, (new_net, CollisionList(new_collisions)))
end

function direct_map_to_collision_offsets(direct_map, collision_ranges)
    first_collision_m1 = first(first(collision_ranges)) - 1
    collision_offsets = zeros(Int, length(direct_map))
    collision_offsets[1:first_collision_m1] .= 1
    for rnge in collision_ranges
        if collision_offsets[first(rnge)] != 0
            @toggleassert all(!iszero, @view(collision_offsets[rnge]))
        else
            for (i, j) in enumerate(rnge)
                @toggleassert collision_offsets[j] == 0
                collision_offsets[j] = i
                k = direct_map[j]
                while k != j
                    @toggleassert collision_offsets[k] == 0
                    collision_offsets[k] = i
                    k = direct_map[k]
                end
            end
        end
    end
    @toggleassert all(!iszero, collision_offsets)
    collision_offsets
end

# ensure that there cannot be any orbit of the form i -> ... -> j -> ... where i and j are in the same collision node
function find_necessary_image(iactual, (direct_map, reverse_map, collisions, to_shrunk))
    first_collision_m1 = first(first(collisions).rnge) - 1
    directactual = direct_map[iactual]
    irnge = directactual < 0 ? -directactual : (to_shrunk[directactual-first_collision_m1] - first_collision_m1)
    currentrnge = to_shrunk[iactual-first_collision_m1] - first_collision_m1
    postprev = iactual
    prev = reverse_map[postprev]
    niter = 0
    nitermax = length(direct_map)
    while niter < nitermax && prev > 0 && to_shrunk[prev-first_collision_m1] != irnge
        niter += 1
        postprev = prev
        prev = reverse_map[postprev]
    end
    niter == nitermax && return -1
    prev > 0 && return prev
    -prev == currentrnge && return postprev
    return 0
end

function attribute_next_available!(vmapprogress, index, iactual, modifystack)
    attribute_next_available_modify!(vmapprogress, 0, index, iactual, 0, modifystack)
end

function attribute_next_available_modify!(vmapprogress, thisindex, returnto, iactual, exclude, modifystack)
    direct_map, reverse_map, collisions, to_shrunk = modifystack
    first_collision_m1 = first(first(collisions).rnge) - 1
    directactual = direct_map[iactual]
    irnge = directactual < 0 ? -directactual : (to_shrunk[directactual-first_collision_m1] - first_collision_m1)
    necessaryimage = find_necessary_image(iactual, modifystack)
    necessaryimage == -1 && return false
    crnge = collisions[irnge].rnge
    rnge = @view crnge[exclude+1:end]
    idxinrnge = if necessaryimage > 0
        necessaryimage in rnge ? necessaryimage-first(rnge)+1 : nothing
    else
        findfirst(<(0), @view reverse_map[rnge])
    end
    idxinrnge isa Nothing && return false
    if thisindex == 0
        push!(vmapprogress, (returnto, iactual => idxinrnge+exclude))
    else
        vmapprogress[thisindex] = (returnto, iactual => idxinrnge+exclude)
    end
    newval = rnge[idxinrnge]
    previous = direct_map[iactual]
    @toggleassert previous != newval
    direct_map[iactual] = newval
    if previous > 0
        reverse_map[previous] = reverse_map[newval]
    end
    reverse_map[newval] = iactual
    return true
end

function backtrack_to!(vmapprogress, index, (modifystack, shrunk_pvmap, reference_direct_map, reference_reverse_map))
    direct_map, reverse_map, collisions, _ = modifystack
    returnto, (before, idx_after) = vmapprogress[index]
    for (_, (i, j)) in @view vmapprogress[index+1:end]
        mi_rnge = direct_map[i] = reference_direct_map[i]
        k = collisions[-mi_rnge].rnge[j]
        reverse_map[k] = reference_reverse_map[k]
    end
    resize!(vmapprogress, index)
    wasavail_backtrack = attribute_next_available_modify!(vmapprogress, index, returnto, before, idx_after, modifystack)
    if wasavail_backtrack # stop backtracking
        false, returnto
    else # continue backtracking
        # k0 = collision_ranges[-reference_direct_map[before]][idx_after]
        # reverse_map[k0] = reference_reverse_map[k0]
        true, index - 1
    end
end

"""
    translation_to_direct_map(shrunk_pvmap, net, collisions::CollisionList, to_shrunk)

Given a `shrunk_pvmap` representing the mapping of the shrunk vertices obtained from a
translation (valid on the shrunk net), return a `direct_map` which is a mapping of the
vertex numbers of the initial `net` obtained from the same translation.

If no such `direct_map` can be found, i.e. if the translation is not valid on the initial
net, return `nothing`.
"""
function translation_to_direct_map(shrunk_pvmap, net, collisions::CollisionList, to_shrunk)
    # @show shrunk_pvmap, collisions, to_shrunk
    @toggleassert !isempty(collisions)
    valid = true # determine if the translation is valid
    first_collision_m1 = first(first(collisions).rnge) - 1
    shrunk_indirect_vmap = zeros(Int, length(shrunk_pvmap) - first_collision_m1)
    for (i, (v, _)) in enumerate(@view shrunk_pvmap[first_collision_m1+1:end])
        shrunk_indirect_vmap[v - first_collision_m1] = i
    end
    n = length(net.pge)
    direct_map = zeros(Int, n)
    reverse_map = zeros(Int, n)
    for i in 1:first_collision_m1
        direct_map[i] = _v = first(shrunk_pvmap[i])
        @toggleassert reverse_map[_v] == 0
        reverse_map[_v] = i
    end
    last_j = first_collision_m1
    targetnodeidx = 0
    cofs = 0
    for (i, node) in enumerate(collisions)
        targetnodeidx = first(shrunk_pvmap[first_collision_m1+i]) - first_collision_m1
        targetnode = collisions[targetnodeidx]
        @toggleassert possibly_equivalent_nodes(node, targetnode)
        cofs = first(node.rnge) - 1
        for subrnge in node.subranges
            ofssubrnge = cofs .+ subrnge
            for j in ((last_j+1):(first(ofssubrnge)-1))
                _v = targetnode.rnge[j - cofs]
                direct_map[j] = _v
                reverse_map[_v] = j
            end
            for j in ofssubrnge
                direct_map[j] = -targetnodeidx
                reverse_map[j] = -shrunk_indirect_vmap[i]
            end
            last_j = last(ofssubrnge)
        end
        for j in (last_j+1):last(node.rnge)
            _v = targetnode.rnge[j - cofs]
            direct_map[j] = _v
            reverse_map[_v] = j
        end
        last_j = last(node.rnge)
    end
    @toggleassert !any(iszero, direct_map) && !any(iszero, reverse_map)
    modifystack = direct_map, reverse_map, collisions, to_shrunk
    backtrackstack = (modifystack, shrunk_pvmap, copy(direct_map), copy(reverse_map))

    for grand_index in (first_collision_m1+1):n
        vgi = direct_map[grand_index]
        vgi > 0 && continue
        vmapprogress = Tuple{Int,Pair{Int,Int}}[]
        grand_wasavail = attribute_next_available!(vmapprogress, 1, grand_index, modifystack)
        @toggleassert grand_wasavail
        index = 1
        valid = false
        backtracking = false

        while index > 0
            yield()
            # @show vmapprogress
            returnto, (before, idx_after) = vmapprogress[index]
            shrunk_before = to_shrunk[before-first_collision_m1] # shrunk node to which "before" belongs
            shrunk_after = first(shrunk_pvmap[shrunk_before]) # shrunk node to which "after" belongs
            rnge = collisions[shrunk_after - first_collision_m1].rnge

            # "before" is currently mapped to "rnge[after]"
            if backtracking
                backtracking, index = backtrack_to!(vmapprogress, index, backtrackstack)
                continue
            end

            after = rnge[idx_after]
            m = degree(net.pge.g, before)
            if m != degree(net.pge.g, after)
                backtracking, index = backtrack_to!(vmapprogress, index, backtrackstack)
                continue
            end

            hadaction = false
            actual_neighbors = neighbors(net.pge.g, before)
            expected_neighbors = [begin
                w = reverse_map[x.v]
                src = to_shrunk[before - first_collision_m1]
                dst = w < 0 ? -w + first_collision_m1 : w ≤ first_collision_m1 ? w : to_shrunk[w - first_collision_m1]
                PeriodicVertex(w, x.ofs + shrunk_pvmap[src].ofs - shrunk_pvmap[dst].ofs)
            end for x in neighbors(net.pge.g, after)]
            # a negative vertex in expected_neighbors indicates a vertex to be assigned part of a given collision node (whose index is the opposite of that given)

            sort!(expected_neighbors; by=x->x.v>0 ? x : PeriodicVertex(typemin(Int) - x.v, x.ofs))
            first_attributed_neighbor_m1 = (@something findfirst(x -> x.v > 0, expected_neighbors) length(expected_neighbors)+1) - 1
            met_attributed = falses(m)
            failure = false
            for i_attributed in (first_attributed_neighbor_m1+1):m
                failure = true
                expected = expected_neighbors[i_attributed]
                for (i_actual, actual) in enumerate(actual_neighbors)
                    met_attributed[i_actual] && continue
                    if actual == expected
                        met_attributed[i_actual] = true
                        failure = false
                        break
                    end
                end
                failure && break
            end
            if failure
                # TODO: improve precision of backtrack_target for performance
                backtracking, index = backtrack_to!(vmapprogress, length(vmapprogress), backtrackstack)
                continue
            end

            nonattributed_progress = falses(first_attributed_neighbor_m1)
            backtrack_target = 0
            for (i_actual, actual) in enumerate(actual_neighbors)
                met_attributed[i_actual] && continue

                # mapping of a colliding vertex
                shrunk_actual = PeriodicVertex(-to_shrunk[actual.v - first_collision_m1] + first_collision_m1, actual.ofs)
                targetrnge = first(shrunk_pvmap[to_shrunk[actual.v - first_collision_m1]]) - first_collision_m1

                # check that there is an available expected neighbor compatible with the actual one
                failure = true
                for i in 1:first_attributed_neighbor_m1
                    nonattributed_progress[i] && continue
                    expected_neighbors[i] == shrunk_actual || continue
                    nonattributed_progress[i] = true
                    failure = false
                    break
                end
                failure && break

                attributed_previously = direct_map[actual.v]
                if attributed_previously > 0
                    # this situation can happen due to interference of a previous iteration
                    # of the "for actual in actual_neighbors" loop. For example, if there
                    # are two neighbors with the same vertex number but different offsets
                    if !(attributed_previously in collisions[targetrnge].rnge)
                        backtrack_target = findlast(x -> first(x[2]) == actual.v, vmapprogress)
                        break
                    end
                else
                    wasavail = attribute_next_available!(vmapprogress, index, actual.v, modifystack)
                    if wasavail
                        hadaction = true
                    else
                        failure = true
                        backtrack_target = length(vmapprogress)
                    end
                end
                failure && break
            end
            if failure
                if backtrack_target == 0
                    backtrack_target = length(vmapprogress)
                end
                backtracking, index = backtrack_to!(vmapprogress, backtrack_target, backtrackstack)
                continue
            end

            if !hadaction && index == length(vmapprogress)
                valid = true
                break
            end
            index += 1
        end

        valid || break
    end
    valid ? direct_map : nothing
end


function find_first_valid_translation_unstable(shrunk_net::CrystalNet{D,T}, (net, collisions)) where {D,T}
    check_symmetry = CheckSymmetryWithCollisions(collisions)
    translations_to_check = Tuple{SVector{D,T},Vector{PeriodicVertex{D}}}[]
    for pt in possible_translations(shrunk_net)
        t = last(pt)
        pvmap = check_symmetry(shrunk_net.pge, t, nothing, shrunk_net.types)
        if pvmap isa Vector{PeriodicVertex{D}}
            push!(translations_to_check, (t, pvmap))
        end
    end
    isempty(translations_to_check) && return nothing
    first_collision_m1 = first(first(collisions).rnge) - 1

    to_shrunk = [i+first_collision_m1 for (i, node) in enumerate(collisions) for _ in node.rnge]
    # to_shrunk[j-first_collision_m1] is the index of the shrunk node corresponding to vertex j

    collision_ranges = [node.rnge for node in collisions]
    for (t, shrunk_pvmap) in translations_to_check
        direct_map = translation_to_direct_map(shrunk_pvmap, net, collisions, to_shrunk)
        if !(direct_map isa Nothing) # found a valid translation!
            @toggleassert isperm(direct_map)
            transformation = find_transformation_matrix(t)
            collision_offsets = direct_map_to_collision_offsets(direct_map, collision_ranges)
            return reduce_unstable_net(shrunk_net, net, collisions, shrunk_pvmap, transformation, collision_offsets)
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
    poscol = [imat * round.(x; digits=12) for x in c.pge.pos] # position in the new reference, will be in [0,1)
    n = length(poscol)
    offset = Vector{SVector{D,Int}}(undef, n) # offsets in the new reference
    for (i, pos) in enumerate(poscol)
        ofs = floor.(Int, pos)
        offset[i] = ofs
        poscol[i] = round.(pos .- ofs; digits=12)
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

