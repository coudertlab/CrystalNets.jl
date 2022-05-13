## Main functions of the algorithm


"""
    check_dimensionality(c::CrystalNet)

Check that the dimensionality of the net (i.e. the number of independent axes along which
it is periodic) is equal to `D`, or throw a DimensionMismatch otherwise.
"""
function check_dimensionality(c::CrystalNet{D}) where {D}
    edgs = [c.pge.pos[dst(e)] .+ PeriodicGraphs.ofs(e) .- c.pge.pos[src(e)] for e in edges(c.pge.g)]
    sort!(edgs)
    unique!(edgs)
    mat = reduce(hcat, edgs)
    if D == 3
        isrank3(mat) || throw(DimensionMismatch("Internal error: the input net does not have expected dimensionality 3."))
    elseif D == 2
        isrank2(mat) || throw(DimensionMismatch("Internal error: the input net does not have expected dimensionality 2."))
    elseif D == 1
        isrank1(mat) || throw(DimensionMismatch("Internal error: the input net does not have expected dimensionality 1."))
    else
        throw(AssertionError("1 ≤ D ≤ 3"))
    end
    nothing
end


"""
    possible_translations(c::CrystalNet)

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
    sort!(sortedpos, by=norm)
    for t in sortedpos
        @toggleassert t == back_to_unit.(t)
        max_den, i_max_den = findmax(denominator.(t))
        numerator(t[i_max_den]) == 1 || continue
        nz = count(iszero, t)
        push!(ts, (nz, i_max_den, max_den, t))
    end
    return sort!(ts; by=(x->(x[1], x[2], x[3])))
end


"""
    find_all_valid_translations(c::CrystalNet{D}, collisions) where D

Return a `D`-tuple of list of tuples `(i_max_den, max_den, t)` (see
[`possible_translations`](@ref) for interpretation) where the `n`-th list contains all
valid translations of the net having exactly `n-1` zeros.

A translation is valid if it maps exactly each vertex to a vertex and each edge to an edge.

See also: [`possible_translations`](@ref), `PeriodicGraphEmbeddings.check_valid_symmetry`
"""
function find_all_valid_translations(shrunk_net::CrystalNet{D,T}, collisions) where {D,T}
    ret = NTuple{D, Vector{Tuple{Int, Int, SVector{D,T}}}}(ntuple(_->[], Val(D)))
    check_symmetry = check_symmetry_with_collisions(collisions)
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
    reduce_with_matrix(c::CrystalNet, mat, shrunk_net, collisions)

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
    poscol = (imat,) .* c.pge.pos
    n = length(poscol)
    offset = Vector{SVector{D,Int}}(undef, n)
    for (i, pos) in enumerate(poscol)
        ofs = floor.(Int, pos)
        offset[i] = ofs
        poscol[i] = pos .- ofs
    end
    I_sort = sort(1:n; by=i->(poscol[i], hash_position(offset[i])))
    _i = popfirst!(I_sort)
    @toggleassert iszero(offset[_i])
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
    end

    sortedcol = SVector{D,Rational{T}}[SVector{D,Rational{T}}(poscol[i]) for i in I_kept]

    vmap = Vector{Int}(undef, n)
    for (i, pos) in enumerate(poscol)
        j = searchsortedfirst(sortedcol, pos)
        @toggleassert j <= length(sortedcol) && sortedcol[j] == pos
        vmap[i] = j
    end

    edges = PeriodicEdge{D}[]
    for i in 1:length(I_kept)
        ofs_i = offset[I_kept[i]]
        for neigh in neighbors(c.pge.g, I_kept[i])
            ofs_x = offset[neigh.v]
            push!(edges, (i, vmap[neigh.v], ofs_x - ofs_i .+ imat*neigh.ofs))
        end
    end

    newcollisions = [CollisionNode(collisions[i], vmap) for i in kept_collisions]
    for newnode in newcollisions
        if !allunique(newnode.neighs)
            newcollisions = nothing
            break
        end
    end

    graph = PeriodicGraph{D}(edges)
    return CrystalNet{D,Rational{T}}(cell, c.types[I_kept], sortedcol, graph, c.options), newcollisions
end


"""
    minimize(net::CrystalNet, [collisions::Vector{CollisionNode}])

Return a CrystalNet representing the same net as the input, but in a unit cell.
If `collisions` is given, also return the corresponding collisions after minimization.

The computed unit cell may depend on the representation of the input, i.e. it is not
topologicallly invariant.
"""
function minimize(net::CrystalNet, collisions::Vector{CollisionNode})
    translations = find_all_valid_translations(net, collisions)
    all(isempty.(translations)) && return net, collisions
    mat = minimal_volume_matrix(translations)
    _net, collisions = reduce_with_matrix(net, mat, collisions)
    @toggleassert all(isempty.(find_all_valid_translations(_net, collisions)))
    return _net, collisions
end
minimize(net::CrystalNet) = minimize(net, CollisionNode[])[1]


function findfirstbasis(offsets)
    newbasis = normal_basis_rational(offsets)
    invbasis = inv(newbasis)
    intcoords = [Int.(invbasis * x) for x in offsets]
    return newbasis, intcoords
end

function findbasis(edges::Vector{Tuple{Int,Int,SVector{D,T}}}) where {D,T}
    m = length(edges)
    positivetrans = SVector{D,T}[]
    map_to_ofs = Vector{Int}(undef, m)
    tmp_map = Int[]
    for i in 1:m
        trans = last(edges[i])
        c = cmp(trans, zero(SVector{D,Int}))
        if iszero(c)
            map_to_ofs[i] = 0
        elseif c < 0
            push!(tmp_map, -i)
            push!(positivetrans, .-trans)
        else
            push!(tmp_map, i)
            push!(positivetrans, trans)
        end
    end

    I_sort = sortperm(positivetrans)
    uniques = SVector{D,T}[]
    last_trans = zero(SVector{D,T})

    for j in I_sort
        i = tmp_map[j]
        trans = positivetrans[j]
        if trans != last_trans
            push!(uniques, trans)
            last_trans = trans
        end
        map_to_ofs[abs(i)] = sign(i) * length(uniques)
    end

    basis, unique_coords = findfirstbasis(uniques)

    newedges = Vector{PeriodicEdge{D}}(undef, m)
    for j in 1:m
        i = map_to_ofs[j]
        src, dst, _ = edges[j]
        if i == 0
            newedges[j] = PeriodicEdge{D}(src, dst, zero(SVector{D,Int}))
        else
            newedges[j] = PeriodicEdge{D}(src, dst, sign(i) .* unique_coords[abs(i)])
        end
    end

    @toggleassert all(z-> begin x, y = z; x == (y.src, y.dst.v, basis*y.dst.ofs) end, zip(edges, newedges))

    return basis, newedges
end


"""
    candidate_key(net::CrystalNet, u, basis, minimal_edgs)

Given the net, a candidate `u => basis` where `u` is the origin and `basis` the triplet of
axes, and `minimal_edgs` the last minimal key (for the pseudo-lexicographical order used),
extract the key corresponding to the current candidate.

The key is the lexicographically ordered list of edges of the graph when its vertices are
numbered according to the candidate. The ordering of keys first compares the list of edges
disregarding the offsets, and then only compares the offsets if the rest is identical.

If the key is larger or equal to `minimal_edgs`, early stop and return two empty lists.
Otherwise, the extracted key is the current best: return the vmap between the initial
vertices and their ordered image in the candidate, as well as the key.

See also: [`find_candidates`](@ref)
"""
function candidate_key(net::CrystalNet{D,T}, u, basis, minimal_edgs) where {D,T}
    n = nv(net.pge.g)
    h = 2 # next node to assign
    origin = net.pge.pos[u]
    newpos = Vector{SVector{D,T}}(undef, n) # positions of the kept representatives
    newpos[1] = zero(SVector{D,T})
    offsets = Vector{SVector{D,Int32}}(undef, n) # offsets of the new representatives with
    # respect to the original one, in the original basis
    offsets[1] = zero(SVector{D,Int32})
    vmap = Vector{Int}(undef, n) # bijection from the old to the new node number
    vmap[1] = u
    rev_vmap = zeros(Int, n) # inverse of vmap
    rev_vmap[u] = 1
    flag_bestedgs = false # marks whether the current list of edges is lexicographically
    # below the best known one. If not by the end of the algorithm, return false
    edgs = Tuple{Int,Int,SVector{D,T}}[]
    bigbasis = T == Rational{BigInt} ? basis : widen(T).(basis)
    mat = T == Rational{BigInt} ? inv(bigbasis) : T.(inv(bigbasis))
    for t in 1:n # t is the node being processed
        neighs = neighbors(net.pge.g, vmap[t])
        ofst = offsets[t]
        pairs = Vector{Tuple{SVector{D,T},Int}}(undef, length(neighs))
        for (i,x) in enumerate(neighs)
            pairs[i] = ((mat*(net.pge.pos[x.v] .+ x.ofs .- origin .+ ofst)), x.v)
        end
        # (x,i) ∈ pairs means that vertex i (in the old numerotation) has position x in the new basis
        order = unique(last.(sort(pairs)))
        inv_order = Vector{Int}(undef, n)
        for i in 1:length(order)
            inv_order[order[i]] = i
        end
        sort!(pairs, lt = ((x, a), (y, b)) -> begin inv_order[a] < inv_order[b] ||
                                                    (inv_order[a] == inv_order[b] && x < y) end)
        # pairs is sorted such that different nodes first appear in increasing order of their position
        # but different representatives of the same node are contiguous and also sorted by position.
        for (coordinate, v) in pairs
            idx = rev_vmap[v]
            if idx == 0 # New node to which h is assigned
                @toggleassert t < h
                vmap[h] = v
                rev_vmap[v] = h
                newpos[h] = coordinate
                offsets[h] = SVector{D,Int32}(bigbasis * coordinate .+ origin .- net.pge.pos[v])
                push!(edgs, (t, h, zero(SVector{D,T})))
                h += 1
            else
                realofs = coordinate .- newpos[idx]
                # offset between this representative of the node and that which was first encountered
                push!(edgs, (t, idx, realofs))
            end
            if !flag_bestedgs
                j = length(edgs)
                c = cmp(minimal_edgs[j], edgs[j])
                c < 0 && return (Int[], Tuple{Int,Int,SVector{D,T}}[]) # early stop
                c > 0 && (flag_bestedgs = true)
            end
        end
    end
    @toggleassert allunique(edgs)
    if !flag_bestedgs # the current list of edges is equal to minimal_edgs
        @toggleassert minimal_edgs == edgs
        return (Int[], Tuple{Int,Int,SVector{D,T}}[])
    end
    return vmap, edgs
end


"""
    partition_by_coordination_sequence(graph, symmetries::AbstractSymmetryGroup=NoSymmetryGroup(graph))

Partition the vertices of the graph into disjoint categories, one for each coordination
sequence. The partition is then sorted by order of coordination sequence.
This partition does not depend on the representation of the graph.
The optional argument `vmaps` is a set of permutations of the vertices that leave the graph
unchanged. In other words, `vmaps` is a set of symmetry operations of the graph.

Return the categories and a list of unique representative for each symmetry class.
"""
function partition_by_coordination_sequence(graph, symmetries::AbstractSymmetryGroup=NoSymmetryGroup(graph))
    n = nv(graph)
    uniques::Vector{Int} = unique(symmetries)
    unique_reprs = [Int[repr] for repr in uniques]
    csequences = [coordination_sequence(graph, repr, 10) for repr in uniques]
    I = sortperm(csequences)
    @toggleassert csequences[I[1]][1] >= 2 # vertices of degree <= 1 should have been removed at input creation

    if symmetries isa NoSymmetryGroup
        categories = [Int[i] for i in 1:n]
    else
        rev_uniques = Vector{Int}(undef, n)
        for (i, repr) in enumerate(uniques)
            rev_uniques[repr] = i
        end
        categories = [Int[] for _ in 1:length(uniques)]
        for i in 1:n
            push!(categories[rev_uniques[symmetries(i)]], i)
        end
    end

    todelete = falses(length(categories))
    last_i = I[1]
    for j in 2:length(I)
        i = I[j]
        if csequences[i] == csequences[last_i]
            todelete[i] = true
            append!(categories[last_i], categories[i])
            push!(unique_reprs[last_i], unique_reprs[i][1])
            # We asserted that unique_reprs[i] only had one element before
            # Note that unique_reprs[last_i] has more than one element now
        else
            last_i = i
        end
    end
    deleteat!(categories, todelete)
    deleteat!(unique_reprs, todelete)
    deleteat!(csequences, todelete)

    num = length(categories)
    numoutgoingedges = Vector{Tuple{Int,Vector{Int}}}(undef, num)
    for i in 1:num
        seq = csequences[i]
        numoutgoingedges[i] = (length(categories[i]) * seq[1], seq)
    end
    sortorder = sortperm(numoutgoingedges)
    # categories are sorted by the total number of outgoing directed edges from each category
    # and by the coordination sequence in case of ex-aequo.
    # categories are thus uniquely determined and ordered independently of the representation of the net

    @toggleassert allunique(csequences)
    # for i in 1:num
    #     @toggleassert all(coordination_sequence(graph, x, 10) == csequences[i] for x in categories[i])
    # end # if enabled, these assertions are somewhat costly (up to ~10% total execution time)
    return categories[sortorder], unique_reprs[sortorder]
end


function check_symmetry_with_collisions(collisions)
    function check_symmetry(pge::PeriodicGraphEmbedding{D,T}, t::SVector{D,T}, r, vtypes) where {D,T}
        vmap = check_valid_symmetry(pge, t, r, vtypes, true)
        vmap isa Nothing && return nothing
        n = length(pge.pos)
        m = n - length(collisions)
        for (i, node) in enumerate(collisions)
            j = vmap[m+i].v
            m + i == j && continue
            j ≤ m && return nothing # a non-collision node is mapped to a collision node
            collisions[j-m] == CollisionNode(node, [x.v for x in vmap]) || return nothing
        end
        return vmap
    end
    check_symmetry
end


"""
    find_candidates(net::CrystalNet{D}, collisions::Vector{CollisionNode}) where D

Return a non-empty set of candidates `u => basis` where `u` is a vertex and `basis` is
matrix whose columns are `D` linearly independent euclidean embeddings of edges.
The returned set is independent of the representation of the graph used in `net`.

Also return a `category_map` linking each vertex to its category number, as defined by
[`partition_by_coordination_sequence`](@ref)

See also: [`candidate_key`](@ref)
"""
function find_candidates(net::CrystalNet{D,T}, collisions::Vector{CollisionNode}) where {D,T}
    L = D*D
    if D == 3
        check_symmetry = check_symmetry_with_collisions(collisions)
        symmetries = find_symmetries(net.pge, net.types, check_symmetry)
        categories, unique_reprs = partition_by_coordination_sequence(net.pge.g, symmetries)
    else
        categories, unique_reprs = partition_by_coordination_sequence(net.pge.g)
    end

    category_map = Vector{Int}(undef, nv(net.pge.g))
    for (i, cat) in enumerate(categories)
        for j in cat
            category_map[j] = i
        end
    end
    # @toggleassert sort.(categories) == sort.(first(partition_by_coordination_sequence(net.pge.g)))
    candidates = Dict{Int,Vector{SMatrix{D,D,T,L}}}()
    for reprs in unique_reprs
        # First, we try to look for triplet of edges all starting from the same vertex within a category
        degree(net.pge.g, first(reprs)) <= D && continue
        candidates = find_candidates_onlyneighbors(net, reprs, category_map)
        !isempty(candidates) && break
    end
    if D >= 3 && isempty(candidates)
        # If we arrive at this point, it means that all vertices only have coplanar neighbors
        # and D >= 3
        # Then we look for all triplets of edges two of which start from a vertex
        # and one does not, stopping at the first pair of categories for which this
        # set of candidates is non-empty.
        if isempty(candidates)
            for reprs in unique_reprs
                candidates = find_candidates_fallback(net, reprs, categories, category_map)
                !isempty(candidates) && break
            end
        end
    end
    if isempty(candidates)
        check_dimensionality(net)
        error("Internal error: no candidate found.")
    end
    if D == 3
        return extract_through_symmetry(candidates, symmetries), category_map
    else
        flattened_candidates = Pair{Int,SMatrix{D,D,T,L}}[]
        for (i, mats) in candidates
            for mat in mats
                push!(flattened_candidates, i => SMatrix{D,D,T,L}(mat))
            end
        end
        return flattened_candidates, category_map
    end
end


"""
    extract_through_symmetry(candidates::Dict{Int,Vector{SMatrix{3,3,T,9}}}, symmetries::AbstractSymmetryGroup) where T

Given the candidates and the list of symmetries of the net, return the flattened list of
candidates after removing candidates that are symmetric images of the kept ones.
"""
function extract_through_symmetry(candidates::Dict{Int,Vector{SMatrix{3,3,T,9}}}, symmetries::AbstractSymmetryGroup) where T
    unique_candidates = Pair{Int,SMatrix{3,3,T,9}}[]
    for (i, mats) in candidates
        @toggleassert i == symmetries(i)
        symms = [symm for symm in symmetries if symm[i] == i] # vmaps for which i is a fixpoint
        min_mats = Set{SVector{9,T}}()
        for mat in mats
            # At this point, `i => mat` is a potential candidate, and so are all its
            # symmetric images. Among those, only the ones of the form `i => mat2` (i.e.
            # with the same origin `i`) have not been eliminated earlier by symmetry.
            # Let's enumerate all such candidates and only keep one of them, for example
            # the lexicographically smallest, since they are all equivalent to `i => mat`.
            min_mat = SVector{9,T}(mat) # flattened to allow lexicographic comparison
            for symm in symms
                new_mat = SVector{9,T}(symm(mat))
                if new_mat < min_mat
                    min_mat = new_mat
                end
            end
            push!(min_mats, min_mat) # min_mats is a Set to eliminate duplicates
        end
        for x in min_mats
            push!(unique_candidates, i => SMatrix{3,3,T,9}(x))
        end
    end
    return unique_candidates
end


"""
    find_initial_candidates(net::CrystalNet{D}, candidates_v, category_map) where D

Given the net, a list of vertices in a given category and the `category_map`, return a
list of pairs `u => (basis, cats)` where `u ∈ candidates_v`, `basis` is a `D`-rank matrix
made by juxtaposing the euclidean embeddings of outgoing edges from `u`, and `cats` are the
categories of the respective neighbors of `u`.

If the `basis` corresponding to vertex `u` is not of rank `D`, it is not included in the
returned list (for instance, if all outgoing edges of a vertex are coplanar with `D == 3`).
"""
function find_initial_candidates(net::CrystalNet{D,T}, candidates_v, category_map) where {D,T}
    @toggleassert 1 ≤ D ≤ 3
    deg = degree(net.pge.g, first(candidates_v)) # The degree is the same for all vertices of the same category
    n = length(candidates_v)
    _initial_candidates = Vector{Pair{Int,Tuple{Matrix{T},Vector{Int}}}}(undef, n)
    valid_initial_candidates = falses(n)

    @threads for i in 1:n
        v = candidates_v[i]
        a = Matrix{T}(undef, D, deg)
        cats = Vector{Int}(undef, deg)
        posi = net.pge.pos[v]
        for (j, x) in enumerate(neighbors(net.pge.g, v))
            a[:,j] .= net.pge.pos[x.v] .+ x.ofs .- posi
            cats[j] = category_map[x.v]
        end
        if (D == 3 && isrank3(a)) || (D == 2 && isrank2(a)) || (D == 1 && ((@toggleassert a[:,1] != 0); true))
            _initial_candidates[i] = (v => (a, cats))
            valid_initial_candidates[i] = true
        end
    end
    return [_initial_candidates[i] for i in 1:n if valid_initial_candidates[i]]
end


"""
    find_candidates_onlyneighbors(net::CrystalNet{D}, candidates_v, category_map) where D

Given the net, a list of vertices in a given category and the `category_map`, return a Dict
whose pairs `u => matv` are such that `u ∈ candidates_v` and `matv` is a list of unique
invertible matrices of size `D` whose columns are euclidean embeddings of outgoing edges
from `u`.
Each such matrix has a category, defined by the `D`-uplet of categories of each
corresponding outneighbor of `u`: the returned Dict is such that all the matrices belonging
to all `matv` share the same category.

The returned Dict is empty iff `find_initial_candidates(net, candidates_v, category_map)`
is empty.
"""
function find_candidates_onlyneighbors end


function find_candidates_onlyneighbors(net::CrystalNet3D{T}, candidates_v, category_map) where T
    initial_candidates = find_initial_candidates(net, candidates_v, category_map)
    candidates = Dict{Int,Vector{SMatrix{3,3,T,9}}}()
    isempty(initial_candidates) && return candidates
    mincats::SVector{3,Int} = SVector{3,Int}((length(category_map), length(category_map), length(category_map)))
    ordertype::Int = 1
    # ordertype designates the kind of ordering of the category of the three edges:
    # ordertype == 1 means c1 == c2 == c3 where ci is the category of edge i
    # ordertype == 2 means c1  < c2 == c3
    # ordertype == 3 means c1 == c2  < c3
    # ordertype == 4 means c1  < c2  < c3

    for (v, (mat, cats)) in initial_candidates
        _, n = size(mat)
        matv = SMatrix{3,3,T,9}[] # the list of bases making a candidate with origin v
        for _i in 1:(n-2), _j in (_i+1):(n-1), _k in (_j+1):n
            m = SMatrix{3,3,T,9}(mat[:,[_i,_j,_k]])
            issingular(m) && continue
            orders = SVector{3,Int}[]
            subcats = cats[SVector{3,Int}(_i, _j, _k)]
            reorder = SVector{3,Int}(begin
                cmp23 = subcats[2] >= subcats[3]
                cmp13 = subcats[1] >= subcats[3]
                if subcats[1] >= subcats[2]
                    cmp23 ? (3,2,1) : (cmp13 ? (2,3,1) : (2,1,3))
                else
                    cmp23 ? (cmp13 ? (3,1,2) : (1,3,2)) : (1,2,3)
                end
            end)
            subcats = subcats[reorder]
            # subcats are the ordered categories of the three edges
            if subcats[1] == subcats[3] # c1 == c2 == c3
                (ordertype != 1 || mincats[1] < subcats[1]) && continue
                if mincats[1] > subcats[1]
                    empty!(matv)
                    empty!(candidates)
                    mincats = subcats
                end
                orders = SVector{3,Int}[(1,2,3), (1,3,2), (2,1,3), (2,3,1), (3,1,2), (3,2,1)]
            elseif subcats[2] == subcats[3] # c1 < c2 == c3
                (ordertype > 2 || (ordertype == 2 && mincats < subcats)) && continue
                if ordertype == 1 || mincats > subcats
                    empty!(matv)
                    empty!(candidates)
                    ordertype = 2
                    mincats = subcats
                end
                orders = SVector{3,Int}[reorder, (reorder[1], reorder[3], reorder[2])]
            elseif subcats[1] == subcats[2] # c1 == c2 < c3
                (ordertype == 4 || (ordertype == 3 && mincats < subcats)) && continue
                if ordertype <= 2 || mincats > subcats
                    empty!(matv)
                    empty!(candidates)
                    ordertype = 3
                    mincats = subcats
                end
                orders = SVector{3,Int}[reorder, (reorder[2], reorder[1], reorder[3])]
            else # c1 < c2 < c3
                ordertype == 4 && mincats < subcats && continue
                if ordertype != 4 || mincats > subcats
                    empty!(matv)
                    empty!(candidates)
                    ordertype = 4
                    mincats = subcats
                end
                orders = SVector{3,Int}[reorder]
            end
            append!(matv, m[:,o] for o in orders)
        end
        candidates[v] = matv
    end

    return candidates
end

function find_candidates_onlyneighbors(net::CrystalNet2D{T}, candidates_v, category_map) where T
    initial_candidates = find_initial_candidates(net, candidates_v, category_map)

    candidates = Dict{Int,Vector{SMatrix{2,2,T,4}}}()
    isempty(initial_candidates) && return candidates
    mincats::SVector{2,Int} = SVector{2,Int}((length(category_map), length(category_map)))
    ordertype::Int = 1

    for (v, (mat, cats)) in initial_candidates
        _, n = size(mat)
        matv = SMatrix{2,2,T,4}[]
        for _i in 1:(n-1), _j in (_i+1):n
            m = SMatrix{2,2,T,4}(mat[:,[_i,_j]])
            issingular(m) && continue
            orders = SVector{2,Int}[]
            subcats = cats[SVector{2,Int}(_i, _j)]
            reorder = SVector{2,Int}(subcats[1] >= subcats[2] ? (2,1) : (1,2))
            subcats = subcats[reorder]
            if subcats[1] == subcats[2]
                (ordertype != 1 || mincats[1] < subcats[1]) && continue
                if mincats[1] > subcats[1]
                    empty!(matv)
                    empty!(candidates)
                    mincats = subcats
                end
                orders = SVector{2,Int}[(1,2), (2,1)]
            else
                ordertype == 2 && mincats < subcats && continue
                if ordertype != 2 || mincats > subcats
                    empty!(matv)
                    empty!(candidates)
                    ordertype = 2
                    mincats = subcats
                end
                orders = SVector{2,Int}[reorder]
            end
            append!(matv, m[:,o] for o in orders)
        end
        candidates[v] = matv
    end

    return candidates
end

function find_candidates_onlyneighbors(net::CrystalNet1D{T}, candidates_v, category_map) where T
    initial_candidates = find_initial_candidates(net, candidates_v, category_map)

    candidates = Dict{Int,Vector{SMatrix{1,1,T,1}}}()
    isempty(initial_candidates) && return candidates
    mincat::Int = length(category_map)

    for (v, (mat, cats)) in initial_candidates
        n = length(mat)
        matv = SMatrix{1,1,T,1}[]
        for i in 1:n
            cat = cats[i]
            mincat < cat && continue
            if mincat > cat
                empty!(matv)
                empty!(candidates)
                mincat = cat
            end
            push!(matv, mat[[i]])
        end
        candidates[v] = matv
    end

    return candidates
end


"""
    find_candidates_fallback(net::CrystalNet3D, reprs, othercats, category_map)

Return candidates in the same form as [`find_candidates_onlyneighbors`](@ref) except that
only two edges start from `u` and one does not.
"""
function find_candidates_fallback(net::CrystalNet3D{T}, reprs, othercats, category_map) where T
    candidates = Dict{Int,Vector{SMatrix{3,3,T,9}}}(u => [] for u in reprs)
    n = length(category_map)
    mincats = SizedVector{3,Int}(fill(n, 3))
    current_cats = SizedVector{3,Int}(fill(0, 3))
    for cat in othercats
        for u in reprs
            mats = SMatrix{3,3,T,9}[]
            posu = net.pge.pos[u]
            neighu = neighbors(net.pge.g, u)
            for (j1, _x1) in enumerate(neighu)
                category_map[_x1.v] > mincats[1] && category_map[_x1.v] > mincats[2] && continue
                for j2 in j1+1:length(neighu)
                    x2 = neighu[j2]
                    x1 = _x1
                    if category_map[x1.v] > category_map[x2.v]
                        x1, x2 = x2, x1
                    end
                    current_cats[1] = category_map[x1.v]
                    current_cats[1] > mincats[1] && continue
                    current_cats[2] = category_map[x2.v]
                    current_cats[1] == mincats[1] && current_cats[2] > mincats[2] && continue
                    vec1 = net.pge.pos[x1.v] .+ x1.ofs .- posu
                    vec2 = net.pge.pos[x2.v] .+ x2.ofs .- posu
                    j = iszero(vec1[1]) ? iszero(vec1[2]) ? 3 : 2 : 1
                    vec1[j] * vec2 == vec2[j] * vec1 && continue # colinear
                    for v in cat
                        posv = net.pge.pos[v]
                        for x3 in neighbors(net.pge.g, v)
                            vec3 = net.pge.pos[x3.v] .+ x3.ofs .- posv
                            mat = SMatrix{3,3,T,9}(hcat(vec1, vec2, vec3))
                            issingular(mat) && continue
                            current_cats[3] = category_map[x3.v]
                            current_cats > mincats && continue
                            if current_cats < mincats
                                mincats .= current_cats
                                foreach(empty!, values(candidates))
                                empty!(mats)
                            end
                            # if d < 0
                            #     mat = hcat(vec2, vec1, vec3)
                            # end
                            push!(mats, mat)
                            if current_cats[1] == current_cats[2]
                                push!(mats, hcat(vec2, vec1, vec3))
                            end
                        end
                    end
                end
            end
            append!(candidates[u], mats)
        end
        # If at this point candidates is not empty, then we have found the first
        # other category cat such that there is a candidate with two edges from
        # maincat and one from cat. This is uniquely determined independently of
        # the representation of the net, so we can simply return the found candidates
        # without needing to check for the remaining categories.
        if !all(isempty, values(candidates))
            for (u, mats) in candidates
                if isempty(mats)
                    delete!(candidates, u)
                end
            end
            return candidates
         end
    end
    @toggleassert all(isempty, values(candidates))
    return Dict{Int,Vector{SMatrix{3,3,T,9}}}()
end


"""
    topological_key(net::CrystalNet)

Return a unique topological key for the net, which is a topological invariant of the net
(i.e. it does not depend on its initial representation).
"""
function topological_key(net::CrystalNet{D}) where D
    isempty(net.pge.pos) && return PeriodicGraph{D}()
    newnet, collisions = collision_nodes(net)
    if collisions isa Nothing
        net.pge.g.width[] = -2 # internal error code
        return net.pge.g
    end
    return topological_key(newnet, collisions::Vector{CollisionNode})
end

function topological_key(net::CrystalNet{D,T}, collisions::Vector{CollisionNode}) where {D,T}
    candidates, category_map = find_candidates(net, collisions)
    v, minimal_basis = popfirst!(candidates)
    n = nv(net.pge.g)
    _edgs = [(n + 1, 0, zero(SVector{D,T}))]
    minimal_vmap, minimal_edgs = candidate_key(net, v, minimal_basis, _edgs)
    for (v, basis) in candidates
        vmap, edgs = candidate_key(net, v, basis, minimal_edgs)
        isempty(vmap) && continue
        @toggleassert edgs < minimal_edgs
        if edgs < minimal_edgs
            minimal_edgs = edgs
            minimal_vmap = vmap
            minimal_basis = basis
        end
    end

    newbasis, edges = findbasis(minimal_edgs)
    graph = PeriodicGraph{D}(n, edges)

    if !isempty(collisions)
        graph = expand_collisions(collisions, graph, minimal_vmap)
    end

    # finalbasis = minimal_basis * newbasis
    # return Int.(finalbasis), minimal_vmap, graph
    return graph
end

