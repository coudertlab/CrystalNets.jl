struct UnstableNetException <: Exception
    g::String
end
UnstableNetException(net::CrystalNet) = UnstableNetException(string("unstable ", net.graph))
function Base.showerror(io::IO, e::UnstableNetException)
    print(io, "The net is unstable, hence it cannot be analyzed.")
end

function check_dimensionality(c::CrystalNet{D,T}) where {D,T}
    edgs = [c.pos[dst(e)] .+ PeriodicGraphs.ofs(e) .- c.pos[src(e)] for e in edges(c.graph)]
    sort!(edgs)
    unique!(edgs)
    mat = reduce(hcat, edgs)
    if D == 3 && !isrank3(mat)
        throw(ArgumentError("the input net does not have dimensionality 3, so it cannot be analyzed."))
    elseif D == 2 && !isrank2(mat)
        throw(ArgumentError("the input net does not have dimensionality 2, so it cannot be analyzed."))
    elseif D == 1 && !isrank1(mat)
        throw(ArgumentError("the input net does not have dimensionality 1, so it cannot be analyzed."))
    else
        throw(AssertionError("1 ≤ D ≤ 3"))
    end
    nothing
end

function check_valid_translation(c::CrystalNet{D,T}, t::SVector{D,<:Rational{<:Integer}}, r=nothing) where {D,T}
    U = soft_widen(T)
    n = length(c.pos)
    vmap = Int[]
    offsets = SVector{D,Int}[]
    for k in 1:n
        curr_pos = c.pos[k]
        transl = U.(isnothing(r) ? curr_pos : (r * curr_pos)) .+ t
        ofs = floor.(Int, transl)
        x = transl .- ofs
        (i, j) = x < curr_pos ? (1, k) : (k, length(c.types)+1)
        while j - i > 1
            m = (j+i)>>1
            if cmp(x, c.pos[m]) < 0
                j = m
            else
                i = m
            end
        end
        (c.pos[i] == x && c.types[i] == c.types[k]) || return nothing
        push!(vmap, i)
        push!(offsets, ofs)
    end
    for e in edges(c.graph)
        src = vmap[e.src]
        dst = vmap[e.dst.v]
        newofs = (isnothing(r) ? e.dst.ofs : r * e.dst.ofs) .+ offsets[e.dst.v] .- offsets[e.src]
        has_edge(c.graph, PeriodicGraphs.unsafe_edge{D}(src, dst, newofs)) || return nothing
    end
    return vmap
end

function possible_translations(c::CrystalNet{D,T}) where {D,T}
    n = length(c.pos)
    ts = Tuple{Int, Int, Int, SVector{D,T}}[]
    sortedpos = copy(c.pos)
    origin = popfirst!(sortedpos)
    @assert iszero(origin)
    sort!(SVector{D,soft_widen(widen(widen(T)))}.(sortedpos), by=norm)
    for t in sortedpos
        @assert t == back_to_unit.(t)
        max_den, i_max_den = findmax(denominator.(t))
        numerator(t[i_max_den]) == 1 || continue
        nz = count(iszero, t)
        push!(ts, (nz, i_max_den, max_den, t))
    end
    return sort!(ts; by=(x->(x[1], x[2], x[3])))
end

# TODO: remove?
function find_first_valid_translations(c::CrystalNet)
    for (nz, i_max_den, max_den, t) in possible_translations(c)
        !isnothing(check_valid_translation(c, t)) && return (nz, i_max_den, max_den, t)
    end
    return nothing
end

function find_all_valid_translations(c::CrystalNet{D,T}) where {D,T}
    ret = NTuple{D, Vector{Tuple{Int, Int, SVector{D,T}}}}(ntuple(_->[], Val(D)))
    for (nz, i_max_den, max_den, t) in possible_translations(c)
        !isnothing(check_valid_translation(c, t)) && push!(ret[nz+1], (i_max_den, max_den, t))
    end
    return ret
end


function minimal_volume_matrix(translations::Tuple{T}) where T
    nz0 = translations[1]
    denmax = 1//1
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
    detmin = 1//1
    best = 0
    @inbounds for i in 1:n
        d = nz0[i][3][1]
        d == 0 && continue
        if d < detmin
            detmin = d
            best = i
        end
    end
    @assert !iszero(best)
    ret = hcat(nz0[best][3])
    if ret[1] < 0
        ret = hcat(.-nz0[best][3])
    end
    return ret
end

function minimal_volume_matrix(translations::Tuple{T,T}) where T
    nz0, nz1 = translations
    denmax = [1//1, 1//1]
    imax = [0, 0]
    for j in 1:length(nz1)
        i, den, _ = nz1[j]
        if den > denmax[i]
            imax[i] = j
            denmax[j] = den
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
    detmin = 1//1
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
    @assert !any(iszero.(best))
    i, j = best
    ret = hcat(all[i][3], all[j][3])
    if det(ret) < 0
        ret = hcat(all[i][3], .-all[j][3])
    end
    return ret
end

function minimal_volume_matrix(translations::Tuple{T,T,T}) where T
    nz0, nz1, nz2 = translations

    denmax = [1//1, 1//1, 1//1]
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
    detmin = 1//1
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
    @assert !any(iszero.(best))
    i, j, k = best
    ret = hcat(all[i][3], all[j][3], all[k][3])
    if det(ret) < 0
        ret = hcat(all[i][3], all[j][3], .-all[k][3])
    end
    return ret
end

function reduce_with_matrix(c::CrystalNet{D,<:Rational{T}}, mat) where {D,T}
    U = widen(T)
    lengths = degree(c.graph)
    if D == 3
        cell = Cell(c.cell, c.cell.mat * mat)
    else
        _mat = SizedMatrix{3,3,BigFloat}(LinearAlgebra.I)
        _mat[1:D,1:D] .= mat
        cell = Cell(c.cell, c.cell.mat * _mat)
    end
    imat = soft_widen(T).(inv(mat)) # The inverse should only have integer coefficients
    poscol = (U.(imat),) .* c.pos

    offset = [floor.(Int, x) for x in poscol]
    for i in 1:length(poscol)
        poscol[i] = poscol[i] .- offset[i]
    end
    I_sort = sort(1:length(poscol); by=i->(poscol[i], hash_position(offset[i])))
    i = popfirst!(I_sort)
    @assert iszero(offset[i])
    I_kept = Int[i]
    sortedcol = SVector{D,Rational{U}}[poscol[i]]
    for i in I_sort
        x = poscol[i]
        if x != sortedcol[end]
            push!(I_kept, i)
            push!(sortedcol, x)
        end
    end
    edges = PeriodicEdge{D}[]
    for i in 1:length(I_kept)
        ofs_i = offset[I_kept[i]]
        for neigh in neighbors(c.graph, I_kept[i])
            x = poscol[neigh.v]
            j = searchsortedfirst(sortedcol, x)
            @assert j <= length(sortedcol) && sortedcol[j] == x
            ofs_x = offset[neigh.v]
            push!(edges, (i, j, ofs_x - ofs_i .+ imat*neigh.ofs))
        end
    end
    graph = PeriodicGraph{D}(edges)
    @assert degree(graph) == lengths[I_kept]
    return CrystalNet(cell, c.types[I_kept], sortedcol, graph)
end

function minimize(net::CrystalNet)
    translations = find_all_valid_translations(net)
    while !all(isempty.(translations))
        mat = minimal_volume_matrix(translations)
        net = reduce_with_matrix(net, mat)
        translations = find_all_valid_translations(net) # TODO don't recompute
    end
    return net
end

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

    for j in 1:length(positivetrans)
        i = tmp_map[I_sort[j]]
        trans = positivetrans[I_sort[j]]
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

    @assert all(z-> begin x, y = z; x == (y.src, y.dst.v, basis*y.dst.ofs) end, zip(edges, newedges))

    return basis, newedges
end


function candidate_key(net::CrystalNet{D,T}, u, basis, minimal_edgs) where {D,T}
    V = soft_widen(widen(T))
    n = nv(net.graph)
    h = 2 # Next node to assign
    origin = net.pos[u]
    newpos = Vector{SVector{D,V}}(undef, n) # Positions of the kept representatives
    newpos[1] = zero(SVector{D,V})
    offsets = Vector{SVector{D,Int}}(undef, n)
    # Offsets of the new representatives with respect to the original one, in the original basis
    offsets[1] = zero(SVector{D,Int})
    vmap = Vector{Int}(undef, n) # Bijection from the old to the new node number
    vmap[1] = u
    rev_vmap = zeros(Int, n) # inverse of vmap
    rev_vmap[u] = 1
    flag_bestedgs = false # Marks whether the current list of edges is lexicographically
    # below the best known one. If not by the end of the algorithm, return false
    edgs = Tuple{Int,Int,SVector{D,V}}[]
    mat = V.(inv(widen(V).(basis)))
    for t in 1:n # t is the node being processed
        neighs = neighbors(net.graph, vmap[t])
        ofst = offsets[t]
        pairs = Vector{Tuple{SVector{D,V},Int}}(undef, length(neighs))
        for (i,x) in enumerate(neighs)
            pairs[i] = (V.(mat*(net.pos[x.v] .+ x.ofs .- origin .+ ofst)), x.v)
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
        bigbasis = widen(V).(basis)
        for (coordinate, v) in pairs
            idx = rev_vmap[v]
            if idx == 0 # New node to which h is assigned
                @assert t < h
                vmap[h] = v
                rev_vmap[v] = h
                newpos[h] = coordinate
                offsets[h] = SVector{D,Int}(bigbasis * coordinate .+ origin .- net.pos[v])
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
                c < 0 && return (Int[], Tuple{Int,Int,SVector{D,V}}[])
                c > 0 && (flag_bestedgs = true)
            end
        end
    end
    @assert allunique(edgs)
    if !flag_bestedgs
        @assert minimal_edgs == edgs
        return (Int[], Tuple{Int,Int,SVector{D,V}}[])
    end
    return vmap, edgs
end


"""
Partition the vertices of the graph into disjoint categories, one for each
coordination sequence. The partition is then sorted by order of coordination sequence.
This partition does not depend on the representation of the graph.
The optional argument vmaps is a set of permutations of the vertices that leave the
graph unchanged. In other words, vmaps is a set of symmetry operations of the graph.

Also returns the map from vertices to an identifier such that all vertices with
the same identifier are symmetric images of one another, as well as a list of
unique representative for each symmetry class.
"""
function partition_by_coordination_sequence(graph, vmaps=nothing)
    # First, union-find on the symmetries to avoid computing useless coordination sequences
    n = nv(graph)
    unionfind = collect(1:n)
    if vmaps isa Nothing || isempty(vmaps) # Absence of symmetry
        unique_reprs = Vector{Int}[Int[i] for i in 1:n]
        cat_map = Int[i for i in 1:n]
        categories = Vector{Int}[[i] for i in 1:n]
    else
        for vmap in vmaps
            for (i,j) in enumerate(vmap)
                i == j && continue
                if i > j
                    i, j = j, i
                end
                repri = i
                while repri != unionfind[repri]
                    repri = unionfind[repri]
                end
                reprj = j
                while reprj != unionfind[reprj]
                    reprj = unionfind[reprj]
                end
                repri == reprj && continue
                if repri < reprj
                    unionfind[j] = repri
                else
                    unionfind[i] = reprj
                end
            end
        end

        cat_map = zeros(Int, n)
        unique_reprs = Vector{Int}[]
        categories = Vector{Int}[]
        for i in 1:n
            repri = i
            if repri != unionfind[repri]
                repri = unionfind[repri]
                while repri != unionfind[repri]
                    repri = unionfind[repri]
                end
                descent = i
                while descent != unionfind[descent]
                    tmp = unionfind[descent]
                    unionfind[descent] = repri
                    descent = tmp
                end
                @assert descent == repri
            end
            cat = cat_map[repri]
            if iszero(cat)
                push!(unique_reprs, [i])
                push!(categories, [i])
                cat = length(categories)
                cat_map[repri] = cat
            else
                push!(categories[cat], i)
            end
            cat_map[i] = cat
        end
    end

    PeriodicGraphs.graph_width!(graph) # setup for the computation of coordination sequences
    csequences = [coordination_sequence(graph, first(cat), 10) for cat in categories]
    I = sortperm(csequences)
    @assert csequences[I[1]][1] >= 2 # vertices of degree <= 1 should have been removed at input creation

    todelete = falses(length(categories))
    last_i = I[1]
    for j in 2:length(I)
        i = I[j]
        if csequences[i] == csequences[last_i]
            todelete[i] = true
            append!(categories[last_i], categories[i])
            push!(unique_reprs[last_i], pop!(unique_reprs[i]))
            empty!(categories[i])
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

    @assert allunique(csequences)
    for i in 1:num
        @assert all(coordination_sequence(graph, x, 10) == csequences[i] for x in categories[i])
    end # TODO comment out these costly asserts
    return categories[sortorder], cat_map, unique_reprs[sortorder]
end


"""
Find a set of pairs (v, basis) that are candidates for finding the key with
candidate_key(net, v, basis).
The returned set is independent of the representation of the graph used in net.
"""
function find_candidates(net::CrystalNet{D,T}) where {D,T}
    L = D*D
    U = soft_widen(T)
    if D == 3
        rotations, vmaps, _ = find_symmetries(net)
        @assert length(rotations) == length(vmaps)
        categories, symmetry_map, unique_reprs = partition_by_coordination_sequence(net.graph, vmaps)
    else
        categories, symmetry_map, unique_reprs = partition_by_coordination_sequence(net.graph)
    end

    category_map = Vector{Int}(undef, nv(net.graph))
    for (i, cat) in enumerate(categories)
        for j in cat
            category_map[j] = i
        end
    end
    @assert sort.(categories) == sort.(first(partition_by_coordination_sequence(net.graph)))
    candidates = Dict{Int,Vector{SMatrix{D,D,U,L}}}()
    for reprs in unique_reprs
        # First, we try to look for triplet of edges all starting from the same vertex within a category
        degree(net.graph, first(reprs)) <= D && continue
        candidates = find_candidates_onlyneighbors(net, reprs, category_map)
        !isempty(candidates) && break
    end
    if D >= 3 && isempty(candidates)
        # If we arrive at this point, it means that all vertices only have coplanar neighbors
        # (in the case D == 3, it should be impossible for D < 3).
        # Then we look for all triplets of edges two of which start from a vertex
        # and one does not, stopping and the first pair of categories for which this
        # set of candidates is non-empty.
        if isempty(candidates)
            for reprs in unique_reprs
                candidates = find_candidates_fallback(net, reprs, categories, category_map)
                !isempty(candidates) && break
            end
        end
    end
    if isempty(candidates)
        if !allunique(net.pos)
            throw(UnstableNetException(net))
        else
            check_dimensionality(net)
            error("Internal error: no candidate found.")
        end
    end
    if D == 3
        return extract_through_symmetry(candidates, vmaps, rotations), category_map
    else
        flattened_candidates = Pair{Int,SMatrix{D,D,U,L}}[]
        for (i, mats) in candidates
            for mat in mats
                push!(flattened_candidates, i => SMatrix{D,D,U,L}(mat))
            end
        end
        return flattened_candidates, category_map
    end
end

function extract_through_symmetry(candidates::Dict{Int,Vector{SMatrix{3,3,T,9}}}, vmaps, rotations) where T
    unique_candidates = Pair{Int,SMatrix{3,3,T,9}}[]
    for (i, mats) in candidates
        @assert i == minimum(vmap[i] for vmap in vmaps)
        indices = [j for j in 1:length(vmaps) if vmaps[j][i] == i]
        min_mats = Set{SVector{9,T}}()
        for mat in mats
            flattenedmat = SVector{9,T}(mat)
            min_mat = flattenedmat
            for j in indices
                new_mat = SVector{9,T}(rotations[j] * mat)
                if new_mat < min_mat
                    min_mat = new_mat
                end
            end
            push!(min_mats, min_mat)
        end
        for x in min_mats
            push!(unique_candidates, i => SMatrix{3,3,T,9}(x))
        end
    end
    return unique_candidates
end

function find_candidates_onlyneighbors(net::CrystalNet3D{T}, candidates_v, category_map) where T
    U = soft_widen(T)
    deg = degree(net.graph, first(candidates_v)) # The degree is the same for all vertices of the same category
    initial_candidates = Pair{Int,Tuple{Matrix{U},Vector{Int}}}[]
    fastlock = SpinLock()
    @threads for v in candidates_v
        a = Matrix{U}(undef, 3, deg)
        cats = Vector{Int}(undef, deg)
        posi = net.pos[v]
        for (j, x) in enumerate(neighbors(net.graph, v))
            a[:,j] .= net.pos[x.v] .+ x.ofs .- posi
            cats[j] = category_map[x.v]
        end
        if isrank3(a)
            pair = v => (a, cats)
            lock(fastlock) do
                push!(initial_candidates, pair)
            end
        end
    end

    candidates = Dict{Int,Vector{SMatrix{3,3,U,9}}}()
    isempty(initial_candidates) && return candidates
    current_cats::SizedVector{3,Int} = SizedVector{3,Int}(fill(length(category_map), 3))
    current_ordertype::Int = 1

    @threads for (v, (mat, cats)) in initial_candidates
        _, n = size(mat)
        matv = SMatrix{3,3,U,9}[]
        lock(fastlock)
        ordertype = current_ordertype
        mincats = copy(current_cats)
        unlock(fastlock)
        for _i in 1:(n-2), _j in (_i+1):(n-1), _k in (_j+1):n
            m = SMatrix{3,3,U,9}(mat[:,[_i,_j,_k]])
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
            if subcats[1] == subcats[3]
                (ordertype != 1 || mincats[1] < subcats[1]) && continue
                if mincats[1] > subcats[1]
                    empty!(matv)
                    mincats = subcats
                end
                orders = SVector{3,Int}[(1,2,3), (1,3,2), (2,1,3), (2,3,1), (3,1,2), (3,2,1)]
            elseif subcats[2] == subcats[3]
                (ordertype > 2 || (ordertype == 2 && mincats < subcats)) && continue
                if ordertype == 1 || mincats > subcats
                    empty!(matv)
                    ordertype = 2
                    mincats = subcats
                end
                orders = SVector{3,Int}[reorder, (reorder[1], reorder[3], reorder[2])]
            elseif subcats[1] == subcats[2]
                (ordertype == 4 || (ordertype == 3 && mincats < subcats)) && continue
                if ordertype <= 2 || mincats > subcats
                    empty!(matv)
                    ordertype = 3
                    mincats = subcats
                end
                orders = SVector{3,Int}[reorder, (reorder[2], reorder[1], reorder[3])]
            else
                ordertype == 4 && mincats < subcats && continue
                if ordertype != 4 || mincats > subcats
                    empty!(matv)
                    ordertype = 4
                    mincats = subcats
                end
                orders = SVector{3,Int}[reorder]
            end
            mats = [m[:,o] for o in orders]
            append!(matv, mats)
        end
        if ordertype < current_ordertype || (ordertype == current_ordertype
                                             && mincats > current_cats)
            continue # fast-path to avoid unecessary locking
        end

        lock(fastlock)
        if ordertype < current_ordertype || (ordertype == current_ordertype
                                             && mincats > current_cats)
            unlock(fastlock)
            continue
        end
        if ordertype > current_ordertype || mincats < current_cats
            empty!(candidates)
            current_ordertype = ordertype
            current_cats = mincats
        end
        candidates[v] = matv
        unlock(fastlock)
    end

    return candidates
end

function find_candidates_onlyneighbors(net::CrystalNet2D{T}, candidates_v, category_map) where T
    U = soft_widen(T)
    deg = degree(net.graph, first(candidates_v)) # The degree is the same for all vertices of the same category
    initial_candidates = Pair{Int,Tuple{Matrix{U},Vector{Int}}}[]
    fastlock = SpinLock()
    @threads for v in candidates_v
        a = Matrix{U}(undef, 2, deg)
        cats = Vector{Int}(undef, deg)
        posi = net.pos[v]
        for (j, x) in enumerate(neighbors(net.graph, v))
            a[:,j] .= net.pos[x.v] .+ x.ofs .- posi
            cats[j] = category_map[x.v]
        end
        if isrank2(a)
            pair = v => (a, cats)
            lock(fastlock) do
                push!(initial_candidates, pair)
            end
        end
    end

    candidates = Dict{Int,Vector{SMatrix{2,2,U,4}}}()
    isempty(initial_candidates) && return candidates
    current_cats::SizedVector{2,Int} = SizedVector{2,Int}(fill(length(category_map), 2))
    current_ordertype::Int = 1

    @threads for (v, (mat, cats)) in initial_candidates
        _, n = size(mat)
        matv = SMatrix{2,2,U,4}[]
        lock(fastlock)
        ordertype = current_ordertype
        mincats = copy(current_cats)
        unlock(fastlock)
        for _i in 1:(n-1), _j in (_i+1):n
            m = SMatrix{2,2,U,4}(mat[:,[_i,_j]])
            issingular(m) && continue
            orders = SVector{2,Int}[]
            subcats = cats[SVector{2,Int}(_i, _j)]
            reorder = SVector{2,Int}(subcats[1] >= subcats[2] ? (2,1) : (1,2))
            subcats = subcats[reorder]
            if subcats[1] == subcats[2]
                (ordertype != 1 || mincats[1] < subcats[1]) && continue
                if mincats[1] > subcats[1]
                    empty!(matv)
                    mincats = subcats
                end
                orders = SVector{2,Int}[(1,2), (2,1)]
            else
                ordertype == 2 && mincats < subcats && continue
                if ordertype != 2 || mincats > subcats
                    empty!(matv)
                    ordertype = 2
                    mincats = subcats
                end
                orders = SVector{2,Int}[reorder]
            end
            mats = [m[:,o] for o in orders]
            append!(matv, mats)
        end
        if ordertype < current_ordertype || (ordertype == current_ordertype
                                             && mincats > current_cats)
            continue # fast-path to avoid unecessary locking
        end

        lock(fastlock)
        if ordertype < current_ordertype || (ordertype == current_ordertype
                                             && mincats > current_cats)
            unlock(fastlock)
            continue
        end
        if ordertype > current_ordertype || mincats < current_cats
            empty!(candidates)
            current_ordertype = ordertype
            current_cats = mincats
        end
        candidates[v] = matv
        unlock(fastlock)
    end

    return candidates
end

function find_candidates_onlyneighbors(net::CrystalNet1D{T}, candidates_v, category_map) where T
    U = soft_widen(T)
    deg = degree(net.graph, first(candidates_v)) # The degree is the same for all vertices of the same category
    initial_candidates = Pair{Int,Tuple{Vector{U},Vector{Int}}}[]
    fastlock = SpinLock()
    @threads for v in candidates_v
        a = Vector{U}(undef, deg)
        cats = Vector{Int}(undef, deg)
        posi = net.pos[v][]
        for (j, x) in enumerate(neighbors(net.graph, v))
            diff = net.pos[x.v][] + x.ofs[] - posi
            @assert !iszero(diff)
            a[j] = diff
            cats[j] = category_map[x.v]
        end
        pair = v => (a, cats)
        lock(fastlock) do
            push!(initial_candidates, pair)
        end
    end

    candidates = Dict{Int,Vector{SMatrix{1,1,U,1}}}()
    isempty(initial_candidates) && return candidates
    current_cat::Int = length(category_map)
    current_ordertype::Int = 1

    @threads for (v, (mat, cats)) in initial_candidates
        n = length(mat)
        matv = SMatrix{1,1,U,1}[]
        lock(fastlock)
        ordertype = current_ordertype
        mincat = current_cat
        unlock(fastlock)
        for i in 1:n
            cat = cats[i]
            mincat < cat && continue
            if mincat > cat
                empty!(matv)
                mincat = cat
            end
            append!(matv, [mat[[i]]])
        end
        if ordertype < current_ordertype || (ordertype == current_ordertype
                                             && mincat > current_cat)
            continue # fast-path to avoid unecessary locking
        end

        lock(fastlock)
        if ordertype < current_ordertype || (ordertype == current_ordertype
                                             && mincat > current_cat)
            unlock(fastlock)
            continue
        end
        if ordertype > current_ordertype || mincat < current_cat
            empty!(candidates)
            current_ordertype = ordertype
            current_cat = mincat
        end
        candidates[v] = matv
        unlock(fastlock)
    end

    return candidates
end

function find_candidates_fallback(net::CrystalNet{3,T}, reprs, othercats, category_map) where T
    U = soft_widen(T)
    candidates = Dict{Int,Vector{SMatrix{3,3,U,9}}}(u => [] for u in reprs)
    n = length(category_map)
    mincats = SizedVector{3,Int}(fill(n, 3))
    current_cats = SizedVector{3,Int}(fill(0, 3))
    for cat in othercats
        for u in reprs
            mats = SMatrix{3,3,T,9}[]
            posu = net.pos[u]
            neighu = neighbors(net.graph, u)
            for (j1, x1) in enumerate(neighu)
                category_map[x1.v] > mincats[1] && category_map[x1.v] > mincats[2] && continue
                for j2 in j1+1:length(neighu)
                    x2 = neighu[j2]
                    if category_map[x1.v] > category_map[x2.v]
                        x1, x2 = x2, x1
                    end
                    current_cats[1] = category_map[x1.v]
                    current_cats[1] > mincats[1] && continue
                    current_cats[2] = category_map[x2.v]
                    current_cats[1] == mincats[1] && current_cats[2] > mincats[2] && continue
                    vec1 = net.pos[x1.v] .+ x1.ofs .- posu
                    vec2 = net.pos[x2.v] .+ x2.ofs .- posu
                    j = iszero(vec1[1]) ? iszero(vec1[2]) ? 3 : 2 : 1
                    vec1[j] * vec2 == vec2[j] .* vec1 && continue # colinear
                    for v in cat
                        posv = net.pos[v]
                        for x3 in neighbors(net.graph, v)
                            vec3 = net.pos[x3.v] .+ x3.ofs .- posv
                            mat = SMatrix{3,3,U,9}(hcat(vec1, vec2, vec3))
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
    @assert all(isempty, values(candidates))
    return Dict{Int,Vector{SMatrix{3,3,U,9}}}()
end

# function parallel_topological_key(net::CrystalNet{D,T}) where {D,T}
#     candidates = find_candidates(net)
#     # @show length(candidates)
#     numthreads = min(nthreads(), length(candidates))
#     minimal_graph = Vector{PeriodicGraph3D}(undef, numthreads)
#     minimal_vmap = Vector{Vector{Int}}(undef, numthreads)
#     minimal_basis = Vector{SMatrix{3,3,T,9}}(undef, numthreads)
#     @threads for i in 1:numthreads
#         v, basis = candidates[end+1-i]
#         minimal_basis[i], minimal_vmap[i], minimal_graph[i] = candidate_key(net, v, basis)
#     end
#     resize!(candidates, length(candidates) - numthreads)
#     @threads for (v, basis) in candidates
#         id = threadid() # necessarily inbounds because otherwise the candidates is empty
#         newbasis, vmap, graph = candidate_key(net, v, basis)
#         if edges(graph) < edges(minimal_graph[id])
#             minimal_graph[id] = graph
#             minimal_vmap[id] = vmap
#             minimal_basis[id] = newbasis
#         end
#     end
#     _, j = findmin(edges.(minimal_graph))
#     return minimal_basis[j], minimal_vmap[j], minimal_graph[j]
# end

function topological_key(net::CrystalNet{D,T}) where {D,T}
    if isempty(net.pos)
        throw(ArgumentError("the net is empty."))
    end
    @assert allunique(net.pos) # FIXME make a more precise check for net stability. Currently fails for sxt
    candidates, category_map = find_candidates(net)
    v, minimal_basis = popfirst!(candidates)
    n = nv(net.graph)
    _edgs = [(n + 1, 0, zero(SVector{D,T}))]
    minimal_vmap, minimal_edgs = candidate_key(net, v, minimal_basis, _edgs)
    for (v, basis) in candidates
        vmap, edgs = candidate_key(net, v, basis, minimal_edgs)
        isempty(vmap) && continue
        @assert edgs < minimal_edgs
        if edgs < minimal_edgs
            minimal_edgs = edgs
            minimal_vmap = vmap
            minimal_basis = basis
        end
    end

    newbasis, coords = findbasis(minimal_edgs)

    return minimal_basis * newbasis, minimal_vmap, PeriodicGraph{D}(n, coords)
end

function find_new_representation(pos, basis, vmap, graph)
    return collect(eachcol(CIFTypes.equilibrium(graph))) # TODO optimize
end


function topological_genome(net::CrystalNet{D,T}, skipminimize=false; forgettypes=true)::String where {D,T}
    if !allunique(net.pos)
        return UnstableNetException(net).g
    end
    if forgettypes
        net = CrystalNet{D,T}(net.cell, fill(Symbol(""), length(net.types)), net.pos, net.graph)
    end
    if !skipminimize
        try
            net = minimize(net)
        catch e
            if T !== BigInt && (e isa OverflowError || e isa InexactError)
                return topological_genome(CrystalNet{D,widen(soft_widen(T))}(net), false; forgettypes=false)
            end
            rethrow()
        end
    end
    genome::String = try
        string(last(topological_key(net)))
    catch e
        if T !== BigInt && (e isa OverflowError || e isa InexactError)
            return topological_genome(CrystalNet{D,widen(soft_widen(T))}(net), true; forgettypes=false)
        end
        if e isa UnstableNetException
            e.g
        else
            rethrow()
        end
    end

    return genome

    # return CrystalNet(graph) # TODO remove this temporary bandaid
    # # FIXME CrystalNet(graph).graph != graph so this is invalid !
    #
    # positions = find_new_representation(net.pos, basis, vmap, graph)
    #
    # ret = CrystalNet(Cell(net.cell, net.cell.mat*basis), net.types[vmap], positions, graph)
    # means = mean(ret.pos)
    # for i in 1:length(ret.pos)
    #     ret.pos[i] = ret.pos[i] .+ [1//2, 1//2, 1//2] .- means
    # end
    # @assert all(ret.pos[i] == mean(ret.pos[x.v] .+ x.ofs for x in neighbors(ret.graph, i)) for i in 1:length(net.pos))
    #
    # return ret
end

function topological_genome(g::PeriodicGraph, skipminimize=true; forgettypes=true)
    net = CrystalNet(g)
    return topological_genome(net, skipminimize; forgettypes)
end

function _loop_group!(ex, id, net, group)
    for j in 1:length(ex.args)
        arg = ex.args[j]
        if arg isa Symbol
            if arg === :id
                ex.args[j] = id
            elseif arg === :net
                ex.args[j] = net
            elseif arg === :group
                ex.args[j] = group
            end
        elseif arg isa Expr
            _loop_group!(arg, id, net, group)
        end
    end
    nothing
end
macro loop_group(ex)
    ret = Expr(:block)
    for i in 1:3
        id = Symbol("id", i)
        net = Symbol("net", i)
        D = Symbol("D", i)
        group = :(group.$D)
        newex = deepcopy(ex)
        _loop_group!(newex, id, net, group)
        push!(ret.args, esc(newex))
    end
    return ret
end

function topological_genome(group::CrystalNetGroup, skipminimize=false; forgettypes=true)
    ret = Tuple{Vector{Int},String}[]
    @loop_group for (id, net) in group
        push!(ret, (id, topological_genome(net, skipminimize; forgettypes)))
    end
    return ret
end

function reckognize_topology(genome::AbstractString, arc=CRYSTAL_NETS_ARCHIVE)
    (startswith(genome, "unstable") || genome == "non-periodic") && return genome
    get(arc, genome, "UNKNOWN")
end


function reckognize_topologies(path; ignore_atoms=(), forgettypes=true)
    ret = Dict{String,String}()
    failed = Dict{String, Tuple{Exception, Vector{Union{Ptr{Nothing}, Base.InterpreterIP}}}}()
    newarc = copy(CRYSTAL_NETS_ARCHIVE)
    @threads for f in readdir(path)
        name = first(splitext(f))
        genomes::Vector{Tuple{Vector{Int},String}} = try
            topological_genome(CrystalNetGroup(parse_chemfile(path*f; ignore_atoms)); forgettypes)
        catch e
            if e isa InterruptException ||
              (e isa TaskFailedException && e.task.result isa InterruptException)
                rethrow()
            end
            if e isa NonPeriodicInputException
                [(Int[], "non-periodic")]
            else
                failed[name] = (e, catch_backtrace())
                Tuple{Vector{Int},String}[]
            end
        end
        for (i, (_, genome)) in enumerate(genomes)
            newname = length(genomes) == 1 ? name : name * '_' * string(i)
            id = reckognize_topology(genome, newarc)
            if id == "UNKNOWN"
                newarc[genome] = newname
                ret[newname] = genome
            else
                ret[newname] = id
            end
        end
    end
    return ret, failed
end
