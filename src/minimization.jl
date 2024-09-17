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
    find_all_valid_translations(c::Union{Crystal,CrystalNet{D}}, collisions) where D

Return a `D`-tuple of list of tuples `(i_max_den, max_den, t)` (see
[`possible_translations`](@ref) for interpretation) where the `n`-th list contains all
valid translations of the net having exactly `n-1` zeros.

A translation is valid if it maps exactly each vertex to a vertex and each edge to an edge.

See also: [`possible_translations`](@ref), `PeriodicGraphEmbeddings.check_valid_symmetry`
"""
function find_all_valid_translations(shrunk_net::CrystalNet{D,T}, collisions) where {D,T}
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
            return c, nothing
        end
    end

    graph = PeriodicGraph{D}(edges)
    opts = permute_mapping!(c.options, vmap)
    return CrystalNet{D,Rational{T}}(cell, c.types[I_kept], sortedcol, graph, opts), newcollisions
end


"""
    minimize(net::CrystalNet, [collisions::Vector{CollisionNode}])

Return a CrystalNet representing the same net as the input, but in a unit cell.
If `collisions` is given, also return the corresponding collisions after minimization.

The computed unit cell may depend on the representation of the input, i.e. it is not
topologicallly invariant.
"""
function minimize(net::CrystalNet, collisions::CollisionList)
    translations = find_all_valid_translations(net, collisions)
    all(isempty.(translations)) && return net, collisions
    mat = minimal_volume_matrix(translations)
    _net, collisions = reduce_with_matrix(net, mat, collisions)
    collisions isa CollisionList && @toggleassert all(isempty.(find_all_valid_translations(_net, collisions)))
    return _net, collisions
end
minimize(net::CrystalNet) = minimize(net, CollisionList(net.pge.g, UnitRange{Int}[]))[1]



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

