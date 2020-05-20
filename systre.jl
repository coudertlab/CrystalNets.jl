include("output.jl")
include("input.jl")
include("arithmetics.jl")

import LinearAlgebra: det
import Statistics: mean
using StaticArrays
using LightGraphs
using .CIFTypes.PeriodicGraphs
import .CIFTypes.PeriodicGraphs: hash_position
using Base.Threads

function back_to_unit(r::Rational)
    return Base.unsafe_rational(mod(numerator(r), denominator(r)), denominator(r))
end

function check_valid_translation(c::CrystalNet{T}, t::SVector{3,<:Rational{<:Integer}}, r=nothing) where T
    n = length(c.pos)
    vmap = Int[]
    offsets = SVector{3,Int}[]
    for k in 1:n
        curr_pos = c.pos[k]
        transl = (isnothing(r) ? curr_pos : (r * curr_pos)) .+ t
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
        has_edge(c.graph, PeriodicGraphs.unsafe_edge{3}(src, dst, newofs)) || return nothing
    end
    return vmap
end

function possible_translations(c::CrystalNet{T}) where T
    n = length(c.pos)
    ts = Tuple{Int, Int, Int, SVector{3,T}}[]
    sortedpos = copy(c.pos)
    origin = popfirst!(sortedpos)
    @assert iszero(origin)
    sort!(sortedpos, by=(T === Rational{Int32} ? norm∘SVector{3,Rational{Int64}} : norm))
    for t in sortedpos
        @assert t == back_to_unit.(t)
        max_den, i_max_den = findmax(denominator.(t))
        numerator(t[i_max_den]) == 1 || continue
        nz = count(iszero, t)
        push!(ts, (nz, i_max_den, max_den, t))
    end
    return sort!(ts; by=(x->(x[1], x[2], x[3])))
end

# TODO remove?
function find_first_valid_translations(c::CrystalNet)
    for (nz, i_max_den, max_den, t) in possible_translations(c)
        !isnothing(check_valid_translation(c, t)) && return (nz, i_max_den, max_den, t)
    end
    return nothing
end

function find_all_valid_translations(c::CrystalNet{T}) where T
    ret = NTuple{3, Vector{Tuple{Int, Int, SVector{3,T}}}}(([], [], []))
    for (nz, i_max_den, max_den, t) in possible_translations(c)
        !isnothing(check_valid_translation(c, t)) && push!(ret[nz+1], (i_max_den, max_den, t))
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
    detmax = 1//1
    best = (0, 0, 0)
    @inbounds for i in 1:n-2
        for j in i+1:n-1
            for k in j+1:n
                d = abs(det(hcat(all[i][3], all[j][3], all[k][3])))
                d == 0 && continue
                if d < detmax
                    detmax = d
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

function reduce_with_matrix(c::CrystalNet{T}, mat) where T
    lengths = degree(c.graph)
    cell = Cell(c.cell, c.cell.mat * mat)
    imat = Int.(inv(mat)) # The inverse should only have integer coefficients
    poscol = (imat,) .* c.pos

    offset = [floor.(Int, x) for x in poscol]
    for i in 1:length(poscol)
        poscol[i] = poscol[i] .- offset[i]
    end
    I_sort = sort(1:length(poscol); by=i->(poscol[i], hash_position(offset[i])))
    i = popfirst!(I_sort)
    @assert iszero(offset[i])
    I_kept = Int[i]
    sortedcol = SVector{3,T}[poscol[i]]
    for i in I_sort
        x = poscol[i]
        if x != sortedcol[end]
            push!(I_kept, i)
            push!(sortedcol, x)
        end
    end
    edges = PeriodicEdge3D[]
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
    graph = PeriodicGraph3D(edges)
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
    newbasis = nf3D(offsets)
    invbasis = inv(newbasis)
    intcoords = [Int.(invbasis * x) for x in offsets]
    return newbasis, intcoords
end

function findbasis(edges::Vector{Tuple{Int,Int,SVector{3,T}}}, nonzerosoffsets) where T
    nzoffsets = [edges[i][3] for i in nonzerosoffsets]
    I_sorted = sortperm(nzoffsets)
    sorted_uniques = SVector{3,T}[nzoffsets[I_sorted[1]]]
    I_uniques = Int[I_sorted[1]]
    perm = Vector{Int}(undef, length(nonzerosoffsets))
    perm[I_sorted[1]] = 1
    for i in 2:length(nonzerosoffsets)
        x = nzoffsets[I_sorted[i]]
        if x != sorted_uniques[end]
            push!(sorted_uniques, x)
            push!(I_uniques, I_sorted[i])
        end
        perm[I_sorted[i]] = length(I_uniques)
    end
    original_order = sortperm(I_uniques)
    uniques = sorted_uniques[original_order]

    basis, unique_coords = findfirstbasis(uniques)

    n = length(edges)
    order = zeros(Int, n)
    reverse_order = invperm(original_order)
    for i in 1:length(nonzerosoffsets)
        order[nonzerosoffsets[i]] = reverse_order[perm[i]]
    end

    newedges = Vector{PeriodicEdge3D}(undef, n)
    for j in 1:n
        i = order[j]
        src, dst, _ = edges[j]
        if i == 0
            newedges[j] = PeriodicEdge3D(src, dst, zero(SVector{3,Int}))
        else
            newedges[j] = PeriodicEdge3D(src, dst, unique_coords[i])
        end
    end

    @assert all(z-> begin x, y = z; x == (y.src, y.dst.v, basis*y.dst.ofs) end, zip(edges, newedges))

    return basis, newedges
end


function candidate_key(net::CrystalNet{T}, u, basis) where T
    U = widen(T)
    V = widen(U)
    n = nv(net.graph)
    h = 2 # Next node to assign
    origin = net.pos[u]
    newpos = Vector{SVector{3,V}}(undef, n) # Positions of the kept representatives
    newpos[1] = zero(SVector{3,V})
    offsets = Vector{SVector{3,Int}}(undef, n)
    # Offsets of the new representatives with respect to the original one, in the original basis
    offsets[1] = zero(SVector{3,Int})
    vmap = Vector{Int}(undef, n) # Bijection from the old to the new node number
    vmap[1] = u
    rev_vmap = zeros(Int, n) # inverse of vmap
    rev_vmap[u] = 1
    edgs = Tuple{Int,Int,SVector{3,V}}[]
    mat = inv(big.(basis))
    nonzerosoffsets = Int[]
    for t in 1:n # t is the node being processed
        neighs = neighbors(net.graph, vmap[t])
        ofst = offsets[t]
        pairs = Vector{Tuple{SVector{3,V},Int}}(undef, length(neighs))
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
        bigbasis = big.(basis)
        for (coordinate, v) in pairs
            idx = rev_vmap[v]
            if idx == 0 # New node to which h is assigned
                vmap[h] = v
                rev_vmap[v] = h
                newpos[h] = coordinate
                offsets[h] = SVector{3,Int}(bigbasis * coordinate .+ origin .- net.pos[v])
                push!(edgs, (t, h, zero(SVector{3,T})))
                h += 1
            else
                realofs = coordinate .- newpos[idx]
                # offset between this representative of the node and that which was first encountered
                push!(edgs, (t, idx, realofs))
                iszero(realofs) || push!(nonzerosoffsets, length(edgs))
            end
        end
    end
    @assert allunique(edgs)

    newbasis, coords = findbasis(edgs, nonzerosoffsets)
    newbasis = basis * newbasis
    @assert abs(det(newbasis)) == 1

    # newedges = collect(edges(PeriodicGraph3D(nv(net.graph), coords))) # TODO improve?
    # @show newbasis
    # @show newedges
    # @show vmap
    # @show collect(edges(net.graph[rev_vmap]))
    # global g = CrystalNet(Cell(net.cell, net.cell.mat), net.types[vmap], newpos, PeriodicGraph3D(n, coords))
    # @assert all(g.pos[i] == mean(g.pos[x.v] .+ x.ofs for x in neighbors(g.graph, i)) for i in 1:n)
    # @assert all(z-> begin x, y = z; (x.src, x.dst.v, (basis*newbasis)*x.dst.ofs) == (y.src, y.dst.v, y.dst.ofs) end, zip(edges(net.graph[vmap]), newedges))

    return newbasis, vmap, PeriodicGraph3D(n, coords)
end


mutable struct SpglibDataset
    spacegroup_number::Cint
    hall_number::Cint
    _international_symbol::NTuple{11,Cchar}
    _hall_symbol::NTuple{17,Cchar}
    _choice::NTuple{6,Cchar}
    _transformation_matrix::NTuple{9,Cdouble}
    _origin_shift::NTuple{3,Cdouble}
    n_operations::Cint
    _rotations::Ptr{Cint}
    _translations::Ptr{Cdouble}
    n_atoms::Cint
    _wyckoffs::Ptr{Cint}
    _site_symmetry_symbols::Ptr{NTuple{7,Cchar}}
    _equivalent_atoms::Ptr{Cint}
    _crystallographic_orbits::Ptr{Cint}
    _primitive_lattice::NTuple{9,Cdouble}
    _mapping_to_primitive::Ptr{Cint}
    n_std_atoms::Cint
    _std_lattice::NTuple{9,Cdouble}
    _std_types::Ptr{Cint}
    _std_positions::Ptr{Cdouble}
    _std_rotation_matrix::NTuple{9,Cdouble}
    _std_mapping_to_primitive::Ptr{Cint}
    _pointgroup_symbol::NTuple{6,Cchar}
end
function Base.getproperty(ds::SpglibDataset, name::Symbol)
    if name === :international_symbol
        letters = Char.(getfield(ds, :_international_symbol))
        fzero = findfirst(==('\0'), collect(letters))
        return join(letters[i] for i in 1:(isnothing(fzero) ? length(letters) : (fzero-1)))
    elseif name === :hall_symbol
        letters = Char.(getfield(ds, :_hall_symbol))
        fzero = findfirst(==('\0'), collect(letters))
        return join(letters[i] for i in 1:(isnothing(fzero) ? length(letters) : (fzero-1)))
    elseif name === :choice
        letters = Char.(getfield(ds, :_choice))
        fzero = findfirst(==('\0'), collect(letters))
        return join(letters[i] for i in 1:(isnothing(fzero) ? length(letters) : (fzero-1)))
    elseif name === :transformation_matrix
        return reshape(collect(getfield(ds, :_transformation_matrix)), (3,3))
    elseif name === :origin_shift
        return reshape(collect(getfield(ds, :_origin_shift)), 3)
    elseif name === :rotations
        return unsafe_wrap(Array, getfield(ds, :_rotations), (3,3,Int(getfield(ds, :n_operations))))
    elseif name === :translations
        return unsafe_wrap(Array, getfield(ds, :_translations), (3,Int(getfield(ds, :n_operations))))
    elseif name === :wyckoffs
        return unsafe_wrap(Array, getfield(ds, :_wyckoffs), Int(getfield(ds, :n_atoms)))
    elseif name === :site_symmetry_symbols
        a = unsafe_wrap(Array, getfield(ds, :_site_symmetry_symbols), getfield(ds, :n_atoms))
        ret = Vector{String}(undef, length(a))
        for i in 1:length(a)
            letters = Char.(a[i])
            fzero = findfirst(==('\0'), collect(letters))
            ret[i] = join(letters[j] for j in 1:(isnothing(fzero) ? length(letters) : (fzero-1)))
        end
        return ret
    elseif name === :equivalent_atoms
        return unsafe_wrap(Array, getfield(ds, :_equivalent_atoms), getfield(ds, :n_atoms))
    elseif name === :crystallographic_orbits
        return unsafe_wrap(Array, getfield(ds, :_crystallographic_orbits), getfield(ds, :n_atoms))
    elseif name === :primitive_lattice
        return reshape(collect(getfield(ds, :_primitive_lattice)), (3,3))
    elseif name === :mapping_to_primitive
        return unsafe_wrap(Array, getfield(ds, :_mapping_to_primitive), getfield(ds, :n_atoms))
    elseif name === :std_lattice
        return reshape(collect(getfield(ds, :_std_lattice)), (3,3))
    elseif name === :std_types
        return unsafe_wrap(Array, getfield(ds, :_std_types), getfield(ds, :n_std_atoms))
    elseif name === :std_positions
        return unsafe_wrap(Array, getfield(ds, :_std_positions), (3,Int(getfield(ds, :n_std_atoms))))
    elseif name === :std_rotation_matrix
        return reshape(collect(getfield(ds, :_std_rotation_matrix)), (3,3))
    elseif name === :std_mapping_to_primitive
        return unsafe_wrap(Array, getfield(ds, :_std_mapping_to_primitive), getfield(ds, :n_std_atoms))
    elseif name === :pointgroup_symbol
        letters = Char.(getfield(ds, :_pointgroup_symbol))
        fzero = findfirst(==('\0'), collect(letters))
        return join(letters[i] for i in 1:(isnothing(fzero) ? length(letters) : (fzero-1)))
    else
        getfield(ds, name)
    end
end

function get_spglib_dataset(net::CrystalNet)
    lattice = Matrix{Cdouble}(net.cell.mat') # transpose to account for row-major operations
    n = nv(net.graph)
    positions = Matrix{Cdouble}(undef, 3, n)
    types = Vector{Cint}(undef, n)
    symb_to_int = Dict{Symbol,Cint}()
    j = 1
    for i in 1:n
        positions[:,i] .= net.cell.mat * net.pos[i]
        j += ((types[i] = get!(symb_to_int, net.types[i], j)) == j)
    end
    ptr = ccall((:spg_get_dataset, :libsymspg), Ptr{SpglibDataset},
                (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
                lattice, positions, types, n, 10*eps(Cdouble))
    dataset = unsafe_load(ptr)
    @assert dataset.n_atoms == n # otherwise the net is not minimal
    return dataset
end

function find_symmetries(net::CrystalNet{Rational{S}}) where S
    lattice = Matrix{Cdouble}(LinearAlgebra.I, 3, 3) # positions are expressed in this basis

    I = sortperm(net.pos)
    uniquepos = SVector{3, Rational{S}}[net.pos[I[1]]]
    symb_to_int = Dict{Vector{Tuple{Symbol,Int}}, Cint}([(net.types[I[1]], 1)] => 1)
    last_types = [(net.types[I[1]], 1)]
    types = Cint[1]
    new_symbol = 2
    for i in 2:length(net.pos)
        j = I[i]
        if net.pos[j] == last(uniquepos)
            flag = true
            for k in 1:length(last_types)
                thistype, x = last_types[k]
                if thistype == net.types[j]
                    last_types[k] = (thistype, x+1)
                    flag = false
                    break
                end
            end
            if flag
                push!(last_types, (net.types[j], 1))
            end
            symb = get!(symb_to_int, sort!(last_types), new_symbol)
            types[end] = symb
            new_symbol += (symb == new_symbol)
        else
            push!(uniquepos, net.pos[j])
            last_types = [(net.types[j], 1)]
            symb = get!(symb_to_int, sort!(last_types), new_symbol)
            push!(types, symb)
            new_symbol += (symb == new_symbol)
        end
    end

    n = length(uniquepos)
    positions = Matrix{Cdouble}(undef, 3, n)
    den = 1
    for i in 1:n
        den = lcm(den, lcm((S == Int32 ? Int64 : S).(denominator.(uniquepos[i]))))
        positions[:,i] .= uniquepos[i]
    end

    len = ccall((:spg_get_multiplicity, :libsymspg), Cint,
                (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
                lattice, positions, types, n, 100*eps(Cdouble))
    @assert len >= 1
    rotations = Array{Cint}(undef, 3, 3, len)
    translations = Array{Cdouble}(undef, 3, len)
    _len = ccall((:spg_get_symmetry, :libsymspg), Cint,
                 (Ptr{Cint}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
                 rotations, translations, len, lattice, positions, types, n, 100*eps(Cdouble))
    @assert _len == len
    # @show len
    den = lcm(len, den)
    # The first symmetry should always be the identity, which we can skip
    @assert isone(rotations[:,:,1])
    @assert iszero(translations[:,1])
    symmetries = Tuple{SMatrix{3,3,Int,9}, SVector{3,Rational{S}}}[(rotations[:,:,1], translations[:,1])]
    vmaps = Vector{Int}[collect(1:nv(net.graph))]
    floatpos = [float(x) for x in net.pos]
    hasmirror = false # whether a mirror plane exists or not. If so, the graph is achiral
    for i in 2:len
        rot = SMatrix{3,3,Int,9}(transpose(rotations[:,:,i]))
        if det(rot) < 0
            @assert det(rot) == -1
            hasmirror = true
            # continue # No need to keep it as its mirror image will be kept
        end
        # @assert det(rot) == 1
        tr = SVector{3,Cdouble}(translations[:,i])
        trans = SVector{3,Rational{S}}(round.(S == Int32 ? Int64 : S, den .* tr) .// den)
        vmap = check_valid_translation(net, trans, rot)
        if isnothing(vmap)
            trans = SVector{3,Rational{S}}(net.pos[last(findmin([norm(x .- tr) for x in floatpos]))])
            vmap = check_valid_translation(net, trans, rot)
            # if isnothing(vmap)
            #     @show translations[:,i]
            #     @show float(trans)
            #     @show trans
            #     @show rot
            #     println("net = CrystalNet(PeriodicGraph3D(\"", net.graph, "\"))")
            #     throw("Invalid symmetry") # This error is for debugging only, false positives without it being an error
            # end
            isnothing(vmap) && continue
        end
        @assert rot != LinearAlgebra.I # Otherwise the symmetry is a pure translation, so the net is not minimal
        push!(symmetries, (rot, trans))
        push!(vmaps, vmap::Vector{Int})
    end
    # @show vmaps
    # @show length(symmetries)
    return (symmetries, vmaps, hasmirror)
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
function partition_by_coordination_sequence(graph, vmaps=Vector{Int}[])
    # First, union-find on the symmetries to avoid computing useless coordination sequences
    n = nv(graph)
    unionfind = collect(1:n)
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

    vsequences = [vertex_sequence(graph, first(cat), 10) for cat in categories]
    I = sortperm(vsequences)
    @assert vsequences[I[1]][1] >= 2 # vertices of degree <= 1 should have been removed at input creation

    todelete = falses(length(categories))
    last_i = I[1]
    for j in 2:length(I)
        i = I[j]
        if vsequences[i] == vsequences[last_i]
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
    deleteat!(vsequences, todelete)

    num = length(categories)
    numoutgoingedges = Vector{Tuple{Int,Vector{Int}}}(undef, num)
    for i in 1:num
        seq = vsequences[i]
        numoutgoingedges[i] = (length(categories[i]) * seq[1], seq)
    end
    sortorder = sortperm(numoutgoingedges)
    # categories are sorted by the total number of outgoing directed edges from each category
    # and by the vertex sequence in case of ex-aequo.
    # categories are thus uniquely determined and ordered independently of the representation of the net

    @assert allunique(vsequences)
    for i in 1:num
        @assert all(vertex_sequence(graph, x, 10) == vsequences[i] for x in categories[i])
    end # TODO comment out these costly asserts
    return categories[sortorder], cat_map, unique_reprs[sortorder]
end


# function find_candidates(net::CrystalNet{T}) where T
#     candidates = Pair{Int,SMatrix{3,3,T,9}}[]
#     n = nv(net.graph)
#     cols = Vector{SVector{3,T}}(undef, 3)
#     for u in 1:n
#         posu = net.pos[u]
#         neighs = neighbors(net.graph, u)
#         for u1 in neighs
#             cols[1] = net.pos[u1.v] .+ u1.ofs .- posu
#             for u2 in neighs
#                 cols[2] = net.pos[u2.v] .+ u2.ofs .- posu
#                 for u3 in neighs
#                     cols[3] = net.pos[u3.v] .+ u3.ofs .- posu
#                     mat = hcat(cols[1], cols[2], cols[3])
#                     d = det(mat)
#                     iszero(d) && continue
#                     i = 1; j = 2; k = 3
#                     if d < 0
#                         j = 3; k = 2
#                     end
#                     for _ in 1:3
#                         mat = hcat(cols[i], cols[j], cols[k])
#                         push!(candidates, u => mat)
#                         i, j, k = j, k, i
#                     end
#                 end
#             end
#         end
#     end
#     return candidates
# end

"""
Find a set of pairs (v, basis) that are candidates for finding the key with
candidate_key(net, v, basis).
The returned set is independent of the representation of the graph used in net.
"""
function find_candidates(net::CrystalNet{T}) where T
    symmetries, vmaps, hasmirror = find_symmetries(net)
    @assert length(symmetries) == length(vmaps)
    categories, symmetry_map, unique_reprs = partition_by_coordination_sequence(net.graph, vmaps)
    # @show unique_reprs
    # @show categories
    # @show first(partition_by_coordination_sequence(net.graph, []))
    @assert sort.(categories) == sort.(first(partition_by_coordination_sequence(net.graph, [])))
    candidates = Dict{Int,Vector{SMatrix{3,3,T,9}}}()
    for reprs in unique_reprs
        # First, we try to look for triplet of edges all starting from the same vertex within a category
        degree(net.graph, first(reprs)) <= 3 && continue
        candidates = find_candidates_onlyneighbors(net, reprs)
        !isempty(candidates) && break
    end
    if isempty(candidates)
        # If we arrive at this point, it means that all vertices only have coplanar neighbors.
        for (i, cat) in enumerate(categories)
            # Then we look for triplet of edges two of which start from a vertex and one does not, but all
            # three belonging to the same category.
            length(cat) == 1 && continue
            for v in unique_reprs[i]
                others = [x for x in cat if x != v]
                mergewith!(append!, candidates, find_candidates_fallback(net, [v], [others]))
            end
            !isempty(candidates) && break
        end
        if isempty(candidates)
            for (i, reprs) in enumerate(unique_reprs)
                # Finally, we look for the triplet of edges two of which start from a vertex and one does not,
                # without the constraint that all three belong to the same category.
                othercats = [categories[j+(j>=i)] for j in 1:length(categories)-1]
                candidates = find_candidates_fallback(net, reprs, othercats)
                !isempty(candidates) && break
            end
        end
    end
    @assert !isempty(candidates)
    # @show candidates
    # isempty(vmaps) && return candidates
    # return candidates
    return extract_through_symmetry(candidates, vmaps, symmetries)
end

function extract_through_symmetry(candidates::Dict{Int,Vector{SMatrix{3,3,T,9}}}, vmaps, symmetries) where T
    # @show sum(length, values(candidates))
    unique_candidates = Pair{Int,SMatrix{3,3,T,9}}[]
    for (i, mats) in candidates
        @assert i == minimum(vmap[i] for vmap in vmaps)
        indices = [j for j in 1:length(vmaps) if vmaps[j][i] == i]
        min_mats = Set{SVector{9,T}}()
        for mat in mats
            # # push!(unique_candidates, i => mat)
            # for j in indices
            #    push!(unique_candidates, i => first(symmetries[j]) * mat)
            # end
            flattenedmat = SVector{9,T}(mat)
            min_mat = flattenedmat
            for j in indices
                new_mat = SVector{9,T}(first(symmetries[j]) * mat)
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
    # @show unique_candidates
    # @show length(collect(unique_candidates))
    return unique_candidates
end

function find_candidates_onlyneighbors(net::CrystalNet{T}, candidates_v) where T
    # throw("UNTESTED")
    deg = degree(net.graph, first(candidates_v)) # The degree is the same for all vertices of the same category
    initial_candidates = Pair{Int,Matrix{T}}[]
    fastlock = SpinLock()
    @threads for v in candidates_v
        a = Matrix{T}(undef, 3, deg)
        posi = net.pos[v]
        for (j, x) in enumerate(neighbors(net.graph, v))
            a[:,j] .= net.pos[x.v] .+ x.ofs .- posi
        end
        if LinearAlgebra.rank(a) >= 3
            lock(fastlock) do
                push!(initial_candidates, v => a)
            end
        end
        # @show hcat(sort(collect(eachcol(a)))...)
    end

    candidates = Dict{Int,Vector{SMatrix{3,3,T,9}}}()
    isempty(initial_candidates) && return candidates
    # @show initial_candidates
    @threads for (v, mat) in initial_candidates
        _, n = size(mat)
        matv = SMatrix{3,3,T,9}[]
        for _i in 1:n, _j in (_i+1):n, _k in (_j+1):n
            m = SMatrix{3,3,T,9}(mat[:,[_i,_j,_k]])
            d = det(big.(m))
            iszero(d) && continue
            i = 1; j = 2; k = 3
            mats = Vector{SMatrix{3,3,T,9}}(undef, 6)
            for l in 1:3
                mats[l] = m[:,[i,j,k]]
                i, j, k = j, k, i
            end
            i, j, k = i, k, j
            for l in 4:6
                mats[l] = m[:,[i,j,k]]
                i, j, k = j, k, i
            end
            append!(matv, mats)
        end
        lock(fastlock) do
            candidates[v] = matv
        end
    end
    return candidates
end

# function find_candidates_onlyneighbors(net::CrystalNet{T}, category, symmetry_map, num_symmetries, hasmirror) where T
#     deg = degree(net.graph, first(category)) # The degree is the same for all vertices of the same category
#     candidates_v = copy(category)
#     neighbormatrices = [Matrix{T}[] for _ in 1:nthreads()]
#     todelete = falses(length(candidates_v))
#     symmetry_done = falses(num_symmetries) # categories that were already explored
#     @threads for (i, v) in collect(enumerate(candidates_v))
#         if symmetry_done[symmetry_map[v]]
#             todelete[i] = true
#             continue
#         end
#         symmetry_done[symmetry_map[v]] = true
#         a = Matrix{T}(undef, 3, deg)
#         posi = net.pos[v]
#         for (j, x) in enumerate(neighbors(net.graph, v))
#             a[:,j] .= net.pos[x.v] .+ x.ofs .- posi
#         end
#         if LinearAlgebra.rank(a) >= 3
#             push!(neighbormatrices[threadid()], a)
#         else
#             todelete[i] = true
#         end
#         # @show hcat(sort(collect(eachcol(a)))...)
#     end
#     neighbormatrix = Iterators.flatten(neighbormatrices)
#     candidates = Pair{Int,SMatrix{3,3,T,9}}[]
#     isempty(neighbormatrix) && return candidates
#     deleteat!(candidates_v, todelete)
#     @assert length(candidates_v) == length(collect(neighbormatrix))
#     for (v, mat) in zip(candidates_v, neighbormatrix)
#         _, n = size(mat)
#         for _i in 1:n, _j in (_i+1):n, _k in (_j+1):n
#             m = SMatrix{3,3,T}(mat[:,[_i,_j,_k]])
#             d = det(m)
#             iszero(d) && continue
#             i = 1; j = 2; k = 3
#             # flag = false
#             # for (rot, _) in symmetries
#             #     cmpmat = rot * m
#             #     for (_, x) in candidates
#             #         if x == m
#             #             flag = true
#             #             break
#             #         end
#             #     end
#             #     flag && break
#             # end
#             # flag && continue
#             if hasmirror
#                 for _ in 1:3
#                     push!(candidates, v => m[:,[i,j,k]])
#                     i, j, k = j, k, i
#                 end
#                 i, j, k = i, k, j
#             else
#                 j += (d < 0)
#                 k -= (d < 0)
#             end
#             for _ in 1:3
#                 push!(candidates, v => m[:,[i,j,k]])
#                 i, j, k = j, k, i
#             end
#         end
#     end
#     return candidates
# end

function find_candidates_fallback(net::CrystalNet{T}, reprs, othercats) where T
    # throw("UNTESTED")
    candidates = Dict{Int,Vector{SMatrix{3,3,T,9}}}(u => [] for u in reprs)
    for cat in othercats
        for u in reprs
            mats = SMatrix{3,3,T,9}[]
            posu = net.pos[u]
            neighu = neighbors(net.graph, u)
            for (j1, x1) in enumerate(neighu), j2 in j1+1:length(neighu)
                x2 = neighu[j2]
                vec1 = net.pos[x1.v] .+ x1.ofs .- posu
                vec2 = net.pos[x2.v] .+ x2.ofs .- posu
                j = iszero(vec1[1]) ? iszero(vec1[2]) ? 3 : 2 : 1
                vec1[j] * vec2 == vec2[j] .* vec1 && continue # colinear
                for v in cat
                    posv = net.pos[v]
                    for x3 in neighbors(net.graph, v)
                        vec3 = net.pos[x3.v] .+ x3.ofs .- posv
                        mat = hcat(vec1, vec2, vec3)
                        d = det(Rational{Int128}.(mat))
                        iszero(d) && continue
                        # if d < 0
                            # mat = hcat(vec2, vec1, vec3)
                        # end
                        push!(mats, mat)
                        push!(mats, hcat(vec2, vec1, vec3))
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
    return Dict{Int,Vector{SMatrix{3,3,T,9}}}()
end

function parallel_systre_key(net::CrystalNet{T}) where T
    candidates = find_candidates(net)
    # @show length(candidates)
    numthreads = min(nthreads(), length(candidates))
    minimal_graph = Vector{PeriodicGraph3D}(undef, numthreads)
    minimal_vmap = Vector{Vector{Int}}(undef, numthreads)
    minimal_basis = Vector{SMatrix{3,3,T,9}}(undef, numthreads)
    @threads for i in 1:numthreads
        v, basis = candidates[end+1-i]
        minimal_basis[i], minimal_vmap[i], minimal_graph[i] = candidate_key(net, v, basis)
    end
    resize!(candidates, length(candidates) - numthreads)
    @threads for (v, basis) in candidates
        id = threadid() # necessarily inbounds because otherwise the candidates is empty
        newbasis, vmap, graph = candidate_key(net, v, basis)
        if edges(graph) < edges(minimal_graph[id])
            minimal_graph[id] = graph
            minimal_vmap[id] = vmap
            minimal_basis[id] = newbasis
        end
    end
    _, j = findmin(edges.(minimal_graph))
    return minimal_basis[j], minimal_vmap[j], minimal_graph[j]
end

function systre_key(net::CrystalNet{T}) where T
    candidates = find_candidates(net)
    # @show length(candidates)
    v, basis = popfirst!(candidates)
    minimal_basis, minimal_vmap, minimal_graph = candidate_key(net, v, basis)
    for (v, basis) in candidates
        newbasis, vmap, graph = candidate_key(net, v, basis)
        if edges(graph) < edges(minimal_graph)
            minimal_graph = graph
            minimal_vmap = vmap
            minimal_basis = newbasis
        end
    end
    return minimal_basis, minimal_vmap, minimal_graph
end

function find_new_representation(pos, basis, vmap, graph)
    return collect(eachcol(CIFTypes.equilibrium(graph))) # TODO optimize
end


function systre(g::PeriodicGraph3D, skipminimize=true)
    net = CrystalNet(g)
    if !skipminimize
        net = minimize(net)
    end
    return last(systre_key(net))
end

function systre(path::AbstractString)
    crystal = Crystal(parse_cif(path))
    return systre(CrystalNet(crystal), false)
    # return Crystal(crystal.cifinfo, net.cell, collect(1:length(net.types)),
    #                net.types, reduce(hcat, net.pos), net.graph)
end

function systre(net::CrystalNet, skipminimize)
    if !skipminimize
        net = minimize(net)
    end
    return systre(net)
end

function systre(net::CrystalNet)
    basis, vmap, graph = systre_key(net)

    return CrystalNet(graph) # TODO remove this temprorary bandaid

    positions = find_new_representation(net.pos, basis, vmap, graph)

    ret = CrystalNet(Cell(net.cell, net.cell.mat*basis), net.types[vmap], positions, graph)
    means = mean(ret.pos)
    for i in 1:length(ret.pos)
        ret.pos[i] = ret.pos[i] .+ [1//2, 1//2, 1//2] .- means
    end
    @assert all(ret.pos[i] == mean(ret.pos[x.v] .+ x.ofs for x in neighbors(ret.graph, i)) for i in 1:length(net.pos))

    return ret
end
