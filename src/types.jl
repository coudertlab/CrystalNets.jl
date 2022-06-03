## Type definitions for intermediate representations of crystals up to a net

import Base: ==

function cell_with_warning(mat::StaticArray{Tuple{3,3},BigFloat})
    if !all(isfinite, mat) || iszero(det(mat))
        @ifwarn @error lazy"Suspicious unit cell of matrix $(Float64.(mat)). Is the input really periodic? Using a cubic unit cell instead."
        return Cell()
    end
    return Cell(Cell(), mat)
end

## CIF

"""
    CIF

Representation of a .cif file.
"""
struct CIF
    cifinfo::Dict{String, Union{String, Vector{String}}}
    cell::Cell{Rational{Int}}
    ids::Vector{Int}
    types::Vector{Symbol}
    pos::Matrix{Float64}
    bonds::Vector{Vector{Tuple{Int,Float32}}}
end

function keepinbonds(bonds::Vector{Vector{Tuple{Int,Float32}}}, keep::Vector{Int})
    @toggleassert issorted(keep)
    (isempty(keep) || isempty(bonds)) && return Vector{Tuple{Int,Float32}}[]
    n = length(bonds)
    idxs =  Vector{Int}(undef, n)
    idx = 0
    keepidx = keep[1]
    for i in 1:n
        if i == keepidx
            idx += 1
            keepidx = idx == length(keep) ? 0 : keep[idx + 1]
        end
        idxs[i] = idx
    end
    ret = [Tuple{Int,Float32}[] for _ in 1:length(keep)]
    for _i in 1:length(keep)
        vec = ret[_i]
        for (j, d) in bonds[keep[_i]]
            idx = idxs[j]
            if idx != 0 && keep[idx] == j
                push!(vec, (idx, d))
            end
        end
    end
    return ret
end

function add_to_bondlist!(bondlist::Vector{Tuple{Int,Float32}}, x::Int, d::Float32)
    addedflag = false
    for (i, (j, _)) in enumerate(bondlist)
        if j ≥ x
            if j == x
                bondlist[i] = (x, d)
            else
                splice!(bondlist, i:i-1, [(x, d)])
            end
            addedflag = true
            break
        end
    end
    addedflag || push!(bondlist, (x, d))
    nothing
end

function get_bondlist(bondlist::Vector{Tuple{Int,Float32}}, x::Int)
    for (j, d) in bondlist
        if j ≥ x
            if j == x
                return d
            end
            return Inf32
        end
    end
    return Inf32
end

function sortprune_bondlist!(bondlist::Vector{Tuple{Int,Float32}})
    sort!(bondlist)
    toremove = Int[]
    k = 0
    for (i, (j, _)) in enumerate(bondlist)
        if j == k
            push!(toremove, i)
        else
            k = j
        end
    end
    isempty(toremove) || deleteat!(bondlist, toremove)
    nothing
end


"""
    remove_partial_occupancy(::CIF)

Only keep one atom per atom site.
"""
function remove_partial_occupancy(cif::CIF)
    points::Vector{SVector{3,Float64}} = collect(eachcol(cif.pos))
    perm = sortperm(points)
    n = length(perm)
    # minimal_length = cbrt(norm(cif.cell.mat))*4e-4 # TODO: find an appropriate value
    last_triplet = points[perm[1]] .+ (1,1,1)
    last_position = 0
    toalias = Vector{Int}[]
    inalias = false
    for i in 1:n
        j = perm[i]
        thispoint = points[j]
        if norm(thispoint .- last_triplet) < 4e-4 # TODO: replace with minimal_length
            if inalias
                push!(toalias[end], j)
            else
                push!(toalias, [last_position, j])
                inalias = true
            end
        else
            inalias = false
            last_position = j
            last_triplet = thispoint
        end
    end
    if isempty(toalias) # Nothing to do, we simply alias the CIF to help the compiler (it may be useless)
        return CIF(cif.cifinfo, cif.cell, cif.ids, cif.types, cif.pos, cif.bonds)
    end
    @ifwarn @warn "This CIF file contains at least one site with multiple atoms. Only one atom will be kept per atom site."
    occupancies = parsestrip.(Float64, get(cif.cifinfo, "atom_site_occupancy",
                                      ["1.0" for _ in 1:length(cif.types)])::Vector{String})
    toremove = Int[]
    bonds = [copy(b) for b in cif.bonds]
    for alias in toalias
        m = length(alias)
        occup = @view occupancies[alias]
        max_occupancy = maximum(occup)
        with_max_occupancy = [alias[i] for i in 1:m if occup[i] == max_occupancy]
        representative = minimum(with_max_occupancy)
        append!(toremove, a for a in alias if a != representative)

        isempty(bonds) && continue
        allbonds = bonds[representative]
        for a in alias
            append!(allbonds, bonds[a])
        end
        sort!(allbonds)

        k = 0
        newbonds = Tuple{Int,Float32}[]
        bonds[representative] = newbonds
        for bond in allbonds
            bond[1] == k && continue
            k = bond[1]
            k == representative && continue
            push!(newbonds, bond)
            bondk = bonds[k]
            add_to_bondlist!(bondk, representative, bond[2])
        end
    end

    sort!(toremove)
    ids = deleteat!(copy(cif.ids), toremove)
    tokeep = deleteat!(collect(1:n), toremove)
    bonds = keepinbonds(bonds, tokeep)
    pos = cif.pos[:,tokeep]
    @toggleassert allunique(collect(eachcol(pos)))
    return CIF(cif.cifinfo, cif.cell, ids, cif.types, pos, bonds)
end

"""
    prune_collisions(::CIF)

For each site where there are atoms suspiciously close to one another, remove all
but one of them. This arises for example when all the possible positions of at atom are
superposed in the CIF file, typically for a solvent which should be disregarded anyway.
"""
function prune_collisions(cif::CIF)
    toremove = Int[]
    points::Vector{SVector{3,Float64}} = collect(eachcol(cif.pos))
    n = length(points)
    smallmat = Float64.(cif.cell.mat)
    buffer, ortho, safemin = prepare_periodic_distance_computations(smallmat)
    for i in 1:n, j in (i+1):n
        buffer .= points[i] .- points[j]
        if periodic_distance!(buffer, smallmat, ortho, safemin) < 0.55
            push!(toremove, j)
        end
    end
    if isempty(toremove) # Nothing to do, we simply alias the CIF to help the compiler (it may be useless)
        return false, CIF(cif.cifinfo, cif.cell, cif.ids, cif.types, cif.pos, cif.bonds)
    end
    @ifwarn @warn "This CIF file contains multiple colliding atoms. Only one atom will be kept per site."
    unique!(sort!(toremove))
    ids = deleteat!(copy(cif.ids), toremove)
    tokeep = deleteat!(collect(1:n), toremove)
    bonds = keepinbonds(cif.bonds, tokeep)
    pos = cif.pos[:,tokeep]
    return true, CIF(cif.cifinfo, cif.cell, ids, cif.types, pos, bonds)
end

"""
    expand_symmetry(::CIF)

Applies all the symmetry operations listed in the CIF file to the atoms and the bonds.
"""
function expand_symmetry(c::CIF)
    cif::CIF = remove_partial_occupancy(c)
    if isempty(cif.cell.equivalents)
        return CIF(cif.cifinfo, deepcopy(cif.cell), cif.ids, cif.types, cif.pos, cif.bonds)
    end

    if isempty(cif.bonds)
        oldbonds = Vector{Tuple{Int,Float32}}[]
        knownbondlengths = true
        bonds = Vector{Tuple{Int,Float32}}[]
    else
        maxid = maximum(cif.ids)
        @toggleassert minimum(cif.ids) ≥ 1
        rev_id = [Int[] for _ in 1:maxid]
        for (i, idi) in enumerate(cif.ids)
            push!(rev_id[idi], i)
        end
        oldbonds = [Tuple{Int,Float32}[] for _ in 1:maxid]
        for (idi, revidi) in enumerate(rev_id)
            bondsidi = oldbonds[idi]
            for i in revidi
                for (j, d) in cif.bonds[i]
                    push!(bondsidi, (cif.ids[j], d))
                end
            end
            sortprune_bondlist!(bondsidi)
        end
        knownbondlengths = !any(any(x -> iszero(x[2]), b) for b in oldbonds)
        bonds = [Tuple{Int,Float32}[] for _ in 1:length(cif.ids)]

        if !knownbondlengths
            # This means that the bonds are determined as symmetric images of bonds of the
            # asymetric unit, which may or may not be the convention used in the CIF files.
            # Most CIF files having geom_bond_atom_site_labels should have geom_bond_distance anyway.
            @ifwarn @warn "Expanding CIF symmetry without knowing the bond lengths: the resulting bonds might be erroneous"
        end
    end

    n = length(cif.ids)
    newids::Vector{Int} = copy(cif.ids)
    newpos::Vector{SVector{3,Float64}} = collect(eachcol(cif.pos))
    smallmat = Float64.(cif.cell.mat)
    buffer, ortho, safemin = prepare_periodic_distance_computations(smallmat)
    #=@inbounds=# for equiv in cif.cell.equivalents
        image = zeros(Int, n)
        for i in 1:length(cif.ids)
            v = newpos[i]
            p = Vector(equiv.mat*v .+ equiv.ofs)
            p .-= floor.(p)
            for j in 1:length(newpos)
                buffer .= newpos[j] .- p
                if periodic_distance!(buffer, smallmat, ortho, safemin) < 0.55
                    image[i] = j
                    break
                end
            end
            if image[i] == 0
                push!(newpos, p)
                push!(newids, cif.ids[i])
                isempty(bonds) || push!(bonds, Tuple{Int,Float32}[])
                image[i] = length(newpos)
            end
        end
        if !knownbondlengths
            for (i, bondi) in enumerate(cif.bonds), (j, _) in bondi
                imi = image[newids[i]]
                imj = image[newids[j]]
                add_to_bondlist!(bonds[imi], imj, 0f0)
                add_to_bondlist!(bonds[imj], imi, 0f0)
            end
        end
    end

    if knownbondlengths && !isempty(bonds)
        m = length(newids)
        for i in 1:m
            for j in (i+1):m
                bondlength = get_bondlist(oldbonds[newids[i]], newids[j])
                bondlength < Inf32 || continue
                buffer .= newpos[i] .- newpos[j]
                if abs(periodic_distance!(buffer, smallmat, ortho, safemin) - bondlength) < 0.55
                    push!(bonds[i], (j, bondlength))
                    push!(bonds[j], (i, bondlength))
                end
            end
        end
    end

    @toggleassert all(issorted, cif.bonds)
    @toggleassert all(allunique, cif.bonds)

    return CIF(cif.cifinfo, deepcopy(cif.cell), newids, copy(cif.types), reduce(hcat, newpos), bonds)
end


"""
    edges_from_bonds(bonds::Vector{Vector{Tuple{Int,Float32}}}, mat, pos)

Given a bond list `bonds` containing triplets `(a, b, dist)` where atoms `a` and `b` are
bonded if their distance is lower than `dist`, the 3×3 matrix of the cell `mat` and the
Vector{SVector{3,Float64}} `pos` whose elements are the fractional positions of the
atoms, extract the list of PeriodicEdge3D corresponding to the bonds.
Since the adjacency matrix wraps bonds across the boundaries of the cell, the edges
are extracted so that the closest representatives are chosen to form bonds.
"""
function edges_from_bonds(bonds::Vector{Vector{Tuple{Int,Float32}}},
                          mat::SMatrix{3,3,Float64,9}, pos::Vector{SVector{3,Float64}})
    n = length(pos)
    edges = PeriodicEdge3D[]
    ref_dst = norm(mat*[1, 1, 1])
    for i in 1:n, (k, maxdist) in bonds[i]
        k < i && continue
        @toggleassert k != i
        iszero(maxdist) && continue
        offset::Vector{SVector{3, Int}} = []
        old_dst = ref_dst
        for ofsx in -1:1, ofsy in -1:1, ofsz in -1:1 # TODO: optimize with the periodic_distance trick?
            dst = norm(mat * (pos[i] .- (pos[k] .+ (ofsx, ofsy, ofsz))))
            if maxdist == -Inf32 # only keep the closest
                if dst < maxdist || abs2(dst - old_dst) < 1e-3
                    push!(offset, (ofsx, ofsy, ofsz))
                    old_dst = (dst + old_dst)/2
                elseif dst < old_dst
                    empty!(offset)
                    push!(offset, (ofsx, ofsy, ofsz))
                    old_dst = dst
                end
            elseif dst < maxdist
                push!(offset, (ofsx, ofsy, ofsz))
            end
        end
        for ofs in offset
            push!(edges, (i, k, ofs))
        end
    end
    return edges
end

## Clusters

"""
    Clusters

Classification of the atoms of a crystalline framework in different clusters.
For simple crystals, every atom is its own cluster.
For a MOF, a cluster is a SBU, which can be either organic or inorganic.
"""
struct Clusters
    sbus::Vector{Vector{PeriodicVertex3D}}
    classes::Vector{Int}
    attributions::Vector{Int}
    offsets::Vector{SVector{3,Int}}
    periodic::BitVector
end

function Clusters(n)
    sbus = [[PeriodicVertex3D(i)] for i in 1:n]
    classes = collect(1:n)
    attributions = collect(1:n)
    offsets = [zero(SVector{3,Int}) for _ in 1:n]
    return Clusters(sbus, classes, attributions, offsets, falses(n))
end

Base.isempty(c::Clusters) = c.attributions == 1:length(c.attributions)

function Base.getindex(c::Clusters, vmap::AbstractVector{<:Integer})
    rev_vmap = zeros(Int, length(c.attributions))
    for (i, j) in enumerate(vmap)
        rev_vmap[j] = i
    end
    sbus = Vector{PeriodicVertex3D}[]
    sbu_vmap = Int[]
    sbu_countermap = Vector{Int}(undef, length(c.sbus))
    for (i, sbu) in enumerate(c.sbus)
        newsbu = PeriodicVertex3D[]
        for x in sbu
            y = rev_vmap[x.v]
            y == 0 && continue
            push!(newsbu, PeriodicVertex3D(y, x.ofs))
        end
        if !isempty(newsbu)
            push!(sbu_vmap, i)
            push!(sbus, newsbu)
        end
        sbu_countermap[i] = length(sbus)
    end
    offsets = [c.offsets[i] for i in vmap]
    attributions = sbu_countermap[c.attributions][vmap]
    classes = c.classes[sbu_vmap]
    periodic = c.periodic[sbu_vmap]
    return Clusters(sbus, classes, attributions, offsets, periodic)
end

## Crystal

"""
    Crystal

Intermediate representation of a crystal, retaining information on the cell, and the
fractional placement of the atoms and their type, as well as the residues which will be used as
vertices for the computation of the underlying topology.
"""
struct Crystal{T<:Union{Nothing,Clusters}}
    pge::PeriodicGraphEmbedding3D{Float64}
    types::Vector{Symbol}
    clusters::T
    options::Options

    function Crystal{Clusters}(pge, types, clusters, options)
        return new{Clusters}(pge, types, clusters, options)
    end

    function Crystal{Nothing}(pge, types, options)
        return new{Nothing}(pge, types, nothing, options)
    end
end

function Crystal{Clusters}(cell, types, clusters, pos, graph, options)
    Crystal{Clusters}(PeriodicGraphEmbedding3D(graph, pos, cell), types, clusters, options)
end
function Crystal{Nothing}(cell, types, pos, graph, options)
    Crystal{Nothing}(PeriodicGraphEmbedding3D(graph, pos, cell), types, options)
end

function Crystal(pge, types, clusters, options)
    if clusters isa Nothing
        return Crystal{Nothing}(pge, types, options)
    end
    return Crystal{Clusters}(pge, types, clusters, options)
end
function Crystal(c::Crystal; kwargs...)
    return Crystal(c.pge, c.types, c.clusters, Options(c.options; kwargs...))
end

function ==(c1::Crystal{T}, c2::Crystal{T}) where T
    c1.pge == c2.pge && c1.types == c2.types && c1.clusters == c2.clusters && c1.options == c2.options
end

function Crystal{Nothing}(c::Crystal{T}; kwargs...) where T
    if isempty(kwargs)
        if T === Nothing
            return c
        end
        return Crystal{Nothing}(c.pge, c.types, c.options)
    end
    return Crystal{Nothing}(c.pge, c.types, Options(c.options; kwargs...))
end
function Crystal{Clusters}(c::Crystal, clusters::Clusters; kwargs...)
    if isempty(kwargs)
        return Crystal{Clusters}(c.pge, c.types, clusters, c.options)
    end
    return Crystal{Clusters}(c.pge, c.types, clusters, Options(c.options; kwargs...))
end
function Crystal{Clusters}(c::Crystal{T}; kwargs...) where T
    if T === Clusters
        return Crystal{Clusters}(c, c.clusters; kwargs...)
    end
    Crystal{Clusters}(c, Clusters(length(c.types)); kwargs...)
end

"""
    trimmed_crystal(c::Crystal{Nothing})

Rebuild the crystal after trimming its graph according to [`trim_topology`](@ref CrystalNets.trim_topology).
"""
function trimmed_crystal(c::Crystal{Nothing})
    g = deepcopy(c.pge.g)
    remove_metal_cluster_bonds!(g, c.types, c.options)
    vmap, graph = trim_topology(g)
    types = c.types[vmap]
    pge = PeriodicGraphEmbedding(graph, c.pge.pos[vmap], c.pge.cell)
    opts = isempty(c.options._pos) ? c.options : Options(c.options; _pos=pge.pos)
    return Crystal{Nothing}(pge, types, opts)
end


function Base.getindex(c::Crystal{T}, vmap::AbstractVector{<:Integer}) where T
    types = c.types[vmap]
    if T === Nothing
        return Crystal{Nothing}(c.pge[vmap], types, c.options)
    else
        return Crystal{Clusters}(c.pge[vmap], types, c.clusters[vmap], c.options)
    end
end

## CrystalNet

# For the remainder of the file, we can work in 1D, 2D or 3D

"""
    equilibrium(g::PeriodicGraph)

Return an equilibrium placement for the vertices of the graph, defined as a set
of positions such that each vertex is at the barycentre of its neighbors.
The returned equilibrium placement is such that the first vertex of the graph
is at the origin of the space.
"""
function equilibrium(g::PeriodicGraph{N}) where N
    n = nv(g)
    iszero(n) && return Matrix{Rational{Int64}}(undef, N, 0)
    Y = Matrix{Int}(undef, n, N)
    A = spzeros(Int, n, n)
    neigh = Vector{Int}(undef, n)
    offset = SizedVector{N,Int}(undef)
    for i in 1:n
        neigh .= 0
        offset .= 0
        count = 0
        for k in neighbors(g, i)
            k.v == i && continue
            count += 1
            neigh[k.v] += 1
            offset .-= k.ofs
        end
        Y[i,:] .= offset
        A[i,:] .= neigh
        A[i,i] = -count
    end

    return dixon_solve(Val(N), A[2:end,2:end], Y[2:end,:])
end

"""
    trim_topology(graph::PeriodicGraph)

Return a pair `(vmap, newgraph)` extracted from the input by removing vertices
of valence lower or equal to 1, and by replacing vertices of valence 2 by edges,
until convergence.
The only exceptions are vertices only bonded to their representatives of another
cell: those will not be replaced by edges even if their valence is 2, since this
latter case indicates an irreducible trivial 1-dimensional topology.

`vmap` maps the vertices of `newgraph` to their counterpart in `graph`.
"""
function trim_topology(graph::PeriodicGraph{N}) where N
    newgraph = PeriodicGraph{N}(nv(graph), collect(edges(graph)))
    remove_idx = Int[]
    flag = any(<=(2), degree(newgraph))
    vmap = collect(1:nv(graph))
    ignorecounter = 0
    while flag # we alternate cycles when we remove valence 1 and valence 2
        # until no such vertex remains
        flag = any(<=(1), degree(newgraph))
        while flag
            n = nv(newgraph)
            for i in 1:n
                if degree(newgraph, i) <= 1
                    push!(remove_idx, i)
                end
            end
            map = rem_vertices!(newgraph, remove_idx)
            vmap = vmap[map]
            #pos = pos[vmap]
            empty!(remove_idx)
            flag = any(isone, degree(newgraph))
        end
        flag = count(==(2), degree(newgraph)) > ignorecounter

        while flag
            n = nv(newgraph)
            ignorecounter = 0
            for i in 1:n
                if degree(newgraph, i) == 2
                    neigh1, neigh2 = neighbors(newgraph, i)
                    if neigh1.v == i
                        ignorecounter += 1
                        continue
                    end
                    add_edge!(newgraph, PeriodicEdge{N}(neigh1.v, neigh2.v, neigh2.ofs .- neigh1.ofs))
                    push!(remove_idx, i)
                end
            end
            map = rem_vertices!(newgraph, remove_idx)
            vmap = vmap[map]
            #pos = pos[vmap]
            empty!(remove_idx)
            flag = count(==(2), degree(newgraph)) > ignorecounter
        end
        flag = any(<=(1), degree(newgraph))
    end
    return vmap, newgraph
end


"""
    CrystalNet{D,T<:Real}

Representation of a net as a topological abstraction of a crystal.

`D` is the dimensionality of the net, which is the number of repeated dimensions
of a single connex component. This dimensionality is not necessarily the dimension
of the space the crystal is embedded into, which would always be 3 for real space.

`T` is the numeric type used to store the exact coordinates of each vertex at the
equilibrium placement.
"""
struct CrystalNet{D,T<:Union{Rational{Int32}, Rational{Int64}, Rational{Int128}, Rational{BigInt}}}
    pge::PeriodicGraphEmbedding{D,T}
    types::Vector{Symbol}
    options::Options
end

CrystalNet{D,T}(net::CrystalNet{D,T}) where {D,T} = net
function CrystalNet{D,T}(net::CrystalNet{N}; kwargs...) where {D,T,N}
    n = length(net.types)
    if N != D
        newpos = Vector{SVector{D,T}}(undef, n)
        for i in 1:n
            if N < D
                _pos = zero(MVector{D,T})
                _pos[1:N] = net.pge.pos[i]
                pos = SVector{D,T}(_pos)
            else
                pos = SVector{D,T}(net.pge.pos[i][1:D])
            end
            newpos[i] = pos
        end
    else
        newpos = net.pge.pos
    end
    pge = PeriodicGraphEmbedding{D,T}(PeriodicGraph{D}(net.pge.g), newpos, net.pge.cell)
    CrystalNet{D,T}(pge, net.types, Options(net.options; kwargs...))
end
CrystalNet{D}(net::CrystalNet{R,T}; kwargs...) where {D,R,T} = CrystalNet{D,T}(net; kwargs...)

const CrystalNet1D = CrystalNet{1}
const CrystalNet2D = CrystalNet{2}
const CrystalNet3D = CrystalNet{3}

Base.ndims(::CrystalNet{D}) where {D} = D

function CrystalNet{D,T}(cell::Cell, types::AbstractVector{Symbol}, graph::PeriodicGraph{D},
                       placement::AbstractMatrix{T}, options::Options) where {D,T<:Real}
    n = nv(graph)
    @toggleassert size(placement) == (D, n)
    pge, s = SortedPeriodicGraphEmbedding{T}(graph, placement, cell)
    types = Symbol[types[s[i]] for i in 1:n]
    # @toggleassert all(pos[i] == mean(pos[x.v] .+ x.ofs for x in neighbors(graph, i)) for i in 1:length(pos))
    return CrystalNet{D,T}(pge, types, options)
end

function CrystalNet{D,T}(cell::Cell, opts::Options) where {D,T<:Real}
    return CrystalNet{D,T}(PeriodicGraphEmbedding{D,T}(cell), Symbol[], opts)
end
CrystalNet{D}(cell::Cell, opts::Options) where {D} = CrystalNet{D,Rational{Int32}}(cell, opts)

function CrystalNet{D,T}(cell::Cell, types::Vector{Symbol}, pos, graph::PeriodicGraph{D}, options::Options) where {D,T}
    CrystalNet{D,T}(PeriodicGraphEmbedding{D,T}(graph, pos, cell), types, options)
end


function CrystalNet{D}(pge::PeriodicGraphEmbedding{D,T}, types::Vector{Symbol}, options::Options) where {D,T}
    CrystalNet{D,T}(pge, types, options)
end
function CrystalNet{D}(cell::Cell, types::Vector{Symbol},
                       graph::PeriodicGraph{D}, options::Options) where D
    placement = equilibrium(graph)
    pge, s = SortedPeriodicGraphEmbedding(graph, placement, cell)
    return CrystalNet{D}(pge, types[s], options)
end

function CrystalNet{D}(graph::PeriodicGraph{D}, options::Options) where D
    CrystalNet{D}(Cell(), fill(Symbol(""), nv(graph)), graph, options)
end
function CrystalNet(graph::PeriodicGraph{D}, options::Options) where D
    CrystalNet{D}(Cell(), fill(Symbol(""), nv(graph)), graph, options)
end


const PseudoGraph{D} = Union{PeriodicGraph{D},AbstractVector{PeriodicEdge{D}}}

function CrystalNet(g::Union{PseudoGraph,AbstractString}; kwargs...)
    CrystalNet(g isa PeriodicGraph ? g : PeriodicGraph(g), Options(; kwargs...))
end
function CrystalNet{D}(g::Union{PseudoGraph,AbstractString}; kwargs...) where {D}
    CrystalNet{D}(g isa PeriodicGraph ? g : PeriodicGraph{D}(g), Options(; kwargs...))
end


function CrystalNet{D}(cell::Cell, types::Vector{Symbol},
                       graph::PeriodicGraph, options::Options) where D
    ne(graph) == 0 && return CrystalNet{D}(cell, types, PeriodicGraph{D}(nv(graph)), options)
    return CrystalNet{D}(cell, types, PeriodicGraph{D}(graph), options)
end


function Base.show(io::IO, x::CrystalNet)
    print(io, typeof(x), " of ", x.options.name, " with ", length(x.types), " vertices and ",
          ne(x.pge.g), " edges")
    if length(x.options.clusterings) == 1
        print(io, " (clustering: ", x.options.clusterings[1], ')')
    end
    if !isempty(x.options.error)
        print(io, " (an error happened: ", x.options.error, ')')
    end
    nothing
end

function separate_components(c::Crystal{T}) where T
    graph = PeriodicGraph3D(c.pge.g)
    dimensions = PeriodicGraphs.dimensionality(graph)
    @ifwarn if haskey(dimensions, 0)
        @warn "Detected structure of dimension 0, possibly solvent residues. It will be ignored for topology computation."
    end
    ret = (Tuple{Vector{Int},Crystal{T}}[], Tuple{Vector{Int},Crystal{T}}[], Tuple{Vector{Int},Crystal{T}}[])
    for i in 1:3
        reti = ret[i]
        for vmap in get(dimensions, i, Vector{Int}[])
            push!(reti, (vmap, c[vmap]))
        end
    end
    return ret
end


function _collect_net!(ret::Vector{<:CrystalNet{D}}, encountered, idx, c, clustering) where D
    vmap, graph = trim_topology(c.pge.g)
    types = c.types[vmap]
    remove_metal_cluster_bonds!(graph, types, c.options)
    remove_homoatomic_bonds!(graph, types, c.options.ignore_homoatomic_bonds, false)
    j = get!(encountered, c.pge.g, idx)
    if j == idx
        export_default(Crystal{Nothing}(c.pge.cell, types, c.pge.pos[vmap], graph, c.options),
            lazy"subnet_$clustering", c.options.name, c.options.export_subnets)
        ret[idx] = try
            CrystalNet{D}(c.pge.cell, types, graph, c.options)
        catch e
            (c.options.throw_error || isinterrupt(e)) && rethrow()
            CrystalNet{D}(c.pge.cell, Options(c.options; error=(string(e)::String)))
        end
    else
        ref = ret[j]
        ret[idx] = typeof(ref)(ref.pge.cell, ref.types, ref.pge.pos, ref.pge.g, c.options)
    end
    nothing
end

function collect_nets(crystals::Vector{Crystal{Nothing}}, ::Val{D}) where D
    ret = Vector{CrystalNet{D}}(undef, length(crystals))
    encountered = Dict{PeriodicGraph3D,Int}()
    idx = 1
    for c in crystals
        clustering = only(c.options.clusterings)
        structure = c.options.structure
        if clustering == Clustering.Auto && (structure == StructureType.MOF || structure == StructureType.Cluster)
            alln = Crystal{Nothing}(c; clusterings=[Clustering.AllNodes])
            pge = PeriodicGraphEmbedding(c.pge.g, c.pge.pos, c.pge.cell)
            singlen = allnodes_to_singlenodes(Crystal{Nothing}(pge, c.types, Options(c.options; clusterings=[Clustering.SingleNodes])))
            resize!(ret, length(ret)+1)
            _collect_net!(ret, encountered, idx, alln, Clustering.AllNodes)
            _collect_net!(ret, encountered, idx+1, singlen, Clustering.SingleNodes)
            idx += 1
        elseif clustering == Clustering.PE
            _collect_net!(ret, encountered, idx, pem_to_pe(c), clustering)
        elseif clustering == Clustering.SingleNodes || clustering == Clustering.Standard # Standard = SingleNode ∘ PEM
            _collect_net!(ret, encountered, idx, allnodes_to_singlenodes(c), clustering)
        else
            _collect_net!(ret, encountered, idx, c, clustering)
            export_default(c, lazy"clusters_$clustering", c.options.name, c.options.export_clusters)
        end
        idx += 1
    end
    return ret
end

"""
    UnderlyingNets

Grouping of the connected components of a structure according to their dimensionality.
"""
struct UnderlyingNets
    D1::Vector{Tuple{Vector{Int},Vector{CrystalNet1D}}}
    D2::Vector{Tuple{Vector{Int},Vector{CrystalNet2D}}}
    D3::Vector{Tuple{Vector{Int},Vector{CrystalNet3D}}}
end
UnderlyingNets() = UnderlyingNets(Tuple{Vector{Int},Vector{CrystalNet1D}}[],
                                  Tuple{Vector{Int},Vector{CrystalNet2D}}[],
                                  Tuple{Vector{Int},Vector{CrystalNet1D}}[],
                                 )

function _repeatgroups!(ex, i)
    for j in 1:length(ex.args)
        arg = ex.args[j]
        if arg isa Symbol
            if arg === :groups
                ex.args[j] = Expr(:., :groups, QuoteNode(Symbol(:D, i)))
            elseif arg === :nets
                ex.args[j] = Symbol(:nets, i)
            elseif arg === :D
                ex.args[j] = i
            # elseif arg === :CrystalNet
            #     ex.args[j] = :(CrystalNet{$i})
            # elseif arg === :UnderlyingNets
            #     ex.args[j] = :groups
            end
        elseif arg isa Expr
            _repeatgroups!(arg, i)
        end
    end
    nothing
end

macro repeatgroups(ex)
    exs = [deepcopy(ex) for _ in 1:3]
    for i in 1:3
        _repeatgroups!(exs[i], i)
    end
    return quote
        $(esc(exs[1]))
        $(esc(exs[2]))
        $(esc(exs[3]))
    end
end

function UnderlyingNets(c::Crystal)
    groups = UnderlyingNets()
    components = separate_components(c)
    if all(isempty, components)
        vmap = collect(1:length(c.types))
        nets = [CrystalNet3D(c.pge.cell, Options(c.options; clusterings=[clust])) for clust in c.options.clusterings]
        push!(groups.D3, (vmap, nets))
        return groups
    end
    @repeatgroups begin
        for (vmap, component) in components[D]
            crystals = collapse_clusters(component)
            nets = collect_nets(crystals, Val(D))
            push!(groups, (vmap, nets))
        end
    end
    return groups
end


function CrystalNet(c::Crystal)
    group = UnderlyingNets(c)
    D = isempty(group.D3) ? isempty(group.D2) ? isempty(group.D1) ? 0 : 1 : 2 : 3
    D == 0 && return CrystalNet{0}(c[1].pge.cell, c[1].options)
    if D == 3
        length(group.D3) > 1 && __throw_interpenetrating(D)
        (isempty(group.D1) && isempty(group.D2)) || __warn_nonunique(D)
        _D3 = last(first(group.D3))
        length(_D3) > 1 && __throw_multiplenets(D)
        return first(_D3)
    elseif D == 2
        length(group.D2) > 1 && __throw_interpenetrating(D)
        isempty(group.D2) || __warn_nonunique(D)
        _D2 = last(first(group.D2))
        length(_D2) > 1 && __throw_multiplenets(D)
        return first(_D2)
    end
    @toggleassert D == 1
    length(group.D1) > 1 && __throw_interpenetrating(D)
    _D1 = last(first(group.D1))
    length(_D1) > 1 && __throw_multiplenets(D)
    return first(_D1)
end
__warn_nonunique(D) = @ifwarn @warn "Presence of periodic structures of different dimensionalities. Only the highest dimensionality ($D here) will be retained."
__throw_interpenetrating(D) = error(ArgumentError("Multiple interpenetrating $D-dimensional structures. Cannot handle this as a single CrystalNet, use UnderlyingNets instead."))
__throw_multiplenets(D) = error(ArgumentError("Found multiple nets of dimension $D, please specify a single `Clustering` option."))


function CrystalNet(x::UnderlyingNets)
    if isempty(x.D3)
        if isempty(x.D2)
            if isempty(x.D1)
                CrystalNet{0}(Cell(), Options(; error="Empty UnderlyingNets cannot be converted to a CrystalNet."))
            end
            if length(x.D1) == 1
                return x.D1[1][2][1]
            end
        end
        if length(x.D2) == 1 && isempty(x.D1)
            return x.D2[1][2][1]
        end
    end
    if length(x.D3) == 1 && isempty(x.D2) && isempty(x.D1)
        return x.D3[1][2][1]
    end
    CrystalNet{0}(Cell(), Options(; error="UnderlyingNets contain multiple nets, cannot be converted to a single CrystalNet."))
end
function CrystalNet{D}(x::UnderlyingNets) where D
    if D == 3
        length(x.D3) == 1 && return x.D3[1][2][1]
        return CrystalNet{D}(Cell(), Options(; error="UnderlyingNets contains $(ifelse(isempty(x.D3), "no net", "multiple nets")) of dimension 3"))
    elseif D == 2
        length(x.D2) == 1 && return x.D2[1][2][1]
        return CrystalNet{D}(Cell(), Options(; error="UnderlyingNets contains $(ifelse(isempty(x.D2), "no net", "multiple nets")) of dimension 2"))
    elseif D == 1
        length(x.D1) == 1 && return x.D1[1][2][1]
        return CrystalNet{D}(Cell(), Options(; error="UnderlyingNets contains $(ifelse(isempty(x.D1), "no net", "multiple nets")) of dimension 1"))
    end
    CrystalNet{D}(Cell(), Options(; error="UnderlyingNets cannot be converted to a CrystalNet of dimension $D."))
end


const SmallPseudoGraph = Union{PseudoGraph{1},PseudoGraph{2},PseudoGraph{3}}

function UnderlyingNets(g::SmallPseudoGraph, options::Options)
    groups = UnderlyingNets()
    graph = g isa PeriodicGraph ? g : PeriodicGraph(g)
    dimensions = PeriodicGraphs.dimensionality(graph)
    @ifwarn if haskey(dimensions, 0)
        @error "Detected substructure of dimension 0 in the input graph. It will be ignored for topology computation."
    end
    cell = Cell()
    @repeatgroups begin
        for vmap in get(dimensions, D, Vector{Int}[])
            nets = PeriodicGraph{D}(graph[vmap])
            types = fill(Symbol(""), nv(nets))
            push!(groups, (vmap, CrystalNet{D}[CrystalNet{D}(cell, types, nets, options)]))
        end
    end
    return groups
end
UnderlyingNets(g::AbstractString, opts::Options) = UnderlyingNets(PeriodicGraph(g), opts)
UnderlyingNets(g::Union{SmallPseudoGraph,AbstractString}; kwargs...) = UnderlyingNets(g, Options(; kwargs...))

const SmallDimPeriodicGraph = Union{PeriodicGraph{0}, PeriodicGraph1D, PeriodicGraph2D, PeriodicGraph3D}

"""
    TopologicalGenome

A topological genome computed by `CrystalNets.jl`.

Store both the actual genome (as a `PeriodicGraph`) and the name of the net, if recognized.

Like for a `PeriodicGraph`, the textual representation of a `TopologicalGenome` can be
parsed back into a `TopologicalGenome`:
```jldoctest
julia> topology = topological_genome(CrystalNet(PeriodicGraph("2  1 2 0 0  2 1 0 1  2 1 1 0")))
hcb

julia> typeof(topology)
TopologicalGenome

julia> PeriodicGraph(topology)  # The actual topological genome, as a PeriodicGraph
PeriodicGraph2D(2, PeriodicEdge2D[(1, 2, (-1,0)), (1, 2, (0,0)), (1, 2, (0,1))])

julia> parse(TopologicalGenome, "hcb") == topology
true
```
"""
struct TopologicalGenome
    genome::SmallDimPeriodicGraph
    name::Union{Nothing,String}
    unstable::Bool
    error::String
end

TopologicalGenome(g, name, unstable=false) = TopologicalGenome(g, name, unstable, "")
TopologicalGenome(x::AbstractString) = TopologicalGenome(PeriodicGraph{0}(), nothing, false, x)
TopologicalGenome() = TopologicalGenome("")

function ==(s1::TopologicalGenome, s2::TopologicalGenome)
    s1.genome == s2.genome || return false
    ndims(s1.genome) == 0 && return s1.error == s2.error
    return true
end
function Base.hash(s::TopologicalGenome, h::UInt)
    isempty(s.error) ? hash(s.error, h) : s.name === nothing ? hash(s.genome, h) : hash(s.name, h)
end

function Base.show(io::IO, x::TopologicalGenome)
    if !isempty(x.error)
        print(io, "FAILED with: ", x.error)
    elseif ndims(x.genome) == 0
        print(io, "non-periodic")
    elseif x.unstable
        print(io, "unstable ", x.genome)
    elseif x.name isa String
        splitcomma = split(x.name, ','; limit=2)
        fstsplit = first(splitcomma)
        if length(fstsplit) ≥ 3 && islowercase(fstsplit[3]) && !startswith(fstsplit, "sqc")
            printstyled(io, fstsplit; bold=true) # RCSR
        else
            print(io, fstsplit)
        end
        length(splitcomma) == 2 && print(io, ',', last(splitcomma))
    else
        print(io, "UNKNOWN ", x.genome)
    end
end

function Base.parse(::Type{TopologicalGenome}, s::AbstractString)
    if startswith(s, "UNKNOWN")
        return TopologicalGenome(PeriodicGraph(s[9:end]), nothing, false)
    end
    s == "non-periodic" && return TopologicalGenome()
    if startswith(s, "unstable")
        return TopologicalGenome(PeriodicGraph(s[10:end]), nothing, true)
    end
    if startswith(s, "FAILED")
        return TopologicalGenome(s[13:end])
    end
    return TopologicalGenome(parse(PeriodicGraph, REVERSE_CRYSTAL_NETS_ARCHIVE[s]), s, false)
end


PeriodicGraphs.PeriodicGraph(t::TopologicalGenome) = t.genome
PeriodicGraphs.PeriodicGraph{D}(t::TopologicalGenome) where D = PeriodicGraph{D}(t.genome)

"""
    TopologyResult

The result of a topology computation on a structure with different [`Clustering`](@ref)
options.

Its representation includes the name of the clustering options along with their
corresponding genome. It is omitted if there is only one clustering option which is
[`Auto`](@ref Clustering).

Like for a [`TopologicalGenome`](@ref) (or a `PeriodicGraph`), the textual representation
of a `TopologyResult` can be parsed back to a `TopologyResult`:
```jldoctest
julia> mof5 = joinpath(dirname(dirname(pathof(CrystalNets))), "test", "cif", "MOF-5.cif");

julia> topologies = determine_topology(mof5, structure=StructureType.MOF, clusterings=[Clustering.Auto, Clustering.Standard, Clustering.PE])
AllNodes, SingleNodes: pcu
Standard: xbh
PE: cab

julia> typeof(topologies)
TopologyResult

julia> parse(TopologyResult, repr(topologies)) == topologies
true
```
"""
struct TopologyResult
    results::SizedVector{8,TopologicalGenome,Vector{TopologicalGenome}}
    attributions::MVector{8,Int8}
    uniques::Vector{Int8}
end

TopologyResult() = TopologyResult(SizedVector{8,TopologicalGenome}(undef),
                                  MVector{8,Int8}(0, 0, 0, 0, 0, 0, 0, 0), Int8[])

function TopologyResult(x::Vector{Tuple{_Clustering,Union{_Clustering,TopologicalGenome}}})
    ret = TopologyResult()
    results = ret.results
    attributions = ret.attributions
    uniques = ret.uniques
    sort!(x; by=t -> t[2] isa _Clustering)

    encountered = Dict{TopologicalGenome,Int}()
    for (_i, result) in x
        if result isa _Clustering
            j = Int(result)
            attributions[Int(_i)] = attributions[j]
        else
            i = Int(_i)
            j = get!(encountered, result, i)
            attributions[i] = j
            if j == i
                push!(uniques, j)
                results[j] = result
            end
        end
    end

    return ret
end

function ==(t1::TopologyResult, t2::TopologyResult)
    for (a1, a2) in zip(t1.attributions, t2.attributions)
        if a1 * a2 == 0
            a1 + a2 == 0 || return false
        else
            t1.results[a1] == t2.results[a2] || return false
        end
    end
    return true
end

function Base.hash(t::TopologyResult, h::UInt)
    for i in 1:8
        a = t.attributions[i]
        a == 0 && continue
        h = hash(t.results[a], hash(i, h))
    end
    return h
end


@static if VERSION < v"1.7-"
    struct Returns{T} <: Function
        x::T
    end
    (r::Returns)(args...; kwargs...) = r.x
end

function Base.get(f::Union{Function,Type}, x::TopologyResult, c::_Clustering)
    if x.attributions[Int(c)] == 0
        f()
    else
        x.results[x.attributions[Int(c)]]
    end
end
Base.get(x::TopologyResult, c::_Clustering, default) = get(Returns(default), x, c)

function Base.getindex(x::TopologyResult, c::_Clustering)
    if x.attributions[Int(c)] == 0
        throw(ArgumentError(lazy"No stored topology result for clustering $c"))
    end
    return x.results[x.attributions[Int(c)]]
end
Base.getindex(x::TopologyResult, c::Symbol) = getindex(x, clustering_from_symb(c))
Base.getindex(x::TopologyResult, i) = x.results[x.uniques[i]]

function Base.setindex!(x::TopologyResult, a::Union{TopologicalGenome,_Clustering}, c::_Clustering)
    i = Int(c)
    notuniqueflag = false
    if a isa _Clustering
        x.attributions[i] = Int(a)
        notuniqueflag = true
    else
        for j in x.uniques
            if x.results[j] == a
                x.attributions[i] = j
                notuniqueflag = true
                break
            end
        end
    end
    k = findfirst(==(i), x.uniques)
    if notuniqueflag
        if k isa Int
            deleteat!(x.uniques, k)
        end
    else
        if k isa Nothing
            push!(x.uniques, i)
        end
        x.attributions[i] = i
        x.results[i] = a::TopologicalGenome
    end
    x
end

TopologyResult(s::String) = setindex!(TopologyResult(), TopologicalGenome(s), Clustering.Auto)

function Base.setindex!(x::TopologyResult, ::Nothing, c::_Clustering)
    i = Int(c)
    if x.attributions[i] == i
        j::Int = findfirst(==(i), x.uniques)
        deleteat!(x.uniques, j)
        for k in 1:8
            if x.attributions[k] == i
                x.attributions[k] = 0
            end
        end
    else
        x.attributions[i] = 0
    end
    x
end


function Base.show(io::IO, x::TopologyResult)
    rev_vmap = zeros(Int, 8)
    for (i, j) in enumerate(x.uniques)
        rev_vmap[j] = i
    end
    samenet = [[k] for k in x.uniques]
    for (i, j) in enumerate(x.attributions)
        (j == 0 || j == i) && continue
        push!(samenet[rev_vmap[j]], i)
    end

    if length(samenet) == 1 && length(samenet[1]) == 1 && x.uniques[1] == 1 # Auto
        print(io, x.results[1])
        return
    end

    compact = get(io, :compact, false)
    for (k, l) in enumerate(samenet)
        join(io, (string(clustering_from_num(k)) for k in l), compact ? ',' : ", ")
        print(io, ": ")
        show(io, x.results[l[1]])
        k == length(samenet) || (compact ? print(io, " | ") : println(io))
    end
end

function Base.iterate(x::TopologyResult, state::Int=1)
    next = iterate(x.uniques, state)
    if next isa Tuple{Int8,Int}
        return x.results[next[1]], next[2]
    end
    nothing
end
Base.eltype(::Type{TopologyResult}) = TopologicalGenome
Base.length(x::TopologyResult) = length(x.uniques)

function Base.parse(::Type{TopologyResult}, s::AbstractString)
    splits = split(s, " | ")
    if length(splits) == 1
        splits = split(s, '\n')
    end
    ret = Tuple{_Clustering,Union{_Clustering,TopologicalGenome}}[]
    if length(splits) == 1 && (startswith(s, "FAILED") || length(split(splits[1], ": ")) == 1)
        push!(ret, (Clustering.Auto, parse(TopologicalGenome, s)))
        return TopologyResult(ret)
    end
    for sp in splits
        _clusterings, result = split(sp, ": ")
        clusterings = split(_clusterings, ',')
        ref_cluster = parse(_Clustering, first(clusterings))
        if length(clusterings) > 1
            for j in 2:length(clusterings)
                clust = clusterings[j]
                push!(ret, (parse(_Clustering, clust[1] == ' ' ? clust[2:end] : clust), ref_cluster))
            end
        end
        push!(ret, (ref_cluster, parse(TopologicalGenome, result)))
    end
    return TopologyResult(ret)
end
