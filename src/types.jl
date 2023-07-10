## Type definitions for intermediate representations of crystals up to a net

import Base: ==

function cell_with_warning(mat::StaticArray{Tuple{3,3},BigFloat})
    if !all(isfinite, mat) || iszero(det(mat))
        throw(ArgumentError(lazy"Suspicious unit cell of matrix $(Float64.(mat)). Please check that the unit cell format is given and readable by Chemfiles."))
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
            push!(toremove, i-1)
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
    maxid = maximum(cif.ids; init=1)

    if isempty(cif.bonds)
        oldbonds = Vector{Tuple{Int,Float32}}[]
        knownbondlengths = true
        bonds = Vector{Tuple{Int,Float32}}[]
    else
        @toggleassert minimum(cif.ids) ≥ 1
        rev_id = [Int[] for _ in 1:maxid]
        for (i, k) in enumerate(cif.ids)
            push!(rev_id[k], i)
        end
        # For each id k, oldbonds[k] is the list of (j, d) where j is the id of an atom
        # bonded to k with distance d.
        oldbonds = [Tuple{Int,Float32}[] for _ in 1:maxid]
        for (k, revidk) in enumerate(rev_id)
            bondsk = oldbonds[k]
            for i in revidk
                append!(bondsk, (cif.ids[j], d) for (j, d) in cif.bonds[i])
            end
            sortprune_bondlist!(bondsk)
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
    @toggleassert allunique(cif.ids)
    newids::Vector{Int} = copy(cif.ids)
    newpos::Vector{SVector{3,Float64}} = collect(eachcol(cif.pos))
    smallmat = Float64.(cif.cell.mat)
    symmetric_aliases = [BitSet(i) for i in 1:maxid]
    buffer, ortho, safemin = prepare_periodic_distance_computations(smallmat)
    #=@inbounds=# for equiv in cif.cell.equivalents
        image = zeros(Int, n) # image[i] is the index of the position of the image of i
        for i in 1:n
            v = newpos[i]
            p = Vector(equiv.mat*v .+ equiv.ofs)
            p .-= floor.(p)
            imi = 0
            for (j, posj) in enumerate(newpos)
                buffer .= posj .- p
                if periodic_distance!(buffer, smallmat, ortho, safemin) < 0.55
                    imi = j
                    break
                end
            end
            if imi == 0
                push!(newpos, p)
                push!(newids, cif.ids[i])
                isempty(bonds) || push!(bonds, Tuple{Int,Float32}[])
                imi = length(newpos)
            elseif imi ≤ n
                bitset = symmetric_aliases[cif.ids[i]]
                id_imi = cif.ids[imi]
                if symmetric_aliases[id_imi] !== bitset
                    union!(bitset, symmetric_aliases[id_imi])
                    symmetric_aliases[id_imi] = bitset
                end
            end
            image[i] = imi
        end
        if !knownbondlengths
            for (i, bondi) in enumerate(cif.bonds), (j, _) in bondi
                imi = image[i]
                imj = image[j]
                add_to_bondlist!(bonds[imi], imj, 0f0)
                add_to_bondlist!(bonds[imj], imi, 0f0)
            end
        end
    end

    if knownbondlengths && !isempty(bonds)
        for i in 1:n
            id_i = cif.ids[i]
            bitset = symmetric_aliases[id_i]
            min_bitset = minimum(bitset)
            if id_i == min_bitset
                new_oldbonds = Tuple{Int,Float32}[]
                for id_j in bitset
                    for (k, d) in oldbonds[id_j]
                        append!(new_oldbonds, (j, d) for j in symmetric_aliases[k])
                    end
                end
                sortprune_bondlist!(new_oldbonds)
                oldbonds[id_i] = new_oldbonds
            else
                oldbonds[id_i] = oldbonds[min_bitset]
            end
        end

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
    return Crystal{Nothing}(pge, types, rev_permute_mapping(opts, vmap, length(c.types)))
end


function Base.getindex(c::Crystal{T}, vmap::AbstractVector{<:Integer}) where T
    types = c.types[vmap]
    opts = rev_permute_mapping(c.options, vmap, length(c.types))
    if T === Nothing
        return Crystal{Nothing}(c.pge[vmap], types, opts)
    else
        return Crystal{Clusters}(c.pge[vmap], types, c.clusters[vmap], opts)
    end
end

function PeriodicGraphs.make_supercell(c::Crystal{Nothing}, t)
    pge = make_supercell(c.pge, t)
    newtypes = repeat(c.types, prod(t))
    return Crystal{Nothing}(pge, newtypes, c.options)
end

## CrystalNet

# For the remainder of the file, we can work in 1D, 2D or 3D

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
    pge, s = SortedPeriodicGraphEmbedding{T}(copy(graph), placement, cell)
    types = Symbol[types[s[i]] for i in 1:n]
    # @toggleassert all(pos[i] == mean(pos[x.v] .+ x.ofs for x in neighbors(graph, i)) for i in 1:length(pos))
    return CrystalNet{D,T}(pge, types, rev_permute_mapping(options, s))
end

function CrystalNet{D,T}(cell::Cell, opts::Options) where {D,T<:Real}
    return CrystalNet{D,T}(PeriodicGraphEmbedding{D,T}(cell), Symbol[], opts)
end
CrystalNet{D}(cell::Cell, opts::Options) where {D} = CrystalNet{D,Rational{Int32}}(cell, opts)

function CrystalNet{D,T}(cell::Cell, types::Vector{Symbol}, pos, graph::PeriodicGraph{D}, options::Options) where {D,T}
    CrystalNet{D,T}(PeriodicGraphEmbedding{D,T}(copy(graph), pos, cell), types, options)
end


function CrystalNet{D}(pge::PeriodicGraphEmbedding{D,T}, types::Vector{Symbol}, options::Options) where {D,T}
    CrystalNet{D,T}(pge, types, options)
end
function CrystalNet{D}(cell::Cell, types::Vector{Symbol},
                       graph::PeriodicGraph{D}, options::Options) where D
    placement = equilibrium(graph) # from PeriodicGraphEquilibriumPlacement.jl
    pge, s = SortedPeriodicGraphEmbedding(copy(graph), placement, cell)
    return CrystalNet{D}(pge, types[s], rev_permute_mapping(options, s))
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

function PeriodicGraphs.make_supercell(net::CrystalNet{N,T}, t) where {N,T}
    N == 0 && return copy(net)
    pge = make_supercell(net.pge, t)
    newtypes = repeat(net.types, prod(t))
    return CrystalNet{N,T}(pge, newtypes, Options(net.options; _pos=SVector{3,Float64}[]))
end

function separate_components(c::Crystal{T}) where T
    graph = PeriodicGraph3D(c.pge.g)
    dimensions = PeriodicGraphs.dimensionality(graph)
    @ifwarn if haskey(dimensions, 0)
        @warn "Detected structure of dimension 0, possibly solvent residues. It will be ignored for topology computation."
    end
    ret = (Tuple{Crystal{T},Int,Vector{Int}}[], Tuple{Crystal{T},Int,Vector{Int}}[], Tuple{Crystal{T},Int,Vector{Int}}[])
    for i in 1:3
        reti = ret[i]
        for (vmap, nfold) in get(dimensions, i, Vector{Int}[])
            push!(reti, (c[vmap], nfold, vmap))
        end
    end
    return ret
end


function _collect_net!(ret::Vector{<:CrystalNet{D}}, encountered, idx, c, clustering) where D
    vmap, graph = trim_topology(c.pge.g)
    types = c.types[vmap]
    opts = rev_permute_mapping(c.options, vmap, length(c.types))
    remove_metal_cluster_bonds!(graph, types, opts)
    remove_homoatomic_bonds!(graph, types, c.options.ignore_homoatomic_bonds, false)
    j = get!(encountered, c.pge.g, idx)
    if j == idx
        export_default(Crystal{Nothing}(c.pge.cell, types, c.pge.pos[vmap], graph, opts),
            lazy"subnet_$clustering", c.options.name, c.options.export_subnets)
        ret[idx] = try
            CrystalNet{D}(c.pge.cell, types, graph, opts)
        catch e
            (c.options.throw_error || isinterrupt(e)) && rethrow()
            CrystalNet{D}(c.pge.cell, Options(opts; error=(string(e)::String)))
        end
    else
        ref = ret[j]
        ret[idx] = typeof(ref)(ref.pge.cell, ref.types, ref.pge.pos, ref.pge.g, opts)
    end
    nothing
end

function collect_nets(crystals::Vector{Crystal{Nothing}}, ::Val{D}) where D
    ret = Vector{CrystalNet{D}}(undef, length(crystals))
    encountered = Dict{PeriodicGraph3D,Int}()
    idx = 1
    clusts = falses(length(instances(Clustering._Clustering)))
    for c in crystals
        clustering = only(c.options.clusterings)
        k = Int(clustering)
        if clusts[k]
            resize!(ret, length(ret)-1)
            continue
        end
        clusts[k] = true
        structure = c.options.structure
        if clustering == Clustering.Auto && (structure == StructureType.MOF || structure == StructureType.Cluster)
            hadalln = clusts[Int(Clustering.AllNodes)]
            hadsinglen = clusts[Int(Clustering.SingleNodes)]
            resize!(ret, length(ret)+1-hadalln-hadsinglen)
            hadalln && hadsinglen && continue
            alln = Crystal{Nothing}(c; clusterings=[Clustering.AllNodes])
            if !hadalln
                clusts[Int(Clustering.AllNodes)] = true
                _collect_net!(ret, encountered, idx, alln, Clustering.AllNodes)
                idx += 1
            end
            if !hadsinglen
                clusts[Int(Clustering.SingleNodes)] = true
                pge = PeriodicGraphEmbedding(c.pge.g, c.pge.pos, c.pge.cell)
                singlen = allnodes_to_singlenodes(Crystal{Nothing}(pge, c.types, Options(c.options; clusterings=[Clustering.SingleNodes])))
                _collect_net!(ret, encountered, idx, singlen, Clustering.SingleNodes)
            end
            !hadalln && hadsinglen && (idx -= 1)
        elseif clustering == Clustering.PE
            _collect_net!(ret, encountered, idx, pem_to_pe(c), clustering)
        elseif clustering == Clustering.SingleNodes || clustering == Clustering.Standard # Standard = SingleNodes ∘ PEM
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
    D1::Vector{Tuple{Vector{CrystalNet1D},Int,Vector{Int}}}
    D2::Vector{Tuple{Vector{CrystalNet2D},Int,Vector{Int}}}
    D3::Vector{Tuple{Vector{CrystalNet3D},Int,Vector{Int}}}
end
UnderlyingNets() = UnderlyingNets(Tuple{Vector{CrystalNet1D},Int,Vector{Int}}[],
                                  Tuple{Vector{CrystalNet2D},Int,Vector{Int}}[],
                                  Tuple{Vector{CrystalNet1D},Int,Vector{Int}}[],
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
        push!(groups.D3, (nets, 1, vmap))
        return groups
    end
    @repeatgroups begin
        for (i, (comp, nfold, vmap)) in enumerate(components[D])
            component = Crystal(comp.pge, comp.types, comp.clusters,
                                Options(comp.options; name=string(comp.options.name,'_',i)))
            crystals = collapse_clusters(component)
            nets = collect_nets(crystals, Val(D))
            push!(groups, (nets, nfold, vmap))
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
        _D3 = first(first(group.D3))
        length(_D3) > 1 && __throw_multiplenets(D)
        return first(_D3)
    elseif D == 2
        length(group.D2) > 1 && __throw_interpenetrating(D)
        isempty(group.D2) || __warn_nonunique(D)
        _D2 = first(first(group.D2))
        length(_D2) > 1 && __throw_multiplenets(D)
        return first(_D2)
    end
    @toggleassert D == 1
    length(group.D1) > 1 && __throw_interpenetrating(D)
    _D1 = first(first(group.D1))
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
        for (vmap, nfold) in get(dimensions, D, Vector{Int}[])
            nets = PeriodicGraph{D}(graph[vmap])
            n = nv(nets)
            types = fill(Symbol(""), n)
            opts = rev_permute_mapping(options, vmap, n)
            push!(groups, (CrystalNet{D}[CrystalNet{D}(cell, types, nets, opts)], nfold, vmap))
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
        return TopologicalGenome(s[14:end])
    end
    return TopologicalGenome(parse(PeriodicGraph, REVERSE_CRYSTALNETS_ARCHIVE[s]), s, false)
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

julia> topologies = only(determine_topology(mof5, structure=StructureType.MOF, clusterings=[Clustering.Auto, Clustering.Standard, Clustering.PE]))[1]
AllNodes, SingleNodes: pcu
Standard: xbh
PE: cab

julia> typeof(topologies)
TopologyResult

julia> parse(TopologyResult, repr(topologies)) == topologies
true
```

See also [`TopologicalGenome`](@ref) and [`InterpenetratedTopologyResult`](@ref).
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


function Base.keys(x::TopologyResult)
    rev_vmap = zeros(Int, 8)
    for (i, j) in enumerate(x.uniques)
        rev_vmap[j] = i
    end
    samenet = [[clustering_from_num(k)] for k in x.uniques]
    for (i, j) in enumerate(x.attributions)
        (j == 0 || j == i) && continue
        push!(samenet[rev_vmap[j]], clustering_from_num(i))
    end
    return samenet
end

function Base.show(io::IO, ::MIME"text/plain", x::TopologyResult)
    samenet = keys(x)
    if length(samenet) == 1 && length(samenet[1]) == 1 && x.uniques[1] == 1 # Auto
        print(io, x.results[1])
        return
    end

    compact = get(io, :compact, false)
    for (k, l) in enumerate(samenet)
        join(io, (string(k) for k in l), compact ? ',' : ", ")
        print(io, ": ")
        show(io, x.results[Int(l[1])])
        k == length(samenet) || (compact ? print(io, " | ") : println(io))
    end
end
Base.show(io::IO, x::TopologyResult) = show(io, MIME"text/plain"(), x)

struct ValuesTopologyResult <: AbstractVector{TopologicalGenome}
    x::TopologyResult
end
Base.values(x::TopologyResult) = ValuesTopologyResult(x)
Base.valtype(::Type{TopologyResult}) = TopologicalGenome
Base.size(vx::ValuesTopologyResult) = (length(vx.x),)
Base.getindex(vx::ValuesTopologyResult, i::Int) = vx.x[i]

function Base.iterate(x::TopologyResult, (p, i)=(pairs(x), nothing))
    next = i isa Nothing ? iterate(p) : iterate(p, i)
    next isa Nothing && return nothing
    return next[1], (p, next[2])
end

Base.eltype(::Type{TopologyResult}) = Pair{Vector{_Clustering},TopologicalGenome}

Base.length(x::TopologyResult) = length(x.uniques)
Base.isempty(x::TopologyResult) = isempty(x.uniques)
Base.firstindex(::TopologyResult) = 1
Base.lastindex(x::TopologyResult) = length(x)

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


"""
    InterpenetratedTopologyResult <: AbstractVector{Tuple{TopologyResult,Int}}

The result of a topology computation on a structure containing possibly several
interpenetrated substructures.

An `InterpenetratedTopologyResult` can be seen as a list of `(topology, n)` pair where
* `topology` is the [`TopologyResult`](@ref) corresponding to the substructures.
* `n` is an integer such that the substructure is composed of an `n`-fold catenated net.

The entire structure can thus be decomposed in a series of substructures, each of them
possibly decomposed into several catenated nets.

!!! info "Vocabulary"
    In this context, *interpenetration* and *catenation* have slightly different meanings:
    - two (or more) substructures are *interpenetrated* if both are present in the unit cell, and
      are composed of vertices that have disjoint numbers. They may or may not all have the
      same topology since they are disjoint and independent subgraphs. For example:
      ```jldoctest
      julia> topological_genome(PeriodicGraph("2   1 1  0 1   2 2  0 1   2 2  1 0"))
      2 interpenetrated substructures:
      ⋅ Subnet 1 → UNKNOWN 1 1 1 1
      ⋅ Subnet 2 → sql
      ```
    - a net is `n`-fold *catenated* if the unit cell of a single connected component of the
      net is `n` times larger than the unit cell of the overall net. In that case, the net
      is actually made of `n` interpenetrating connected components, which all have the
      same topology. For example:
      ```jldoctest
      julia> topological_genome(PeriodicGraph("3   1 1  2 0 0   1 1  0 1 0   1 1  0 0 1"))
      (2-fold) pcu
      ```
    Both may occur inside a single structure, for example:
    ```jldoctest
    julia> topological_genome(PeriodicGraph("2   1 1  0 2   2 2  0 1   2 2  1 0"))
    2 interpenetrated substructures:
    ⋅ Subnet 1 → (2-fold) UNKNOWN 1 1 1 1
    ⋅ Subnet 2 → sql
    ```

# Example
```jldoctest
julia> mof14 = joinpath(dirname(dirname(pathof(CrystalNets))), "test", "cif", "MOFs", "MOF-14.cif");

julia> topologies = determine_topology(mof14, structure=StructureType.MOF, clusterings=[Clustering.Auto, Clustering.Standard, Clustering.PE])
2 interpenetrated substructures:
⋅ Subnet 1 → AllNodes,SingleNodes,Standard: pto | PE: sqc11259
⋅ Subnet 2 → AllNodes,SingleNodes,Standard: pto | PE: sqc11259

julia> typeof(topologies)
InterpenetratedTopologyResult

julia> parse(InterpenetratedTopologyResult, repr(topologies)) == topologies
true

julia> topologies[2]
(AllNodes, SingleNodes, Standard: pto
PE: sqc11259, 1)

julia> topology, n = topologies[2]; # second subnet

julia> n # catenation multiplicity
1

julia> topology
AllNodes, SingleNodes, Standard: pto
PE: sqc11259

julia> typeof(topology)
TopologyResult
```
"""
struct InterpenetratedTopologyResult <: AbstractVector{Tuple{TopologyResult,Int}}
    data::Vector{Tuple{TopologyResult,Int,Vector{Int}}}
end
InterpenetratedTopologyResult() = InterpenetratedTopologyResult(Tuple{TopologyResult,Int,Vector{Int}}[])
InterpenetratedTopologyResult(e::AbstractString) = InterpenetratedTopologyResult([(TopologyResult(string(e)), 1, Int[])])
Base.size(x::InterpenetratedTopologyResult) = (length(x.data),)
Base.getindex(x::InterpenetratedTopologyResult, i) = (y = x.data[i]; (y[1], y[2]))

function Base.show(io::IO, ::MIME"text/plain", x::InterpenetratedTopologyResult)
    compact = length(x) > 1
    if compact
        print(io, length(x), " interpenetrated substructures:")
    elseif length(x) == 0
        print(io, "non-periodic")
    end
    for (i, (topology, nfold)) in enumerate(x)
        if compact
            print(io, "\n⋅ Subnet ", i, " → ")
        end
        hasnfold = nfold > 1
        if hasnfold
            @static if VERSION < v"1.10-"
                printstyled(io, '(', nfold, "-fold) ", color=:yellow)
            else
                printstyled(io, '(', nfold, "-fold) ", italic=true)
            end
        end
        print(IOContext(io, :compact=>(compact|hasnfold)), topology)
    end
end
Base.show(io::IO, x::InterpenetratedTopologyResult) = show(io, MIME("text/plain"), x)

function parse_nfold_topologyresult(x::AbstractString)
    nfold = 1
    num_digits = -8
    if x[1] == '('
        nfold = parse(Int, first(split(@view(x[2:end]), !isnumeric; limit=2)))
        num_digits = ndigits(nfold)
        @assert @view(x[(2+num_digits):(8+num_digits)]) == "-fold) "
    end
    parse(TopologyResult, @view(x[(9+num_digits):end])), nfold
end

function Base.parse(::Type{InterpenetratedTopologyResult}, x::AbstractString)
    s = split(x; limit=4)
    length(s) == 1 && x == "non-periodic" && return InterpenetratedTopologyResult()
    if length(s) > 3 && s[2] == "interpenetrated" && s[3] == "substructures:"
        lines = split(s[4], '\n')
        data = Vector{Tuple{TopologyResult,Int,Vector{Int}}}(undef, length(lines))
        for l in lines
            @assert @view(l[1:11]) == "⋅ Subnet "
            splits = split(@view(l[11:end]); limit=3)
            i = parse(Int, splits[1])
            @assert splits[2] == "→"
            topo, nfold = parse_nfold_topologyresult(splits[3])
            data[i] = (topo, nfold, Int[])
        end
        return InterpenetratedTopologyResult(data)
    end
    topo1, nfold1 = parse_nfold_topologyresult(x)
    return InterpenetratedTopologyResult([(topo1, nfold1, Int[])])
end
