## Type definitions for intermediate representations of crystals up to a net

include("specialsolver.jl")
include("options.jl")

import Base: ==
using Tokenize



## EquivalentPosition

"""
    EquivalentPosition

Representation of a symmetry operation in 3D, defined by an affine function.
"""
struct EquivalentPosition
    mat::SMatrix{3,3,Rational{Int}, 9}
    ofs::SVector{3,Rational{Int}}
end

"""
Find the reference identifiers for the three dimensions for the CIF group called
symmetry_equiv_pos_as_xyz or space_group_symop_operation_xyz.
Usually this is simply ("x", "y", "z").
"""
function find_refid(eqs)
    isempty(eqs) && return ("x", "y", "z")
    for eq in eqs # Normally it should be the first one but I'm not sure of the specs here
        ts = collect(tokenize(eq))
        any(Tokenize.Tokens.isoperator∘Tokenize.Tokens.kind, ts) && continue
        refid = [""]
        not_at_the_end = true
        for t in ts
            @assert not_at_the_end
            k = Tokenize.Tokens.kind(t)
            if k === Tokenize.Tokens.IDENTIFIER
                @assert refid[end] == ""
                refid[end] = Tokenize.Tokens.untokenize(t)
            elseif k === Tokenize.Tokens.COMMA || k === Tokenize.Tokens.SEMICOLON
                @assert refid[end] != ""
                push!(refid, "")
            elseif k === Tokenize.Tokens.ENDMARKER
                not_at_the_end = false
            elseif k !== Tokenize.Tokens.WHITESPACE
                error("Input string {$eq} is not a valid symmetry equivalent")
            end
        end
        if not_at_the_end
            error("Unknown end of line marker for symmetry equivalent {$eq}")
        end
        if length(refid) != 3 || refid[end] == ""
            error("Input string {$eq} is not a valid symmetry equivalent")
        end
        return tuple(refid[1], refid[2], refid[3])
    end
    return ("x", "y", "z")
end

function Base.parse(::Type{EquivalentPosition}, s::AbstractString, refid=("x", "y", "z"))
    const_dict = Dict{String, Int}(refid[1]=>1, refid[2]=>2, refid[3]=>3)
    mat = zeros(Rational{Int}, 3, 3)
    ofs = zeros(Rational{Int}, 3)
    curr_num::Union{Int, Nothing} = nothing
    curr_val::Union{Rational{Int}, Nothing} = nothing
    curr_sign::Union{Bool, Nothing} = nothing
    encountered_div::Bool = false
    i = 1
    something_written = false
    for x in tokenize(s)
        k = Tokenize.Tokens.kind(x)
        k === Tokenize.Tokens.WHITESPACE && continue
        if k === Tokenize.Tokens.INTEGER
            @assert isnothing(curr_val)
            if encountered_div
                curr_val = Rational{Int}(Int(curr_num), parse(Int, x.val))
                curr_num = nothing
                encountered_div = false
            else
                @assert isnothing(curr_num)
                curr_num = parse(Int, x.val)
            end
        else
            @assert !encountered_div
            if k === Tokenize.Tokens.IDENTIFIER
                if !isnothing(curr_num)
                    @assert isnothing(curr_val)
                    curr_val = curr_num
                    curr_num = nothing
                end
                sign = isnothing(curr_sign) ? 1 : 2*curr_sign - 1
                val = isnothing(curr_val)  ? Rational{Int}(1) : curr_val
                j = const_dict[Tokenize.Tokens.untokenize(x)]
                mat[i,j] += sign * val
                curr_val = nothing
                curr_sign = nothing
                something_written = true
            else
                if x.kind === Tokenize.Tokens.FWD_SLASH
                    @assert isnothing(curr_val)
                    @assert !isnothing(curr_num)
                    encountered_div = true
                    continue
                end
                if !isnothing(curr_num)
                    @assert isnothing(curr_val)
                    curr_val = curr_num
                    curr_num = nothing
                end
                if !isnothing(curr_val)
                    sign = isnothing(curr_sign) ? 1 : 2*curr_sign - 1
                    @ifwarn begin if !iszero(ofs[i])
                            @warn "Existing offset already existing for position $i in {$s}"
                        end
                    end
                    ofs[i] += sign * Rational{Int}(curr_val)
                    curr_val = nothing
                    curr_sign = nothing
                else
                    @assert isnothing(curr_sign)
                end
                if x.kind === Tokenize.Tokens.PLUS
                    curr_sign = true
                elseif x.kind === Tokenize.Tokens.MINUS
                    curr_sign = false
                elseif k === Tokenize.Tokens.COMMA || k === Tokenize.Tokens.SEMICOLON
                    i > 2 && error("Too many dimensions specified for symmetry equivalent {$s}")
                    something_written || error("{$s} is not a symmetry equivalent (no dependency expressed in position $i)")
                    something_written = false
                    i += 1
                else
                    k !== Tokenize.Tokens.ENDMARKER && error("Unknown end of line marker for symmetry equivalent {$s}")
                    i!= 3 && error("Input string {$s} is not a valid symmetry equivalent")
                end
            end
        end
    end
    EquivalentPosition(SMatrix{3, 3, Rational{Int}, 9}(mat), SVector{3, Rational{Int}}(ofs))
end

function Base.show(io::IO, eq::EquivalentPosition)
    function rationaltostring(x::Rational{<:Integer}, notofs::Bool, first::Bool)
        if notofs && (x == 1 || x == -1)
            return x < 0 ? '-' : first ? "" : "+"
        end
        sign = x < 0 || first ? "" : "+"
        sign * (x.den == 1 ? string(x.num) : string(x.num)*'/'*string(x.den))
    end
    xyz = ('x', 'y', 'z')
    for i in 1:3
        first = true
        for j in 1:3
            if eq.mat[i,j] != 0
                coeff = rationaltostring(eq.mat[i,j], true, first)
                first = false
                print(io, coeff)
                print(io, xyz[j])
            end
        end
        if eq.ofs[i] != 0
            print(io, rationaltostring(eq.ofs[i], false, first))
        end
        i < 3 && print(io, ',')
    end
    nothing
end

## Cell

"""
    Cell

Representation of a periodic cell in 3D. Contains information about the cell
(axes lengths and angles) and its symmetry group, through its Hall number.

See SPACE_GROUP_HALL, SPACE_GROUP_FULL, SPACE_GROUP_HM et SPACE_GROUP_IT for
the correspondance between Hall number and usual symbolic representations.
"""
struct Cell
    hall::Int
    mat::SMatrix{3,3,BigFloat,9} # Cartesian coordinates of a, b and c
    equivalents::Vector{EquivalentPosition}

    Cell(hall, mat::AbstractMatrix, equivalents) = new(hall, mat, equivalents)
end

function ==(c1::Cell, c2::Cell)
    c1.hall == c2.hall && c1.mat == c2.mat
end

function Cell(hall, (a, b, c), (α, β, γ), eqs=EquivalentPosition[])
    cosα = cosd(α); cosβ = cosd(β); cosγ = cosd(γ); sinγ = sind(γ)
    ω = sqrt(1 - cosα^2 - cosβ^2 - cosγ^2 + 2*cosα*cosβ*cosγ)
    mat = SMatrix{3,3,BigFloat,9}([a  b*cosγ  c*cosβ ;
                                  0   b*sinγ  c*(cosα - cosβ*cosγ)/sinγ ;
                                  0   0       c*ω/sinγ ])
    if isempty(eqs)
        eqs = get_symmetry_equivalents(hall)
    end
    return Cell(hall, mat, eqs)
end
Cell() = Cell(1, (10, 10, 10), (90, 90, 90), EquivalentPosition[])

function cell_parameters(mat::StaticArray{Tuple{3,3},BigFloat})
    _a, _b, _c = eachcol(mat)
    a = norm(_a)
    b = norm(_b)
    c = norm(_c)
    α = acosd(_b'_c/(b*c))
    β = acosd(_c'_a/(c*a))
    γ = acosd(_a'_b/(a*b))
    return a, b, c, α, β, γ
end
cell_parameters(cell::Cell) = cell_parameters(cell.mat)

function Cell(cell::Cell, mat::StaticArray{Tuple{3,3},BigFloat})
    return Cell(cell.hall, mat, cell.equivalents)
end

function Cell(mat::StaticArray{Tuple{3,3},BigFloat})
    if !all(isfinite, mat) || iszero(det(mat))
        @ifwarn @error "Suspicious unit cell of matrix $(Float64.(mat)). Is the input really periodic? Using a cubic unit cell instead."
        return Cell()
    end
    return Cell(Cell(), mat)
end

function Base.show(io::IO, cell::Cell)
    a, b, c, α, β, γ = Float64.(cell_parameters(cell))
    print(io, Cell, "($(cell.hall), ($a, $b, $c), ($α, $β, $γ))")
end
function Base.show(io::IO, ::MIME"text/plain", cell::Cell)
    (a, b, c), (α, β, γ) = Float64.(cell_parameters(cell))
    hall_symbol, crystal_system = HALL_SYMBOLS[cell.hall]
    print(io, Cell, " with Hall symbol $hall_symbol ($crystal_system) and parameters")
    print(io, "a=$a, b=$b, c=$c, α=$α, β=$β, γ=$γ")
end

## CIF

"""
    CIF

Representation of a .cif file.
"""
struct CIF
    cifinfo::Dict{String, Union{String, Vector{String}}}
    cell::Cell
    ids::Vector{Int}
    types::Vector{Symbol}
    pos::Matrix{Float64}
    bonds::Matrix{Float32}
end


function keep_atoms(cif::CIF, kept)
    kept_ids = sort!([cif.ids[i] for i in kept])
    unique!(kept_ids)
    idmap = Vector{Int}(undef, length(cif.types)) # upper bound on maximum(kept_ids)
    for (i,x) in enumerate(kept_ids)
        idmap[x] = i
    end
    return CIF(cif.cifinfo, cif.cell, [idmap[cif.ids[i]] for i in kept],
               cif.types[kept_ids], cif.pos[:, kept], cif.bonds[kept, kept])
end


"""
    periodic_distance(u, v, mat)

Distance between points u and v, given as triplets of fractional coordinates, in
a repeating unit cell of matrix mat.
The distance is the shortest between all equivalents of u and v.
"""
function periodic_distance(u, v, mat)
    x = similar(u)
    #=@inbounds=# for i in 1:3
        diff = u[i] - v[i] + 0.5
        x[i] = diff - floor(diff) - 0.5
    end
    return norm(mat*x)
end


"""
    remove_partial_occupancy(::CIF)

Only keep one atom per atom site.
"""
function remove_partial_occupancy(cif::CIF)
    points::Vector{SVector{3,Float64}} = collect(eachcol(cif.pos))
    perm = sortperm(points)
    n = length(perm)
    minimal_length = cbrt(norm(cif.cell.mat))*4e-4 # TODO: find an appropriate value
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
    occupancies = parse.(Float64, get(cif.cifinfo, "atom_site_occupancy", 
                                      ["1.0" for _ in 1:length(cif.types)])::Vector{String})
    toremove = Int[]
    bonds = copy(cif.bonds)
    for alias in toalias
        m = length(alias)
        occup = @view occupancies[alias]
        max_occupancy = maximum(occup)
        with_max_occupancy = [alias[i] for i in 1:m if occup[i] == max_occupancy]
        representative = minimum(with_max_occupancy)
        addtoremove = Int[]
        for a in alias
            a == representative && continue
            push!(addtoremove, a)
            bonds[representative,:] .= min.((@view bonds[a,:]), @view bonds[representative,:])
            bonds[:,representative] .= min.((@view bonds[:,a]), @view bonds[:,representative])
        end
        bonds[representative,representative] = Inf32
        append!(toremove, addtoremove)
    end
    sort!(toremove)
    ids = deleteat!(copy(cif.ids), toremove)
    tokeep = deleteat!(collect(1:n), toremove)
    bonds = bonds[tokeep, tokeep]
    pos = cif.pos[:,tokeep]
    @assert allunique(collect(eachcol(pos)))
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
    for i in 1:n, j in (i+1):n
        if periodic_distance(points[i], points[j], smallmat) < 0.55
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
    bonds = cif.bonds[tokeep, tokeep]
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
    n = length(cif.ids)
    _oldbonds = Dict{Tuple{Int64, Int64}, Float32}((cif.ids[i],cif.ids[j]) => cif.bonds[i,j]
                    for i in 1:n for j in (i+1):n if cif.bonds[i,j] < Inf32)
    knownbondlengths = !any(iszero, values(_oldbonds))
    
    if !isempty(_oldbonds) && !knownbondlengths
        # This means that the bonds are determined as symmetric images of bonds of the
        # asymetric unit, which may or may not be the convention used in the CIF files.
        # Most CIF files having geom_bond_atom_site_labels should have geom_bond_distance anyway.
        @ifwarn @warn "Expanding CIF symmetry without knowing the bond lengths: the resulting bonds might be erroneous"
    end
    newbonds = Set{Tuple{Int,Int}}(keys(_oldbonds))
    oldbonds = copy(_oldbonds)
    for ((i,j), bond) in _oldbonds
        oldbonds[(j,i)] = bond
    end
    newids::Vector{Int} = copy(cif.ids)
    newpos::Vector{SVector{3,Float64}} = collect(eachcol(cif.pos))
    smallmat = Float64.(cif.cell.mat)
    #=@inbounds=# for equiv in cif.cell.equivalents
        image = zeros(Int, n)
        for i in 1:length(cif.ids)
            v = newpos[i]
            p = Vector(equiv.mat*v .+ equiv.ofs)
            p .-= floor.(p)
            for j in 1:length(newpos)
                if periodic_distance(newpos[j], p, smallmat) < 0.55
                    image[i] = j
                    break
                end
            end
            if image[i] == 0
                push!(newpos, p)
                push!(newids, cif.ids[i])
                image[i] = length(newpos)
            end
        end
        if !knownbondlengths
            for (i,j) in keys(oldbonds)
                push!(newbonds, minmax(image[i], image[j]))
            end
        end
    end
    m = length(newids)
    bonds = fill(Inf32, m, m)
    if knownbondlengths
        for i in 1:m
            for j in (i+1):m
                bondlength = get(oldbonds, (newids[i], newids[j]), Inf32)
                bondlength < Inf32 || continue
                if abs(periodic_distance(newpos[i], newpos[j], smallmat) - bondlength) < 0.55
                    bonds[i,j] = bonds[j,i] = bondlength
                end
            end
        end
    else
        for (i,j) in newbonds
            bonds[i,j] = bonds[j,i] = zero(Float32)
        end
    end

    return CIF(cif.cifinfo, deepcopy(cif.cell), newids, copy(cif.types), reduce(hcat, newpos), bonds)
end

"""
    edges_from_bonds(bonds, mat, pos)

Given a bond list `bonds` containing triplets `(a, b, dist)` where atoms `a` and `b` are
bonded if their distance is lower than `dist`, the 3×3 matrix of the cell `mat` and the
Vector{SVector{3,Float64}} `pos` whose elements are the fractional positions of the
atoms, extract the list of PeriodicEdge3D corresponding to the bonds.
Since the adjacency matrix wraps bonds across the boundaries of the cell, the edges
are extracted so that the closest representatives are chosen to form bonds.
"""
function edges_from_bonds(bonds, mat, pos)
    n = length(pos)
    adjacency = zeros(n, n)
    for (a, b, dist) in bonds
        adjacency[a, b] = dist
        adjacency[b, a] = dist
    end
    edges = PeriodicEdge3D[]
    ref_dst = norm(mat*[1, 1, 1])
    for i in 1:n, k in (i+1):n
        maxdist = adjacency[k,i]
        iszero(maxdist) && continue
        offset::Vector{SVector{3, Int}} = []
        old_dst = ref_dst
        for ofsx in -1:1, ofsy in -1:1, ofsz in -1:1 # TODO: optimize with the periodic_distance trick?
            dst = norm(mat * (pos[i] .- (pos[k] .+ (ofsx, ofsy, ofsz))))
            if maxdist == -1f0
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
end

function Clusters(n)
    sbus = [[PeriodicVertex3D(i)] for i in 1:n]
    classes = collect(1:n)
    attributions = collect(1:n)
    offsets = [zero(SVector{3,Int}) for _ in 1:n]
    return Clusters(sbus, classes, attributions, offsets)
end

Base.isempty(c::Clusters) = c.attributions == 1:length(c.attributions)

function Base.getindex(c::Clusters, vmap::AbstractVector{<:Integer})
    rev_vmap = zeros(Int, length(c.attributions))
    for (i, j) in enumerate(vmap)
        rev_vmap[j] = i
    end
    sbus = Vector{PeriodicVertex3D}[]
    sbu_vmap = Int[]
    for (i, sbu) in enumerate(c.sbus)
        newsbu = PeriodicVertex3D[]
        for x in sbu
            y = rev_vmap[x.v]
            y == 0 && continue
            push!(newsbu, PeriodicVertex3D(y, x.ofs))
        end
        if !isempty(sbu)
            push!(sbu_vmap, i)
            push!(sbus, newsbu)
        end
    end
    offsets = [c.offsets[i] for i in vmap]
    attributions = sbu_vmap[c.attributions[vmap]]
    classes = c.classes[sbu_vmap]
    return Clusters(sbus, classes, attributions, offsets)
end

## Crystal

"""
    Crystal

Intermediate representation of a crystal, retaining information on the cell, and the
fractional placement of the atoms and their type, as well as the residues which will be used as
vertices for the computation of the underlying topology.
"""
struct Crystal{T<:Union{Nothing,Clusters}}
    cell::Cell
    types::Vector{Symbol}
    clusters::T
    pos::Vector{SVector{3,Float64}}
    graph::PeriodicGraph3D
    options::Options

    function Crystal{Clusters}(cell, types, clusters, pos, graph, options)
        return new{Clusters}(cell, types, clusters, pos, graph, options)
    end

    function Crystal{Nothing}(cell, types, pos, graph, options)
        return new{Nothing}(cell, types, nothing, pos, graph, options)
    end
end




function Crystal(cell, types, clusters, pos, graph, options)
    if clusters isa Nothing
        return Crystal{Nothing}(cell, types, pos, graph, options)
    end
    return Crystal{Clusters}(cell, types, clusters, pos, graph, options)
end
function Crystal(c::Crystal; kwargs...)
    return Crystal(c.cell, c.types, c.clusters, c.pos, c.graph, Options(c.options; kwargs...))
end

function ==(c1::Crystal{T}, c2::Crystal{T}) where T
    c1.cell == c2.cell && c1.types == c2.types && c1.clusters == c2.clusters &&
    c1.pos == c2.pos && c1.graph == c2.graph && c1.options == c2.options
end

function Crystal{Nothing}(c::Crystal{T}; kwargs...) where T
    if isempty(kwargs)
        if T === Nothing
            return c
        end
        return Crystal{Nothing}(c.cell, c.types, c.pos, c.graph, c.options)
    end
    return Crystal{Nothing}(c.cell, c.types, c.pos, c.graph,
                            Options(c.options; kwargs...))
end
function Crystal{Clusters}(c::Crystal, clusters::Clusters; kwargs...)
    if isempty(kwargs)
        return Crystal{Clusters}(c.cell, c.types, clusters, c.pos, c.graph, c.options)
    end
    return Crystal{Nothing}(c.cell, c.types, clusters, c.pos, c.graph,
                            Options(c.options; kwargs...))
end

trimmed_crystal(c::Crystal{Clusters}) = trimmed_crystal(coalesce_sbus(c, c.clusters))
function trimmed_crystal(c::Crystal{Nothing})
    vmap, graph = trim_topology(c.graph)
    return Crystal{Nothing}(c.cell, c.types[vmap], c.pos[vmap], graph, c.options)
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
    iszero(n) && return Matrix{Rational{Int128}}(undef, N, 0)
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
struct CrystalNet{D,T<:Real}
    cell::Cell
    types::Vector{Symbol}
    pos::Vector{SVector{D,T}}
    graph::PeriodicGraph{D}
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
                _pos[1:N] = net.pos[i]
                pos = SVector{D,T}(_pos)
            else
                pos = SVector{D,T}(net.pos[i][1:D])
            end
            newpos[i] = pos
        end
    else
        newpos = net.pos
    end
    CrystalNet{D,T}(net.cell, net.types, newpos,
                    PeriodicGraphs.change_dimension(PeriodicGraph{D}, net.graph),
                    Options(net.options; kwargs...))
end
CrystalNet{D}(net::CrystalNet{R,T}; kwargs...) where {D,R,T} = CrystalNet{D,T}(net; kwargs...)

const CrystalNet1D = CrystalNet{1}
const CrystalNet2D = CrystalNet{2}
const CrystalNet3D = CrystalNet{3}

Base.ndims(::CrystalNet{D}) where {D} = D

function CrystalNet{D,T}(cell::Cell, types::AbstractVector{Symbol}, graph::PeriodicGraph{D},
                       placement::AbstractMatrix{T}, options::Options) where {D,T<:Real}
    n = nv(graph)
    @assert size(placement) == (D, n)
    pos = Vector{SVector{D,T}}(undef, n)
    offsets = Vector{SVector{D,Int}}(undef, n)
    @inbounds for (i, x) in enumerate(eachcol(placement))
        offsets[i] = floor.(Int, x)
        pos[i] = x .- offsets[i]
    end
    s = sortperm(pos)
    pos = pos[s]
    types = Symbol[types[s[i]] for i in 1:n]
    graph = offset_representatives!(graph, .-offsets)[s]
    # @assert all(pos[i] == mean(pos[x.v] .+ x.ofs for x in neighbors(graph, i)) for i in 1:length(pos))
    return CrystalNet{D,T}(cell, types, pos, graph, options)
end

function trim_crystalnet!(graph, types, tohandle, keep)
    sort!(tohandle)
    toremove = keep ? deleteat!(collect(1:length(types)), tohandle) : tohandle
    vmap = rem_vertices!(graph, toremove)
    return vmap
end

struct CrystalNetGroup
    D1::Vector{Tuple{Vector{Int},CrystalNet1D}}
    D2::Vector{Tuple{Vector{Int},CrystalNet2D}}
    D3::Vector{Tuple{Vector{Int},CrystalNet3D}}
end
CrystalNetGroup() = CrystalNetGroup(Tuple{Vector{Int},CrystalNet1D}[], Tuple{Vector{Int},CrystalNet2D}[], Tuple{Vector{Int},CrystalNet3D}[])

function _separategroups!(ex, groups, i)
    for j in 1:length(ex.args)
        arg = ex.args[j]
        if arg isa Symbol
            if arg === :CrystalNetGroups
                ex.args[j] = groups
            elseif arg === :CrystalNet
                ex.args[j] = :(CrystalNet{$i})
            elseif arg === groups
                ex.args[j] = Expr(:., groups, QuoteNode(Symbol(:D, i)))
            end
        elseif arg isa Expr
            _separategroups!(arg, groups, i)
        end
    end
    nothing
end

macro separategroups(D, groups, ex)
    exs = [deepcopy(ex) for _ in 1:3]
    for i in 1:3
        _separategroups!(exs[i], groups, i)
    end
    return quote
        if $(esc(D)) == 1
            $(esc(exs[1]))
        elseif $(esc(D)) == 2
            $(esc(exs[2]))
        elseif $(esc(D)) == 3
            $(esc(exs[3]))
        else
            throw(AssertionError("1 ≤ D ≤ 3"))
        end
    end
end

struct NonPeriodicInputException <: Exception end
Base.showerror(io::IO, ::NonPeriodicInputException) = print(io, "Non-periodic input")

function CrystalNetGroup(cell::Cell, types::AbstractVector{Symbol},
                         graph::PeriodicGraph, opts::Options)
    initialvmap, graph = trim_topology(PeriodicGraphs.change_dimension(PeriodicGraph3D, graph))
    types = types[initialvmap]
    remove_homoatomic_bonds!(graph, types, opts.ignore_homoatomic_bonds)
    if !isempty(opts.export_net) && !isempty(opts._pos)
        pos = opts._pos[initialvmap]
        export_default(Crystal{Nothing}(cell, types, pos, graph, opts), 
                       "net", opts.name, opts.export_net; repeats=2)
    end
    dimensions = PeriodicGraphs.dimensionality(graph)

    if haskey(dimensions, 0)
        @ifwarn @warn "Removing complex structure of dimension 0, possibly solvent residues."
        # Note that the above warning will not appear for small solvent molecules, because
        # those will have been removed by trim_topology beforehand.
        dim0::Vector{Int} = reduce(vcat, @inbounds dimensions[0]; init=Int[])
        vmap0 = trim_crystalnet!(graph, types, dim0, false)
        types = types[vmap0]
        dimensions = PeriodicGraphs.dimensionality(graph)
        initialvmap = initialvmap[vmap0]
    end
    isempty(dimensions) && throw(NonPeriodicInputException())

    groups = CrystalNetGroup()

    for D in sort(collect(keys(dimensions)))
        dimD = @inbounds dimensions[D]
        if length(dimensions) == 1
            graphD = graph
            typesD = types
            vmapD = initialvmap
        else
            graphD = deepcopy(graph)
            catdimD::Vector{Int} = reduce(vcat, dimD; init=Int[])
            vmapD = trim_crystalnet!(graphD, types, catdimD, true)
            invvmapD = zeros(Int, length(types))
            for i in eachindex(vmapD)
                invvmapD[vmapD[i]] = i
            end
            for l in dimD
                for i in eachindex(l)
                    l[i] = invvmapD[l[i]]
                end
            end
            typesD = types[vmapD]
            vmapD = initialvmap[vmapD]
        end
        for l in dimD # l is a connex component of dimensionality D
            if length(dimD) == 1
                g = graphD
                t = typesD
                vmap = vmapD
            else
                g = deepcopy(graphD)
                tokeep::Vector{Int} = reduce(vcat, l; init=Int[])
                vmap = trim_crystalnet!(g, typesD, tokeep, true)
                t = typesD[vmap]
                vmap = vmapD[vmap]
            end
            @separategroups D groups push!(groups, (vmap, CrystalNet(cell, t, g, opts)))
        end
    end

    return groups
end

function CrystalNet(cell::Cell, types::AbstractVector{Symbol},
                    graph::PeriodicGraph, options::Options)
    group = CrystalNetGroup(cell, types, graph, options)
    D = isempty(group.D3) ? (isempty(group.D2) ? 1 : 2) : 3
    nonuniqueD = D != 1 && (!isempty(group.D1) || (D == 3 && !isempty(group.D2)))
    if nonuniqueD
        @ifwarn @warn "Presence of periodic structures of different dimensionalities. Only the highest dimensionality ($D here) will be retained."
    end
    intertwinnedD = false
    @separategroups D group begin
        intertwinnedD = length(group) > 1
    end
    if intertwinnedD
        error(ArgumentError("Multiple intertwinned $D-dimensional structures. Cannot handle this as a single CrystalNet, use CrystalNetGroup instead."))
    end
    @separategroups D group begin
        return last(only(group))
    end
end

function CrystalNet{D}(cell::Cell, types::AbstractVector{Symbol},
                       graph::PeriodicGraph, options::Options) where D
    g = change_dimension(PeriodicGraph{D}, graph)
    return CrystalNet{D}(cell, types, g, options)
end



macro tryinttype(T)
    tmin = :(typemin($T))
    tmax = :(typemax($T))
    return esc(quote
        if (($tmin <= m) & (M <= $tmax))
            net = CrystalNet{D,Rational{$T}}(cell, types, graph,
                        Rational{$T}.(placement), options)
            return net
        end
    end)
end

function CrystalNet{D}(cell::Cell, types::AbstractVector{Symbol},
                       graph::PeriodicGraph{D}, options::Options) where D
    placement = equilibrium(graph)
    m = min(minimum(numerator.(placement)), minimum(denominator.(placement)))
    M = max(maximum(numerator.(placement)), maximum(denominator.(placement)))
    @tryinttype Int8
    @tryinttype Int16
    @tryinttype Int32
    @tryinttype Int64
    @tryinttype Int128
    return CrystalNet{D,Rational{BigInt}}(cell, types, graph, placement, options)
    # Type-unstable function, but yields better performance than always falling back to BigInt
end

struct EmptyGraphException end
function Base.showerror(io::IO, ::EmptyGraphException)
    print(io, "Empty graph. This probably means that the bonds have not been correctly attributed. Please switch to an input file containing explicit bonds.")
end



function CrystalNet(crystal::Crystal)
    c = coalesce_sbus(crystal)
    CrystalNet(c.cell, c.types, c.graph, c.options)
end

function CrystalNetGroup(crystal::Crystal)
    c = coalesce_sbus(crystal)
    CrystalNetGroup(c.cell, c.types, c.graph, c.options)
end


find_clusters(c::Crystal) = find_clusters(c, c.options.clustering_mode)
function find_clusters(c::Crystal{T}, mode::_ClusteringMode)::Tuple{Clusters, _ClusteringMode} where T
    if mode == ClusteringMode.EachVertex
        return Clusters(length(c.types)), ClusteringMode.EachVertex
    elseif mode == ClusteringMode.Auto
        if T === Clusters
            return c.clusters, ClusteringMode.Auto
        else
            return find_clusters(c, ClusteringMode.EachVertex)
        end
    elseif mode == ClusteringMode.Input
        if T === Nothing
            throw(ArgumentError("Cannot use input residues as clusters: the input does not have residues."))
        else
            return c.clusters, ClusteringMode.Input
        end
    elseif mode == ClusteringMode.MOF
        clusters = find_sbus(c)
        if length(clusters.sbus) <= 1
            return find_clusters(c, ClusteringMode.MOFWiderOrganicSBUs)
        end
        return clusters, ClusteringMode.MOF
    elseif mode == ClusteringMode.MOFWiderOrganicSBUs
        sbus2 = SBUKinds([
            [:metal, :actinide, :lanthanide, :metalloid], [:C, :halogen, :nonmetal],
        ])
        clusters = find_sbus(c, sbus2)
        if length(clusters.sbus) <= 1
            return find_clusters(c, ClusteringMode.MOFMetalloidIsMetal)
        end
        return clusters, ClusteringMode.MOF
    elseif mode == ClusteringMode.MOFMetalloidIsMetal
        sbus3 = SBUKinds([
            [:metal, :actinide, :lanthanide, :metalloid], [:C, :halogen], [:nonmetal]
        ], Set{Int}(3))
        clusters = find_sbus(c, sbus3)
        if length(clusters.sbus) <= 1
            throw(InvalidSBU("ClusteringMode.MOF leads to a single cluster, choose a different clustering mode."))
        end
        return clusters, ClusteringMode.MOF
    elseif mode == ClusteringMode.Guess
        crystal = Crystal{Nothing}(c)
        uniquetypes = unique!(sort(c.types))
        if (:C ∈ uniquetypes) && any(x -> ismetal[atomic_numbers[x]], uniquetypes)
            try
                return find_clusters(crystal, ClusteringMode.MOF)
            catch e
                if !(e isa ClusteringError)
                    rethrow()
                end
            end
        end
        return find_clusters(c, ClusteringMode.Auto)
    end
end

function CrystalNet(g::Union{PeriodicGraph,AbstractString,AbstractVector{PeriodicEdge{D}} where D},
                    options::Options)
    graph = PeriodicGraph(g)
    cell = Cell()
    n = nv(graph)
    types = [Symbol("") for _ in 1:n]
    return CrystalNet(cell, types, graph, options)
end
function CrystalNet(g::Union{PeriodicGraph,AbstractString,AbstractVector{PeriodicEdge{D}} where D};
                    kwargs...)
    return CrystalNet(g, Options(; kwargs...))
end
