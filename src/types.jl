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
    mat::SMatrix{3,3,Int,9}
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
        refid = [""]
        not_at_the_end = true
        immediate_continue = false
        for t in ts
            @toggleassert not_at_the_end
            k = Tokenize.Tokens.kind(t)
            if k === Tokenize.Tokens.IDENTIFIER
                @toggleassert refid[end] == ""
                refid[end] = Tokenize.Tokens.untokenize(t)
            elseif k === Tokenize.Tokens.COMMA || k === Tokenize.Tokens.SEMICOLON
                @toggleassert refid[end] != ""
                push!(refid, "")
            elseif k === Tokenize.Tokens.ENDMARKER
                not_at_the_end = false
            elseif k !== Tokenize.Tokens.WHITESPACE && t.kind !== Tokenize.Tokens.PLUS
                immediate_continue = true
                break
            end
        end
        immediate_continue && continue
        if not_at_the_end
            error("Unknown end of line marker for symmetry equivalent {$eq}")
        end
        if length(refid) != 3 || refid[end] == ""
            error("Input string {$eq} is not a valid symmetry equivalent")
        end
        return tuple(lowercase(refid[1]), lowercase(refid[2]), lowercase(refid[3]))
    end
    return ("x", "y", "z")
end

function Base.parse(::Type{EquivalentPosition}, s::AbstractString, refid=("x", "y", "z"))
    const_dict = Dict{String, Int}(refid[1]=>1, refid[2]=>2, refid[3]=>3)
    mat = zeros(Int, 3, 3)
    ofs = zeros(Rational{Int}, 3)
    curr_num::Union{Int, Nothing} = nothing
    curr_val::Union{Rational{Int}, Nothing} = nothing
    curr_sign::Union{Bool, Nothing} = nothing
    encountered_div::Bool = false
    i = 1
    something_written = false
    for x in tokenize(lowercase(s))
        k = Tokenize.Tokens.kind(x)
        k === Tokenize.Tokens.WHITESPACE && continue
        if k === Tokenize.Tokens.INTEGER
            @toggleassert isnothing(curr_val)
            if encountered_div
                curr_val = Rational{Int}(Int(curr_num), parse(Int, x.val))
                curr_num = nothing
                encountered_div = false
            else
                @toggleassert isnothing(curr_num)
                curr_num = parse(Int, x.val)
            end
        else
            @toggleassert !encountered_div
            if k === Tokenize.Tokens.IDENTIFIER
                if !isnothing(curr_num)
                    @toggleassert isnothing(curr_val)
                    curr_val = curr_num
                    curr_num = nothing
                end
                sign = isnothing(curr_sign) ? 1 : 2*curr_sign - 1
                val = isnothing(curr_val)  ? 1 : curr_val
                j = const_dict[Tokenize.Tokens.untokenize(x)]
                mat[i,j] += sign * val
                curr_val = nothing
                curr_sign = nothing
                something_written = true
            else
                if x.kind === Tokenize.Tokens.FWD_SLASH
                    @toggleassert isnothing(curr_val)
                    @toggleassert !isnothing(curr_num)
                    encountered_div = true
                    continue
                end
                if x.kind === Tokenize.Tokens.FLOAT
                    @toggleassert isnothing(curr_val)
                    @toggleassert isnothing(curr_num)
                    @toggleassert !encountered_div
                    curr_val = Rational{Int}(rationalize(Int8, parse(Float16, x.val)))
                    continue
                end
                if !isnothing(curr_num)
                    @toggleassert isnothing(curr_val)
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
                    @toggleassert isnothing(curr_sign)
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
                    i != 3 && error("Input string \"$s\" is not a valid symmetry equivalent")
                end
            end
        end
    end
    EquivalentPosition(SMatrix{3,3,Int,9}(mat), SVector{3,Rational{Int}}(ofs))
end

function Base.show(io::IO, eq::EquivalentPosition)
    function rationaltostring(x, notofs::Bool, first::Bool)
        if notofs && (x == 1 || x == -1)
            return x < 0 ? '-' : first ? "" : "+"
        end
        sign = x < 0 || first ? "" : "+"
        sign * (denominator(x) == 1 ? string(numerator(x)) : string(numerator(x))*'/'*string(denominator(x)))
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

See [`SPACE_GROUP_HALL`](@ref), [`SPACE_GROUP_FULL`](@ref), [`SPACE_GROUP_HM`](@ref)`
and [`SPACE_GROUP_IT`](@ref) for the correspondance between Hall number and usual symbolic
representations.
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
        id = popfirst!(eqs)
        @toggleassert isone(id.mat) && iszero(id.ofs)
    end
    return Cell(hall, mat, eqs)
end
Cell() = Cell(1, (10, 10, 10), (90, 90, 90), EquivalentPosition[])

function cell_parameters(mat::AbstractMatrix)
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
    a, b, c, α, β, γ = Float64.(cell_parameters(cell))
    hall_symbol, crystal_system = HALL_SYMBOLS[cell.hall]
    print(io, Cell, " with Hall symbol $hall_symbol ($crystal_system) and parameters ")
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
    bonds::Vector{Vector{Tuple{Int,Float32}}}
end

function keepinbonds(bonds, keep)
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

function add_to_bondlist!(bondlist, x, d)
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

function get_bondlist(bondlist, x)
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

function sortprune_bondlist!(bondlist)
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

function keep_atoms(cif::CIF, kept)
    kept_ids = sort!([cif.ids[i] for i in kept])
    unique!(kept_ids)
    idmap = Vector{Int}(undef, length(cif.types)) # upper bound on maximum(kept_ids)
    for (i,x) in enumerate(kept_ids)
        idmap[x] = i
    end
    return CIF(cif.cifinfo, cif.cell, [idmap[cif.ids[i]] for i in kept],
               cif.types[kept_ids], cif.pos[:, kept], keepinbonds(cif.bonds, kept))
end


function prepare_periodic_distance_computations(mat)
    a, b, c, α, β, γ = cell_parameters(mat)
    ortho = all(x -> isapprox(Float16(x), 90; rtol=0.02), (α, β, γ))
    _a, _b, _c = eachcol(mat)
    safemin = min(dot(cross(_b, _c), _a)/(b*c),
                  dot(cross(_c, _a), _b)/(a*c),
                  dot(cross(_a, _b), _c)/(a*b))/2
    # safemin is the half-distance between opposite planes of the unit cell
    return MVector{3,Float64}(undef), ortho, safemin
end

function periodic_distance!(buffer, u, mat, ortho, safemin)
    @simd for i in 1:3
        diff = u[i] + 0.5
        buffer[i] = diff - floor(diff) - 0.5
    end
    ref = norm(mat*buffer)
    (ortho || ref ≤ safemin) && return ref
    @inbounds for i in 1:3
        buffer[i] += 1
        newnorm = norm(mat*buffer)
        newnorm < ref && return newnorm # in a reduced lattice, there should be at most one
        buffer[i] -= 2
        newnorm = norm(mat*buffer)
        newnorm < ref && return newnorm
        buffer[i] += 1
    end
    return ref
end

"""
    periodic_distance(u, mat, ortho=nothing, safemin=nothing)

Distance between point `u` and the origin, given as a triplet of fractional coordinates, in
a repeating unit cell of matrix `mat`.
The distance is the shortest between all equivalents of `u` and the origin.
If `ortho` is set to `true`, the angles α, β and γ of the cell are assumed right, which
accelerates the computation by up to 7 times.
If a distance lower than `safemin` is computed, stop trying to find a periodic image of `u`
closer to the origin.
If unspecified, both `ortho` and `safemin` are automatically determined from `mat`.

This implementation assumes that the cell corresponds to a reduced lattice. It may be
invalid for some edge cases otherwise.

For optimal performance, use `periodic_distance!` with `buffer`, `ortho` and `safemin`
obtained from `prepare_periodic_distance_computations`.
"""
function periodic_distance(u, mat, ortho=nothing, safemin=nothing)
    if ortho === nothing || safemin === nothing
        _, ortho, safemin = prepare_periodic_distance_computations(mat)
    end
    periodic_distance!(similar(u), u, mat, ortho, safemin)
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
        if periodic_distance!(buffer, points[i] .- points[j], smallmat, ortho, safemin) < 0.55
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
                if periodic_distance!(buffer, newpos[j] .- p, smallmat, ortho, safemin) < 0.55
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
                if abs(periodic_distance!(buffer, newpos[i] .- newpos[j], smallmat, ortho, safemin) - bondlength) < 0.55
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
    @toggleassert size(placement) == (D, n)
    pos = Vector{SVector{D,T}}(undef, n)
    offsets = Vector{SVector{D,Int}}(undef, n)
    @inbounds for (i, x) in enumerate(eachcol(placement))
        offsets[i] = floor.(Int, x)
        pos[i] = Base.unsafe_rational.(getfield.(x, :num) .- getfield.(x, :den).*offsets[i], getfield.(x, :den))
    end
    s = sortperm(pos)
    pos = pos[s]
    types = Symbol[types[s[i]] for i in 1:n]
    graph = offset_representatives!(graph, .-offsets)[s]
    # @toggleassert all(pos[i] == mean(pos[x.v] .+ x.ofs for x in neighbors(graph, i)) for i in 1:length(pos))
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


function CrystalNetGroup(cell::Cell, types::AbstractVector{Symbol},
                         graph::PeriodicGraph, opts::Options)
    initialvmap, graph = trim_topology(PeriodicGraphs.change_dimension(PeriodicGraph3D, graph))
    types = types[initialvmap]
    remove_homoatomic_bonds!(graph, types, opts.ignore_homoatomic_bonds)
    if !isempty(opts.export_net) && !isempty(opts._pos)
        pos = opts._pos[initialvmap]
        export_default(Crystal{Nothing}(cell, types, pos, graph, opts), 
                       "net", opts.name, opts.export_net)
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
                    graph::PeriodicGraph, options::Options)::CrystalNet
    group = CrystalNetGroup(cell, types, graph, options)
    D = isempty(group.D3) ? isempty(group.D2) ? isempty(group.D1) ? 0 : 1 : 2 : 3
    D == 0 && return CrystalNet{0,Int8}(cell, Symbol[], SVector{0,Int8}[], PeriodicGraph{0}(), options)
    nonuniqueD = D != 1 && (!isempty(group.D1) || (D == 3 && !isempty(group.D2)))
    if nonuniqueD
        @ifwarn @warn "Presence of periodic structures of different dimensionalities. Only the highest dimensionality ($D here) will be retained."
    end
    interpenetratingD = false
    @separategroups D group begin
        interpenetratingD = length(group) > 1
    end
    if interpenetratingD
        error(ArgumentError("Multiple interpenetrating $D-dimensional structures. Cannot handle this as a single CrystalNet, use CrystalNetGroup instead."))
    end
    @separategroups D group begin
        return last(first(group))
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
    isempty(placement) && return CrystalNet{D,Rational{Int}}(cell, types, graph, placement, options)
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


function CrystalNet(crystal::Crystal)
    c = coalesce_sbus(crystal)
    CrystalNet(c.cell, c.types, c.graph, c.options)
end

function CrystalNetGroup(crystal::Crystal)
    c = coalesce_sbus(crystal)
    CrystalNetGroup(c.cell, c.types, c.graph, c.options)
end


find_clusters(c::Crystal) = find_clusters(c, c.options.structure)
function find_clusters(c::Crystal{T}, mode::_StructureType)::Tuple{Clusters, _StructureType} where T
    if mode == StructureType.EachVertex || mode == StructureType.Zeolite
        return Clusters(length(c.types)), mode
    elseif mode == StructureType.Auto
        if T === Clusters
            return c.clusters, StructureType.Auto
        else
            return find_clusters(c, StructureType.EachVertex)
        end
    elseif mode == StructureType.Input
        if T === Nothing
            throw(ArgumentError("Cannot use input residues as clusters: the input does not have residues."))
        else
            return c.clusters, StructureType.Input
        end
    elseif mode == StructureType.MOF || mode == StructureType.Cluster
        clusters = find_sbus(c, c.options.cluster_kinds)
        if length(clusters.sbus) <= 1 && c.options.cluster_kinds == default_sbus
            return find_clusters(c, StructureType.MOFWiderOrganicSBUs)
        end
        return clusters, StructureType.MOF
    elseif mode == StructureType.MOFWiderOrganicSBUs
        sbus2 = ClusterKinds([
            [:metal, :actinide, :lanthanide, :metalloid], [:C, :halogen, :nonmetal],
        ])
        clusters = find_sbus(c, sbus2)
        if length(clusters.sbus) <= 1
            return find_clusters(c, StructureType.MOFMetalloidIsMetal)
        end
        return clusters, StructureType.MOF
    elseif mode == StructureType.MOFMetalloidIsMetal
        sbus3 = ClusterKinds([
            [:metal, :actinide, :lanthanide, :metalloid], [:C, :halogen], [:nonmetal]
        ], Int[3])
        clusters = find_sbus(c, sbus3)
        return clusters, StructureType.MOF
    elseif mode == StructureType.Guess
        crystal = Crystal{Nothing}(c)
        try
            return find_clusters(crystal, StructureType.Cluster)
        catch e
            if !(e isa ClusteringError)
                rethrow()
            end
        end
        return find_clusters(c, StructureType.Auto)
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
