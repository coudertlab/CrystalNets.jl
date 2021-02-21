## Functions to handle input files.

using Chemfiles
using Statistics: std

include("clustering.jl")


"""
    parse_cif(file_path)

Parse a CIF file and return a dictionary where each identifier (without the
starting '_' character) is linked to its value.
Values are either a string or a vector of string (if defined in a loop).
"""
function parse_cif(file)
    all_data = Dict{String, Union{String, Vector{String}}}()
    inloop = false
    loopisspecified = false
    loopspec = String[]
    loop_n = 0

    l = read(file, String)
    i, j, x = nextword(l, 0)
    while i != 0

        if inloop
            if !loopisspecified
                if l[i] != '_' # This indicates the start of the values
                    loop_n = length(loopspec)
                    loopisspecified = true
                else # The identifier is part of the loop specification
                    push!(loopspec, l[i+1:j])
                    all_data[l[i+1:j]] = String[]
                    i, j, x = nextword(l, x); continue
                end
            end

            # From this point, the loop has been specified
            @assert loopisspecified
            if l[i] != '_' && l[i:j] != "loop_"
                for k in 1:loop_n
                    push!(all_data[loopspec[k]], l[i:j])
                    i, j, x = nextword(l, x)
                end
                continue
            end
            if l[i:j] == "loop_"
                loopisspecified = false
                loopspec = String[]
                i, j, x = nextword(l, x); continue
            end

            # This point can only be reached if we just quitted a loop
            inloop = false
        end

        @assert !inloop
        if l[i] == '_' # Simple identifier definition
            next_i, next_j, next_x = nextword(l, x)
            @assert next_i != 0
            all_data[l[i+1:j]] = l[next_i:next_j]
            x = next_x
        else
            if l[i:j] == "loop_"
                inloop = true
                loopisspecified = false
                loopspec = String[]
            elseif j-i > 4 && l[i:i+4] == "data_"
                @assert !haskey(all_data, "data")
                all_data["data"] = l[i+5:j]
            else
                k::Int = findprev(isequal('\n'), l, i)
                n::Int = count("\n", l[1:k])
                error("Unkown word \"$(l[i:j])\" at line $(n+1), position $(i-k):$(j-k)")
            end
        end

        i, j, x = nextword(l, x)
    end

    return all_data
end

function parsestrip(s)
    s = s[end] == ')' ? s[1:prevind(s, findlast('(', s))] : s
    return parse(BigFloat, s)
end


"""
    CIF(file_path::AbstractString)

Make a CIF object out of the parsed file.
"""
CIF(file::AbstractString) = CIF(parse_cif(file))
function CIF(parsed::Dict{String, Union{Vector{String},String}})
    natoms = length(parsed["atom_site_label"])
    equivalentpositions = pop!(parsed,
        haskey(parsed, "symmetry_equiv_pos_as_xyz") ?
            "symmetry_equiv_pos_as_xyz" : "space_group_symop_operation_xyz")
    refid = find_refid(equivalentpositions)
    cell = Cell(Symbol(haskey(parsed, "symmetry_cell_setting") ?
                            pop!(parsed, "symmetry_cell_setting") : ""),
                pop!(parsed, "symmetry_space_group_name_H-M"),
                haskey(parsed, "symmetry_Int_Tables_number") ?
                    parse(Int, pop!(parsed, "symmetry_Int_Tables_number")) : 0,
                parsestrip(pop!(parsed, "cell_length_a")),
                parsestrip(pop!(parsed, "cell_length_b")),
                parsestrip(pop!(parsed, "cell_length_c")),
                parsestrip(pop!(parsed, "cell_angle_alpha")),
                parsestrip(pop!(parsed, "cell_angle_beta")),
                parsestrip(pop!(parsed, "cell_angle_gamma")),
                parse.(EquivalentPosition, equivalentpositions, Ref(refid)))

    haskey(parsed, "symmetry_equiv_pos_site_id") && pop!(parsed, "symmetry_equiv_pos_site_id")
    removed_identity = false
    for i in eachindex(cell.equivalents)
        eq = cell.equivalents[i]
        if isone(eq.mat) && iszero(eq.ofs)
            deleteat!(cell.equivalents, i)
            removed_identity = true
            break
        end
    end
    @assert removed_identity

    labels =  pop!(parsed, "atom_site_label")
    _symbols = haskey(parsed, "atom_site_type_symbol") ?
                    pop!(parsed, "atom_site_type_symbol") : copy(labels)
    symbols = String[]
    @inbounds for x in _symbols
        i = findfirst(!isletter, x)
        push!(symbols, isnothing(i) ? x : x[1:i-1])
    end
    pos_x = pop!(parsed, "atom_site_fract_x")
    pos_y = pop!(parsed, "atom_site_fract_y")
    pos_z = pop!(parsed, "atom_site_fract_z")


    types = Symbol[]
    pos = Matrix{Float64}(undef, 3, natoms)
    correspondence = Dict{String, Int}()
    for i in 1:natoms
        @assert !haskey(correspondence, labels[i])
        correspondence[labels[i]] = i
        push!(types, Symbol(symbols[i]))
        pos[:,i] = parsestrip.([pos_x[i], pos_y[i], pos_z[i]])
        pos[:,i] .-= floor.(Int, pos[:,i])
    end

    invids = sortperm(types)
    types = types[invids]
    ids = invperm(invids)
    bonds = falses(natoms, natoms)
    if haskey(parsed, "geom_bond_atom_site_label_1") &&
       haskey(parsed, "geom_bond_atom_site_label_2")
        bond_a = pop!(parsed, "geom_bond_atom_site_label_1")
        bond_b = pop!(parsed, "geom_bond_atom_site_label_2")
        for i in 1:length(bond_a)
            try
                x = correspondence[bond_a[i]]
                y = correspondence[bond_b[i]]
                bonds[x,y] = bonds[y,x] = 1
            catch e
                if e isa KeyError
                    @ifwarn @warn "Atom $(e.key) will be ignored since it has no placement in the CIF file."
                else
                    rethrow()
                end
            end
        end
    end

    return CIF(parsed, cell, ids, types, pos, bonds)
end



"""
    parse_arc(file)

Parse a .arc Systre archive such as the one used by the RCSR.
Return a pair `(flag, pairs)`.

`flag` is set if the archive corresponds to one generated by a compatible release
of CrystalNets. If unset, the genomes of the archive may not be the same as those
computed by CrystalNets for the same nets.
`pairs` is a `Dict{String,String}` whose entries have the form `genome => id` where `id`
is the name of the net and `genome` is the topological genome corresponding to this net
(given as a string of whitespace-separated values parseable by `PeriodicGraph`).
"""
function parse_arc(file)
    pairs = Dict{String,String}()
    curr_key = ""
    counter = 1
    firstline = readline(file)
    flag = startswith(firstline, "Made by CrystalNets.jl v")
    if flag
        flag = CRYSTAL_NETS_VERSION == VersionNumber(firstline[25:end])
    end
    for l in eachline(file)
        if length(l) > 3 && l[1:3] == "key"
            @assert isempty(curr_key)
            i = 4
            while isspace(l[i])
                i += 1
            end
            curr_key = l[i:end]
            @assert !haskey(pairs, curr_key)
        elseif length(l) > 2 && l[1:2] == "id"
            @assert !isempty(curr_key)
            i = 3
            while isspace(l[i])
                i += 1
            end
            pairs[curr_key] = l[i:end]
            curr_key = ""
        end
    end
    return flag, pairs
end

"""
    parse_arcs(file)

Parse a folder containing .arc Systre archives such as the one used by the RCSR.
Return a pair `(flag, pairs)` with the same convention than `parse_arc`
"""
function parse_arcs(path)
    combine(x, y) = x * ", " * y
    dict = Dict{String,String}()
    flag = true
    for f in readdir(path; sort=true)
        _flag, _dict = parse_arc(path*f)
        flag &= _flag
        mergewith!(combine, dict, _dict)
    end
    return flag, dict
end

function parse_atom_name(name::AbstractString)
    firstsep = findfirst(x -> ispunct(x) || isspace(x) || isnumeric(x), name)
    if firstsep isa Nothing
        return String(name)
    else
        firstsep::Int
        if firstsep != 1 && !any(isuppercase, @view name[nextind(name, firstsep):end])
            return String(name[1:prevind(name, firstsep)])
        else
            return String(name)
        end
    end
end

function parse_atom(name)
    atom = Chemfiles.Atom(name)
    set_type!(atom, parse_atom_name(name))
    return atom
end


function set_unique_bond_type!(frame::Frame, types, bond_length, bonded_atoms::Tuple{Symbol, Symbol}, onlykeep, tol)
    @ifwarn @info "To avoid guessing bonds, use a file format that contains the bonds."
    n = Int(size(frame))
    pos = Chemfiles.positions(frame)
    mat = Chemfiles.matrix(UnitCell(frame))'
    invmat = inv(mat)
    indices = [i for i in 1:n if types[i] ∈ onlykeep]
    invpos = [invmat*pos[:,i] for i in indices]
    #=@inbounds=# for _i in 1:length(indices)
        i = indices[_i]
        for _j in _i+1:length(indices)
            j = indices[_j]
            if minmax(types[i], types[j]) == bonded_atoms
                bonded = abs2(periodic_distance(invpos[_i], invpos[_j], mat) - bond_length) <= tol
                if bonded
                    Chemfiles.add_bond!(frame, i-1, j-1)
                end
            end
        end
    end
    nothing
end

function try_guess_bonds!(frame::Frame, types)
    unique_types = unique!(sort(types))
    if unique_types == [:C] || unique_types == [:C, :H]
        @ifwarn @warn "Guessing bonds. The structure seems to be made of only carbons (and possibly hydrogens): using an interatomic distance of 1.54±0.3 Å to assign edges between C atoms."
        set_unique_bond_type!(frame, types, 1.54, (:C, :C), (:C,), 0.3)
    elseif unique_types == [:O, :Si] || unique_types == [:Si]
        @ifwarn @warn "Guessing bonds. The structure seems to be a zeolite: using an interatomic distance of 3.1±0.3 Å to assign edges between Si atoms."
        set_unique_bond_type!(frame, types, 1.63, (:O, :Si), (:O, :Si), 0.15)
    else
        @ifwarn begin
            @warn "Guessing bonds through Chemfiles. This may take a while for big structures and may be inexact."
            @info "To avoid guessing bonds, use a file format that contains the bonds."
        end
        guess_bonds!(frame)
    end
end


function attribute_residues(topology, assert_use_existing_residues)
    m = Int(count_residues(topology))
    n = Int(size(topology))
    atoms_in_residues = m == 0 ? 0 : sum(length(atoms(Residue(topology, i))) for i in 0:(m-1))
    @assert atoms_in_residues <= n

    if atoms_in_residues < n
        if assert_use_existing_residues
            throw(ArgumentError("""
            Cannot use existing residues as vertices because not all atoms have an associated residue.
            To fix this, either assign a residue to each atom or provide another way to detect the vertices.
            """))
        end
        @ifwarn begin
            if atoms_in_residues > 0
                @warn "Some but not all atoms have an associated residue, so we cannot rely on existing residues"
            end
        end
        attributions = Int[]
    else
        attributions = zeros(Int, n)
        for i_residue in 1:m
            for atom in atoms(Residue(topology, i_residue-1))
                attributions[atom+1] = i_residue
            end
        end
    end
    return attributions
end

@static if !isdefined(Chemfiles, :atoms) # up to 0.9.3 included
    function atoms(residue::Residue)
        count = size(residue)
        result = Array{UInt64}(undef, count)
        Chemfiles.__check(Chemfiles.lib.chfl_residue_atoms(Chemfiles.__const_ptr(residue), pointer(result), count))
        return result
    end
end


function check_collision(pos, mat)
    _n, n = size(pos)
    @assert _n == 3
    invmat = inv(mat)
    toremove = Int[]
    for i in 1:n, j in (i+1):n
        if periodic_distance(invmat*pos[:,i], invmat*pos[:,j], mat) < 0.55
            push!(toremove, i)
            push!(toremove, j)
        end
    end
    tokeep = collect(1:n)
    if !isempty(toremove)
        @ifwarn @warn "This file contains multiple colliding atoms. All colliding atoms will be removed."
        unique!(sort!(toremove))
    end
    return toremove
end


# macro checkvalence(type, degree, mode)
#     msg = Expr(:call, :string, "In ", mode, " mode, found ", type,
#                " with incorrect valence (", :d, ").")
#     return esc(quote
#         if t === $type
#             d == $degree && continue
#             CrystalNets.@ifwarn @warn $msg
#             return true
#         end
#     end)
# end
# @checkvalence :O 2 specialmode
# @checkvalence :Si 4 specialmode
# if specialmode === :zeolite
#     throw(ArgumentError("Found atom $t in zeolite mode, where only Si and O are accepted. Choose a different mode."))
# end

"""
    least_plausible_neighbours(Δs, n, δ)

Find the positions of the n least probable neighbours of an atom, given the list
Δs of the distance between their position and that of the atom.

This function is highly empirical and should not be considered utterly reliable.
"""
function least_plausible_neighbours(Δs, n)
    m = length(Δs)
    p::Vector{Int} = sortperm(Δs)
    n == 1 && return [p[1]]
    cycle = Δs[p[1:n]]
    avg = mean(cycle)
    means = Float64[avg]
    stds = Float64[std(cycle; mean=avg)]
    for i in n+1:m
        popfirst!(cycle)
        push!(cycle, Δs[p[i]])
        avg += (-Δs[p[i-n]] + Δs[p[i]])/n
        push!(means, avg)
        push!(stds, std(cycle; mean=avg))
    end
    if m == n+1
        _min, _max = minmax(stds[1], stds[2])
        if _max - _min > 0.5/n
            return stds[1] < stds[2] ? [p[1]] : [p[2]]
        end
        return [p[1]]
    else
        j = argmin(collect(stds[i] + abs(means[i] - 1.65)/2 + means[i]/3 for i in 1:(m-n)))
        ret = p[1:j-1]
        append!(ret, @view p[n+j:end])
        return ret
    end
end

macro reduce_valence(n)
    return esc(quote
        neighs = neighbors(graph, i)
        m = length(neighs)
        m == $n && continue
        (!dofix || m < $n) && push!(invalidatoms, t)
        if dofix && m > $n
            posi = pos[i]
            Δs = [norm(pos[x.v] .+ mat * x.ofs - posi) for x in neighs]
            toremove = least_plausible_neighbours(Δs, $n)
            neighs = copy(neighs) # otherwise the list is modified by rem_edge!
            for v in toremove
                rem_edge!(graph, PeriodicEdge{N}(i, neighs[v]))
            end
        end
    end)
end

macro reduce_valence(n1, n2)
    return esc(quote
        neighs = neighbors(graph, i)
        m = length(neighs)
        $n1 ≤ m ≤ $n2 && continue
        (!dofix || m < $n1) && push!(invalidatoms, t)
        if dofix && m > $n2
            posi = pos[i]
            Δs = [norm(pos[x.v] .+ mat * x.ofs - posi) for x in neighs]
            toremove = least_plausible_neighbours(Δs, $n2)
            neighs = copy(neighs) # otherwise the list is modified by rem_edge!
            for v in toremove
                rem_edge!(graph, PeriodicEdge{N}(i, neighs[v]))
            end
        end
    end)
end
#
# macro loop_reduce_valence(typs...)
#     if length(typs) == 1
#         return esc(quote
#             if $(typs[1].args[2]) ∈ uniquetypes
#                 for i in 1:n
#                     t = types[i]
#                     if t === $(typs[1].args[2])
#                         @reduce_valence $(typs[1].args[3])
#                     end
#                 end
#             end
#         end)
#     else
#         x1 = pop!(typs)
#         x2 = pop!(typs)
#         cond = Expr(:||, Expr(:call, :∈, QuoteNode(x1.args[2]), :uniquetypes),
#                          Expr(:call, :∈, QuoteNode(x2.args[2]), :uniquetypes))
#         subcond = Expr(:||, Expr(:call, :===, :t, QuoteNode(x1.args[2])),
#                             Expr(:call, :===, :t, QuoteNode(x2.args[2])))
#         for x in typs
#             cond = Expr(:||, Expr(:call, :∈, QuoteNode(x.args[2]), :uniquetypes), cond)
#             subcond = Expr(:||, Expr(:call, :===, :t, QuoteNode(x.args[2])))
#         end
#     return quote
#         if $cond
#             for i in 1:n
#                 t = types[i]
#                 if $subcond
#                     @reduce_valence
#     end
# end

#hasHneighbor(types, graph, i) = any(x -> types[x.v] === :H, neighbors(graph,i))
function fix_valence!(graph::PeriodicGraph{N}, pos, types, mat, ::Val{dofix}) where {N,dofix}
    n = length(types)
    ## Small atoms valence check
    invalidatoms = Set{Symbol}()
    uniquetypes = unique!(sort(types))
    if :H ∈ uniquetypes || :Na ∈ uniquetypes
        for i in 1:n
            t = types[i]
            if t === :H || t === :Na
                @reduce_valence 1
            end
        end
    end
    if :O ∈ uniquetypes
        for i in 1:n
            t = types[i]
            if t === :O
                @reduce_valence 2
            end
        end
    end
    if :N ∈ uniquetypes || :C ∈ uniquetypes
        for i in 1:n
            t = types[i]
            if t === :N
                @reduce_valence 2 4
            elseif t === :C # sp1 carbon is not a common occurence
                @reduce_valence 3 4
            end
        end
    end
    if !isempty(invalidatoms)
        s = String.(collect(invalidatoms))
        @ifwarn @warn (dofix ? "After fix, f" : "F")*"ound $(join(s, ',', " and ")) with invalid number of bonds."
        return true
    end
    return false
end


"""
    @enum BondingMode

Selection mode for the detection of bonds. The choices are:
-   `InputBonds`: use the input bonds. Fail if those are not specified.
-   `ChemfilesBonds`: use chemfiles built-in bond detection mechanism.
-   `AutoBonds`: if the input specifies bonds, use them unless they look suspicious (too or too
    large according to a heuristic). Otherwise, fall back to `ChemfilesBonds`.
"""
@enum BondingMode begin
    InputBonds
    ChemfilesBonds
    AutoBonds
end


function sanity_checks!(graph, pos, types, mat, bondingmode)
    ## Bond length check
    removeedges = PeriodicEdge3D[]
    for e in edges(graph)
        s, d = e.src, e.dst.v
        bondlength = norm(pos[d] .+ (mat * e.dst.ofs) .- pos[s])
        if (bondlength < 0.85 && types[s] !== :H && types[d] !== :H)
            push!(removeedges, e)
        elseif bondlength > 3
            @ifwarn @warn "Suspiciously large bond found: $bondlength pm between $(types[s]) and $(types[d])."
            return true
        end
    end
    if bondingmode == AutoBonds
        @ifwarn begin
            if !isempty(removeedges)
                @warn "Suspicious bonds lengths found. Such bonds are probably spurious and will be deleted."
                @info "To force retaining these bonds, use --bond-detect=input or --bond-detect=chemfiles"
            end
        end
        for e in removeedges
            rem_edge!(graph, e)
        end
    end
    return false
end


"""
       parse_chemfiles(path)

Parse a file given in any reckognised chemical format and extract the topological
information.
Such format can be .cif or any file format reckognised by Chemfiles.jl that
contains all the necessary topological information.
"""
function parse_chemfile(_path, exportto=tempdir(), bondingmode::BondingMode=AutoBonds, assert_use_existing_residues=false; ignore_atoms=[])
    # Separate the cases unhandled by Chemfiles from the others
    path = expanduser(_path)
    cell = Cell()
    guessed_bonds = false
    frame::Frame = if lowercase(last(splitext(path))) == ".cif"
        global framecif = Frame()
        cif = expand_symmetry(CIF(path))
        a, b, c, α, β, γ = cell_parameters(cif.cell)
        set_cell!(framecif, UnitCell(a, b, c, α, β, γ))
        n = length(cif.ids)
        ignored = 0
        types = Symbol[]
        for i in 1:n
            typ = cif.types[cif.ids[i]]
            if typ ∈ ignore_atoms
                ignored += 1
                continue
            end
            push!(types, Symbol(parse_atom_name(String(typ))))
            pos = cif.cell.mat * cif.pos[:,i]
            atom = parse_atom(string(typ))
            add_atom!(framecif, atom, Vector{Float64}(pos))
        end
        n -= ignored
        if iszero(cif.bonds)
            guessed_bonds = true
            try_guess_bonds!(framecif, types)
        else
            for i in 1:n, j in (i+1):n
                if cif.bonds[i,j]
                    add_bond!(framecif, i-1, j-1)
                end
            end
        end
        attributions = Int[]


        framecif
    else # The rest is handled by Chemfiles

        frameelse = read(Trajectory(path))
        for i in Int(size(frameelse))-1:-1:0
            if Symbol(type(Chemfiles.Atom(frameelse, i))) ∈ ignore_atoms
                Chemfiles.remove_atom!(frameelse, i)
            end
        end

        toremove = reverse(check_collision(positions(frameelse), matrix(UnitCell(frame))))
        for j in toremove
            Chemfiles.remove_atom!(frameelse, j-1)
        end

        topology = Topology(frameelse)
        n = Int(size(topology))
        types = [Symbol(type(Chemfiles.Atom(frameelse, i))) for i in 0:(n-1)]

        if bonds_count(topology) == 0
            guessed_bonds = true
            try_guess_bonds!(frameelse, types)
            topology = Topology(frameelse) # safer but useless since the underlying pointer is the same
        end

        attributions = attribute_residues(topology, assert_use_existing_residues)

        frameelse
    end

    cell = Cell(Cell(), SMatrix{3,3,BigFloat}(matrix(UnitCell(frame)))')
    if !all(isfinite, cell.mat) || iszero(det(cell.mat))
        @ifwarn @warn "Suspicious unit cell of matrix $(Float64.(cell.mat)). Is the input really periodic? Using a cubic unit cell instead"
        cell = Cell()
    end

    poss = copy(positions(frame))
    colpos::Vector{SVector{3,Float64}} = collect(eachcol(poss))

    adjacency = zeros(Bool, n, n)
    topology = Topology(frame)
    for (a,b) in eachcol(bonds(topology))
        adjacency[a+1,b+1] = true
        adjacency[b+1,a+1] = true
    end
    mat = Float64.(cell.mat)
    graph = PeriodicGraph3D(n, edges_from_bonds(adjacency, mat, colpos))

    if bondingmode != InputBonds
        bad_valence = fix_valence!(graph, colpos, types, mat, Val(false))
        recompute_bonds = bad_valence || sanity_checks!(graph, colpos, types, mat, bondingmode)
        if recompute_bonds
            if !guessed_bonds
                @ifwarn begin
                    @warn "Disregarding all bonds from the input file."
                    @info "To force retaining the initial bonds, use --bond-detect=input"
                end
                try_guess_bonds!(frame, types)
                topology = Topology(frame)
                adjacency = zeros(Bool, n, n)
                for (a,b) in eachcol(bonds(topology))
                    adjacency[a+1,b+1] = true
                    adjacency[b+1,a+1] = true
                end
                graph = PeriodicGraph3D(n, edges_from_bonds(adjacency, cell.mat, colpos))
            end
            remaining_not_fixed = fix_valence!(graph, colpos, types, mat, Val(true))
            @ifwarn if remaining_not_fixed
                @warn "Remaining atoms with invalid valence. Proceeding anyway."
            end
            recompute_bonds = sanity_checks!(graph, colpos, types, cell.mat, bondingmode)
            @ifwarn if recompute_bonds
                 @warn "Remaining bonds of suspicious lengths. Proceeding anyway."
            end
            attributions = attribute_residues(topology, assert_use_existing_residues)
        end
    end

    if isempty(attributions)
        crystalnothing = Crystal{Nothing}(cell, types, nothing, poss, graph)
        ifexport(crystalnothing, splitext(splitdir(path)[2])[1], exportto)
        return crystalnothing
    else
        crystalclusters = Crystal{Clusters}(cell, types, regroup_sbus(graph, attributions), poss, graph)
        ifexport(crystalclusters, splitext(splitdir(path)[2])[1], exportto)
        return crystalclusters
    end
end
