## Functions to handle input files.

using Chemfiles
using Statistics: std

include("clustering.jl")
include("guessbonds.jl")


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
    bonds = fill(Inf32, natoms, natoms)
    if haskey(parsed, "geom_bond_atom_site_label_1") &&
       haskey(parsed, "geom_bond_atom_site_label_2")
        bond_a = pop!(parsed, "geom_bond_atom_site_label_1")
        bond_b = pop!(parsed, "geom_bond_atom_site_label_2")
        dists = haskey(parsed, "geom_bond_distance") ? 
                    parse.(Float32, pop!(parsed, "geom_bond_distance")) :
                    fill(zero(Float32), length(bond_a))
        for i in 1:length(bond_a)
            try
                x = correspondence[bond_a[i]]
                y = correspondence[bond_b[i]]
                bonds[x,y] = bonds[y,x] = dists[i]
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



function attribute_residues(residues, n, assert_use_existing_residues)
    m = length(residues)
    atoms_in_residues = m == 0 ? 0 : sum(length(atoms(r)) for r in residues)
    @assert atoms_in_residues <= n

    if atoms_in_residues < n
        if assert_use_existing_residues
            throw(ArgumentError("""
            Cannot use existing residues as vertices because not all atoms have an associated residue.
            To fix this, either assign a residue to each atom or provide another way to detect the vertices.
            """))
        end
        if atoms_in_residues > 0
            @ifwarn @warn "Some but not all atoms have an associated residue, so we cannot rely on existing residues"
        end
        return Int[]
    end
    attributions = zeros(Int, n)
    for i_residue in 1:m
        for atom in atoms(residues[i_residue])
            attributions[atom+1] = i_residue
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
    n = length(pos)
    toremove = Int[]
    for i in 1:n
        posi = pos[i]
        for j in (i+1):n
            if periodic_distance(posi, pos[j], mat) < 0.55
                push!(toremove, j)
            end
        end
    end
    if !isempty(toremove)
        @ifwarn @warn "This file contains multiple colliding atoms. Only one atom will be kept per site."
        sort!(toremove)
        unique!(toremove)
    end

    return toremove
end



"""
    least_plausible_neighbours(Δs, n)

Find the positions of the n least probable neighbours of an atom, given the list
Δs of the distance between their position and that of the atom.

This function is highly empirical and should not be considered utterly reliable.
"""
function least_plausible_neighbours(Δs, n)
    m = length(Δs)
    m ≤ n && return collect(1:m) # may happen because H bonds cannot be removed
    p::Vector{Int} = sortperm(Δs)
    # If we should only remove n neighbours, remove those which minimize
    # the dispersion among the remaining others
    m == n+1 && return p[2:end] # only one remaining neighbour: keep the closest
    if std(Δs[p[n+1:end]]) > 0.9std(Δs[p[1:end-n]])
        return p[end-n+1:end]
    end
    return p[1:n]
end

macro reduce_valence(n)
    return esc(quote
        neighs = neighbors(graph, i)
        m = length(neighs)
        m == $n && continue
        (!dofix || m < $n) && push!(invalidatoms, t)
        if dofix && m > $n
            posi = pos[i]
            Δs = [norm(mat * (pos[x.v] .+ x.ofs - posi)) for x in neighs]
            toremove = least_plausible_neighbours(Δs, m - $n)
            neighs = copy(neighs) # otherwise the list is modified by rem_edge!
            for v in toremove
                rem_edge!(graph, PeriodicEdge{N}(i, neighs[v]))
            end
        end
    end)
end

macro reduce_valence(n1, n2)
    comparison = n1 == 0 ? :(m ≤ $n2) : :($n1 ≤ m ≤ $n2)
    invalidcond = n1 == 0 ? :(!dofix) : :(!dofix || m < $n1)
    return esc(quote
        neighs = neighbors(graph, i)
        m = length(neighs)
        $comparison && continue
        ($invalidcond) && push!(invalidatoms, t)
        if dofix && m > $n2
            posi = pos[i]
            noHatoms = [x for x in neighs if types[x.v] !== :H]
            Δs = [norm(mat * (pos[x.v] .+ x.ofs - posi)) for x in noHatoms]
            toremove = least_plausible_neighbours(Δs, m - $n2)
            for v in toremove
                rem_edge!(graph, PeriodicEdge{N}(i, noHatoms[v]))
            end
        end
    end)
end


#hasHneighbor(types, graph, i) = any(x -> types[x.v] === :H, neighbors(graph,i))
function fix_valence!(graph::PeriodicGraph{N}, pos, types, mat, ::Val{dofix}) where {N,dofix}
    # Small atoms valence check
    n = length(types)
    invalidatoms = Set{Symbol}()
    # First pass over H, since those are likely bonded to their closest neighbor
    for i in 1:n
        t = types[i]
        if t === :H
            @reduce_valence 1
        end
    end
    monovalent = Set{Symbol}([:Li, :Na, :K, :F, :Br, :Cl, :I])
    for i in 1:n
        t = types[i]
        if t ∈ monovalent
            @reduce_valence 0 1
        elseif t === :O
            @reduce_valence 0 2
        elseif t === :N
            @reduce_valence 2 4
        elseif t === :C  # sp1 carbon is not a common occurence
            @reduce_valence 3 4
        end
    end
    if !isempty(invalidatoms)
        s = String.(collect(invalidatoms))
        @ifwarn @warn (dofix ? "After attempted fix, f" : "F")*"ound $(join(s, ", ", " and ")) with invalid number of bonds."
    end
    return invalidatoms
end



function sanity_checks!(graph, pos, types, mat, options)
    ## Bond length check
    removeedges = PeriodicEdge3D[]
    for e in edges(graph)
        s, d = e.src, e.dst.v
        bondlength = norm(mat *(pos[d] .+ e.dst.ofs .- pos[s]))
        if (bondlength < 0.65 && types[s] !== :H && types[d] !== :H)
            push!(removeedges, e)
        elseif bondlength > 3 && options.cutoff_coeff ≤ 0.85
            @ifwarn @warn "Suspiciously large bond found: $bondlength pm between $(types[s]) and $(types[d])."
            return true
        end
    end
    if options.bonding_mode == AutoBonds
        if !isempty(removeedges)
            @ifwarn begin
                @warn "Suspicious small bond lengths found. Such bonds are probably spurious and will be deleted."
                @info "To force retaining these bonds, use --bond-detect=input or --bond-detect=chemfiles"
            end
            if !(options.dryrun isa Nothing)
                options.dryrun[:try_noAutoBonds] = nothing
                options.dryrun[:suspect_smallbonds] = union(
                    Set{Symbol}([types[src(e)] for e in removeedges]),
                    Set{Symbol}([types[dst(e)] for e in removeedges]))
            end
        end
        for e in removeedges
            rem_edge!(graph, e)
        end
    end
    return false
end

function remove_monotonic_bonds!(graph::PeriodicGraph{D}, types, targets) where D
    isempty(targets) && return
    n = length(types)
    for i in 1:n
        t = types[i]
        t ∈ targets || continue
        rem_edges = PeriodicVertex{D}[]
        for x in neighbors(graph, i)
            types[x.v] == t || continue
            push!(rem_edges, x)
        end
        for x in rem_edges
            rem_edge!(graph, i, x)
        end
    end
    nothing
end


function parse_as_cif(cif::CIF, options, name)
    guessed_bonds = false
    n = length(cif.ids)
    ignored = Int[]
    types = Symbol[]
    pos = SVector{3,Float64}[]
    occupancies = options.ignore_low_occupancy ?
                        get(cif.cifinfo, :atom_site_occupancy, ones(n)) : ones(n)
    for i in 1:n
        typ = cif.types[cif.ids[i]]
        if typ ∈ options.ignore_atoms || occupancies[i] < 0.5
            push!(ignored, i)
            continue
        end
        push!(types, Symbol(parse_atom_name(String(typ))))
        push!(pos, cif.pos[:,i])
        #push!(pos, cif.cell.mat * cif.pos[:,i])
    end
    if all(==(Inf32), cif.bonds)
        if options.bonding_mode == InputBonds
            throw(ArgumentError("Cannot use input bonds since there are none. Use another option for --bonds-detect or provide bonds in the CIF file."))
        end
        guessed_bonds = true

        bonds = guess_bonds(pos, types, Float64.(cif.cell.mat), options.cutoff_coeff)
    else
        bonds = Tuple{Int,Int}[]
        _i = 1 # correction to i based on ignored atoms
        current_ignored_i = get(ignored, 1, 0)
        for i in 1:n
            if i == current_ignored_i
                _i += 1
                current_ignored_i = get(ignored, _i, 0)
                continue
            end
            _j = _i # correction to j based on ignored atoms
            current_ignored_j = current_ignored_i
            for j in (i+1):n
                if j == current_ignored_j
                    _j += 1
                    current_ignored_j = get(ignored, _j, 0)
                    continue
                end
                if cif.bonds[i,j] < Inf32
                    push!(bonds, (i - _i + 1, j - _j + 1))
                end
            end
        end
    end
    n -= length(ignored)

    cell = Cell(Cell(), cif.cell.mat)
    return finalize_checks(cell, pos, types, Int[], bonds, guessed_bonds, options, name)
end


function parse_as_chemfile(frame, options, name)
    types = Symbol[]
    for i in Int(size(frame))-1:-1:0
        typ = Symbol(type(Chemfiles.Atom(frame, i)))
        if typ ∈ options.ignore_atoms
            Chemfiles.remove_atom!(frame, i)
        else
            push!(types, typ)
        end
    end

    cell = Cell(Cell(), SMatrix{3,3,BigFloat}(matrix(UnitCell(frame)))')
    pos::Vector{SVector{3,Float64}} = Ref(inv(cell.mat)) .* eachcol(positions(frame))

    toremove = check_collision(pos, cell.mat)
    if !isempty(toremove)
        deleteat!(types, toremove)
        for j in Iterators.reverse(toremove)
            Chemfiles.remove_atom!(frame, j-1)
        end
        if !(options.dryrun isa Nothing)
            options.dryrun[:collisions] = nothing
        end
    end

    n = length(pos)
    guessed_bonds = false
    topology = Topology(frame)
    if bonds_count(topology) == 0
        guessed_bonds = true
        bonds = guess_bonds(pos, types, Float64.(cell.mat), options.cutoff_coeff)
        for (u,v) in bonds
            add_bond!(frame, u-1, v-1)
        end
    else
        bonds = [(a+1, b+1) for (a,b) in eachcol(Chemfiles.bonds(topology))]
        if !(options.dryrun isa Nothing) && options.bonding_mode == AutoBonds
            options.dryrun[:try_InputBonds] = nothing
        end
    end

    topology = Topology(frame) # Just a precaution since frame was possibly modified
    m = Int(count_residues(topology))
    residues = [Residue(topology, i) for i in 0:(m-1)]

    attributions = attribute_residues(residues, n, options.clustering == InputClustering)
    return finalize_checks(cell, pos, types, attributions, bonds, guessed_bonds, options, name)
end


function finalize_checks(cell, pos, types, attributions, bonds, guessed_bonds, options, name)
    if !all(isfinite, cell.mat) || iszero(det(cell.mat))
        @ifwarn @error "Suspicious unit cell of matrix $(Float64.(cell.mat)). Is the input really periodic? Using a cubic unit cell instead."
        cell = Cell()
    end

    n = length(pos)

    adjacency = falses(n, n)
    for (a,b) in bonds
        adjacency[a, b] = true
        adjacency[b, a] = true
    end
    mat = Float64.(cell.mat)
    graph = PeriodicGraph3D(n, edges_from_bonds(adjacency, mat, pos))

    if options.bonding_mode != InputBonds
        bad_valence = !isempty(fix_valence!(graph, pos, types, mat, Val(false)))
        recompute_bonds = bad_valence || sanity_checks!(graph, pos, types, mat, options)
        if recompute_bonds
            if !guessed_bonds
                @ifwarn begin
                    @warn "Disregarding all bonds from the input file."
                    @info "To force retaining the initial bonds, use --bond-detect=input"
                end
                bonds = guess_bonds(pos, types, mat, options.cutoff_coeff)
                adjacency = falses(n, n)
                for (a,b) in bonds
                    adjacency[a, b] = true
                    adjacency[b, a] = true
                end
                graph = PeriodicGraph3D(n, edges_from_bonds(adjacency, cell.mat, pos))
            end
            invalidatoms = fix_valence!(graph, pos, types, mat, Val(true))
            remaining_not_fixed = !isempty(invalidatoms)
            @ifwarn if remaining_not_fixed
                @warn "Remaining atoms with invalid valence. Proceeding anyway."
            end
            if !(options.dryrun isa Nothing)
                options.dryrun[:invalidatoms] = invalidatoms
            end
            recompute_bonds = sanity_checks!(graph, pos, types, cell.mat, options)
            @ifwarn if recompute_bonds
                @warn "Remaining bonds of suspicious lengths. Proceeding anyway."
            end
        end
    end

    remove_monotonic_bonds!(graph, types, options.ignore_monotonic_bonds)

    if isempty(attributions)
        crystalnothing = Crystal{Nothing}(cell, types, pos, graph, options)
        export_input(crystalnothing, name, options.export_input)
        return crystalnothing
    else
        crystalclusters = Crystal{Clusters}(cell, types, regroup_sbus(graph, attributions),
                                            pos, graph, options)
        export_input(crystalclusters, name, options.export_input)
        return crystalclusters
    end
end


"""
       parse_chemfile(path, options::Options)

Parse a file given in any reckognised chemical format and extract the topological
information.
Such format can be .cif or any file format reckognised by Chemfiles.jl that
contains all the necessary topological information.
"""
function parse_chemfile(_path, options::Options)
    # Separate the cases unhandled by Chemfiles from the others
    path = expanduser(_path)
    name = splitext(splitdir(path)[2])[1]
    if options.name == "unnamed"
        options = Options(options; name)
    end
    if lowercase(last(splitext(path))) == ".cif"
        cif = expand_symmetry(CIF(path))
        if options.authorize_pruning
            neededprune, cif = prune_collisions(cif)
            if neededprune && !(options.dryrun isa Nothing)
                options.dryrun[:collisions] = nothing
            end
        end
        return parse_as_cif(cif, options, name)
    end
    frame = read(Trajectory(path))
    return parse_as_chemfile(frame, options, name)
end
parse_chemfile(path; kwargs...) = parse_chemfile(path, Options(; kwargs...))
