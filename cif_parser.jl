

include("types.jl")
using .CIFTypes

"""
    nextword(l, k)

Return the pair of indices (i, j) such that l[i:j] is the next word in the
string l after position k.
Use k = j to get the following word, and so forth.

(0, 0) is returned if there is no next word.
"""
function nextword(l, i)
    n = lastindex(l)
    i == n && return (0, 0)

    i = nextind(l, i)
    start = 0
    if l[i] == '\'' || l[i] == '"'
        if l[nextind(l, i)] == '#'
            i = findnext(isequal('\n'), l, i)
            i == nothing && return (0, 0)
        else
            @assert isspace(l[nextind(l, i)])
            i = nextind(l, nextind(l, i))
        end
    end
    while i <= n
        if isspace(l[i])
            i = nextind(l, i)
        else
            if l[i] == '#'
                i = findnext(isequal('\n'), l, i)
                i == nothing && return (0, 0)
            else
                start = i
                break
            end
        end
    end
    start == 0 && return (0, 0)

    inquote = false
    quotesymb = "'"
    while i <= n
        c = l[i]
        if !isspace(c)
            c == '#' && return (start, prevind(l, i))
            if c == '\'' || c == '"'
                inquote = true
                quotesymb = c
                break
            end
            i = nextind(l, i)
        else
            return (start, prevind(l, i))
        end
    end
    if inquote
        while i < n
            c = l[i]
            if c == quotesymb && (isspace(l[nextind(l, i)]) ||
                                 l[nextind(l, i)] == '#')
                return (start+1, prevind(l, i))
            end
            i = nextind(l, i)
        end
        if !l[n] == quotesymb
            throw("Invalid syntax: opening quote $quotesymb at position $start is not closed")
        end
        return (start+1, prevind(l, i))
    end
    return (start, n)
end

"""
    parse_cif(file)

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
    i, j = nextword(l, 0)
    while i != 0

        if inloop
            if !loopisspecified
                if l[i] != '_' # This indicates the start of the values
                    loop_n = length(loopspec)
                    loopisspecified = true
                else # The identifier is part of the loop specification
                    push!(loopspec, l[i+1:j])
                    all_data[l[i+1:j]] = String[]
                    i, j = nextword(l, j); continue
                end
            end

            # From this point, the loop has been specified
            @assert loopisspecified
            if l[i] != '_' && l[i:j] != "loop_"
                for k in 1:loop_n
                    push!(all_data[loopspec[k]], l[i:j])
                    i, j = nextword(l, j)
                end
                continue
            end
            if l[i:j] == "loop_"
                loopisspecified = false
                loopspec = String[]
                i, j = nextword(l, j); continue
            end

            # This point can only be reached if we just quitted a loop
            inloop = false
        end

        @assert !inloop
        if l[i] == '_' # Simple identifier definition
            next_i, next_j = nextword(l, j)
            @assert next_i != 0
            all_data[l[i+1:j]] = l[next_i:next_j]
            j = next_j
        else
            if l[i:j] == "loop_"
                inloop = true
                loopisspecified = false
                loopspec = String[]
            elseif j-i > 4 && l[i:5] == "data_"
                # simply skip
            else
                k = findprev(isequal('\n'), l, i)
                n = count("\n", l[1:k])
                throw("Unkown word \"$(l[i:j])\" at line $n, position $(i-k):$(j-k)")
            end
        end

        i, j = nextword(l, j)
    end

    return all_data
end

"""
    CIF(parsed_data)

Make a CIF object out of the parsed file.
"""
function CIF(parsed)
    natoms = length(parsed["atom_site_label"])
    cell = Cell(Symbol(parsed["symmetry_cell_setting"]),
                parsed["symmetry_space_group_name_H-M"],
                parse(Int, parsed["symmetry_Int_Tables_number"]),
                parse(Float64, parsed["cell_length_a"]),
                parse(Float64, parsed["cell_length_b"]),
                parse(Float64, parsed["cell_length_c"]),
                parse(Float64, parsed["cell_angle_alpha"]),
                parse(Float64, parsed["cell_angle_beta"]),
                parse(Float64, parsed["cell_angle_gamma"]))

    labels = parsed["atom_site_label"]
    symbols = parsed["atom_site_type_symbol"]
    pos_x = parsed["atom_site_fract_x"]
    pos_y = parsed["atom_site_fract_y"]
    pos_z = parsed["atom_site_fract_z"]

    atoms = Symbol[]
    pos = Matrix{Float64}(undef, 3, natoms)
    correspondence = Dict{String, Int}()
    atom_types = Symbol[]
    for i in 1:natoms
        @assert !haskey(correspondence, labels[i])
        correspondence[labels[i]] = i
        push!(atoms, Symbol(symbols[i]))
        pos[:,i] = parse.(Float64, [pos_x[i], pos_y[i], pos_z[i]])
    end
    bond_a = parsed["geom_bond_atom_site_label_1"]
    bond_b = parsed["geom_bond_atom_site_label_2"]
    bonds = zero(BitMatrix(undef, natoms, natoms))
    for i in 1:length(bond_a)
        x = correspondence[bond_a[i]]
        y = correspondence[bond_b[i]]
        bonds[x,y] = bonds[y,x] = 1
    end
    return CIF(natoms, cell, atoms, bonds, pos)
end

Crystal(parsed::Dict) = Crystal(CIF(parsed))
