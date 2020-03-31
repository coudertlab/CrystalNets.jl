includet("types.jl")
using .CIFTypes

function export_dataline(f, x)
    inbetween = '\t'
    if any(isspace, x)
        if any(==('\n'), x)
            x = ';'*x*"\n;"
            inbetween = '\n'
        else
            @assert x[1]!= '\''
            x = '\''*x*'\''
        end
    end
    println(f, inbetween*x)
end

function export_cif(file, c::Union{Crystal, CIF})
    info = copy(c.cifinfo)
    loops = Dict{Int, Tuple{Vector{String}, Vector{Vector{String}}}}()
    open(file, write=true) do f
        print(f, "data_"); println(f, pop!(info, "data")*"\n\n")
        for (id, data) in info
            if data isa String
                print(f, '_'); print(f, id)
                export_dataline(f, data)
            else
                @assert data isa Vector{String}
                n = length(data)
                if haskey(loops, n)
                    push!(loops[n][1], id)
                    push!(loops[n][2], data)
                else
                    loops[n] = [String[id], Vector{String}[data]]
                end
            end
        end

        println(f, """

        _cell_length_a\t\t$(c.geometry.a/1.2)
        _cell_length_b\t\t$(c.geometry.b/1.2)
        _cell_length_c\t\t$(c.geometry.c/1.2)
        _cell_angle_alpha\t\t$(c.geometry.α)
        _cell_angle_beta\t\t$(c.geometry.β)
        _cell_angle_gamma\t\t$(c.geometry.γ)

        _symmetry_space_group_name_H-M\t'$(c.geometry.spacegroup)'
        _symmetry_Int_Tables_number\t$(c.geometry.tablenumber)
        _symmetry_cell_setting\t$(c.geometry.latticesystem)

        loop_
        _symmetry_equiv_pos_as_xyz
        x, y, z""")
        for eq in c.geometry.equivalents
            println(f, eq)
        end
        println(f)

        if !haskey(loops, c.natoms)
            loops[c.natoms] = (String[], Vector{String}[])
        end
        append!(loops[c.natoms][1],
                ["_atom_site_label", "atom_site_type_symbol",
                 "atom_site_fract_x", "atom_site_fract_y", "atom_site_fract_z"])
        labels = String[string(c.atoms[i])*string(i) for i in 1:c.natoms]
        pos = string.(round.(c.pos; sigdigits=6))
        append!(loops[c.natoms][2], [labels, string.(c.atoms), pos[1,:], pos[2,:], pos[3,:]])

        if c isa Crystal
            bonds = edges(c.graph)
            n = length(bonds)
            if !haskey(loops, n)
                loops[n] = (String[], Vector{String}[])
            end
            append!(loops[n][1], ["_geom_bond_atom_site_label_1", "_geom_bond_atom_site_label_2"])
            append!(loops[n][2], [String[], String[]])
            for e in bonds
                push!(loops[n][2][end-1], string(labels[e.src]))
                push!(loops[n][2][end], string(labels[e.dst.v]))
            end
        end

        for (n, (ids, datas)) in sort(loops)
            println(f, "loop_")
            for id in ids
                print(f, '_'); println(f, id)
            end
            for i in 1:n
                for j in 1:length(datas)
                    print(f, datas[j][i]); print(f, '\t')
                end
                println(f)
            end
            println(f)
        end
        println(f, "\n#END")
    end
    nothing
end
