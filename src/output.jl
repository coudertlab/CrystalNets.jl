import PeriodicGraphs

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

function export_vtf(file, c::CrystalNet, repeatedges=1, colorname=false)
    mkpath(splitdir(file)[1])
    n = length(c.types)
    @assert length(c.pos) == n
    invcorres = [PeriodicVertex3D(i) for i in 1:n]
    corres = Dict{PeriodicVertex3D,Int}([invcorres[i]=>i for i in 1:n])
    encounteredtypes = Dict{Symbol,String}()
    numencounteredtypes = 0
    open(file, write=true) do f
        println(f, """
        ###############################
        # written by PeriodicGraphs.jl
        ###############################
        """)

        for i in 1:n
            ty = c.types[i]
            sty = ty === Symbol("") ? string(i) : string(ty)
            name = colorname ? get!(encounteredtypes, ty) do
                string(numencounteredtypes+=1)
            end : string(i)
            println(f, "atom $(i-1) type $sty name $name resid $(i-1)")
        end
        j = n + 1
        for _ in 1:repeatedges
            jmax = j - 1
            for i in 1:jmax
                vertex = invcorres[i]
                for x in neighbors(c.graph, vertex.v)
                    y = PeriodicVertex3D(x.v, x.ofs .+ vertex.ofs)
                    if get!(corres, y, j) == j
                        j += 1
                        push!(invcorres, y)
                    end
                end
            end
        end
        for i in n+1:length(invcorres)
            v = invcorres[i].v
            ofs = invcorres[i].ofs
            ty = c.types[v]
            sty = ty === Symbol("") ? string(i) : string(ty)
            name = colorname ? encounteredtypes[ty] : string(v)
            println(f, "atom $(i-1) type $sty name $name resid $v")
        end
        println(f)

        for (i,x) in enumerate(invcorres)
            for neigh in neighbors(c.graph, x.v)
                j = get(corres, PeriodicVertex3D(neigh.v, neigh.ofs .+ x.ofs), nothing)
                isnothing(j) && continue
                if i < j
                    println(f, "bond ", i - 1, ':', j - 1)
                end
            end
        end
        println(f)

        _a, _b, _c, α, β, γ = cell_parameters(c.cell)
        println(f, "pbc $_a $_b $_c $α $β $γ\n")

        println(f, "ordered")
        for x in invcorres
            coord = c.cell.mat * (widen.(c.pos[x.v]) .+ x.ofs)
            println(f, join(round.(Float64.(coord); digits=15), ' '))
        end
    end
end

function export_vtf(file, cif::CIF, repeatedges=1)
    export_vtf(file, CrystalNet(cif), repeatedges, true)
end

function export_cif(file, c::Union{Crystal, CIF})
    mkpath(splitdir(file)[1])
    info = c isa CIF ? copy(c.cifinfo) : Dict{String,Union{String,Vector{String}}}("data"=>last(splitdir(file)))
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
                    loops[n] = (String[id], Vector{String}[data])
                end
            end
        end

        #scale_factor::Float64 = unique!(sort(c.types)) == [:Si] && c.pos[:,1] != [0,0,0] ? 1.2 : 1.0
        _a, _b, _c, α, β, γ = Float64.(cell_parameters(c.cell))

        println(f, """

        _cell_length_a\t\t$_a
        _cell_length_b\t\t$_b
        _cell_length_c\t\t$_c
        _cell_angle_alpha\t\t$α
        _cell_angle_beta\t\t$β
        _cell_angle_gamma\t\t$γ

        _symmetry_space_group_name_H-M\t'$(c.cell.spacegroup)'""")
        if c.cell.tablenumber != 0
            println(f, "_symmetry_Int_Tables_number\t$(c.cell.tablenumber)")
        end
        if c.cell.latticesystem !== Symbol("")
            println(f, "_symmetry_cell_setting\t$(c.cell.latticesystem)")
        end

        println(f, """

        loop_
        _symmetry_equiv_pos_as_xyz
        x,y,z""")
        for eq in c.cell.equivalents
            println(f, eq)
        end
        println(f)

        n = c isa CIF ? length(c.ids) : length(c.types)
        if !haskey(loops, n)
            loops[n] = (String[], Vector{String}[])
        end
        append!(loops[n][1],
                ["atom_site_label", "atom_site_type_symbol",
                 "atom_site_fract_x", "atom_site_fract_y", "atom_site_fract_z"])
        labels = String[string(c.types[c isa CIF ? c.ids[i] : i])*string(i) for i in 1:n]
        pos = string.(round.(c.pos; sigdigits=6))
        append!(loops[n][2], [labels, string.(c.types[x] for x in (c isa CIF ? c.ids : collect(1:n))),
                pos[1,:], pos[2,:], pos[3,:]])

        if c isa Crystal
            bonds = edges(c.graph)
            src = String[]
            dst = String[]
            mult = Int[]
            last_bond = (0, 0)
            n = 0
            for e in bonds
                if (e.src, e.dst.v) == last_bond
                    mult[end] += 1
                else
                    push!(src, string(labels[e.src]))
                    push!(dst, string(labels[e.dst.v]))
                    push!(mult, 1)
                    last_bond = (e.src, e.dst.v)
                    n += 1
                end
            end
            if !haskey(loops, n)
                loops[n] = (String[], Vector{String}[])
            end
            append!(loops[n][1], ["geom_bond_atom_site_label_1", "geom_bond_atom_site_label_2", "geom_bond_multiplicity"])
            append!(loops[n][2], [src, dst, string.(mult)])
        end

        for (n, (ids, datas)) in loops
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

function export_cgd(file, c::Crystal)
    mkpath(splitdir(file)[1])
    open(file, write=true) do f
        println(f, "CRYSTAL\n")
        margin = 1
        println(f, "\tNAME\t", last(splitdir(file)))
        #scale_factor::Float64 = unique!(sort(c.types)) == [:Si] && c.pos[:,1] != [0,0,0] ? 1.2 : 1.0
        _a, _b, _c, α, β, γ = Float64.(cell_parameters(c.cell))
        println(f, "\tGROUP\t\"", join(split(c.cell.spacegroup, ' ')), "\"")
        println(f, "\tCELL\t", _a, ' ', _b, ' ', _c, ' ', α, ' ', β, ' ', γ, ' ')
        println(f, "\tATOM")
        to_revisit = Int[]
        for i in 1:length(c.ids)
            # c.ids[i] ∈ (@view c.ids[to_revisit]) && continue
            push!(to_revisit, i)
            println(f, "\t\t", i, ' ', degree(c.graph, i), ' ', c.pos[1,i], ' ',
                    c.pos[2,i], ' ', c.pos[3,i])
        end
        println(f, "\tEDGE")
        for i in to_revisit
            for e in neighbors(c.graph, i)
                e.v < i && continue
                dest = c.pos[:,e.v] .+ e.ofs
                println(f, "\t\t", i, '\t', dest[1], ' ', dest[2], ' ', dest[3])
            end
        end
        println(f, "\nEND")
    end
end

function export_cgd(file, g::PeriodicGraph)
    mkpath(splitdir(file)[1])
    open(file, write=true) do f
        println(f, "PERIODIC_GRAPH\n")
        println(f, "ID ", basename(splitext(file)[1]), '\n')
        println(f, "EDGES")
        repr = reverse(split(string(g)))
        n = parse(Int, pop!(repr))
        m = length(repr) ÷ (n+2)
        @assert iszero(length(repr) % (n+2))
        for _ in 1:m
            src = pop!(repr)
            dst = pop!(repr)
            ofs = Vector{String}(undef, n)
            for i in 1:n
                ofs[i] = pop!(repr)
            end
            println(f, '\t', src, ' ', dst, ' ', join(ofs, ' '))
        end
        println(f, "END\n")
    end
end


function export_clusters(crystal::Crystal{Clusters}, path=joinpath(tempdir(),tempname()))
    frame = Chemfiles.Frame()
    m = length(crystal.clusters.classes)
    # residues = [Chemfiles.Residue(string(i)) for i in 1:m]
    # for i in 1:m
    #     Chemfiles.set_property!(residues[i], "chainname", string(i))
    #     Chemfiles.set_property!(residues[i], "chainid", string(i))
    #     Chemfiles.add_residue!(frame, residues[i])
    # end
    _a, _b, _c, α, β, γ = cell_parameters(crystal.cell)
    Chemfiles.set_cell!(frame, Chemfiles.UnitCell(_a, _b, _c, α, β, γ))
    recenter::SVector{3,Float64} = minimum(crystal.pos; dims=2)
    for i in 1:length(crystal.types)
        # resid = crystal.clusters.attributions[i]
        atom = Chemfiles.Atom(string(crystal.clusters.classes[crystal.clusters.attributions[i]]))
        Chemfiles.set_type!(atom, string(crystal.types[i]))
        Chemfiles.add_atom!(frame, atom, collect(Float64.(crystal.cell.mat * (crystal.pos[:,i] .- recenter))))
        # Chemfiles.add_atom!(residues[resid], i)
    end
    for e in edges(crystal.graph)
        iszero(PeriodicGraphs.ofs(e)) || continue
        Chemfiles.add_bond!(frame, src(e)-1, dst(e)-1)
    end
    target = isempty(last(splitext(path))) ? path*".pdb" : path
    output = Chemfiles.Trajectory(target, 'w')
    write(output, frame)
    close(output)
    return target
end


"""
    export_arc(path)

Export the current archive to the specified path.
"""
function export_arc(path, empty=false, arc=CRYSTAL_NETS_ARCHIVE)
    major = CRYSTAL_NETS_VERSION.major
    minor = CRYSTAL_NETS_VERSION.minor
    open(path, "w") do f
        println(f, "Made by CrystalNets.jl v$major.$minor")
        if !empty
            println(f)
            pairs = sort(collect(arc); by=last)
            for (genome, id) in pairs
                println(f, "key\t", genome)
                println(f, "id\t", id)
                println(f)
            end
        end
    end
end
