import PeriodicGraphs

function export_dataline(f, x)
    inbetween = '\t'
    if any(isspace, x)
        if any(==('\n'), x)
            x = ';'*x*"\n;"
            inbetween = '\n'
        else
            @toggleassert x[1]!= '\''
            x = '\''*x*'\''
        end
    end
    println(f, inbetween*x)
end

function export_vtf(file, __c::Union{Crystal,CrystalNet}, repeatedges=6, colorname=false)
    mkpath(splitdir(file)[1])
    n = length(__c.types)
    @toggleassert length(__c.pge.pos) == n
    c = __c isa CrystalNet ? CrystalNet3D(__c) : __c
    open(file, write=true) do f
        invcorres = [PeriodicVertex3D(i) for i in 1:n]
        corres = Dict{PeriodicVertex3D,Int}([invcorres[i]=>i for i in 1:n])
        encounteredtypes = Dict{Symbol,String}()
        atomnums = Int[]
        numencounteredtypes = 0
        println(f, """
        ###############################
        # written by CrystalNets.jl
        ###############################
        """)

        for i in 1:n
            ty = c.types[i]
            sty = ty === Symbol("") ? string(i) : string_atomtype(ty)
            if length(sty) > 16
                sty = sty[1:13]*"etc" # otherwise VMD fails to load the .vtf
            end
            name = if colorname
                _name = get(encounteredtypes, ty, missing)
                if _name isa String
                    _name
                else
                    numencounteredtypes += 1
                    encounteredtypes[ty] = string(" name ", numencounteredtypes)
                end
            else
                n ≥ 32768 ? "" : string(" name ", i)
            end
            atomnum = ty === Symbol("") ? 0 : representative_atom(ty, i)[2]
            push!(atomnums, atomnum)
            resid = colorname ? i : 0
            println(f, "atom $(i-1) type $sty$name resid $resid atomicnumber $atomnum")
        end
        j = n + 1
        for _ in 1:repeatedges
            jmax = j - 1
            for i in 1:jmax
                vertex = invcorres[i]
                for x in neighbors(c.pge.g, vertex.v)
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
            sty = ty === Symbol("") ? string(i) : string_atomtype(ty)
            if length(sty) > 16
                sty = sty[1:13]*"etc"
            end
            name = colorname ? string(" name ", encounteredtypes[ty]) :
                    n ≥ 32768 ? "" : string(" name ", v)
            atomnum = atomnums[v]
            resid = colorname ? i : PeriodicGraphs.hash_position(ofs)
            println(f, "atom $(i-1) type $sty$name resid $resid atomicnumber $atomnum")
        end
        println(f)

        for (i,x) in enumerate(invcorres)
            for neigh in neighbors(c.pge.g, x.v)
                j = get(corres, PeriodicVertex3D(neigh.v, neigh.ofs .+ x.ofs), nothing)
                isnothing(j) && continue
                if i < j
                    println(f, "bond ", i - 1, ':', j - 1)
                end
            end
        end
        println(f)

        ((_a, _b, _c), (_α, _β, _γ)), mat = cell_parameters(c.pge.cell)
        println(f, "pbc $_a $_b $_c $_α $_β $_γ\n")

        println(f, "ordered")
        for x in invcorres
            coord = mat * (widen.(c.pge.pos[x.v]) .+ x.ofs)
            join(f, round.(Float64.(coord); digits=15), ' ')
            println(f)
        end
    end
end


function export_cif(file, __c::Union{Crystal, CIF})
    mkpath(splitdir(file)[1])
    c, info = if __c isa CIF
        ___c = expand_symmetry(__c)
        ___c, copy(__c.cifinfo)
    else
        __c, Dict{String,Union{String,Vector{String}}}("data"=>last(splitdir(file)))
    end
    loops = Dict{Int, Tuple{Vector{String}, Vector{Vector{String}}}}()
    open(file, write=true) do f
        print(f, "data_"); println(f, popstring!(info, "data")*"\n\n")
        for (id, data) in info
            if data isa String
                print(f, '_'); print(f, id)
                export_dataline(f, data)
            else
                @toggleassert data isa Vector{String}
                n = length(data)
                if haskey(loops, n)
                    push!(loops[n][1], id)
                    push!(loops[n][2], data)
                else
                    loops[n] = (String[id], Vector{String}[data])
                end
            end
        end

        cell::Cell = c isa CIF ? c.cell : c.pge.cell
        #scale_factor::Float64 = unique!(sort(c.types)) == [:Si] && c.pge.pos[:,1] != [0,0,0] ? 1.2 : 1.0
        ((__a, __b, __c), (__α, __β, __γ)), _ = cell_parameters(cell)
        _a, _b, _c, _α, _β, _γ = Float64.((__a, __b, __c, __α, __β, __γ))

        println(f, """

        _cell_length_a\t\t$_a
        _cell_length_b\t\t$_b
        _cell_length_c\t\t$_c
        _cell_angle_alpha\t\t$_α
        _cell_angle_beta\t\t$_β
        _cell_angle_gamma\t\t$_γ

        _symmetry_space_group_name_H-M\t'P 1'
        _symmetry_Int_Tables_number\t1
        _symmetry_cell_setting\ttriclinic
        """)

        println(f, """
        loop_
        _symmetry_equiv_pos_as_xyz
        x,y,z""")

        for eq in cell.equivalents
            println(f, eq)
        end
        println(f)

        n::Int = c isa CIF ? length(c.ids) : length(c.types)
        if !haskey(loops, n)
            loops[n] = (String[], Vector{String}[])
        end
        append!(loops[n][1],
                ["atom_site_label", "atom_site_type_symbol",
                 "atom_site_fract_x", "atom_site_fract_y", "atom_site_fract_z"])
        labels = String[string_atomtype(c.types[c isa CIF ? c.ids[i] : i]::Symbol)*string(i) for i in 1:n]
        pos = string.(round.(c isa CIF ? c.pos::Matrix{Float64} : reduce(hcat, c.pge.pos)::Matrix{Float64}; sigdigits=6))
        append!(loops[n][2], [labels, string_atomtype.((c.types[x])::Symbol for x in (c isa CIF ? c.ids : collect(1:n))::Vector{Int}),
                pos[1,:], pos[2,:], pos[3,:]])

        if c isa Union{Crystal{Clusters}, Crystal{Nothing}}
            bonds = edges((c isa CIF ? c.graph : c.pge.g)::PeriodicGraph3D)
            src = String[]
            dst = String[]
            mult = Int[]
            last_bond = (0, 0)
            n = 0
            for e in bonds
                if (e.src, e.dst.v) == last_bond
                    mult[end] += 1
                else
                    push!(src, labels[e.src])
                    push!(dst, labels[e.dst.v])
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
        #scale_factor::Float64 = unique!(sort(c.types)) == [:Si] && c.pge.pos[:,1] != [0,0,0] ? 1.2 : 1.0
        ((__a, __b, __c), (__α, __β, __γ)), _ = cell_parameters(c.pge.cell)
        _a, _b, _c, _α, _β, _γ = Float64.((__a, __b, __c, __α, __β, __γ))
        print(f, "\tGROUP\t\"")
        print(f, RAW_SYMMETRY_DATA[c.pge.cell.hall][4])
        println("\"")
        println(f, "\tCELL\t", _a, ' ', _b, ' ', _c, ' ', _α, ' ', _β, ' ', _γ, ' ')
        println(f, "\tATOM")
        to_revisit = Int[]
        for i in 1:length(c.types)
            push!(to_revisit, i)
            pos = c.pge.pos[i]
            println(f, "\t\t", i, ' ', degree(c.pge.g, i), ' ', pos[1], ' ',
                    pos[2], ' ', pos[3])
        end
        println(f, "\tEDGE")
        for i in to_revisit
            for e in neighbors(c.pge.g, i)
                e.v < i && continue
                dest = c.pge.pos[e.v] .+ e.ofs
                println(f, "\t\t", i, '\t', dest[1], ' ', dest[2], ' ', dest[3])
            end
        end
        println(f, "\nEND")
    end
end

function export_cgd(file, g::PeriodicGraph)
    mkpath(splitdir(file)[1])
    f = open("/tmp/foobar", "w")#open(file, write=true) do f
        println(f, "PERIODIC_GRAPH\n")
        println(f, "ID ", basename(splitext(file)[1]), '\n')
        println(f, "EDGES")
        repr = reverse(split(string(g)))
        n = parse(Int, pop!(repr))
        m = length(repr) ÷ (n+2)
        @toggleassert iszero(length(repr) % (n+2))
        for _ in 1:m
            src = pop!(repr)
            dst = pop!(repr)
            ofs = Vector{String}(undef, n)
            for i in 1:n
                ofs[i] = pop!(repr)
            end
            print(f, '\t', src, ' ', dst, ' ')
            join(f, ofs, ' ')
            println(f)
        end
        println(f, "END\n")
    #end
end
export_cgd(file, c::CrystalNet) = export_cgd(file, change_dimension(PeriodicGraph3D, c.pge.g))


function export_attributions(crystal::Crystal{Clusters}, path=joinpath(tempdir(),tempname()))
    frame = Chemfiles.Frame()
    m = length(crystal.clusters.classes)
    # residues = [Chemfiles.Residue(string(i)) for i in 1:m]
    # for i in 1:m
    #     Chemfiles.set_property!(residues[i], "chainname", string(i))
    #     Chemfiles.set_property!(residues[i], "chainid", string(i))
    #     Chemfiles.add_residue!(frame, residues[i])
    # end
    ((_a, _b, _c), (_α, _β, _γ)), mat = cell_parameters(c.pge.cell)
    Chemfiles.set_cell!(frame, Chemfiles.UnitCell(Float64[_a, _b, _c], Float64[_α, _β, _γ]))
    recenter::SVector{3,Float64} = minimum(reduce(hcat, crystal.pge.pos); dims=2)
    for i in 1:length(crystal.types)
        # resid = crystal.clusters.attributions[i]
        atom = Chemfiles.Atom(string(crystal.clusters.classes[crystal.clusters.attributions[i]]))
        Chemfiles.set_type!(atom, string_atomtype(crystal.types[i]))
        Chemfiles.add_atom!(frame, atom, collect(Float64.(mat * (crystal.pge.pos[i] .- recenter))))
        # Chemfiles.add_atom!(residues[resid], i)
    end
    for e in edges(crystal.pge.g)
        iszero(PeriodicGraphs.ofs(e)) || continue
        Chemfiles.add_bond!(frame, src(e)-1, dst(e)-1)
    end
    target = isempty(last(splitext(path))) ? path*".pdb" : path
    output = Chemfiles.Trajectory(target, 'w')
    write(output, frame)
    close(output)
end


"""
    export_arc(path)

Export the current archive to the specified path.
"""
function export_arc(path, empty=false, arc=CRYSTAL_NETS_ARCHIVE)
    mkpath(splitdir(path)[1])
    open(path, "w") do f
        println(f, "Made by CrystalNets.jl v$CRYSTAL_NETS_VERSION")
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
