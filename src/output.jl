import PeriodicGraphEmbeddings: export_vtf

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

function perturbate_positions(pge::PeriodicGraphEmbedding{D}) where {D}
    lens, _ = cell_parameters(pge.cell.mat)
    newpos = [pos .+ rand((-1,1), 3) .* rand(3) ./ (8 .* lens) for pos in pge.pos]
    return PeriodicGraphEmbedding{D}(copy(pge.g), newpos, pge.cell)
end

_representative_atom(ty, i) = last(representative_atom(ty, i))
function export_vtf(file, c::Union{CrystalNet,Crystal}, repeatedges=6, perturb=false, colorname=false)
    if c isa CrystalNet
        net = CrystalNet3D(c)
        pge = perturb ? perturbate_positions(net.pge) : net.pge
        return export_vtf(file, pge, net.types, repeatedges, colorname, string_atomtype, _representative_atom)
    end
    pge = perturb ? perturbate_positions(c.pge) : c.pge
    return export_vtf(file, pge, c.types, repeatedges, colorname, string_atomtype, _representative_atom)
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

        cell::Cell{Rational{Int}} = c isa CIF ? c.cell : c.pge.cell
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

PeriodicGraphEmbeddings.export_cgd(file, c::Crystal) = export_cgd(file, c.pge)
PeriodicGraphEmbeddings.export_cgd(file, c::CrystalNet) = export_cgd(file, PeriodicGraph3D(c.pge.g))


function export_attributions(crystal::Crystal{Clusters}, path=joinpath(tempdir(),tempname()))
    frame = Chemfiles.Frame()
    m = length(crystal.clusters.classes)
    # residues = [Chemfiles.Residue(string(i)) for i in 1:m]
    # for i in 1:m
    #     Chemfiles.set_property!(residues[i], "chainname", string(i))
    #     Chemfiles.set_property!(residues[i], "chainid", string(i))
    #     Chemfiles.add_residue!(frame, residues[i])
    # end
    ((_a, _b, _c), (_α, _β, _γ)), mat = cell_parameters(crystal.pge.cell)
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
        iszero(last(last(e))) || continue
        Chemfiles.add_bond!(frame, first(e)-1, first(last(e))-1)
    end
    target = isempty(last(splitext(path))) ? path*".pdb" : path
    output = Chemfiles.Trajectory(target, 'w')
    write(output, frame)
    close(output)
end

function _export_trimmed_and_attributions(crystal::Crystal{Nothing}, clusters::Clusters)
    export_default(crystal, "trimmed", crystal.options.name, crystal.options.export_trimmed)
    if !isempty(crystal.options.export_attributions)
        path = tmpexportname(crystal.options.export_attributions, "attribution_", crystal.options.name, ".pdb")
        export_attributions(Crystal{Clusters}(crystal, clusters), path)
        println("Attributions of atoms into SBUs represented represented at ", replace(path, ('\\'=>'/')))
    end
    nothing
end

"""
    export_arc(path, arc=CRYSTAL_NETS_ARCHIVE)

Export archive `arc` to the specified `path`. If unspecified, the exported archive is the
current one.
"""
function export_arc(path, arc=CRYSTAL_NETS_ARCHIVE)
    mkpath(splitdir(path)[1])
    open(path, "w") do f
        println(f, "Made by CrystalNets.jl v$CRYSTAL_NETS_VERSION")
        if !isnothing(arc)
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
