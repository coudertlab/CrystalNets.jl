include("types.jl")
using .CIFTypes
using .CIFTypes.PeriodicGraphs
import LinearAlgebra: norm

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

function export_vtf(file, c::CrystalNet, repeatedges=1)
    mkpath(splitdir(file)[1])
    n = length(c.types)
    @assert length(c.pos) == n
    invcorres = [PeriodicVertex3D(i) for i in 1:n]
    corres = Dict{PeriodicVertex3D,Int}([invcorres[i]=>i for i in 1:n])
    open(file, write=true) do f
        println(f, """
        ###############################
        # written by PeriodicGraphs.jl
        ###############################
        """)

        for i in 1:n
            ty = c.types[i]
            println(f, "atom $(i-1) name $i type ",
                    ty == Symbol("") ? string(i) : string(ty), " resid $(i-1)")
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
            println(f, "atom $(i-1) name $v type ",
                    ty == Symbol("") ? string(i) : string(ty), " resid $(i-1)")
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
        #
        # for e in edges(c.graph)
        #     print(f, "bond ", e.src - 1, ':', corres[e.dst] - 1)
        #     if iszero(e.dst.ofs)
        #         println(f)
        #     else
        #         println(f, ", ", corres[PeriodicVertex3D(e.src, .-e.dst.ofs)] - 1, ':', e.dst.v - 1)
        #     end
        # end
        println(f)

        _a, _b, _c, α, β, γ = CIFTypes.cell_parameters(c.cell)
        println(f, "pbc $_a $_b $_c $α $β $γ\n")

        # axis1 = (c.cell.mat[:,1] .+ [_a, 0, 0]) ./ 2
        # if norm(axis1) <= 1e-8 # The origin is at the middle point
        #     axis1 .= [0, 1, 0]
        # end
        # rotation1 = 2*axis1*axis1'/(axis1'axis1) - LinearAlgebra.I(3)
        # @assert inv(rotation1) ≈ rotation1 ≈ rotation1'
        # mat = rotation1 * c.cell.mat
        # if abs(mat[3,2]) > 1e-8
        #     θ = sign(mat[3,2]*mat[2,2])*π/2 + atan(mat[2,2] / mat[3,2])
        #     mat = [1 0 0; 0 cos(θ) -sin(θ); 0 sin(θ) cos(θ)] * mat
        # end
        mat = c.cell.mat

        println(f, "ordered")
        for x in invcorres
            coord = mat * (c.pos[x.v] .+ x.ofs)
            println(f, join(round.(Float64.(coord); digits=15), ' '))
        end
    end
end

function export_cif(file, c::Union{Crystal, CIF})
    mkpath(splitdir(file)[1])
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
                    loops[n] = (String[id], Vector{String}[data])
                end
            end
        end

        #scale_factor::Float64 = unique!(sort(c.types)) == [:Si] && c.pos[:,1] != [0,0,0] ? 1.2 : 1.0
        _a, _b, _c, α, β, γ = Float64.(CIFTypes.cell_parameters(c.cell))

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

        n = length(c.ids)
        if !haskey(loops, n)
            loops[n] = (String[], Vector{String}[])
        end
        append!(loops[n][1],
                ["atom_site_label", "atom_site_type_symbol",
                 "atom_site_fract_x", "atom_site_fract_y", "atom_site_fract_z"])
        labels = String[string(c.types[c.ids[i]])*string(i) for i in 1:n]
        pos = string.(round.(c.pos; sigdigits=6))
        append!(loops[n][2], [labels, string.(c.types[x] for x in c.ids),
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

function export_cgd(file, c::Crystal)
    mkpath(splitdir(file)[1])
    open(file, write=true) do f
        println(f, "CRYSTAL\n")
        margin = 1
        println(f, "\tNAME\t", c.cifinfo["data"])
        #scale_factor::Float64 = unique!(sort(c.types)) == [:Si] && c.pos[:,1] != [0,0,0] ? 1.2 : 1.0
        _a, _b, _c, α, β, γ = Float64.(CIFTypes.cell_parameters(c.cell))
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
