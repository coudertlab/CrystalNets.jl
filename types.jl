module CIFTypes

include("./PeriodicGraphs.jl")
include("./rationallu.jl")
using Serialization
using StaticArrays, Tokenize
using .PeriodicGraphs
import LightGraphs: nv, neighbors, is_connected
import LinearAlgebra: norm, issuccess, det
import Statistics: mean
export EquivalentPosition, Cell, CIF, Crystal, CrystalNet

struct EquivalentPosition
    mat::SMatrix{3,3,Rational{Int}, 9}
    ofs::SVector{3,Rational{Int}}
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
                    if !iszero(ofs[i])
                        @warn "Existing offset already existing for position $i in {$s}"
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
                    i > 2 && throw("Too many dimensions specified for symmetry equivalent {$s}")
                    something_written || throw("{$s} is not a symmetry equivalent (no dependency expressed in position $i)")
                    something_written = false
                    i += 1
                else
                    k !== Tokenize.Tokens.ENDMARKER && throw("Unknown end of line marker for symmetry equivalent {$s}")
                    i!= 3 && throw("Input string {$s} is not a valid symmetry equivalent")
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


struct Cell
    latticesystem::Symbol
    spacegroup::String
    tablenumber::Int
    mat::SMatrix{3,3,BigFloat,9} # Cartesian coordinates of a, b and c
    equivalents::Vector{EquivalentPosition}

    function Cell(lattice, space, table, a, b, c, α, β, γ, eq)
        cosα = cosd(α); cosβ = cosd(β); cosγ = cosd(γ); sinγ = sind(γ)
        ω = sqrt(1 - cosα^2 - cosβ^2 - cosγ^2 + 2*cosα*cosβ*cosγ)
        mat = SMatrix{3,3,BigFloat,9}([a  b*cosγ  c*cosβ ;
                                      0   b*sinγ  c*(cosα - cosβ*cosγ)/sinγ ;
                                      0   0       c*ω/sinγ ])
        return new(lattice, space, table, mat, eq)
    end

    # function Cell(c::Cell, mat::StaticArray{Tuple{3,3},BigFloat})
    #     return new(c.latticesystem, c.spacegroup, c.tablenumber, mat, c.equivalents)
    # end

    function Cell(c::Cell, eqs::Vector{EquivalentPosition})
        return new(c.latticesystem, c.spacegroup, c.tablenumber, c.mat, eqs)
    end
end
function cell_parameters(mat::StaticArray{Tuple{3,3},BigFloat})
    _a, _b, _c = eachcol(mat)
    a = norm(_a)
    b = norm(_b)
    c = norm(_c)
    α = acosd(_b'_c/(b*c))
    β = acosd(_c'_a/(c*a))
    γ = acosd(_a'_b/(a*b))
    return (a, b, c, α, β, γ)
end
cell_parameters(cell::Cell) = cell_parameters(cell.mat)
function Cell(cell::Cell, mat::StaticArray{Tuple{3,3},BigFloat})
    a, b, c, α, β, γ = cell_parameters(mat)
    return Cell(cell.latticesystem, cell.spacegroup, cell.tablenumber, a, b, c, α, β, γ, cell.equivalents)
end
function Base.show(io::IO, cell::Cell)
    a, b, c, α, β, γ = Float64.(cell_parameters(cell))
    print(io, "Cell(\"$(cell.spacegroup)\", ($a, $b, $c), ($α, $β, $γ))")
end

struct CIF
    cifinfo::Dict{String, Union{String, Vector{String}}}
    cell::Cell
    ids::Vector{Int}
    types::Vector{Symbol}
    pos::Matrix{Float64}
    bonds::Matrix{Bool}
end

function keep_atoms(cif::CIF, kept)
    kept_ids = sort!([cif.ids[i] for i in kept])
    unique!(kept_ids)
    idmap = Vector{Int}(undef, length(cif.types))
    for (i,x) in enumerate(kept_ids)
        idmap[x] = i
    end
    return CIF(cif.cifinfo, cif.cell, [idmap[cif.ids[i]] for i in kept],
               cif.types[kept_ids], cif.pos[:, kept], cif.bonds[kept, kept])
end

function strip_atoms(cif::CIF, atoms)
    vmap = [i for i in cif.ids if cif.types[i] ∉ atoms]
    return CIF(cif.cifinfo, cif.cell, collect(1:length(vmap)),
               cif.types[vmap], cif.pos[:, vmap], cif.bonds[vmap, vmap])
end

# function strip_atoms(cif::CIF, atoms::Union{Tuple{Vararg{Symbol}},AbstractVector{Symbol}})
#     kept = [i for i in 1:length(cif.ids) if cif.types[import[]] ∉ atoms]
#     keep_atoms(cif, kept)
# end

function periodic_distance(u, v)
    dst = 0.0
    #=@inbounds=# for i in 1:3
        x = abs2(u[i] - v[i])
        if x > 0.25
            x = (1 - sqrt(x))^2
        end
        dst += x
    end
    return sqrt(dst)
end

function periodic_distance(u, v, mat)
    dst = 0.0
    x = similar(u)
    uu = copy(u)
    vv = copy(v)
    #=@inbounds=# for i in 1:3
        diff = u[i] - v[i]
        if diff > 0.5
            x[i] = diff - 1
        elseif - diff > 0.5
            x[i] = diff + 1
        else
            x[i] = diff
        end
    end
    return norm(mat*x)
end

function expand_symmetry(cif::CIF)
    @assert iszero(cif.bonds)
    newids = copy(cif.ids)
    newpos::Vector{Vector{Float64}} = collect(eachcol(cif.pos))
    ret = Vector{Vector{Int}}
    #=@inbounds=# for equiv in cif.cell.equivalents, i in 1:length(cif.ids)
        v = newpos[i]
        p = Vector(equiv.mat*v + equiv.ofs)
        @. p = p - floor(p)
        already_present = false
        for j in 1:length(newpos)
            if periodic_distance(newpos[j], p) < 4e-4
                already_present = true
                break
            end
        end
        if !already_present
            push!(newpos, p)
            push!(newids, cif.ids[i])
        end
    end
    return CIF(cif.cifinfo, cif.cell, newids, cif.types, reduce(hcat, newpos),
               Matrix{Bool}(undef, length(newids), length(newids)))
end

function set_unique_bond_type!(cif::CIF, bond_length, bonded_atoms::Tuple{Symbol, Symbol}, tol=0.1)
    #=@inbounds=# for i in 1:length(cif.ids)
    cif.bonds[i,i] = false
        Threads.@threads for j in i+1:length(cif.ids)
            bonded = abs2(periodic_distance(cif.pos[:,i], cif.pos[:,j], cif.cell.mat) - bond_length) <= tol
            cif.bonds[i,j] = bonded
            cif.bonds[j,i] = bonded
            if bonded && minmax(cif.types[cif.ids[i]], cif.types[cif.ids[j]]) != bonded_atoms
                throw("$(cif.types[cif.ids[i]])-$(cif.types[cif.ids[j]]) is not an $(bonded_atoms[1])-$(bonded_atoms[2]) bond")
            end
        end
    end
    nothing
end

struct Crystal
    cifinfo::Dict{String, Union{String, Vector{String}}}
    cell::Cell
    ids::Vector{Int}
    types::Vector{Symbol}
    pos::Matrix{Float64}
    graph::PeriodicGraph3D
end

function Crystal(cif::CIF)
    c1 = Ref(cif)
    if iszero(cif.bonds)
        atoms::Vector{Symbol} = unique!(sort(c1[].types))
        if atoms == [:O, :Si] || atoms == [:Si] # zeolite
            c1[] = expand_symmetry(strip_atoms(cif, (:O,)))
            set_unique_bond_type!(c1[], 3.1, (:Si, :Si))
        elseif atoms == [:C] || atoms == [:C, :H]
            c1[] = expand_symmetry(strip_atoms(cif, (:H,)))
            set_unique_bond_type!(c1[], 1.54, (:C, :C), 0.3)
        elseif atoms == [:Cl, :Na]
            c1[] = expand_symmetry(cif)
            set_unique_bond_type!(c1[], 2.82, (:Cl, :Na), 0.1)
        else
            @show atoms
            throw("Missing bonds on CIF object $(cif.cifinfo["data"])")
        end
    end
    c2 = c1[]
    @assert !iszero(c2.bonds)
    bondedatoms::Vector{Int} = [i for i in 1:length(c2.ids) if count(c2.bonds[:,i]) > 1]
    c = keep_atoms(c2, bondedatoms) # Only keep those that matter for the topology
    n = length(c.ids)
    edges = PeriodicEdge3D[]
    for i in 1:n, k in findall(@view c.bonds[:,i])
        k < i && continue
        offset::Vector{SVector{3, Int}} = []
        old_dst = norm(c.cell.mat*[1, 1, 1])
        for ofsx in -1:1, ofsy in -1:1, ofsz in -1:1
            # dst = norm(c.cell.mat*(c.pos[:,i] .- (c.pos[:,k] + [ofsx; ofsy; ofsz])))
            dst = norm(c.cell.mat*(c.pos[:,i] .- (c.pos[:,k] .+ [ofsx; ofsy; ofsz])))
            if abs2(dst - old_dst) < 1e-4
                push!(offset, (ofsx, ofsy, ofsz))
                old_dst = (dst + old_dst)/2
            elseif dst < old_dst
                offset = [(ofsx, ofsy, ofsz)]
                old_dst = dst
            end
        end
        for ofs in offset
            push!(edges, (i, k, ofs))
        end
    end
    @assert !isempty(edges)
    Crystal(c.cifinfo, c.cell, c.ids, c.types, c.pos, PeriodicGraph3D(n, edges))
end


function equilibrium(g::PeriodicGraph3D)
    n = nv(g)
    iszero(n) && return Matrix{Rational{Int}}(undef, 3, 0)
    Y = Array{Rational{BigInt}}(undef, n, 3)
    A = spzeros(Int, n, n)
    neigh = Array{Int}(undef, n)
    offset = SizedVector{3,Int}(undef)
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

    B = rational_lu(A[2:end,2:end], false)
    if !issuccess(B)
        throw(is_connected(g) ? "Singular exception while equilibrating" :
                                "Cannot equilibrate a disconnected graph")
    end
    Y = linsolve!(B, Y[2:end,:])
    ret = hcat(zeros(Rational{Int128}, 3), Rational{Int128}.(Y)')
    # Rational{Int64} is not enough for tep for instance.
    return ret
end

struct CrystalNet{T<:Real}
    cell::Cell
    types::Vector{Symbol}
    pos::Vector{SVector{3,T}}
    graph::PeriodicGraph3D
end

function CrystalNet{T}(cell::Cell, types::AbstractVector{Symbol}, ids::AbstractVector{<:Integer},
                    graph::PeriodicGraph3D, eq::AbstractMatrix{T}) where T
    n = nv(graph)
    pos = Vector{SVector{3,T}}(undef, n)
    offsets = Vector{SVector{3,Int}}(undef, n)
    for (i, x) in enumerate(eachcol(eq))
        offsets[i] = floor.(Int, x)
        pos[i] = x .- offsets[i]
    end
    s = sortperm(pos)
    pos = pos[s]
    types = Symbol[types[ids[s[i]]] for i in 1:n]
    cell = Cell(cell, EquivalentPosition[])
    graph = offset_representatives!(graph, .-offsets)[s]
    @assert all(pos[i] == mean(pos[x.v] .+ x.ofs for x in neighbors(graph, i)) for i in 1:length(pos))
    return CrystalNet(cell, types, pos, graph)
end

function CrystalNet(cell::Cell, types::AbstractVector{Symbol}, ids::AbstractVector{<:Integer},
                    graph::PeriodicGraph3D, eq::Matrix{T}) where T
    CrystalNet{T}(cell, types, ids, graph, eq)
end

function CrystalNet(cell::Cell, types::AbstractVector{Symbol}, ids::AbstractVector{Int}, graph::PeriodicGraph3D)
    eq = equilibrium(graph)
    minnum = minimum(numerator.(eq))
    maxnum = maximum(numerator.(eq))
    minden = minimum(denominator.(eq))
    maxden = maximum(denominator.(eq))
    for T in (Int32, Int64)
        if ((typemin(T) < min(minnum, minden)) & (max(maxnum, maxden) < typemax(T)))
            return CrystalNet{Rational{T}}(cell, types, ids, graph, Rational{T}.(eq))
        end
    end
    return CrystalNet{Rational{Int128}}(cell, types, ids, graph, eq)
    # Type-unstable function, but yields better performance than always falling back on Int128
end
function CrystalNet{T}(cell::Cell, types::AbstractVector{Symbol}, ids::AbstractVector{Int}, graph::PeriodicGraph3D) where T
    CrystalNet{T}(cell, types, ids, graph, T.(equilibrium(graph)))
end

CrystalNet(c::Crystal) = CrystalNet(c.cell, c.types, c.ids, c.graph)

function CrystalNet(cell::Cell, types::AbstractVector{Symbol}, graph::PeriodicGraph3D)
    n = nv(graph)
    return CrystalNet(cell, types, 1:n, graph)
end

function CrystalNet(g::Union{PeriodicGraph3D,AbstractString,AbstractVector{PeriodicEdge3D}})
    graph = PeriodicGraph3D(g)
    cell = Cell(Symbol(""), "P 1", 0, 10, 10, 10, 90, 90, 90, EquivalentPosition[])
    n = nv(g)
    types = [Symbol("") for _ in 1:n]
    return CrystalNet(cell, types, graph)
end


include("symmetries.jl")

end
