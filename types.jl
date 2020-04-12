module CIFTypes

include("./PeriodicGraphs.jl")
using Serialization
using StaticArrays, Tokenize
using .PeriodicGraphs
import LinearAlgebra: norm
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
        cosα = cosd(α); cosβ = cosd(β); cosγ = cosd(γ); sinβ = sind(β)
        ω = sqrt(1 - cosα^2 - cosβ^2 - cosγ^2 + 2*cosα*cosβ*cosγ)
        mat = SMatrix{3,3,BigFloat,9}([a   b*cosγ                      c*cosβ ;
                                      0   b*ω/sinβ                    0      ;
                                      0   b*(cosα - cosβ*cosγ)/sinβ   c*sinβ ])
        return new(lattice, space, table, mat, eq)
    end

    function Cell(c::Cell, mat::SMatrix{3,3,BigFloat,9})
        return new(c.latticesystem, c.spacegroup, c.tablenumber, mat, c.equivalents)
    end

    function Cell(c::Cell, eqs::Vector{EquivalentPosition})
        return new(c.latticesystem, c.spacegroup, c.tablenumber, c.mat, eqs)
    end
end
function cell_parameters(cell::Cell)
    a, b, c = eachcol(cell.mat)# ./ scale_factor)
    α = Float64(acosd(b'c/(norm(b)*norm(c))))
    β = Float64(acosd(c'a/(norm(c)*norm(a))))
    γ = Float64(acosd(a'b/(norm(a)*norm(b))))
    return (Float64(norm(a)), Float64(norm(b)), Float64(norm(c)), α, β, γ)
end
function Base.show(io::IO, cell::Cell)
    a, b, c, α, β, γ = cell_parameters(cell)
    print(io, "Cell(\"$(cell.spacegroup)\", ($a, $b, $c), ($α, $β, $γ)")
end

struct CIF
    cifinfo::Dict{String, Union{String, Vector{String}}}
    cell::Cell
    ids::Vector{Int}
    types::Vector{Symbol}
    pos::Matrix{Float64}
    bonds::BitMatrix
end

function strip_atoms(cif::CIF, atoms)
    vmap = [i for i in cif.ids if cif.types[i] ∉ atoms]
    return CIF(cif.cifinfo, cif.cell, collect(1:length(vmap)),
               cif.types[vmap], cif.pos[:, vmap], cif.bonds[vmap, vmap])
end

function periodic_distance(u, v)
    dst = 0.0
    @inbounds for i in 1:3
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
    @inbounds for i in 1:3
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
    @inbounds for equiv in cif.cell.equivalents, i in 1:length(cif.ids)
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
               BitMatrix(undef, length(newids), length(newids)))
end

function set_unique_bond_type!(cif::CIF, bond_length, bonded_atoms::Tuple{Symbol, Symbol}, tol=0.1)
    @inbounds for i in 1:length(cif.ids)
        @simd for j in 1:length(cif.ids)
            bonded = abs2(periodic_distance(cif.pos[:,i], cif.pos[:,j], cif.cell.mat) - bond_length) <= tol
            cif.bonds[i,j] = bonded
            if bonded && (cif.types[cif.ids[i]], cif.types[cif.ids[j]]) != bonded_atoms
                throw("Not an $(bonded_atoms[1])-$(bonded_atoms[2]) bond")
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

function Crystal(c::CIF)
    if iszero(c.bonds)
        atoms = unique!(sort(c.types))
        if atoms == [:O, :Si] || atoms == [:Si] # zeolite
            c = expand_symmetry(strip_atoms(c, (:O,)))
            set_unique_bond_type!(c, 3.1, (:Si, :Si))
        elseif atoms == [:C]
            c = expand_symmetry(c)
            set_unique_bond_type!(c, 1.54, (:C, :C), 0.3)
        else
            throw("Missing bonds on CIF object $(c.cifinfo["data"])")
        end
    end
    edges = PeriodicEdge3D[]
    n = length(c.ids)
    for i in 1:n, k in findall(@view c.bonds[i,:])
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
    Crystal(c.cifinfo, c.cell, c.ids, c.types, c.pos, PeriodicGraph3D(n, edges))
end

struct CrystalNet
    cell::Cell
    types::Vector{Symbol}
    pos::Matrix{Rational{Int}}
    graph::PeriodicGraph3D
end


include("symmetries.jl")

end
