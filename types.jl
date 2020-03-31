module CIFTypes

include("./PeriodicGraphs.jl")
using .PeriodicGraphs
using StaticArrays, Tokenize
import LinearAlgebra: norm, I
export Cell, CIF, Crystal, EquivalentPosition

struct EquivalentPosition
    mat::SMatrix{3, 3, Rational{Int}, 9}
    ofs::SVector{3, Rational{Int}}
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
                val = isnothing(curr_val)  ? Rational{Int}(1) : val
                j = const_dict[Tokenize.Tokens.untokenize(x)]
                mat[i,j] += sign * val
                curr_val = nothing
                curr_sign = nothing
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
                        @warn "Existing offset already existing for in \"$s\""
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
                    @assert 1 <= i <= 2
                    i += 1
                else
                    @assert k === Tokenize.Tokens.ENDMARKER
                    # @assert i == 3
                end
            end
        end
    end
    EquivalentPosition(SMatrix{3, 3, Rational{Int}, 9}(mat), SVector{3, Rational{Int}}(ofs))
end

function Base.show(io::IO, eq::EquivalentPosition)
    function rationaltostring(x::Rational{<:Integer})
        sign = x > 0 ? '+' : '-'
        (x == 1 || x == -1) && return string(sign)
        sign * (x.den == 1 ? string(x.num) : string(x.num)*'/'*string(x.den))
    end
    xyz = ('x', 'y', 'z')
    for i in 1:3
        eq.ofs[i] == 0 || print(io, rationaltostring(eq.ofs[i]))
        for j in 1:3
            if eq.mat[i,j] != 0
                print(io, rationaltostring(eq.mat[i,j]))
                print(io, xyz[j])
            end
        end
        i < 3 && print(io, ',')
    end
    nothing
end

struct Cell
    latticesystem::Symbol
    spacegroup::String
    tablenumber::Int
    a::Float64
    b::Float64
    c::Float64
    α::Float64
    β::Float64
    γ::Float64
    equivalents::Vector{EquivalentPosition}
    mat::SMatrix{3, 3, Float64, 9} # Conversion between fractional and cartesian coordinates

    function Cell(lattice, space, table, a, b, c, α, β, γ, eq)
        cosα = cosd(α); cosβ = cosd(β); cosγ = cosd(γ); sinβ = sind(β)
        ω = sqrt(1 - cosα^2 - cosβ^2 - cosγ^2 + 2*cosα*cosβ*cosγ)
        mat = SMatrix{3,3,Float64,9}([a   b*cosγ                      c*cosβ ;
                                      0   b*ω/sinβ                    0      ;
                                      0   b*(cosα - cosβ*cosγ)/sinβ   c*sinβ ])
        return new(lattice, space, table, a, b, c, α, β, γ, eq, mat)
    end
end
function show(io::IO, c::Cell)
    print(io, "Cell($(c.spacegroup), ($(c.a), $(c.b), $(c.c)), ($(c.α), $(c.β), $(c.γ))")
end

struct CIF
    natoms::Int
    cifinfo::Dict{String, Union{String, Vector{String}}}
    geometry::Cell
    atoms::Vector{Symbol}
    pos::Matrix{Float64}
    bonds::BitMatrix
end

function strip_atoms(cif::CIF, atoms)
    tokeep = cif.atoms .∉ Ref(atoms)
    natoms = count(tokeep)
    return CIF(natoms, cif.cifinfo, cif.geometry, cif.atoms[tokeep],
               cif.pos[:, tokeep], cif.bonds[tokeep, tokeep])
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
    newatoms = copy(cif.atoms)
    newpos::Vector{Vector{Float64}} = collect(eachcol(cif.pos))
    ret = Vector{Vector{Int}}
    @inbounds for equiv in cif.geometry.equivalents, i in 1:cif.natoms
        v = newpos[i]
        p = Vector(equiv.mat*v + equiv.ofs)
        @. p = p - floor(p)
        already_present = false
        for j in 1:length(newpos)
            if periodic_distance(newpos[j], p) < 1e-4
                already_present = true
                break
            end
        end
        if !already_present
            push!(newpos, p)
            push!(newatoms, cif.atoms[i])
        end
    end
    natoms = length(newatoms)
    return CIF(natoms, cif.cifinfo, cif.geometry, newatoms, reduce(hcat, newpos),
               BitMatrix(undef, natoms, natoms))
end

function set_unique_bond_type!(cif::CIF, bond_length, bonded_atoms::Tuple{Symbol, Symbol}, tol=0.1)
    @inbounds for i in 1:cif.natoms
        @simd for j in 1:cif.natoms
            bonded = abs2(periodic_distance(cif.pos[:,i], cif.pos[:,j], cif.geometry.mat) - bond_length) <= tol
            cif.bonds[i,j] = bonded
            if bonded && minmax(cif.atoms[i], cif.atoms[j]) != bonded_atoms
                throw("Not an $(bonded_atoms[1])-$(bonded_atoms[2]) bond")
            end
        end
    end
    nothing
end

struct Crystal
    natoms::Int
    cifinfo::Dict{String, Union{String, Vector{String}}}
    geometry::Cell
    atoms::Vector{Symbol}
    pos::Matrix{Float64}
    graph::PeriodicGraph3D
end

function Crystal(c::CIF)
    if iszero(c.bonds)
        if unique!(sort(c.atoms)) == [:O, :Si] # zeolite
            c = expand_symmetry(strip_atoms(c, (:O,)))
            set_unique_bond_type!(c, 3.1, (:Si, :Si))
        else
            throw("Missing bonds on CIF object $(c.name)")
        end
    end
    edges = PeriodicEdge3D[]
    n = c.natoms
    for i in 1:n, k in findall(@view c.bonds[i,:])
        k < i && continue
        offset::Vector{Tuple{Int, Int, Int}} = []
        old_dst = norm(c.geometry.mat*[1, 1, 1])
        for ofsx in -1:1, ofsy in -1:1, ofsz in -1:1
            # dst = norm(c.geometry.mat*(c.pos[:,i] .- (c.pos[:,k] + [ofsx; ofsy; ofsz])))
            dst = norm(c.geometry.mat*(c.pos[:,i] .- (c.pos[:,k] .+ [ofsx; ofsy; ofsz])))
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
    Crystal(n, c.cifinfo, c.geometry, c.atoms, c.pos, PeriodicGraph3D(n, edges))
end


end
