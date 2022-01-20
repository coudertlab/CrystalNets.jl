# Wrapper around the spglib library

include("symmetry_groups.jl")

import spglib_jll: libsymspg

mutable struct SpglibDataset
    spacegroup_number::Cint
    hall_number::Cint
    _international_symbol::NTuple{11,Cchar}
    _hall_symbol::NTuple{17,Cchar}
    _choice::NTuple{6,Cchar}
    _transformation_matrix::NTuple{9,Cdouble}
    _origin_shift::NTuple{3,Cdouble}
    n_operations::Cint
    _rotations::Ptr{Cint}
    _translations::Ptr{Cdouble}
    n_atoms::Cint
    _wyckoffs::Ptr{Cint}
    _site_symmetry_symbols::Ptr{NTuple{7,Cchar}}
    _equivalent_atoms::Ptr{Cint}
    _crystallographic_orbits::Ptr{Cint}
    _primitive_lattice::NTuple{9,Cdouble}
    _mapping_to_primitive::Ptr{Cint}
    n_std_atoms::Cint
    _std_lattice::NTuple{9,Cdouble}
    _std_types::Ptr{Cint}
    _std_positions::Ptr{Cdouble}
    _std_rotation_matrix::NTuple{9,Cdouble}
    _std_mapping_to_primitive::Ptr{Cint}
    _pointgroup_symbol::NTuple{6,Cchar}
end
function Base.getproperty(ds::SpglibDataset, name::Symbol)
    if name === :international_symbol
        letters = Char.(getfield(ds, :_international_symbol))
        fzero = findfirst(==('\0'), collect(letters))
        return join(letters[i] for i in 1:(isnothing(fzero) ? length(letters) : (fzero-1)))
    elseif name === :hall_symbol
        letters = Char.(getfield(ds, :_hall_symbol))
        fzero = findfirst(==('\0'), collect(letters))
        return join(letters[i] for i in 1:(isnothing(fzero) ? length(letters) : (fzero-1)))
    elseif name === :choice
        letters = Char.(getfield(ds, :_choice))
        fzero = findfirst(==('\0'), collect(letters))
        return join(letters[i] for i in 1:(isnothing(fzero) ? length(letters) : (fzero-1)))
    elseif name === :transformation_matrix
        return reshape(collect(getfield(ds, :_transformation_matrix)), (3,3))
    elseif name === :origin_shift
        return reshape(collect(getfield(ds, :_origin_shift)), 3)
    elseif name === :rotations
        return unsafe_wrap(Array, getfield(ds, :_rotations), (3,3,Int(getfield(ds, :n_operations))))
    elseif name === :translations
        return unsafe_wrap(Array, getfield(ds, :_translations), (3,Int(getfield(ds, :n_operations))))
    elseif name === :wyckoffs
        return unsafe_wrap(Array, getfield(ds, :_wyckoffs), Int(getfield(ds, :n_atoms)))
    elseif name === :site_symmetry_symbols
        a = unsafe_wrap(Array, getfield(ds, :_site_symmetry_symbols), getfield(ds, :n_atoms))
        ret = Vector{String}(undef, length(a))
        for i in 1:length(a)
            letters = Char.(a[i])
            fzero = findfirst(==('\0'), collect(letters))
            ret[i] = join(letters[j] for j in 1:(isnothing(fzero) ? length(letters) : (fzero-1)))
        end
        return ret
    elseif name === :equivalent_atoms
        return unsafe_wrap(Array, getfield(ds, :_equivalent_atoms), getfield(ds, :n_atoms))
    elseif name === :crystallographic_orbits
        return unsafe_wrap(Array, getfield(ds, :_crystallographic_orbits), getfield(ds, :n_atoms))
    elseif name === :primitive_lattice
        return reshape(collect(getfield(ds, :_primitive_lattice)), (3,3))
    elseif name === :mapping_to_primitive
        return unsafe_wrap(Array, getfield(ds, :_mapping_to_primitive), getfield(ds, :n_atoms))
    elseif name === :std_lattice
        return reshape(collect(getfield(ds, :_std_lattice)), (3,3))
    elseif name === :std_types
        return unsafe_wrap(Array, getfield(ds, :_std_types), getfield(ds, :n_std_atoms))
    elseif name === :std_positions
        return unsafe_wrap(Array, getfield(ds, :_std_positions), (3,Int(getfield(ds, :n_std_atoms))))
    elseif name === :std_rotation_matrix
        return reshape(collect(getfield(ds, :_std_rotation_matrix)), (3,3))
    elseif name === :std_mapping_to_primitive
        return unsafe_wrap(Array, getfield(ds, :_std_mapping_to_primitive), getfield(ds, :n_std_atoms))
    elseif name === :pointgroup_symbol
        letters = Char.(getfield(ds, :_pointgroup_symbol))
        fzero = findfirst(==('\0'), collect(letters))
        return join(letters[i] for i in 1:(isnothing(fzero) ? length(letters) : (fzero-1)))
    else
        getfield(ds, name)
    end
end


function find_hall_number(hallsymbol, hm, it)
    if hallsymbol != ""
        hsymbol = lowercase(hallsymbol)
        hall = get(SPACE_GROUP_HALL, replace(hsymbol, ('_' => ' ')), 0)
        if hall ∈ (322, 326, 330)
            if hall == 322
                if endswith(hm, 'b')
                    hall = 324
                end
            elseif hall == 326
                if occursin('c', hm)
                    hall = 328
                end
            elseif occursin('a', hm)
                hall = 332
            end
        end
        hall == 0 || return hall
        @ifwarn @warn "Hall symbol provided but not recognised: $hallsymbol"
    end
    if hm != ""
        hm = hm[1]*lowercase(@view hm[2:end])
        dense_hm = replace(join(split(hm)), ('_'=>""))
        hall = get(SPACE_GROUP_HM, dense_hm, 0)
        hall == 0 || return hall
        hall = get(SPACE_GROUP_FULL, dense_hm, 0)
        hall == 0 || return hall
        @ifwarn begin
            @warn "H-M symbol provided but not recognised : $hm"
            keysHM = collect(keys(SPACE_GROUP_HM))
            suggestions = keysHM[findall(startswith(dense_hm), keysHM)]
            keysFULL = collect(keys(SPACE_GROUP_FULL))
            append!(suggestions, keysFULL[findall(startswith(dense_hm), keysFULL)])
            if !isempty(suggestions)
                msg = join(suggestions, ", ", " or ")
                @info "Did you perhaps mean to use $msg?"
            end
        end
    end
    if it != 0
        if it < 1 || it > 230
            @ifwarn @error "International Table number provided is outside of the allowed range (1-230) : $it"
            return 1
        end
        return SPACE_GROUP_IT[it]
    end
    @ifwarn @error "Could not determine the space group of this file. It will be considered P1."
    return 1
end

function get_symmetry_equivalents(hall)
    rotations = Array{Cint}(undef, 3, 3, 192)
    translations = Array{Cdouble}(undef, 3, 192)
    len = ccall((:spg_get_symmetry_from_database, libsymspg), Cint,
                (Ptr{Cint}, Ptr{Cdouble}, Cint), rotations, translations, hall)
    @toggleassert 0 < len ≤ 192
    eqs = EquivalentPosition[]
    for i in 1:len
        rot = SMatrix{3,3,Rational{Int},9}(transpose(@inbounds rotations[:,:,i]))
        tr = SVector{3,Cdouble}(@inbounds translations[:,i])
        trans = SVector{3,Rational{Int}}(round.(Int, 360 .* tr) .// 360)
        push!(eqs, EquivalentPosition(rot, trans))
    end
    return eqs
end


function get_spglib_dataset(net::CrystalNet3D)
    lattice = Matrix{Cdouble}(adjoint(net.cell.mat)) # transpose to account for row-major operations
    n = nv(net.graph)
    positions = Matrix{Cdouble}(undef, 3, n)
    types = Vector{Cint}(undef, n)
    symb_to_int = Dict{Symbol,Cint}()
    j = 1
    for i in 1:n
        positions[:,i] .= net.cell.mat * net.pos[i]
        j += ((types[i] = get!(symb_to_int, net.types[i], j)) == j)
    end
    ptr = ccall((:spg_get_dataset, libsymspg), Ptr{SpglibDataset},
                (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
                lattice, positions, types, n, 10*eps(Cdouble))
    dataset = unsafe_load(ptr)
    @toggleassert dataset.n_atoms == n # otherwise the net is not minimal
    return dataset
end


function find_symmetries(net::CrystalNet3D{Rational{S}}, collisions) where S
    T = soft_widen(S)
    U = soft_widen(T)
    lattice = Matrix{Cdouble}(LinearAlgebra.I, 3, 3) # positions are expressed in this basis

    I = sortperm(net.pos)
    uniquepos = SVector{3, Rational{S}}[net.pos[I[1]]]
    symb_to_int = Dict{Vector{Tuple{Symbol,Int}}, Cint}([(net.types[I[1]], 1)] => 1)
    last_types = [(net.types[I[1]], 1)]
    types = Cint[1]
    new_symbol = 2
    for i in 2:length(net.pos)
        j = I[i]
        if net.pos[j] == last(uniquepos)
            flag = true
            for k in 1:length(last_types)
                thistype, x = last_types[k]
                if thistype == net.types[j]
                    last_types[k] = (thistype, x+1)
                    flag = false
                    break
                end
            end
            if flag
                push!(last_types, (net.types[j], 1))
            end
            symb = get!(symb_to_int, sort!(last_types), new_symbol)
            types[end] = symb
            new_symbol += (symb == new_symbol)
        else
            push!(uniquepos, net.pos[j])
            last_types = [(net.types[j], 1)]
            symb = get!(symb_to_int, sort!(last_types), new_symbol)
            push!(types, symb)
            new_symbol += (symb == new_symbol)
        end
    end

    n = length(uniquepos)
    positions = Matrix{Cdouble}(undef, 3, n)
    den = 1
    for i in 1:n
        den = lcm(den, lcm((U∘denominator).(uniquepos[i])))
        positions[:,i] .= uniquepos[i]
    end

    len = ccall((:spg_get_multiplicity, libsymspg), Cint,
                (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
                lattice, positions, types, n, 100*eps(Cdouble))
    @toggleassert len >= 1
    rotations = Array{Cint}(undef, 3, 3, len)
    translations = Array{Cdouble}(undef, 3, len)
    _len = ccall((:spg_get_symmetry, libsymspg), Cint,
                 (Ptr{Cint}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
                 rotations, translations, len, lattice, positions, types, n, 100*eps(Cdouble))
    @toggleassert _len == len
    # @show len
    den = lcm(len, den)
    # The first symmetry should always be the identity, which we can skip
    @toggleassert isone(rotations[:,:,1])
    @toggleassert iszero(translations[:,1])
    symmetries = SMatrix{3,3,Int,9}[rotations[:,:,1]]
    vmaps = Vector{Int}[collect(1:nv(net.graph))]
    floatpos = [float(x) for x in net.pos]
    hasmirror = false # whether a mirror plane exists or not. If so, the graph is achiral
    for i in 2:len
        rot = SMatrix{3,3,T,9}(transpose(rotations[:,:,i]))
        if det(rot) < 0
            @toggleassert det(rot) == -1
            hasmirror = true
            # continue # No need to keep it as its mirror image will be kept
        end
        # @toggleassert det(rot) == 1
        tr = SVector{3,Cdouble}(translations[:,i])
        trans = SVector{3,Rational{S}}(round.(U, den .* tr) .// den)
        vmap = check_valid_symmetry(net, trans, collisions, rot)
        if isnothing(vmap)
            trans = SVector{3,Rational{S}}(net.pos[last(findmin([norm(x .- tr) for x in floatpos]))])
            vmap = check_valid_symmetry(net, trans, collisions, rot)
            # if isnothing(vmap)
            #     @show translations[:,i]
            #     @show float(trans)
            #     @show trans
            #     @show rot
            #     println("net = CrystalNet(PeriodicGraph3D(\"", net.graph, "\"))")
            #     throw("Invalid symmetry") # This error is for debugging only, false positives without it being an error
            # end
            isnothing(vmap) && continue
        end
        @toggleassert rot != LinearAlgebra.I # Otherwise the symmetry is a pure translation, so the net is not minimal
        push!(symmetries, rot)
        push!(vmaps, vmap::Vector{Int})
    end
    # @show vmaps
    # @show length(symmetries)
    return (symmetries, vmaps, hasmirror)
end
