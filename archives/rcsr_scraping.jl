using LinearAlgebra
using HTTP, Gumbo, Cascadia, StaticArrays, PeriodicGraphs, Graphs, CrystalNets
using Base.Threads
CrystalNets.toggle_warning(false)
CrystalNets.toggle_export(false)

import CrystalNets: periodic_distance, expand_symmetry, Cell, CIF

function getrcsr(n)
    rcsr = String(HTTP.get("http://rcsr.anu.edu.au/data/$(n)dall.txt").body)
    @assert startswith(rcsr, "start")
    return split(rcsr, r"(\r\n|\n)")
end

const rcsr3D = getrcsr(3)

function precise_round(x, n)
    u = 10^(n+1)
    y = floor(Int, u*x)
    r = rem(y+5, 10)
    ret = (y - r + ifelse(r == 0, 0, 5)) / u
    ret == 1.0 && return precise_round(x, n+1)
    ret
end

"""
Given a position `u`, return a vector of offsets `ofs` such that `u .+ ofs` are the
periodic images of `u` closest to the origin.
"""
function periodic_neighbor!(buffer, u, mat, ortho, safemin, ε)
    ofs = MVector{3,Int}(undef)
    @simd for i in 1:3
        diff = u[i] + 0.5
        ofs[i] = -floor(Int, diff)
        buffer[i] = diff + ofs[i] - 0.5
    end

    ref = norm(mat*buffer)
    ofss = SVector{3,Int}[ofs]
    ref ≤ safemin && return ofss, ref

    totry = if ortho
        _totry = SVector{3,Int}[]
        for i in 1:3
            if isapprox(buffer[i], -0.5; rtol=0.01)
                __totry = MVector{3,Int}(ofs)
                __totry[i] += 1
                push!(_totry, __totry)
                __totry[i] -= 2
                push!(_totry, __totry)
            end
        end
        _totry
    else
        [ofs .+ x for x in SVector{3,Int}[(1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1)]]
    end
    # totry = [ofs .+ x for x in SVector{3,Int}[(1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1)]]

    for newofs in totry
        buffer .= u .+ newofs
        newnorm = norm(mat*buffer)
        if newnorm < ref - ε
            ofss[1] = newofs
            resize!(ofss, 1)
            ref = newnorm
        elseif newnorm ≤ ref + ε
            push!(ofss, newofs)
        end
    end
    return ofss, ref
end

"""
Check whether the given graph has the correct list of coordination sequences
"""
function check_graph(g, seqs)
    n = length(seqs)
    nv(g) != n && return false
    degree(g) != first.(seqs) && return false
    coord_flag = false
    for (i, seq) in enumerate(seqs)
        if coordination_sequence(g, i, 10) != seq
            coord_flag = true
            break
        end
    end
    coord_flag && return false
    return true
end

"""
Return three lists: ids, dists and closest_pos, where for each index i, i is the index of
a potential closest neighbour of x such that:
- ids[i] is the index of the corresponding vertex in list poss
- dists[i] is the distance between that vertex and x
- closest_pos[i] is the position of the vertex
The number of closest neighbours returned is at most maxsize.
"""
function closest_positions(x, mat, poss, ε, maxsize)
    _, ortho, _safemin = CrystalNets.prepare_periodic_distance_computations(mat)
    safemin = 6*_safemin/7
    buffer = MVector{3,Float32}(undef)
    closest_pos = SVector{3,Float32}[]
    ids = Int[]
    dists = Float32[]
    for (i, y) in enumerate(poss)
        ofss, ref = periodic_neighbor!(buffer, x .- y, mat, ortho, safemin, ε)
        append!(ids, i for _ in ofss)
        append!(dists, ref for _ in ofss)
        append!(closest_pos, [y .+ ofs for ofs in ofss])
    end
    sizeI = min(maxsize, length(dists))
    I = sizeI == length(dists) ? sortperm(dists) : collect(partialsortperm(dists, 1:sizeI))
    if ε != Inf
        d1 = dists[I[1]]
        for i in 2:sizeI
            if dists[I[i]] > d1 + ε
                resize!(I,  i-1)
                break
            end
        end
    end
    return ids[I], dists[I], closest_pos[I]
end

function nearest_neighbors(x, poss, dpos, mat, ε, round_i, maxsize=6)
    ids, dists, closest_pos = closest_positions(x, mat, poss, ε, maxsize)
    ret = PeriodicEdge3D[]
    for (d, pos) in zip(dists, closest_pos)
        otherpos = x .- (pos .- x) # so that x is the midpoint of pos and otherpos
        ofspos = floor.(Int, pos)
        ofsother = floor.(Int, otherpos)
        src = get(dpos, round.(pos .- ofspos; digits=round_i), nothing)
        src === nothing && continue
        dst = get(dpos, round.(otherpos .- ofsother; digits=round_i), nothing)
        dst === nothing && continue
        edg = PeriodicEdge3D(src, dst, ofsother .- ofspos)
        if PeriodicGraphs.isindirectedge(edg)
            edg = reverse(edg)
        end
        push!(ret, edg)
    end
    sort!(ret); unique!(ret)
    return ret
end

function try_from_edges(_cife, mat, round_i, pos, dpos, ε, seqs)
    cife = expand_symmetry(_cife)
    pose = [SVector{3,Float32}(round.(x; digits=round_i+1)) for x in eachcol(cife.pos)]
    @assert allunique(pose)
    nns = [nearest_neighbors(x, pos, dpos, mat, ε, round_i) for x in pose]
    progression = [length(nn) for nn in nns]
    progression[end] = 0
    progression_made = true
    edgs = PeriodicEdge3D[]
    counter = 0
    while progression_made
        progression_made = false
        empty!(edgs)
        counter ≥ 8192 && break
        counter += 1
        for (i, nn) in enumerate(nns)
            if !progression_made
                if progression[i] == length(nn)
                    progression[i] = isempty(nn) ? 0 : 1
                else
                    progression[i] += 1
                    progression_made = true
                end
            end
            isempty(nn) && continue
            push!(edgs, nn[progression[i]])
        end
        g = PeriodicGraph(edgs)
        check_graph(g, seqs) && return g
    end
    return nothing
end

function try_from_closest(mat, poss, seqs)
    edgs = PeriodicEdge3D[]
    for (i, pos) in enumerate(poss)
        idxs, _, closest_pos = closest_positions(pos, mat, poss, Inf, seqs[i][1] + 1)
        length(closest_pos) == seqs[i][1] + 1 || return nothing
        @assert idxs[1] == i && iszero(floor.(Int, closest_pos[1]))
        popfirst!(idxs); popfirst!(closest_pos)
        append!(edgs, PeriodicEdge3D(i, idx, floor.(Int, pos)) for (idx, pos) in zip(idxs, closest_pos))
    end
    g = PeriodicGraph(edgs)
    check_graph(g, seqs) && return g
    nothing
end


function try_from_closest_OLD(_cifv::CIF, cifv::CIF, mat, round_i, poss, dpos, ε, seqs)
    n = length(_cifv.ids)
    m = length(cifv.ids)
    num_symm = length(cifv.cell.equivalents)
    @assert num_symm == length(_cifv.cell.equivalents)

    symmetries = Vector{Tuple{Vector{PeriodicVertex3D},SMatrix{3,3,Int,9}}}(undef, num_symm)
    posbuffer = MVector{3,Float32}(undef)
    ofsbuffer = MVector{3,Int}(undef)

    equivalents = cifv.cell.equivalents
    #=@inbounds=#for (j, equiv) in enumerate(equivalents)
        newpos = Vector{PeriodicVertex3D}(undef, m)
        for i in 1:m
            posbuffer .= equiv.mat * cifv.pos[:,i] .+ equiv.ofs
            ofsbuffer .= floor.(Int, posbuffer)
            v = dpos[round.(posbuffer .- ofsbuffer; digits=round_i)]
            newpos[i] = PeriodicVertex3D(v, ofsbuffer)
        end
        symmetries[j] = (newpos, equiv.mat)
    end
    I = sortperm(symmetries; by=x->(x[1], reshape(x[2], 9)))
    todelete = Int[]
    for i in 2:length(I)
        if symmetries[I[i]] == symmetries[I[i-1]]
            push!(todelete, i)
        end
    end
    if !isempty(todelete)
        deleteat!(I, todelete)
        equivalents = cifv.cell.equivalents[I]
        symmetries = symmetries[I]
    end

    closest_neighbours = Vector{Vector{PeriodicVertex3D}}(undef, n)
    for (i, pos) in enumerate(eachcol(_cifv.pos))
        idxs, _, closest_pos = closest_positions(pos, mat, poss, ε, seqs[i][1] + 3)
        @assert length(closest_pos) ≥ seqs[i][1] + 1
        closest_neighbours[i] = [PeriodicVertex3D(idx, floor.(Int, pos)) for (idx, pos) in zip(idxs, closest_pos)]
        @assert closest_neighbours[i][1] == PeriodicVertex3D(i)
        popfirst!(closest_neighbours[i])
    end

    orbits = [[symm[1][i] for symm in symmetries] for i in 1:m]
    rots = [symm[2] for symm in symmetries]
    progression = ones(Int, n)
    progression[1] = 0
    progression_made = true
    edgs = PeriodicEdge3D[]
    g = PeriodicGraph3D()
    while progression_made
        progression_made = false
        skip_check_g = false
        empty!(edgs)
        for (i, closest) in enumerate(closest_neighbours)
            progression_made_now = !progression_made
            if progression_made_now
                if progression[i] == length(closest)
                    progression[i] = 1
                    progression_made_now = false
                    skip_check_g = true
                else
                    progression[i] += 1
                    progression_made = true
                end
            end
            neigh = closest[progression[i]]
            if progression_made_now && !skip_check_g
                initial_level = progression[i] - 1
                while has_edge(g, PeriodicEdge3D(i, neigh))
                    if progression[i] == length(closest)
                        progression[i] = initial_level
                        neigh = closest[progression[i]]
                        progression_made = false
                        break
                    end
                    progression[i] += 1
                    neigh = closest[progression[i]]
                end
            end

            push!(edgs, PeriodicEdge3D(i, neigh))
            for (ui, uneigh, rot) in zip(orbits[i], orbits[neigh.v], rots)
                push!(edgs, PeriodicEdge3D(ui.v, uneigh.v, rot*neigh.ofs .+ uneigh.ofs .- ui.ofs))
            end
        end
        g = PeriodicGraph(edgs)
        check_graph(g, seqs) && return g
    end
    return nothing
end

function determine_graph(_cifv, _cife, cell, _seqs)
    cifv = expand_symmetry(_cifv)
    mat = Float64.(cell.mat)
    seqs = _seqs[cifv.ids]
    _pos = Float32.(cifv.pos)
    pos = _pos
    round_i = 1
    while round_i ≤ 8
        pos = [SVector{3,Float32}(x) for x in eachcol(round.(_pos; digits=round_i))]
        round_i += 1
        allunique(pos) && break
    end
    pos = [SVector{3,Float32}(x) for x in eachcol(precise_round.(_pos, round_i))]
    @assert allunique(pos)
    dpos = Dict{SVector{3,Float32},Int}([p => j for (j,p) in enumerate(pos)])
    ε = cbrt(det(mat) / length(cifv.ids))

    # g1 = try_from_closest(_cifv, cifv, mat, round_i, pos, dpos, ε, seqs)
    g1 = try_from_closest(mat, pos, seqs)
    g1 === nothing || return g1

    g2 = try_from_edges(_cife, mat, round_i, pos, dpos, ε, seqs)
    g2 === nothing || return g2

    return nothing
end


"""
Organizes the RCSR web data into a list of names, CIF (for vertices and edges), cells and
coordination sequences. Also records invalid symmetries.
"""
function extract_rcsr_data(rcsr=rcsr3D)
    nrcsr = length(rcsr)
    symmetry_issues = Tuple{String,String,Int}[]
    names = String[]
    cifvs = CIF[]
    cifes = CIF[]
    cells = Cell[]
    seqss = Vector{Vector{Int}}[]
    i = 2
    while i < nrcsr
        last_i = i
        strip(rcsr[i]) == "-1" && break
        i += 1
        name = strip(rcsr[i])
        i += 2
        try
            weaving = false
            for j in 1:5
                name ∈ ("cdz-e", "ssf-e", "bor-y", "pok", "qok") || @assert rcsr[i][1] == ' '
                splits = split(rcsr[i])
                num = parse(Int, first(splits))
                if j == 4
                    @assert last(splits) == "keywords"
                    if num != 0
                        for _ in 1:num
                            i += 1
                            if strip(rcsr[i]) == "weaving"
                                weaving = true
                            end
                        end
                    end
                    i += 1
                else
                    i += num + 1
                end
            end

            _symb, _spgroup = split(rcsr[i])
            spgroup = parse(Int, _spgroup)
            symb = filter(x -> x != '(' && x != ')', _symb)
            hall = get(CrystalNets.SPACE_GROUP_HM, symb, 0)
            if hall == 0
                hall = CrystalNets.SPACE_GROUP_IT[spgroup]
                push!(symmetry_issues, (name, _symb, spgroup))
                # @info "Invalid symmetry $symb for $name (defaulting to $hall from spgroup $spgroup)"
            elseif hall != CrystalNets.SPACE_GROUP_IT[spgroup]
                @warn "symb $symb leads to hall number $hall but spgroup $spgroup leads to $(CrystalNets.SPACE_GROUP_IT[spgroup])"
            end

            i += 1
            a, b, c, α, β, γ = parse.(Float64, split(rcsr[i]))
            cell = Cell(hall, (a, b, c), (α, β, γ))

            i += 1
            numv = parse(Int, rcsr[i])
            posv = Matrix{Float64}(undef, 3, numv)
            coordination = Vector{Int}(undef, numv)

            i += 1
            for v in 1:numv
                splits = split(rcsr[i])
                name ∈ ("odf-d", "gwe-a", "qtz-t", "moo", "cot-a", "fnh-b", "zaz", "lwa-d") || @assert splits[1] == "V"*string(v)
                coordination[v] = parse(Int, splits[2])
                i += 1
                posv[:, v] .= parse.(Float64, split(rcsr[i]))
                i += 5
            end

            nume = parse(Int, rcsr[i])
            pose = Matrix{Float64}(undef, 3, nume)
            i += 1
            for e in 1:nume
                splits = split(rcsr[i], x -> isspace(x) || !isprint(x); keepempty=false)
                @assert splits[2] == "2"
                name ∈ ("ntb", "bcu-dia-c", "oku") || splits[1] == "E"*string(e) || @warn "$(splits[1]) != E$e"
                i += 1
                pose[:, e] .= parse.(Float64, split(rcsr[i]))
                i += 4
            end

            i += 5
            seqs = Vector{Vector{Int}}(undef, numv)
            for v in 1:numv
                seq = parse.(Int, split(rcsr[i]))
                @assert length(seq) == 11
                name ∈ ("xxv", "rpa", "qyc") || @assert seq[1] == coordination[v]
                pop!(seq)
                seqs[v] = seq
                i += 1
            end
            if name == "xxv"
                @assert length(seqs) == 2
                seqs[1], seqs[2] = seqs[2], seqs[1]
                @assert seqs[1][1] == coordination[1]
                @assert seqs[2][1] == coordination[2]
            elseif name == "qyc"
                @assert seqs[3][1] == 3
                seqs[3][1] = 4
            end

            if !weaving # store the net
                push!(names, name)
                push!(cifvs, CIF(Dict{String, Union{String, Vector{String}}}(), cell,
                    collect(1:numv), [Symbol("") for _ in 1:numv], posv, Vector{Tuple{Int,Float32}}[]))
                push!(cifes, CIF(Dict{String, Union{String, Vector{String}}}(), cell,
                    collect(1:nume), [Symbol("") for _ in 1:nume], pose, Vector{Tuple{Int,Float32}}[]))
                push!(cells, cell)
                push!(seqss, seqs)
            end

            while i < length(rcsr) && rcsr[i] != "start"
                i += 1
            end
            i += 1
        catch e
            @show name, last_i
            rethrow()
        end
    end
    return names, cifvs, cifes, cells, seqss, symmetry_issues
end


"""
Automatically collect the periodic graph corresponding to the nets in the RCSR.
"""
function extract_graphs(rcsr=rcsr3D)
    names, cifvs, cifes, cells, seqss, symmetry_issues = extract_rcsr_data(rcsr)
    n = length(names)
    ret = Vector{Pair{String,PeriodicGraph3D}}(undef, n)
    errored = [Int[] for _ in 1:nthreads()]
    failed = [Int[] for _ in 1:nthreads()]
    @threads for i in 1:n
        name = names[i]
        cifv = cifvs[i]
        cife = cifes[i]
        cell = cells[i]
        seqs = seqss[i]
        graph = try
            determine_graph(cifv, cife, cell, seqs)
        catch e
            push!(errored[threadid()], i)
            continue
        end
        if graph isa PeriodicGraph3D
            ret[i] = name => graph
        else
            push!(failed[threadid()], i)
        end
    end
    _failed = reduce(vcat, failed)
    _errored = reduce(vcat, errored)
    toremove = vcat(_failed, _errored)
    sort!(toremove)
    deleteat!(ret, toremove)

    return Dict(ret), names[_failed], names[_errored], symmetry_issues
end



# Utility for comparing with EPINET

function altnamesrcsr(rcsr=rcsr3D)
    ret = Pair{String,Vector{String}}[]
    @assert rcsr[1] == "start"
    i = 1
    flag = true
    while flag
        i += 1
        num = rcsr[i]
        num == "-1" && break
        i += 1
        net = rcsr[i]
        @assert !haskey(Dict(ret), net)
        while true
            i += 1
            if i > length(rcsr)
                flag = false
                break
            end
            l = rcsr[i]
            if endswith(l, "number of names")
                num = parse(Int, split(l)[1])
                altnames = String[]
                for _ in 1:num
                    i += 1
                    l = rcsr[i]
                    push!(altnames, strip(l))
                end
                push!(ret, (strip(net) => altnames))
                break
            end
        end
        while true
            i += 1
            if i > length(rcsr)
                flag = false
                break
            end
            l = rcsr[i]
            l == "start" && break
        end
    end
    d = Dict(ret)
    @assert length(d) == length(ret)
    return d
end

function epinet_comparison(altnames=altnamesrcsr())
    ret = Pair{String,Int}[]
    for (k, v) in altnames
        for name in v
            if startswith(name, "sqc") && all(isnumeric, name[4:end])
                push!(ret, (k => parse(Int, name[4:end])))
            end
        end
    end
    d = Dict(ret)
    @assert length(d) == length(ret)
    return d
end



