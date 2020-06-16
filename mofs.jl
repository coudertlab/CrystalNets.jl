using Statistics: mean

function find_sbus(net::CrystalNet)
    n = nv(net.graph)
    classes = Vector{Int}(undef, n)
    for i in 1:n
        typ = net.types[i]
        if typ === :C || typ === :H
            classes[i] = 1
        elseif typ === :O || typ === :Cu
            classes[i] = 2
        else
            error("Unknown atom type")
        end
    end
    sbus = Vector{PeriodicVertex3D}[]
    sbu_classes = Int[]
    attributions = zeros(Int, n)
    offsets = Vector{SVector{3,Int}}(undef, n)
    for i in 1:n
        iszero(attributions[i]) || continue
        class = classes[i]
        push!(sbus, [PeriodicVertex3D(i)])
        push!(sbu_classes, class)
        attr = length(sbus)
        attributions[i] = attr
        Q = PeriodicVertex3D[PeriodicVertex3D(i)]
        for u in Q
            @assert attributions[u.v] == attr
            for neigh in neighbors(net.graph, u.v)
                x = neigh.v
                if classes[x] == class
                    attrx = attributions[x]
                    if iszero(attrx)
                        ofs = neigh.ofs .+ u.ofs
                        offsets[x] = ofs
                        attributions[x] = attr
                        y = PeriodicVertex3D(x, ofs)
                        push!(sbus[end], y)
                        push!(Q, y)
                    else
                        @assert attrx == attr
                    end
                end
            end
        end
    end
    return sbus, sbu_classes, attributions, offsets
end

function coalesce_sbus(net::CrystalNet{T}, (sbus, sbu_classes, attributions, offsets)) where T
    n = length(sbus)
    pos = Vector{SVector{3,T}}(undef, n)
    types = Vector{Symbol}(undef, n)
    for (i, sbu) in enumerate(sbus)
        pos[i] = mean(net.pos[x.v] .+ x.ofs for x in sbu)
        types[i] = Symbol(sbu_classes[i]) #Symbol(join(sort!([net.types[x.v] for x in sbu])))
    end
    edgs = PeriodicEdge3D[]
    for e in edges(net.graph)
        s, d, of = src(e), dst(e), ofs(e)
        atts = attributions[s]
        attd = attributions[d]
        attributions[s] == attributions[d] && continue
        push!(edgs, PeriodicEdge3D(atts, attd, of .+ offsets[s] .- offsets[d]))
    end
    return CrystalNet(net.cell, types, pos, PeriodicGraph3D(edgs))
end
