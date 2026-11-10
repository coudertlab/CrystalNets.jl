# Paired collisions with a single global exchange invisible to the placement.
# Larger kernels, local independent exchanges, typed nets and mapping tracking
# need separate treatment; this helper does not claim general canonicalization.
function paired_translation_variants(g::PeriodicGraph{D}) where D
    groups = UnderlyingNets(g; ignore_types=true, authorize_pruning=false, trim=false, throw_error=true)
    nets = getproperty(groups, Symbol("D", D))
    length(nets) == 1 || return nothing
    net = only(first(only(nets)))
    _, (equiv, collisions) = collision_nodes(net)
    isempty(collisions) && return nothing
    first(first(collisions).rnge) == 1 || return nothing
    all(c -> length(c.rnge) == 2, collisions) || return nothing
    native = equiv.pge.g
    n, m = nv(native), length(collisions)
    flip, block = collect(1:n), zeros(Int, n)
    for (i, c) in enumerate(collisions)
        a, b = c.rnge
        flip[a], flip[b] = b, a
        block[a] = block[b] = i
    end
    all(has_edge(native, flip[e.src], PeriodicVertex{D}(flip[e.dst.v], e.dst.ofs))
        for e in edges(native)) || return nothing

    # Pair-flip constraints must connect the entire infinite graph. A connected
    # quotient alone would incorrectly admit independent flips in different cells.
    constraints = PeriodicEdge{D}[]
    for u in 1:n, v in neighbors(native, u)
        a, b = block[u], block[v.v]
        a == b && iszero(v.ofs) && continue
        if !has_edge(native, flip[u], PeriodicVertex{D}(v.v, v.ofs))
            push!(constraints, PeriodicEdge{D}(a, b, v.ofs))
        end
    end
    components = dimensionality(PeriodicGraph{D}(m, constraints))
    length(components) == 1 && haskey(components, D) || return nothing
    length(components[D]) == 1 || return nothing
    vertices, multiplicity = only(components[D])
    length(vertices) == m && multiplicity == 1 || return nothing

    variants = PeriodicGraph{D}[]
    for mask in 0:((1 << D)-1)
        edgs = PeriodicEdge{D}[]
        # Transform both directions before choosing undirected representatives.
        for u in 1:n, v in neighbors(native, u)
            parity = false
            for axis in 1:D
                iszero(mask & (1 << (axis-1))) || (parity ⊻= isodd(v.ofs[axis]))
            end
            target = parity ? flip[v.v] : v.v
            push!(edgs, PeriodicEdge{D}(u, target, v.ofs))
        end
        other = PeriodicGraph{D}(n, edgs)
        @toggleassert nv(other) == n && ne(other) == ne(native)
        @toggleassert degree(other) == degree(native)
        push!(variants, other)
    end
    variants
end

function normalize_paired_genome(original::PeriodicGraph{D}) where D
    pending = [string(original)]
    orbit = Set(pending)
    while !isempty(pending)
        variants = paired_translation_variants(PeriodicGraph{D}(pop!(pending)))
        isnothing(variants) && continue
        for graph in variants
            groups = UnderlyingNets(graph; ignore_types=true, authorize_pruning=false, trim=false, throw_error=true)
            nets = getproperty(groups, Symbol("D", D))
            length(nets) == 1 || error("Paired translation changed component multiplicity")
            net = only(first(only(nets)))
            result = topological_genome(net; _normalize_paired=false)
            isempty(result.error) || error(result.error)
            nv(result.genome) == nv(original) || error("Unsupported paired orbit: minimized cell size changes")
            candidate = string(result.genome)
            if candidate ∉ orbit
                push!(orbit, candidate)
                push!(pending, candidate)
                length(orbit) <= 64 || error("Paired translation orbit limit exceeded")
            end
        end
    end
    key = first(sort!(collect(orbit); by=s -> Tuple(parse.(Int, split(s)))))
    PeriodicGraph{D}(key)
end
