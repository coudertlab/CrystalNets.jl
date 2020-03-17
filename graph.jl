include("types.jl")
using .CIFTypes
using LightGraphs

using GraphPlot

const ORGANIC_ELEMENTS = (:C, :H, :O, :N)

function find_identical_vertex(bonds::BitMatrix)
end

function coordination_sequence(bonds::BitMatrix)

end

function identify_subgroups(c::Crystal, x::Symbol)
    @warn "Not finished"
    subgroups = Vector{Int}[]
    xgroup = ===(x).(c.atoms)
    n = length(xgroup)
    for _i in 1:n, _j in (i+1):n
        i = xgroup[_i]; j = xgroup[_j]
        # TODO
    end
end

function find_clusters(c::Crystal, tol=5)
    g = c.graph
    components = LightGraphs.connected_components(g)
    if length(components) == 1
        return find_clusters(c, c.graph; tol=tol)
    end
    ret = Vector{Int}[]
    for comp in components
        sg, vmap = LightGraphs.induced_subgraph(g, comp)
        subcrystal = Crystal(length(vmap), c.atoms[vmap], c.bonds[vmap, vmap],
                             c.pos[:,vmap], sg)
        subdecomposition = find_clusters(subcrystal; tol=tol)
        for cluster in subdecomposition
            push!(ret, [vmap[i] for i in cluster])
        end
    end
    return ret
end

"""
    find_clusters_connected(g::AbstractGraph)

Separate the set of vertices into clusters that are repeated
"""
function find_clusters(c::Crystal; tol)
    @warn "Not finished"
    g = c.graph
    if nv(g) <= tol
        return 1:nv(g)
    end
    elements = unique(c.atoms)
    metallic = [x for x ∈ elements if x ∉ ORGANIC_ELEMENTS]
    for metal in metallic
        subgroups = identify_subgroups(c, metal)
        # TODO
    end
end

function find_equilibrium
    #TODO
end
