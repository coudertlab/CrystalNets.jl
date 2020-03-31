module PeriodicGraphs

using Reexport
using LinearAlgebra
@reexport using LightGraphs
export PeriodicVertex, PeriodicEdge3D, PeriodicGraph3D,
       cellgraph, periodiccellgraph, equilibrium
import Base:(==)

struct PeriodicVertex
    v::Int
    ofs::Tuple{Int, Int, Int}
end
PeriodicVertex(n::Integer) = PeriodicVertex(n, (0, 0, 0))
function Base.show(io::IO, x::PeriodicVertex)
    print(io, "PeriodicVertex($(x.v), $(x.ofs))")
end
function Base.convert(::Type{<:PeriodicVertex}, (dst, offset)::Tuple{Int, Tuple{Int, Int, Int}})
    PeriodicVertex(dst, offset)
end
function Base.isless(x::PeriodicVertex, y::PeriodicVertex)
    isless(x.v, y.v) || (isequal(x.v, y.v) && isless(x.ofs, y.ofs))
end
function ==(x::PeriodicVertex, y::PeriodicVertex)
    isequal(x.v, y.v) && isequal(x.ofs, y.ofs)
end
ZtoN(x::Int) = -(x<0) + 2*abs(x)

"""
    hash_position((x1,x2,x3)::Tuple{Int, Int, Int})

Bijection from Z^3 -> N.
For n = max(abs.(x1, x2, x3)), m = max(abs.(y1, y2, y3)), the hash is such that
n < m implies hash_position((x1, x2, x3)) < hash_position((y1, y2, y3))

"""
function hash_position((x1,x2,x3)::Tuple{Int, Int, Int})
    _b = x1 >= x2
    b1 = _b & (x1 >= x3)
    b2 = (!_b) & (x2 >= x3)
    b3 = (!b1) & (!b2)
    return b1*(x2*(x1 + 1) + x3 + x1^3) +
           b2*((x2 + 1)*(x1 + x2 + 1) + x3 + x2^3) +
           b3*((x3 + 1)*(2x3 + 1) + x3*(x1 + x3^2) + x2)
end
function hash_position(x::PeriodicVertex, n::Int)
    v, t::Tuple{Int, Int, Int} = x.v, ZtoN.(x.ofs)
    return v + n*hash_position(t)
end

struct PeriodicEdge3D <: LightGraphs.SimpleGraphs.AbstractSimpleEdge{Int}
    src::Int
    dst::PeriodicVertex
    function PeriodicEdge3D(src::Int, dst::Int, offset::Tuple{Int, Int, Int}, check)
        if check
            src == dst && offset == (0, 0, 0) && throw("Loops are forbidden edges : PeriodicEdge3D($((src, dst, offset))) is invalid")
            src > dst && begin src, dst, offset = dst, src, .-offset end#throw("Order bonds such that src <= dst : PeriodicEdge3D($((src, dst, offset))) is invalid")
            src == dst && offset > .-offset && (offset = .-offset)
        end
        return new(src, PeriodicVertex(dst, offset))
    end
end
function PeriodicEdge3D(src::Int, dst::Int, offset::Tuple{Int, Int, Int})
    PeriodicEdge3D(src, dst, offset, true)
end
function PeriodicEdge3D(src::Int, dst::PeriodicVertex)
    PeriodicEdge3D(src, dst.v, dst.ofs, false)
end
function Base.convert(::Type{<:PeriodicEdge3D}, (src, dst, offset)::Tuple{T,T,Tuple{T,T,T}}) where {T<:Integer}
    PeriodicEdge3D(src, dst, offset, true)
end
function LightGraphs.reverse(e::PeriodicEdge3D)
    if e.src == e.dst
        return e
    end
    return PeriodicEdge3D(e.dst.v, e.src, .-e.dst.ofs, false)
end
LightGraphs.src(e::PeriodicEdge3D) = e.src
LightGraphs.dst(e::PeriodicEdge3D) = e.dst.v
function Base.isless(x::PeriodicEdge3D, y::PeriodicEdge3D)
    return isless(x.src, y.src) || (isequal(x.src, y.src) && isless(x.dst, y.dst))
end
function ==(x::PeriodicEdge3D, y::PeriodicEdge3D)
    return isequal(x.src, y.src) && isequal(x.dst, y.dst)
end
function Base.show(io::IO, x::PeriodicEdge3D)
    print(io, "PeriodicEdge3D($(x.src), $(x.dst))")
end

struct PeriodicGraph3D <: AbstractGraph{Int}
    edges::Vector{PeriodicVertex}
    indices::Vector{Int}
end

function PeriodicGraph3D(n::Integer)
    @assert n >= 0
    return PeriodicGraph3D(PeriodicVertex[], [1 for _ in 1:n+1])
end
function PeriodicGraph3D(nv::Int, t::Vector{PeriodicEdge3D})
    sort!(t); unique!(t)
    ne = length(t)
    g = PeriodicGraph3D(nv)
    for e in t
        add_edge!(g, e)
    end
    return g
end
function PeriodicGraph3D(t::Vector{PeriodicEdge3D})
    return PeriodicGraph3D(maximum(max(e.src, e.dst.v) for e in t), t)
end

function ==(g1::PeriodicGraph3D, g2::PeriodicGraph3D)
    return g1.indices == g2.indices && g1.edges == g2.edges
end

"""
    search_dst_between(l, i, j, d, offset)

For a sorted list of edges l and given starting and index indices i and j,
return the index t at which PeriodicEdge(d, offset) is present in the list if
so, or at which to add the edge if absent.
"""
@inline function search_dst_between(l, i, j, d, offset)
    i == j && return i
    while j-i > 1
        m = div(i+j, 2)
        if l[m].v < d || (l[m].v == d && l[m].ofs <= offset)
            i = m
        else
            j = m
        end
    end
    (l[i].v < d || (l[i].v == d && l[i].ofs < offset)) && return i+1
    return i
end

@inline function search_dst_between(l, i, j, d)
    i > length(l) && return i
    while j-i > 1
        m = div(i+j, 2)
        if l[m].v <= d
            i = m
        else
            j = m
        end
    end
    l[i].v < d && return i+1
    return i
end

function _add_edge!(g::PeriodicGraph3D, e::PeriodicEdge3D)
    s, dst = e.src, e.dst
    rev = reverse(e)
    @inbounds i = search_dst_between(g.edges, g.indices[s], g.indices[s+1], dst.v, dst.ofs)
    i <= length(g.edges) && g.edges[i] == dst && return false
    insert!(g.edges, i, dst)
    @inbounds for k in (s+1):length(g.indices)
        g.indices[k] += 1
    end
    return true
end
function LightGraphs.add_edge!(g::PeriodicGraph3D, e::PeriodicEdge3D)
    success = _add_edge!(g, e) && _add_edge!(g, reverse(e))
    return success
end
function LightGraphs.add_edge!(g::PeriodicGraph3D, src::Int, dst::Int, offset::Tuple{Int, Int, Int})
    add_edge!(g, PeriodicEdge3D(src, dst, offset))
end

LightGraphs.ne(g::PeriodicGraph3D) = length(g.edges) >> 1
LightGraphs.nv(g::PeriodicGraph3D) = length(g.indices) - 1
LightGraphs.vertices(g::PeriodicGraph3D) = Base.OneTo(nv(g))
function LightGraphs.edges(g::PeriodicGraph3D)
    ret = PeriodicEdge3D[]
    for i in vertices(g)
        for j in g.indices[i]:g.indices[i+1]-1
            push!(ret, PeriodicEdge3D(i, g.edges[j]))
        end
    end
    return ret
 end
Base.eltype(g::PeriodicGraph3D) = PeriodicVertex
LightGraphs.edgetype(::PeriodicGraph3D) = PeriodicEdge3D
function LightGraphs.has_edge(g::PeriodicGraph3D, s::Int, d::Int)
    (s < 1 || s >= length(g.indices)) && return false
    @inbounds i = search_dst_between(g.edges, g.indices[s], g.indices[s+1], d)
    return i <= length(g.edges) && @inbounds g.edges[i].v == d
end
function LightGraphs.has_edge(g::PeriodicGraph3D, e::PeriodicEdge3D)
    s, d, ofs = e.src, e.dst.v, e.dst.ofs
    (s < 1 || s >= length(g.indices)) && return false
    @inbounds i = search_dst_between(g.edges, g.indices[s], g.indices[s+1], d, ofs)
    return i <= length(g.edges) && @inbounds g.edges[i] == e.dst
end
function LightGraphs.outneighbors(g::PeriodicGraph3D, v::Int)
    return @inbounds @view g.edges[g.indices[v]:g.indices[v+1]-1]
end
@inline function num_neighbors(g::PeriodicGraph3D, v::Int)
    return g.indices[v+1] - g.indices[v]
end
LightGraphs.inneighbors(g::PeriodicGraph3D, v::Integer) = outneighbors(g, v)
Base.zero(::Type{PeriodicGraph3D}) = PeriodicGraph3D(0)
LightGraphs.is_directed(::Type{<:PeriodicGraph3D}) = false
@inline has_contiguous_vertices(::Type{<:PeriodicGraph3D}) = true
LightGraphs.has_vertex(g::PeriodicGraph3D, v::Integer) = 1 <= v <= nv(g)
function LightGraphs.SimpleGraphs.add_vertices!(g::PeriodicGraph3D, n::Integer)
    append!(g.indices, fill(length(g.edges) + 1, n))
    return n
end
LightGraphs.SimpleGraphs.add_vertex!(g::PeriodicGraph3D) = (add_vertices!(g, 1); true)

function make_local_clique!(g::PeriodicGraph3D, v::Integer)
    for i in g.indices[v]:g.indices[v+1]-1
        u1 = g.edges[i]
        for j in (i+1):g.indices[v+1]-1
            u2 = g.edges[j]
            add_edge!(g, u1.v, u2.v, u2.ofs .- u1.ofs)
        end
    end
    nothing
end

# function LightGraphs.SimpleGraphs.rem_vertex!(g::PeriodicGraph3D, v::Integer)
#     n = g.indices[v+1] - g.indices[v]
#     for i in g.indices[v]:g.indices[v+1]-1
#         u = g.edges[i]
#         j = search_dst_between(g.edges, g.indices[u.v]-i+1, g.indices[u.v+1]-i+1, v, .-u.ofs)
#         @assert g.edges[j] == PeriodicVertex(v, .-u.ofs)
#         deleteat!(g.edges, j) # et corriger les voisins d'après
#         throw("À compléter")
#     end
#
#     deleteat!(g.edges, g.indices[v]:g.indices[v+1]-1)
#     deleteat!(g.indices, v)
#     for i in v:length(g.indices)
#         g.indices[i] -= n
#     end
#     throw("À compléter")
# end

function LightGraphs.SimpleGraphs.rem_vertices!(g::PeriodicGraph3D, t::Vector{T} where T<:Integer, keep_order=true)
    @assert keep_order
    sort!(t)
    while !isempty(t) && first(t) < 1
        popfirst!(t)
    end
    while !isempty(t) && last(t) > nv(g)
        pop!(t)
    end
    isempty(t) && return collect(1:nv(g))
    n = length(t)
    edges = copy(g.edges)
    indices = copy(g.indices)
    empty!(g.edges)
    empty!(g.indices)
    i_next_vertex_to_del = 1
    next_vertex_to_del = t[i_next_vertex_to_del]

    vmap = Int[]
    rev_vmap = Int[]
    for i in 1:length(indices)-1
        if i == next_vertex_to_del
            if i_next_vertex_to_del == n
                next_vertex_to_del = length(indices)
            else
                i_next_vertex_to_del += 1
                next_vertex_to_del = t[i_next_vertex_to_del]
            end
            push!(rev_vmap, 0)
            continue
        end
        push!(vmap, i)
        push!(rev_vmap, length(vmap))
    end

    for i in vmap
        push!(g.indices, length(g.edges)+1)
        for j in indices[i]:indices[i+1]-1
            e = edges[j]
            rev_vmap[e.v] == 0 && continue
            push!(g.edges, PeriodicVertex(rev_vmap[e.v], e.ofs))
        end
    end
    push!(g.indices, length(g.edges)+1)
    return vmap
end

function LightGraphs.SimpleGraphs.rem_vertex!(g::PeriodicGraph3D, v::Integer)
    n = length(g.indices)
    return length(rem_vertices!(g, [v])) == n - 2
end

function LightGraphs.connected_components(g::PeriodicGraph3D)
    nvg = nv(g)
    label = zeros(Int, nvg)

    for u in vertices(g)
        label[u] != 0 && continue
        label[u] = u
        Q = Int[]
        push!(Q, u)
        @inbounds while !isempty(Q)
            src = popfirst!(Q)
            for dst in outneighbors(g, src)
                vertex = dst.v
                if label[vertex] == 0
                    push!(Q, vertex)
                    label[vertex] = u
                end
            end
        end
    end
    c, d = LightGraphs.components(label)
    return c
end

function LightGraphs._neighborhood(g::PeriodicGraph3D, v::Integer, d::Real, distmx::AbstractMatrix{U}, neighborfn::Function) where U <: Real
    @assert typeof(neighborfn) === typeof(outneighbors)
    Q = Vector{Tuple{PeriodicVertex, U}}()
    d < zero(U) && return Q
    start_vertex = PeriodicVertex(v, (0, 0, 0))
    push!(Q, (start_vertex, zero(U),) )
    n = nv(g)
    seen = BitArray(0 for _ in 1:(n*(2d+1)^3))
    seen[hash_position(start_vertex, n)] = true
    @inbounds for (src, currdist) in Q
        currdist >= d && continue
        @simd for dst in outneighbors(g, src.v)
            dst = PeriodicVertex(dst.v, dst.ofs .+ src.ofs)
            position = hash_position(dst, n)
            if !seen[position]
                seen[position] = true
                if currdist + distmx[src.v, dst.v] <= d
                    push!(Q, (dst , currdist + distmx[src.v, dst.v],))
                end
            end
        end
    end
    return Q
end

function cellgraph(g::PeriodicGraph3D)
    edges = Edge{Int}[]
    for i in vertices(g)
        for j in outneighbors(g, i)
            j.ofs == (0, 0, 0) && push!(edges, Edge{Int}(i, j.v))
        end
    end
    ret = SimpleGraph(edges)
    add_vertices!(ret, nv(g) - nv(ret))
    return ret
end

function periodiccellgraph(g::PeriodicGraph3D)
    ret = SimpleGraph([Edge{Int}(i, j.v) for i in vertices(g) for j in outneighbors(g, i)])
    add_vertices!(ret, nv(g) - nv(ret))
    return ret
end

function equilibrium(g::PeriodicGraph3D)
    n = nv(g)
    Y = Array{Rational{Int128}}(undef, 3, n)
    A = Array{Rational{Int128}}(undef, n, n)
    neigh = Array{Int}(undef, n)
    offset = Array{Int}(undef, 3)
    for i in 1:n
        m = num_neighbors(g, i)
        neigh .= 0
        offset .= 0
        for k in neighbors(g, i)
            neigh[k.v] += 1
            offset .+= k.ofs
        end
        A[i,:] = neigh ./ m
        Y[:,i] = offset ./ m
    end
    @show factorize(A)
    @show Y
    return factorize(A)\Y'
end

end
