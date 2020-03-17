using LightGraphs


module PeriodicGraphs

using LightGraphs
export PeriodicEdge3D, PeriodicGraph3D, LightGraphs, cellgraph, periodiccellgraph
import Base:(==)

struct PeriodicEdge3D <: LightGraphs.SimpleGraphs.AbstractSimpleEdge{Int}
    src::Int
    dst::Int
    offset::Tuple{Int, Int, Int}
    function PeriodicEdge3D(src::Int, dst::Int, offset::Tuple{Int, Int, Int}, check)
        if check
            src == dst && offset == (0, 0, 0) && throw("Loops are forbidden edges : PeriodicEdge3D($((src, dst, offset))) is invalid")
            src > dst && begin src, dst, offset = dst, src, .-offset end#throw("Order bonds such that src <= dst : PeriodicEdge3D($((src, dst, offset))) is invalid")
            src == dst && offset > .-offset && (offset = .-offset)
        end
        return new(src, dst, offset)
    end
end
function PeriodicEdge3D(src::Int, dst::Int, offset::Tuple{Int, Int, Int})
    PeriodicEdge3D(src, dst, offset, true)
end
function Base.convert(::Type{<:PeriodicEdge3D}, (src, dst, offset)::Tuple{T,T,Tuple{T,T,T}}) where {T<:Integer}
    PeriodicEdge3D(src, dst, offset, true)
end
function LightGraphs.reverse(e::PeriodicEdge3D)
    if e.src == e.dst
        return e
    end
    return PeriodicEdge3D(e.dst, e.src, .-e.offset, false)
end
LightGraphs.dst(e::PeriodicEdge3D) = e.dst
LightGraphs.src(e::PeriodicEdge3D) = e.src
function Base.isless(x::PeriodicEdge3D, y::PeriodicEdge3D)
    return isless(x.src, y.src) || (isequal(x.src, y.src) && isless(x.dst, y.dst))
end
function ==(x::PeriodicEdge3D, y::PeriodicEdge3D)
    return isequal(x.src, y.src) && isequal(x.dst, y.dst) && isequal(x.offset, y.offset)
end
function Base.show(io::IO, x::PeriodicEdge3D)
    print(io, "PeriodicEdge3D($(x.src), $(x.dst), $(x.offset))")
end


struct PeriodicGraph3D <: AbstractGraph{Int}
    ne::Ref{Int}
    edges::Vector{PeriodicEdge3D}
    indices::Vector{Int}
end

function PeriodicGraph3D(n::Integer)
    @assert n >= 0
    return PeriodicGraph3D(Ref(0), PeriodicEdge3D[], [1 for _ in 1:n+1])
end
function PeriodicGraph3D(t::Vector{PeriodicEdge3D})
    sort!(t); unique!(t)
    nv = maximum(e.dst for e in t)
    ne = length(t)
    g = PeriodicGraph3D(nv)
    for e in t
        add_edge!(g, e)
    end
    return g
end

@inline function search_dst_between(l, i, j, d, offset)
    i == j && return i
    while j-i > 1
        m = div(i+j, 2)
        if l[m].dst < d || (l[m].dst == d && l[m].offset <= offset)
            i = m
        else
            j = m
        end
    end
    (l[i].dst < d || (l[i].dst == d && l[i].offset < offset)) && return i+1
    return i
end

@inline function search_dst_between(l, i, j, d)
    i > length(l) && return i
    while j-i > 1
        m = div(i+j, 2)
        if l[m].dst <= d
            i = m
        else
            j = m
        end
    end
    l[i].dst < d && return i+1
    return i
end

function _add_edge!(g::PeriodicGraph3D, e::PeriodicEdge3D)
    s, d, ofs = e.src, e.dst, e.offset
    rev = reverse(e)
    @inbounds i = search_dst_between(g.edges, g.indices[s], g.indices[s+1], d, ofs)
    i <= length(g.edges) && g.edges[i] == e && return false
    insert!(g.edges, i, e)
    @inbounds for k in (s+1):length(g.indices)
        g.indices[k] += 1
    end
    return true
end
function LightGraphs.add_edge!(g::PeriodicGraph3D, e::PeriodicEdge3D)
    success = _add_edge!(g, e) && (e.src == e.dst || _add_edge!(g, reverse(e)))
    g.ne[] += success
    return success
end
function LightGraphs.add_edge!(g::PeriodicGraph3D, src::Int, dst::Int, offset::Tuple{Int, Int, Int})
    add_edge!(g, PeriodicEdge3D(src, dst, offset))
end

LightGraphs.ne(g::PeriodicGraph3D) = g.ne[]
LightGraphs.nv(g::PeriodicGraph3D) = length(g.indices)-1
LightGraphs.edges(g::PeriodicGraph3D) = g.edges
LightGraphs.vertices(g::PeriodicGraph3D) = Base.OneTo(nv(g))
Base.eltype(g::PeriodicGraph3D) = Int
LightGraphs.edgetype(::PeriodicGraph3D) = PeriodicEdge3D
function LightGraphs.has_edge(g::PeriodicGraph3D, s::Int, d::Int)
    (s < 1 || s >= length(g.indices)) && return false
    @inbounds i = search_dst_between(g.edges, g.indices[s], g.indices[s+1], d)
    return i <= length(g.edges) && @inbounds g.edges[i].dst == d
end
function LightGraphs.has_edge(g::PeriodicGraph3D, e::PeriodicEdge3D)
    s, d, ofs = e.src, e.dst, e.offset
    (s < 1 || s >= length(g.indices)) && return false
    @inbounds i = search_dst_between(g.edges, g.indices[s], g.indices[s+1], d, ofs)
    return i <= length(g.edges) && @inbounds g.edges[i].dst == d && g.edges[i].offset == ofs
end
function LightGraphs.outneighbors(g::PeriodicGraph3D, (v, offset)::Tuple{Int, Tuple{Int, Int, Int}})
    ret = Tuple{Int, Tuple{Int, Int, Int}}[]
    self_neighbours = Tuple{Int, Tuple{Int, Int, Int}}[]
    @inbounds for i in g.indices[v]:(g.indices[v+1]-1)
        e = g.edges[i]
        push!(ret, (e.dst, e.offset .+ offset))
        if e.src == e.dst
            push!(self_neighbours, (e.dst, offset .- e.offset))
        else !isempty(self_neighbours)
            for _ in 1:length(self_neighbours)
                push!(ret, pop!(self_neighbours))
            end
        end
    end
    return ret
end
LightGraphs.outneighbors(g::PeriodicGraph3D, v::Integer) = outneighbors(g, (v, (0, 0, 0)))
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
            for (vertex, _) in all_neighbors(g, src)
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
    Q = Vector{Tuple{Tuple{Int, Tuple{Int, Int, Int}}, U}}()
    d < zero(U) && return Q
    push!(Q, ((v, (0, 0, 0)), zero(U),) )
    seen = Set{Tuple{Int, Tuple{Int, Int, Int}}}([(v, (0, 0, 0))])
    @inbounds for (src, currdist) in Q
        currdist >= d && continue
        for dst::Tuple{Int, Tuple{Int, Int, Int}} in neighborfn(g, src)
            if dst ∉ seen
                push!(seen, dst)
                if currdist + distmx[first(src), first(dst)] <= d
                    push!(Q, (dst , currdist + distmx[first(src), first(dst)],))
                end
            end
        end
    end
    return Q
end

function cellgraph(g::PeriodicGraph3D)
    ret = SimpleGraph([Edge{Int}(e.src, e.dst) for e in g.edges if e.offset == (0, 0, 0)])
    add_vertices!(ret, nv(g) - nv(ret))
    return ret
end

function periodiccellgraph(g::PeriodicGraph3D)
    ret = SimpleGraph([Edge{Int}(e.src, e.dst) for e in g.edges])
    add_vertices!(ret, nv(g) - nv(ret))
    return ret
end

end
