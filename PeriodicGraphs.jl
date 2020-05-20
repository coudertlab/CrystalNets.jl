module PeriodicGraphs

using LinearAlgebra
using StaticArrays
using LightGraphs
export PeriodicVertex, PeriodicEdge, PeriodicGraph, LoopException,
       PeriodicVertex3D, PeriodicEdge3D, PeriodicGraph3D, vertex_sequence,
       ofs, cellgraph, periodiccellgraph, equilibrium, offset_representatives!,
       swap_axes!, find_edges
import Base: (==), isless, convert, show, showerror, eltype, iterate, zero,
             length, in, ndims, print, cmp
import Base.Order: Forward, Lt


struct PeriodicVertex{N}
    v::Int
    ofs::SVector{N,Int}
end
PeriodicVertex{N}(n::Integer) where {N} = PeriodicVertex{N}(n, zero(SVector{N,Int}))

function show(io::IO, x::PeriodicVertex{N}) where N
    if get(io, :typeinfo, Any) != PeriodicVertex{N}
        print(io, PeriodicVertex{N})
    end
    print(io, '(', x.v, ", ", x.ofs, ')')
end
function convert(::Type{<:PeriodicVertex{N}}, (dst, offset)::Tuple{Integer,Any}) where N
    PeriodicVertex{N}(dst, offset)
end
function cmp(x::PeriodicVertex{N}, y::PeriodicVertex{N}) where N
    c = cmp(x.v, y.v)
    iszero(c) || return c
    return cmp(x.ofs, y.ofs)
end
isless(x::PeriodicVertex{N}, y::PeriodicVertex{N}) where {N} = cmp(x, y) < 0
==(x::PeriodicVertex{N}, y::PeriodicVertex{N}) where {N} = iszero(cmp(x, y))

const PeriodicVertex2D = PeriodicVertex{2}
const PeriodicVertex3D = PeriodicVertex{3}


ZtoN(x::Signed) = -(x<0) + 2*abs(x)
function hash_position((x1,x2,x3)::SVector{3,<:Integer})
    x1 = ZtoN(x1); x2 = ZtoN(x2); x3 = ZtoN(x3)
    _b = x1 >= x2
    b1 = _b & (x1 >= x3)
    b2 = (!_b) & (x2 >= x3)
    b3 = (!b1) & (!b2)
    return b1*(x2*(x1 + 1) + x3 + x1^3) +
           b2*((x2 + 1)*(x1 + x2 + 1) + x3 + x2^3) +
           b3*((x3 + 1)*(2x3 + 1) + x3*(x1 + x3^2) + x2)
end

function hash_position((x1,x2)::SVector{2,<:Integer})
    x1 = ZtoN(x1); x2 = ZtoN(x2)
    b = x1 >= x2
    return b*(x2 + x1^2) + (!b)*(x1 + x2*(x2 + 1) + 1)
end

"""
    hash_position(x::PeriodicVertex{N}, n::Integer) where {N}

Given a PeriodicVertex{N} and the number n of vertices in a graph, compute a
unique hash for the given vertex.
This hash function is a bijection between the set of all possible vertices and
N*. Its value is an integer between 1+n*(2d-1)^N (or 1 if d = 0) and n*(2d+1)^N,
where d = max.(abs.(x.ofs)). This means that when one unit cell A is further
than another B (for the Manhattan distance), all vertices in A will have a
larger hash than all vertices in B.
"""
function hash_position(x::PeriodicVertex, n::Integer)
    return x.v + n*hash_position(x.ofs)
end


struct LoopException <: Exception
    src::Int
end
showerror(io::IO, e::LoopException) = println("LoopException: a loop from vertex $(e.src) to itself in the same unit cell is a forbidden edges. Maybe the offset is wrong?")
@noinline __throw_loopexception(src) = throw(LoopException(src))

struct unsafe_edge{N} end

struct PeriodicEdge{N} <: LightGraphs.SimpleGraphs.AbstractSimpleEdge{Int}
    src::Int
    dst::PeriodicVertex{N}
    function (::Type{unsafe_edge{N}})(src, dst) where N
        return new{N}(src, dst)
    end
end
unsafe_edge{N}(src, dst, ofs) where {N} = unsafe_edge{N}(src, PeriodicVertex{N}(dst, ofs))
function PeriodicEdge{N}(src, dst::PeriodicVertex{N}) where N
    src == dst.v && all(iszero.(dst.ofs)) && __throw_loopexception(src)
    return unsafe_edge{N}(src, dst)
end
function PeriodicEdge{N}(src, dst, offset) where N
    return PeriodicEdge{N}(src, PeriodicVertex{N}(dst, offset))
end
function PeriodicEdge{N}((src, dst, offset)::Tuple{T,T,Union{SVector{3,T},NTuple{N,T}}}) where {N,T<:Integer}
    PeriodicEdge{N}(src, dst, offset)
end
function convert(::Type{PeriodicEdge{N}}, (src, dst, offset)::Tuple{T,T,Union{SVector{3,T},NTuple{N,T}}}) where {N,T<:Integer}
    PeriodicEdge{N}(src, dst, offset)
end

function show(io::IO, x::PeriodicEdge{N}) where N
    if get(io, :typeinfo, Any) != PeriodicEdge{N}
        print(io, PeriodicEdge{N})
    end
    print(io, '(', x.src, ", ", x.dst.v, ", (", join(x.dst.ofs, ", "), ')', ')')
end

const PeriodicEdge2D = PeriodicEdge{2}
const PeriodicEdge3D = PeriodicEdge{3}

function LightGraphs.reverse(e::PeriodicEdge{N}) where N
    if e.src == e.dst
        return e
    end
    return unsafe_edge{N}(e.dst.v, e.src, .-e.dst.ofs)
end
LightGraphs.src(e::PeriodicEdge) = e.src
LightGraphs.dst(e::PeriodicEdge) = e.dst.v
ofs(e::PeriodicEdge) = e.dst.ofs
isindirectedge(e::PeriodicEdge) = e.src > e.dst.v || (e.src == e.dst.v && e.dst.ofs < zero(e.dst.ofs))
function cmp(x::PeriodicEdge{N}, y::PeriodicEdge{N}) where N
    c = cmp(x.src, y.src)
    iszero(c) || return c
    return cmp(x.dst, y.dst)
end
isless(x::PeriodicEdge{N}, y::PeriodicEdge{N}) where {N} = cmp(x, y) < 0
==(x::PeriodicEdge{N}, y::PeriodicEdge{N}) where {N} = iszero(cmp(x, y))

struct PeriodicGraph{N} <: AbstractGraph{Int}
    ne::Base.RefValue{Int}
    nlist::Vector{Vector{PeriodicVertex{N}}} # For each vertex, the sorted list of its neighbors
    directedgestart::Vector{Int}
    width::Base.RefValue{Rational{Int}}
end

const PeriodicGraph2D = PeriodicGraph{2}
const PeriodicGraph3D = PeriodicGraph{3}

function PeriodicGraph{N}(nv::Integer, t, s) where N
    return PeriodicGraph{N}(Ref(nv), t, s, Ref(-1//1))
end
function PeriodicGraph{N}(n::Integer) where N
    @assert n >= 0
    return PeriodicGraph{N}(0, [PeriodicVertex{N}[] for _ in 1:n], [1 for _ in 1:n])
end
function PeriodicGraph{N}(nv::Integer, t::Vector{PeriodicEdge{N}}) where N
    sort!(t); unique!(t)
    ne = length(t)
    g = PeriodicGraph{N}(nv)
    for e in t
        add_edge!(g, e)
    end
    return g
end
PeriodicGraph(nv::Integer, t::AbstractVector{PeriodicEdge{N}}) where {N} = PeriodicGraph{N}(nv, t)
function PeriodicGraph{N}(t::AbstractVector{PeriodicEdge{N}}) where N
    return PeriodicGraph(maximum(max(e.src, e.dst.v) for e in t), t)
end
PeriodicGraph(t::AbstractVector{PeriodicEdge{N}}) where {N} = PeriodicGraph{N}(t)
function PeriodicGraph(s::AbstractString)
    key = parse.(Int, split(s, ' '))
    N = popfirst!(key)
    edgs = PeriodicEdge{N}[]
    len = N+2
    for i in 0:(length(key) ÷ len)-1
        push!(edgs, PeriodicEdge{N}(key[i*len + 1], key[i*len + 2], key[i*len+3:(i+1)*len]))
    end
    return PeriodicGraph{N}(edgs)
end
function PeriodicGraph{N}(s::AbstractString) where N
    g = PeriodicGraph(s)
    ndims(g) == N || throw("Cannot construct a $N-dimensional graph from a $(ndims(g))-dimensional key")
    return g
end
function PeriodicGraph{N}(g::PeriodicGraph{N}) where N
    PeriodicGraph{N}(g.ne[], [copy(x) for x in g.nlist], copy(g.directedgestart))
end
PeriodicGraph(g::PeriodicGraph{N}) where {N} = PeriodicGraph{N}(g)


function show(io::IO, g::PeriodicGraph{N}) where N
    if get(io, :typeinfo, Any) != PeriodicGraph{N}
        print(io, PeriodicGraph{N})
    end
    print(io, '(', nv(g), ',', ' ', collect(edges(g)), ')')
end
function print(io::IO, g::PeriodicGraph{N}) where N
    print(io, N)
    for e in edges(g)
        print(io, ' ', e.src, ' ', e.dst.v, ' ', join(e.dst.ofs, ' '))
    end
end

function ==(g1::PeriodicGraph{N}, g2::PeriodicGraph{M}) where {N,M}
    return N == M && nv(g1) == nv(g2) && edges(g1) == edges(g2)
end

ndims(::PeriodicGraph{N}) where {N} = N

LightGraphs.ne(g::PeriodicGraph) = g.ne[]
LightGraphs.nv(g::PeriodicGraph) = length(g.nlist)
LightGraphs.vertices(g::PeriodicGraph) = Base.OneTo(nv(g))
LightGraphs.edges(g::PeriodicGraph{N}) where {N} = PeriodicEdgeIter{N}(g)
eltype(g::PeriodicGraph{N}) where {N} = PeriodicVertex{N}
LightGraphs.edgetype(::PeriodicGraph{N}) where {N} = PeriodicEdge{N}
function LightGraphs.has_edge(g::PeriodicGraph, s::Int, d::Int)
    ((s < 1) | (s > nv(g))) && return false
    #=@inbounds=# begin
        start = g.directedgestart[s]
        lo, hi = s > d ? (1, start-1) : (start, lastindex(g.nlist[s]))
        i = searchsortedfirst(g.nlist[s], d, lo, hi, Lt((x,y)->isless(x.v, y)))
        return i <= length(g.nlist[s]) && g.nlist[s][i].v == d
    end
end
function find_edges(g::PeriodicGraph, s::Int, d::Int)
    ((s < 1) | (s > nv(g))) && return false
    #=@inbounds=# begin
        start = g.directedgestart[s]
        lo, hi = s > d ? (1, start-1) : (start, lastindex(g.nlist[s]))
        rng = searchsorted(g.nlist[s], d, lo, hi, Lt((x,y)->isless(x.v, y)))
        if s == d.v
            rng = (2*first(rng) - last(rng) - 1):last(rng)
        end
        return g.nlists[rng]
    end
end
function LightGraphs.has_edge(g::PeriodicGraph, e::PeriodicEdge)
    s, d = e.src, e.dst
    ((s < 1) | (s > nv(g))) && return false
    #=@inbounds=# begin
        start = g.directedgestart[s]
        lo, hi = isindirectedge(e) ? (1, start-1) : (start, lastindex(g.nlist[s]))
        i = searchsortedfirst(g.nlist[s], d, lo, hi, Forward)
        return i <= length(g.nlist[s]) && g.nlist[s][i] == d
    end
end
LightGraphs.outneighbors(g::PeriodicGraph, v::Integer) = g.nlist[v]
LightGraphs.inneighbors(g::PeriodicGraph, v::Integer) = outneighbors(g, v)
zero(::Type{PeriodicGraph{N}}) where N = PeriodicGraph{N}(0)
LightGraphs.is_directed(::Type{<:PeriodicGraph}) = false
@assert !isdefined(LightGraphs, :has_contiguous_vertices) # FIXME when LightGraphs v1.3.1 is tagged
@inline has_contiguous_vertices(::Type{<:PeriodicGraph}) = true # TODO use LightGraphs.has_contiguous_vertices (LightGraphs version > 1.3.1)
LightGraphs.has_vertex(g::PeriodicGraph, v::Integer) = 1 <= v <= nv(g)
function LightGraphs.SimpleGraphs.add_vertices!(g::PeriodicGraph{N}, n::Integer) where N
    append!(g.nlist, [PeriodicVertex{N}[] for _ in 1:n])
    append!(g.directedgestart, [1 for _ in 1:n])
    g.width[] = -1
    return n
end
function LightGraphs.SimpleGraphs.add_vertex!(g::PeriodicGraph{N}) where N
    push!(g.nlist, PeriodicVertex{N}[])
    push!(g.directedgestart, 1)
    g.width[] = -1
    true
end

function _add_edge!(g::PeriodicGraph, e::PeriodicEdge, ::Val{check}) where check
    #=@inbounds=# begin
        s, dst = e.src, e.dst
        neigh = g.nlist[s]
        start = g.directedgestart[s]
        indirectedge = isindirectedge(e)
        lo, hi = indirectedge ? (1, start-1) : (start, lastindex(neigh))
        i = searchsortedfirst(neigh, dst, lo, hi, Forward)
        if check
            i <= length(neigh) && neigh[i] == dst && return false
        end
        g.directedgestart[s] += indirectedge
        insert!(neigh, i, dst)
        return true
    end
end
function LightGraphs.add_edge!(g::PeriodicGraph, e::PeriodicEdge)
    (src(e) < 1 || src(e) > nv(g) || dst(e) < 1 || dst(e) > nv(g)) && return false
    success = _add_edge!(g, e, Val(true)) && _add_edge!(g, reverse(e), Val(false))
    if success
        g.ne[] += 1
        g.width[] = -1
    end
    return success
end

function _rem_edge!(g::PeriodicGraph, e::PeriodicEdge, ::Val{check}) where check
    #=@inbounds=# begin
        s, dst = e.src, e.dst
        neigh = g.nlist[s]
        start = g.directedgestart[s]
        indirectedge = isindirectedge(e)
        lo, hi = indirectedge ? (1, start-1) : (start, lastindex(neigh))
        i = searchsortedfirst(neigh, dst, lo, hi, Forward)
        if check
            i <= length(neigh) && neigh[i] == dst || return false
        end
        g.directedgestart[s] -= indirectedge
        deleteat!(neigh, i)
        return true
    end
end
function LightGraphs.rem_edge!(g::PeriodicGraph, e::PeriodicEdge)
    (src(e) < 1 || src(e) > nv(g) || dst(e) < 1 || dst(e) > nv(g)) && return false
    success = _rem_edge!(g, e, Val(true)) && _rem_edge!(g, reverse(e), Val(false))
    if success
        g.ne[] -= 1
        g.width[] = -1
    end
    return success
end


# function make_local_clique!(g::PeriodicGraph3D, v::Integer)
#     for i in g.indices[v]:g.indices[v+1]-1
#         u1 = g.nlist[i]
#         for j in (i+1):g.indices[v+1]-1
#             u2 = g.nlist[j]
#             add_edge!(g, u1.v, u2.v, u2.ofs .- u1.ofs)
#         end
#     end
#     nothing
# end

# function LightGraphs.SimpleGraphs.rem_vertex!(g::PeriodicGraph3D, v::Integer)
#     n = g.indices[v+1] - g.indices[v]
#     for i in g.indices[v]:g.indices[v+1]-1
#         u = g.nlist[i]
#         j = search_dst_between(g.nlist, g.indices[u.v]-i+1, g.indices[u.v+1]-i+1, v, .-u.ofs)
#         @assert g.nlist[j] == PeriodicVertex(v, .-u.ofs)
#         deleteat!(g.nlist, j) # et corriger les voisins d'après
#         throw("À compléter")
#     end
#
#     deleteat!(g.nlist, g.indices[v]:g.indices[v+1]-1)
#     deleteat!(g.indices, v)
#     for i in v:length(g.indices)
#         g.indices[i] -= n
#     end
#     throw("À compléter")
# end

function LightGraphs.SimpleGraphs.rem_vertices!(g::PeriodicGraph{N}, t::AbstractVector{<:Integer}, keep_order::Bool=false) where N
    isempty(t) && return collect(1:nv(g))
    sort!(t)
    (first(t) < 1 || last(t) > nv(g)) && throw(ArgumentError("Vertices to be removed must be in the range 1:nv(g)."))

    bt = falses(nv(g))
    bt[t] .= true

    for i in t
        g.ne[] -= count((i < x.v || (i == x.v && x.ofs > zero(x.ofs))) for x in g.nlist[i])
    end

    vmap = Int[]
    rev_vmap = zeros(Int, nv(g))

    if keep_order
        append!(vmap, collect(1:nv(g)))
        deleteat!(vmap, bt)
        rev_vmap[vmap] .= 1:length(vmap)
        deleteat!(g.nlist, bt)
    else
        i_next_vertex_to_del = 1
        next_vertex_to_del = t[i_next_vertex_to_del]
        t_end = length(t)
        g_end = length(g.nlist)
        i = 1
        while i <= g_end
            if i == next_vertex_to_del
                i_next_vertex_to_del += 1
                if i_next_vertex_to_del > t_end
                    next_vertex_to_del = nv(g) + 1
                else
                    next_vertex_to_del = t[i_next_vertex_to_del]
                end
                while t_end >= i_next_vertex_to_del && g_end == t[t_end]
                    t_end -= 1
                    g_end -= 1
                end
                if i < g_end
                    g.nlist[i], g.nlist[g_end] = g.nlist[g_end], g.nlist[i]
                    push!(vmap, g_end)
                    rev_vmap[g_end] = length(vmap)
                end
                g_end -= 1
            else
                push!(vmap, i)
                rev_vmap[i] = length(vmap)
            end
            i += 1
        end
    end
    resize!(g.nlist, g_end)


    for i in vertices(g)
        neighbors = g.nlist[i]
        remove_edges = falses(length(neighbors))
        startoffset = 1
        for (k, x) in enumerate(neighbors)
            if bt[x.v]
                remove_edges[k] = true
            else
                neigh = PeriodicVertex{N}(rev_vmap[x.v], x.ofs)
                neighbors[k] = neigh
                startoffset += isindirectedge(PeriodicEdge{N}(i, neigh))
            end
        end
        deleteat!(neighbors, remove_edges)
        sort!(neighbors)
        g.directedgestart[i] = startoffset
    end

    g.width[] = -1

    return vmap
end

function LightGraphs.SimpleGraphs.rem_vertex!(g::PeriodicGraph, v::Integer)
    n = nv(g)
    return length(rem_vertices!(g, [v])) == n - 1
end



struct PeriodicEdgeIter{N} <: AbstractEdgeIter
    g::PeriodicGraph{N}
end

eltype(::Type{PeriodicEdgeIter{N}}) where {N} = PeriodicEdge{N}
length(iter::PeriodicEdgeIter) = iter.g.ne[]

function iterate(iter::PeriodicEdgeIter{N}, (vertex, neigh)) where N
    nlists = iter.g.nlist
    n = length(nlists)
    #=@inbounds=# while vertex <= n
        if iszero(neigh)
            neigh = iter.g.directedgestart[vertex]
        end
        neighbors = nlists[vertex]
        if neigh > length(neighbors)
            vertex += 1
            neigh = 0
            continue
        end
        return (unsafe_edge{N}(vertex, neighbors[neigh]), (vertex, neigh+1))
    end
    return nothing
end

iterate(iter::PeriodicEdgeIter) = iterate(iter, (1, 0))

function in(edge, iter::PeriodicEdgeIter{N}) where N
    has_edge(iter.g, edge)
end

function cmp(it1::PeriodicEdgeIter{N}, it2::PeriodicEdgeIter{N}) where N
    n = length(it1)
    m = length(it2)
    n == m || return cmp(n,m)
    n == 0 && return 0
    e1::PeriodicEdge{N}, st1::Tuple{Int,Int} = iterate(it1)
    e2::PeriodicEdge{N}, st2::Tuple{Int,Int} = iterate(it2)
    for _ in 1:n-1
        c = cmp(e1, e2)
        if iszero(c)
            e1, st1 = iterate(it1, st1)
            e2, st2 = iterate(it2, st2)
            continue
        end
        return c
    end
    return cmp(e1, e2)
end

function isless(it1::PeriodicEdgeIter{N}, it2::PeriodicEdgeIter{N}) where N
    return cmp(it1, it2) < 0
end
function ==(it1::PeriodicEdgeIter{N}, it2::PeriodicEdgeIter{N}) where N
    return iszero(cmp(it1, it2))
end

function vertex_permutation(g::PeriodicGraph{N}, vlist) where N
    n = length(vlist)
    newvid = Vector{Int}(undef, n)
    for i in 1:n
        newvid[vlist[i]] = i
    end
    edges = Vector{Vector{PeriodicVertex{N}}}(undef, n)
    startoffsets = [1 for _ in 1:n]
    #=@inbounds=# for i in 1:n
        startoffset = 1
        neighs = copy(g.nlist[vlist[i]])
        edges[i] = neighs
        for j in 1:length(neighs)
            dst = neighs[j]
            neigh = PeriodicVertex{N}(newvid[dst.v], dst.ofs)
            neighs[j] = neigh
            startoffsets[i] += isindirectedge(PeriodicEdge{N}(i, neigh))
        end
        sort!(neighs)
    end
    return PeriodicGraph{N}(g.ne[], edges, startoffsets)
end

function LightGraphs.induced_subgraph(g::PeriodicGraph{N}, vlist::AbstractVector{U}) where {N, U<:Integer}
    allunique(vlist) || __throw_unique_vlist()
    n = length(vlist)
    n == nv(g) && return (vertex_permutation(g, vlist), vlist)
    newvid = zeros(Int, nv(g))
    for i in 1:n
        newvid[vlist[i]] = i
    end

    ne = 0
    edges = Vector{Vector{PeriodicVertex{N}}}(undef, n)
    startoffsets = [1 for _ in 1:n]
    for i in 1:n
        edges[i] = PeriodicVertex{N}[]
        startne = ne
        for dst in g.nlist[vlist[i]]
            v = newvid[dst.v]
            iszero(v) && continue
            neigh = PeriodicVertex{N}(v, dst.ofs)
            push!(edges[i], neigh)
            ne += isindirectedge(PeriodicEdge{N}(i, neigh))
        end
        startoffsets[i] = 1 + ne - startne
        sort!(edges[i])

    end
    return (PeriodicGraph{N}(ne, edges, startoffsets), vlist)
end
@noinline __throw_unique_vlist() = throw(ArgumentError("Vertices in subgraph list must be unique"))

function offset_representatives!(g::PeriodicGraph{N}, offsets) where N
    n = nv(g)
    length(offsets) == n || __throw_invalid_offsets()
    for i in 1:n
        neighs = g.nlist[i]
        ofsi = offsets[i]
        startoffset = 1
        for j in 1:length(neighs)
            x = neighs[j]
            neigh = PeriodicVertex{N}(x.v, x.ofs .+ ofsi .- offsets[x.v])
            neighs[j] = neigh
            startoffset += isindirectedge(PeriodicEdge{N}(i, neigh))
        end
        sort!(neighs)
        g.directedgestart[i] = startoffset
        g.width[] = -1
    end
    g
end
@noinline __throw_invalid_offsets() = throw(ArgumentError("The size of offsets does not match the number of vertices"))

function swap_axes!(g::PeriodicGraph{N}, t::AbstractVector{<:Integer}) where N
    length(t) == N || __throw_invalid_axesswap()
    for i in vertices(g)
        neighs = g.nlist[i]
        for j in 1:length(neighs)
            x = neighs[j]
            neigh = PeriodicVertex{N}(x.v, x.ofs[t])
            neighs[j] = neigh
        end
        sort!(neighs)
        # Note: g.width[] is unchanged
    end
    g
end
@noinline __throw_invalid_axesswap() = throw(ArgumentError("The number of axes must match the dimension of the graph"))


function LightGraphs.connected_components(g::PeriodicGraph)
    nvg = nv(g)
    label = zeros(Int, nvg)

    for u in vertices(g)
        label[u] != 0 && continue
        label[u] = u
        Q = Int[]
        push!(Q, u)
        #=@inbounds=# while !isempty(Q)
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


# function LightGraphs._neighborhood(g::Union{PeriodicGraph{2}, PeriodicGraph{3}}, v::Integer, d::Real, distmx::AbstractMatrix{U}, neighborfn::Function) where {U <: Real}
#     @assert typeof(neighborfn) === typeof(outneighbors)
#     N = ndims(g)
#     Qs::Vector{Vector{Tuple{PeriodicVertex{N},U}}} = [Tuple{PeriodicVertex{N},U}[] for _ in 1:nthreads()]
#     d < zero(U) && return Tuple{PeriodicVertex{N},U}[]
#     start_vertex = PeriodicVertex{N}(v)
#
#     push!(Qs[1], (start_vertex, zero(U)) )
#     waiting = collect(2:nthreads())
#     numbusy = Atomic{Int}(1)
#     waitnow = Atomic{Bool}[Atomic{Bool}(true) for _ in 1:nthreads()]
#     waitnow[1][] = false
#     waitinglock = SpinLock()
#
#     n = nv(g)
#     minoffset = 0
#     maxoffset = 0
#     for e in edges(g)
#         min_e, max_e = extrema(e.dst.ofs)
#         minoffset = min(min_e, minoffset)
#         maxoffset = max(max_e, maxoffset)
#     end
#     seen_size = n*(max(2*(-minoffset)*d, 2*maxoffset*d+1))^N # TODO find a closer upper bound
#     seen = [Atomic{Bool}(false) for _ in 1:seen_size]
#     seen[hash_position(start_vertex, n)][] = true
#
#     @threads for i in 1:nthreads()
#         id = threadid()
#
#         while true
#             # @show i
#
#             while !iszero(numbusy[]) && atomic_or!(waitnow[id], true) end
#             iszero(numbusy[]) && break
#             # print("ACCESS for ", id, " (", i, "); NUMBUSY IS ", numbusy[]); println()
#
#             #=@inbounds=# for (src, currdist) in Qs[id]
#                 # print("ID ", id, " ACCESSING NODE ", src); println()
#                 currdist == d && continue # should be in Q but all its neighbours are too far
#                 # print("ID ", id, " ACCESSING NODE ", src, " : GRANTED"); println()
#                 for dst in outneighbors(g, src.v)
#                     dst = PeriodicVertex{N}(dst.v, dst.ofs .+ src.ofs)
#                     position = hash_position(dst, n)
#                     # print("ID ", id, "  (", dst, ")"); println()
#                     if !seen[position][] && !Base.Threads.atomic_or!(seen[position], true)
#                         distance = currdist + distmx[src.v, dst.v]
#                         # print("ID ", id, "  [", dst, "] not seen before, proceeding"); println()
#                         if distance <= d
#                             # print("ID ", id, "[", dst, "] at the right distance; lock is ", islocked(waitinglock)); println()
#                             j = isempty(waiting) ? id : lock(waitinglock) do
#                                 isempty(waiting) ? id : pop!(waiting)
#                             end
#                             push!(Qs[j], (dst, distance))
#                             if j != id
#                                 flag = atomic_and!(waitnow[j], false)
#                                 @assert flag
#                                 atomic_add!(numbusy, 1)
#                                 # @show id, j
#                             end
#                         end
#                     end
#                 end
#             end
#             atomic_sub!(numbusy, 1)
#             # print("WAITING for ", id, " (", i, "); NUMBUSY IS ", numbusy[]); println()
#         end
#     end
#     return collect(Iterators.flatten(Qs))
# end

"""
Sets the width field of the graph so that the n-th neighbor or any vertex v of
the unit cell (0,0,0) is in a unit cell (i,j,k) such that
max∘abs.((i,k,k)) ≤ 1 + fld((n - 1), width)
"""
function graph_width!(g::PeriodicGraph{N}) where N
    distances = LightGraphs.Parallel.floyd_warshall_shortest_paths(cellgraph(g)).dists
    extremalpoints = NTuple{N,NTuple{2,Vector{Tuple{Int,Int}}}}((([],[]),([],[]),([],[])))
    # a, x ∈ extremalpoints[i][j] where i ∈ ⟦1,N⟧ and j ∈ ⟦1,2⟧ means that
    # vertex x has a neighbor whose offset is a*(-1)^(j-1) along dimension i
    maxa = 1
    for e in edges(g)
        iszero(ofs(e)) && continue
        offset = ofs(e)
        for i in 1:N
            iszero(offset[i]) && continue
            j = signbit(offset[i]) + 1
            a = abs(offset[i])
            if a > maxa
                maxa = a
            end
            push!(extremalpoints[i][j], (a, src(e)))
            push!(extremalpoints[i][3-j], (a, dst(e)))
        end
    end

    width = Rational(nv(g)+1)
    for i in 1:N
        @assert all((!isempty).(extremalpoints[i]))
        for (a1, x1) in extremalpoints[i][1], (a2, x2) in extremalpoints[i][2]
            dist = distances[x1, x2]
            if dist == typemax(Int)
                dist = 1
            end
            d = (dist + 1) // (a1 + a2)
            if d < width
                width = d
            end
        end
    end
    g.width[] = width == nv(g)+1 ? Rational(maxa) : width
end

function LightGraphs._neighborhood(g::Union{PeriodicGraph{2}, PeriodicGraph{3}}, v::Integer, d::Real, distmx::AbstractMatrix{U}, neighborfn::Function) where U <: Real
    @assert typeof(neighborfn) === typeof(outneighbors)
    N = ndims(g)
    Q = Tuple{PeriodicVertex{N}, U}[]
    d < zero(U) && return Q
    start_vertex = PeriodicVertex{N}(v)
    push!(Q, (start_vertex, zero(U),) )
    n = nv(g)
    width = g.width[]
    if width == -1
        width = graph_width!(g)
    end
    seen_size = n*(2*(1 + fld(d-1, width)) + 1)^N
    seen = falses(seen_size)
    seen[hash_position(start_vertex, n)] = true
    #=@inbounds=# for (src, currdist) in Q
        currdist == d && continue # should be in Q but all its neighbours are too far
        for dst in outneighbors(g, src.v)
            dst = PeriodicVertex{N}(dst.v, dst.ofs .+ src.ofs)
            position = hash_position(dst, n)
            if !seen[position]
                seen[position] = true
                distance = currdist + distmx[src.v, dst.v]
                if distance <= d
                    push!(Q, (dst, distance))
                end
            end
        end
    end
    return Q
end

function LightGraphs._neighborhood(g::PeriodicGraph{N}, v::Integer, d::Real, distmx::AbstractMatrix{U}, neighborfn::Function) where {N,U <: Real}
    @assert typeof(neighborfn) === typeof(outneighbors)
    Q = Tuple{PeriodicVertex, U}[]
    d < zero(U) && return Q
    start_vertex = PeriodicVertex{N}(v)
    push!(Q, (start_vertex, zero(U),) )
    n = nv(g)
    seen = Set{PeriodicVertex{N}}()
    push!(seen, start_vertex)
    #=@inbounds=# for (src, currdist) in Q
        currdist == d && continue # should be in Q but all its neighbours are too far
        @simd for dst in outneighbors(g, src.v)
            dst = PeriodicVertex(dst.v, dst.ofs .+ src.ofs)
            if dst ∉ seen
                push!(seen, dst)
                distance = currdist + distmx[src.v, dst.v]
                if distance <= d
                    push!(Q, (dst, distance))
                end
            end
        end
    end
    return Q
end

function vertex_sequence(g::PeriodicGraph, v::Integer, dmax)
    Q = LightGraphs._neighborhood(g, v, dmax, weights(g), outneighbors)
    popfirst!(Q)
    ret = zeros(Int, dmax)
    for (_, d) in Q
        ret[d] += 1
    end
    return ret
end

function cellgraph(g::PeriodicGraph)
    edgs = [Edge{Int}(src(x), dst(x)) for x in edges(g) if iszero(ofs(x))]
    ret = SimpleGraph(edgs)
    add_vertices!(ret, nv(g) - nv(ret))
    return ret
end

function periodiccellgraph(g::PeriodicGraph)
    ret = SimpleGraph([Edge{Int}(i, j.v) for i in vertices(g) for j in outneighbors(g, i)])
    add_vertices!(ret, nv(g) - nv(ret))
    return ret
end

end
