includet("./PeriodicGraphs.jl")
using Test, .PeriodicGraphs, Random, ProgressMeter

Random.seed!(12)

@testset "Edge construction and reverse" begin
    @test PeriodicEdge3D(1, 1, (0, 0, 1)) == PeriodicEdge3D(1, 1, (0, 0, -1))
    @test PeriodicEdge3D(2, 1, (1, 0, -1)) == PeriodicEdge3D(1, 2, (-1, 0, 1))
    @test reverse(PeriodicEdge3D(1, 2, (0, 0, -1))) == PeriodicEdge3D(2, 1, (0, 0, 1), false)
    @test reverse(PeriodicEdge3D(1, 1, (0, 0, -1))) == PeriodicEdge3D(1, 1, (0, 0, -1), false)
end

global graph = PeriodicGraph3D(0)

function gooindices(g::PeriodicGraph3D)
    checked_vertex = 0
    for k in 1:length(g.edges)
        e = g.edges[k]
        if e.src > checked_vertex
            for j in checked_vertex+1:e.src
                g.indices[j] == k || return false
            end
            checked_vertex = e.src
        end
    end
    for j in checked_vertex+1:length(g.indices)
        g.indices[j] == length(g.edges) + 1 || return false
    end
    return true
end

@testset "Graph construction" begin
    g::PeriodicGraph3D = PeriodicGraph3D(PeriodicEdge3D[(1, 2, (0, 0, 0)),
                                                        (2, 3, (0, 1, 0)),
                                                        (2, 3, (0, -1, 0)),
                                                        (3, 1, (0, 0, 1))])
    global graph = g
    @test neighbors(g, 1) == [(2, (0, 0, 0)), (3, (0, 0, -1))]
    @test_throws "Loops are forbidden edges : PeriodicEdge3D((1, 1, (0, 0, 0))) is invalid" add_edge!(g, 1, 1, (0, 0, 0))
    add_vertices!(g, 997)
    @test nv(g) == 1000
    count = ne(g)
    @test count == 4
    @progress for i in 1:10000
        _e = (rand(1:1000), rand(1:1000), (rand(-1:1), rand(-1:1), rand(-1:1)))
        _e[1] == _e[2] && _e[3] == (0, 0, 0) && continue
        e = PeriodicEdge3D(_e...)
        shouldsucceed = !has_edge(g, e)
        @test shouldsucceed == !has_edge(g, reverse(e))
        success = add_edge!(g, e)
        @test success == shouldsucceed
        @test has_edge(g, e)
        @test has_edge(g, reverse(e))
        @test gooindices(g)
        count += success
    end
    @test issorted(g.edges)
    for i in 1:nv(g)
        @test g.indices[i] == g.indices[i+1] || g.edges[g.indices[i]].src == i
        @test g.indices[i] == 1 || g.edges[g.indices[i]-1].src != i
    end
    @test ne(g) == count
    @test length(connected_components(g)) == 1
    @test only(connected_components(g)) == collect(1:nv(g))
    add_vertex!(g)
    @test length(connected_components(g)) == 2
    @test add_edge!(g, 1000, 1001, (0, 0, 0))
    @test has_edge(g, 1000, 1001)
    @test has_edge(g, 1001, 1000)
    @test length(connected_components(g)) == 1
    @test gooindices(g)
    @test !add_edge!(g, 1000, 1001, (0, 0, 0))
    @test has_edge(g, 1000, 1001)
    @test has_edge(g, 1001, 1000)
    @test gooindices(g)
    for e in edges(g)
        @test has_edge(g, e)
        @test has_edge(g, reverse(e))
    end
    @test allunique(g.edges)
end

@testset "Neighbors" begin
    g::PeriodicGraph3D = graph
    @test add_vertices!(g, 2) == 2
    @test add_edge!(g, 1002, 1003, (0, 0, 1))
    @test only(neighbors(g, 1002)) == (1003, (0, 0, 1))
    @test only(neighbors(g, 1003)) == (1002, (0, 0, -1))
    @test !has_edge(g, 1002, 1002)
    @test add_edge!(g, 1002, 1002, (1, 0, 0))
    @test only(neighbors(g, 1003)) == (1002, (0, 0, -1))
    @test length(neighbors(g, 1002)) == 3
    @test has_edge(g, 1002, 1002)
    @test neighborhood(g, 1, 1) == pushfirst!(neighbors(g, 1), (1, (0, 0, 0)))
end

@testset "Graph reduction" begin
    g::PeriodicGraph3D = graph
    c = cellgraph(g)
    p = periodiccellgraph(g)
    @test ne(g) >= ne(p) >= ne(c)
    @test nv(g) == nv(p) == nv(c)
    @test connected_components(p) == connected_components(g)
    @test connected_components(p) >= connected_components(g)
end

nothing
