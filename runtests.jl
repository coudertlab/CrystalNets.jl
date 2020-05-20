include("PeriodicGraphs.jl")
using Test, .PeriodicGraphs, Random, LightGraphs
import StaticArrays: SVector

Random.seed!(12)

global graph = PeriodicGraph3D(0)

@testset "Graph construction" begin
    g::PeriodicGraph3D = PeriodicGraph3D([PeriodicEdge3D(1, PeriodicVertex3D(2, (0,0,0))),
                                          PeriodicEdge3D(2, PeriodicVertex3D(3, (0,1,0))),
                                          PeriodicEdge3D(2, PeriodicVertex3D(3, (0,-1,0))),
                                          PeriodicEdge3D(3, PeriodicVertex3D(1, (0,0,1)))])
    global graph = g
    @test neighbors(g, 1) == PeriodicVertex3D[(2, (0, 0, 0)), (3, (0, 0, -1))]
    @test_throws LoopException add_edge!(g, PeriodicEdge3D(1, 1, (0, 0, 0)))
    add_vertices!(g, 997)
    @test nv(g) == 1000
    count = ne(g)
    @test count == 4
    neigh_1 = 2
    x = rand(2:999)
    neigh_x = 0
    neigh_1000 = 0
    @progress for i in 1:10000
        _e = (rand(1:1000), PeriodicVertex3D(rand(1:1000), (rand(-1:1), rand(-1:1), rand(-1:1))))
        _e[1] == _e[2].v && _e[2].ofs == (0, 0, 0) && continue
        e = PeriodicEdge3D(_e...)
        shouldsucceed = !has_edge(g, e)
        @test shouldsucceed == !has_edge(g, reverse(e))
        success = add_edge!(g, e)
        if success != shouldsucceed
            @show e
            throw("!")
        end
        @test success == shouldsucceed
        @test has_edge(g, e)
        @test has_edge(g, reverse(e))
        count += success
        if e.src == 1 || e.dst.v == 1
            neigh_1 += success
        elseif e.src == x || e.dst.v == x
            neigh_x += success
        elseif e.src == 1000 || e.dst.v == 1000
            neigh_1000 += success
        end
    end
    @test ne(g) == count
    @test length(connected_components(g)) == 1
    @test only(connected_components(g)) == collect(1:nv(g))
    add_vertex!(g)
    @test length(connected_components(g)) == 2
    @test add_edge!(g, PeriodicEdge3D(1000, 1001, (0, 0, 0)))
    @test has_edge(g, 1000, 1001)
    @test has_edge(g, 1001, 1000)
    @test length(connected_components(g)) == 1
    @test !add_edge!(g, PeriodicEdge3D(1000, 1001, (0, 0, 0)))
    @test has_edge(g, 1000, 1001)
    @test has_edge(g, 1001, 1000)
    for e in edges(g)
        @test has_edge(g, e)
        @test has_edge(g, reverse(e))
    end
end

@testset "Neighbors" begin
    g::PeriodicGraph3D = graph
    @test add_vertices!(g, 2) == 2
    @test add_edge!(g, PeriodicEdge3D(1002, 1003, (0, 0, 1)))
    @test only(neighbors(g, 1002)) == PeriodicVertex3D(1003, (0, 0, 1))
    @test only(neighbors(g, 1003)) == PeriodicVertex3D(1002, (0, 0, -1))
    @test !has_edge(g, 1002, 1002)
    @test add_edge!(g, PeriodicEdge3D(1002, 1002, (1, 0, 0)))
    @test only(neighbors(g, 1003)) == PeriodicVertex3D(1002, (0, 0, -1))
    @test length(neighbors(g, 1002)) == 3
    @test has_edge(g, 1002, 1002)
    @test neighborhood(g, 1, 1) == pushfirst!(collect(neighbors(g, 1)), PeriodicVertex3D(1, (0, 0, 0)))
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

@testset "Periodic vertex hash" begin
    n::Int = 4
    seen = falses((2n-1)^3)
    for i in -n+1:n-1, j in -n+1:n-1, k in -n+1:n-1
        x = PeriodicGraphs.hash_position(SVector{3,Int}(i,j,k))+1
        @test !seen[x]
        seen[x] = true
    end
    @test all(seen)
    append!(seen, falses((2n+1)^3-(2n-1)^3))
    for i in (-n,n), j in -n:n, k in -n:n
        x = PeriodicGraphs.hash_position(SVector{3,Int}(i,j,k))+1
        @test !seen[x]
        seen[x] = true
    end
    for i in -n+1:n-1, j in (-n,n), k in -n:n
        x = PeriodicGraphs.hash_position(SVector{3,Int}(i,j,k))+1
        @test !seen[x]
        seen[x] = true
    end
    for i in -n+1:n-1, j in -n+1:n-1, k in (-n,n)
        x = PeriodicGraphs.hash_position(SVector{3,Int}(i,j,k))+1
        @test !seen[x]
        seen[x] = true
    end

    for i in 1:length(seen)
        if !seen[i]
            @show i
        end
    end
    @test all(seen)
    g = PeriodicGraph3D([PeriodicEdge3D(1, 1, (0, 0, 1)),
                         PeriodicEdge3D(1, 1, (1, 1, 0))])
    @test neighborhood(g, 1, 1) == PeriodicVertex3D[(1, (0, 0, 0)),
            (1, (-1, -1, 0)), (1, (0, 0, -1)), (1, (0, 0, 1)), (1, (1, 1, 0))]
    @test neighborhood(g, 1, 2) == PeriodicVertex3D[(1, (0, 0, 0)), (1, (-1, -1, 0)),
            (1, (0, 0, -1)), (1, (0, 0, 1)), (1, (1, 1, 0)), (1, (-2, -2, 0)),
            (1, (-1, -1, -1)), (1, (-1, -1, 1)), (1, (0, 0, -2)), (1, (1, 1, -1)),
            (1, (0, 0, 2)), (1, (1, 1, 1)), (1, (2, 2, 0))]
end

@testset "Edge iteration" begin
    g::PeriodicGraph3D = PeriodicGraph3D([PeriodicEdge3D(2, 3, (-1, 0, 0)),
                                          PeriodicEdge3D(2, 3, (0, 1, 0)),
                                          PeriodicEdge3D(4, 4, (-1, 1, 0)),
                                          PeriodicEdge3D(1, 3, (0, 0, 1)),
                                          PeriodicEdge3D(3, 2, (0, 0, 0))])
    expected = PeriodicEdge3D[(1, 3, (0, 0, 1)), (2, 3, (-1, 0, 0)), 
                              (2, 3, (0, 0, 0)), (2, 3, (0, 1, 0)),
                              (4, 4, (1, -1, 0))]
    @test collect(edges(g)) == expected
    @test PeriodicEdge3D(4, 4, (1, -1, 0)) ∈ edges(g)
    add_vertex!(g)
    @test collect(edges(g)) == expected
    @test length(edges(g)) == ne(g)
end

@testset "Vertex removal" begin
    g::PeriodicGraph3D = PeriodicGraph3D([PeriodicEdge3D(1, 3, (0, 0, 1)),
                                          PeriodicEdge3D(2, 3, (0, 1, 0)),
                                          PeriodicEdge3D(2, 3, (0, -1, 0)),
                                          PeriodicEdge3D(3, 2, (0, 0, 0)),
                                          PeriodicEdge3D(4, 4, (1, 1, 0))])
    global gg = deepcopy(g)
    @test rem_vertices!(g, Int[]) == [1, 2, 3, 4]
    @test g == gg
    @test rem_vertices!(g, [2, 4], false) == [1, 3]
    @test nv(g) == 2
    @test ne(g) == 1
    @test neighbors(g, 1) == PeriodicVertex3D[(2, (0, 0, 1))]
    @test neighbors(g, 2) == PeriodicVertex3D[(1, (0, 0, -1))]
    @test rem_vertex!(gg, 2)
    @test rem_vertex!(gg, 2)
    @test g == gg
end

nothing
