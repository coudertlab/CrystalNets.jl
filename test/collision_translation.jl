using CrystalNets, Test, Random, PeriodicGraphs, StaticArrays, Graphs
import CrystalNets as CNets

@testset "Collision translation regressions" begin
    fixtures = joinpath(@__DIR__, "fixtures", "collision_translation")
    loadgraph(name) = PeriodicGraph(read(joinpath(fixtures, name), String))
    classify(g) = CNets.one_topology(topological_genome(g;
        ignore_types=true, authorize_pruning=false, trim=false, throw_error=true))

    @testset "Forced image must be available" begin
        collisions = CNets.CollisionList([
            CNets.CollisionNode(1:6, Vector{Int}[], Int[], [1:3, 4:6]),
            CNets.CollisionNode(7:12, Vector{Int}[], Int[], [1:3, 4:6]),
        ])
        direct = [9, 8, 7, 12, 11, 10, 3, 2, 1, 6, 4, -1]
        reverse = [9, 8, 7, 11, -2, 10, 3, 2, 1, 6, 5, 4]
        state = (direct, reverse, collisions, [fill(1, 6); fill(2, 6)])
        saved = (copy(direct), copy(reverse))
        @test CNets.find_necessary_image(12, state) == 4
        @test !CNets.attribute_next_available_modify!(
            Tuple{Int,Pair{Int,Int}}[], 0, 1, 12, 0, state)
        @test (direct, reverse) == saved
    end

    @testset "Index origin with fixed vertices" begin
        groups = CNets.UnderlyingNets(loadgraph("small_index_origin.pgraph");
            ignore_types=true, authorize_pruning=false, trim=false, throw_error=true)
        net = only(first(only(groups.D1)))
        _, (_, collisions) = CNets.collision_nodes(net)
        @test [c.rnge for c in collisions] == [3:4, 5:6]
        state = ([2, 1, 5, -2, -1, 3], [2, 1, 6, -2, 3, -1],
                 collisions, [3, 3, 4, 4])
        @test CNets.find_necessary_image(5, state) == 3
    end

    @testset "Small public-API crash reproducers" begin
        # These threw isperm(direct_map) and a negative-index BoundsError,
        # respectively. throw_error=true prevents an error sentinel from passing.
        for file in ("small_permutation.pgraph", "small_bounds.pgraph")
            @testset "$file" begin
                genome = classify(loadgraph(file))
                @test isempty(genome.error)
                @test ndims(genome.genome) == 1
                @test nv(genome.genome) > 0
            end
        end
    end

    function check_translation_maps(shrunk, equiv, collisions)
        D = ndims(shrunk.pge.g)
        fixed = first(first(collisions).rnge) - 1
        to_shrunk = [i + fixed for (i, node) in enumerate(collisions) for _ in node.rnge]
        shrunk_vertex(v) = v <= fixed ? v : to_shrunk[v - fixed]
        symmetry = CNets.CheckSymmetryWithCollisions(collisions)
        for pt in CNets.possible_translations(shrunk)
            pvmap = symmetry(shrunk.pge, last(pt), nothing, shrunk.types)
            pvmap isa Vector{PeriodicVertex{D}} || continue
            direct = CNets.translation_to_direct_map(pvmap, equiv, collisions, to_shrunk)
            isnothing(direct) && continue
            @test isperm(direct)
            for e in edges(equiv.pge.g)
                u, v = e.src, e.dst.v
                ofs = e.dst.ofs + pvmap[shrunk_vertex(v)].ofs - pvmap[shrunk_vertex(u)].ofs
                @test has_edge(equiv.pge.g, direct[u], PeriodicVertex{D}(direct[v], ofs))
            end
        end
    end

    recovered = PeriodicGraph[]
    for line in eachline(joinpath(fixtures, "cases.tsv"))
        startswith(line, '#') && continue
        id, dimension, file, expected = split(line, '\t')
        @testset "$id original error input" begin
            groups = CNets.UnderlyingNets(loadgraph(file);
                ignore_types=true, authorize_pruning=false, trim=false, throw_error=true)
            nets = getproperty(groups, Symbol("D", dimension))
            @test length(nets) == 1
            net = only(first(only(nets)))
            shrunk, (equiv, collisions) = CNets.collision_nodes(net)
            @test !isempty(collisions)
            for step in 1:10
                check_translation_maps(shrunk, equiv, collisions)
                reduced = CNets.find_first_valid_translation_unstable(shrunk, (equiv, collisions))
                isnothing(reduced) && break
                @test step < 10
                shrunk, (equiv, collisions) = reduced
            end
            genome = CNets.topological_genome(net)
            @test isempty(genome.error)
            @test genome.genome == PeriodicGraph(expected)
            push!(recovered, genome.genome)
            for seed in 1:4
                rng = MersenneTwister(20261109 + seed)
                native = net.pge.g
                n, D = nv(native), ndims(native)
                offsets = [SVector{D,Int}(rand(rng, -3:3, D)) for _ in 1:n]
                transformed = swap_axes!(
                    offset_representatives!(native[randperm(rng, n)], offsets), randperm(rng, D))
                other = classify(transformed)
                @test isempty(other.error)
                @test other.genome == genome.genome
            end
        end
    end
    @test length(recovered) == 8
    @test length(unique(recovered[1:5])) == 1
end
