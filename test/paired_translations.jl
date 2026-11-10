using CrystalNets, Test, Random, PeriodicGraphs, StaticArrays, Graphs
import CrystalNets as CNets

@testset "Paired translation fingerprints" begin
    fixtures = joinpath(@__DIR__, "fixtures", "paired_translations")
    loadgraph(name) = PeriodicGraph(read(joinpath(fixtures, name), String))
    function genome(g)
        result = CNets.one_topology(topological_genome(g;
            ignore_types=true, authorize_pruning=false, trim=false, throw_error=true))
        @test isempty(result.error)
        @test ndims(result.genome) > 0
        @test nv(result.genome) > 0
        result.genome
    end

    @testset "Small equivalent periodic representations" begin
        examples = Set{String}()
        reductions = collect(Iterators.drop(eachline(joinpath(fixtures, "reductions.tsv")), 1))
        @test length(reductions) == 23
        for line in reductions
            _, example = split(line, '\t')
            push!(examples, example)
        end
        @test length(examples) == 12
        for example in sort!(collect(examples))
            @testset "$example" begin
                a = genome(loadgraph(example * "_a.pgraph"))
                b = genome(loadgraph(example * "_b.pgraph"))
                @test a == b
                @test genome(a) == a
            end
        end
    end

    @testset "All original source and Systre-key inputs" begin
        count = 0
        for line in eachline(joinpath(fixtures, "cases.tsv"))
            startswith(line, '#') && continue
            id, dimension, source, systre, expected = split(line, '\t')
            @testset "$id" begin
                a, b = genome(loadgraph(source)), genome(loadgraph(systre))
                @test ndims(a) == parse(Int, dimension)
                @test a == b == PeriodicGraph(expected)
                # Deterministic local RNG: do not change the rest of the test suite.
                rng = MersenneTwister(20261109 + count)
                for file in (source, systre)
                    g = loadgraph(file)
                    n, D = nv(g), ndims(g)
                    offsets = [SVector{D,Int}(rand(rng, -3:3, D)) for _ in 1:n]
                    changed = swap_axes!(
                        offset_representatives!(g[randperm(rng, n)], offsets), randperm(rng, D))
                    @test genome(changed) == a
                end
                # Exercise a doubled cell in each intrinsic periodic direction.
                for axis in 1:ndims(a)
                    factors = ntuple(i -> i == axis ? 2 : 1, ndims(a))
                    @test genome(make_supercell(a, factors)) == a
                end
                @test genome(a) == a
            end
            count += 1
        end
        @test count == 23
    end

    @testset "Original primitive-cell comparisons" begin
        count = 0
        for line in eachline(joinpath(fixtures, "primitive_pairs.tsv"))
            startswith(line, '#') && continue
            id, a, b = split(line, '\t')
            @testset "$id" begin
                @test genome(PeriodicGraph(a)) == genome(PeriodicGraph(b))
            end
            count += 1
        end
        @test count == 23
    end

    @testset "Three-dimensional paired translations" begin
        g = loadgraph("synthetic_3d.pgraph")
        expected = genome(g)
        variants = CNets.paired_translation_variants(g)
        @test !isnothing(variants)
        @test length(variants) == 8
        for variant in variants
            @test nv(variant) == nv(g)
            @test ne(variant) == ne(g)
            @test degree(variant) == degree(g)
            @test genome(variant) == expected
        end
    end

    @testset "Scope and graph preservation" begin
        ladder = loadgraph("R01_a.pgraph")
        variants = CNets.paired_translation_variants(ladder)
        @test !isnothing(variants)
        @test length(variants) == 2
        for g in variants
            @test nv(g) == nv(ladder)
            @test ne(g) == ne(ladder)
            @test degree(g) == degree(ladder)
        end
        # Quotient connectivity alone is insufficient: these pair flips are
        # independent in different cells, outside the supported finite kernel.
        independent_flips = PeriodicGraph("1 1 1 1 1 2 -1 1 2 1 2 2 1")
        @test isnothing(CNets.paired_translation_variants(independent_flips))
        @test isnothing(CNets.paired_translation_variants(PeriodicGraph("1 1 1 1")))
        @test isnothing(CNets.paired_translation_variants(PeriodicGraph("2 1 1 0 1 1 1 1 0")))
        @test genome(ladder) != genome(PeriodicGraph("1 1 1 1"))
        @test genome(PeriodicGraph("2 1 1 0 1 1 1 1 0")) !=
              genome(PeriodicGraph("2 1 1 0 1 1 1 1 0 1 1 1 1"))
    end
end
