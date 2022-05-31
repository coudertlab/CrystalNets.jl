using CrystalNets
using Test, Random
using PeriodicGraphs
using StaticArrays
using Graphs
using Combinatorics
import Base.Threads

CrystalNets.toggle_export(false)

function _finddirs()
    root = dirname(dirname(pathof(CrystalNets)))
    return joinpath(root, "test", "cif"), root
end

function capture_out(name)
    result = open(name, "w") do out
        redirect_stderr(devnull) do
            redirect_stdout(CrystalNets.julia_main, out)
        end
    end
    written = readlines(name)
    return result, written
end

const safeARCHIVE = deepcopy(CrystalNets.CRYSTAL_NETS_ARCHIVE)
const safeREVERSE = deepcopy(CrystalNets.REVERSE_CRYSTAL_NETS_ARCHIVE)
function __reset_archive!(safeARCHIVE, safeREVERSE)
    empty!(CrystalNets.CRYSTAL_NETS_ARCHIVE)
    empty!(CrystalNets.REVERSE_CRYSTAL_NETS_ARCHIVE)
    merge!(CrystalNets.CRYSTAL_NETS_ARCHIVE, safeARCHIVE)
    merge!(CrystalNets.REVERSE_CRYSTAL_NETS_ARCHIVE, safeREVERSE)
    nothing
end

import CrystalNets.Clustering: SingleNodes, AllNodes, Standard, PE, PEM


CrystalNets.toggle_export(false)

@testset "MOF examples" begin
    cifs, crystalnetsdir = _finddirs()
    kwargs = (; structure=StructureType.MOF, clusterings=[Clustering.Auto,Standard,PE,PEM])
    println(stderr, "The following warning about warnings altering performance is expected")
    mofdataset = determine_topology_dataset(joinpath(cifs, "MOFs"), false; kwargs...)

    @testset "Dataset analysis" begin
        @test length(mofdataset) == 14

        hkust1 = mofdataset["HKUST-1.cif"]
        @test hkust1[SingleNodes] == hkust1[AllNodes] == hkust1[Standard]
        @test hkust1[SingleNodes].name == "tbo"

        jxust1 = mofdataset["JXUST-1.cif"]
        @test jxust1[SingleNodes] == jxust1[AllNodes]
        @test PeriodicGraph(jxust1[SingleNodes]) == PeriodicGraph(REVERSE_CRYSTAL_NETS_ARCHIVE["pcu"])

        mil53 = mofdataset["MIL-53.cif"]
        @test string(mil53[SingleNodes]) == mil53[Standard].name == "bpq"
        @test mil53[AllNodes].name == "rna"
        @test mil53[PE].name == "sra, ABW"

        mil100 = mofdataset["MIL-100.cif"]
        @test mil100[SingleNodes].name == mil100[AllNodes].name == "moo"

        mil101 = mofdataset["MIL-101.cif"]
        @test mil101[SingleNodes].name == mil101[AllNodes].name == "mtn-e"
        @test mil101[PE].name == "mtn-e-a"

        mof5 = mofdataset["MOF-5.cif"]
        @test mof5[SingleNodes].name == mof5[AllNodes].name == mof5[Standard].name == "tbo"

        mof14 = mofdataset["MOF-14.cif/1"]
        @test mof14 == mofdataset["MOF-14.cif/2"]
        @test mof14[SingleNodes].name == mof14[AllNodes].name == mof14[Standard].name == "pto"
        @test mof14[PE] == parse(TopologicalGenome, "sqc11259")

        mof801 = mofdataset["MOF-801.cif"]
        @test mof801[SingleNodes].name == mof801[AllNodes].name == "fcu"
        @test mof801[Standard].name == "xbi"
        @test mof801[PE].name == "ubt"

        pcn700 = mofdataset["PCN-700.cif"]
        @test pcn700[SingleNodes].name == pcn700[AllNodes].name == "bcu"
        @test pcn700[PE].name == "pcb, ACO"

        uio66 = mofdataset["UiO-66.cif"]
        @test uio66[SingleNodes].name == uio66[AllNodes].name == "fcu"
        @test uio66[Standard].name == "xbi"
        @test uio66[PE].name == "ubt"

        zif8 = mofdataset["ZIF-8.cif"]
        @test mofdataset["ZIF-8.cif"] == parse(TopologyResult, "AllNodes, SingleNodes, Standard, PEM: sod, SOD | PE: sod-e")

        zif67 = mofdataset["ZIF-67.cif"]
        @test string(zif67) == "AllNodes, SingleNodes, Standard, PEM: sod, SOD\nPE: sod-e"
    end

    CrystalNets.toggle_warning(false)
    @test mofdataset["MIL-53.cif"] == determine_topology(joinpath(cifs, "MIL-53.cif"); kwargs...)
    @test_throws ArgumentError determine_topology(joinpath(cifs, "MIL-53.cif"); kwargs..., bonding=Bonding.Input)

    @test mofdataset["UiO-66.cif"] == determine_topology(joinpath(cifs, "UiO-66.cif"); kwargs...)

    juc101 = determine_topology(joinpath(cifs, "JUC-101.cif"); kwargs..., bonding=Bonding.Guess)
    @test juc101 == determine_topology(joinpath(cifs, "JUC-101.cif"); kwargs..., bonding=Bonding.Input)
    @test juc101[SingleNodes].name == "nia"
    @test juc101[AllNodes].name == "jjt"
    @test string(juc101[Standard].genome) == "3 1 2 0 0 0 1 2 0 1 1 1 3 0 0 0 1 3 0 1 0 1 4 0 0 0 1 4 0 0 1 1 5 0 0 0 1 5 0 1 1 1 6 0 0 0 1 6 0 1 0 1 7 0 0 0 1 7 0 0 1 2 3 0 0 -1 2 4 0 0 0 2 8 0 0 0 2 8 0 1 1 3 4 0 0 1 3 8 0 0 1 3 8 0 1 1 4 8 0 1 0 4 8 0 1 1 5 6 0 0 0 5 7 0 -1 0 5 8 1 0 0 5 8 1 1 1 6 7 0 -1 0 6 8 1 0 1 6 8 1 1 1 7 8 1 1 0 7 8 1 1 1"
    @test string(juc101[PE]) == "UNKNOWN 3 1 2 0 0 0 1 3 0 0 0 1 4 0 0 0 1 5 0 0 0 2 3 0 0 0 2 6 0 0 0 2 7 0 0 0 3 8 0 0 0 3 9 0 0 0 4 7 0 0 0 4 8 0 0 0 4 10 0 0 0 5 11 0 0 0 5 12 0 0 0 6 11 1 0 0 6 13 0 0 0 7 8 0 0 0 7 14 0 0 0 8 15 0 0 0 9 11 0 1 0 9 16 0 0 0 10 17 0 0 0 10 18 0 0 0 12 13 -1 1 0 12 16 -1 0 0 12 18 0 0 -1 13 16 0 -1 0 13 19 0 0 0 14 17 1 0 0 14 19 0 0 1 15 17 0 1 0 15 20 0 0 0 16 20 0 0 -1 18 19 -1 1 1 18 20 -1 0 0 19 20 0 -1 -1"
    @test juc101[PEM].genome == PeriodicGraph("3 1 2 0 0 0 1 3 0 0 0 1 4 0 0 0 1 5 0 0 0 1 6 0 0 0 1 7 0 0 0 2 4 0 0 0 2 8 0 0 0 3 4 0 0 0 3 9 0 0 0 4 6 0 0 0 4 10 0 0 0 4 11 0 0 0 5 6 0 0 0 5 12 0 0 0 6 7 0 0 0 6 10 0 0 0 6 11 0 0 0 7 13 0 0 0 8 14 0 0 0 8 15 0 0 0 9 16 0 0 0 9 17 0 0 0 10 18 0 0 0 11 19 0 0 0 12 15 -1 0 0 12 20 0 0 0 13 16 -1 0 0 13 21 0 0 0 14 22 0 0 0 14 23 0 0 0 15 18 1 0 -1 16 19 1 0 -1 17 22 -1 1 0 17 23 -1 1 0 18 24 0 0 0 19 25 0 0 0 20 22 -1 0 1 20 26 0 0 0 21 22 -2 1 1 21 26 -1 1 0 22 23 0 0 0 22 26 1 0 -1 23 24 0 0 -1 23 25 1 -1 -1 23 26 1 0 -1 24 26 1 0 0 25 26 0 1 0")

    ewetuw = determine_topology(joinpath(cifs, "EWETUW_clean.cif"); kwargs...)
    @test allunique(ewetuw)
    @test ewetuw[SingleNodes].name == "fit"
    @test unique!(sort!(degree(ewetuw[PE].genome))) == [3, 5, 6]
    @assert allunique(unique!(sort!(degree(x.genome))) for x in ewetuw)

    wemfif = determine_topology(joinpath(cifs, "WEMFIF_clean.cif"); kwargs...)
    @test wemfif[AllNodes] == wemfif[SingleNodes] == wemfif[Standard] == wemfif[PEM]
    @test wemfif[AllNodes].name == "dia"
    @test wemfif[PE].name == "crs"
    CrystalNets.toggle_warning(true)
end

@testset "Archive" begin
    @info "Checking that all known topologies are recognized (this can take a few minutes)."
    if Threads.nthreads() == 1
        @info "Use multiple threads to reduce this time"
    end
    reverse_archive = collect(CrystalNets.CRYSTAL_NETS_ARCHIVE)
    failurelock = ReentrantLock()
    failures = 0
    Threads.@threads for (genome, id) in reverse_archive
        test = try
            topological_genome(CrystalNet(PeriodicGraph(genome))).name == id
        catch e
            CrystalNets.isinterrupt(e) && rethrow()
            false
        end
        if !test
            lock(failurelock) do
                failures += 1
                # The following will throw as non-boolean, hence printing the failing test
                @test "$id failed (Archive)"
            end
        end
    end
    Test.get_testset().n_passed += length(reverse_archive) - failures
end

@testset "Module" begin
    targets = ["pcu", "afy, AFY", "apc, APC", "bam", "bcf", "cdp", "cnd", "ecb", "fiv",
    "ftd", "ftj", "ins", "kgt", "mot", "moz", "muh", "pbz", "qom", "sig",
    "sma", "sod-f", "sod-h", "utj", "utp"]
    failurelock = ReentrantLock()
    failures = 0
    Threads.@threads for target in targets
        @info "Testing $target"
        graph = PeriodicGraph(CrystalNets.REVERSE_CRYSTAL_NETS_ARCHIVE[target])
        n = PeriodicGraphs.nv(graph)
        for k in 1:50
            r = randperm(n)
            offsets = [SVector{3,Int}([rand(-3:3) for _ in 1:3]) for _ in 1:n]
            graph = swap_axes!(offset_representatives!(graph[r], offsets), randperm(3))
            if topological_genome(CrystalNet(graph)).name != target
                lock(failurelock) do
                    failures += 1
                    # The following will throw as non-boolean, hence printing the failing test
                    @test "$target failed (Module) with g = $(string(graph))"
                end
            end
        end
    end

    Test.get_testset().n_passed += length(targets) - failures
    cifs, crystalnetsdir = _finddirs()
    @test topological_genome(CrystalNet(redirect_stderr(devnull) do;
            parse_chemfile(joinpath(cifs, "Moganite.cif"))
          end)).name == "mog"
end

# # The following testset is too long to be run on CI
# @testset "Full-randomization test 3D" begin
#     reverse_archive3D = Tuple{String,String}[(g, id) for (g, id) in CrystalNets.CRYSTAL_NETS_ARCHIVE if g[1] == '3']
#     failurelock = ReentrantLock()
#     failures = 0
#     Threads.@threads for (genome, id) in reverse_archive3D
#         graph = PeriodicGraph3D(genome)
#         n = PeriodicGraphs.nv(graph)
#         for k in 1:30
#             r = randperm(n)
#             offsets = [SVector{3,Int}(rand(-3:3, 3)) for _ in 1:n]
#             graph = swap_axes!(offset_representatives!(graph[r], offsets), randperm(3))
#             topresult = topological_genome(CrystalNet(graph))
#             if topresult.name != id
#                 lock(failurelock) do
#                     failures += 1
#                     # The following will throw as non-boolean, hence printing the failing test
#                     @error """$id failed (full-random$(isempty(topresult.error) ? "" : "; error: "*
#                             Base._truncate_at_width_or_chars(true, topresult.error, 100)))"""
#                 end
#                 break
#             end
#         end
#     end
#     @info "Passed: $(length(reverse_archive3D) - failures)"
# end


@testset "Executable" begin
    cifs, crystalnetsdir = _finddirs()
    safeARGS = deepcopy(ARGS)

    out = tempname()
    empty!(ARGS)

    push!(ARGS, "-k", "3   1 2  0 0 0   1 2  0 0 1   1 2  0 1 0   1 2  1 0 0")
    result, written = capture_out(out)
    @test result == 0
    @test written == ["dia"]

    empty!(ARGS)
    path = joinpath(cifs, "ABW.cif")
    push!(ARGS, path)
    result, written = capture_out(out)
    @test result == 0
    @test written == ["sra, ABW"]

    empty!(ARGS)
    path = joinpath(cifs, "ABW.cif")
    push!(ARGS, "-a", joinpath(CrystalNets.arc_location, "rcsr.arc"), path)
    result, written = capture_out(out)
    @test result == 0
    @test written == ["sra"]
    __reset_archive!(safeARCHIVE, safeREVERSE)

    empty!(ARGS)
    path = joinpath(cifs, "RRO.cif")
    push!(ARGS, path)
    result, written = capture_out(out)
    @test result == 0
    @test written == ["RRO"]

    empty!(ARGS)
    path = joinpath(cifs, "RRO.cif")
    push!(ARGS, "-a", joinpath(CrystalNets.arc_location, "rcsr.arc"), path)
    result, written = capture_out(out)
    @test result == 1
    @test startswith(only(written), "UNKNOWN")
    __reset_archive!(safeARCHIVE, safeREVERSE)

    empty!(ARGS)
    path = joinpath(cifs, "HKUST-1.cif")
    push!(ARGS, "-s", "mof", path)
    result, written = capture_out(out)
    @test result == 0
    @test written == ["AllNodes, SingleNodes: tbo"]

    empty!(ARGS)
    path = joinpath(cifs, "HKUST-1_sym.cif")
    push!(ARGS, "-s", "mof", path)
    result, written = capture_out(out)
    @test result == 0
    @test written == ["AllNodes, SingleNodes: tbo"]

    empty!(ARGS)
    path = joinpath(cifs, "Diamond.cif")
    push!(ARGS, "-c", "atom", path)
    result, written = capture_out(out)
    @test result == 0
    @test written == ["EachVertex: dia"]

    empty!(ARGS)
    path = joinpath(cifs, "Diamond.cif")
    push!(ARGS, path)
    result, written = capture_out(out)
    @test result == 0
    @test written == ["dia"]

    empty!(ARGS)
    path = joinpath(cifs, "ALPO-3.10.7.163.001.cif")
    push!(ARGS, "-s", "guess", "-b", "input", path)
    result, written = capture_out(out)
    @test result == 0
    @test written == ["AllNodes, SingleNodes: gme, GME"]

    empty!(ARGS)
    path = joinpath(cifs, "ALPO-3.10.7.163.001.cif")
    push!(ARGS, "-s", "guess", "-b", "guess", path)
    result, written = capture_out(out)
    @test result == 0
    @test written == ["AllNodes, SingleNodes: gme, GME"]

    empty!(ARGS)
    path = joinpath(cifs, "ALPO-3.1.1.37.001.cif")
    push!(ARGS, "-s", "guess", "-b", "input", path)
    result, written = capture_out(out)
    @test result == 1 # Unknown topology with the input bonds

    empty!(ARGS)
    path = joinpath(cifs, "ALPO-3.1.1.37.001.cif")
    push!(ARGS, "-b", "auto", path)
    result, written = capture_out(out)
    @test_broken result == 0
    @test_broken written == ["afi, AFI"]

    # Test automatic removal of solvent residues and sites with multiple atoms
    empty!(ARGS)
    path = joinpath(cifs, "ALPO-3.1.1.128.001.cif")
    push!(ARGS, path)
    result, written = capture_out(out)
    @test result == 0
    @test written == ["sas, SAS"]

    empty!(ARGS)
    path = joinpath(cifs, "ALPO-3.1.2.20.002.cif")
    push!(ARGS, path)
    result, written = capture_out(out)
    @test result == 0
    @test written == ["gis, GIS"]

    empty!(ARGS)
    push!(ARGS, "--help")
    result, written = capture_out(out)
    @test result == 0
    @test startswith(popfirst!(written), "usage: CrystalNets")
    @test !isempty(popfirst!(written))
    @test occursin("CRYSTAL_FILE", popfirst!(written))
    # @test occursin("Form B", popfirst!(written))
    # @test occursin("Form C", popfirst!(written))
    @test isempty(popfirst!(written))
    @test isempty(popfirst!(written))
    @test popfirst!(written) == "Automatic recognition of crystal net topologies."

    empty!(ARGS)
    append!(ARGS, safeARGS)
    if basename(@__DIR__) != "test" # if used with include("runtests.jl")
        CrystalNets._reset_archive!()
    end
end


@testset "Unstable nets" begin
    for n in 2:4
        for m in 0:div(n*(n-1), 2)
            @info "Testing unstable node canonicalization for n = $n; m = $m"
            seen = SimpleGraph[]
            for _ in 1:500
                g = SimpleGraph(n, m; seed=-1)
                any(==(g), seen) && continue
                seencolors = Set{Vector{Int}}()
                seensubnodes = Set{Vector{Vector{Int}}}()
                for _ in 1:300
                    colors = rand(1:n, n)
                    sort!(colors) # allowed since the order of vertices is random
                    colors ∈ seencolors && continue
                    push!(seencolors, colors)
                    subnodes = [Int[] for _ in 1:n]
                    for i in 1:n
                        push!(subnodes[colors[i]], i)
                    end
                    filter!(!isempty, subnodes)
                    length(subnodes) == n && continue
                    subnodes ∈ seensubnodes && continue
                    push!(seensubnodes, subnodes)
                    sig = g[CrystalNets._order_collision(g, subnodes)]
                    for r in permutations(1:n)
                        newcolors = colors[r]
                        newsubnodes = [Int[] for _ in 1:n]
                        for i in 1:n
                            push!(newsubnodes[newcolors[i]], i)
                        end
                        filter!(!isempty, newsubnodes)
                        g2 = g[r]
                        @test sig == g2[CrystalNets._order_collision(g2, newsubnodes)]
                    end
                end
            end
        end
    end
end
