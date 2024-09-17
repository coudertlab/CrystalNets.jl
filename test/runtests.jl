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

const safeARCHIVE = deepcopy(CrystalNets.CRYSTALNETS_ARCHIVE)
const safeREVERSE = deepcopy(CrystalNets.REVERSE_CRYSTALNETS_ARCHIVE)
function __reset_archive!(safeARCHIVE, safeREVERSE)
    empty!(CrystalNets.CRYSTALNETS_ARCHIVE)
    empty!(CrystalNets.REVERSE_CRYSTALNETS_ARCHIVE)
    merge!(CrystalNets.CRYSTALNETS_ARCHIVE, safeARCHIVE)
    merge!(CrystalNets.REVERSE_CRYSTALNETS_ARCHIVE, safeREVERSE)
    nothing
end

function extract1(x)
    topo, nfold = only(x)
    @test nfold == 1
    topo
end

import CrystalNets.Clustering: SingleNodes, AllNodes, Standard, PE, PEM

@testset "MOF examples" begin
    cifs, crystalnetsdir = _finddirs()
    kwargs = (; structure=StructureType.MOF, clusterings=[Clustering.Auto,Standard,PE,PEM])
    println(stderr, "The following warning about altering performance, the two messages on warning toggling and the three error statements are expected.")
    mofdataset = determine_topology_dataset(joinpath(cifs, "MOFs"), save=false, showprogress=false; kwargs...)

    @testset "Dataset analysis" begin
        @test length(mofdataset) == 15

        hkust1 = extract1(mofdataset["HKUST-1.cif"])
        @test hkust1[SingleNodes] == hkust1[AllNodes] == hkust1[Standard]
        @test hkust1[SingleNodes].name == "tbo"

        jxust1, nfoldjxust1 = only(mofdataset["JXUST-1.cif"])
        @test nfoldjxust1 == 2
        @test jxust1[SingleNodes] == jxust1[AllNodes]
        @test PeriodicGraph(jxust1[SingleNodes]) == PeriodicGraph(REVERSE_CRYSTALNETS_ARCHIVE["pcu"])

        mil53 = extract1(mofdataset["MIL-53.cif"])
        @test string(mil53[SingleNodes]) == mil53[Standard].name == "bpq"
        @test mil53[AllNodes].name == "rna"
        @test mil53[PE].name == "sra, ABW"

        mil100 = extract1(mofdataset["MIL-100.cif"])
        @test mil100[SingleNodes].name == mil100[AllNodes].name == "moo"

        mil101 = extract1(mofdataset["MIL-101.cif"])
        @test mil101[SingleNodes].name == mil101[AllNodes].name == "mtn-e"
        @test mil101[PE].name == "mtn-e-a"

        mof5 = extract1(mofdataset["MOF-5.cif"])
        @test first(mof5) == ([AllNodes, SingleNodes, Standard, PEM] => parse(TopologicalGenome, "tbo"))
        @test last(collect(mof5)) == ([PE] => last(mof5))
        @test mof5[SingleNodes].name == mof5[AllNodes].name == mof5[Standard].name == "tbo"

        mof11 = extract1(mofdataset["MOF-11.cif"])
        @test mof11[PE].name == "pts-f"
        @test_broken first(mof11) == ([AllNodes, SingleNodes, Standard, PEM] => parse(TopologicalGenome, "pts"))

        topologies_mof14 = mofdataset["MOF-14.cif"]
        @test topologies_mof14[1] == topologies_mof14[2]
        mof14, nfoldmof14 = topologies_mof14[1]
        @test nfoldmof14 == 1
        @test mof14[SingleNodes].name == mof14[AllNodes].name == mof14[Standard].name == "pto"
        @test mof14[PE] == parse(TopologicalGenome, "sqc11259")

        mof801 = extract1(mofdataset["MOF-801.cif"])
        @test mof801[SingleNodes].name == mof801[AllNodes].name == "fcu"
        @test mof801[Standard].name == "xbi"
        @test mof801[PE].name == "ubt"

        pcn700 = extract1(mofdataset["PCN-700.cif"])
        @test pcn700[SingleNodes].name == pcn700[AllNodes].name == "bcu"
        @test pcn700[PE].name == "pcb, ACO"

        uio66 = extract1(mofdataset["UiO-66.cif"])
        @test uio66[SingleNodes].name == uio66[AllNodes].name == "fcu"
        @test uio66[Standard].name == "xbi"
        @test uio66[PE].name == "ubt"

        zif8 = extract1(mofdataset["ZIF-8.cif"])
        @test mofdataset["ZIF-8.cif"] == parse(InterpenetratedTopologyResult, "AllNodes, SingleNodes, Standard, PEM: sod, SOD | PE: sod-e")

        zif67 = extract1(mofdataset["ZIF-67.cif"])
        @test string(zif67) == "AllNodes, SingleNodes, Standard, PEM: sod, SOD\nPE: sod-e"

        webzek = mofdataset["WEBZEK.cif"]
        @test length(webzek) == 2
        webzekA, webzekB = webzek
        @test webzekA[2] == webzekB[2] == 1
        @test webzekA[1][:AllNodes] == webzekA[1][Clustering.SingleNodes] == parse(TopologicalGenome, "bcu")
        @test string(webzekB[1]) == "AllNodes, SingleNodes, Standard, PEM: dia\nPE: dia-a"
    end

    CrystalNets.toggle_warning(false)

    @test mofdataset["MIL-53.cif"] == determine_topology(joinpath(cifs, "MIL-53.cif"); kwargs...)
    @test_throws ArgumentError determine_topology(joinpath(cifs, "MIL-53.cif"); kwargs..., bonding=Bonding.Input)

    @test mofdataset["UiO-66.cif"] == determine_topology(joinpath(cifs, "UiO-66.cif"); kwargs...)

    juc101 = extract1(determine_topology(joinpath(cifs, "JUC-101.cif"); kwargs..., bonding=Bonding.Guess))
    @test juc101 == extract1(determine_topology(joinpath(cifs, "JUC-101.cif"); kwargs..., bonding=Bonding.Input))
    @test juc101[SingleNodes].name == "nia"
    @test juc101[AllNodes].name == "jjt"
    @test string(juc101[Standard].genome) == "3 1 2 0 0 0 1 2 0 1 1 1 3 0 0 0 1 3 0 1 0 1 4 0 0 0 1 4 0 0 1 1 5 0 0 0 1 5 0 1 1 1 6 0 0 0 1 6 0 1 0 1 7 0 0 0 1 7 0 0 1 2 3 0 0 -1 2 4 0 0 0 2 8 0 0 0 2 8 0 1 1 3 4 0 0 1 3 8 0 0 1 3 8 0 1 1 4 8 0 1 0 4 8 0 1 1 5 6 0 0 0 5 7 0 -1 0 5 8 1 0 0 5 8 1 1 1 6 7 0 -1 0 6 8 1 0 1 6 8 1 1 1 7 8 1 1 0 7 8 1 1 1"
    @test string(juc101[PE]) == "UNKNOWN 3 1 2 0 0 0 1 3 0 0 0 1 4 0 0 0 1 5 0 0 0 2 3 0 0 0 2 6 0 0 0 2 7 0 0 0 3 8 0 0 0 3 9 0 0 0 4 7 0 0 0 4 8 0 0 0 4 10 0 0 0 5 11 0 0 0 5 12 0 0 0 6 11 1 0 0 6 13 0 0 0 7 8 0 0 0 7 14 0 0 0 8 15 0 0 0 9 11 0 1 0 9 16 0 0 0 10 17 0 0 0 10 18 0 0 0 12 13 -1 1 0 12 16 -1 0 0 12 18 0 0 -1 13 16 0 -1 0 13 19 0 0 0 14 17 1 0 0 14 19 0 0 1 15 17 0 1 0 15 20 0 0 0 16 20 0 0 -1 18 19 -1 1 1 18 20 -1 0 0 19 20 0 -1 -1"
    @test juc101[PEM].genome == PeriodicGraph("3 1 2 0 0 0 1 3 0 0 0 1 4 0 0 0 1 5 0 0 0 1 6 0 0 0 1 7 0 0 0 2 4 0 0 0 2 8 0 0 0 3 4 0 0 0 3 9 0 0 0 4 6 0 0 0 4 10 0 0 0 4 11 0 0 0 5 6 0 0 0 5 12 0 0 0 6 7 0 0 0 6 10 0 0 0 6 11 0 0 0 7 13 0 0 0 8 14 0 0 0 8 15 0 0 0 9 16 0 0 0 9 17 0 0 0 10 18 0 0 0 11 19 0 0 0 12 15 -1 0 0 12 20 0 0 0 13 16 -1 0 0 13 21 0 0 0 14 22 0 0 0 14 23 0 0 0 15 18 1 0 -1 16 19 1 0 -1 17 22 -1 1 0 17 23 -1 1 0 18 24 0 0 0 19 25 0 0 0 20 22 -1 0 1 20 26 0 0 0 21 22 -2 1 1 21 26 -1 1 0 22 23 0 0 0 22 26 1 0 -1 23 24 0 0 -1 23 25 1 -1 -1 23 26 1 0 -1 24 26 1 0 0 25 26 0 1 0")

    ewetuw = extract1(determine_topology(joinpath(cifs, "EWETUW_clean.cif"); kwargs...))
    @test allunique(ewetuw)
    @test ewetuw[SingleNodes].name == "fit"
    @test unique!(sort!(degree(ewetuw[PE].genome))) == [3, 5, 6]
    @assert allunique(unique!(sort!(degree(last(x).genome))) for x in ewetuw)

    wemfif = extract1(determine_topology(joinpath(cifs, "WEMFIF_clean.cif"); kwargs...))
    @test wemfif[AllNodes] == wemfif[SingleNodes] == wemfif[Standard] == wemfif[PEM]
    @test wemfif[AllNodes].name == "dia"
    @test wemfif[PE].name == "crs"
    CrystalNets.toggle_warning(true)

    println(stderr, "The following warning about symmetry and the three warnings about 0-dimensional structures are expected.")

    # test cell minimization with collision nodes
    nott112 = extract1(determine_topology(joinpath(cifs, "NOTT-112.cif"); kwargs..., bonding=Bonding.Input))
    @test startswith(string(nott112), "AllNodes, PEM: ntt\nSingleNodes, Standard: nts\nPE: ")

    # test input bonding when different symmetric images of the same atoms have different bonds
    fowwar = extract1(determine_topology(joinpath(cifs, "FOWWAR.cif"); kwargs..., bonding=Bonding.Input, clusterings=[Clustering.Standard]))
    @test string(fowwar) == "Standard: UNKNOWN 3 1 1 0 0 1 1 2 0 0 0 1 3 0 0 0 1 4 0 0 0 2 5 0 0 0 2 6 0 0 0 3 6 0 0 0 3 7 0 0 0 4 5 1 0 0 4 7 0 -1 0 5 5 0 1 1 5 8 0 0 0 6 6 0 0 1 6 8 0 1 0 7 7 0 1 1 7 8 1 1 0"

    # Test 0-dimensional input
    calfig = extract1(determine_topology(joinpath(cifs, "CALFIG.cif"); kwargs..., clusterings=[Clustering.Auto]))
    @test string(calfig) == "0-dimensional"

    println(stderr, "The following @error about carbon cycle disorder is expected.")
    # Test carbon cycle disorder detection
    cizpos = extract1(determine_topology(joinpath(cifs, "CIZPOS.cif"); kwargs..., clusterings=[Clustering.Auto], bonding=Bonding.Guess))
    @test string(cizpos) == "AllNodes: fof\nSingleNodes: nbo"

    @test_broken determine_topology(joinpath(cifs, "c8ce00653a2.cif"); kwargs..., clusterings=[Clustering.PE], bonding=Bonding.Guess)

    bexvad = determine_topology(joinpath(cifs, "BEXVAD.cif"); kwargs..., clusterings=[Clustering.PEM], bonding=Bonding.Guess)
    @test_broken bexvad[1][1] == bexvad[2][1]
end

@testset "Other kinds of structures" begin
    cifs, crystalnetsdir = _finddirs()
    @test string(determine_topology(joinpath(cifs, "Clathrate_hydrate.cif"); Hbonds=true)) == "ict"
    @test string(determine_topology(joinpath(cifs, "Clathrate_hydrate.cif"); Hbonds=true, structure=StructureType.Guess)) == "AllNodes, SingleNodes: ict"
    @test string(determine_topology(joinpath(cifs, "Lithosite.cif"); structure=StructureType.Zeolite, ignore_atoms=(:K,))) == "-LIT"
end

@testset "Archive" begin
    @info "Checking that all known topologies are recognized (this can take a few minutes)."
    Threads.nthreads() == 1 && @info "Use multiple threads to reduce this time"
    reverse_archive = collect(CrystalNets.CRYSTALNETS_ARCHIVE)
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
                @error "$id failed (Archive)"
            end
        end
    end
    Test.get_testset().n_passed += length(reverse_archive) - failures
    @test failures == 0
end

@testset "Module" begin
    targets = ["pcu", "afy, AFY", "apc, APC", "bam", "bcf", "cdp", "cnd", "ecb", "fiv",
    "ftd", "ftj", "ins", "kgt", "mot", "moz", "muh", "pbz", "qom", "sig",
    "sma", "sod-f", "sod-h", "utj", "utp", "nts", "lth"]
    failurelock = ReentrantLock()
    failures = 0
    Threads.@threads for target in targets
        @info "Testing $target"
        graph = PeriodicGraph(CrystalNets.REVERSE_CRYSTALNETS_ARCHIVE[target])
        for k in 1:10
            superg = make_supercell(graph, clamp.(rand(0:fld(30, nv(graph)), 3), 1, 5))
            n = nv(superg)
            r = randperm(n)
            offsets = [SVector{3,Int}([rand(-3:3) for _ in 1:3]) for _ in 1:n]
            newgraph = swap_axes!(offset_representatives!(superg[r], offsets), randperm(3))
            if topological_genome(CrystalNet(newgraph)).name != target
                lock(failurelock) do
                    failures += 1
                    @error "$target failed (Module) with g = $(string(graph))"
                end
            end
        end
    end

    Test.get_testset().n_passed += length(targets) - failures
    @test failures == 0
    cifs, crystalnetsdir = _finddirs()
    @test redirect_stderr(devnull) do;
        c = parse_chemfile(joinpath(cifs, "Moganite.cif"))
        superc = make_supercell(c, (3, 1, 2))
        topological_genome(CrystalNet(c)).name == topological_genome(CrystalNet(superc)).name == "mog"
    end
end

# # The following testset is too long to be run on CI
# @testset "Full-randomization test 3D" begin
#     reverse_archive3D = Tuple{String,String}[(g, id) for (g, id) in CrystalNets.CRYSTALNETS_ARCHIVE if g[1] == '3']
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
    @test result == 0
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
    @test result == 0
    @test startswith(only(written), "AllNodes, SingleNodes: UNKNOWN 3") # Unknown topology with the input bonds

    empty!(ARGS)
    path = joinpath(cifs, "ALPO-3.1.1.37.001.cif")
    push!(ARGS, "-b", "auto", path)
    result, written = capture_out(out)
    @test result == 0
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

    #= FIXME: https://github.com/carlobaldassi/ArgParse.jl/pull/128
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
    =#

    empty!(ARGS)
    append!(ARGS, safeARGS)
    if basename(@__DIR__) != "test" # if used with include("runtests.jl")
        CrystalNets._reset_archive!()
    end
end

@testset "Unstable nets" begin
    minimize_to_unstable = PeriodicGraph("2 1 1 0 1 1 3 0 0 1 4 0 0 1 5 0 0 1 6 0 0 2 2 0 1 2 3 1 0 2 4 1 0 2 5 0 0 2 6 0 0")
    net_minimize_to_unstable = topological_genome(CrystalNet(minimize_to_unstable))
    @test string(net_minimize_to_unstable) == "unstable 2"


    mini2 = PeriodicGraph("2  1 4 0 0  1 2 0 0  1 3 0 0  2 3 0 1  4 5 0 0  4 6 0 0  5 6 1 0")
    mini3_3 = PeriodicGraph("3 1 2 0 0 0 1 3 0 0 0 1 4 0 0 0 1 7 0 0 0 2 3 0 0 1 4 5 0 0 0 4 6 0 0 0 4 7 0 0 0 5 6 0 1 0 7 8 0 0 0 7 9 0 0 0 8 9 1 0 0")
    mini3_2 = PeriodicGraph("3 1 2 0 0 0 1 3 0 0 0 1 7 0 0 0 2 3 0 0 1 4 5 0 0 0 4 6 0 0 0 4 7 0 0 0 5 6 0 1 0 7 8 0 0 0 7 9 0 0 0 8 9 1 0 0")
    small = PeriodicGraph("3 1 2 1 0 0 1 3 0 0 0 1 4 0 0 0 1 5 0 0 0 1 6 0 0 0 1 7 0 0 0 2 3 0 0 0 2 4 0 0 0 2 5 0 0 0 2 6 0 0 0 2 7 0 0 0 3 4 0 1 0 3 5 0 0 0 3 6 0 0 0 3 8 0 0 0 4 5 0 0 0 4 6 0 0 0 4 8 0 0 0 5 6 0 0 1 5 9 0 0 0 6 9 0 0 0 7 8 0 0 0 8 9 0 0 0")

    unstabletry = Union{PeriodicGraph2D,PeriodicGraph3D}[mini2, mini3_2, mini3_3, small]

    failurelock = ReentrantLock()
    failures = 0
    Threads.@threads for graph in unstabletry
        genome = topological_genome(CrystalNet(graph))
        @test !genome.unstable
        N = ndims(graph)
        for k in 1:40
            supercell = make_supercell(graph, rand(1:3, N))
            n = nv(supercell)
            r = randperm(n)
            offsets = [SVector{N,Int}([rand(-3:3) for _ in 1:N]) for _ in 1:n]
            newgraph = swap_axes!(offset_representatives!(supercell[r], offsets), randperm(N))
            if topological_genome(CrystalNet(newgraph)) != genome
                lock(failurelock) do
                    failures += 1
                    @error "Unstable graph $graph failed (Module) with g = $(string(newgraph))"
                end
                break
            end
        end
    end
    Test.get_testset().n_passed += length(unstabletry) - failures
    @test failures == 0
end

@testset "Non-CIF Files" begin
    testcase = joinpath(last(_finddirs()), "test", "Moganite.xyz")
    mogtopo = extract1(determine_topology(testcase))
    @test string(mogtopo) == "mog"
    # Test ordering of species is correct while parsing with chemfiles
    @test parse_chemfile(testcase).types == [repeat([:Si], 12); repeat([:O], 24)]
end
