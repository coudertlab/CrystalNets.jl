using CrystalNets
import CrystalNets as CNets
using Test, Random
using PeriodicGraphs
using StaticArrays
using Graphs
using Combinatorics
import Base.Threads

CNets.toggle_export(false)

function _finddirs()
    root = dirname(dirname(pathof(CrystalNets)))
    return joinpath(root, "test", "cif"), root
end

function capture_out(name)
    result = open(name, "w") do out
        redirect_stderr(devnull) do
            redirect_stdout(CNets.julia_main, out)
        end
    end
    written = readlines(name)
    return result, written
end

const safeARCHIVE = deepcopy(CNets.CRYSTALNETS_ARCHIVE)
const safeREVERSE = deepcopy(CNets.REVERSE_CRYSTALNETS_ARCHIVE)
function _reset_archive!(safeARCHIVE, safeREVERSE)
    empty!(CNets.CRYSTALNETS_ARCHIVE)
    empty!(CNets.REVERSE_CRYSTALNETS_ARCHIVE)
    merge!(CNets.CRYSTALNETS_ARCHIVE, safeARCHIVE)
    merge!(CNets.REVERSE_CRYSTALNETS_ARCHIVE, safeREVERSE)
    nothing
end

function extract1(x)
    topo, nfold = only(x)
    @test nfold == 1
    topo
end

function count_valid_tests(n, failures)
    # accounting the correct number of tests
    for _ in 1:(n-failures)
        @test true
    end
    @test failures == 0
end

import CrystalNets.Clustering: SingleNodes, AllNodes, Standard, PE, PEM

@testset "MOF examples" begin
    cifs, crystalnetsdir = _finddirs()
    kwargs = (; structure=StructureType.MOF, clusterings=[Clustering.Auto,Standard,PE,PEM])
    println(stderr, "The following warning about altering performance, the two messages on warning toggling and the three error statements are expected.")
    mofdataset = determine_topology_dataset(joinpath(cifs, "MOFs"), save=false, showprogress=false; kwargs...)

    @testset "Dataset analysis" begin
        @test length(mofdataset) == 16

        afixew, afixew_nfold = only(mofdataset["AFIXEW_ASR.cif"])
        @test afixew_nfold == 4
        @test afixew[AllNodes] == afixew[SingleNodes] == parse(TopologicalGenome, "bex")

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

    CNets.toggle_warning(false)

    @test mofdataset["MIL-53.cif"] == determine_topology(joinpath(cifs, "MIL-53.cif"); kwargs...)
    @test_throws ArgumentError determine_topology(joinpath(cifs, "MIL-53.cif"); kwargs..., bonding=Bonding.Input)

    @test mofdataset["UiO-66.cif"] == determine_topology(joinpath(cifs, "UiO-66.cif"); kwargs...)

    juc101 = extract1(determine_topology(joinpath(cifs, "JUC-101.cif"); kwargs..., bonding=Bonding.Guess))
    @test juc101 == extract1(determine_topology(joinpath(cifs, "JUC-101.cif"); kwargs..., bonding=Bonding.Input))
    @test juc101[SingleNodes].name == "nia"
    @test juc101[AllNodes].name == "jjt"
    @test string(juc101[Standard].genome) == "3 1 2 0 0 0 1 2 0 1 1 1 3 0 0 0 1 3 0 1 0 1 4 0 0 0 1 4 0 0 1 1 5 0 0 0 1 5 0 1 1 1 6 0 0 0 1 6 0 1 0 1 7 0 0 0 1 7 0 0 1 2 3 0 0 -1 2 4 0 0 0 2 8 0 0 0 2 8 0 1 1 3 4 0 0 1 3 8 0 0 1 3 8 0 1 1 4 8 0 1 0 4 8 0 1 1 5 6 0 0 0 5 7 0 -1 0 5 8 1 0 0 5 8 1 1 1 6 7 0 -1 0 6 8 1 0 1 6 8 1 1 1 7 8 1 1 0 7 8 1 1 1"
    @test string(juc101[PE]) == "3-jqifjikzyy (3 1 2 0 0 0 1 3 0 0 0 1 4 0 0 0 1 5 0 0 0 2 3 0 0 0 2 6 0 0 0 2 7 0 0 0 3 8 0 0 0 3 9 0 0 0 4 7 0 0 0 4 8 0 0 0 4 10 0 0 0 5 11 0 0 0 5 12 0 0 0 6 11 1 0 0 6 13 0 0 0 7 8 0 0 0 7 14 0 0 0 8 15 0 0 0 9 11 0 1 0 9 16 0 0 0 10 17 0 0 0 10 18 0 0 0 12 13 -1 1 0 12 16 -1 0 0 12 18 0 0 -1 13 16 0 -1 0 13 19 0 0 0 14 17 1 0 0 14 19 0 0 1 15 17 0 1 0 15 20 0 0 0 16 20 0 0 -1 18 19 -1 1 1 18 20 -1 0 0 19 20 0 -1 -1)"
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
    CNets.toggle_warning(true)

    println(stderr, "The following warning about symmetry and the three warnings about 0-dimensional structures are expected.")

    # test cell minimization with collision nodes
    nott112 = extract1(determine_topology(joinpath(cifs, "NOTT-112.cif"); kwargs..., bonding=Bonding.Input))
    @test startswith(string(nott112), "AllNodes, PEM: ntt\nSingleNodes, Standard: nts\nPE: 3-pxakgyqzyg")

    # test input bonding when different symmetric images of the same atoms have different bonds
    fowwar = extract1(determine_topology(joinpath(cifs, "FOWWAR.cif"); kwargs..., bonding=Bonding.Input, clusterings=[Clustering.Standard]))
    @test string(fowwar) == "Standard: 3-pmyqvetluv (3 1 1 0 0 1 1 2 0 0 0 1 3 0 0 0 1 4 0 0 0 2 5 0 0 0 2 6 0 0 0 3 6 0 0 0 3 7 0 0 0 4 5 1 0 0 4 7 0 -1 0 5 5 0 1 1 5 8 0 0 0 6 6 0 0 1 6 8 0 1 0 7 7 0 1 1 7 8 1 1 0)"

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
    reverse_archive = collect(CNets.CRYSTALNETS_ARCHIVE)
    failurelock = ReentrantLock()
    failures = 0
    Threads.@threads for (genome, id) in reverse_archive
        test = try
            topological_genome(CrystalNet(PeriodicGraph(genome))).name == id
        catch e
            CNets.isinterrupt(e) && rethrow()
            false
        end
        if !test
            lock(failurelock) do
                failures += 1
                @error "$id failed (Archive)"
            end
        end
    end
    count_valid_tests(length(reverse_archive), failures)
end

@testset "Module" begin
    targets = ["pcu", "afy, AFY", "apc, APC", "bam", "bcf", "cdp", "cnd", "ecb", "fiv",
    "ftd", "ftj", "ins", "kgt", "mot", "moz", "muh", "pbz", "qom", "sig",
    "sma", "sod-f", "sod-h", "utj", "utp", "nts", "lth", "cys", "llw-z", "nts"]
    failurelock = ReentrantLock()
    failures = 0
    Threads.@threads for target in targets
        @info "Testing $target"
        graph = PeriodicGraph(CNets.REVERSE_CRYSTALNETS_ARCHIVE[target])
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

    count_valid_tests(length(targets), failures)

    cifs, crystalnetsdir = _finddirs()
    @test redirect_stderr(devnull) do;
        c = parse_chemfile(joinpath(cifs, "Moganite.cif"))
        superc = make_supercell(c, (3, 1, 2))
        topological_genome(CrystalNet(c)).name == topological_genome(CrystalNet(superc)).name == "mog"
    end
end

# # The following testset is too long to be run on CI
# @testset "Full-randomization test 3D" begin
#     reverse_archive3D = Tuple{String,String}[(g, id) for (g, id) in CNets.CRYSTALNETS_ARCHIVE if g[1] == '3']
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
    push!(ARGS, "-a", joinpath(CNets.arc_location, "rcsr.arc"), path)
    result, written = capture_out(out)
    @test result == 0
    @test written == ["sra"]
    _reset_archive!(safeARCHIVE, safeREVERSE)

    empty!(ARGS)
    path = joinpath(cifs, "RRO.cif")
    push!(ARGS, path)
    result, written = capture_out(out)
    @test result == 0
    @test written == ["RRO"]

    empty!(ARGS)
    path = joinpath(cifs, "RRO.cif")
    push!(ARGS, "-a", joinpath(CNets.arc_location, "rcsr.arc"), path)
    result, written = capture_out(out)
    @test result == 0
    @test startswith(only(written), "3-fvyjfoutao")
    _reset_archive!(safeARCHIVE, safeREVERSE)

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
    @test startswith(only(written), "AllNodes, SingleNodes: 3-wtacqyopex") # Unknown topology with the input bonds

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
        CNets._reset_archive!()
    end
end

@testset "Unstable nets" begin
    mini2 = PeriodicGraph("2  1 4 0 0  1 2 0 0  1 3 0 0  2 3 0 1  4 5 0 0  4 6 0 0  5 6 1 0")
    mini3_2 = PeriodicGraph("3 1 2 0 0 0 1 3 0 0 0 1 7 0 0 0 2 3 0 0 1 4 5 0 0 0 4 6 0 0 0 4 7 0 0 0 5 6 0 1 0 7 8 0 0 0 7 9 0 0 0 8 9 1 0 0")
    mini3_3 = PeriodicGraph("3 1 2 0 0 0 1 3 0 0 0 1 4 0 0 0 1 7 0 0 0 2 3 0 0 1 4 5 0 0 0 4 6 0 0 0 4 7 0 0 0 5 6 0 1 0 7 8 0 0 0 7 9 0 0 0 8 9 1 0 0")
    small = PeriodicGraph("3 1 2 1 0 0 1 3 0 0 0 1 4 0 0 0 1 5 0 0 0 1 6 0 0 0 1 7 0 0 0 2 3 0 0 0 2 4 0 0 0 2 5 0 0 0 2 6 0 0 0 2 7 0 0 0 3 4 0 1 0 3 5 0 0 0 3 6 0 0 0 3 8 0 0 0 4 5 0 0 0 4 6 0 0 0 4 8 0 0 0 5 6 0 0 1 5 9 0 0 0 6 9 0 0 0 7 8 0 0 0 8 9 0 0 0")

    u1A = PeriodicGraph("1  1 1 1  1 2 0  1 3 0  2 3 0")
    u1B = PeriodicGraph("1  1 1 1  1 2 0  1 3 0  2 3 0  4 5 0  5 4 1  1 4 0  1 5 0")
    u2A = PeriodicGraph("2 1 2 0 0 1 4 -1 0 1 4 0 -1 2 3 -1 -1 2 3 0 0 3 4 0 0")
    u2B = PeriodicGraph("2 1 4 -1 0 1 4 0 0 1 5 0 -1 2 3 -1 0 2 3 0 0 2 6 0 -1 3 5 -1 1 3 5 0 1 4 6 -1 1 4 6 0 1")
    u2C = PeriodicGraph("2 1 3 -1 0 1 4 -1 1 1 4 0 -1 2 3 -1 1 2 3 0 -1 2 4 -1 0")
    u3A = PeriodicGraph("3 1 2 0 0 0 1 3 -1 -1 -1 1 3 -1 0 0 1 3 0 -1 0 1 3 0 0 1 2 3 -1 -1 0 2 3 0 0 0")
    u3B = PeriodicGraph("3 1 2 0 0 0 1 3 0 0 0 1 4 0 0 0 1 5 0 0 0 1 6 -1 0 0 1 7 0 0 0 2 3 0 0 0 2 4 0 0 0 2 5 0 -1 0 2 6 0 0 0 2 9 0 0 0 3 4 0 0 -1 3 5 0 0 0 3 6 0 0 0 3 8 0 0 0 4 5 0 0 0 4 6 0 0 0 4 8 0 0 0 5 6 0 0 0 5 9 0 0 0 6 7 0 0 0 7 9 0 0 0 8 9 0 0 1")

    gm = PeriodicGraph("2  1 1 1 0  1 2 0 0  1 3 0 0  1 4 0 0  1 2 0 1  1 3 0 1  1 4 0 1  2 3 0 0");
    gmx3 = PeriodicGraph("2  1 2 0 0  1 2 0 -1  1 3 0 0  1 3 0 -1  1 4 0 0  1 4 0 -1  1 5 0 0  1 9 -1 0
                             5 7 0 0  5 7 0 -1  5 8 0 0  5 8 0 -1  5 6 0 0  5 6 0 -1  5 9 0 0
                             9 10 0 0  9 10 0 -1  9 11 0 0  9 11 0 -1  9 12 0 0  9 12 0 -1
                             2 3 0 0  7 8 0 0  10 12 0 0");
    gmx5 = PeriodicGraph("2  1 6 0 0  1 7 0 0  1 8 0 -1  1 5 0 0  1 2 -1 0  1 3 -1 0  1 4 0 0
                             4 8 0 0  4 9 0 0  4 10 0 0  4 11 0 0  4 6 0 1  4 7 0 1  4 12 0 1
                             12 10 0 0  12 11 0 0  12 9 0 0  12 15 0 -1  12 14 0 -1  12 16 0 -1  12 13 0 0
                             13 14 0 0  13 15 0 0  13 16 0 0  13 17 0 0  13 18 0 0  13 19 0 0  13 20 0 0
                             20 2 0 0  20 3 0 0  20 1 1 1  20 5 1 0  20 18 0 1  20 17 0 1  20 19 0 1
                             2 5 1 0  7 8 0 -1  10 11 0 0  14 16 0 0  18 19 0 0");

    gn = PeriodicGraph("3  1 1 1 0 0  1 2 0 0 0  1 2 0 1 0  1 2 0 0 1  1 2 0 1 1  1 3 0 0 0  1 4 0 0 0  1 5 0 0 0  1 3 0 1 0  1 4 0 1 0  1 5 0 1 0  3 4 0 0 0  2 2 1 0 0  2 6 0 0 0  2 7 0 0 0  2 8 0 0 0  2 6 0 0 1  2 7 0 0 1  2 8 0 0 1  6 7 0 0 0");

    # Nets with invalid translations
    u1D = PeriodicGraph("1  1 2 0  2 1 1  1 3 0  3 4 0  3 5 0  4 6 0  5 6 0  4 7 0  5 7 0   2 8 0  8 9 0  8 10 0  9 10 0  9 11 0  10 12 0  11 12 0")
    u2D = PeriodicGraph("2  1 5 0 0  1 8 0 0  1 10 -1 0  1 9 -1 0  2 6 0 0  2 7 0 0  2 6 0 1  2 7 0 1  3 11 0 0  3 12 0 0  3 11 0 1  3 12 0 1  4 5 0 0  4 8 0 0  4 9 0 0  4 10 0 0  5 8 0 0  5 7 0 0  6 8 0 0  9 11 0 0  10 12 0 0  11 12 0 0")

    u3Droot = PeriodicGraph("3  1 1 0 0 1  1 1 0 1 0  1 1 1 0 0  1 2 0 0 0  2 3 0 0 0  2 4 0 0 0  3 5 0 0 0  3 6 0 0 0  4 5 0 0 0  4 6 0 0 0")
    gu3Droot = CNets.one_topology(topological_genome(u3Droot)).genome
    u3D = make_supercell(u3Droot, (2, 2, 2))
    coordinationsu3D = unique!(sort!([coordination_sequence(u3D, i, 10) for i in 1:nv(u3D)]))
    _t = rem_edge!(u3D, 3, PeriodicVertex3D(6)) &&
         rem_edge!(u3D, 4, PeriodicVertex3D(5)) &&
         add_edge!(u3D, 3, PeriodicVertex3D(4)) &&
         add_edge!(u3D, 5, PeriodicVertex3D(6))
    @test _t
    @test coordinationsu3D == unique!(sort!([coordination_sequence(u3D, i, 10) for i in 1:nv(u3D)]))
    gu3D = CNets.one_topology(topological_genome(u3D)).genome
    @test gu3D != gu3Droot
    @test nv(gu3D) == nv(u3D)

    _netgm3 = CrystalNet(make_supercell(gm, (3, 1))[[8, 11, 6, 12, 3, 2, 7, 9, 10, 4, 1, 5]])
    shrunk_net3, (equiv_net3, collisions3) = CNets.collision_nodes(_netgm3)
    t3 = last(CNets.possible_translations(shrunk_net3)[2])
    @test t3 == [1//3, 0]
    pvmap3 = CNets.CheckSymmetryWithCollisions(collisions3)(shrunk_net3.pge, t3, nothing, shrunk_net3.types)
    @test pvmap3 == PeriodicVertex2D[(2, (0,0)), (3, (0,0)), (1, (1,0)), (5, (0,0)), (6, (0,0)), (4, (1,0))]
    collisions_ranges3 = [node.rnge for node in collisions3]
    to_shrunk3 = [i+3 for (i, rnge) in enumerate(collisions_ranges3) for _ in rnge]
    direct_map3 = CNets.translation_to_direct_map(pvmap3, equiv_net3, collisions3, to_shrunk3)
    @test direct_map3 == [2, 3, 1, 7, 8, 9, 10, 11, 12, 4, 5, 6]
    transformation3 = CNets.find_transformation_matrix(t3)
    collision_offsets3 = CNets.direct_map_to_collision_offsets(direct_map3, collisions_ranges3)
    new_shrunk_net3, (new_net3, new_collision_ranges3) = CNets.reduce_unstable_net(shrunk_net3, equiv_net3, collisions3, pvmap3, transformation3, collision_offsets3)
    collision_offsets3_alt = [1, 1, 1, 1, 2, 3, 2, 1, 3, 3, 2, 1]
    new_shrunk_net3_alt, (new_net3_alt, new_collision_ranges3_alt) = CNets.reduce_unstable_net(shrunk_net3, equiv_net3, collisions3, pvmap3, transformation3, collision_offsets3)
    @test topological_genome(new_net3) == topological_genome(new_net3_alt) == topological_genome(CrystalNet(gm))

    _netgm4 = CrystalNet(make_supercell(gm, (2, 2))[[12, 5, 15, 16, 13, 1, 6, 2, 8, 3, 11, 10, 9, 7, 14, 4]])
    shrunk_net4, (equiv_net4, collisions4) = CNets.collision_nodes(_netgm4)
    collision_ranges4 = [5:7, 8:10, 11:13, 14:16]
    t4 = last(CNets.possible_translations(shrunk_net4)[1])
    @test t4 == [1//2, 1//2]
    check_symmetry4 = CNets.CheckSymmetryWithCollisions(collisions4)
    pvmap4 = check_symmetry4(shrunk_net4.pge, t4, nothing, shrunk_net4.types)
    @test pvmap4 == PeriodicVertex2D[(4, (0,0)), (3, (0,1)), (2, (1,0)), (1, (1,1)), (8, (0,1)), (7, (0,0)), (6, (1,1)), (5, (1,0))]
    collision_offsets4 = [1, 1, 1, 1, 1, 2, 3, 3, 2, 1, 2, 1, 3, 3, 1, 2]
    transformation4 = CNets.find_transformation_matrix(t4)
    @test CNets.direct_map_to_collision_offsets([4, 3, 2, 1, 10, 9, 8, 13, 11, 12, 16, 15, 14, 7, 5, 6], [5:7,8:10,11:13,14:16]) == collision_offsets4
    new_shrunk_net4, (new_net4, new_collisions4) = CNets.reduce_unstable_net(shrunk_net4, equiv_net4, collisions4, pvmap4, transformation4, collision_offsets4)

    t2 = last(CNets.possible_translations(new_shrunk_net4)[2])
    @test t2 == [0, 1//2]
    pvmap2 = CNets.CheckSymmetryWithCollisions(new_collisions4)(new_shrunk_net4.pge, t2, nothing, new_shrunk_net4.types)
    @test pvmap2 == PeriodicVertex2D[(2, (0,0)), (1, (0,1)), (4, (0,0)), (3, (0,1))]
    collision_offsets2 = [1, 1, 1, 2, 3, 1, 3, 2]
    transformation2 = CNets.find_transformation_matrix(t2)
    new_shrunk_net2, (new_net2, new_collision_ranges2) = CNets.reduce_unstable_net(new_shrunk_net4, new_net4, new_collisions4, pvmap2, transformation2, collision_offsets2)
    @test topological_genome(new_net2) == topological_genome(CrystalNet(gm))

    @test nv(CNets.one_topology(topological_genome(u1D)).genome) == 12

    unstabletry = Union{PeriodicGraph1D,PeriodicGraph2D,PeriodicGraph3D}[
        mini2, mini3_2, mini3_3, small, u1A, u1B, u2A, u2B, u2C, u3A, u3B, gm, gn, u1D, u2D, u3Droot, u3D,
    ]
    unstablematerials = Union{PeriodicGraph1D,PeriodicGraph2D,PeriodicGraph3D}[
        PeriodicGraph("3 1 2 0 0 0 1 4 0 0 -1 1 5 0 0 0 1 13 0 -1 0 1 17 -1 0 0 1 17 0 0 0 2 3 0 0 0 2 6 0 0 0 2 14 0 -1 0 2 18 -1 0 0 2 18 0 0 0 3 4 0 0 0 3 7 0 0 0 3 15 0 -1 0 3 19 -1 0 0 3 19 0 0 0 4 8 0 0 0 4 16 0 -1 0 4 20 -1 0 0 4 20 0 0 0 5 6 0 0 0 5 8 0 0 -1 5 9 0 0 0 5 21 -1 0 0 5 21 0 0 0 6 7 0 0 0 6 10 0 0 0 6 22 -1 0 0 6 22 0 0 0 7 8 0 0 0 7 11 0 0 0 7 23 -1 0 0 7 23 0 0 0 8 12 0 0 0 8 24 -1 0 0 8 24 0 0 0 9 10 0 0 0 9 12 0 0 -1 9 13 0 0 0 9 25 -1 0 0 9 25 0 0 0 10 11 0 0 0 10 14 0 0 0 10 26 -1 0 0 10 26 0 0 0 11 12 0 0 0 11 15 0 0 0 11 27 -1 0 0 11 27 0 0 0 12 16 0 0 0 12 28 -1 0 0 12 28 0 0 0 13 14 0 0 0 13 16 0 0 -1 13 29 -1 0 0 13 29 0 0 0 14 15 0 0 0 14 30 -1 0 0 14 30 0 0 0 15 16 0 0 0 15 31 -1 0 0 15 31 0 0 0 16 32 -1 0 0 16 32 0 0 0 17 18 0 0 0 17 20 0 0 -1 17 21 0 0 0 17 29 0 -1 0 18 19 0 0 0 18 22 0 0 0 18 30 0 -1 0 19 20 0 0 0 19 23 0 0 0 19 31 0 -1 0 20 24 0 0 0 20 32 0 -1 0 21 22 0 0 0 21 24 0 0 -1 21 25 0 0 0 22 23 0 0 0 22 26 0 0 0 23 24 0 0 0 23 27 0 0 0 24 28 0 0 0 25 26 0 0 0 25 28 0 0 -1 25 29 0 0 0 26 27 0 0 0 26 30 0 0 0 27 28 0 0 0 27 31 0 0 0 28 32 0 0 0 29 30 0 0 0 29 32 0 0 -1 30 31 0 0 0 31 32 0 0 0"),
        PeriodicGraph("3 1 4 0 0 -1 1 8 0 0 0 1 11 0 -1 -1 1 51 -1 0 0 1 52 -1 -1 -1 1 58 -1 -1 0 2 7 0 0 0 2 11 0 0 0 2 16 0 0 0 2 48 -1 0 0 2 52 -1 0 0 2 59 -1 0 0 3 5 0 -1 0 3 10 0 0 0 3 60 -1 0 0 4 13 0 0 0 4 14 0 -1 0 4 55 -1 0 0 5 9 0 0 0 5 12 0 0 0 5 18 0 0 0 5 53 -1 0 0 5 58 -1 0 0 6 10 0 0 0 6 19 0 0 0 6 56 -1 0 0 6 60 -1 0 0 7 14 0 0 0 7 15 0 0 0 7 54 -1 0 0 8 13 0 0 -1 8 21 0 -1 0 8 60 -1 0 0 9 15 0 0 0 9 19 0 0 0 9 56 -1 0 0 10 12 0 -1 0 10 17 0 0 0 10 22 0 0 0 10 57 -1 0 0 11 14 0 0 0 11 20 0 0 1 12 15 0 0 0 13 16 0 0 0 13 23 0 0 0 13 27 0 0 0 13 59 -1 0 0 14 17 0 1 0 14 23 0 1 0 14 24 0 0 0 15 24 0 0 0 15 25 0 0 0 16 20 0 0 1 16 28 0 0 0 17 29 0 0 0 18 20 0 0 0 18 21 0 0 0 18 34 0 0 0 19 22 0 0 0 19 25 0 0 0 19 26 0 0 0 19 31 0 0 0 20 26 0 0 0 20 30 0 0 -1 21 22 0 1 0 21 33 0 0 0 22 37 0 -1 0 23 29 0 0 0 23 38 0 -1 0 24 28 0 0 0 24 39 0 0 0 25 34 0 0 0 25 41 0 0 0 26 34 0 0 0 26 36 0 0 0 27 28 0 0 0 27 36 0 0 1 27 38 0 -1 0 28 30 0 0 0 28 32 0 0 0 28 45 0 0 0 29 31 0 0 0 29 32 0 0 0 29 35 0 -1 0 29 42 0 -1 0 30 47 0 0 0 31 37 0 -1 0 31 41 0 0 0 32 41 0 0 0 33 38 0 0 -1 33 46 0 0 0 34 43 0 0 0 34 44 0 0 0 34 50 0 0 0 35 37 0 0 0 35 39 0 0 0 36 40 0 0 0 36 51 0 0 0 37 40 0 1 0 37 44 0 0 0 37 53 0 0 0 38 42 0 0 0 38 48 0 0 0 38 52 0 0 0 39 48 0 0 0 39 49 0 0 0 40 46 0 -1 0 41 49 0 0 0 41 50 0 0 0 41 57 0 0 0 42 54 0 0 0 43 46 0 0 0 43 47 0 0 -1 44 46 0 0 0 44 56 0 0 0 45 47 0 0 0 45 55 0 0 0 46 52 0 0 -1 46 58 0 0 0 47 48 0 0 0 47 51 0 0 1 47 59 0 0 0 49 55 0 0 0 49 56 0 0 0 50 56 0 0 0 50 60 0 0 0 51 60 0 0 0 53 54 0 0 0 53 56 0 0 0 54 57 0 1 0 55 57 0 0 0 55 59 0 0 0"),
        PeriodicGraph("2 1 2 -1 0 1 2 0 -1 1 2 0 0 2 3 0 0 3 4 -1 0 3 4 0 -1 3 4 0 0"),
        PeriodicGraph("2 1 2 0 0 1 21 -1 -1 1 21 -1 0 1 21 0 -1 2 11 -1 0 2 11 0 -1 2 11 0 0 3 4 0 0 3 13 -1 0 3 13 0 -1 3 13 0 0 3 24 -1 -1 3 24 -1 0 3 24 0 -1 4 14 -1 0 4 14 0 -1 4 14 0 0 5 6 0 0 5 26 -1 -1 5 26 -1 0 5 26 0 -1 6 7 0 0 6 16 -1 0 6 16 0 -1 6 16 0 0 6 27 -1 -1 6 27 -1 0 6 27 0 -1 7 17 -1 0 7 17 0 -1 7 17 0 0 8 9 0 0 8 29 -1 -1 8 29 -1 0 8 29 0 -1 9 10 0 0 9 19 -1 0 9 19 0 -1 9 19 0 0 9 30 -1 -1 9 30 -1 0 9 30 0 -1 10 20 -1 0 10 20 0 -1 10 20 0 0 11 12 0 0 11 22 -1 0 11 22 0 -1 11 22 0 0 12 23 -1 0 12 23 0 -1 12 23 0 0 13 14 0 0 14 15 0 0 14 25 -1 0 14 25 0 -1 14 25 0 0 15 26 -1 0 15 26 0 -1 15 26 0 0 16 17 0 0 17 18 0 0 17 28 -1 0 17 28 0 -1 17 28 0 0 18 29 -1 0 18 29 0 -1 18 29 0 0 19 20 0 0 22 23 0 0 23 24 0 0 25 26 0 0 26 27 0 0 28 29 0 0 29 30 0 0"),
        PeriodicGraph("2 1 13 0 -1 1 18 -1 0 1 18 0 0 1 24 -1 0 1 24 0 0 2 8 0 -1 2 8 0 0 2 19 -1 0 2 19 0 0 3 14 0 -1 3 16 0 -1 3 25 -1 0 3 25 0 0 3 27 -1 0 3 27 0 0 3 46 -1 -1 3 46 0 -1 4 9 0 -1 4 9 0 0 4 20 -1 0 4 20 0 0 5 15 0 -1 5 26 -1 0 5 26 0 0 6 11 0 -1 6 11 0 0 6 21 -1 0 6 21 0 0 7 12 0 -1 7 12 0 0 7 23 -1 0 7 23 0 0 8 29 -1 0 8 29 0 0 8 34 -1 0 8 34 0 0 8 40 -1 0 8 40 0 0 9 30 -1 0 9 30 0 0 9 36 -1 0 9 36 0 0 9 43 -1 0 9 43 0 0 10 10 1 0 10 17 0 0 10 28 -1 0 10 28 0 0 10 35 -1 0 10 35 0 0 11 31 -1 0 11 31 0 0 11 38 -1 0 11 38 0 0 11 44 -1 0 11 44 0 0 12 32 -1 0 12 32 0 0 12 39 -1 0 12 39 0 0 12 45 -1 0 12 45 0 0 13 33 -1 0 13 33 0 0 14 33 -1 0 14 33 0 0 15 37 -1 0 15 37 0 0 15 42 -1 0 15 42 0 0 16 41 -1 0 16 41 0 0 17 22 -1 1 17 22 0 1 18 19 0 0 19 34 0 -1 19 34 0 0 20 29 0 0 20 36 0 -1 20 36 0 0 20 40 0 -1 21 30 0 0 21 38 0 -1 21 38 0 0 21 43 0 -1 22 28 0 0 22 32 0 0 22 35 0 -1 22 35 0 0 22 45 0 -1 23 31 0 0 23 39 0 -1 23 39 0 0 23 44 0 -1 24 33 0 0 25 33 0 0 26 26 0 1 26 37 0 0 27 41 0 0 33 33 0 1 41 41 0 1 41 42 0 0 41 46 0 0"),
        PeriodicGraph("2 1 7 0 -1 1 7 0 0 1 8 0 -1 1 8 0 0 1 13 -1 0 1 13 0 0 1 14 -1 0 1 14 0 0 1 20 -1 -1 1 20 -1 0 1 20 0 -1 1 20 0 0 2 7 0 -1 2 7 0 0 2 9 0 -1 2 9 0 0 2 13 -1 0 2 13 0 0 2 15 -1 0 2 15 0 0 2 21 -1 -1 2 21 -1 0 2 21 0 -1 2 21 0 0 3 8 0 -1 3 8 0 0 3 10 0 -1 3 10 0 0 3 14 -1 0 3 14 0 0 4 11 0 -1 4 11 0 0 4 16 -1 0 4 16 0 0 4 17 -1 0 4 17 0 0 5 11 0 -1 5 11 0 0 5 12 0 -1 5 12 0 0 5 17 -1 0 5 17 0 0 5 18 -1 0 5 18 0 0 5 24 -1 -1 5 24 -1 0 5 24 0 -1 5 24 0 0 6 12 0 -1 6 12 0 0 6 18 -1 0 6 18 0 0 6 19 -1 0 6 19 0 0 7 13 -1 0 7 13 -1 1 7 13 0 0 7 13 0 1 7 20 -1 0 7 20 0 0 7 21 -1 0 7 21 0 0 8 14 -1 0 8 14 -1 1 8 14 0 0 8 14 0 1 8 20 -1 0 8 20 0 0 8 22 -1 0 8 22 0 0 9 15 -1 0 9 15 -1 1 9 15 0 0 9 15 0 1 9 21 -1 0 9 21 0 0 10 22 -1 0 10 22 0 0 11 17 -1 0 11 17 -1 1 11 17 0 0 11 17 0 1 11 23 -1 0 11 23 0 0 11 24 -1 0 11 24 0 0 12 18 -1 0 12 18 -1 1 12 18 0 0 12 18 0 1 12 24 -1 0 12 24 0 0 12 25 -1 0 12 25 0 0 13 20 0 -1 13 20 0 0 13 21 0 -1 13 21 0 0 14 20 0 -1 14 20 0 0 14 22 0 -1 14 22 0 0 15 16 0 0 15 21 0 -1 15 21 0 0 16 23 0 -1 16 23 0 0 17 23 0 -1 17 23 0 0 17 24 0 -1 17 24 0 0 18 24 0 -1 18 24 0 0 18 25 0 -1 18 25 0 0 19 25 0 -1 19 25 0 0"),
        PeriodicGraph("2 1 1 1 0 1 8 0 0 1 28 -1 -1 1 28 0 -1 2 5 0 0 2 19 -1 0 2 19 0 0 3 4 0 0 3 20 -1 0 3 20 0 0 4 4 1 0 4 9 0 0 4 29 -1 -1 4 29 0 -1 5 5 1 0 5 15 0 -1 5 22 -1 0 5 22 0 0 6 7 0 0 6 21 -1 0 6 21 0 0 7 7 1 0 7 10 0 0 7 30 -1 -1 7 30 0 -1 8 25 -1 0 8 25 0 0 9 26 -1 0 9 26 0 0 10 27 -1 0 10 27 0 0 11 12 0 0 11 25 -1 0 11 25 0 0 12 12 1 0 12 16 0 0 12 23 -1 0 12 23 0 0 13 14 0 0 13 26 -1 0 13 26 0 0 14 14 1 0 14 17 0 0 14 24 -1 0 14 24 0 0 15 27 -1 0 15 27 0 0 16 20 -1 1 16 20 0 1 17 21 -1 1 17 21 0 1 18 18 1 0 18 19 0 0 19 19 1 0 20 23 0 0 21 24 0 0 22 27 0 0 25 28 0 0 26 29 0 0 27 30 0 0"),
        PeriodicGraph("2 1 22 -1 -1 1 22 0 -1 1 24 -1 0 1 24 0 0 1 42 -1 -1 2 3 0 0 2 23 -1 0 2 23 0 0 2 25 -1 0 2 25 0 0 3 24 -1 0 3 24 0 0 3 26 -1 0 3 26 0 0 4 5 0 0 4 25 -1 0 4 25 0 0 4 27 -1 0 4 27 0 0 5 5 1 0 6 7 0 0 6 27 -1 0 6 27 0 0 6 29 -1 0 6 29 0 0 7 28 -1 0 7 28 0 0 7 30 -1 0 7 30 0 0 8 9 0 0 8 29 -1 0 8 29 0 0 9 9 1 0 9 30 -1 0 9 30 0 0 10 19 0 0 10 31 -1 0 10 31 0 0 10 32 -1 0 10 32 0 0 11 18 0 0 11 32 -1 0 11 32 0 0 11 33 -1 0 11 33 0 0 12 12 1 0 12 13 0 0 13 38 -1 0 13 38 0 0 14 15 0 0 14 34 -1 0 14 34 0 0 14 35 -1 0 14 35 0 0 15 15 1 0 15 37 -1 0 15 37 0 0 16 17 0 0 16 35 -1 0 16 35 0 0 16 38 -1 0 16 38 0 0 17 36 -1 0 17 36 0 0 17 37 -1 0 17 37 0 0 18 34 -1 0 18 34 0 0 18 39 -1 0 18 39 0 0 19 39 -1 0 19 39 0 0 19 40 -1 0 19 40 0 0 20 21 0 0 20 40 -1 0 20 40 0 0 21 22 0 0 21 41 -1 0 21 41 0 0 21 42 -1 0 21 42 0 0 23 24 0 0 23 42 -1 -1 23 42 0 -1 25 26 0 0 26 27 0 0 27 28 0 0 27 31 0 0 29 30 0 0 31 40 0 0 32 39 0 0 33 34 0 0 35 37 0 0 36 38 0 0 40 41 0 0"),
        PeriodicGraph("2 1 3 -1 0 1 3 0 -1 1 3 0 0 1 5 -1 -1 1 5 -1 0 1 5 0 -1 2 4 -1 0 2 4 0 -1 2 4 0 0 4 5 -1 0 4 5 0 -1 4 5 0 0"),
        PeriodicGraph("1 1 2 0 1 3 0 1 4 0 1 5 0 1 13 -1 1 13 0 1 14 -1 1 14 0 2 6 0 2 15 -1 2 15 0 3 7 0 3 16 -1 3 16 0 4 7 0 4 17 -1 4 17 0 5 6 0 5 18 -1 5 18 0 6 19 -1 6 19 0 7 20 -1 7 20 0 7 21 -1 7 21 0 8 10 0 8 11 0 8 17 -1 8 17 0 9 11 0 9 18 -1 9 18 0 10 12 0 10 22 -1 10 22 0 11 23 -1 11 23 0 12 24 -1 12 24 0 13 15 0 13 16 0 14 17 0 14 18 0 15 19 0 16 20 0 17 21 0 17 22 0 17 23 0 18 23 0 21 24 0 22 24 0"),
        PeriodicGraph("2 1 17 -1 -1 1 17 -1 0 1 17 0 0 2 3 0 0 2 10 -1 -1 2 10 0 -1 2 10 0 0 3 11 -1 -1 3 11 0 -1 3 11 0 0 4 5 0 0 4 20 -1 -1 4 20 -1 0 4 20 0 0 5 13 -1 -1 5 13 0 -1 5 13 0 0 6 7 0 0 6 14 -1 -1 6 14 0 -1 6 14 0 0 7 21 -1 -1 7 21 -1 0 7 21 0 0 8 16 -1 -1 8 16 0 -1 8 16 0 0 9 10 0 0 9 18 -1 0 9 18 0 0 9 18 0 1 11 12 0 0 12 19 -1 0 12 19 0 0 12 19 0 1 13 14 0 0 15 16 0 0 15 22 -1 0 15 22 0 0 15 22 0 1 17 18 0 0 19 20 0 0 21 22 0 0"),
        PeriodicGraph("2 1 2 0 0 1 10 -1 0 1 10 0 -1 1 10 0 0 1 11 -1 0 1 11 0 -1 1 11 0 0 1 20 -1 -1 1 20 -1 0 1 20 0 -1 2 21 -1 -1 2 21 -1 0 2 21 0 -1 3 12 -1 0 3 12 0 -1 3 12 0 0 4 5 0 0 4 14 -1 0 4 14 0 -1 4 14 0 0 5 6 0 0 5 15 -1 0 5 15 0 -1 5 15 0 0 5 24 -1 -1 5 24 -1 0 5 24 0 -1 6 25 -1 -1 6 25 -1 0 6 25 0 -1 7 8 0 0 7 17 -1 0 7 17 0 -1 7 17 0 0 8 9 0 0 8 18 -1 0 8 18 0 -1 8 18 0 0 8 27 -1 -1 8 27 -1 0 8 27 0 -1 9 28 -1 -1 9 28 -1 0 9 28 0 -1 10 12 0 0 12 23 -1 0 12 23 0 -1 12 23 0 0 13 14 0 0 13 21 -1 0 13 21 0 -1 13 21 0 0 14 15 0 0 14 22 -1 0 14 22 0 -1 14 22 0 0 16 17 0 0 16 25 -1 0 16 25 0 -1 16 25 0 0 17 18 0 0 17 26 -1 0 17 26 0 -1 17 26 0 0 19 28 -1 0 19 28 0 -1 19 28 0 0 21 22 0 0 24 25 0 0 25 26 0 0 27 28 0 0"),
        PeriodicGraph("1 1 2 0 1 3 0 1 22 -1 1 22 0 2 23 -1 2 23 0 2 24 -1 2 24 0 3 25 -1 3 25 0 4 5 0 4 22 -1 4 22 0 4 27 -1 4 27 0 5 31 -1 5 31 0 6 9 0 6 24 -1 6 24 0 7 11 0 7 26 -1 7 26 0 8 13 0 8 28 -1 8 28 0 8 34 -1 8 34 0 9 9 1 10 15 0 10 29 -1 10 29 0 11 35 -1 11 35 0 12 16 0 12 31 -1 12 31 0 13 36 -1 13 36 0 14 19 0 14 32 -1 14 32 0 15 37 -1 15 37 0 16 33 -1 16 33 0 16 39 -1 16 39 0 17 20 0 17 36 -1 17 36 0 18 21 0 18 37 -1 18 37 0 19 38 -1 19 38 0 19 42 -1 19 42 0 20 20 1 21 40 -1 21 40 0 22 23 0 23 26 0 24 28 0 25 29 0 27 30 0 30 30 1 31 33 0 32 42 0 34 36 0 35 38 0 37 40 0 39 41 0 41 41 1"),
        PeriodicGraph("2 1 2 0 0 1 4 -1 -1 1 4 -1 0 1 4 0 -1 2 3 -1 0 2 3 0 -1 2 3 0 0"),
        PeriodicGraph("2 1 7 -1 0 1 7 0 -1 1 7 0 0 2 18 -1 -1 2 18 -1 0 2 18 0 -1 3 4 0 0 3 14 -1 -1 3 14 -1 0 3 14 0 -1 4 9 -1 0 4 9 0 -1 4 9 0 0 5 6 0 0 5 16 -1 -1 5 16 -1 0 5 16 0 -1 6 11 -1 0 6 11 0 -1 6 11 0 0 7 8 0 0 8 13 -1 0 8 13 0 -1 8 13 0 0 9 10 0 0 10 15 -1 0 10 15 0 -1 10 15 0 0 11 12 0 0 12 17 -1 0 12 17 0 -1 12 17 0 0 13 14 0 0 15 16 0 0 17 18 0 0"),
        PeriodicGraph("2 1 6 -1 -1 1 6 -1 0 1 6 0 -1 2 4 -1 0 2 4 0 -1 2 4 0 0 2 7 -1 -1 2 7 -1 0 2 7 0 -1 3 5 -1 0 3 5 0 -1 3 5 0 0 4 6 -1 0 4 6 0 -1 4 6 0 0 5 7 -1 0 5 7 0 -1 5 7 0 0"),
        PeriodicGraph("2 1 15 -1 0 1 15 0 0 1 21 -1 0 1 21 0 -1 2 5 0 0 2 10 0 -1 2 10 0 0 3 3 0 1 3 14 -1 0 3 14 0 0 4 11 0 -1 4 11 0 0 4 16 -1 0 4 16 0 0 5 24 -1 0 5 24 0 -1 6 12 0 -1 6 12 0 0 6 17 -1 0 6 17 0 0 7 18 -1 0 7 18 0 0 7 22 -1 0 7 22 0 -1 8 9 0 0 8 13 0 -1 8 13 0 0 9 23 -1 0 9 23 0 -1 10 20 -1 1 10 20 0 0 11 14 -1 1 11 14 0 0 12 16 -1 1 12 16 0 0 13 17 -1 1 13 17 0 0 15 22 0 -1 15 22 0 0 18 24 0 -1 18 24 0 0 19 20 0 0 19 23 0 -1 19 23 0 0 21 21 0 1"),
        PeriodicGraph("2 1 3 -1 -1 1 3 0 -1 1 3 0 0 2 4 -1 -1 2 4 0 -1 2 4 0 0 2 5 -1 -1 2 5 -1 0 2 5 0 0 3 5 -1 0 3 5 0 0 3 5 0 1"),
        PeriodicGraph("1 1 2 0 1 3 0 1 13 -1 1 13 0 1 14 -1 1 14 0 2 16 -1 2 16 0 3 16 -1 3 16 0 4 6 0 4 15 -1 4 15 0 5 6 0 5 15 -1 5 15 0 6 19 -1 6 19 0 6 20 -1 6 20 0 7 8 0 7 9 0 7 17 -1 7 17 0 7 18 -1 7 18 0 8 22 -1 8 22 0 9 22 -1 9 22 0 10 11 0 10 23 -1 10 23 0 11 12 0 11 24 -1 11 24 0 11 25 -1 11 25 0 12 23 -1 12 23 0 13 15 0 14 15 0 16 17 0 16 18 0 19 21 0 20 21 0 21 21 1 22 24 0 22 25 0"),
    ]

    failurelock = ReentrantLock()
    failures = 0
    nlimit = length(unstabletry)
    Nunstable = (length(unstabletry)+length(unstablematerials))
    #=Threads.@threads=# for i in 1:Nunstable
        j = i > nlimit ? i - nlimit : i
        graph = (i > nlimit ? unstablematerials : unstabletry)[j]
        listname = i > nlimit ? :unstablematerials : :unstabletry
        genome = CNets.one_topology(topological_genome(graph))
        N = ndims(graph)
        #=Threads.nthreads() == 1 &&=# (println(i, " / ", Nunstable))
        for k in 1:(i > nlimit ? 50 : 100)
            #=Threads.nthreads() == 1 &&=# print(".")
            local super, supercell, n, r, offsets, axesperm, newgraph, newgenome
            super = rand(1:(3 - (i>nlimit)), N)
            if i > nlimit
                super[rand(1:N)] += 1
            end
            supercell = make_supercell(graph, super)
            n = nv(supercell)
            r = randperm(n)
            offsets = [SVector{N,Int}([rand(-3:3) for _ in 1:N]) for _ in 1:n]
            axesperm = randperm(N)
            newgraph = swap_axes!(offset_representatives!(supercell[r], offsets), axesperm)
            newgenome = CNets.one_topology(topological_genome(newgraph))
            if newgenome != genome
                @lock failurelock begin
                    failures += 1
                    @error "Unstable graph failed ($listname[$j]): g1 = PeriodicGraph(\"$graph\"); g2 = PeriodicGraph(\"$newgraph\");\ngen1 = $genome; gen2 = $newgenome;\nsupercell = $super; offsets = $offsets; axesperm = $axesperm; r = $r"
                end
                break
            end
        end
        #=Threads.nthreads() == 1 &&=# println()
    end

    @test topological_genome(gm) == topological_genome(gmx3) == topological_genome(gmx5)

    g1 = PeriodicGraph("2 1 4 -1 0 1 4 0 0 1 5 0 -1 2 3 -1 0 2 3 0 0 2 6 0 -1 3 5 -1 1 3 5 0 1 4 6 -1 1 4 6 0 1");
    g2 = PeriodicGraph("2 1 2 2 1 1 6 3 -1 1 17 -2 1 2 9 0 -1 2 11 4 1 3 5 -7 1 3 6 -1 -2 3 17 -6 0 4 6 2 2 4 7 3 4 4 12 -2 6 5 8 4 -3 5 11 5 1 6 14 -2 2 7 13 -2 2 7 14 -3 0 7 18 -4 -3 8 10 -2 1 8 15 -4 -1 8 16 -5 4 9 12 1 2 9 15 -1 -2 9 16 -2 3 10 11 3 3 10 13 -1 3 11 12 -3 0 13 17 -4 -1 14 16 -1 2 15 18 3 -1 17 18 2 -4");
    @test topological_genome(g1) == topological_genome(g2)

    count_valid_tests(length(unstabletry), failures)
end

@testset "Non-CIF Files" begin
    testcase = joinpath(last(_finddirs()), "test", "Moganite.xyz")
    mogtopo = extract1(determine_topology(testcase))
    @test string(mogtopo) == "mog"
    # Test ordering of species is correct while parsing with chemfiles
    @test parse_chemfile(testcase).types == [repeat([:Si], 12); repeat([:O], 24)]
end
