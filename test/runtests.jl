using CrystalNets
using Test, Random
using PeriodicGraphs
using StaticArrays
import Base.Threads

CrystalNets.toggle_export(false)

function _finddirs()
    curr = last(splitdir(@__DIR__))
    root = curr == "CrystalNets" ? normpath(@__DIR__) : normpath(@__DIR__, "..")
    return joinpath(root, "test", "cif"), root
end

const known_unstable_nets = ("sxt", "llw-z") # special case for these known unstable nets

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


@testset "Archive" begin
    @info "Checking that all known topologies are reckognized (this can take a few minutes)."
    tests = Dict{String,Bool}([x=>false for x in values(CrystalNets.CRYSTAL_NETS_ARCHIVE)
                               if x ∉ known_unstable_nets])
    Threads.@threads for (genome, id) in collect(CrystalNets.CRYSTAL_NETS_ARCHIVE)
        if id ∈ known_unstable_nets
            @test_broken reckognize_topology(topological_genome(PeriodicGraph(CrystalNets.REVERSE_CRYSTAL_NETS_ARCHIVE[id]))) == id
            continue
        end
        tests[id] = reckognize_topology(topological_genome(PeriodicGraph(genome))) == id
    end
    for (id, b) in tests
        if !b
            @show "Failed for $id (Archive)"
        end
        @test b
    end
    #=
    failures = String[]
    Threads.@threads for (genome, id) in collect(CrystalNets.CRYSTAL_NETS_ARCHIVE)
        (id == "sxt" || id == "llw-z") && continue # special case for these known unstable nets
        if reckognize_topology(topological_genome(PeriodicGraph(genome))) != id
            push!(failures, id)
        end
    end
    if !isempty(failures)
        @show "Failed: $failures"
    end
    @test isempty(failures)
    =#
end

@testset "Module" begin
    targets = ["pcu", "afy, AFY", "apc, APC", "bam", "bcf", "cdp", "cnd", "ecb", "fiv",
    "ftd", "ftj", "ins", "kgt", "mot", "moz", "muh", "pbz", "qom", "sig",
    "sma", "sod-f", "sod-h", "utj", "utp"]
    tests = Dict{String,Bool}([x=>true for x in targets])
    Threads.@threads for target in targets
        @info "Testing $target"
        graph = PeriodicGraph(CrystalNets.REVERSE_CRYSTAL_NETS_ARCHIVE[target])
        n = PeriodicGraphs.nv(graph)
        for k in 1:50
            r = randperm(n)
            offsets = [SVector{3,Int}([rand(-3:3) for _ in 1:3]) for _ in 1:n]
            graph = swap_axes!(offset_representatives!(graph[r], offsets), randperm(3))
            tests[target] &= reckognize_topology(topological_genome(graph)) == target
        end
    end
    for (id, b) in tests
        if !b
            @show "Failed for $id (Module)"
        end
        @test b
    end

    cifs, crystalnetsdir = _finddirs()
    @test reckognize_topology(topological_genome(CrystalNet(
        redirect_stderr(devnull) do; parse_chemfile(joinpath(cifs, "Moganite.cif")) end))) == "mog"
end

@testset "Executable" begin
    cifs, crystalnetsdir = _finddirs()
    safeARGS = deepcopy(ARGS)

    out = tempname()
    empty!(ARGS)

    push!(ARGS, "-g", "3   1 2  0 0 0   1 2  0 0 1   1 2  0 1 0   1 2  1 0 0")
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
    push!(ARGS, "-a", CrystalNets.arc_location*"rcsr.arc", path)
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
    push!(ARGS, "-a", CrystalNets.arc_location*"rcsr.arc", path)
    result, written = capture_out(out)
    @test result == 1
    @test written == ["UNKNOWN"]
    __reset_archive!(safeARCHIVE, safeREVERSE)

    empty!(ARGS)
    path = joinpath(cifs, "HKUST-1.cif")
    push!(ARGS, "-c", "mof", path)
    result, written = capture_out(out)
    @test result == 0
    @test length(written) == 2
    @test last(written) == "tbo"

    empty!(ARGS)
    path = joinpath(cifs, "HKUST-1_sym.cif")
    push!(ARGS, "-c", "mof", path)
    result, written = capture_out(out)
    @test result == 0
    @test length(written) == 2
    @test last(written) == "tbo"

    empty!(ARGS)
    path = joinpath(cifs, "Diamond.cif")
    push!(ARGS, "-c", "atom", path)
    result, written = capture_out(out)
    @test result == 0
    @test written == ["dia"]

    empty!(ARGS)
    path = joinpath(cifs, "ALPO-3.10.7.163.001.cif")
    push!(ARGS, "-c", "guess", "-b", "input", path)
    result, written = capture_out(out)
    @test result == 0
    @test written == ["gme, GME"]

    empty!(ARGS)
    path = joinpath(cifs, "ALPO-3.10.7.163.001.cif")
    push!(ARGS, "-c", "guess", "-b", "chemfiles", path)
    result, written = capture_out(out)
    @test result == 0
    @test written == ["gme, GME"]

    empty!(ARGS)
    path = joinpath(cifs, "ALPO-3.1.1.37.001.cif")
    push!(ARGS, "-c", "guess", "-b", "input", path)
    result, written = capture_out(out)
    @test result == 1 # Unknown topology with the input bonds

    empty!(ARGS)
    path = joinpath(cifs, "ALPO-3.1.1.37.001.cif")
    push!(ARGS, "-b", "auto", path)
    result, written = capture_out(out)
    @test result == 0
    @test written == ["afi, AFI"]

    # Test automatic removal of solvent residues and sites with multiple atoms
    empty!(ARGS)
    path = joinpath(cifs, "ALPO-3.1.1.128.001.cif")
    push!(ARGS, path, "-b", "input")
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
    append!(ARGS, safeARGS)
    if splitdir(@__DIR__) != "test" # if used with include("runtests.jl")
        CrystalNets._reset_archive!()
    end
end

@testset "Executable help" begin
    safeARGS = deepcopy(ARGS)
    out = tempname()
    empty!(ARGS)
    push!(ARGS, "--help")
    result, written = capture_out(out)
    @test result == 0
    @test startswith(popfirst!(written), "usage: CrystalNets")
    @test !isempty(popfirst!(written))
    @test occursin("CRYSTAL_FILE", popfirst!(written))
    @test occursin("Form B", popfirst!(written))
    @test occursin("Form C", popfirst!(written))
    @test isempty(popfirst!(written))
    @test isempty(popfirst!(written))
    @test popfirst!(written) == "Automatic reckognition of crystal net topologies."

    empty!(ARGS)
    append!(ARGS, safeARGS)
end
