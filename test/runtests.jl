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


@testset "Archive" begin
    @info "Checking that all known topologies are recognized (this can take a few minutes)."
    if Threads.nthreads() == 1
        @info "Use multiple threads to reduce this time"
    end
    tests = Dict{String,Bool}([x=>false for x in values(CrystalNets.CRYSTAL_NETS_ARCHIVE)])
    reverse_archive = collect(CrystalNets.CRYSTAL_NETS_ARCHIVE)
    Threads.@threads for (genome, id) in reverse_archive
        tests[id] = try
            topological_genome(CrystalNet(PeriodicGraph(genome))).name == id
        catch
            false
        end
    end
    for (id, b) in tests
        if !b
            @info "Failed for $id (Archive)"
        end
        @test b
    end
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
            tests[target] &= topological_genome(CrystalNet(graph)).name == target
        end
    end
    for (id, b) in tests
        if !b
            @show "Failed for $id (Module)"
        end
        @test b
    end

    cifs, crystalnetsdir = _finddirs()
    @test topological_genome(CrystalNet(redirect_stderr(devnull) do;
            parse_chemfile(joinpath(cifs, "Moganite.cif"))
          end)).name == "mog"
end


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

