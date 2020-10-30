using CrystalNets
using Test, Random
using PeriodicGraphs
using StaticArrays

const cifs = joinpath(@__DIR__, "cif")
const crystalnetsdir = normpath(@__DIR__, "..")
const exename = `$(Base.julia_cmd()) --project=$crystalnetsdir
                 $(joinpath(crystalnetsdir, "src", "CrystalNets.jl"))`

@testset "Module" begin
    for target in ["pbc", "afy", "apc", "bam", "bcf", "cdp", "cnd", "ecb", "fiv",
                   "ftd", "ftj", "ins", "kgt", "llw-z", "mot", "moz", "muh", "pbz",
                   "qom", "sig", "sma", "sod-f", "sod-h", "sxt", "utj", "utp"]
       graph = PeriodicGraph(CrystalNets.REVERSE_CRYSTAL_NETS_ARCHIVE[target])
       n = PeriodicGraphs.nv(graph)
       for k in 1:50
           r = randperm(n)
           offsets = [SVector{3,Int}([rand(-3:3) for _ in 1:3]) for _ in 1:n]
           graph = swap_axes!(offset_representatives!(graph[r], offsets), randperm(3))
           @test reckognize_topology(topological_genome(graph)) == target
       end
   end

   @info """The following warnings about guessing bonds are expected."""
   @test reckognize_topology(topological_genome(CrystalNet(parse_chemfile(joinpath(cifs, "Moganite.cif"))))) == "mog"
end


@testset "Executable" begin
    ans = read(`$exename -g "3   1 2  0 0 0   1 2  0 0 1   1 2  0 1 0   1 2  1 0 0"`, String)
    @test ans == "dia\n"

    path = joinpath(cifs, "ABW.cif")
    ans = read(`$exename $path`, String)
    @test split(ans, "\n")[end-1] == "sra"

    path = joinpath(cifs, "RRO.cif")
    ans = read(`$exename $path`, String)
    @test split(ans, "\n")[end-1] == "Unknown topology."

    path = joinpath(cifs, "HKUST-1.cif")
    ans = read(`$exename -c mof $path`, String)
    @test split(ans, "\n")[end-1] == "tbo"

    path = joinpath(cifs, "Diamond.cif")
    ans = read(`$exename -c atom $path`, String)
    @test split(ans, "\n")[end-1] == "dia"

    help = split(read(`$exename --help`, String), "\n")[1:7]
    @test startswith(help[1], "usage: CrystalNets")
    @test occursin("CRYSTAL_FILE", help[2])
    @test occursin("Form A", help[2])
    @test occursin("Form B", help[3])
    @test occursin("Form C", help[4])
    @test isempty(help[5]) && isempty(help[6])
    @test help[7] == "Automatic reckognition of crystal net topologies."
end
