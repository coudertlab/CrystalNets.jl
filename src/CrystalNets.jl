"""
    CrystalNets

Module for automatic recognition of crystal net topologies.
To use as an executable, run the source file in a shell:
```bash
julia --project=/path/to/CrystalNets/ /path/to/CrystalNets/src/CrystalNets.jl
```
Otherwise, as a module, to try to recognize the net underlying a crystal given in a
chemical file format called FILE, the entry point is the following execution:
```julia
julia> using CrystalNets

julia> determine_topology(FILE)
```
"""
module CrystalNets

export CrystalNet,
       CrystalNet1D, CrystalNet2D, CrystalNet3D,
       UnderlyingNets,
       TopologicalGenome,
       TopologyResult,
       determine_topology,
       determine_topology_dataset,
       guess_topology,
       guess_topology_dataset,
       parse_chemfile,
       topological_genome,
       recognize_topology,
       StructureType,
       Bonding,
       Clustering,
       ClusterKinds

import LinearAlgebra: det, norm, rank, cross
using Base.Threads
import Serialization

using PeriodicGraphEmbeddings
using PeriodicGraphs
using PeriodicGraphs: hash_position, change_dimension
using StaticArrays
using Graphs

import Logging

const DOWARN = Base.RefValue{Bool}(false)
const DOEXPORT = Base.RefValue{Bool}(false)

"""
    toggle_warning(to=nothing)

Toggle warnings on (if `to == true`) or off (if `to == false`).
Without an argument, toggle on and off repeatedly at each call.
"""
function toggle_warning(to=nothing)
    global DOWARN
    DOWARN[] = to isa Nothing ? !DOWARN[] : to
end

"""
    toggle_export(to=nothing)

Toggle exports on (if `to == true`) or off (if `to == false`).
Without an argument, toggle on and off repeatedly at each call.
"""
function toggle_export(to=nothing)
    global DOEXPORT
    DOEXPORT[] = to isa Nothing ? !DOEXPORT[] : to
end

function __init__()
    toggle_warning("--no-warn" ∉ ARGS)
    toggle_export("--no-export" ∉ ARGS)
    @static if abspath(PROGRAM_FILE) == @__FILE__
        global DOWARN
        if DOWARN[]
            @info "Proceed to a full installation (see README.md) for better performance."
        end
        @show "EXECUTING JULIA MAIN"
        ret_code = CrystalNets.julia_main()
        exit(ret_code)
    end
    nothing
end

include("utils.jl")
__precompile__(true)
include("specialsolver.jl") # Exact sparse symmetric integer matrix factorization
include("options.jl") # Computation options
include("types.jl") # Main internal type definitions used to represent topologies
include("input.jl") # Crystal file parsing and conversion to an internal type
include("archive.jl") # Manipulation of the topological archive
include("output.jl") # Crystal file exports
include("arithmetics.jl") # Handling of (possibly sparse) integer/rational matrices
include("stability.jl") # Functions related to unstable nets
include("topology.jl") # Main functions of the algorithm
include("query.jl") # Entry point for the user-facing functions
include("executable.jl") # Entry point for the argument parsing of the executable
include("precompile.jl")

end # module CrystalNets

nothing
