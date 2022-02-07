"""
    CrystalNets

Module for automatic recognition of crystal net topologies.
To use as an executable, run the source file in a shell:
```bash
julia --project=$(normpath(@__DIR__, "..")) $(@__FILE__)
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
       UnderlyingNets,
       parse_chemfile,
       topological_genome,
       recognize_topology,
       determine_topology,
       determine_topologies,
       guess_topology,
       guess_topologies,
       topologies_dataset,
       guess_dataset,
       StructureType,
       Bonding,
       Clustering

import LinearAlgebra: det, norm, rank
using Base.Threads
import Serialization

using PeriodicGraphs
import PeriodicGraphs: hash_position, change_dimension
using StaticArrays
using Graphs

import Logging

const DOWARN = Base.RefValue{Bool}(false)
const DOEXPORT = Base.RefValue{Bool}(false)

function toggle_warning(to=nothing)
    global DOWARN
    DOWARN[] = to isa Nothing ? !DOWARN[] : to
end
function toggle_export(to=nothing)
    global DOEXPORT
    DOEXPORT[] = to isa Nothing ? !DOEXPORT[] : to
end

function __init__()
    toggle_warning("--no-warn" ∉ ARGS)
    toggle_export("--no-export" ∉ ARGS)
    if abspath(PROGRAM_FILE) == @__FILE__
        global DOWARN
        if DOWARN[]
            @info "Proceed to a full installation (see README.md) for better performance."
        end
        ret_code = CrystalNets.julia_main()
        exit(ret_code)
    end
    nothing
end

include("utils.jl")
__precompile__(true)
include("types.jl") # Main internal type definitions used to represent topologies
include("input.jl") # Crystal file parsing and conversion to an internal type
include("archive.jl") # Manipulation of the topological archive
include("output.jl") # Crystal file exports
include("arithmetics.jl") # Handling of (possibly sparse) integer/rational matrices
include("symmetries.jl") # Extraction of symmetries
include("stability.jl") # Functions related to unstable nets
include("topology.jl") # Main functions of the algorithm
include("query.jl") # Entry point for the user-facing functions
include("executable.jl") # Entry point for the argument parsing of the executable
include("precompile.jl")

end # module CrystalNets

nothing
