"""
    CrystalNets

Module for automatic reckognition of crystal net topologies.
To use as an executable, run the source file in a shell:
```bash
julia --project=$(normpath(@__DIR__, "..")) $(@__FILE__)
```
Otherwise, as a module, to try to reckognize the net underlying a crystal given in a
chemical file format called FILE, the entry point is the following execution:
```julia
julia> using CrystalNets

julia> reckognize_topology(topological_genome(CrystalNet(parse_chemfile(FILE))))
```
"""
module CrystalNets

export CrystalNet, topological_genome, parse_chemfile, reckognize_topology

import LinearAlgebra: det, norm, rank
using Base.Threads

using PeriodicGraphs
import PeriodicGraphs: hash_position, change_dimension
using StaticArrays
using Graphs

import Logging
import Logging: Warn, Info, @logmsg

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
    nothing
end

include("utils.jl")
__precompile__(true)
include("types.jl") # Main internal type definitions used to represent topologies
include("input.jl") # Crystal file parsing and conversion to an internal type
include("archive.jl") # Manipulation of the topological archive
include("output.jl")
include("arithmetics.jl")
include("symmetries.jl")
include("topology.jl") # Entry point for the main algorithm
include("executable.jl") # Entry point for the argument parsing of the executable
include("precompile.jl")

end # module CrystalNets

if abspath(PROGRAM_FILE) == @__FILE__
    CrystalNets.julia_main()
end
