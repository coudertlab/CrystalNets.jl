module CIFTypes

using LightGraphs
export LightGraphs
export Cell, CIF, Crystal

struct Cell
    latticesystem::Symbol
    spacegroup::String
    tablenumber::Int
    a::Float64
    b::Float64
    c::Float64
    α::Float64
    β::Float64
    γ::Float64
end

# TODO fill in equivalent positions
struct EquivalentPositions end


struct CIF
    natoms::Int
    geometry::Cell
    atoms::Vector{Symbol}
    bonds::BitMatrix
    pos::Matrix{Float64}
end


struct Crystal
    natoms::Int
    atoms::Vector{Symbol}
    bonds::BitMatrix
    pos::Matrix{Float64}
    graph::SimpleGraph
end
Crystal(c::CIF) = Crystal(c.natoms, c.atoms, c.bonds, c.pos, SimpleGraph(c.bonds))


end
