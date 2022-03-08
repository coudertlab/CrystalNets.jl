## Types
```@docs
CrystalNets.Crystal
CrystalNets.Clusters
CrystalNets.CollisionNode
CrystalNets.Cell
CrystalNets.CIF
CrystalNets.EquivalentPosition
```

## Core topology functions
```@docs
CrystalNets.topological_key
CrystalNets.CRYSTAL_NETS_ARCHIVE
CrystalNets.minimize
CrystalNets.candidate_key
CrystalNets.possible_translations
CrystalNets.find_all_valid_translations
CrystalNets.minimal_volume_matrix
CrystalNets.reduce_with_matrix
CrystalNets.partition_by_coordination_sequence
CrystalNets.find_candidates
CrystalNets.extract_through_symmetry
CrystalNets.find_initial_candidates
CrystalNets.find_candidates_onlyneighbors
CrystalNets.find_candidates_fallback
```

## Input
```@docs
CrystalNets.parse_cif
CrystalNets.CIF(file_path::AbstractString)
CrystalNets.find_refid
Base.parse(::Type{CrystalNets.EquivalentPosition}, s::AbstractString, refid=("x", "y", "z"))
CrystalNets.check_collision
CrystalNets.fix_atom_on_a_bond!
CrystalNets.least_plausible_neighbours
CrystalNets.fix_valence!
CrystalNets.sanitize_removeatoms!
CrystalNets.remove_triangles!
CrystalNets.remove_homoatomic_bonds!
CrystalNets.sanity_checks!
```

## Crystal and CIF handling
```@docs
CrystalNets.cell_parameters
CrystalNets.periodic_distance
CrystalNets.remove_partial_occupancy
CrystalNets.prune_collisions
CrystalNets.expand_symmetry
CrystalNets.trim_monovalent
CrystalNets.trimmed_crystal
CrystalNets.trim_topology
CrystalNets.equilibrium
```

## Bond guessing
```@docs
CrystalNets.guess_bonds
CrystalNets.edges_from_bonds
```

## Clustering algorithm
```@docs
CrystalNets.find_sbus
CrystalNets.regroup_sbus
CrystalNets.regroup_paddlewheel!
CrystalNets.split_sbu!
CrystalNets.reclassify!
CrystalNets.add_to_newclass!
CrystalNets.in_small_cycles_around
CrystalNets.group_cycle
CrystalNets.collapse_clusters
CrystalNets.pem_to_pe
CrystalNets.allnodes_to_singlenodes
```

## Symmetry handling
```@docs
CrystalNets.SPACE_GROUP_HALL
CrystalNets.SPACE_GROUP_HM
CrystalNets.SPACE_GROUP_FULL
CrystalNets.SPACE_GROUP_IT
CrystalNets.HALL_SYMBOLS
```

## Unstable nets
```@docs
CrystalNets.shrink_collisions
CrystalNets.order_collision
CrystalNets.expand_collisions
CrystalNets.CollisionNode(graph::PeriodicGraph, node::UnitRange{Int}, vmap=nothing)
CrystalNets.collision_nodes
```

## Utils
```@docs
CrystalNets.@toggleassert
CrystalNets.check_dimensionality
CrystalNets.check_valid_symmetry
```