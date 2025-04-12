# Internal types and functions

See also the documentations of [`PeriodicGraphs.jl`](https://liozou.github.io/PeriodicGraphs.jl/)
and of [`PeriodicGraphEmbeddings.jl`](https://liozou.github.io/PeriodicGraphEmbeddings.jl)

## Types

```@docs
CrystalNets.Crystal
CrystalNets.Clusters
CrystalNets.CollisionNode
CrystalNets.CIF
```

## Core topology functions

```@docs
CrystalNets.topological_key
CrystalNets.CRYSTALNETS_ARCHIVE
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
CrystalNets.remove_partial_occupancy
CrystalNets.prune_collisions
CrystalNets.expand_symmetry
CrystalNets.trim_monovalent
CrystalNets.trimmed_crystal
CrystalNets.trim_topology
```

## Bond guessing

```@docs
CrystalNets.guess_bonds
CrystalNets.edges_from_bonds
```

## Clustering algorithm

```@docs
CrystalNets.find_sbus!
CrystalNets.regroup_sbus
CrystalNets.regroup_paddlewheel!
CrystalNets.split_sbu!
CrystalNets.reclassify!
CrystalNets.add_to_newclass!
CrystalNets.group_cycle
CrystalNets.collapse_clusters
CrystalNets.pem_to_pe
CrystalNets.allnodes_to_singlenodes
```

## Unstable nets

```@docs
CrystalNets.collision_nodes
CrystalNets.CollisionNode
```

## Archives

```@docs
CrystalNets.make_archive
```

## Utils

```@docs
CrystalNets.@toggleassert
CrystalNets.check_dimensionality
```

## Other

```@docs
CrystalNets.guess_topology
CrystalNets.guess_topology_dataset
CrystalNets.recognize_topology
CrystalNets.total_interpenetration
```
