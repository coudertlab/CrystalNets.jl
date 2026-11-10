# Collision-translation regression inputs

`E01` through `E08` are all eight exact original internal-error graphs from a
Materials Project 2026.09.28 elemental-sublattice search. These preserve the
original ambient dimension 3; their connected nets have intrinsic dimension 1 or 2.
Each `.pgraph` is a native PeriodicGraphs serialization. The matching `.edges.txt`
file spells out one undirected periodic edge per line (one-based vertices;
reverse edge implicit). Per-case `.txt` files record the original observations.
No atom coordinates, API access, downloaded metadata, or Systre installation is needed.

`cases.tsv` maps each original graph to its intrinsic dimension and successful
post-fix CrystalNets genome. Those expected genomes are regression snapshots of
this branch, rather than independent Systre classifications: Systre did not
successfully classify these inputs. The tests additionally verify actual edge
preservation for accepted translation permutations and seeded representation invariance.

The three `small_*.pgraph` files are the reduced eight-vertex assertion example,
24-vertex bounds-error example, and six-vertex index-origin example from the PR
notes. The index-origin example isolates a wrong helper result, not an independent
third public-API crash. They are locally reduced/constructed examples, with no
claim of globally smallest size.

`test/collision_translation.jl` is included by `test/runtests.jl`, so the ordinary
`julia --project -e 'using Pkg; Pkg.test()'` runs these regressions in CI.
The file can also be included directly in an environment with CrystalNets and its
test dependencies. It asserts success, exact output, and mapping validity; expected
exceptions and zero-dimensional error sentinels are not treated as successful outputs.
