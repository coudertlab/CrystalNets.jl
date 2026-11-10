# Paired-translation fingerprint regression inputs

This folder contains every original disagreement pair `D01` through `D23`:
23 raw source `.pgraph` files and 23 companion graphs decoded from recorded Systre
keys. They preserve ambient dimension 3 and every original periodic edge. Each
has a matching human-readable `.edges.txt` file and a per-case `.txt` observation
record. Native serialized graphs are loaded with `PeriodicGraph(read(path, String))`.
No structures, atom coordinates, API access, or Systre installation is needed.

`cases.tsv` identifies the source/companion files, intrinsic dimension, and expected
normalized CrystalNets genome. Equality of the underlying infinite graphs is the
requirement; literal agreement with Systre's numbering convention is not required.
`primitive_pairs.tsv` contains the 23 topology-preserving primitive-cell comparisons
from the PR description, including both recorded fingerprints before normalization.

`R01` through `R12` are the twelve distinct smaller comparison pairs in the PR
notes; `reductions.tsv` maps all 23 cases to them. The reduced examples have 2–6
vertices and need not retain the topology of their original material graph.
`synthetic_3d.pgraph` adds a constructed 3D paired net to exercise all eight
translation-generator choices. The reduced examples are locally minimal under the
documented pair-compatible reduction moves,
not claimed globally smallest.

`test/paired_translations.jl` is included by `test/runtests.jl`, so the normal
`julia --project -e 'using Pkg; Pkg.test()'` runs the regressions in CI. It checks
all original pairs and primitive pairs, the reduced pairs, relabeling, representative
offsets, axis permutations, doubled cells, idempotence, preservation of edges/degrees,
and exclusion of independent local flips. Distinct chain/ladder and square/triangular
controls guard against trivial collapse to a shared fingerprint. Successful, nonzero
genomes are required before equality is checked. No new Julia dependencies are added.

## Independent equivalence certificates

The `proof/` directory includes integer-edge isomorphism witnesses for all 23
original pairs, the D13-versus-D14 equivalence, and all 23 reductions. The optional
Python verifier uses only the standard library and does not invoke CrystalNets,
Systre, floating-point placement, or any graph-matching library. Run:

```sh
python3 test/fixtures/paired_translations/proof/verify_input_certificates.py
```

The original-pair witnesses cover intrinsic lattice extraction, cell covers,
vertex offsets, pair exchanges and unimodular basis changes. These supporting
proofs do not add Python as a dependency of the Julia test suite.
