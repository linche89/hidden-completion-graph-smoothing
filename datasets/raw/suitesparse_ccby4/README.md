# SuiteSparse Matrix Collection subset

Retrieved on 2026-09-13 from the official SuiteSparse Matrix Collection.
The original Matrix Market archives are retained without modification.

The Collection states that its matrix data are licensed under CC BY 4.0.
`LICENSE_CC-BY-4.0.txt` is the Creative Commons legal code downloaded from
https://creativecommons.org/licenses/by/4.0/legalcode.txt. Matrix-specific
provenance and citation metadata remain embedded in the original archives.

Selected matrices:

- `usroads.tar.gz`: sparse undirected road graph, for anchored graph-Laplacian
  and pairwise relative-measurement WLS experiments.
- `wathen100.tar.gz`: real symmetric positive-definite matrix, for generic
  sparse information-system tests.
- `bcsstk18.tar.gz`: real symmetric positive-definite structural matrix with
  severe conditioning, for numerical robustness tests.

These last two matrices are solver/information-matrix tests, not claimed to be
literal sensor-network measurement models. A sparse SPD matrix does not imply
that a sparse local measurement factorization is supplied with the dataset.

