# Source data attribution and licences

These are unchanged published benchmark inputs, not new caller runs. Files were
renamed locally for clarity; compressed bytes remain identical to their sources.
Exact paths, checksums, versions and limits are recorded in `SOURCE_MANIFEST.json`.

## Locityper

Timofey Prodanov, **Locityper loci database, 1KGP genotypes and benchmarking data**,
14 July 2025, pinned Zenodo record [15879519](https://doi.org/10.5281/zenodo.15879519).
The record metadata declares [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/).
Files under `locityper/` are reproduced from this record without content changes.
The record README reports benchmark v0.17.3; embedded table commands refer to
v0.17.0 summary filenames. That discrepancy remains unresolved. The current
website's older version statement must not override the pinned source evidence.

## COSIGT paper repository

Davide Bolognini and collaborators,
[COSIGT paper repository](https://github.com/davidebolo1993/cosigt_paper),
commit `3d84bf8490bbff5ec33cf72706e713d1a41730ed`.
Files under `cosigt/` are licensed under MIT. The complete original licence and
copyright notice are preserved as `cosigt/LICENSE`. Cite the accompanying paper:
Bolognini et al., *COSIGT: population-scalable genotyping of complex loci from
low-coverage sequencing data using pangenome graphs*, Genome Biology 27, 286
(2026), [doi:10.1186/s13059-026-04242-4](https://doi.org/10.1186/s13059-026-04242-4).

`locus_groups.tsv` is our derivative of the original CMRG BED: it adds overlap
component identifiers while preserving locus names and coordinates. The source
MIT licence and attribution apply to the underlying BED data.

External reads and assemblies are not redistributed in this snapshot. Their
source terms are separate from these benchmark-table licences.
