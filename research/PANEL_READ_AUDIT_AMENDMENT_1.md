# Read-audit source and acquisition amendment

10 October 2026, before any new cohort alignment/event classification or
junction-read-support analysis. This amends the source provisions of the
[frozen protocol](PANEL_READ_AUDIT_PROTOCOL.md), first published at
`986336d1308ef903123fcbc52f0084e1659c972d`. Case identities, controls, marker
rules, thresholds, exclusions and interpretation gates remain as specified.

## Missing published alignment files

The complete current database archive passed its official MD5 check. Its
SHA256 is `0f3082cd6a032f52fd4758c7f1bab2a795ececf97de7d3e1b695e5a2fe79a992`.
It contains FASTAs, BEDs and alias files but **no PAF alignment files**, despite
their description in the upstream documentation. All 40 selected case/control
records' truth and predicted haplotype identifiers resolve in the verified
sequence/alias data. Numerical alignment compatibility remains untested until
an alignment source is available.

We fix this source-recovery decision tree before observing structural or read
support outcomes:

1. Inspect the older official database from Zenodo record
   [12607584](https://zenodo.org/records/12607584), `db_327.tar.gz`,
   934,483,715 bytes, official MD5 `7d5982b60529c8457fd1cc72fcc8ffac`.
   Reuse published PAFs only if the selected loci have identical normalized
   FASTA identifier-to-sequence mappings, equivalent alias partitions and
   identical reference intervals to the primary verified archive. Merely
   having the same locus names is insufficient. Record all source/member hashes
   and the original source of each recovered alignment.
2. If compatible archived PAFs are unavailable, regenerate alignments with the
   **release-matched Locityper 0.17.3** executable described below, using the
   verified primary FASTAs. This branch is selected by source availability,
   never by which alignment parameters give more favourable results.
3. Apply the already specified archive-divergence compatibility gate to every
   case/control. Keep failures unresolved and report their numbers. Do not
   search alternative parameters, binaries or alignment paths to improve
   agreement after seeing these results.

Recovered or regenerated CIGARs describe one alignment representation. Even
matching scalar discrepancies cannot prove that a regenerated path is the
identical path used in the original benchmark, or establish normalized indel
equivalence across repeats.

## Fixed release-matched regeneration

The Bioconda package `locityper-0.17.3-h919a2d8_0.tar.bz2` has SHA256
`8a9d28e7ba8a9605fbab6a18d370dc41948822365eadbd64413a7aadfff3db0a`.
The downloaded executable reports v0.17.3; its hash, package metadata, recipe
source and synthetic-only smoke checks are recorded under
`studies/panel_read_audit/legacy/`. It runs locally without an installation into
system directories. The package recipe did not pin its WFA2 dependency to a
commit. Therefore it is **release-matched regeneration**, not proof of using
the authors' identical benchmark executable.

The official `locityper-bench` README at commit
`ed5c5b4c6f7495287a9fff6b62eb7666409ef2e2` (SHA256
`363f616250840f94fe062c7f01b57763c3d20be7869b536616aea540992b1f02`)
specifies all-pairs alignment with `-A -D 1`. Make the other historical defaults
explicit, without substituting current-version defaults:

```sh
locityper align -i haplotypes.fa.gz -o haplotypes.paf.gz \
    -A -D 1 -m 15 15 -k 25,51,101 -g 500 -a 9 \
    -M 4 -O 6 -E 1 -@ 8
```

Run this fixed command for all 15 selected loci if regeneration is needed.
Keep source FASTAs unchanged, record input/binary/output hashes and commands,
and preserve stderr/exit status. Apply the same rule to cases and controls.
The runner limits each locus to 1,800 seconds and 16 GiB address space; timeouts
and resource failures remain unresolved and are reported, never silently dropped.
The alignment method uses heuristics and is not asserted to be globally exact.
The structural engine's lexical nearest-representative tie rule and 1e-8
divergence compatibility tolerance remain unchanged; independently reproduce
the upstream evaluator's original tie semantics for comparison.

## Failed remote-reference acquisition

The first read-acquisition attempt suffered remote reference decoding failures.
Some subsequent requests on failed handles returned empty outputs. Every BAM
from that attempt, including nonempty ones, was set aside; none is eligible
for the evidence analysis. Preserve its logs and an explicitly unusable
manifest as failure provenance. Empty files are not evidence of absent reads
or alleles.

Reacquire all 40 cases/controls using complete required GRCh38 reference
chromosomes downloaded through indexed byte ranges. Verify the actual sequence
of each chromosome against the official dictionary's MD5, record source ranges
and the local FASTA SHA256, and retain normal CRAM reference checks. This
changes transport and verification, not coordinates or biological selection.

Use a fresh CRAM handle for each locus. A usable slice requires positive record
count, successful complete local BAM decoding, matching acquisition/scan record
counts, a verified BAM SHA256 and explicit verified status in the manifest.
Per-locus failures stay visible without discarding other completed records.
The scanner must refuse unverified or stale outputs. Whole source CRAM MD5s
are publisher-side checksums, not locally recomputed whole-file hashes: only
regional reads are retrieved. Exact historical library equivalence remains
unproven unless supported by additional recorded metadata.

## Interpretation clarifications

The controls have low overall **sequence discrepancy** under the historical
metric; they are not proven structurally represented negative controls. A
private 50–1,000 bp event in a long locus can pass the 0.5% cutoff. Finding
junction support in a control therefore cannot be called a false positive.

A junction word can repeat at multiple positions in a truth sequence or occur
at both insertion boundaries. Supporting both word sets does not then localize
two independent breakpoints. Inspect multiplicity in positive cases and
describe the result as sequence-word evidence unless localization is
independently established. Keep the pre-specified support rule unchanged.

Record indel-base burdens in each size bin and mismatch counts, in addition to
event counts. A supported target-sized gap might not explain most of a
candidate's original discrepancy if a much larger event dominates it. These
are descriptive safeguards, not outcome-dependent eligibility changes.

## Source gate outcome before analytical execution

The older archive passed its official MD5 check; its SHA256 is
`72b839ac883a1fdf92a0349a87cac32e819a5a19f2e85ce8c389bd6e7fde2bf7`.
It contains the original PAFs. All 15 selected loci passed exact normalized
sequence, alias and reference-interval equivalence checks. The unmodified
upstream v0.17.3 evaluator reproduced all 120 full-panel, leave-one-out and
available-panel comparisons across the 40 selected units: all 480 edit,
alignment-size, divergence and quality fields matched. The audit therefore
uses the recovered published PAFs; the regeneration fallback is unused.

All 40 replacement read slices passed verification and are nonempty, with
993,855 alignment records in total. Their input manifest SHA256 is
`ed282f12f9d417ac12495169035233b31a7b095404b0cf1b35b857cbad62a86a`.
These availability and provenance checks preceded structural classification
and junction-support analysis.
