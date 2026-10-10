# Independent checks of the read-level case audit

10 October 2026. This document records checks performed independently of the
case-selection and structural/read-analysis implementations. It will distinguish
completed checks from checks that remain unavailable; it is not a certification
of assembly truth or a prospective accuracy study.

## Case selection: passed

The reviewer reconstructed the candidate set directly from the committed
`heldout_predictions.csv.gz` and `heldout_truth.csv.gz`, without importing the
selection program. The reconstruction required a leave-one-out call accepted
by the matched-95% unexplained-read-fraction filter, prediction divergence
greater than 0.005, available-panel divergence greater than 0.005, and paired
full-panel divergence no greater than 0.005.

All 22 candidate donor–locus keys matched `cases.tsv`: ten donors and 15 loci.
The same-locus controls were independently reconstructed from finite,
nonnegative full and leave-one-out divergence at most 0.005 and positive
comparison size. Candidate processing order, prior-control-use counts, and
the salted SHA-256 tie rule reproduced every match. There are 18 unique
controls for 21 candidate matches, three reused matches, and one unmatched
candidate, HG03492 at MUC1. Reused controls are not independent observations.

These are truth-selected cases from an already inspected held-out dataset.
They are now development material. Neither their ascertainment nor these
checks permit estimation of prospective detection accuracy or new-algorithm
improvement.

The controls have low overall sequence discrepancy; they are not proven
structurally represented negatives. A private 50–1,000 bp event can contribute
less than 0.5% across a sufficiently long diploid locus. Junction evidence in
a control therefore cannot be labeled an assay false positive.

## Interpretation and implementation checks required

The audit distinguishes three different assertions: a structural difference
between two haplotypes, absence of an equivalent structural allele from the
available panel, and distinguishability using real sequencing reads. A single
pairwise CIGAR, a panel-absent sequence word, or read support alone does not
establish all three.

Independent inspection of upstream Locityper `v0.17.3:extra/gt_dist.py` identified
two implementation requirements. First, discarded haplotype names are aliases
of identical sequences. Removing a donor's two haplotype names must retain a
representative sequence if any non-donor alias remains. Truth names may also
exist only as discarded aliases. Second, some historical PAF records contain
an extra `+` before the match-count fields; the upstream parser explicitly
handles this layout. CIGAR interpretation must additionally preserve alignment
direction and query/target coordinate consumption.

Breakpoint-spanning 31-mers require explicit junction coordinates, reverse-
complement canonicalization and a sequence-word specificity check. Markers
absent from the non-donor panel can reflect nearby substitutions even when a
structurally equivalent indel exists. Specificity against a collection of
locus panels is not whole-genome uniqueness.

Read counts must collapse mates and repeated regional records to fragment
names, preserve the role of mapping ambiguity, and separately identify
duplicate, secondary, supplementary and QC-failed records. Support thresholds
are exploratory evidence rules, not calibrated genotype likelihoods. No
support can reflect limited coverage, sequencing error, repeats, or read
recruitment rather than absence of the allele. Regional CRAM extraction also
misses molecules mapped elsewhere or wholly unmapped.

## Structural and read evidence

Pending independent review of completed outputs. The case-selection check
above does not certify these downstream stages.

## Checks performed before outcome execution

The rolling canonical 31-mer implementation was compared with an independent
naive string-window calculation on 200 randomized sequences, including
ambiguous bases and the Phred-19/20 boundary. All comparisons passed. The
review also checked PAF direction reversal and gap-operation swapping against
upstream conventions and the source implementation's synthetic tests.

Several issues were identified and corrected before analytical execution:
the comparator tie rule now follows the frozen lexical representative rule;
archived comparisons use the stated divergence tolerance; incomplete
cross-locus backgrounds fail closed; fixed flank words have explicit
availability states; and paired read starts collapse to one fragment
signature. Read consumption verifies each regional BAM against its recorded
checksum. These are implementation corrections before outcomes, not choices
made to improve an observed result.

During acquisition review, truncated remote reads were followed by nominally
successful empty BAM files from a reused remote handle. These files cannot
represent biological absence or valid negative evidence. Acquisition must be
repaired or those cases retained as unresolved before read-support outcomes
are interpreted. All first-attempt BAMs were excluded, and acquisition was
repeated with verified local reference chromosomes and fresh remote handles.

The replacement acquisition independently passed exact 40-case/control key
coverage, every BAM's SHA256, sequential BAM decoding and recorded counts,
donor read-group identity, reference chromosome MD5 metadata linkage, and
the full local reference FASTA SHA256. All 40 slices were nonempty and marked
verified: 993,855 alignment records and 70,545,821 BAM bytes. The local FASTA
checksum matches the source manifest that records actual sequence-MD5
verification for all ten selected chromosomes. The independent check verified
that linkage; it did not separately repeat each chromosome MD5 calculation.
See `results/panel_read_audit/acquisition_independent_audit.json`. No biological
marker-support outcome was calculated during this acquisition check.

Even a marker occurring at both insertion junctions may be the same repeated
sequence word, potentially also present inside the truth insertion. Support
for both junction word sets therefore need not independently localize two
breakpoints. Positive results require direct witness inspection and retain
the protocol's narrower sequence-evidence interpretation.
