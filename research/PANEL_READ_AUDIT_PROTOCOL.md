# Exploratory read-level audit of panel-sensitive discrepancies

10 October 2026. This protocol fixes case selection and evidence rules before
new structural/read-support outcomes are calculated. The historical baseline
at commit `9b33cd9af12c1750e94e299c34ab44cbb739265e` remains unchanged. This is
an **outcome-selected exploratory investigation**, not a new held-out accuracy
study or an evaluation of a novel algorithm.

## Question and scope

Which of the 22 candidate discrepancies contain a 50–1,000 bp pairwise indel,
and can public short reads support a diagnostic junction sequence that is
absent from the other donors' locus-panel sequences?

Three distinct claims must remain separate: a structural difference from the
closest supplied sequence, absence of equivalent structural alleles from the
whole panel, and observability using reads. This audit directly examines the
first and a bounded sequence-marker assay for the third. It does not establish
the second solely from a missing k-mer. A represented indel on another sequence
background may have different junction words. Repeat copy-number changes can
also be real yet yield no diagnostic 31-mer.

## Frozen cases and controls

Cases are the 22 held-out historical LOO calls which were accepted by the
calibration-selected fraction filter at matched 95% retention, have
`div > 0.005` and `avail_div > 0.005`, and have paired full-panel `div <= 0.005`.
They span 15 loci and ten donors. Their exact identities and source hashes are
in `data/panel_read_audit/cases.tsv` and `SELECTION_MANIFEST.json`.

Controls must be a different held-out donor at the same locus, with finite,
scorable full and LOO comparisons both at or below 0.005. They need not pass
the baseline confidence filter. Process cases by `(locus, sample)`, preferring
the least-used eligible control, then the smallest SHA256 of
`panel-read-audit-v1|control|<case_sample>|<locus>|<control_sample>`.
No native confidence metric, discrepancy magnitude or new read outcome ranks
controls. This produces 18 distinct controls for 21 matched cases; three
SMOC2 matches reuse a control. HG03492/MUC1 is unmatched. Do not enlarge the
control pool after seeing outcomes. The analysis unit for inventory is 40
unique donor–locus cases, not 43 independent case/control matches.

Attempt sequence assessment and read acquisition for every selected case and
distinct control. Keep unsupported, ambiguous and unavailable cases visible in
the denominator. No prospective sensitivity, specificity or general population
prevalence can be estimated from this truth-selected cohort.

## Sources and compatibility checks

Use the historical Locityper database `db.tar.gz` from pinned
[Zenodo record 15879519](https://doi.org/10.5281/zenodo.15879519), 785,474,256
bytes, official MD5 `f782f6e72426b7a818c7e1fd42f50277`. Verify the whole archive
and record its SHA256 plus every extracted member's hash before analysis.
Required members are the original `ref.bed`, unique haplotype FASTA, discarded
haplotype aliases, and all-v-all PAF alignments for the 15 selected loci.
Acquiring immutable inputs and examining schema/coordinates is allowed before
the analytical freeze; structural and read-support outcomes are not.

Resolve truth IDs through the discarded-haplotype alias groups. A sequence
remains available after excluding one donor if **any** alias belongs to another
donor. Never remove it solely because the retained FASTA identifier names the
excluded donor. Include the reference sequence if the historical evaluation
included it. Check the original evaluation semantics against official source
commit `146cfd42e9179bafd63d1f41db7fc449940989de` (v0.17.3), while preserving the
benchmark's unresolved v0.17.0/v0.17.3 provenance.

For both truth haplotypes, select the non-donor sequence with smallest upstream
`edit / alignment_size`; resolve exact ties lexically by representative ID and
report their count. Reproduce the archived closest-available diploid divergence
to its nine-decimal rounding precision (absolute tolerance 1e-8). Verify full
and predicted comparison compatibility where source records permit. If source
or alias incompatibility prevents reproduction, label the case unresolved
rather than silently changing its comparator. The closest-per-haplotype
aggregate is not asserted to be a globally optimal diploid ratio. If equally
minimal haplotype ratios have different alignment sizes, the lexical tie rule
can change the diploid aggregate; identify a tie-dependent incompatibility
explicitly rather than calling it source corruption.

Use public Illumina reads for the named donors, preferring the HPRC source of
the original benchmark. Record exact source URLs, indices, official object
checksums where supplied, alignment headers/read-group identifiers and library
provenance. A second sequencing library from the same donor is independent
support, not a rerun of the original benchmark; label it accordingly.

Extract the historical database's GRCh38 intervals, with 1,000 bp added to each
side and bounds clipped to the reference contig. Verify reference identities
against CRAM headers and reference MD5 information. The COSIGT refined BED is
not a substitute. Include every interval specified by the original locus BED.
No source may be sliced using coordinates from a different assembly.

Acquire indexed mapped-window reads with original sequence and qualities,
without MAPQ or proper-pair selection. Where fetching mates outside the window
is practical, do so under a documented deterministic rule; otherwise explicitly
record missing mates. Do not imply that mapped windows include all globally
unmapped reads, reads mapped elsewhere or paralogous copies. Persist source and
slice hashes, commands, windows, record counts and any acquisition failures.
Large external inputs may be reproducibly retrieved from their pinned sources;
do not commit whole-genome CRAMs or the complete database to this repository.

## Structural description

Inspect both truth-haplotype versus closest-available alignments in every
selected case/control. Parse PAF/CIGAR using its actual historical schema;
verify orientation, coordinate bounds and query/target consumption. Do not
silently treat a partial alignment as an end-to-end comparison or invent gap
events for unaligned sequence. Unparseable or incomplete comparisons remain
explicitly unresolved.

Record all CIGAR insertion/deletion runs, grouping lengths into <50 bp,
50–1,000 bp inclusive, and >1,000 bp. A qualifying single gap is **isolated**
only if no other positive-length I/D lies within 31 query bases of either
junction, including the boundary. Handle zero-query-length deletion junctions
and intervals consistently. Record compound differences separately. Absence of
an isolated qualifying CIGAR gap does not establish absence of a meaningful
repeat or structural change. These are alignment-representation descriptions,
not normalized structural allele classes or independent assembly validation.

## Frozen diagnostic marker assay

Use exact canonical 31-mers (minimum of the word and its reverse complement)
without ambiguous bases. For each qualifying gap, enumerate junction-spanning
words with at least ten bases on either side. A truth insertion has left and
right junctions; a truth deletion has one junction. Keep the junctions separate.

Keep raw marker counts for audit, but the primary junction assay requires
**SNP-clean aligned context**: bases outside the focal insertion must match
their aligned alternative bases (`=` or verified identical bases in `M`);
reject markers whose aligned context includes `X` or another indel. For a
deletion, all 31 query bases must satisfy this requirement around the gap.
The inserted bases themselves need not align to the alternative. Apply the
symmetric guard to alternative-junction markers. This reduces nearby-SNV
confounding without proving structural-allele absence.

A primary truth marker must be absent from every available non-donor sequence
at the same locus, respecting aliases, and from the panels of the other 14
audited loci. Record both filtering stages. This is a limited background
screen, not proof of genome-wide uniqueness. Alternative-junction markers must
be absent from both of the donor's truth haplotype sequences and pass the same
cross-locus background screen. Their presence does not itself quantify a
genotype likelihood.

For descriptive exposure, use the fixed truth-flank words at
`[q_start-62:q_start-31]` and `[q_end+31:q_end+62]`, if in bounds and present in
the closest alternative. Record unavailable words rather than move them to a
favourable position. These are local read-exposure checks, not calibrated depth
estimates or genome-wide origin certificates.

## Read evidence and interpretation

Count a marker only when all 31 read bases have Phred quality >=20. Primary
counts exclude secondary, supplementary, QC-failed and duplicate-flagged
records (flags 0x100, 0x800, 0x200, 0x400), with exclusions reported. Do not
require a proper pair. Retain MAPQ-zero primary records and separately report
the MAPQ>=20 subset. A locus-selection filter must not predetermine support.

Collapse fragment names across mates and overlapping windows. Record distinct
fragment names and fragment-start signatures per marker and junction. A
signature is the sorted set of mapped `(contig, alignment_start, strand)` and
mapped-mate `(next_contig, next_start, mate_strand)` positions, using available
SAM fields; both mates should produce the same signature. Mark incomplete
signatures when a mapped mate is unavailable. Multiple k-mers or both mates of
one fragment do not create independent molecules. Preserve witness identifiers
and relevant alignment/quality information for inspection. Distinct signatures
remain a heuristic against unflagged PCR duplicates, not molecular barcodes.

An isolated junction receives an **exploratory support flag** with at least
three distinct fragment names and at least two distinct fragment-start signatures.
A fully marked insertion requires this at both junctions; a deletion requires
it at its one junction. One-sided insertion evidence is partial. If either
required junction has no eligible marker, report it as unassayable by this
test, not unsupported biology. Report complex-gap evidence separately from the
primary isolated-gap flag. No threshold is calibrated as an error probability.

Zero support can arise from low exposure, sequencing errors, repeat ambiguity,
different sequencing libraries or mapped-window ascertainment. Three positive
fragments do not independently certify their locus of origin or the full
structural allele. Report all support counts, exposure and background-screen
limitations; do not turn unavailable or uninformative evidence into negatives.

## Outputs, audit and decision

Publish the protocol and case/control definitions before analyzing new
outcomes; freeze the executable analysis and acquired source manifests before
running it. Record freeze commit(s), code hashes, environment and commands.
Save resolved haplotypes, chosen comparisons, event descriptions, eligible
markers, per-junction/fragment support and a complete case-status table.
Independent review must check alias exclusion, CIGAR orientation/consumption,
marker specificity, fragment deduplication, every reported denominator and a
selection of positive witnesses directly from reads.

This is a bounded mechanistic follow-up. A useful result identifies concrete
cases with a meaningful structural difference and inspectable read evidence
surviving ordinary confidence filtering. It still does not show that a new
algorithm beats current Locityper, COSIGT or an appropriate local-assembly
method. Current-caller reruns require valid database/background/recruitment
inputs and a separately recorded experiment; binary version checks do not
count. No association-study benefit is evaluated here.

Recommend algorithm development only after structural equivalence and read
identifiability are adequately established against current baselines. Otherwise
report whether the limiting factor is missing truth, unsuitable event sizes,
repeat ambiguity, source mismatch, inadequate read evidence or existing tools.
Cases inspected here are development material for any later project; an
independent future test must use newly frozen cases.
