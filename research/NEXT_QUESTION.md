# Proposed next question: detecting inadequate reference haplotype panels

Literature assessment: 9 October 2026. Status: research recommendation, not an
implemented method, established novelty claim, or completed experiment.
The v6 alignment prototype remains a documented negative result.

## Recommendation

**Can DNA read alignments detect when a sample carries an indel or structural
haplotype that is inadequately represented in a reference panel, and reduce
confidently wrong genotype calls at the same call retention?**

Start with targeted diploid loci containing insertions/deletions and repeats.
The intended output distinguishes adequate panel support, evidence of panel
inadequacy, and insufficient information. It does not reconstruct every missing
allele. This is an alignment-evidence problem at the interface with genotyping,
with a potential application to genotype–phenotype studies.

The key distinction is between ranking the supplied candidates and testing
whether any supplied candidate adequately explains the observations. If the
panel offers A and B but the sample carries C, a large A-versus-B score margin
does not establish that A is correct. The actual difficulty is distinguishing
a missing allele from sequencing errors, low coverage, contamination, incorrect
locus mapping and read-recruitment failure.

## Evidence for the gap and closest prior art

| Source | What it already does / remaining limitation |
| --- | --- |
| [COSIGT, September 2026](https://doi.org/10.1186/s13059-026-04242-4) | Selects complex-locus haplotypes; explicitly describes assignment of novel structural variants to the closest available haplotype. It already evaluates excluded samples/panels and uses a cosine-score margin for filtering. A new holdout benchmark or another margin threshold alone is not novel. |
| [Locityper official documentation](https://locityper.vercel.app/genotype) | Provides genotype quality, weighted distance, unexplained-read counts and warnings. The authors give no general filtering rule. Its unexplained-read statistic concerns reads fitting some panel haplotypes but not the selected genotype; it is not itself a test that all panel haplotypes fit poorly. |
| [Orthanq, July 2024](https://doi.org/10.1186/s12859-024-05832-4) | Already supplies probabilistic haplotype uncertainty and examination of mutations beyond reported HLA types. Generic uncertainty-aware genotyping is established prior art. |
| [Kourami, 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC5804087/) | Already uses graph-guided assembly to discover novel HLA alleles. Missing-allele discovery in general is not new; include an appropriate assembly alternative where it is applicable. |

The proposed residual question is whether an explicit adequacy model improves
error-versus-retention performance across unseen samples and loci, beyond these
existing uncertainty measures and simple filters. The search supports an
unresolved practical limitation. It does not prove that nobody has studied the
exact proposed method. No new method has yet been specified or tested.

## Why prioritize this question

It links alignment algorithms to the quality of a genotype matrix, and it can
be tested against mature tools without building another whole-genome aligner.
The potential value is avoiding unsupported genotype assignments and identifying
samples needing additional analysis. Improved association power or false-positive
control would require a separate downstream study. Selective missingness can
itself bias association analyses, so simply issuing more no-calls is not success.

Several other directions remain scientifically important but have closer prior
art or less clearly isolated opportunities:

| Direction considered | Reason not selected for this first project |
| --- | --- |
| Repeat sequence/interruption genotyping | A [2026 ONT benchmark](https://doi.org/10.1186/s13059-026-04210-y) identifies residual sequence-level errors, especially with long, heterozygous and nonreference alleles. But [LongTR](https://doi.org/10.1186/s13059-024-03319-2) already uses candidate haplotypes and probabilistic realignment, and [nPoRe](https://doi.org/10.1186/s12859-023-05193-4) already models repeat-dependent alignment penalties. A new alignment model needs a better isolated failure. |
| General mapping-confidence calibration | [Qtip](https://doi.org/10.1186/s13059-017-1290-3) already addresses mapper-specific calibration. MAPQ concerns locus placement, not whether local indel evidence is correctly represented. |
| General reference-bias correction | [biastools](https://doi.org/10.1186/s13059-024-03240-8) documents remaining indel-related assignment bias, but local realignment, allele-aware mapping and probabilistic calling are established approaches. A generic proposal would be insufficiently distinct. |

## First milestone: test the existing methods

This is a proposed study outline. Freeze the exact loci, versions, event
definitions, exclusions and decision rule before examining new test results.

1. **Reuse public benchmark material.** Start from the
   [COSIGT paper repository](https://github.com/davidebolo1993/cosigt_paper),
   with its published predictions and scripts, and the linked HPRC/HGSVC
   assemblies and matching 1000 Genomes reads. Use public benchmark data;
   the controlled-access Moli-sani cohort is unnecessary. Begin with a small
   fixed locus subset and paired-end short reads shared by COSIGT and Locityper.
   Assess long-read generalization only after the first result.
2. **Define meaningful panel inadequacy.** Exact nucleotide absence is too
   strict: a private SNV need not invalidate the structural genotype. A sensible
   pilot focuses on independently assembly-supported, scorable 50–1,000 bp
   insertion/deletion allele classes at bounded loci. Define allele equivalence
   and tolerances before testing, accounting for shifted repeat breakpoints.
   Retain uncertain truth as unscorable, rather than making either caller truth.
3. **Use both natural and controlled exclusions.** Remove whole donors and
   relatives from the panel. Separately remove all representatives of specified
   structural allele families. Include represented-allele controls. Development,
   calibration and final testing must separate both donors and loci; duplicate
   allele copies must not leak across the intended exclusions.
4. **Establish observability.** Check for informative breakpoint-spanning or
   other discriminating fragments. Assembly truth gives the genotype, not the
   origin of every real read in a repeat. Use simulations only for controlled
   origin/identifiability checks. Measure natural donor exclusions separately
   from deliberately damaged panels. Breakpoint evidence may reveal an
   unsupported indel without determining its full sequence or distinguishing
   all allele families; evaluate those observability requirements separately.
5. **Make the baselines difficult to beat.** Compare current pinned Locityper
   and COSIGT, their documented filters, absolute fit as well as score margins,
   and simple error/clipping/breakpoint-consistency filters. Give ordinary
   post-hoc calibration of existing outputs the same training data. Include
   an applicable local-assembly caller on the selected loci; Orthanq/Kourami
   are additional comparators if HLA enters the study. Count recruitment,
   preprocessing and filtering costs.
6. **Measure error at matched retention.** Plot structurally incorrect accepted
   calls against the fraction of all eligible sample–locus cases retained.
   Report missing-allele detection separately from generic failure detection.
   Check coverage, repeat class, allele frequency and panel-representation
   strata; report differential missingness. Summarize uncertainty across both
   donors and loci, not as independent read-level observations.

An example later method-evaluation gate is at least **20% fewer structurally wrong accepted
calls at 95% call retention**, versus the strongest calibrated baseline, with
an uncertainty interval excluding no improvement and acceptable measured cost.
These percentages are proposed project criteria, not literature-established
clinical thresholds. Finalize them after feasibility checks, before held-out
outcomes. If too few errors exist to assess the gate reliably, do not claim a win.
The first milestone only establishes whether sufficient observable errors
survive existing calibrated filters to justify developing a method.

**Stop before algorithm development** if existing filters already remove nearly
all avoidable errors, if reference/truth uncertainty explains the apparent gap,
or if the available reads do not distinguish missing from represented alleles.
Require an effect on natural donor exclusions; a gain only on artificial panel
deletions is insufficient.

## Possible algorithm, only after demonstrating headroom

An inexpensive local model could compare the fit of all panel haplotype pairs
with an error-aware alternative that allows consistent unsupported indel
evidence. Repeated agreement across independent molecules, flanking placement,
and residual alignment structure could help separate an unrepresented allele
from isolated sequencing noise. Reads poorly aligned to every panel haplotype
must remain available; examining only the already accepted reads can hide the
failure being tested.

The model must not win merely by adding unconstrained alternatives that overfit
sequencing errors. Compare it against simple absolute-fit statistics and
ordinary calibration before adding complexity. Report empirically calibrated
uncertainty under the evaluated conditions, not a universal mathematical
guarantee under arbitrary population or technology shifts. Novel-allele
reconstruction and partial genotypes are later questions, outside the pilot.

The first deliverable should therefore be a baseline-and-headroom report.
Implementation is justified only if that report identifies a material,
observable failure surviving the strongest existing approaches.
