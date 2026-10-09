# Frozen external comparison: v6 evaluation

Fixed before candidate measurements on 9 October 2026. The v5 prototype is
unchanged; this tests its practical competitiveness under its existing global
unit-cost scoring model. The project and subsequent updates live at
https://github.com/comparativechrono/certified-alignment.

## Candidate, baselines and decision

The single candidate is the linear band predictor with sparse certification.
External comparators are Edlib and WFA2 in high and ultralow memory modes,
each configured for exact global Levenshtein alignment and a complete path.
Both WFA2 configurations are fixed in advance; no per-pair winner selection.
Prototype linear/length and linear/qgram are diagnostic internal baselines.

Edlib uses an unknown distance threshold (`k=-1`), global NW and PATH. WFA2
uses the edit metric, end-to-end span, compute-alignment scope, no heuristic
pruning, and one thread. No method receives the oracle score or generating
alignment. Matching is literal and substitutions, insertions and deletions each
cost one. Different optimal paths are accepted if they replay correctly.

The development gate requires candidate runtime at least 10% below EACH of
Edlib, WFA2-high and WFA2-ultralow on the primary genomic-background suite,
with each paired locus-bootstrap 95% upper ratio below one, and no correctness
or budget failures. Failure against any external configuration fails the gate.
A large loss supports stopping optimization of this prototype. Success would
justify further work, not establish biological or production benefits.

## Frozen workload

The primary 192 pairs are semisynthetic: actual E. coli K-12 MG1655 reference
background with controlled edits. This is not a real-read dataset and the noise
is not calibrated to a sequencing platform. Reference accession NC_000913.3,
assembly GCF_000005845.2; sequence length 4,641,652. Source and checksums are
recorded under `data/reference/`.

Select eight nonoverlapping 8,192-base reference windows reproducibly, seed
2026101001, before alignment. Each locus supplies centered nested source lengths
1,024, 4,096 and 8,192. For each length use:

- One interior insertion or deletion, positioned between one third and one half
  of the source; gap length `max(1,floor(n*fraction))`, fraction 0.02 or 0.08.
- Low noise: 0.4% substitution probability, no additional indels.
- Mixed noise: 2% substitution probability and 0.2% single-base indel-event
  probability per base of the target after the large edit, equally split between
  insertion and deletion events.

This gives 24 pairs per locus. A secondary set of 192 random-background pairs
shares the nested-window geometry, large-edit position/length, substitution
masks and nonzero replacement offsets, small-indel event masks/directions and
inserted bases. Random controls are analyzed separately. Metadata specifies the
exact edit-order convention and enough information to reconstruct every pair.

All generator choices, exact sequences, source/dependency hashes, seeds and
configuration are frozen before candidate output is inspected. Preserve all
cases, including failures and slow cases. No tuning on these outcomes.

## Fair native measurement

Use native C/C++ adapters for all methods, with preloaded input and full paths.
The prototype wrapper links the unchanged source and calls the same prediction
and alignment functions; it changes only how timings are collected. All are
compiled with `-O3 -DNDEBUG`, without `-march=native`. Pin official dependency
commits and record compiler and system information.

Allow documented WFA2 workspace reuse within a batch. Measure per-pair required
prediction/alignment, allocation, traceback and conversion to the shared CIGAR
format. Include workspace setup and teardown once per batch. Input parsing,
validation and output serialization are outside all method timings. Retained
output containers are handled consistently. Report any unavoidable differences.

Run each 192-pair suite as its own batch: seven timed repetitions and one
excluded warmup per method. Rotate/reverse method order and shuffle cases with
the same order for all methods within each repetition. No competing build,
oracle or correctness jobs during timing.

The main estimator sums per-case median alignment times and adds median batch
setup/teardown overhead for that suite. Also report each observed complete batch
total. For paired bootstrap, amortize that fixed overhead equally across the
192 cases and resample the eight complete loci, keeping all 24 variants together.
Use 9,999 resamples and seed 2026101002. Intervals describe synthetic-locus
variation conditional on this machine, not cross-machine uncertainty. Report
length/noise/direction/gap strata descriptively without replacing the endpoint.

Whole-process peak resident memory, if available, includes preloaded sequences
and retained outputs and is distinct from the prototype's bound-payload metric.
Do not compare payload bytes with another program's process RSS.

## Correctness, budgets and scope

Compute independent unrestricted two-row Levenshtein distances for all 384 pairs
before timing. Check the helper against the separate Python oracle. Replay every
returned path, check every score and finite prototype certificate interval, and
preserve explicit failure records. Compare adapter results with the frozen CLI
on a small independent validation set before benchmarking.

Keep prototype parameters unchanged: width 8, predictor block 64, seed k=11,
epsilon zero, order block 8, sparse work 4,000,000, bound payload 64 MiB, main DP
50,000,000 cells. External methods use their exact unrestricted algorithms;
timeouts and allocation failures are reported. Do not increase a cap to rescue
a case after seeing results. Time-to-error is never credited as a speedup.

The frozen v5 engine binary identity is
`4a492041691bd03ab0577152cfe0aa056a19f8f4be18551e43538d4a8d476a25`.
The native adapter binary is separately hashed because it changes the driver;
hashes of every linked prototype source must remain identical to v5.

No affine-gap, mapping, real-read accuracy, variant-calling or association-study
improvement follows from this experiment. A public reference background improves
sequence realism without removing the limitations of an injected-error model.
