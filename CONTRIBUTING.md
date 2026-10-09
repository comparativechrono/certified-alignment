# Project workflow

The canonical project is https://github.com/comparativechrono/certified-alignment.
Keep future source changes, experimental protocols and reviewed results here.
The initial import preserves the v5 source and data. Large recorded result files
are stored losslessly as gzip to keep the repository manageable. Run
`python3 scripts/restore_recorded_results.py` to restore and verify their original
bytes. `results/COMPRESSED_RESULTS.json` records both checksums;
`MANIFEST.sha256` identifies the complete uncompressed v5 snapshot.

Build with `make`. `make test` checks global unit-cost alignment against an
independent oracle; `make test-sparse` audits the sparse lower bound and complete
certificates. See README.md for the earlier regression targets and benchmarks.
Build products and local reproduction directories are excluded from git.

For an experiment, state the candidate, baselines, scoring, workload, budgets
and decision rule before timing. Pin input and dependency identities. Include
preprocessing, certification and traceback costs consistently, and preserve
failures and unfavorable cases. Keep independent score and script checks.

Commit new results alongside their reproducible harness and protocol. Preserve
published historical runs: use a new result directory for a rerun. Distinguish
primary endpoints from exploratory subsets and synthetic data from actual
sequencing reads. A speedup between variants of this prototype is not an
external production-aligner comparison.

Keep changes to the alignment engine separate from evaluations of a frozen
engine. Do not tune a tested candidate on the data used to claim its performance.
