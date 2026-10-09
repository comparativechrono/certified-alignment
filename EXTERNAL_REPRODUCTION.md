# Reproduce the external comparison

Use Linux, Git, GCC/G++ with C++17 support, GNU Make, ar, and Python 3.9+.
No Python packages, CMake or Meson are needed. The reference FASTA is included;
network access is needed to fetch the pinned Edlib and WFA2 sources. See
[the dependency and timing details](benchmarks/EXTERNAL_DEPENDENCIES.md) and
[the fixed experimental protocol](EXTERNAL_BENCHMARK_PLAN.md).

Use the repository revision accompanying the published result. A newer revision
may change the experiment. The commands below preserve the published result
directory and freeze newly built binaries for a separate local run; machine and
compiler differences can change timings and binary hashes.

## Run a new measurement

From a fresh checkout:

```sh
git clone https://github.com/comparativechrono/certified-alignment.git
cd certified-alignment
```

Run the following in Bash from the repository root. Change the final directory
name for each additional run. The output path must be absolute and inside the
checkout because the runtime manifest stores paths relative to the repository.

```sh
set -eu
EXTERNAL_RUN="$PWD/results/external_reproduction"
test ! -e "$EXTERNAL_RUN"

make build/certalign build/reference_distance
python3 benchmarks/build_external.py
python3 tests/test_external_adapter.py \
  --output "$EXTERNAL_RUN/adapter_validation.json"

python3 benchmarks/run_external_bench.py \
  --output "$EXTERNAL_RUN" --generate-only

python3 tests/validate_repeat_oracle.py \
  --pairs "$EXTERNAL_RUN/pairs.tsv" \
  --metadata "$EXTERNAL_RUN/cases.json" \
  --output "$EXTERNAL_RUN/oracle_scores.json" \
  --expected-count 384 --expected-analytic 0

python3 benchmarks/run_external_bench.py \
  --output "$EXTERNAL_RUN" --freeze-runtime

python3 benchmarks/run_external_bench.py \
  --output "$EXTERNAL_RUN" --run
```

The existing oracle helper is reused with the explicit external-dataset counts;
its default counts belong to the earlier repeat experiment. It computes exact
scores independently before timing. The generator produces the same 384 pairs
from the committed reference and fixed seeds. The adapter test also checks the
wrapper against the unchanged prototype CLI.

After `--freeze-runtime`, do not rebuild or edit any frozen source, protocol,
input, oracle output or runtime file. Finish other builds and tests before
`--run`; competing CPU jobs contaminate timing. The harness uses one excluded
warmup and seven timed repetitions of all six methods, checks every returned
score/path and finite certificate interval, then writes its analysis automatically.
Failures remain failures and are not credited as speedups. Partial timing files
cannot be resumed or overwritten; use another fresh directory for a new run.

Inspect `validation.json`, `summary.json`, `per_case.json`, and `raw.jsonl` in
the reproduction directory. No published result is overwritten. Reproduction
directories matching `results/*reproduction*/` are ignored by Git.

## Reanalyze the published measurements without timing

The restoration helper verifies and decompresses the `.gz` files registered in
`results/COMPRESSED_RESULTS.json`. It refuses to overwrite different existing
bytes. Restore first, then analyze a separate copy so published summaries remain
untouched:

```sh
set -eu
python3 scripts/restore_recorded_results.py

EXTERNAL_ANALYSIS="$PWD/results/external_analysis_reproduction"
test ! -e "$EXTERNAL_ANALYSIS"
cp -a results/external_benchmark "$EXTERNAL_ANALYSIS"
python3 benchmarks/run_external_bench.py \
  --output "$EXTERNAL_ANALYSIS" --analyze
```

Reanalysis needs no external library build and performs no alignments. It uses
the recorded validated outputs and input hashes, so the checkout's engine source
must still match the published input freeze. It regenerates the copied summaries
and bootstrap intervals; it does not rerun score/path validation. For new native
measurements, use the complete build, validation and freeze sequence above.

`MANIFEST.sha256` identifies the original uncompressed v5 snapshot, not every
file in the evolving repository. Current experiment provenance is recorded in
its input/runtime freezes and build manifest.
