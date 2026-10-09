# Native external alignment baselines

The comparison uses the official native APIs, compiled from pinned public source.
Dependency sources and binaries are build inputs, not vendored project files.

| Library | Official source and pinned commit | Version description | License |
|---|---|---|---|
| Edlib | [Martinsos/edlib, `0ddc23ea06abd95ab6227442acb4c27bd7607b68`](https://github.com/Martinsos/edlib/tree/0ddc23ea06abd95ab6227442acb4c27bd7607b68) | `1.3.9.post1-8-g0ddc23e` | [MIT](https://github.com/Martinsos/edlib/blob/0ddc23ea06abd95ab6227442acb4c27bd7607b68/LICENSE) |
| WFA2 | [smarco/WFA2-lib, `3282f834713cb9767884949143fd0c803acae697`](https://github.com/smarco/WFA2-lib/tree/3282f834713cb9767884949143fd0c803acae697) | `v2.3.6-21-g3282f83` | [MIT](https://github.com/smarco/WFA2-lib/blob/3282f834713cb9767884949143fd0c803acae697/LICENSE) |

The version descriptions are `git describe` output, not claims that these commits
are tagged releases. `build/external_build.json` records compiler versions, build
commands, license hashes, dependency commits, frozen prototype source hashes,
adapter source hash and adapter binary hash.

## Build and validate

```sh
make
python3 benchmarks/build_external.py
python3 tests/test_external_adapter.py
```

The build script needs Git, GCC/G++, GNU Make, ar, and Python 3. It can use an
existing dependency directory with `--deps /absolute/path`. Existing tracked
source modifications or a different checked-out commit are rejected. It builds
Edlib's single source translation unit directly, avoiding a CMake/Meson dependency.
WFA2 uses its own Makefile. WFA2 is forcibly rebuilt to avoid silently reusing
objects made with different compiler flags.

All algorithm code is compiled with `-O3 -DNDEBUG`, without `-march=native`, LTO,
or OpenMP. C++ uses C++17. WFA2's make build receives `-fPIE` as well. Execution
is single-threaded. Optimized production libraries can use hardware features
available under their normal generic compiler target; the scalar prototype has
not been rewritten or vectorized for this comparison.

## Alignment settings

All methods compute exact end-to-end unit-cost Levenshtein alignment and return
a complete path. No oracle distance, error-rate estimate, fixed score threshold,
approximation tolerance, free end gap, or heuristic pruning is supplied.

* Edlib uses `edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0)`.
  Its internal threshold selection is enabled by `k=-1`.
* WFA2 uses `distance_metric=edit`, `alignment_scope=compute_alignment`,
  `alignment_form.span=alignment_end2end`, `extension=false`, and
  `heuristic.strategy=wf_heuristic_none`. Both `wavefront_memory_high` and
  `wavefront_memory_ultralow` are measured separately. The latter uses WFA2's
  bidirectional algorithm. Other system settings are the upstream defaults;
  `max_num_threads=1` is explicit. Any non-complete status is a hard error.
* Prototype methods use `cba::predict` and `cba::align` from the frozen source,
  a linear predictor and the length, q-gram, or sparse bound. All options retain
  their v5 defaults: initial width 8, exact tolerance zero, cumulative DP limit
  50,000,000 cells, bound payload cap 64 MiB, ordered block size 8 and sparse
  proof-work cap 4,000,000. A prototype work-limit failure remains an explicit
  uncertified output; the harness must not count it as a fast successful alignment.

The build only renames the existing prototype CLI entry point using the compiler
definition `-Dmain=certalign_cli_unused_main` on its translation unit. Engine and
predictor files are neither edited nor copied into a new implementation.

## Batch interface and measurement

```sh
build/external_adapter --method edlib --input examples/pairs.tsv
build/external_adapter --method wfa2-ultralow --input examples/pairs.tsv
build/external_adapter --method prototype-sparse --input examples/pairs.tsv
```

Methods are `edlib`, `wfa2-high`, `wfa2-ultralow`, `prototype-length`,
`prototype-qgram`, and `prototype-sparse`. Optional `--repeats N` and `--warmups N`
default to 1 and 0. Every TSV row is `id<TAB>x<TAB>y`; empty sequences are allowed.

The adapter preloads the entire TSV batch, then executes all alignments without
printing or reading input between them. It creates one WFA2 workspace and reuses
it for every pair and repeat, following the upstream recommendation in the
[WFA2 README](https://github.com/smarco/WFA2-lib/blob/3282f834713cb9767884949143fd0c803acae697/README.md).
Edlib uses its ordinary per-call API. Prototype prediction and bounds are built
per pair, just as in the frozen implementation.

Per-pair `total_ms` includes native alignment, traceback, conversion to a compact
`=XID` CIGAR string, and release of per-call work structures. It excludes input
parsing, JSON serialization, the retained returned CIGAR/diagnostic objects, and
their eventual destruction. The prototype's diagnostic trace is retained for
independent certificate auditing. WFA2 workspace allocations performed during
alignment are naturally included; its one-time construction and deletion are
reported separately as `batch_setup_ms` and `batch_teardown_ms`.

For a batch with one measured repeat and no internal warmup, the primary compute
time is `batch_setup_ms + alignment_sum_ms + batch_teardown_ms`. `wallcompute_ms`
also reports actual elapsed time across that whole region, including loop/result
storage overhead. With internal warmups, it includes warmup work as well;
`warmup_alignment_ms` makes that explicit. Process startup and file I/O are
excluded from these compute measures for every method. Any external process RSS
measurement includes preloaded inputs, retained outputs, linked libraries and
runtime overhead, and is not an algorithm-only memory figure.

The first stdout JSON line is a `type=batch` metadata object. Remaining lines are
`type=alignment` records with method, ID, repeat, lengths, distance, lower bound,
CIGAR, exactness, timing, prototype diagnostics and all certificate round bounds.
Unavailable prototype bounds use `-1`; libraries return `lower_bound=distance`
and an empty trace. A failure exits nonzero, except the frozen prototype's
explicit per-record work-limit/error result, which preserves its error field.

The CIGAR convention is deletion `D` consumes x, insertion `I` consumes y.
Edlib's raw `EDLIB_EDOP_INSERT` consumes its query x and is therefore mapped to
`D`; its `EDLIB_EDOP_DELETE` maps to `I`. WFA2's native I/D convention already
matches. Its `M` is an exact match and is normalized to `=`. Independent replay
checks every consumption, match/mismatch, and total path score.

Two empty-input API conventions are normalized explicitly. Edlib's empty-input
fast path returns the correct distance but omits operations; the adapter creates
the unique all-insertion or all-deletion path. WFA2 ultralow can leave its raw
`cigar->score` unset for empty inputs even when the returned CIGAR is valid. The
adapter uses WFA2's public `cigar_score_edit()` for every WFA2 result, matching
the edit-distance scoring in upstream benchmark utilities. Normalization and
scoring are included in measured time; no expected/oracle score is supplied.

## API validation

`tests/test_external_adapter.py` compares every result with an independent
unrestricted scalar edit-distance oracle. It covers all 961 pairs of binary
strings of length at most four, empty sequences, mixed-case/ambiguous symbols,
byte values including NUL, and lengths around 64/128/256. Larger divergent and
noisy inputs exercise nontrivial ultralow-memory work. Cases run in mixed order
with one warmup and three retained repeats to test reused workspace state.

The same test checks every path, every finite prototype certificate bound,
repeat determinism, and equivalence of wrapper outputs to the unchanged prototype
CLI (distance, path, rounds, work, fallbacks and traces). The production benchmark
separately validates each of its outputs; passing API tests alone does not
establish accuracy or speed on that workload.
