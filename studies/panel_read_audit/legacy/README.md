# Release-matched alignment regeneration

The verified Zenodo v0.17 database advertises haplotype PAFs but does not
contain them. These scripts prepare a bounded fallback; they do not establish
that regenerated alignments equal the original benchmark alignments.

The Bioconda package `locityper-0.17.3-h919a2d8_0.tar.bz2` is 2,904,704 bytes,
SHA-256 `8a9d28e7ba8a9605fbab6a18d370dc41948822365eadbd64413a7aadfff3db0a`.
Its executable runs as `locityper v0.17.3`, SHA-256
`cc73d9e285584e1969f68637ce95f789e77c985fbd76c9c818cce632fc7b371d`.
The source release resolves to commit
`146cfd42e9179bafd63d1f41db7fc449940989de`.

Preparation requires no Conda installation, Docker daemon, or system changes:

```bash
python studies/panel_read_audit/legacy/prepare_legacy_aligner.py \
  --destination /absolute/scratch/path/locityper_legacy_0.17.3 \
  --manifest-output /absolute/scratch/path/runtime_verification.json
```

Preparation performs only synthetic sequence alignment. It verifies all six
pairs of four synthetic sequences, including a 60 bp insertion, an 80 bp
deletion and a one-base substitution. Checks cover actual aligned bases,
query/reference consumption, pair count, edit count, score and divergence.
The committed `runtime_verification.json` records that completed smoke check.
The all-pairs runner was also exercised on these synthetic inputs.

The published benchmark command is:

```bash
locityper align -i loci/$locus/haplotypes.fa.gz \
  -o loci/$locus/haplotypes.paf.gz -A -D 1
```

Source: [`locityper-bench/README.md`](https://github.com/tprodanov/locityper-bench/blob/ed5c5b4c6f7495287a9fff6b62eb7666409ef2e2/README.md),
SHA-256 `363f616250840f94fe062c7f01b57763c3d20be7869b536616aea540992b1f02`.
Our command explicitly supplies the v0.17.3 defaults as well:
`-m 15 15 -k 25,51,101 -g 500 -a 9 -M 4 -O 6 -E 1 -@ 8`.
The published `-D 1` is retained; it overrides the v0.17.3 default `-D 0.5`.

After the protocol amendment is frozen, the runner can be invoked as:

```bash
python studies/panel_read_audit/legacy/regenerate_alignments.py \
  --binary /absolute/scratch/path/locityper_legacy_0.17.3/prefix/bin/locityper \
  --database /absolute/scratch/path/extracted/db \
  --panel-manifest data/panel_read_audit/panel_source_manifest.json \
  --output-database /absolute/scratch/path/extracted/db \
  --manifest-output /absolute/scratch/path/regeneration_manifest.json \
  --timeout-seconds 1800 --max-address-space-gb 16
```

The runner verifies every FASTA against the source manifest before launching
any alignment. It runs loci sequentially, uses all pairs, records commands and
limits, validates complete unique pair coverage and CIGAR presence, and writes
deterministic gzip output atomically. It refuses to overwrite existing PAFs.
Failures remain explicit; parameters never change automatically.

The Bioconda recipe used the v0.17.3 release archive but cloned WFA2 without a
pinned commit. Thus this is a fixed, release-matched build, not proven identical
to the authors' original build. The current 1.7.4 alignment implementation has
material changes beyond default flags and is not an interchangeable fallback.
Accuracy 9 does not certify globally optimal alignments: the historical
backbone and max-gap procedure remains heuristic. Archived numeric outcomes
must pass the study's compatibility checks before regenerated alignments can
support primary read-level interpretation.
