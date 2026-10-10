# Regional read acquisition

`resolve_sources.py` resolves the donor-specific `Illumina/child/*.final.cram`
objects by listing the public HPRC bucket linked in the Locityper paper. It records
object size, ETag, last-modified metadata, the published whole-CRAM MD5 text, and
computed hashes of the downloaded CRAI and MD5 sidecars. The whole CRAM is not
downloaded or independently hashed.

`acquire_slices.py` uses pysam 0.24.1 / HTSlib 1.24 to fetch indexed CRAM records
from each historical database `ref.bed` interval, expanded by 1,000 bases on each
side and clamped to the chromosome. Overlapping windows are merged. The original
GRCh38 full analysis set plus decoy/HLA reference is supplied explicitly as a
local FASTA containing the required full chromosomes. Each chromosome sequence
is checked against the official dictionary M5; normal CRAM reference checksum
checks remain enabled. A fresh remote CRAM handle is opened for each locus.

No mapping-quality, proper-pair, duplicate, QC, secondary or supplementary filter
is applied during acquisition. Bases, qualities, alignment flags and mapping
qualities are retained in the BAM. Downstream evidence counting applies its
pre-specified filters. Coordinate-sorted BAMs and indexes have local SHA-256
checksums; each record in the manifest includes reference sequence names,
lengths and the CRAM header M5 values. Read-group donor names are checked.

These are donor-matched public Illumina data from the paper's cited source. Exact
raw-library equivalence to the historical caller inputs is not asserted.

Reads whose stored coordinates fall outside these windows are not recruited.
In particular, this does not search globally unmapped reads or recover mates at
other coordinates. The metadata count primary paired query names lacking both
ends in the extracted records. Positive sequence evidence may therefore be
informative; zero observations do not establish biological absence.

Example from the repository root, with paths adjusted to the local cache:

```bash
python studies/panel_read_audit/reads/resolve_sources.py \
  --cases data/panel_read_audit/cases.tsv \
  --controls data/panel_read_audit/controls.tsv \
  --output data/panel_read_audit/read_sources.json
python studies/panel_read_audit/reads/acquire_slices.py \
  --cases data/panel_read_audit/cases.tsv \
  --controls data/panel_read_audit/controls.tsv \
  --sources data/panel_read_audit/read_sources.json \
  --db /path/to/verified/db \
  --reference-fasta /path/to/verified/selected_GRCh38.fa \
  --output /path/to/read_cache/bams \
  --cache /path/to/read_cache/indexes \
  --manifest results/panel_read_audit/read_manifest.json \
  --flank 1000 --workers 4 \
  --protocol-commit 986336d1308ef903123fcbc52f0084e1659c972d
```

The BAM paths in the manifest identify the original execution workspace. Re-run
acquisition to populate paths in another workspace; source URLs and interval
hashes remain the reproducible identifiers.

The first acquisition attempt used a remote reference URL. HTSlib fetched whole
chromosomes and a reference transport failure left poisoned CRAM handles that
subsequently yielded empty records. All BAMs from that attempt were set aside;
none enter read-evidence counting. The corrected attempt uses a verified local
reference, a fresh CRAM handle per locus, rejects every zero-record extraction as
unresolved, and independently re-reads every local BAM to check its record count.
The original failure log and an explicitly unusable attempt manifest are kept.
Final successes carry `status=verified` and `acquisition_verified=true`.
