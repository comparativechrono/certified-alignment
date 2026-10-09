# Native baseline runtime and public read-data feasibility

Status: preparation completed on 9 October 2026. **No fresh genotyping experiment
was performed.** These preparations are separate from retrospective analysis of
published benchmark results and do not establish accuracy or calibration of the
current software.

## A small public matched read/panel dataset is available

The [official Locityper workflow](https://locityper.vercel.app/workflow) uses the
[version 2 test dataset](https://zenodo.org/records/21903460), published on
12 August 2026. Its archive was downloaded and verified against the publisher's
MD5 checksum. Local inventory, sequence headers, and FASTQ lengths were checked.

| Item | Verified content |
| --- | --- |
| Archive | `data.tar.xz`, 107,320,572 bytes |
| Publisher MD5 | `505077b0fea95e5054f0a7ba0b981721` |
| SHA-256 | `f5c199e31723b361a5200e409a35b65e0df98bc7cc3e0fb84322ed48f72fe495` |
| Read donor | HG00733 |
| Paired reads | 878,832 records per mate file; every read is 150 bases |
| NDUFA6 panel | 586 FASTA entries; lengths 25,220–25,281 bases |
| NDUFB1 panel | 588 FASTA entries; lengths 25,550–25,898 bases |
| HG00733 truth sequences | Both named haplotypes occur in each panel |
| NDUFA6 truth lengths | HG00733.1: 25,226; HG00733.2: 25,251 bases |
| NDUFB1 truth lengths | Both HG00733 haplotypes: 25,589 bases |
| NDUFA6 target | CHM13 chr22:42,555,116–42,580,369, BED coordinates |
| NDUFB1 target | CHM13 chr14:86,335,661–86,361,545, BED coordinates |

The archive contains the paired FASTQs, Brotli-compressed local panels, target
BED files, and assembly aliases. The documentation identifies the panel sources
as HPRC2 and HGSVC3. It contains no ready-made Locityper target database, reference
k-mer counts, or preprocessed background distributions. Presence of named truth
haplotypes is not independent validation of their assembly quality, and FASTQ
inventory alone does not establish recruitment completeness.

This is an engineering dataset with one donor and two loci. It is unsuitable for
estimating generalization across donors and loci, and was not selected for
50–1,000 bp indel allele absence. The published workflow already illustrates
very close leave-one-out predictions. We should not select this demonstration
as evidence of biological headroom.

The machine-readable inventory is
[`../studies/panel_adequacy/runtime/locityper_test_v2_inventory.json`](../studies/panel_adequacy/runtime/locityper_test_v2_inventory.json).

## Current native executables are runnable

Locityper's [official installation guide](https://locityper.vercel.app/install)
provides a Docker image. The image was resolved to its Linux/amd64 manifest,
downloaded, and every layer size and SHA-256 digest verified. Binaries and shared
libraries were extracted into a scratch directory. Dynamically linked binaries
were launched through the bundled musl loader; the statically linked Strobealign
binary runs directly. No Docker daemon, privilege escalation, or global software
installation was needed.

Pinned platform manifest:

`docker.io/tprodanov/locityper@sha256:a9575354804c47be19b4cb5da74f551d3f018ebaaef86da501808ae20d4abbae`

| Executable | Observed version | Verification |
| --- | --- | --- |
| Locityper | 1.7.4 | Version and genotype help completed |
| Strobealign | 0.18.0 | Version completed |
| Samtools / HTSlib | 1.19 / 1.19 | Version completed |
| Jellyfish | 2.3.1 | Version completed |
| Minimap2 | 2.31-r1302 | Version completed |

The compressed image layers total 37,087,565 bytes. Binary SHA-256 values,
resolved manifest, release-tag commit and observed version outputs are saved in
[`../studies/panel_adequacy/runtime/locityper_current_runtime.json`](../studies/panel_adequacy/runtime/locityper_current_runtime.json).
The Locityper release tag `v1.7.4` resolves to
`7881f454517552cae8d8d62e6f6f91bb7f06375d`. The container is pinned by its own digest;
the tag's commit alone is not proof of the exact source used in that image.

Version checks establish runtime compatibility, not successful end-to-end
genotyping. The bundled Samtools reports no curl/S3/GCS support, so remote
indexed BAM/CRAM slicing cannot be assumed available through that binary.

Reproduce the verified preparation on Linux x86_64 with Python 3.12 or newer:

```bash
python3 studies/panel_adequacy/runtime/prepare_locityper_runtime.py \
    --destination /absolute/path/outside/the/repository
```

This checks immutable layer/archive hashes, prepares executable wrappers, runs
version checks, and extracts the demonstration data. It never calls genotypes.
Use `--skip-data` to prepare only the runtime. The compressed download is about
145 MB before extraction.

## Remaining prerequisites for a fresh biological run

The official workflow constructs whole-reference Jellyfish k-mer counts and
background depth/error distributions. It downloads a roughly 936 MB compressed
CHM13 reference, then runs:

```bash
samtools faidx CHM13.fa
jellyfish count --canonical --lower-count 2 --out-counter-len 2 --mer-len 25 \
    --threads 8 --size 3G --output CHM13.jf CHM13.fa
locityper target -d db -r CHM13.fa -j CHM13.jf -L data/panels/targets.bed
locityper augment -d db
locityper preproc -i data/HG00733.1.fastq.gz data/HG00733.2.fastq.gz \
    -o preproc/HG00733 -r CHM13.fa -j CHM13.jf
```

The provided local panels avoid needing the separate 3.1 GB assembly archive
used to demonstrate haplotype extraction. A future full-versus-excluded-panel
engineering check can then run `locityper genotype` on the same reads/database
with fixed seed and thread count, once a separate protocol is frozen. The
documented `--leave-out HG00733` removes that donor's named candidates; it does
not by itself remove relatives or equivalent allele copies under different
names. The program's own extended help explicitly warns that identical
haplotypes with other names can remain. Natural donor exclusion and deliberate
allele-family removal therefore require different controls.

The current preparation intentionally did not substitute locus-only k-mer
counts or guessed depth/error models for the required whole-reference and
background information. Those substitutions could change confidence and
recruitment behavior, which is the object of this study.

The [official benchmark-data page](https://locityper.vercel.app/available_data)
also links [Zenodo record 15879519](https://zenodo.org/records/15879519):

| File | Compressed bytes | Potential use |
| --- | ---: | --- |
| `locityper-benchmarking.tar.gz` | 14,709,179 | Published native confidence and truth evaluations; acquired separately for the retrospective baseline |
| `bg-1kgp.tar.gz` | 7,962,713 | Published background models |
| `db.tar.gz` | 785,474,256 | Prepared locus databases and haplotype alignments |
| `1kgp_genotypes.csv.gz` | 13,511,857 | Published cohort predictions |

Both demo loci and HG00733 occur in the published full/leave-one-out Illumina
benchmark keys. However, reusing that older database/background with the newer
demonstration FASTQs would require checking panel boundaries, reference,
coverage, preprocessing provenance, and format compatibility. Their mere
availability does not establish that the combination is a valid fresh run.

## COSIGT feasibility

The official [COSIGT source](https://github.com/davidebolo1993/cosigt) was inspected
at `3b8a905483128f1a0aed718430f2491267c33214`. Its core scorer is a Go program,
separate from the Snakemake workflow. It consumes observed node-coverage vectors,
haplotype traversal vectors, and cluster assignments, and writes absolute cosine
scores and ranked candidate combinations. The source identifies its scorer as
version 0.2 and notes behavior changes from 0.1.x, so published and current
outputs must not be silently pooled.

The [official example](https://davidebolo1993.github.io/cosigtdoc/usecases/usecases.html)
uses two 1000 Genomes samples and C4/GYPA loci. It starts from whole-genome CRAMs,
a complete reference, and chromosome-scale HPRC graphs. The full workflow needs
graph construction/alignment and coverage tools in addition to the small scorer.
Published truth-comparison tables cannot be treated as substitutes for its
native cosine scores or read-coverage evidence. The present preparation does
not include a fresh COSIGT run or an assembly-caller comparison.

## Consequence for this milestone

The downloaded published Locityper benchmark offers a substantially broader
immediate baseline than a two-locus rerun, because it preserves native deployable
confidence columns alongside separately identifiable truth-derived outcomes.
Current native reruns are technically feasible but remain a separate follow-up.
Neither sequence-error filtering in the published archive nor the prepared
two-locus demonstration alone resolves missing structural allele detection.
