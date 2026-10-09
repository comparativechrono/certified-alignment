# Benchmarking

Benchmark regions, exclusion masks, HLA resources, and prediction tables used to
evaluate COSIGT in the manuscript.

## Structure

- **[cmrgs/](cmrgs)** - Clinically-relevant medically-actionable gene regions
- **[svs/](svs)** - Structural variant regions
- **[flagger/](flagger)** - Regions excluded from benchmarking due to assembly errors
- **[hla_typing/](hla_typing)** - HLA typing comparisons and regions
- **[leave_zero_out/](leave_zero_out)** - Leave-zero-out prediction results (COSIGT and Locityper where applicable)
- **[leave_all_out/](leave_all_out)** - Leave-all-out prediction results (COSIGT)

Each subdirectory contains a dedicated `README.md` with detailed information.

## Benchmarking Strategies

### Leave-zero-out

Samples are genotyped against graphs that include matching haplotypes. These
tables cover modern-read conditions at 1X, 2X, 5X, and 30X, plus simulated
ancient-DNA conditions at 1X and 2X with 0% or 10% contamination.

### Leave-all-out

Samples are genotyped against graphs built from an alternate pangenome, so the
exact target haplotypes are absent. Tables include the COSIGT prediction and the
best available alternate-pangenome match used as the lower-bound comparator.
