# CMRGs

Clinically relevant medically actionable gene regions used for COSIGT
benchmarking.

## Files

- `cmrgs.hprcy1.bed` - Original 326-gene CMRG target set adapted from the Locityper paper/database
- `cmrgs_refined.hprcy1.bed` - Same 326 genes after the COSIGT region-refinement step

## Description

The refined BED is the region set used by the COSIGT benchmark tables and plots.
Both files use reference-style chromosome coordinates and gene names.

**Format**: BED file with 4 columns (chromosome, start, end, gene name)
