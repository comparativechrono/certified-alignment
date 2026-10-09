# SVs

Structural variant-containing gene regions used for benchmarking.

## Files

- `svs_refined.hgsvcv3.bed` - 265 refined SV-containing gene regions
- `svs_refined.annotation.bed` - SVtype annotations

## Description

This target set was defined for the COSIGT study and refined with the COSIGT
region-refinement workflow. It is used for the HGSVCv3-based SV benchmark tables
and figure panels.

**Format**: BED file with 4 columns (chromosome, start, end, gene name) and BED file with (chromomosome,start,end,gene name, SVtype, SVtype of the longest variant, length of the longest variant, number of samples with the longest variant, SVtype of the 2nd longest variant, length of the 2nd longest variant, number of samples with the 2nd longest variant)
