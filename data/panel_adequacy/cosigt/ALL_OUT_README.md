# Leave-All-Out

Prediction results for samples whose matching haplotypes are absent from the
pangenome graph.

## Structure

- **cmrgs_predictions/** - HPRC Year 1 samples genotyped on HGSVCv3-based graphs
- **svs_predictions/** - HGSVCv3 samples genotyped on HPRC Year 1-based graphs

## Files

- `cosigt_cmrgs_hprcy1_samples_on_hgsvcv3_assemblies.tsv.gz` - HPRC Year 1 CMRG samples genotyped against HGSVCv3 assemblies
- `cosigt_cmrgs_hprcy1_AMRnoPUR_samples_on_hgsvcv3_nonAMR_assemblies.tsv.gz` - Population-subset CMRG table used for cross-population analysis
- `cosigt_svs_hgsvcv3_samples_on_hprcy1_assemblies.tsv.gz` - HGSVCv3 SV samples genotyped against HPRC Year 1 assemblies

## File Format

Tab-separated files with the following columns:

- `sample`, `region`, `gene_name` - Sample and region identifiers
- `hap_1_pred`, `hap_2_pred` - COSIGT predicted haplotypes from the alternate pangenome
- `QV_1_pred`, `QV_2_pred`, `QV_sum_pred` - Phred-scaled prediction quality values
- `hap_1_best`, `hap_2_best` - Best available alternate-pangenome haplotypes used as the lower-bound comparator
- `QV_1_best`, `QV_2_best`, `QV_sum_best` - Quality values for those best matches
- `error_rate_*_pred` - Error rates for predictions
- `error_rate_*_best` - Error rates for best matches
- `error_rate_diff_pred_minus_best` - Difference between predicted and best error rates
- `QV_diff_pred_minus_best` - Difference between predicted and best QV

Haplotype predictions are written as `SAMPLE#HAPLOTYPE#CONTIG:START-END`.
