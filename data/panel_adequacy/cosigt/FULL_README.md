# Leave-Zero-Out

Prediction results for samples whose matching haplotypes are present in the
pangenome graph.

## Structure

- **cmrgs_predictions/modern/** - CMRG COSIGT-vs-Locityper comparison tables at 1X, 2X, 5X, and 30X, plus a COSIGT VCF at 30X
- **cmrgs_predictions/aDNA/** - Simulated ancient-DNA CMRG predictions at 1X and 2X with 0% or 10% contamination
- **svs_predictions/** - SV-region predictions for HGSVCv3 samples at 30X

The aDNA directory contains COSIGT-vs-Locityper comparison tables and
COSIGT-only tables for the `anc5k_split20k` and `anc20k_split100k` simulation
settings used in the paper.

## File Naming

Comparison tables follow this broad pattern:

- `cosigt_locityper_<target>_<graph>_<condition>.tsv.gz` for COSIGT and Locityper comparisons
- `cosigt_<target>_<graph>_<condition>.tsv.gz` for COSIGT-only aDNA summaries
- `cosigt_<target>_<graph>_<coverage>.vcf.gz` for COSIGT VCF output where included

Common condition tokens include coverage (`1X`, `2X`, `5X`, `30X`),
contamination (`cont0pct`, `cont10pct`), and aDNA demographic labels such as
`anc5k_split20k` or `anc20k_split100k`.

## File Format

Comparison TSVs contain these core columns:

- `sample`, `region`, `gene_name` - Sample and region identifiers
- `hap_1_pred_locityper`, `hap_2_pred_locityper` - Locityper haplotype predictions, where available
- `hap_1_pred_cosigt`, `hap_2_pred_cosigt` - COSIGT haplotype predictions
- `QV_1_*`, `QV_2_*`, `QV_sum_*` - Phred-scaled quality values
- `error_rate_1_*`, `error_rate_2_*`, `avg_error_rate_*` - Per-haplotype and average error rates
- `error_rate_diff_cosigt_minus_locityper` - COSIGT minus Locityper error-rate difference, where available
- `QV_diff_cosigt_minus_locityper` - COSIGT minus Locityper QV difference, where available

Haplotype predictions are written as `SAMPLE#HAPLOTYPE#CONTIG:START-END`.
