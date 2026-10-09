# Locityper loci database, 1KGP genotypes and benchmarking data

This collection contains [Locityper](https://github.com/tprodanov/locityper) `v0.17`-`v0.18` results.

## Locityper genotypes for the 1000 genomes project (1KGP)

Genotyping results can be found in the `1kgp_genotypes.csv.gz` file with the following columns:
- `sample`, `locus`,
- `genotype`: Primary genotype prediction,
- `quality`: Phred-scale probability of a non-primary genotype (possibly very similar) being correct,
- `total_reads`: How many reads were recruited for genotyping this locus,
- `unexpl_reads`: How many reads were not explained by the primary genotype,
- `weight_dist`: Weighted distance between the primary genotype prediction and other genotypes.
    Here, distance was calculated as the Jaccard distance between multisets of (15,15)-minimizers.
    Distances are weighted according to the secondary genotype probabilities.
- `warnings`: Other heuristic warnings.

Genotypes were calculated using Locityper `v0.17.3` based on 3202 Illumina WGS datasets.

## Target loci database

A database of 326 target loci can be found in `db.tar.gz` (challenging medically relevant genes + CFH + MUC + CYP2 genes).
For each locus, the following files can be useful:
- `ref.bed`: reference coordinates of the locus,
- `discarded_haplotypes.txt`: haplotypes, identical to others and discarded from futher analysis,
- `haplotypes.fa.gz`: all unique HPRC haplotypes for the locus,
- `haplotypes.paf.gz`: pairwise alignments of all haplotypes to each other.

## Background parameters for the full 1KGP cohort

Locityper-preprocessed background parameters for 3202 Illumina WGS datasets can be found in `bg-1kgp.tar.gz`.
After decompressing, corresponding directories used during genotyping as `locityper genotype ... -p bg-1kgp/SAMPLE`.

## Locityper benchmarking

Locityper benchmarking results can be found in the `locityper-benchmarking.tar.gz` archive.
Except for the trio comparison, benchmarking was performed based on the HPRC samples using Locityper `v0.17.3`.

`QV` subdirectory contains 9 files, of them 4 `full_*` files with Locityper evaluation using full reference panel,
`loo_*` files for Locityper in the leave-one-out (LOO) configuration,
and `nygc.csv.gz` for evaluating NYGC 1KGP call set haplotypes.

Each file has three lines for each sample and locus: based on the `query_type` column, two `hap` lines represent haplotype evaluation and `gt` line shows diploid genotype evaluation.
In addition to previously explained columns, the files may have following columns:
- `edit`, `size`: edit distance and alignment size in the haplotype/genotype alignment,
- `div = edit / size`,
- `qv`: Phred-like transformation of `div`,
- `avail_*`: in the LOO-setting, what is the similarity between the actual and the closest available haplotypes/genotypes in the LOO database.

### Trio concordance

Trio concordance based on the 1KGP samples can be found in the `trio_conc.csv.gz` file.
For each child (`indiv` column) and each locus the file shows concordance between child-parent haplotypes/genotypes.
Last (`expl`) column explains which child-parent haplotype pairs were matched together.

### Benchmarking HLA/KIR calling

HLA/KIR results can be found in the `HLA` subfolder.
`immuannot.csv.gz` contains baseline HLA/KIR annotation for the HPRC samples.
`alleles1` and `alleles2` represent alleles on the first and second haplotype of the corresponding samples.
If a haplotype contains a novel allele, it can be annotated with several close alleles, combined via `|`.

Three other files evaluate HLA/KIR calling accuracy, where each sample and locus can be written several times, according to the number of haplotypes.

Accuracy column shows how many allele fields were matched, while `base_length` shows how many fields there are in the baseline allele.
Finally, `avail` column shows if this the first two baseline fields are available in the LOO setting.

### Variant calling comparison

Variant calling comparison between Locityper `v0.15.1`, 1KGP and Pangenie `v3.02` can be found in the `varcall.csv.gz` file.
There, number of true/false positive/negative variant calls are shown for each combination of locus,tool,sample, variant type and quality threshold (0/1/10).

## Simulated data

File `sim_seeds.txt.gz` contains ART illumina v2.5.8 seeds (in hex form) for each of the pangenome haplotypes using the following command:
```bash
art_illumina -q -ss HS25 -m 500 -s 20 -l 150 -f 1 -na -rs $seed -i $in_haplotype -o $out_fastq
```
Each haplotype is run 15 times with different seeds with coverage 1, for diploid samples this produces a total coverage of 30.

