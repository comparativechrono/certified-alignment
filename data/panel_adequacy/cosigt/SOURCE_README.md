# COSIGT Paper

Companion repository for the COSIGT manuscript. It contains the resource
manifests, benchmark outputs, plotting inputs, and plotting scripts used to
generate the manuscript figures.

The COSIGT implementation itself lives in the separate tool repository linked
below. This repository is mainly a compact reproducibility bundle for the paper:
it records which references, assemblies, regions, prediction tables, and figure
scripts were used.

## Study Context

The manuscript evaluates COSIGT for graph-based genotyping of difficult genomic
regions. The bundled results cover:

- Clinically relevant medically actionable gene regions (CMRGs)
- Structural variant-containing gene regions (SVs)
- HLA typing summaries
- Leave-zero-out experiments, where target samples are represented in the graph
- Leave-all-out experiments, where target samples are genotyped against an
  alternate pangenome
- Modern-read, downsampled, and simulated ancient-DNA conditions

Where available, prediction tables compare COSIGT with Locityper. HLA summaries
compare COSIGT and T1K accuracy across the HLA loci included in the paper.

## Repository Structure

- **[plot/](plot)** - Scripts and data for generating all figures (main and supplementary)
- **[resources/](resources)** - Sample lists, assembly identifiers, benchmarking regions, and prediction results

Each directory contains a dedicated `README.md` with detailed information.

## Citation

If you use this code or data, please cite:

XXX

## License

See [LICENSE](LICENSE) for details.

## Related

- [COSIGT tool repository](https://github.com/davidebolo1993/cosigt)
