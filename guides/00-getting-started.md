---
title: "Getting started"
subtitle: "Software, reference genomes, conventions and example data for every guide"
---

Read this once before any guide. It covers the set-up every guide shares.

::: {.callout-tip}
## How to read these guides
Each guide follows one workflow from raw reads to interpretation: design, QC, processing,
statistics, plots and inference. Where there is a common alternative, a collapsed box names it
and its trade-offs. Treat each guide as a documented starting point, not the only correct way.
Read the linked tool documentation, and adapt thresholds and parameters to your own data.
:::

| Guide | Question | Processing | Analysis in R |
|---|---|---|---|
| [RNA-seq](01-rna-seq.md) | Which genes change expression? | `pipelines/rnaseq` | DESeq2, enrichment |
| [ATAC-seq](02-atac-seq.md) | Which regulatory elements open or close? | `pipelines/chip_atac` (`assay: atac`) | Consensus peaks, DESeq2, ChIPseeker, motifs |
| [ChIP-seq](03-chip-seq.md) | Where does a protein bind, or a histone mark sit, and does it change? | `pipelines/chip_atac` (`assay: chip`) | DiffBind, ChIPseeker |
| [Hi-C](04-hic.md) | How is the genome folded, and does folding change? | `pipelines/hic` | Contact matrices, APA, plotgardener |
| [scRNA-seq](05-scrna-seq.md) | Which cell types exist, and how do they differ? | `pipelines/scrna` | OSCA/Bioconductor, pseudobulk |
| [Comparative genomics](appendix-comparative-genomics.md) | How do you compare species? | `pipelines/consensus_genome` | liftOver, orthologs, outgroups |

## Hardware

| Task | Minimum sensible machine |
|---|---|
| R analysis of processed data (all guides) | Laptop, 16 GB RAM |
| RNA-seq / ATAC / ChIP processing | Workstation: 16 cores, 64 GB RAM, 1 TB disk (STAR needs ~32 GB for human) |
| Hi-C processing, Cell Ranger | 32+ cores, 64–128 GB RAM, several TB of scratch, or an HPC cluster |

All pipelines run the same way on a workstation or a cluster; see
[Processing pipelines](../pipelines/README.md).

## Software

```bash
# 1. Miniforge (conda + mamba): https://github.com/conda-forge/miniforge
# 2. Snakemake, in its own environment
mamba create -n snakemake -c conda-forge -c bioconda snakemake=9 && conda activate snakemake
# 3. The guides' repository
git clone https://github.com/Delta-43/comparative-genomics-guides.git
# 4. R + Bioconductor for the analysis sections (or install R from CRAN and use BiocManager)
mamba env create -n guides-r -f comparative-genomics-guides/pipelines/envs/r_analysis.yaml
```

Pipelines create their own tool environments on first run (`--use-conda`). If you prefer a
plain R installation:

```r
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "apeglm", "airway", "clusterProfiler", "org.Hs.eg.db",
                       "ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene", "DiffBind", "Rsubread",
                       "scater", "scran", "scDblFinder", "SingleR", "celldex", "TENxPBMCData",
                       "batchelor", "speckle", "DropletUtils", "plotgardener", "biomaRt"))
install.packages(c("pheatmap", "ggplot2", "data.table", "harmony"))
```

The R code in these guides was run with R 4.5 and Bioconductor 3.22 (e.g. DESeq2 1.50,
DiffBind 3.20, scran 1.38). Older releases mostly work, but function arguments occasionally change.

## Reference genomes and annotation

Use a genome FASTA and gene annotation (GTF) **from the same source and release**:

| Source | Human example | Chromosome names |
|---|---|---|
| Ensembl | `Homo_sapiens.GRCh38.dna.primary_assembly.fa` + `Homo_sapiens.GRCh38.<release>.gtf` | `1`, `2`, … `MT` |
| GENCODE | `GRCh38.primary_assembly.genome.fa` + `gencode.v<N>.primary_assembly.annotation.gtf` | `chr1`, … `chrM` |
| UCSC | `hg38.fa` + UCSC/RefSeq GTF, chain files for liftOver | `chr1`, … `chrM` |
| 10x Genomics | `refdata-gex-GRCh38-2024-A` (prebuilt for Cell Ranger) | `chr1`, … |

**Chromosome naming is the most common silent failure.** The FASTA, index, GTF, blacklist,
chain files and chrom.sizes must all use the same style. Check before every run:

```bash
grep '^>' genome.fa | head -3                     # FASTA
grep -v '^#' genes.gtf | cut -f1 | sort -u | head -3    # GTF
cut -f1 blacklist.bed | sort -u | head -3         # blacklist
```

Other useful resources: [ENCODE blacklists](https://github.com/Boyle-Lab/Blacklist) (Amemiya
*et al.* 2019), [UCSC liftOver chains](https://hgdownload.soe.ucsc.edu/downloads.html), and
[deepTools effective genome sizes](https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html).

## Conventions used in the guides

- Code blocks marked `bash` run in a terminal; blocks marked `r` run in R.
- R examples use **public example datasets** that Bioconductor downloads for you (airway,
  PBMC 3k, DiffBind's tamoxifen data), so every block runs as written. Where a block reads your
  own pipeline output instead, the file paths match what the pipelines write (`results/...`).
- Collapsed "Alternatives" boxes describe other standard tools and their trade-offs.
- Each guide ends with a **reproducibility checklist** and **exercises**; answers are in a
  collapsed box.

## Reproducibility checklist (every guide)

1. **Record versions:** pipeline commit, conda environment files, genome and annotation releases,
   chain files. Save `snakemake --report report.html` and R's `sessionInfo()`.
2. **Set seeds** before every random step: random background regions, UMAP/t-SNE, k-means,
   `quickCluster`, doublet simulation, permutation tests (`set.seed(1)`).
3. **Dry-run first:** `snakemake -n` lists every job and command.
4. **Keep raw counts raw.** Normalise inside the statistical method. TPM, CPM and RPGC values are
   for plots and browsers only.
5. **Write down every exclusion:** which samples, cells or regions were removed, and why.
6. **Use relative paths** in scripts (e.g. `here::here("results", ...)`), never
   machine-specific absolute paths.

## Example data

The processing pipelines can be tried on a published, openly available study of primate
astrocyte evolution: **Ciuba *et al.* 2025**, "Molecular signature of primate astrocytes reveals
pathways and regulatory changes contributing to human brain evolution", *Cell Stem Cell*
32:426–444 ([doi:10.1016/j.stem.2024.12.011](https://doi.org/10.1016/j.stem.2024.12.011)). It has
bulk RNA-seq, ATAC-seq, H3K27ac and H3K4me3 ChIP-seq from human, chimpanzee and rhesus macaque
iPSC-derived astrocytes, and in situ Hi-C from human and chimpanzee iPSC-derived astrocytes and
human primary fetal astrocytes:

| Assay | ArrayExpress | ENA project |
|---|---|---|
| RNA-seq | E-MTAB-13252 | ERP150223 |
| ATAC-seq | E-MTAB-13253 | ERP150225 |
| ChIP-seq H3K27ac | E-MTAB-13254 | ERP150231 |
| ChIP-seq H3K4me3 | E-MTAB-13255 | ERP150220 |
| Hi-C | E-MTAB-13259 | ERP150235 |

The run tables are in `pipelines/demo_data/`. One command turns a table into a sample sheet and a
checksum-verified download script:

```bash
cd comparative-genomics-guides/pipelines
python demo_data/make_samples.py demo_data/E-MTAB-13253_ATAC-seq.ena.tsv --pipeline chip_atac --outdir ~/data/atac
bash ~/data/atac/download.sh
```

Check the size before downloading: the RNA-seq, ATAC and ChIP runs are 1–6 GB each; the Hi-C runs
are 150+ GB each. The guides quote a few of the study's published results as examples of what each
analysis produces. They are there for illustration, not as numbers you should expect to reproduce
exactly.

::: {.callout-note}
## Example configs
Each pipeline's `config.yaml` is a generic single-species default. For this study's
cross-species set-up (RNA-seq on a consensus genome; ATAC, ChIP and Hi-C on each species' own
genome), copy the complete configs in `pipelines/<assay>/examples/`, as described in
[Processing pipelines](../pipelines/README.md#example-profiles).
:::
