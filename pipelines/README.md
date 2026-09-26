---
title: "Processing pipelines"
---

Portable [Snakemake](https://snakemake.readthedocs.io) workflows that turn raw FASTQ files into
the processed files every guide starts its analysis from. Each one runs the same way on a
laptop (small data), a workstation or an HPC cluster.

| Pipeline | Input | Produces | Guide |
|---|---|---|---|
| [`rnaseq`](https://github.com/Delta-43/comparative-genomics-guides/tree/main/pipelines/rnaseq) | Paired-end RNA-seq FASTQ | Gene count matrix, strandedness call, BPM bigWigs, MultiQC report | [RNA-seq](../guides/01-rna-seq.md) |
| [`chip_atac`](https://github.com/Delta-43/comparative-genomics-guides/tree/main/pipelines/chip_atac) | Paired-end ChIP-seq or ATAC-seq FASTQ | Filtered BAMs, MACS2 narrow/broad peaks (per sample and per group), RPGC bigWigs, FRiP / fingerprint / TSS-enrichment QC | [ATAC-seq](../guides/02-atac-seq.md), [ChIP-seq](../guides/03-chip-seq.md) |
| [`hic`](https://github.com/Delta-43/comparative-genomics-guides/tree/main/pipelines/hic) | Paired-end Hi-C FASTQ | `.hic` file, per-chromosome contact dumps, TopDom domains, loops | [Hi-C](../guides/04-hic.md) |
| [`scrna`](https://github.com/Delta-43/comparative-genomics-guides/tree/main/pipelines/scrna) | 10x Genomics FASTQ (+ antibody capture) | Cell Ranger reference and filtered feature-barcode matrix | [scRNA-seq](../guides/05-scrna-seq.md) |
| [`consensus_genome`](https://github.com/Delta-43/comparative-genomics-guides/tree/main/pipelines/consensus_genome) | UCSC assemblies + chain files | A reference with cross-species divergent sites masked | [Comparative genomics](../guides/appendix-comparative-genomics.md) |

::: {.callout-note}
## Testing status
- **`rnaseq` and `chip_atac` have been run end to end** (both `atac` and `chip` modes, with an input
  control) on small simulated datasets, using the conda environments in `envs/`. On a spliced,
  reverse-stranded simulation, the RNA-seq pipeline inferred `-s 2`, and DESeq2 on its counts
  recovered exactly the genes simulated as changed. The ATAC run removed all mitochondrial
  reads and duplicates, and recovered 40 of 40 simulated peaks in the pooled call, with none
  called elsewhere.
- `hic`, `scrna` and `consensus_genome` are dry-run tested in all their modes (their environments
  build). The masking script in `consensus_genome` has unit tests on synthetic alignments. They need
  real-size inputs (Juicer, Cell Ranger, UCSC genomes), which weren't run here.
- None of the pipelines has yet been run on a full-size real dataset. If you do, please report
  problems through the repository's issue tracker.
:::

## Install once

```bash
# Miniforge gives you conda + mamba: https://github.com/conda-forge/miniforge
mamba create -n snakemake -c conda-forge -c bioconda snakemake=9
conda activate snakemake
git clone https://github.com/Delta-43/comparative-genomics-guides.git
cd comparative-genomics-guides/pipelines
```

Every rule declares its own software environment in [`envs/`](https://github.com/Delta-43/comparative-genomics-guides/tree/main/pipelines/envs).
With `--use-conda`, Snakemake builds these on first use, so you never install the tools by hand.
Two tools can't come from conda and are handled separately:

- **Cell Ranger** needs a free download from 10x Genomics (licence click-through). Set its path
  in `scrna/config.yaml`.
- **Juicer**'s scripts are cloned from GitHub at a pinned commit by the `hic` pipeline itself.

## Run

```bash
cd rnaseq                              # or chip_atac, hic, scrna, consensus_genome
$EDITOR config.yaml samples.tsv        # reference files, options, your samples
snakemake --use-conda --cores 16 -n    # dry run: lists every job and command. Always do this first.
snakemake --use-conda --cores 16       # run
snakemake --report report.html         # provenance report: software versions, commands, runtimes
```

Snakemake only reruns steps whose inputs changed, so after an interruption just run the same
command again.

### Workstation or cluster?

`--cores` caps total CPU use; per-rule `threads` and `mem_mb` come from `config.yaml`.
Approximate memory needs:

| Step | Memory | Notes |
|---|---|---|
| STAR genome index (human) | 32–40 GB | Built once per genome |
| STAR alignment | ~32 GB per job | Use `--resources mem_mb=<total RAM>` so parallel jobs fit |
| BBSplit (two mammalian genomes) | ~96 GB | Mixed-species samples only |
| bowtie2, MACS2, deepTools | 4–8 GB | |
| Juicer (per sample) | 16–32 GB | Plus ~150 GB scratch per billion read pairs |
| Cell Ranger | ≥ 64 GB | |

On an HPC cluster, install an executor plugin and submit through it; nothing in the Snakefiles
changes:

```bash
pip install snakemake-executor-plugin-slurm
snakemake --use-conda --executor slurm --jobs 100 --default-resources slurm_account=<acct> runtime=720
```

## Using the example data

[`demo_data/`](https://github.com/Delta-43/comparative-genomics-guides/tree/main/pipelines/demo_data)
holds the ENA run tables for the example study used throughout the guides (Ciuba *et al.* 2025,
ArrayExpress E-MTAB-13252 to 13259). `make_samples.py` turns a run table into a ready-to-use
`samples.tsv` plus a resumable, checksum-verified download script:

```bash
python demo_data/make_samples.py demo_data/E-MTAB-13253_ATAC-seq.ena.tsv \
    --pipeline chip_atac --outdir ~/data/example_atac
bash ~/data/example_atac/download.sh           # downloads FASTQ and checks md5 sums
```

For your own ENA/SRA study, download its run table with the ENA portal API (fields
`run_accession,sample_alias,library_layout,read_count,fastq_md5,fastq_ftp`) and provide a
`--species-map` if it has several species.

## Design choices and defaults

Defaults follow widely used community practice (ENCODE, nf-core, HBC training materials). Every
choice below is a config option; the guides discuss the pros and cons of each.

**All pipelines**

- Trimming with Trim Galore (Cutadapt) at its default settings (Phred ≥ 20, reads ≥ 20 bp),
  followed by FastQC before and after, and MultiQC to collect everything into one report.
- Relative paths in `samples.tsv` are resolved against the sample sheet's own folder.
- Chromosome naming must match across FASTA, GTF, blacklist and chain files (`chr1` vs `1`).
  The mitochondrial and unplaced-contig filter recognises both styles.

**RNA-seq**

- Strandedness is **inferred from the data** (STAR `--quantMode GeneCounts` gives counts for
  both strand orientations), not assumed. A mismatch with your config triggers a warning.
- featureCounts counts read pairs (`-p --countReadPairs`) on exons, summed per `gene_id`.
  Multi-mapping reads and reads overlapping two genes are not counted by default; `-O` / `-M`
  change that at the cost of some ambiguity.

**ChIP-seq / ATAC-seq**

- Kept: properly paired reads with MAPQ ≥ 30, primary alignments only. PCR duplicates are
  removed with `samtools markdup` (pair-aware), as are mitochondrial and unplaced contigs and
  ENCODE blacklist regions (if a blacklist is given).
- ATAC: reads are shifted +4/−5 bp for the Tn5 insertion offset (`alignmentSieve --ATACshift`).
- MACS2 runs in paired-end mode (`-f BAMPE`), so fragment sizes come from the read pairs, not a
  model. ChIP samples use their input/IgG control if one is listed. Peaks are called per sample
  and on the pooled BAM of each `group`.
- QC: FRiP, deepTools fingerprint, insert-size histogram, TSS enrichment (if a GTF is given).

**Mixed-species samples** (co-cultures, xenografts, two species pooled in one experiment)

- List two genomes in `config.yaml`. Reads are first assigned to a species with BBSplit, then
  each species is processed on its own genome. Reads matching both genomes equally well are
  discarded (`ambiguous2=toss`). That is the safe choice for closely related species, where
  many reads are identical in both genomes.
