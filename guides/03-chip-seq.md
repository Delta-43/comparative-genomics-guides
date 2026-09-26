---
title: "ChIP-seq: transcription factor binding and histone marks"
subtitle: "Peak calling, quality control and differential binding"
---

::: {.callout-tip}
## Learning objectives
- Design a ChIP-seq experiment with the right controls, depth and replication.
- Process reads to filtered alignments and narrow or broad peaks.
- Judge enrichment quality (FRiP, fingerprint, cross-correlation, replicate agreement).
- Test for differential binding with DiffBind, annotate peaks and find motifs.
- Integrate binding with expression, and compare binding across species.

**Time:** about 3 hours. **You need:** the set-up from [Getting started](00-getting-started.md).
:::

## 1. What ChIP-seq measures

Chromatin immunoprecipitation (ChIP) uses an antibody to pull down DNA fragments bound by a
protein (a transcription factor such as CTCF) or wrapped in nucleosomes carrying a histone mark
(e.g. H3K27ac at active enhancers and promoters, H3K4me3 at promoters, H3K27me3 over repressed
domains). Sequencing the pulled-down DNA shows where that protein or mark sits across the genome.

Two peak shapes need different treatment:

| Type | Examples | Peak shape | Peak calling |
|---|---|---|---|
| Point-source / narrow | TFs, CTCF, H3K4me3, H3K27ac | Sharp, hundreds of bp | MACS2 default (narrowPeak) |
| Broad | H3K27me3, H3K36me3, H3K9me3, H3K4me1 | Diffuse, kb to hundreds of kb | MACS2 `--broad`, or domain callers (SICER, epic2) |

**It can tell you:** where the protein or mark is enriched, and whether enrichment differs
between conditions. **It can't tell you:** that binding is functional, the absolute occupancy
(without spike-ins), or anything the antibody doesn't specifically recognise.

## 2. Experimental design

Standards below are from ENCODE ([TF](https://www.encodeproject.org/chip-seq/transcription_factor/),
[histone](https://www.encodeproject.org/chip-seq/histone/)).

| Factor | Recommendation |
|---|---|
| Replicates | ≥ 2 biological replicates; ≥ 3 for differential binding |
| Control | An **input** (sonicated chromatin, no antibody) per condition, matched in read length and layout. IgG is an alternative, but usually yields little DNA |
| Depth | TFs and narrow marks: 20 M usable fragments per replicate. Broad marks: 45 M |
| Antibody | Validated for ChIP (knockout/knockdown, or ENCODE characterisation). The antibody is the largest single source of failure |
| Global changes | If a treatment changes the mark globally (e.g. EZH2 inhibitor and H3K27me3), add a spike-in: ChIP-Rx, exogenous *Drosophila* chromatin ([Orlando *et al.* 2014](https://doi.org/10.1016/j.celrep.2014.10.018)) |

::: {.callout-note collapse="true"}
## Alternatives: CUT&RUN and CUT&Tag
CUT&RUN ([Skene & Henikoff 2017](https://doi.org/10.7554/eLife.21856)) and CUT&Tag
([Kaya-Okur *et al.* 2019](https://doi.org/10.1038/s41467-019-09982-5)) tether a nuclease or Tn5
to the antibody in intact nuclei. They need far fewer cells and reads (3–8 M), and have lower
background, so no input is needed. The same pipeline works: use `assay: chip` with no control,
and keep duplicates if the library is low complexity (`remove_duplicates: false`), because
enzyme cutting creates genuine identical fragments.
:::

## 3. Raw-read QC

As in the [RNA-seq guide](01-rna-seq.md#3-quality-control-of-raw-reads). ChIP libraries often
show moderate duplication. High duplication (> 50 %) usually means little starting material.

## 4. From FASTQ to peaks

```bash
cd pipelines/chip_atac
$EDITOR config.yaml samples.tsv    # assay: "chip"; each ChIP row names its input in `control`
snakemake --use-conda --cores 16 -n && snakemake --use-conda --cores 16
```

The same filtering as the [ATAC-seq guide](02-atac-seq.md#4-from-fastq-to-peaks) applies:
MAPQ ≥ 30, duplicates removed, mitochondrial and unplaced contigs removed, blacklist.
Tn5 shifting is off for ChIP. MACS2 then compares each ChIP to its input:

```bash
# what the pipeline runs for each sample (narrow; the broad call adds --broad)
macs2 callpeak -f BAMPE -t chip.bam -c input.bam -g <genome size> -n sample --keep-dup all
```

Outputs per sample: `*_peaks.narrowPeak` (peaks with summit position, fold enrichment,
−log10 p and q) and `*_peaks.broadPeak`, plus peaks from the pooled replicates of each `group`.

::: {.callout-warning}
## Peaks without a control
Without an input, MACS2 uses a flat local background. Open chromatin, copy-number gains and
"hyper-ChIPable" regions then turn up as false peaks. If there is no input, use the blacklist,
be stricter on q-values, and never compare peak *counts* between samples with and without
controls.
:::

## 5. Quality control

| Metric | What it measures | Guide value |
|---|---|---|
| Usable fragments | Depth after filtering | ≥ 20 M (narrow) / ≥ 45 M (broad), ENCODE |
| NRF / PBC1 / PBC2 | Library complexity (fraction of distinct fragments) | > 0.9 / > 0.9 / > 10 (ENCODE) |
| FRiP | Fraction of reads in peaks: signal-to-noise | No fixed ENCODE cut-off. TFs often 1–10 %, H3K4me3/H3K27ac 10–50 %. Compare like with like |
| NSC / RSC | Strand cross-correlation (phantompeakqualtools) | NSC > 1.05, RSC > 0.8 ([Landt *et al.* 2012](https://doi.org/10.1101/gr.136184.111)) |
| Fingerprint | Cumulative coverage (`results/qc/fingerprint/`) | ChIP curves bend sharply away from the input's diagonal |
| Replicate agreement | Peak overlap, IDR, correlation of signal | Most strong peaks shared. IDR for TFs (ENCODE) |

Replicate overlap in R:

```r
library(GenomicRanges)
read_peaks <- function(f) {
  x <- read.table(f)[, 1:3]; GRanges(x$V1, IRanges(x$V2 + 1, x$V3))
}
r1 <- read_peaks("results/peaks/narrow/treat_rep1.GRCh38_peaks.narrowPeak")
r2 <- read_peaks("results/peaks/narrow/treat_rep2.GRCh38_peaks.narrowPeak")
c(rep1 = length(r1), rep2 = length(r2),
  rep1_in_rep2 = mean(overlapsAny(r1, r2)), rep2_in_rep1 = mean(overlapsAny(r2, r1)))
```

A strong experiment shows most of the *stronger* replicate's top peaks in the other replicate.
Weak peaks are always less reproducible, so rank by significance before comparing.

## 6. Differential binding with DiffBind

[DiffBind](https://bioconductor.org/packages/DiffBind/) builds a consensus peak set, counts reads
(input-subtracted), normalises and runs DESeq2 or edgeR. It starts from a sample sheet:

```text
SampleID,Condition,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller
ctrl_1,ctrl,1,results/filtered/ctrl_rep1.GRCh38.bam,in_ctrl,results/filtered/input_ctrl.GRCh38.bam,results/peaks/narrow/ctrl_rep1.GRCh38_peaks.narrowPeak,narrow
...
```

```r
library(DiffBind)
db <- dba(sampleSheet = "diffbind_samples.csv")
db <- dba.count(db, summits = 200)           # re-centre peaks on summits, ±200 bp
db <- dba.normalize(db)                      # default: library-size (full-library) normalisation
db <- dba.contrast(db, design = "~ Condition", contrast = c("Condition", "treat", "ctrl"),
                   minMembers = 2)            # default minMembers = 3: with 2 replicates per group
                                              # no contrast is added unless you lower it
db <- dba.analyze(db)                         # DESeq2 by default; applies blacklist + input greylist.
                                              # If greylisting fails (very sparse input): bGreylist = FALSE
dba.show(db, bContrasts = TRUE)
sites <- dba.report(db, th = 0.05)            # GRanges of differentially bound sites
dba.plotPCA(db, attributes = DBA_CONDITION, label = DBA_ID)
dba.plotMA(db); dba.plotVolcano(db)
```

To practise without data, DiffBind includes a pre-counted example (ER ChIP-seq in breast-cancer
lines with differing tamoxifen response):

```r
data(tamoxifen_counts)                        # loads a counted DBA object named `tamoxifen`
tam <- dba.contrast(tamoxifen, design = "~ Tissue + Condition")
tam <- dba.analyze(tam)
dba.show(tam, bContrasts = TRUE)
```

::: {.callout-note collapse="true"}
## Choices inside differential binding, with pros and cons
- **Library-size vs. RLE normalisation** (`dba.normalize(normalize = DBA_NORM_RLE)`): RLE
  (DESeq2-style) assumes most sites don't change. Library size is safer when a condition gains
  or loses binding globally.
- **Background normalisation** (`background = TRUE`, large bins) is robust to global shifts
  without a spike-in.
- **Spike-in normalisation** (`spikein = TRUE`) is the only way to measure genuinely global
  changes.
- **Consensus rule** (`minOverlap`): peaks present in ≥ 2 samples by default. Stricter rules
  gain confidence but lose condition-specific sites.
:::

## 7. Visualisation

```bash
computeMatrix reference-point --referencePoint center -a 2000 -b 2000 \
    -R differential_sites.bed -S results/bigwig/*.RPGC.bw -o sites.mat.gz
plotHeatmap -m sites.mat.gz -o sites_heatmap.png --kmeans 2
plotProfile -m sites.mat.gz -o sites_profile.png --perGroup
```

Check top sites in IGV with the ChIP **and** input tracks. A real peak is absent from the
input. Compare treated and control replicates at the same y-axis scale.

## 8. Annotation and motifs

```r
library(ChIPseeker); library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
anno <- annotatePeak(sites, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = "org.Hs.eg.db")
plotAnnoBar(anno); plotDistToTSS(anno)
covplot(sites, weightCol = "Fold")            # peak positions along chromosomes
```

For a TF, the expected motif should be strongly enriched and **centred** on peak summits. If
it isn't, suspect the antibody or an indirect-binding mode. Use ±50–100 bp around summits:

```bash
awk 'BEGIN{OFS="\t"}{s=$2+$10; print $1, s-50, s+50, $4}' sample_peaks.narrowPeak > summits100.bed
findMotifsGenome.pl summits100.bed hg38 homer_out/ -size given -p 8
# or MEME-ChIP on summit sequences, which also runs CentriMo to test motif centrality
```

## 9. Integrating with gene expression

Binding near a gene doesn't mean regulation. Test the association instead of assuming it.

```r
# Are differentially expressed genes more often near a differential peak than other genes?
# de_genes / all_genes: gene IDs from the RNA-seq guide; genes_gr: gene GRanges from the same annotation
near <- overlapsAny(promoters(genes_gr, upstream = 50000, downstream = 50000), sites)
tab <- table(DE = names(genes_gr) %in% de_genes, near_peak = near)
fisher.test(tab)
```

Or use GREAT-style region-to-gene enrichment ([rGREAT](https://bioconductor.org/packages/rGREAT/)),
which accounts for how much genomic territory each gene regulates.

## 10. What you can and can't conclude

- **Peak ≠ function.** Many binding events have no measurable effect on nearby genes.
- **Differential peak height ≠ occupancy change** if normalisation is wrong. Global changes
  need spike-ins.
- **Peak number depends on depth, antibody and thresholds.** Don't compare raw peak counts
  across experiments. Compare signal over a common consensus set.
- **Input matters.** Copy-number differences between samples (e.g. cancer lines) create
  false differential peaks unless each ChIP has its own matched input.
- **Broad marks are domains, not peaks.** Summarise signal over genes or domains rather than
  summit windows.

## 11. Extension: comparing species

Comparing binding of the same factor across species (e.g. CTCF in human and macaque;
[Schmidt *et al.* 2012](https://doi.org/10.1016/j.cell.2011.11.058)) adds three problems:
different coordinates, different mappability, and making normalisation fair.

1. **Call peaks per species on its own genome.** Centre on summits and use a fixed width
   (e.g. 50–200 bp for a TF).
2. **liftOver each species' peaks to the other genome** (`-minMatch=0.5` is a common,
   tolerant setting for short intervals). Classify each peak as:
   - *conserved*: lifts over **and** overlaps a peak called independently in the other species;
   - *species-specific*: lifts over, but no peak in the other species (a candidate gain or loss);
   - *unmappable*: doesn't lift over, so absent sequence rather than absent binding.
3. **Count each species' reads on its own coordinates** for a shared set of regions, and test
   with DESeq2 or DiffBind. Include a **random background set** of matched regions. If random
   regions come out "differential", normalisation or mappability is biased, not biology.
4. **Check sensitivity to the counting window** (e.g. ±50 bp vs ±100 bp around summits). A
   conclusion that flips with the window size isn't robust.

::: {.callout-note collapse="true"}
## Pooling two species in one immunoprecipitation
Chromatin from two species can be mixed before the IP, so both experience identical antibody,
washing and library conditions. Reads are then assigned to species computationally
(`pipelines/chip_atac` with two genomes, via BBSplit). Two consequences follow:

- Reads identical in both genomes can't be assigned and are discarded (`ambiguous2=toss`), so
  highly conserved regions lose coverage in *both* species equally.
- Both species' samples from one IP share its technical efficiency. A natural normalisation is
  therefore **one size factor per IP**, estimated on the two species' pooled counts and applied
  to both:

```r
# counts: regions x samples; ip: factor saying which IP each sample came from
pooled <- sapply(split(seq_len(ncol(counts)), ip), function(j) rowSums(counts[, j, drop = FALSE]))
sf_ip  <- DESeq2::estimateSizeFactorsForMatrix(pooled)
sizeFactors(dds) <- sf_ip[as.character(ip)]          # same factor for both species in one IP
dds <- DESeq(dds)                                    # design e.g. ~ ip + species
```
This removes IP-to-IP efficiency differences without letting each species' own binding pattern
set its normalisation.
:::

## 12. Exercises

1. A TF ChIP has FRiP 0.4 % and NSC 1.02. What would you conclude before looking at any
   differential result?
2. You find 5,000 more peaks in treated than control samples, but treated samples were
   sequenced twice as deep. Is binding increased?
3. In a cross-species comparison, 30 % of human peaks don't lift over to macaque. Are these
   human-specific binding sites?

::: {.callout-note collapse="true"}
## Answers
1. Very weak enrichment: likely a failed IP or a poor antibody. Check the fingerprint and
   the motif centrality before investing more time.
2. Not necessarily. Peak calling sensitivity rises with depth. Compare signal over a
   consensus set with proper normalisation, or subsample to equal depth and re-call.
3. No. Unmappable peaks sit on sequence without a clear ortholog (insertions, deletions,
   rearrangements, or poor assembly). Keep them as a separate "structurally different" class.
:::

**Further reading:** HBC's [Investigating chromatin biology using ChIP-seq](https://hbctraining.github.io/Investigating-chromatin-biology-ChIPseq/)
and [Intro-to-peak-analysis](https://hbctraining.github.io/Intro-to-peak-analysis/), the
[DiffBind vignette](https://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf),
and [ChIPseeker](https://doi.org/10.1093/bioinformatics/btv145).
