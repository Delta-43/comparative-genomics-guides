---
title: "ATAC-seq: chromatin accessibility and regulatory elements"
subtitle: "From FASTQ to differential accessibility, annotation, motifs and footprints"
---

::: {.callout-tip}
## Learning objectives
- Plan an ATAC-seq experiment and know the QC metrics that show it worked.
- Process paired-end reads into filtered, Tn5-corrected alignments and peak calls.
- Build a consensus peak set, count reads and test for differential accessibility in R.
- Annotate peaks, look for enriched motifs, and understand what footprinting can add.
- Interpret accessibility changes cautiously, and extend the analysis across species.

**Time:** about 3 hours. **You need:** the set-up from [Getting started](00-getting-started.md).
:::

## 1. What ATAC-seq measures

The Tn5 transposase cuts and tags DNA where chromatin is open: nucleosome-free regions such as
active promoters, enhancers and insulators. Sequencing the tagged fragments gives a genome-wide
map of **accessibility** ([Buenrostro *et al.* 2013](https://doi.org/10.1038/nmeth.2688)).
Fragment lengths carry information too. Short fragments (< 100 bp) come from nucleosome-free
DNA, and fragments of ~200, ~400… bp span one, two… nucleosomes.

**It can tell you:** where regulatory elements are open, and which ones open or close between
conditions. **It can't tell you** by itself which factor binds there, whether an element is
active or poised (combine it with H3K27ac ChIP-seq, [next guide](03-chip-seq.md)), or which gene
it regulates (see the [Hi-C guide](04-hic.md)).

## 2. Experimental design

| Factor | Recommendation |
|---|---|
| Replicates | ≥ 2 biological replicates (ENCODE); ≥ 3 for differential analysis |
| Depth | ≥ 25 M non-duplicate, non-mitochondrial fragments per replicate (ENCODE: 50 M paired-end reads) |
| Input material | 50,000 viable cells is the classic protocol. Dead cells release open, unprotected DNA, so check viability |
| Mitochondrial reads | Expect 5–50 % depending on cell type. The Omni-ATAC protocol ([Corces *et al.* 2017](https://doi.org/10.1038/nmeth.4396)) reduces them |
| Controls | No input control is needed; the peak caller models background from the data |

## 3. Raw-read QC

The same FastQC/MultiQC checks as in the [RNA-seq guide](01-rna-seq.md#3-quality-control-of-raw-reads)
apply. Two ATAC-specific features are normal. Nextera adapter read-through appears in short
fragments (Trim Galore detects and removes it), and the first bases show Tn5's sequence
preference.

## 4. From FASTQ to peaks

```bash
cd pipelines/chip_atac
$EDITOR config.yaml samples.tsv       # set assay: "atac"; group replicates in the `group` column
snakemake --use-conda --cores 16 -n && snakemake --use-conda --cores 16
```

| Step | Tool / setting | Why | Trade-off |
|---|---|---|---|
| Align | bowtie2 `--very-sensitive -X 2000` | `-X 2000` allows fragments up to 2 kb, keeping multi-nucleosome fragments | Slower than the default |
| Keep good pairs | Proper pairs, MAPQ ≥ 30 | Removes multi-mappers, which pile up falsely in repeats | Loses signal in recent duplications and repeats |
| Remove duplicates | `samtools markdup -r` | PCR duplicates inflate peaks from low-complexity libraries | Very deep libraries have some true duplicates |
| Remove mito and contigs | Name filter (`chrM`/`MT`, `_random`, `chrUn`) | mtDNA isn't chromatinised and can take half the reads | None |
| Blacklist | ENCODE blacklist ([Amemiya *et al.* 2019](https://doi.org/10.1038/s41598-019-45839-z)) | Artefact regions with signal in every experiment | Needs a blacklist for your genome build |
| Tn5 shift | `alignmentSieve --ATACshift` (+4 / −5 bp) | Tn5 binds as a dimer and inserts 9 bp apart. The shift centres reads on the cut site | Essential for footprinting; small effect on peaks |
| Peaks | MACS2 `-f BAMPE` | Uses each pair's real fragment extent; no shift model needed | Calls fragment-level peaks, not cut-site peaks |

::: {.callout-note collapse="true"}
## Alternatives worth knowing
- **Cut-site peaks:** treat each read end as an insertion,
  `macs2 callpeak -f BAM --nomodel --shift -50 --extsize 100 --keep-dup all`. This gives sharper
  peaks and is often preferred before footprinting.
- **Genrich** (ATAC mode, handles replicates jointly) and **HMMRATAC** (models nucleosome
  positions) are alternative peak callers.
- **nf-core/atacseq** is a full community pipeline with extensive QC, if you'd rather not
  maintain your own.
:::

**Checkpoints** ([ENCODE ATAC-seq standards](https://www.encodeproject.org/atac-seq/)):

| Metric | Where | Ideal | Acceptable |
|---|---|---|---|
| Alignment rate | bowtie2 log (MultiQC) | > 95 % | > 80 % |
| Non-dup, non-mito fragments | `samtools view -c -f 64` on the filtered BAM | ≥ 25 M | lower is usable, but has less power |
| FRiP (fraction of reads in peaks) | `results/qc/frip/*.txt` | > 0.3 | > 0.2 |
| TSS enrichment (GRCh38) | `results/qc/tss/*.png` | > 7 | 5–7 (mm10: > 15 ideal, 10–15 acceptable) |
| Library complexity | NRF / PBC1 / PBC2 | > 0.9 / > 0.9 / > 3 | |
| Fragment sizes | `results/qc/fragment_size/*.txt` | A nucleosome-free peak (< 100 bp) and a mono-nucleosome peak (~200 bp) | Required |

Plot the fragment-size distribution in R:

```r
library(ggplot2)
fs <- do.call(rbind, lapply(list.files("results/qc/fragment_size", full.names = TRUE), function(f) {
  d <- read.table(f, col.names = c("size", "n"))
  data.frame(sample = sub("\\..*", "", basename(f)), size = d$size, frac = d$n / sum(d$n))
}))
ggplot(subset(fs, size <= 1000), aes(size, frac, colour = sample)) +
  geom_line() + scale_y_log10() + labs(x = "fragment length (bp)", y = "fraction of fragments") + theme_bw()
```

A good library shows a saw-tooth pattern with a ~10 bp period (the DNA helical repeat) below
150 bp, then peaks at ~200 and ~400 bp. A smooth curve with no nucleosome peaks suggests
over-digestion or dead cells.

## 5. A consensus peak set

Differential analysis needs one set of regions shared by all samples. A common, robust recipe is
to use fixed-width peaks around each summit, and keep regions called in at least two samples:

```r
library(GenomicRanges)
read_narrowpeak <- function(f, width = 500) {
  x <- read.table(f, col.names = c("chr","start","end","name","score","strand","fc","p","q","summit"))
  gr <- GRanges(x$chr, IRanges(x$start + x$summit + 1, width = 1))   # summit, 1-based
  resize(gr, width, fix = "center")
}
files <- list.files("results/peaks/narrow", "narrowPeak$", full.names = TRUE)
peaks <- lapply(files, read_narrowpeak)
merged <- reduce(unlist(GRangesList(peaks)))
n_samples <- rowSums(sapply(peaks, function(p) overlapsAny(merged, p)))
consensus <- merged[n_samples >= 2]
length(consensus)
```

::: {.callout-note collapse="true"}
## Alternatives
The pipeline also calls peaks on the pooled BAM of each group (`results/peaks_merged/`): higher
sensitivity, one call per condition. [IDR](https://github.com/nboley/idr) between replicates is
ENCODE's reproducibility standard. [DiffBind](https://bioconductor.org/packages/DiffBind/)
builds consensus peaks, counts and runs DESeq2/edgeR in one package (`dba()`, `dba.count()`,
`dba.analyze()`). It is convenient, but hides the steps shown here.
:::

## 6. Differential accessibility

Count fragments per consensus region, then use DESeq2 exactly as for RNA-seq genes:

```r
library(Rsubread); library(DESeq2)
saf <- data.frame(GeneID = paste0("peak", seq_along(consensus)), Chr = seqnames(consensus),
                  Start = start(consensus), End = end(consensus), Strand = "+")
bams <- list.files("results/filtered", "\\.bam$", full.names = TRUE)
fc <- featureCounts(bams, annot.ext = saf, isPairedEnd = TRUE, nthreads = 8)
counts <- fc$counts
colnames(counts) <- sub("\\..*", "", basename(colnames(counts)))

coldata <- data.frame(condition = factor(c("ctrl","ctrl","treat","treat")), row.names = colnames(counts))
dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "treat", "ctrl"), alpha = 0.05)
summary(res)
vsd <- vst(dds, nsub = min(1000, nrow(dds)))   # nsub: vst needs it lowered for < 1000 regions
plotPCA(vsd, intgroup = "condition")
```

::: {.callout-warning}
## Normalisation assumes most regions don't change
DESeq2's default size factors assume most peaks are unchanged. If a treatment opens or closes
chromatin globally (for example a chromatin-remodeller knockout), that assumption fails. Then
normalise on reads in background bins (`csaw::normFactors` on large windows), or on a
spike-in.
:::

## 7. Visualising accessibility

```r
plotMA(res, ylim = c(-4, 4))
da <- consensus[which(res$padj < 0.05)]
rtracklayer::export(da, "differential_peaks.bed")
```

```bash
# Heatmap and profile of signal around differential peaks (deepTools)
computeMatrix reference-point -R differential_peaks.bed -S results/bigwig/*.RPGC.bw \
    --referencePoint center -a 2000 -b 2000 -o da.mat.gz
plotHeatmap -m da.mat.gz -o da_heatmap.png
```

Always look at a handful of top regions in a genome browser (IGV) with the bigWigs from every
replicate. A believable change is consistent across replicates and isn't a single spike.

## 8. Annotation: where are the peaks?

```r
library(ChIPseeker); library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
anno <- annotatePeak(consensus, TxDb = txdb, tssRegion = c(-1000, 1000), annoDb = "org.Hs.eg.db")
plotAnnoPie(anno)                     # promoter / intron / distal intergenic ...
plotDistToTSS(anno)
head(as.data.frame(anno)[, c("seqnames","start","annotation","SYMBOL","distanceToTSS")])
```

Promoters are usually open in every cell type. Cell-type- and condition-specific signal sits
mostly in **distal** elements (introns, intergenic). A "nearest gene" is only a guess at the
target: enhancers can act over hundreds of kilobases.

## 9. Motifs and footprints

**Motif enrichment** asks which transcription-factor motifs are over-represented in a set of
regions compared with a background (e.g. up-regulated vs all consensus peaks):

```bash
# HOMER: GC-matched background chosen automatically
findMotifsGenome.pl differential_up.bed hg38 homer_up/ -size 200 -p 8
# MEME suite SEA: explicit background sequences
sea --p up.fa --n background.fa --m JASPAR2024_CORE_vertebrates.meme --o sea_up
```

**chromVAR** ([Schep *et al.* 2017](https://doi.org/10.1038/nmeth.4401)) scores per-sample
accessibility deviations for every motif. It's useful for many samples and for single-cell ATAC.

**Footprinting** looks for a dip in Tn5 cuts where a bound factor protects DNA, inside an open
peak. It needs Tn5-shifted reads, deep data (tens of millions of fragments in peaks; deeper is
better) and correction for Tn5's own sequence bias, e.g. [TOBIAS](https://github.com/loosolab/TOBIAS) (`ATACorrect` → `ScoreBigwig` →
`BINDetect`) or [HINT-ATAC](https://reg-gen.readthedocs.io/). Footprints of factors with short
residence times (many nuclear receptors) are often invisible, so a missing footprint doesn't
show that a factor is absent.

## 10. What you can and can't conclude

| Claim | Supported by ATAC-seq alone? |
|---|---|
| "Region X is more accessible in condition B" | Yes, if replicated and normalisation assumptions hold |
| "Region X is an active enhancer" | No. Add H3K27ac (active) / H3K4me1 (primed), or a reporter assay |
| "Factor Y binds more in B" | Suggestive, if Y's motif is enriched and footprinted. Confirm by ChIP/CUT&RUN |
| "Region X regulates gene Z" | No. Needs 3D contact (Hi-C, HiChIP), eQTL or perturbation (CRISPRi) evidence |
| "Accessibility change causes expression change" | Correlation only, unless tested by perturbation |

Common pitfalls: skipping the TSS-enrichment check, comparing samples with very different
FRiP (normalisation then reflects quality, not biology), and reading every open promoter
as condition-specific.

## 11. Extension: comparing species

Peaks from different species live in different coordinate systems. A robust approach:

1. **Call peaks per species** on its own genome, then fix peak width around summits (e.g. 500 bp).
2. **liftOver** every peak set into every other species with UCSC chain files, requiring a
   minimum fraction of bases to map (`liftOver -minMatch=0.5`). Keep only peaks that map to all
   species being compared ("alignable" peaks). Peaks that don't map are usually deletions,
   insertions or rearrangements, not accessibility differences.
3. **Merge** alignable peaks into one universe with coordinates in *every* genome, then count
   each species' reads **on its own genome's coordinates** and join the counts by region ID.
4. Test with DESeq2 per species pair. With a third species as outgroup, call a region
   **gained on one lineage** only if it differs in the same direction from both other species
   (same logic as the [RNA-seq guide](01-rna-seq.md#11-extension-comparing-species)).

```bash
liftOver -minMatch=0.5 human_peaks.bed hg38ToPanTro6.over.chain.gz human_on_chimp.bed human_unmapped.bed
```

```r
# After counting each species on its own coordinates (count matrices with the same row IDs):
counts <- cbind(counts_human[ids, ], counts_chimp[ids, ], counts_macaque[ids, ])
```

::: {.callout-note}
## Example of results
Ciuba *et al.* (2025) applied this to iPSC-derived astrocytes (ATAC-seq: ArrayExpress
E-MTAB-13253; H3K27ac: E-MTAB-13254; H3K4me3: E-MTAB-13255). From ~224,000 regions alignable
across human, chimpanzee and macaque, they classified enhancers (non-promoter, non-H3K4me3
ATAC peaks) as gained or lost on the human lineage. They then linked the gained ones to
up-regulated genes within 500 kb. Human-gained enhancers near up-regulated genes were enriched
for sequence changes in motifs of "stripe" transcription factors, a family of GC-rich-motif
zinc-finger factors (SP, KLF, EGR, ZBTB) that make chromatin accessible for co-binding partners ([Zhao *et al.* 2022](https://doi.org/10.1016/j.molcel.2022.06.029)).
:::

**Testing causality.** Accessibility differences are correlative. Massively parallel reporter
assays (MPRA) test whether a specific sequence change drives activity. Clone each element in
several versions (derived, "ancestralised" at the changed positions, and randomly mutated at
the same number of positions), measure RNA/DNA barcode ratios, and compare the versions with
[MPRAnalyze](https://bioconductor.org/packages/MPRAnalyze/) or
[mpra](https://bioconductor.org/packages/mpra/).

## 12. Exercises

1. Two samples have FRiP 0.35 and 0.12. What happens to their DESeq2 size factors, and how could
   that bias a differential test?
2. Your TSS enrichment is 3. Name two likely causes and what you'd check.
3. Why must reads be Tn5-shifted for footprinting but not necessarily for peak calling?

::: {.callout-note collapse="true"}
## Answers
1. The low-FRiP sample has fewer reads in peaks, so peak-based size factors scale it up. Any
   residual quality difference (a flatter signal) then looks like a biological change. Check
   PCA for a quality axis, and consider background-bin normalisation.
2. Dead or dying cells (free DNA gives uniform background), over-digestion, or a wrong TSS
   annotation (chromosome naming mismatch). Check viability, the fragment-size plot, and that
   the GTF matches the genome.
3. Peaks are hundreds of bp wide, so a 4–5 bp offset barely moves them. Footprints are 10–20 bp
   features, and an unshifted read end sits beside the true cut site, blurring the dip.
:::

**Further reading:** [ENCODE ATAC-seq standards](https://www.encodeproject.org/atac-seq/),
[Yan *et al.* 2020, "From reads to insight: a hitchhiker's guide to ATAC-seq data analysis"](https://doi.org/10.1186/s13059-020-1929-3),
HBC's [Intro-to-peak-analysis](https://hbctraining.github.io/Intro-to-peak-analysis/).
